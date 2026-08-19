# Assembles the real, server-computed analysis results into the exact JSON
# shape the browser-side visualization (www/ptracker_viz.js) expects. Kept
# as pure functions (no Shiny reactives) so they're directly unit-testable.

.IDENTITY_PALETTE_R <- c("#1f8a70", "#d9730d", "#5b6470", "#7a5cff", "#2a9d8f", "#e07a5f", "#3d5a80")

#' Tally MS2 fragment-ion tiers (common/partial/unique/neutral) across
#' both ion series for the stat-tile summary panel. "Informative" here
#' means "unique" specifically -- a fragment mass no OTHER checked
#' proteoform/confounder shares, i.e. one that on its own would confirm
#' this proteoform's identity if observed. "partial" bonds are shared with
#' SOME but not all others, so they're not double-counted as either
#' fully-informative or fully-uninformative.
#'
#' @param tier_b,tier_y character vectors ("common"/"partial"/"unique"/"neutral")
#' @return list(total, unique, partial, common, neutral)
tally_ms2_tiers <- function(tier_b, tier_y) {
  all_tiers <- c(tier_b, tier_y)
  list(
    total = length(all_tiers),
    unique = sum(all_tiers == "unique"),
    partial = sum(all_tiers == "partial"),
    common = sum(all_tiers == "common"),
    neutral = sum(all_tiers == "neutral")
  )
}

#' Categorical stringency bucket for a propensity score, mirroring the
#' client's own filter thresholds (www/ptracker_viz.js TIER_THRESHOLDS) --
#' used for the peak-data TSV export (R/viz_json.R's flatten_*_peaks()),
#' NOT for any chart rendering decision itself (those stay purely
#' threshold-comparison client-side, this is just a readable label for a
#' downloaded file). The two scoring modes are on different numeric scales
#' (see R/fragmentation_propensity.R vs R/fragmentation_propensity_rf.R) so
#' scoring_mode must always be passed through, never assumed.
#'
#' @param score numeric vector of propensity scores
#' @param scoring_mode "glm" or "rf"
#' @return character vector, one of "Baseline"/"Elevated"/"High"
propensity_group <- function(score, scoring_mode = c("glm", "rf")) {
  scoring_mode <- match.arg(scoring_mode)
  thresholds <- if (scoring_mode == "rf") {
    c(elevated = RF_ELEVATED_THRESHOLD, high = RF_HIGH_THRESHOLD)
  } else {
    c(elevated = GLM_ELEVATED_THRESHOLD, high = GLM_HIGH_THRESHOLD)
  }
  ifelse(score >= thresholds[["high"]], "High",
    ifelse(score > thresholds[["elevated"]], "Elevated", "Baseline"))
}

#' Converts predict_ms1_peaks()'s R-native output (R/isotope_envelope.R)
#' into the nested list-of-lists shape jsonlite::toJSON(auto_unbox = TRUE)
#' needs for the MS1 chart: one entry per populated charge state, each
#' carrying its own array of {mz, rel} points -- either the real discrete
#' isotope comb (resolved = TRUE) or a handful of samples along a smooth
#' envelope curve (resolved = FALSE), decided per charge state against the
#' instrument's own resolving power. See predict_ms1_peaks()'s doc comment
#' for the physics.
#'
#' @param ms1_peaks output of predict_ms1_peaks()
#' @return list of list(z=, resolved=, points = list(list(mz=, rel=), ...))
ms1_peaks_json <- function(ms1_peaks) {
  lapply(ms1_peaks$charge_states, function(cs) {
    list(
      z = cs$z,
      resolved = cs$resolved,
      points = Map(function(mz, rel) list(mz = round(mz, 4), rel = round(rel, 4)),
                    cs$points$mz, cs$points$rel)
    )
  })
}


#' Build the JSON-ready payload for section 1 (relevant-proteoform
#' comparison): one entry per checked proteoform, each with its own ladder,
#' propensity, tier colors, PTM markers, and its mapping onto the shared
#' exon axis (so proteoforms from the same gene visually align).
#'
#' @param pf_list named list of proteoform objects (checked proteoforms)
#' @param iso_key_of named character vector, proteoform id -> its underlying
#'   transcript id (PTM variants share their bare parent's transcript id)
#' @param masses named numeric vector, proteoform id -> total mass
#' @param tiers output of compute_ladder_tiers(pf_list, ...)
#' @param exon_table real per-transcript exon structure (from
#'   get_exon_structure()), used to build the shared axis
#' @param residue_offset named list/vector, proteoform id -> integer offset
#'   into its iso_key parent's full-length numbering (0 for an ordinary
#'   full-length proteoform; `start - 1` for a middle-down peptide, whose own
#'   sequence is 1-based in ITS OWN coordinates but needs shifting back onto
#'   its parent transcript's shared axis before lookup). exon_blocks are left
#'   alone -- axis space is shared across the whole transcript regardless of
#'   which sub-range of it any one row happens to display.
#' @param scoring_mode "glm" or "rf" -- which scoring mode produced `tiers`
#'   (see compute_ladder_tiers()); sent through so the client knows which
#'   tier thresholds apply (www/ptracker_viz.js, RF_TIER_THRESHOLDS_JS-style
#'   constants -- the two modes are NOT on the same numeric scale).
#' @param ms1_stats optional list(total_peaks, crowded_peaks, clean_peaks)
#'   from envelope_crowding_check() (server.R), for the stat-tile panel;
#'   NULL if not computed (e.g. fewer than 2 proteoforms checked)
#' @param gene_symbol optional gene symbol covering this whole checked set
#'   (Option 1/2/3 all resolve to proteoforms of one gene at a time
#'   currently); carried through only for the peak-data TSV export
#'   (flatten_section1_peaks()), NA if unknown
#' @return list ready for jsonlite::toJSON(auto_unbox = TRUE)
build_section1_payload <- function(pf_list, iso_key_of, masses, tiers, exon_table, residue_offset = list(),
                                    scoring_mode = "glm", ms1_stats = NULL, gene_symbol = NA_character_) {
  ids <- names(pf_list)
  unique_iso_keys <- unique(iso_key_of[ids])
  axis <- if (!is.null(exon_table) && length(unique_iso_keys) > 0) {
    build_shared_exon_axis(exon_table, unique_iso_keys)
  } else NULL

  proteoforms <- lapply(seq_along(ids), function(i) {
    id <- ids[i]
    pf <- pf_list[[id]]
    tt <- tiers[[id]]
    iso_key <- iso_key_of[[id]]
    offset <- residue_offset[[id]] %||% 0L
    segs <- if (!is.null(axis)) axis$segments[[iso_key]] else NULL
    axis_pos <- if (!is.null(segs)) {
      vapply(seq_len(nchar(pf$sequence) - 1), function(p) {
        v <- map_local_to_axis(segs, p + offset)
        if (is.na(v)) -1L else as.integer(v)
      }, integer(1))
    } else {
      seq_len(nchar(pf$sequence) - 1)
    }
    exon_blocks <- if (!is.null(segs)) {
      lapply(seq_len(nrow(segs)), function(r) list(start = segs$axis_start[r], end = segs$axis_end[r]))
    } else list()

    list(
      id = id,
      label = pf$id,
      transcript_id = iso_key %||% NA_character_,
      color = .IDENTITY_PALETTE_R[((i - 1) %% length(.IDENTITY_PALETTE_R)) + 1],
      mass = round(masses[[id]], 2),
      len = nchar(pf$sequence),
      sequence = pf$sequence,
      ptms = lapply(pf$ptms, function(p) list(pos = p$site, axis_pos = if (!is.null(segs)) {
        v <- map_local_to_axis(segs, p$site + offset); if (is.na(v)) p$site else as.integer(v)
      } else p$site, name = p$name %||% "PTM")),
      b_mass = round(tt$ladder$b_mass, 2),
      y_mass = round(tt$ladder$y_mass, 2),
      propensity = round(tt$propensity, 2),
      axis_pos = axis_pos,
      exon_blocks = exon_blocks,
      tier_b = tt$tier_b,
      tier_y = tt$tier_y,
      ms2_stats = tally_ms2_tiers(tt$tier_b, tt$tier_y)
    )
  })

  list(
    proteoforms = proteoforms,
    gene_symbol = gene_symbol,
    axis_length = if (!is.null(axis)) axis$axis_length else max(vapply(pf_list, function(p) nchar(p$sequence), integer(1))),
    scoring_mode = scoring_mode,
    ms1_stats = ms1_stats
  )
}

#' Build the JSON-ready payload for section 2 (confounding-protein search):
#' the target's own ladder/tiers/exon structure PLUS every selected
#' confounder's own ladder/tiers -- so the MS2 fragment ladder can show every
#' peak, target and confounders alike, the same way the MS1 chart already
#' does (renderSection2() in www/ptracker_viz.js already draws one MS1 curve
#' per confounder; the MS2 ladder used to draw only the target's row, which
#' is why this now takes a combined `tiers` covering both instead of a
#' target-only `target_tiers`).
#'
#' @param target_pf proteoform object (search target)
#' @param target_mass total mass (Da)
#' @param tiers output of compute_ladder_tiers() run on the COMBINED
#'   named list c(setNames(list(target_pf), target_id), confounder_pfs) --
#'   gives target and every confounder N-way tiers (common/partial/unique)
#'   against that whole set, using the exact same tiering semantics section 1
#'   uses for its checked-proteoform comparison (a confounder tiered
#'   "unique" here means unique against the target AND every other selected
#'   confounder, not just the target alone)
#' @param target_id the key `target_pf` is stored under in `tiers`
#' @param target_exon_table this transcript's own exon rows (single
#'   transcript, no cross-transcript alignment needed)
#' @param confounder_search result of search_confounding_proteins_real()
#' @param confounder_envs named list, confounder id -> predict_ms1_peaks()
#'   output (R/isotope_envelope.R)
#' @param scoring_mode "glm" or "rf" -- see build_section1_payload()
#' @param ms1_stats optional list(total_peaks, colliding_peaks, clean_peaks)
#'   -- how many of the TARGET's own charge-state peaks collide with at
#'   least one real confounder, from confounder_search$mz_collision_detail
#'   (server.R); NULL if not computed
#' @param gene_symbol optional gene symbol for the TARGET only (each
#'   confounder gets its own, from confounder_search$candidates$gene_symbol
#'   when present); carried through only for the peak-data TSV export
#'   (flatten_section2_peaks())
#' @return list ready for jsonlite::toJSON(auto_unbox = TRUE)
build_section2_payload <- function(target_pf, target_mass, tiers, target_id, target_exon_table,
                                    confounder_search, confounder_envs, scoring_mode = "glm", ms1_stats = NULL,
                                    gene_symbol = NA_character_) {
  target_tiers <- tiers[[target_id]]
  exon_blocks <- if (!is.null(target_exon_table) && nrow(target_exon_table) > 0) {
    te <- target_exon_table[order(target_exon_table$residue_start), ]
    lapply(seq_len(nrow(te)), function(r) list(start = te$residue_start[r], end = te$residue_end[r]))
  } else list()

  confounders <- list()
  if (!is.null(confounder_search) && nrow(confounder_search$candidates) > 0) {
    cands <- confounder_search$candidates
    has_gene_col <- "gene_symbol" %in% names(cands)
    confounders <- lapply(seq_len(nrow(cands)), function(r) {
      cid <- cands$id[r]
      env <- confounder_envs[[cid]]
      tt <- tiers[[cid]]
      list(
        id = cid, mass = round(cands$mass[r], 2), len = cands$length[r],
        sequence = if ("sequence" %in% names(cands)) cands$sequence[r] else NA_character_,
        gene_symbol = if (has_gene_col) (cands$gene_symbol[r] %||% NA_character_) else NA_character_,
        found_via = cands$found_via[r] %||% "mass",
        n_colliding_peaks = cands$n_colliding_peaks[r] %||% 0L,
        env = if (!is.null(env)) ms1_peaks_json(env) else list(),
        b_mass = if (!is.null(tt)) round(tt$ladder$b_mass, 2) else list(),
        y_mass = if (!is.null(tt)) round(tt$ladder$y_mass, 2) else list(),
        propensity = if (!is.null(tt)) round(tt$propensity, 2) else list(),
        tier_b = if (!is.null(tt)) tt$tier_b else list(),
        tier_y = if (!is.null(tt)) tt$tier_y else list(),
        ms2_stats = if (!is.null(tt)) tally_ms2_tiers(tt$tier_b, tt$tier_y) else NULL
      )
    })
  }

  list(
    target = list(
      id = target_pf$id, gene_symbol = gene_symbol, mass = round(target_mass, 2), len = nchar(target_pf$sequence),
      sequence = target_pf$sequence,
      ptms = lapply(target_pf$ptms, function(p) list(pos = p$site, name = p$name %||% "PTM")),
      b_mass = round(target_tiers$ladder$b_mass, 2),
      y_mass = round(target_tiers$ladder$y_mass, 2),
      propensity = round(target_tiers$propensity, 2),
      exon_blocks = exon_blocks,
      tier_b = target_tiers$tier_b,
      tier_y = target_tiers$tier_y,
      ms2_stats = tally_ms2_tiers(target_tiers$tier_b, target_tiers$tier_y)
    ),
    confounders = confounders,
    window_da = if (!is.null(confounder_search)) round(confounder_search$window_da, 3) else NA,
    best_charge_state = if (!is.null(confounder_search)) confounder_search$best_charge_state else NA,
    scoring_mode = scoring_mode,
    ms1_stats = ms1_stats,
    n_candidates_before_cap = if (!is.null(confounder_search)) confounder_search$n_candidates_before_cap %||% length(confounders) else NA
  )
}

#' Flatten section 1's payload into a long-format data.frame -- one row per
#' MS1 isotope/envelope point and one row per MS2 backbone bond (split into
#' its b-ion and y-ion rows), for the "Download peak data" button. Operates
#' on the exact same R list build_section1_payload() returns (before
#' jsonlite::toJSON), so what's downloaded always matches what's drawn.
#'
#' @param payload1 output of build_section1_payload(), with `env` already
#'   attached per proteoform (server.R does this after the initial build)
#' @return data.frame, columns: proteoform_id, label, transcript_id,
#'   gene_id, mass_da, len, peak_type (MS1/MS2), charge_state, mz,
#'   relative_intensity, resolved, ion_type, position, fragment_mass_da,
#'   propensity_score, propensity_group, uniqueness_tier, scoring_mode
#'   (NA on MS1 rows -- MS1 peaks don't depend on scoring mode at all, only
#'   the MS2 propensity/tier columns do)
flatten_section1_peaks <- function(payload1) {
  gene_id <- payload1$gene_symbol %||% NA_character_
  scoring_mode <- payload1$scoring_mode %||% "glm"
  rows <- lapply(payload1$proteoforms, function(pf) .flatten_one_proteoform_peaks(pf, gene_id, scoring_mode))
  do.call(rbind, rows)
}

#' Same idea as flatten_section1_peaks() but for section 2: the target plus
#' every confounder currently in the payload (i.e. whichever ones the user
#' left checked in "Compare selected confounders"), with an extra `role`
#' column ("target"/"confounder") and the confounder-search-specific
#' `found_via`/`n_colliding_peaks` columns (NA for target rows, since those
#' describe how a CANDIDATE was found relative to the target, not a
#' property of the target itself).
#'
#' @param payload2 output of build_section2_payload()
#' @return data.frame, same columns as flatten_section1_peaks() plus role,
#'   found_via, n_colliding_peaks
flatten_section2_peaks <- function(payload2) {
  scoring_mode <- payload2$scoring_mode %||% "glm"
  t <- payload2$target
  target_row <- .flatten_one_proteoform_peaks(
    c(t, list(transcript_id = NA_character_)), t$gene_symbol %||% NA_character_, scoring_mode
  )
  target_row$role <- "target"
  target_row$found_via <- NA_character_
  target_row$n_colliding_peaks <- NA_integer_

  conf_rows <- lapply(payload2$confounders, function(cand) {
    r <- .flatten_one_proteoform_peaks(
      c(cand, list(label = cand$id, transcript_id = NA_character_)),
      cand$gene_symbol %||% NA_character_, scoring_mode
    )
    r$role <- "confounder"
    r$found_via <- cand$found_via %||% NA_character_
    r$n_colliding_peaks <- cand$n_colliding_peaks %||% NA_integer_
    r
  })
  do.call(rbind, c(list(target_row), conf_rows))
}

#' Shared row-builder for one proteoform-shaped list (as found in either
#' payload's `proteoforms`/`target`/`confounders` entries) -- long-format:
#' one row per MS1 envelope point, one row per MS2 bond per ion type (b/y).
#' Internal helper, not exported/documented for direct use.
.flatten_one_proteoform_peaks <- function(pf, gene_id, scoring_mode) {
  # scoring_mode deliberately NOT included here -- it only applies to MS2
  # rows (propensity_score/propensity_group/uniqueness_tier all depend on
  # it; MS1 isotope points don't depend on it at all, see compute_ladder_
  # tiers()'s doc comment). Set to NA for MS1 rows and the real value for
  # MS2 rows below, rather than stamping every row with a mode that had no
  # bearing on half of them.
  base <- data.frame(
    proteoform_id = pf$id, label = pf$label %||% pf$id, transcript_id = pf$transcript_id %||% NA_character_,
    gene_id = gene_id %||% NA_character_, mass_da = pf$mass %||% NA_real_, len = pf$len %||% NA_integer_,
    stringsAsFactors = FALSE
  )
  empty_extra <- data.frame(
    peak_type = character(0), charge_state = integer(0), mz = numeric(0), relative_intensity = numeric(0),
    resolved = logical(0), ion_type = character(0), position = integer(0), fragment_mass_da = numeric(0),
    propensity_score = numeric(0), propensity_group = character(0), uniqueness_tier = character(0),
    scoring_mode = character(0), stringsAsFactors = FALSE
  )

  ms1_rows <- lapply(pf$env %||% list(), function(cs) {
    pts <- cs$points %||% list()
    if (length(pts) == 0) return(NULL)
    data.frame(
      peak_type = "MS1", charge_state = cs$z, mz = vapply(pts, function(p) p$mz, numeric(1)),
      relative_intensity = vapply(pts, function(p) p$rel, numeric(1)), resolved = isTRUE(cs$resolved),
      ion_type = NA_character_, position = NA_integer_, fragment_mass_da = NA_real_,
      propensity_score = NA_real_, propensity_group = NA_character_, uniqueness_tier = NA_character_,
      scoring_mode = NA_character_, stringsAsFactors = FALSE
    )
  })
  ms1_df <- if (length(ms1_rows) > 0) do.call(rbind, ms1_rows) else empty_extra

  n_bonds <- length(pf$b_mass %||% list())
  ms2_df <- if (n_bonds > 0) {
    grp <- propensity_group(unlist(pf$propensity), scoring_mode)
    rbind(
      data.frame(
        peak_type = "MS2", charge_state = NA_integer_, mz = NA_real_, relative_intensity = NA_real_,
        resolved = NA, ion_type = "b", position = seq_len(n_bonds), fragment_mass_da = unlist(pf$b_mass),
        propensity_score = unlist(pf$propensity), propensity_group = grp, uniqueness_tier = unlist(pf$tier_b),
        scoring_mode = scoring_mode, stringsAsFactors = FALSE
      ),
      data.frame(
        peak_type = "MS2", charge_state = NA_integer_, mz = NA_real_, relative_intensity = NA_real_,
        resolved = NA, ion_type = "y", position = seq_len(n_bonds), fragment_mass_da = unlist(pf$y_mass),
        propensity_score = unlist(pf$propensity), propensity_group = grp, uniqueness_tier = unlist(pf$tier_y),
        scoring_mode = scoring_mode, stringsAsFactors = FALSE
      )
    )
  } else empty_extra

  cbind(base, rbind(ms1_df, ms2_df))
}
