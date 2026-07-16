# Assembles the real, server-computed analysis results into the exact JSON
# shape the browser-side visualization (www/ptracker_viz.js) expects. Kept
# as pure functions (no Shiny reactives) so they're directly unit-testable.

.IDENTITY_PALETTE_R <- c("#1f8a70", "#d9730d", "#5b6470", "#7a5cff", "#2a9d8f", "#e07a5f", "#3d5a80")

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
#' @return list ready for jsonlite::toJSON(auto_unbox = TRUE)
build_section1_payload <- function(pf_list, iso_key_of, masses, tiers, exon_table, residue_offset = list()) {
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
      tier_y = tt$tier_y
    )
  })

  list(
    proteoforms = proteoforms,
    axis_length = if (!is.null(axis)) axis$axis_length else max(vapply(pf_list, function(p) nchar(p$sequence), integer(1)))
  )
}

#' Build the JSON-ready payload for section 2 (confounding-protein search):
#' the target's own ladder/tiers/exon structure plus real confounder masses
#' and envelopes.
#'
#' @param target_pf proteoform object (search target)
#' @param target_mass total mass (Da)
#' @param target_tiers output of compute_confounder_tiers()
#' @param target_exon_table this transcript's own exon rows (single
#'   transcript, no cross-transcript alignment needed)
#' @param confounder_search result of search_confounding_proteins_real()
#' @param confounder_envs named list, confounder id -> data.frame(z, mz,
#'   relative_intensity)
#' @return list ready for jsonlite::toJSON(auto_unbox = TRUE)
build_section2_payload <- function(target_pf, target_mass, target_tiers, target_exon_table,
                                    confounder_search, confounder_envs) {
  exon_blocks <- if (!is.null(target_exon_table) && nrow(target_exon_table) > 0) {
    te <- target_exon_table[order(target_exon_table$residue_start), ]
    lapply(seq_len(nrow(te)), function(r) list(start = te$residue_start[r], end = te$residue_end[r]))
  } else list()

  confounders <- list()
  if (!is.null(confounder_search) && nrow(confounder_search$candidates) > 0) {
    cands <- confounder_search$candidates
    confounders <- lapply(seq_len(nrow(cands)), function(r) {
      cid <- cands$id[r]
      env <- confounder_envs[[cid]]
      list(
        id = cid, mass = round(cands$mass[r], 2), length = cands$length[r],
        env = if (!is.null(env)) Map(function(z, mz, ri) list(z = z, mz = round(mz, 2), rel = round(ri, 4)),
                                      env$z, env$mz, env$relative_intensity) else list()
      )
    })
  }

  list(
    target = list(
      id = target_pf$id, mass = round(target_mass, 2), len = nchar(target_pf$sequence),
      sequence = target_pf$sequence,
      ptms = lapply(target_pf$ptms, function(p) list(pos = p$site, name = p$name %||% "PTM")),
      b_mass = round(target_tiers$ladder$b_mass, 2),
      y_mass = round(target_tiers$ladder$y_mass, 2),
      propensity = round(target_tiers$propensity, 2),
      exon_blocks = exon_blocks,
      tier_b = target_tiers$tier_b,
      tier_y = target_tiers$tier_y
    ),
    confounders = confounders,
    window_da = if (!is.null(confounder_search)) round(confounder_search$window_da, 3) else NA,
    best_charge_state = if (!is.null(confounder_search)) confounder_search$best_charge_state else NA
  )
}
