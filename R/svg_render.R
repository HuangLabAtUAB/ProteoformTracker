# Server-rendered (static, non-interactive) SVG views for the real Shiny
# app -- MS1 charge-envelope overlay and MS2 fragment-ladder tier view.
#
# This is a first real pass: it reuses the exact tier-coloring logic already
# validated in the standalone prototype (common/partial/unique via mass
# collision), but renders one static SVG per Run click rather than the
# prototype's drag/zoom/click JS. Each proteoform is plotted against its own
# residue axis (0..length) rather than a shared cross-isoform exon axis --
# aligning arbitrary genes' isoforms onto one axis needs real per-gene exon
# matching (feasible via reference_exon_index, deferred as a follow-up).

.TIER_COLOR <- c(common = "#eda100", partial = "#8952e0", unique = "#2a78d6", neutral = "#8a8a86")
.PTM_COLOR_DEFAULT <- "#e0433d"
.IDENTITY_PALETTE <- c("#1f8a70", "#d9730d", "#5b6470", "#7a5cff", "#2a9d8f", "#e07a5f", "#3d5a80")

#' For each proteoform in pf_list, compute per-bond tiers (common / partial /
#' unique) by checking its own fragment ladder against every *other*
#' proteoform in the set, reusing fragment_mass_collision_check() --
#' symmetric N-way comparison via K one-target-vs-(K-1)-candidates calls.
#'
#' @param pf_list named list of proteoform objects (id -> proteoform)
#' @param r_ref,mz_ref,safety_margin passed through to fragment_mass_collision_check()
#' @param scoring_mode "glm" (default, calibrated + length-aware) or "rf"
#'   (length-free ranking mode, see R/fragmentation_propensity_rf.R -- use
#'   when a proteoform's own length is suppressing every bond below the
#'   GLM's tier thresholds and within-proteoform ranking is what's needed)
#' @return named list, id -> list(ladder, tier_b, tier_y) where tier_b/tier_y
#'   are character vectors ("common"/"partial"/"unique"/"neutral"), one per
#'   cleavage position
compute_ladder_tiers <- function(pf_list, r_ref = 120000, mz_ref = 200, safety_margin = DEFAULT_SAFETY_MARGIN,
                                  mode = "denatured", scoring_mode = c("glm", "rf")) {
  scoring_mode <- match.arg(scoring_mode)
  ids <- names(pf_list)
  result <- list()
  for (id in ids) {
    others <- pf_list[setdiff(ids, id)]
    ladder <- generate_fragment_ladder(pf_list[[id]])
    prop <- if (scoring_mode == "rf") {
      fragmentation_propensity_rf(pf_list[[id]])$propensity_score
    } else {
      fragmentation_propensity(pf_list[[id]], mode = mode, method = "HCD")$propensity_score
    }
    n_bonds <- length(ladder$b_mass)
    if (length(others) == 0) {
      result[[id]] <- list(ladder = ladder, propensity = prop, tier_b = rep("neutral", n_bonds), tier_y = rep("neutral", n_bonds))
      next
    }
    cc <- fragment_mass_collision_check(pf_list[[id]], others, r_ref = r_ref, mz_ref = mz_ref, safety_margin = safety_margin)
    n_others <- length(others)
    match_b <- Reduce(`+`, lapply(cc$per_candidate, function(d) as.integer(d$matched[d$ion_type == "b"])))
    match_y <- Reduce(`+`, lapply(cc$per_candidate, function(d) as.integer(d$matched[d$ion_type == "y"])))
    tier_of <- function(mc) ifelse(mc == n_others, "common", ifelse(mc == 0, "unique", "partial"))
    result[[id]] <- list(ladder = ladder, propensity = prop, tier_b = tier_of(match_b), tier_y = tier_of(match_y))
  }
  result
}

#' Tiers for the section-2 confounder view: same flat id -> list(ladder,
#' propensity, tier_b, tier_y) shape as compute_ladder_tiers() (drop-in
#' replacement at the call site), but deliberately NOT the same algorithm.
#' compute_ladder_tiers() is full N-way (every member's own bonds checked
#' against every OTHER member) -- fine for section 1's typically-small
#' checked set, but confirmed directly to take 1-3+ minutes for a target +
#' 30 confounders (the search step's own max_candidates cap), since
#' fragment_mass_collision_check() is O(bonds x bonds) per pair and N-way
#' means O(N) such calls each against O(N) others -- ~930 pairwise ladder
#' comparisons at N=31 vs. the 30 this function does instead:
#'   - the TARGET's own tier_b/tier_y is still multi-way (common only if a
#'     bond matches EVERY shown confounder, unique if none, partial
#'     otherwise) -- this is the question that actually matters ("is this
#'     target bond confounded by anything currently shown") and stays cheap:
#'     one batched fragment_mass_collision_check(target, all confounders) call.
#'   - each CONFOUNDER's own tier_b/tier_y is pairwise against the target
#'     ONLY (binary common/unique, no "partial" -- confounder-vs-confounder
#'     comparison isn't the question this view answers), one cheap call per
#'     confounder rather than one call per confounder against every other
#'     confounder too.
#' Net effect: O(N) pairwise comparisons total instead of O(N^2).
#'
#' @param target_pf,target_id proteoform object and the id to key its own
#'   entry under in the returned list (matches whatever key build_section2_
#'   payload()'s `tiers[[target_id]]` lookup expects)
#' @param confounder_pfs named list of confounder proteoform objects
#' @param scoring_mode "glm" or "rf" -- see compute_ladder_tiers()
#' @return named list, target_id/confounder id -> list(ladder, propensity,
#'   tier_b, tier_y) -- same shape compute_ladder_tiers() returns
compute_confounder_tiers <- function(target_pf, target_id, confounder_pfs, r_ref = 120000, mz_ref = 200,
                                      safety_margin = DEFAULT_SAFETY_MARGIN, mode = "denatured",
                                      scoring_mode = c("glm", "rf")) {
  scoring_mode <- match.arg(scoring_mode)
  score_of <- function(pf) {
    if (scoring_mode == "rf") fragmentation_propensity_rf(pf)$propensity_score
    else fragmentation_propensity(pf, mode = mode, method = "HCD")$propensity_score
  }

  target_ladder <- generate_fragment_ladder(target_pf)
  target_prop <- score_of(target_pf)
  n_bonds_t <- length(target_ladder$b_mass)

  if (length(confounder_pfs) == 0) {
    result <- list(list(ladder = target_ladder, propensity = target_prop,
                         tier_b = rep("neutral", n_bonds_t), tier_y = rep("neutral", n_bonds_t)))
    names(result) <- target_id
    return(result)
  }

  # Target vs. ALL confounders, one batched (cheap) call -- unchanged from
  # the original compute_confounder_tiers()'s approach.
  cc <- fragment_mass_collision_check(target_pf, confounder_pfs, r_ref = r_ref, mz_ref = mz_ref, safety_margin = safety_margin)
  n_others <- length(confounder_pfs)
  match_b <- Reduce(`+`, lapply(cc$per_candidate, function(d) as.integer(d$matched[d$ion_type == "b"])))
  match_y <- Reduce(`+`, lapply(cc$per_candidate, function(d) as.integer(d$matched[d$ion_type == "y"])))
  tier_of_multi <- function(mc) ifelse(mc == n_others, "common", ifelse(mc == 0, "unique", "partial"))
  target_entry <- list(ladder = target_ladder, propensity = target_prop,
                        tier_b = tier_of_multi(match_b), tier_y = tier_of_multi(match_y))
  target_all_masses <- c(target_ladder$b_mass, target_ladder$y_mass)

  # Each confounder vs. the target ONLY -- pairwise, binary. Computed
  # directly here (mirroring fragment_mass_collision_check()'s own matching
  # logic) rather than by calling it a second time per confounder: that
  # function recomputes generate_fragment_ladder() for BOTH arguments on
  # every call, so calling it once per confounder would regenerate the
  # TARGET's own ladder -- a real pyteomics/reticulate round-trip -- 30
  # extra times over. Confirmed directly: that redundancy alone was the
  # dominant cost at a 30-confounder candidate list (turned a ~10s job into
  # 35+s and climbing). target_all_masses (above) is already computed once;
  # reusing it here keeps this loop to one new pyteomics call per confounder
  # (its own ladder, unavoidable) plus cheap vectorized R arithmetic.
  conf_entries <- lapply(confounder_pfs, function(cpf) {
    cladder <- generate_fragment_ladder(cpf)
    cprop <- score_of(cpf)
    tol_b <- safety_margin * fwhm_mass(cladder$b_mass, mz_for_charge(cladder$b_mass, 1), r_ref, mz_ref)
    tol_y <- safety_margin * fwhm_mass(cladder$y_mass, mz_for_charge(cladder$y_mass, 1), r_ref, mz_ref)
    matched_b <- vapply(seq_along(cladder$b_mass), function(i) any(abs(target_all_masses - cladder$b_mass[i]) <= tol_b[i]), logical(1))
    matched_y <- vapply(seq_along(cladder$y_mass), function(i) any(abs(target_all_masses - cladder$y_mass[i]) <= tol_y[i]), logical(1))
    list(
      ladder = cladder, propensity = cprop,
      tier_b = ifelse(matched_b, "common", "unique"),
      tier_y = ifelse(matched_y, "common", "unique")
    )
  })

  result <- c(setNames(list(target_entry), target_id), conf_entries)
  result
}


#' Cheap (no isotope-pattern) PAIRWISE MS2 fragment-ion overlap between the
#' target and EACH confounder candidate independently -- one number per
#' candidate, target vs. that candidate alone, regardless of any other
#' candidate. Used at the confounder-SEARCH step (server.R) to help the
#' user decide which candidates are worth including in the full multi-way
#' comparison BEFORE paying the isotope-envelope cost of "Compare selected
#' confounders": fragment_mass_collision_check() only needs each candidate's
#' own deterministic fragment ladder (arithmetic + one batched pyteomics
#' call per candidate), not an isotope pattern, so this is affordable for
#' the whole candidate list up front.
#'
#' Deliberately NOT the same number as compute_confounder_tiers()'s
#' "unique"/"common" tiers: those are inherently SET-dependent (a fragment
#' only counts as "unique" if NO confounder in the currently-selected set
#' collides with it), so they can't be known until the user has actually
#' picked which confounders to compare against. This is just "does THIS one
#' candidate, alone, share a fragment with the target" -- a simpler,
#' selection-independent number that's still useful for spotting which
#' candidates barely overlap the target at MS2 (safe to deselect if you
#' want a smaller comparison) vs. which ones matter most.
#'
#' @param target_pf proteoform object
#' @param candidate_pfs named list of candidate proteoform objects
#' @param r_ref,mz_ref,safety_margin passed to fragment_mass_collision_check()
#' @return named integer vector, candidate id -> number of the target's b/y
#'   ions (out of 2*(n-1) total) that collide with this ONE candidate's own ladder
confounder_ms2_overlap_counts <- function(target_pf, candidate_pfs, r_ref = 120000, mz_ref = 200,
                                           safety_margin = DEFAULT_SAFETY_MARGIN) {
  if (length(candidate_pfs) == 0) return(integer(0))
  cc <- fragment_mass_collision_check(target_pf, candidate_pfs, r_ref = r_ref, mz_ref = mz_ref, safety_margin = safety_margin)
  counts <- vapply(cc$per_candidate, function(d) sum(d$matched), integer(1))
  names(counts) <- names(candidate_pfs)
  counts
}

#' Render an MS1 charge-state envelope overlay as a static SVG (one stick
#' plot per proteoform), real relative-intensity heights from
#' predict_charge_envelope().
#'
#' @param pf_list named list of proteoform objects
#' @param masses named numeric vector, id -> total mass (Da), same names as pf_list
#' @param mode "denatured" or "native"
#' @return character(1) SVG markup (viewBox 0 0 640 130)
render_ms1_svg <- function(pf_list, masses, mode = "denatured") {
  ids <- names(pf_list)
  if (length(ids) == 0) return("<svg viewBox='0 0 640 40'><text x='10' y='20'>No proteoforms checked.</text></svg>")

  envs <- lapply(ids, function(id) predict_charge_envelope(pf_list[[id]]$sequence, masses[[id]], mode = mode))
  names(envs) <- ids
  all_mz <- unlist(lapply(envs, function(e) e$mz))
  mz_min <- min(all_mz) * 0.95
  mz_max <- max(all_mz) * 1.05
  mzx <- function(mz) 40 + (mz - mz_min) / (mz_max - mz_min) * 560
  base_y <- 86

  els <- sprintf('<line x1="40" y1="%d" x2="600" y2="%d" stroke="#ccc" stroke-width="1"/>', base_y, base_y)
  for (v in c(0, 0.5, 1)) {
    y <- base_y - v * 70
    els <- c(els, sprintf('<line x1="36" y1="%.1f" x2="600" y2="%.1f" stroke="#ddd" stroke-width="0.5" stroke-dasharray="2,2"/><text x="4" y="%.1f" font-size="9" fill="#666">%.1f</text>', y, y, y + 3, v))
  }
  for (i in seq_along(ids)) {
    id <- ids[i]
    color <- .IDENTITY_PALETTE[((i - 1) %% length(.IDENTITY_PALETTE)) + 1]
    env <- envs[[id]]
    for (r in seq_len(nrow(env))) {
      x <- mzx(env$mz[r]); h <- env$relative_intensity[r] * 70
      els <- c(els, sprintf('<line x1="%.1f" y1="%d" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2" stroke-opacity="0.85"/>', x, base_y, x, base_y - h, color))
    }
    els <- c(els, sprintf('<text x="44" y="%d" font-size="10" fill="%s" font-weight="600">%s: %.1f Da</text>', 12 + (i - 1) * 12, color, id, masses[[id]]))
  }
  ticks <- c(round(mz_min), round((mz_min + mz_max) / 2), round(mz_max))
  anchors <- c("start", "middle", "end")
  for (k in 1:3) {
    els <- c(els, sprintf('<text x="%.1f" y="98" font-size="9" fill="#666" text-anchor="%s">%d m/z</text>', mzx(ticks[k]), anchors[k], ticks[k]))
  }
  sprintf('<svg viewBox="0 0 640 130" class="pt-viz">%s</svg>', paste(els, collapse = ""))
}

#' Render the MS2 fragment-ladder tier view as a static SVG: one row per
#' proteoform (its own b-ion row, its own y-ion row, PTM markers), colored by
#' the tiers from compute_ladder_tiers().
#'
#' @param pf_list named list of proteoform objects
#' @param tiers output of compute_ladder_tiers(pf_list, ...)
#' @return character(1) SVG markup
render_ladder_svg <- function(pf_list, tiers) {
  ids <- names(pf_list)
  if (length(ids) == 0) return("<svg viewBox='0 0 640 40'><text x='10' y='20'>No proteoforms checked.</text></svg>")

  row_h <- 78
  px0 <- 90; px1 <- 620
  height <- length(ids) * row_h + 10
  els <- character(0)

  for (i in seq_along(ids)) {
    id <- ids[i]
    pf <- pf_list[[id]]
    tt <- tiers[[id]]
    len <- nchar(pf$sequence)
    top <- (i - 1) * row_h + 4
    color <- .IDENTITY_PALETTE[((i - 1) %% length(.IDENTITY_PALETTE)) + 1]
    xs <- function(pos) px0 + (pos / len) * (px1 - px0)

    els <- c(els, sprintf('<text x="4" y="%d" font-size="11" font-weight="600" fill="%s">%s</text>', top + 9, color, id))
    els <- c(els, sprintf('<text x="4" y="%d" font-size="9" fill="#666">%d aa</text>', top + 20, len))
    els <- c(els, sprintf('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#ddd" stroke-width="1"/>', px0, top + 24, px1, top + 24))

    for (ptm in pf$ptms) {
      x <- xs(ptm$site)
      pcol <- .PTM_COLOR_DEFAULT
      els <- c(els, sprintf('<line x1="%.1f" y1="24" x2="%.1f" y2="13" stroke="%s" stroke-width="1.5" transform="translate(0,%d)"/><circle cx="%.1f" cy="11" r="3" fill="%s" transform="translate(0,%d)"/><text x="%.1f" y="9" font-size="7" text-anchor="middle" fill="#555" transform="translate(0,%d)">%s</text>',
        x, x, pcol, top, x, pcol, top, x, top, ptm$name %||% "PTM"))
    }

    n_bonds <- length(tt$ladder$b_mass)
    for (p in seq_len(n_bonds)) {
      x <- xs(p)
      bc <- .TIER_COLOR[[tt$tier_b[p]]]
      yc <- .TIER_COLOR[[tt$tier_y[p]]]
      els <- c(els, sprintf('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke-width="1.2" stroke="%s" stroke-opacity="0.85"/>', x, top + 30, x, top + 42, bc))
      els <- c(els, sprintf('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke-width="1.2" stroke="%s" stroke-opacity="0.85"/>', x, top + 48, x, top + 60, yc))
    }
    els <- c(els, sprintf('<text x="4" y="%d" font-size="9" fill="#999">b</text>', top + 38))
    els <- c(els, sprintf('<text x="4" y="%d" font-size="9" fill="#999">y</text>', top + 56))
  }
  sprintf('<svg viewBox="0 0 640 %d" class="pt-viz">%s</svg>', height, paste(els, collapse = ""))
}
