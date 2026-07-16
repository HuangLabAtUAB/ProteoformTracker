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
#' @return named list, id -> list(ladder, tier_b, tier_y) where tier_b/tier_y
#'   are character vectors ("common"/"partial"/"unique"/"neutral"), one per
#'   cleavage position
compute_ladder_tiers <- function(pf_list, r_ref = 120000, mz_ref = 200, safety_margin = DEFAULT_SAFETY_MARGIN, mode = "denatured") {
  ids <- names(pf_list)
  result <- list()
  for (id in ids) {
    others <- pf_list[setdiff(ids, id)]
    ladder <- generate_fragment_ladder(pf_list[[id]])
    prop <- fragmentation_propensity(pf_list[[id]], mode = mode, method = "HCD")$propensity_score
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

#' Same idea as compute_ladder_tiers() but for the single-target confounder
#' view: tiers a target's own ladder against a set of confounder proteoforms
#' rather than other checked proteoforms.
#'
#' @param target_pf proteoform object (the confounder-search target)
#' @param confounder_pfs named list of proteoform objects (real confounders,
#'   from build_confounder_proteoforms())
#' @return list(ladder, propensity, tier_b, tier_y) for the target
compute_confounder_tiers <- function(target_pf, confounder_pfs, r_ref = 120000, mz_ref = 200,
                                      safety_margin = DEFAULT_SAFETY_MARGIN, mode = "denatured") {
  ladder <- generate_fragment_ladder(target_pf)
  prop <- fragmentation_propensity(target_pf, mode = mode, method = "HCD")$propensity_score
  n_bonds <- length(ladder$b_mass)
  if (length(confounder_pfs) == 0) {
    return(list(ladder = ladder, propensity = prop, tier_b = rep("neutral", n_bonds), tier_y = rep("neutral", n_bonds)))
  }
  cc <- fragment_mass_collision_check(target_pf, confounder_pfs, r_ref = r_ref, mz_ref = mz_ref, safety_margin = safety_margin)
  n_others <- length(confounder_pfs)
  match_b <- Reduce(`+`, lapply(cc$per_candidate, function(d) as.integer(d$matched[d$ion_type == "b"])))
  match_y <- Reduce(`+`, lapply(cc$per_candidate, function(d) as.integer(d$matched[d$ion_type == "y"])))
  tier_of <- function(mc) ifelse(mc == n_others, "common", ifelse(mc == 0, "unique", "partial"))
  list(ladder = ladder, propensity = prop, tier_b = tier_of(match_b), tier_y = tier_of(match_y))
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
