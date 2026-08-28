# Fragmentation propensity model (design spec section 5): predicts, per
# backbone cleavage position, the probability a b/y fragment ion would
# actually be observed in real HCD/CID top-down MS2 data.
#
# v1 scope: HCD/CID only (b/y ions). ETD/ECD (c/z ions, flatter
# residue-selectivity, structural-accessibility term via IUPred3/AlphaFold
# DB) and UVPD are deferred.
#
# This is a single joint logistic regression (data/fragmentation_
# propensity_glm.rds, built by scripts/build_propensity_glm_model.R),
# loaded once at startup into `propensity_glm_model` (global.R) and scored
# via predict() here -- the SAME load-once-and-predict() pattern
# fragmentation_propensity_rf() already uses for RF mode. Output is a
# genuine probability (0-1), not the old hand-multiplied score.
#
# This replaces an earlier version of this file that assembled the score
# from several SEPARATELY-calibrated pieces multiplied together (residue-
# pair/terminal-distance/charge-density from one joint fit with no length
# term at all; length from a separate per-dataset GAM exploration;
# phospho-proximity from yet another separate pass) under the assumption
# that independently-calibrated effects combine multiplicatively without
# interaction. That assumption was never validated and turned out to be
# wrong: cross-study holdout AUROC for the assembled version was ~0.64,
# vs. ~0.78-0.82 for one joint fit on the exact same feature categories
# (see chat history / scripts/build_propensity_glm_model.R's header for the
# full investigation). The features themselves are unchanged from before --
# same residue-pair indicators, same 6-level terminal-distance cap, same
# 3-level charge-density buckets, same phospho-proximity window, same
# length-clamping philosophy -- only HOW they're combined changed.
#
# Fit on real matched b/y ions from two independent public top-down
# datasets (MassIVE MSV000094311 and MSV000098558, TDPortal-searched, both
# HCD, ~2.7M pooled confident per-bond observations). `dataset` and
# `ion_type` were included as nuisance covariates during fitting (the two
# datasets have different baseline hit rates; there's a real, consistently
# -replicated b/y detection asymmetry -- y-ions ~1.8x more likely matched
# than b-ions at the same bond) so the terms of interest are correctly
# adjusted for both, but neither is a real property of a bond in isolation
# -- predict_glm_propensity() below marginalizes both out by averaging
# predictions across all four (ion_type x dataset) combinations, preserving
# the existing one-score-per-bond contract (shared between the resulting b
# and y ions) rather than returning per-ion-type scores.
#
# No native-mode-specific calibration exists (no native-condition dataset
# has been validated) -- native mode currently reuses this same
# denatured-calibrated model as the best available stand-in, same
# unresolved gap the previous placeholder constants had, just now on a
# consistent 0-1 scale instead of a separate, never-validated multiplicative
# guess.
#
# PTM effects ARE partially captured, at the one site-resolved level that
# had enough data to calibrate: phosphorylation specifically. Every OTHER
# PTM type is still only used as a crude proteoform-level confound control
# during validation, not modeled at the site level -- phosphorylation was
# the only type with enough confident instances to fit reliably.
#
# Proteoform length IS captured, deliberately, despite being a whole-protein
# property rather than a per-bond chemistry effect: longer proteins get
# systematically worse per-bond MS2 coverage under a fixed instrument duty
# cycle. This is real, large, and confirmed directly (scripts/
# build_propensity_glm_model.R's validation) to still leave most bonds in
# proteoforms over ~300-400 residues below the Elevated tier even under the
# corrected joint fit -- that's an honest reflection of the real effect,
# not a modeling bug, and is exactly why the length-free RF ranking mode
# (R/fragmentation_propensity_rf.R) exists as a complementary option for
# long proteoforms.

# Stringency-tier cutpoints for the calibrated (GLM) propensity score --
# mirrored client-side in www/ptracker_viz.js's TIER_THRESHOLDS.glm (must
# stay in sync; see RF_ELEVATED_THRESHOLD etc. in
# R/fragmentation_propensity_rf.R for the length-free mode's own, numerically
# unrelated cutpoints -- both are now on a comparable 0-1 probability scale,
# but the two models are still fit independently and must never share
# thresholds). Selected via the same fold-enrichment sweep methodology as
# RF's thresholds: ~3.1x/~4.7x fold-enrichment for matched fragment ions
# over the ~6.3% pooled baseline match rate.
GLM_ELEVATED_THRESHOLD <- 0.10
GLM_HIGH_THRESHOLD <- 0.20

# local basic-residue ("mobile proton") density: fraction of K/R within a
# +/-5-residue window centered on the bond. Counter to the classic
# tryptic-peptide mobile-proton intuition (more nearby K/R -> more
# cleavage), whole-protein HCD top-down data shows the OPPOSITE: bonds
# near locally K/R-dense stretches are matched less often. Real, replicated
# across both validation datasets independently -- not yet mechanistically
# explained, but consistent enough to calibrate.
CHARGE_DENSITY_WINDOW <- 5

# phosphorylation suppresses local backbone fragmentation (established
# phosphoproteomics MS literature: neutral-loss-of-phosphate pathways and
# charge sequestration compete with standard b/y formation under CID/HCD).
# Detected via UNIMOD:21 in the proteoform's own ptms list (see
# proteoform_schema.R), not the TDPortal-specific PSI-MOD ids used only
# during validation-data extraction.
PHOSPHO_PROXIMITY_WINDOW <- 5
PHOSPHO_UNIMOD_ID <- "UNIMOD:21"

#' Raw local basic-residue (K/R) density in a +/-window around the bond,
#' as a fraction -- the underlying feature both this file's GLM model and
#' the RF ranking mode (R/fragmentation_propensity_rf.R) are calibrated
#' against; shared here so the two scoring paths can never silently drift
#' apart in how this feature is computed.
#'
#' @param residues full sequence as a character vector of single residues
#' @param cleavage_position 1-based interior cleavage position
#' @param window residues on each side to include (default CHARGE_DENSITY_WINDOW)
local_basic_density <- function(residues, cleavage_position, window = CHARGE_DENSITY_WINDOW) {
  n <- length(residues)
  lo <- max(1, cleavage_position - window + 1)
  hi <- min(n, cleavage_position + window)
  mean(residues[lo:hi] %in% c("K", "R"))
}

#' Raw yes/no: is any phosphorylated residue (UNIMOD:21) within +/-window
#' of the bond? The underlying feature both this file's GLM model and the
#' RF ranking mode (R/fragmentation_propensity_rf.R) are calibrated
#' against; shared here for the same reason as local_basic_density(). Only
#' numeric (internal) PTM sites are considered -- "N-term"/"C-term" ptm()
#' entries have no residue position to measure proximity from.
#'
#' @param ptms list of ptm() objects (proteoform$ptms)
#' @param cleavage_position 1-based interior cleavage position
#' @param window residues on each side to include (default PHOSPHO_PROXIMITY_WINDOW)
is_near_phosphosite <- function(ptms, cleavage_position, window = PHOSPHO_PROXIMITY_WINDOW) {
  if (length(ptms) == 0) return(FALSE)
  is_phospho <- vapply(ptms, function(p) identical(p$unimod_id, PHOSPHO_UNIMOD_ID) || identical(p$name, "Phospho"), logical(1))
  sites <- vapply(ptms[is_phospho], function(p) if (is.numeric(p$site)) p$site else NA_integer_, integer(1))
  sites <- sites[!is.na(sites)]
  if (length(sites) == 0) return(FALSE)
  any(abs(sites - cleavage_position) <= window)
}

#' Fragmentation propensity score at every backbone cleavage position of a
#' proteoform, for HCD/CID (b/y ions): a genuine predicted probability
#' (0-1) that a b/y ion at that bond would be observed in real top-down
#' MS2 data, from the joint GLM model (data/fragmentation_propensity_glm.rds).
#'
#' @param proteoform a proteoform object (see proteoform_schema.R)
#' @param mode "denatured" or "native" -- native mode currently reuses the
#'   same denatured-calibrated model (no native-condition dataset has been
#'   validated yet), accepted for interface compatibility/future use
#' @param method dissociation method; only "HCD"/"CID" implemented (spec
#'   groups both under "collision-based fragmentation")
#' @param glm_model a fitted glm object (typically the app's global
#'   `propensity_glm_model`); NULL is treated as "mode unavailable" and
#'   errors clearly rather than silently falling back
#' @return data.frame(cleavage_position, residue_before, residue_after,
#'   propensity_score) -- propensity_score is a probability, not a
#'   multiplicative score; do not compare directly against RF mode's
#'   thresholds (R/fragmentation_propensity_rf.R) even though both are now
#'   0-1 scales, since the two models are fit independently
fragmentation_propensity <- function(proteoform, mode = c("denatured", "native"), method = "HCD",
                                      glm_model = propensity_glm_model) {
  mode <- match.arg(mode)
  if (!inherits(proteoform, "proteoform")) {
    stop("fragmentation_propensity() requires a proteoform object")
  }
  if (!toupper(method) %in% c("HCD", "CID")) {
    stop("only HCD/CID are implemented for v1 -- ETD/ECD/UVPD are deferred")
  }
  if (is.null(glm_model)) {
    stop("Calibrated (GLM) scoring mode unavailable -- run scripts/build_propensity_glm_model.R")
  }

  sequence <- proteoform$sequence
  n <- nchar(sequence)
  if (n < 2) {
    stop("sequence must have at least 2 residues to score fragmentation propensity")
  }

  cleavage_positions <- seq_len(n - 1)
  residues <- strsplit(sequence, "")[[1]]
  residue_before <- residues[cleavage_positions]
  residue_after <- residues[cleavage_positions + 1]

  has_pro <- residue_after == "P"
  has_asp <- residue_before == "D"
  terminal_distance <- pmin(cleavage_positions, n - cleavage_positions)
  td_bucket <- factor(pmin(terminal_distance, 6), levels = 1:6)
  basic_density <- vapply(cleavage_positions, local_basic_density, numeric(1), residues = residues)
  basic_bucket <- cut(basic_density, breaks = c(-Inf, 0, 0.15, Inf), labels = c("none", "low", "high"))
  near_phospho <- vapply(cleavage_positions, is_near_phosphosite, logical(1), ptms = proteoform$ptms)

  length_clamp_min <- attr(glm_model, "length_clamp_min") %||% 20
  length_clamp_max <- attr(glm_model, "length_clamp_max") %||% 320
  clamped_len <- min(max(n, length_clamp_min), length_clamp_max)

  newdata_base <- data.frame(
    has_pro = has_pro, has_asp = has_asp, td_bucket = td_bucket,
    basic_bucket = basic_bucket, near_phospho = factor(near_phospho, levels = c(FALSE, TRUE)),
    clamped_len = clamped_len
  )
  # Marginalize the two nuisance covariates the model was fit with
  # (ion_type, dataset -- see file header) by averaging predictions across
  # all four combinations, rather than picking one arbitrarily.
  nuisance_combos <- expand.grid(
    ion_type = factor(c("B", "Y"), levels = c("B", "Y")),
    dataset = factor(c("MSV000094311", "MSV000098558"), levels = c("MSV000094311", "MSV000098558"))
  )
  preds <- vapply(seq_len(nrow(nuisance_combos)), function(i) {
    nd <- newdata_base
    nd$ion_type <- nuisance_combos$ion_type[i]
    nd$dataset <- nuisance_combos$dataset[i]
    predict(glm_model, newdata = nd, type = "response")
  }, numeric(length(cleavage_positions)))
  propensity_score <- if (is.matrix(preds)) rowMeans(preds) else mean(preds)

  data.frame(
    cleavage_position = cleavage_positions,
    residue_before = residue_before,
    residue_after = residue_after,
    propensity_score = propensity_score
  )
}
