# Fragmentation propensity model (design spec section 5): a literature-
# grounded, rule-based per-bond cleavage propensity score -- not a trained
# intensity-prediction model (none exists for intact-protein fragmentation).
#
# v1 scope: HCD/CID only (b/y ions). ETD/ECD (c/z ions, flatter
# residue-selectivity, structural-accessibility term via IUPred3/AlphaFold
# DB) and UVPD are deferred -- see ptracker_vibe.md.
#
# Four multiplicative effects apply to HCD/CID:
#   - residue-pair effects (Pro/Asp enhancement)
#   - positional (terminal-vs-internal) bias
#   - local basic-residue ("charge density") effect
#   - structural accessibility -- ETD/ECD only, so always a 1.0 no-op here
# combined multiplicatively into a per-bond propensity score.
#
# denatured-mode constants below are CALIBRATED against real matched b/y
# ions from two independent public top-down datasets (MassIVE MSV000094311
# and MSV000098558, TDPortal-searched, both HCD, ~1.83M + 0.89M confident
# per-bond observations), via a pooled logistic regression of `observed`
# on the same features this model uses, with `dataset` as a nuisance
# covariate -- every effect below replicated in sign and rough magnitude
# across both datasets independently before being folded in here. See
# scripts/validate_propensity_combined.R. native-mode constants are still
# the original literature-plausible placeholders (spec direction only, not
# fit against native-condition spectra -- no native dataset validated yet).
#
# One thing this calibration deliberately does NOT capture (see chat/
# DEVLOG for why): a real, consistently-replicated b/y detection asymmetry
# (y-ions ~1.8x more likely to be matched than b-ions at the same bond) --
# folding this in would mean returning per-ion-type scores instead of one
# score per bond, a bigger schema change deferred for now.
#
# PTM effects ARE partially captured, at the one site-resolved level that
# had enough data to calibrate: phosphorylation specifically (see
# PHOSPHO_SUPPRESSION comment). Every OTHER PTM type is still only used as
# a crude proteoform-level confound control during validation, not modeled
# at the site level -- phosphorylation was the only type with enough
# confident instances (>4,000 proteoforms across both datasets, ~44% of
# all PTM-site instances) to fit reliably; the rest have far too few
# observations per type to trust a per-type effect.
#
# One thing it DOES capture, deliberately, despite being a whole-protein
# property rather than a per-bond chemistry effect: overall proteoform
# length. Longer proteins get systematically worse per-bond MS2 coverage
# under a fixed instrument duty cycle (the same limited fragment-ion
# "attention" spreads over more candidate bonds) -- this is real, large,
# and directly relevant to what a user needs to know when comparing
# proteoforms of different length. It's a proteoform-level multiplier
# (identical for every bond in one proteoform), so it never changes which
# bond within a single proteoform ranks highest -- it shifts absolute tier
# thresholds and cross-proteoform comparisons, which is exactly the
# intended effect.

PROLINE_ENHANCEMENT <- list(denatured = 2.1, native = 5.0)
ASPARTATE_ENHANCEMENT <- list(denatured = 1.4, native = 5.0)
# indexed by min(terminal_distance, 6); distance-1 (touching the very
# terminus) is empirically the WEAKEST near-terminal position, not the
# strongest -- propensity peaks around distance 4-5, then relaxes to an
# intermediate plateau (still above the distance-1 floor) further in.
TERMINAL_PROXIMITY_WEIGHTS <- c(1.0, 1.2, 2.5, 4.1, 4.2, 1.5)
BASELINE_WEIGHT <- 1.0
# Stringency-tier cutpoints for the calibrated (GLM) propensity score --
# mirrored client-side in www/ptracker_viz.js's TIER_THRESHOLDS.glm (must
# stay in sync; see RF_ELEVATED_THRESHOLD etc. in
# R/fragmentation_propensity_rf.R for the length-free mode's own, numerically
# unrelated cutpoints). Fold-enrichment for matched fragment ions over
# baseline at these two thresholds: ~2.7x/~4.5x. Elevated lowered from an
# original 1.5 cutpoint: confirmed directly (scripts/validate_propensity_
# stringency_result.rds) that at 1.5+, any proteoform over ~300 residues has
# its max score across the ENTIRE ladder fall below the threshold -- the
# length term alone pushes everything under it, so realistic-length top-down
# targets cleared zero bonds at every tier. 1.2 recovers partial signal for
# short/medium proteoforms without changing the underlying model. A third
# "Very high" tier (was 10, fold ~5.9x) was dropped for the same reason: it
# almost never fired for realistic-length proteoforms either way.
GLM_ELEVATED_THRESHOLD <- 1.2
GLM_HIGH_THRESHOLD <- 4
# local basic-residue ("mobile proton") density: fraction of K/R within a
# +/-5-residue window centered on the bond. Counter to the classic
# tryptic-peptide mobile-proton intuition (more nearby K/R -> more
# cleavage), whole-protein HCD top-down data shows the OPPOSITE: bonds
# near locally K/R-dense stretches are matched less often. Real, replicated
# across both validation datasets independently -- not yet mechanistically
# explained, but consistent enough to calibrate.
CHARGE_DENSITY_WINDOW <- 5
CHARGE_DENSITY_WEIGHTS <- list(none = 1.0, low = 0.76, high = 0.33) # 0, 1, or 2+ K/R nearby
# proteoform length: log-linear in log(seq_len), coefficient fit adjusted
# for every other term above (odds ratio 0.38 per doubling of length,
# p < 1e-300, n=2.7M). Reference length 100 residues -> multiplier 1.0.
# Input is clamped to [LENGTH_CLAMP_MIN, LENGTH_CLAMP_MAX] before applying
# the power law: real proteins run from single digits to tens of thousands
# of residues, and an unclamped power law would produce runaway multipliers
# (>40x) for very short middle-down digestion peptides or >0 but vanishing
# multipliers for very long proteins -- clamping keeps the term inside the
# range the data can actually support (5th-95th percentile of validated
# proteoform lengths), extrapolating the calibrated trend smoothly within
# that band and holding flat beyond it, same spirit as the terminal-
# proximity table's clamp at distance 6+.
LENGTH_REFERENCE <- 100
LENGTH_EXPONENT <- -1.4
LENGTH_CLAMP_MIN <- 30
LENGTH_CLAMP_MAX <- 400
# phosphorylation suppresses local backbone fragmentation (established
# phosphoproteomics MS literature: neutral-loss-of-phosphate pathways and
# charge sequestration compete with standard b/y formation under CID/HCD).
# Confirmed against real data, adjusted for every other term in this file
# including length (odds ratio 0.23, p < 1e-300, n=2.7M): a bond within
# PHOSPHO_PROXIMITY_WINDOW residues of a phosphorylated site is
# substantially less likely to be matched. Detected via UNIMOD:21 in the
# proteoform's own ptms list (see proteoform_schema.R), not the
# TDPortal-specific PSI-MOD ids used only during validation-data extraction.
PHOSPHO_PROXIMITY_WINDOW <- 5
PHOSPHO_SUPPRESSION <- 0.23
PHOSPHO_UNIMOD_ID <- "UNIMOD:21"

#' Residue-pair propensity multiplier for one bond: enhanced when cleavage
#' is N-terminal to proline (residue_after == "P") and/or C-terminal to
#' aspartate (residue_before == "D"); the two enhancements multiply if both
#' apply at the same bond. Stronger under native than denatured conditions,
#' per spec.
#'
#' @param residue_before residue immediately N-terminal to the cleaved bond
#' @param residue_after residue immediately C-terminal to the cleaved bond
#' @param mode "denatured" or "native"
residue_pair_propensity <- function(residue_before, residue_after, mode = c("denatured", "native")) {
  mode <- match.arg(mode)
  score <- BASELINE_WEIGHT
  if (residue_after == "P") {
    score <- score * PROLINE_ENHANCEMENT[[mode]]
  }
  if (residue_before == "D") {
    score <- score * ASPARTATE_ENHANCEMENT[[mode]]
  }
  score
}

#' Positional (terminal-vs-internal) propensity multiplier, calibrated
#' against real matched-ion data (see TERMINAL_PROXIMITY_WEIGHTS comment):
#' peaks a few residues in from either terminus rather than right at it.
#'
#' @param cleavage_position 1-based interior cleavage position
#' @param sequence_length total sequence length (n)
positional_propensity <- function(cleavage_position, sequence_length) {
  terminal_distance <- min(cleavage_position, sequence_length - cleavage_position)
  idx <- min(terminal_distance, length(TERMINAL_PROXIMITY_WEIGHTS))
  TERMINAL_PROXIMITY_WEIGHTS[idx]
}

#' Raw local basic-residue (K/R) density in a +/-window around the bond,
#' as a fraction -- the underlying feature both charge_density_propensity()
#' (GLM multiplier) and the RF ranking mode (R/fragmentation_propensity_rf.R)
#' are calibrated against; shared here so the two scoring paths can never
#' silently drift apart in how this feature is computed.
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

#' Local basic-residue ("charge density") propensity multiplier: counts
#' K/R residues in a +/-CHARGE_DENSITY_WINDOW window around the bond.
#' HCD/CID only -- see CHARGE_DENSITY_WEIGHTS comment for why this is
#' currently HCD/CID-specific and unvalidated for ETD/ECD.
#'
#' @param residues full sequence as a character vector of single residues
#' @param cleavage_position 1-based interior cleavage position
charge_density_propensity <- function(residues, cleavage_position) {
  n <- length(residues)
  lo <- max(1, cleavage_position - CHARGE_DENSITY_WINDOW + 1)
  hi <- min(n, cleavage_position + CHARGE_DENSITY_WINDOW)
  n_basic <- sum(residues[lo:hi] %in% c("K", "R"))
  if (n_basic == 0) {
    CHARGE_DENSITY_WEIGHTS$none
  } else if (n_basic == 1) {
    CHARGE_DENSITY_WEIGHTS$low
  } else {
    CHARGE_DENSITY_WEIGHTS$high
  }
}

#' Proteoform-length propensity multiplier: identical for every bond in a
#' given proteoform (see LENGTH_REFERENCE/LENGTH_EXPONENT comment for why
#' this is calibrated and clamped). Shorter proteoforms get a boost,
#' longer ones a penalty, relative to a 100-residue reference.
#'
#' @param sequence_length total sequence length (n)
length_propensity <- function(sequence_length) {
  clamped <- min(max(sequence_length, LENGTH_CLAMP_MIN), LENGTH_CLAMP_MAX)
  (clamped / LENGTH_REFERENCE) ^ LENGTH_EXPONENT
}

#' Raw yes/no: is any phosphorylated residue (UNIMOD:21) within +/-window
#' of the bond? The underlying feature both phospho_proximity_propensity()
#' (GLM multiplier) and the RF ranking mode (R/fragmentation_propensity_rf.R)
#' are calibrated against; shared here for the same reason as
#' local_basic_density(). Only numeric (internal) PTM sites are considered
#' -- "N-term"/"C-term" ptm() entries have no residue position to measure
#' proximity from.
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

#' Phosphorylation-proximity propensity multiplier: suppressed when any
#' phosphorylated residue (UNIMOD:21) sits within
#' PHOSPHO_PROXIMITY_WINDOW residues of the bond. Only numeric (internal)
#' PTM sites are considered -- "N-term"/"C-term" ptm() entries have no
#' residue position to measure proximity from.
#'
#' @param ptms list of ptm() objects (proteoform$ptms)
#' @param cleavage_position 1-based interior cleavage position
phospho_proximity_propensity <- function(ptms, cleavage_position) {
  if (is_near_phosphosite(ptms, cleavage_position)) PHOSPHO_SUPPRESSION else BASELINE_WEIGHT
}

#' Structural-accessibility propensity multiplier. Only meaningful for
#' ETD/ECD (gas-phase charge density / higher-order structure effects);
#' always a 1.0 no-op for HCD/CID. Deferred, not yet implemented for
#' ETD/ECD -- see ptracker_vibe.md.
#'
#' @param method dissociation method
structural_accessibility_propensity <- function(method) {
  BASELINE_WEIGHT
}

#' Fragmentation propensity score at every backbone cleavage position of a
#' proteoform, for HCD/CID (b/y ions). Combines residue-pair, positional,
#' and structural-accessibility (no-op here) effects multiplicatively.
#'
#' @param proteoform a proteoform object (see proteoform_schema.R)
#' @param mode "denatured" or "native"
#' @param method dissociation method; only "HCD"/"CID" implemented (spec
#'   groups both under "collision-based fragmentation")
#' @return data.frame(cleavage_position, residue_before, residue_after,
#'   residue_pair_score, positional_score, charge_density_score,
#'   length_score, phospho_score, accessibility_score, propensity_score)
fragmentation_propensity <- function(proteoform, mode = c("denatured", "native"), method = "HCD") {
  mode <- match.arg(mode)
  if (!inherits(proteoform, "proteoform")) {
    stop("fragmentation_propensity() requires a proteoform object")
  }
  if (!toupper(method) %in% c("HCD", "CID")) {
    stop("only HCD/CID are implemented for v1 -- ETD/ECD/UVPD are deferred (see ptracker_vibe.md)")
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

  residue_pair_score <- mapply(
    residue_pair_propensity, residue_before, residue_after,
    MoreArgs = list(mode = mode)
  )
  positional_score <- vapply(cleavage_positions, positional_propensity, numeric(1), sequence_length = n)
  charge_density_score <- vapply(cleavage_positions, charge_density_propensity, numeric(1), residues = residues)
  length_score <- rep(length_propensity(n), length(cleavage_positions))
  phospho_score <- vapply(cleavage_positions, phospho_proximity_propensity, numeric(1), ptms = proteoform$ptms)
  accessibility_score <- rep(structural_accessibility_propensity(method), length(cleavage_positions))

  data.frame(
    cleavage_position = cleavage_positions,
    residue_before = residue_before,
    residue_after = residue_after,
    residue_pair_score = unname(residue_pair_score),
    positional_score = positional_score,
    charge_density_score = charge_density_score,
    length_score = length_score,
    phospho_score = phospho_score,
    accessibility_score = accessibility_score,
    propensity_score = unname(residue_pair_score) * positional_score * charge_density_score *
      length_score * phospho_score * accessibility_score
  )
}
