# Fragmentation propensity model (design spec section 5): a literature-
# grounded, rule-based per-bond cleavage propensity score -- not a trained
# intensity-prediction model (none exists for intact-protein fragmentation).
#
# v1 scope: HCD/CID only (b/y ions). ETD/ECD (c/z ions, flatter
# residue-selectivity, structural-accessibility term via IUPred3/AlphaFold
# DB) and UVPD are deferred -- see ptracker_vibe.md.
#
# Three of the spec's four multiplicative effects apply to HCD/CID:
#   - residue-pair effects (Pro/Asp enhancement)
#   - positional (terminal-vs-internal) bias
#   - structural accessibility -- ETD/ECD only, so always a 1.0 no-op here
# combined multiplicatively into a per-bond propensity score.
#
# All specific multiplier values below are illustrative, not calibrated
# against real spectra -- the spec states the *direction* of each effect
# (stronger under native, terminal ions dominate for ~3 events) but not
# universal magnitudes. Phase 4 validates against real example pairs.

PROLINE_ENHANCEMENT <- list(denatured = 3.0, native = 5.0)
ASPARTATE_ENHANCEMENT <- list(denatured = 3.0, native = 5.0)
TERMINAL_PROXIMITY_EVENTS <- 3
TERMINAL_PROXIMITY_WEIGHTS <- c(4.0, 3.5, 3.0) # for distance-from-terminus 1, 2, 3
BASELINE_WEIGHT <- 1.0

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

#' Positional (terminal-vs-internal) propensity multiplier: terminal
#' fragment ions are statistically favored to carry the majority of
#' fragment-ion current for the first ~3 backbone cleavage events from
#' either terminus, regardless of protein size; flat baseline beyond that.
#'
#' @param cleavage_position 1-based interior cleavage position
#' @param sequence_length total sequence length (n)
positional_propensity <- function(cleavage_position, sequence_length) {
  terminal_distance <- min(cleavage_position, sequence_length - cleavage_position)
  if (terminal_distance <= TERMINAL_PROXIMITY_EVENTS) {
    TERMINAL_PROXIMITY_WEIGHTS[terminal_distance]
  } else {
    BASELINE_WEIGHT
  }
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
#'   residue_pair_score, positional_score, accessibility_score, propensity_score)
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
  accessibility_score <- rep(structural_accessibility_propensity(method), length(cleavage_positions))

  data.frame(
    cleavage_position = cleavage_positions,
    residue_before = residue_before,
    residue_after = residue_after,
    residue_pair_score = unname(residue_pair_score),
    positional_score = positional_score,
    accessibility_score = accessibility_score,
    propensity_score = unname(residue_pair_score) * positional_score * accessibility_score
  )
}
