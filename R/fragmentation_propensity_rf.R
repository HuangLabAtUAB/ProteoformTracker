# Optional secondary scoring mode: a length-free RF ranking model, sitting
# alongside (never replacing) the primary calibrated GLM formula in
# R/fragmentation_propensity.R. See scripts/build_propensity_rf_model.R for
# the full rationale, training data, and held-out validation numbers.
#
# In short: the GLM's length term is a real, validated effect, but it's a
# per-proteoform constant multiplier that can suppress every bond in one
# long proteoform below the tier thresholds, even though the RELATIVE
# ranking of bonds within that proteoform is what a user actually needs.
# This mode deliberately excludes length so long proteoforms stay usable.
# It requires data/fragmentation_propensity_rf.rds (built by the script
# above) and the `ranger` package; both loaded once at startup in global.R
# into `propensity_rf_model` (NULL if unavailable -- callers should check).
#
# Output is a genuine probability (0-1), not a multiplicative score like
# the GLM's -- do not compare RF and GLM values directly, and do not reuse
# GLM tier thresholds here. See RF_TIER_* constants below, calibrated
# separately against the same validation data (scripts/
# validate_propensity_rf_deploy_candidate.rds sweep).

RF_ELEVATED_THRESHOLD <- 0.08
RF_HIGH_THRESHOLD <- 0.15
RF_VERY_HIGH_THRESHOLD <- 0.30
# isotope-panel computation gate (see R/isotope_envelope.R callers):
# reuses the Elevated cutoff, since RF has no natural "baseline = 1.0"
# concept the way the GLM's multiplicative formula does.
RF_ISOTOPE_GATE <- RF_ELEVATED_THRESHOLD

#' Length-free RF-ranked fragmentation propensity at every backbone
#' cleavage position of a proteoform. Same bond-level granularity and
#' column-naming spirit as fragmentation_propensity(), but `propensity_score`
#' here is a model-predicted probability (0-1), not a multiplicative score
#' -- see file header before wiring this into anything that also consumes
#' fragmentation_propensity()'s output on the same scale.
#'
#' @param proteoform a proteoform object (see proteoform_schema.R)
#' @param rf_model a fitted ranger model (typically the app's global
#'   `propensity_rf_model`); NULL is treated as "mode unavailable" and
#'   errors clearly rather than silently falling back
#' @return data.frame(cleavage_position, residue_before, residue_after,
#'   terminal_distance, basic_density, near_phospho, propensity_score)
fragmentation_propensity_rf <- function(proteoform, rf_model = propensity_rf_model) {
  if (!inherits(proteoform, "proteoform")) {
    stop("fragmentation_propensity_rf() requires a proteoform object")
  }
  if (is.null(rf_model)) {
    stop("RF ranking mode unavailable -- run scripts/build_propensity_rf_model.R ",
         "and ensure the 'ranger' package is installed")
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
  basic_density <- vapply(cleavage_positions, local_basic_density, numeric(1), residues = residues)
  near_phospho <- vapply(cleavage_positions, is_near_phosphosite, logical(1), ptms = proteoform$ptms)

  newdata <- data.frame(
    has_pro = has_pro, has_asp = has_asp, terminal_distance = terminal_distance,
    basic_density = basic_density, near_phospho = factor(near_phospho, levels = c(FALSE, TRUE))
  )
  propensity_score <- predict(rf_model, data = newdata)$predictions

  data.frame(
    cleavage_position = cleavage_positions,
    residue_before = residue_before,
    residue_after = residue_after,
    terminal_distance = terminal_distance,
    basic_density = basic_density,
    near_phospho = near_phospho,
    propensity_score = propensity_score
  )
}
