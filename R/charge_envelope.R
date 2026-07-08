# Charge-state envelope and signal-to-noise modeling.
# No trained forward-simulation model is wired in yet (UniDec is the
# upgrade path, see Phase 4 in the design spec) -- these are literature-
# grounded heuristics, clearly approximate, sufficient to drive Phase 2
# scoring and visualization.

#' Count sequence-derivable protonatable basic sites (Arg, Lys, plus the
#' free N-terminal amine). Used as the charge-state ceiling for small
#' denatured proteoforms, where the spec calls for a ceiling that "tracks
#' basic-residue count."
#'
#' @param sequence amino acid sequence
#' @param include_n_terminal_amine count the free N-terminal amine as an
#'   additional protonatable site (default TRUE)
basic_residue_count <- function(sequence, include_n_terminal_amine = TRUE) {
  chars <- strsplit(toupper(sequence), "")[[1]]
  count <- sum(chars %in% c("R", "K"))
  if (include_n_terminal_amine) count <- count + 1
  max(count, 1)
}

#' Rayleigh-limit charge-state ceiling: z ~= 0.0778 * sqrt(mass), the
#' de la Mora (2000) scaling relating maximum native ESI charge to a
#' droplet-like surface-area/mass relationship. Applies above ~40 kDa
#' (denatured) or generally under native-like conditions, per spec.
#'
#' @param mass neutral mass, Da
rayleigh_charge_ceiling <- function(mass) {
  max(1, floor(0.0778 * sqrt(mass)))
}

#' Charge-state ceiling for a proteoform, following the spec's rule:
#' denaturing and <=~40 kDa -> basic-residue-count-derived; above that, or
#' native-like generally -> Rayleigh-limit relationship.
#'
#' @param sequence amino acid sequence (used only for the small-denatured case)
#' @param mass neutral mass, Da
#' @param mode "denatured" or "native"
charge_state_ceiling <- function(sequence, mass, mode = c("denatured", "native")) {
  mode <- match.arg(mode)
  if (mode == "denatured" && mass <= 40000) {
    basic_residue_count(sequence)
  } else {
    rayleigh_charge_ceiling(mass)
  }
}

#' Predict a plausible charge-state envelope for a proteoform: which charge
#' states are populated, and their relative (not absolute) intensities.
#'
#' Heuristic Gaussian-shaped envelope over [z_min, ceiling], not a forward
#' simulation. Denatured envelopes are modeled broader and peaked lower
#' relative to their ceiling (more uniform accessibility when unfolded);
#' native envelopes narrower and peaked closer to their (lower) ceiling.
#' Flag as approximate in the UI -- Phase 4 validates against UniDec.
#'
#' @param sequence amino acid sequence
#' @param mass neutral mass, Da (mono or average, caller's choice -- must be
#'   consistent with how R(m/z)/FWHM are computed downstream)
#' @param mode "denatured" or "native"
#' @return data.frame(z, mz, relative_intensity), relative_intensity peak-normalized to 1
predict_charge_envelope <- function(sequence, mass, mode = c("denatured", "native")) {
  mode <- match.arg(mode)
  ceiling <- charge_state_ceiling(sequence, mass, mode)

  if (mode == "denatured") {
    z_min <- max(1, round(0.35 * ceiling))
    z_peak <- 0.75 * ceiling
    spread <- max(0.25 * ceiling, 1)
  } else {
    z_min <- max(1, round(0.65 * ceiling))
    z_peak <- 0.9 * ceiling
    spread <- max(0.12 * ceiling, 1)
  }
  z_min <- min(z_min, ceiling)

  z <- z_min:ceiling
  intensity <- exp(-0.5 * ((z - z_peak) / spread)^2)
  intensity <- intensity / max(intensity)

  data.frame(z = z, mz = mz_for_charge(mass, z), relative_intensity = intensity)
}

#' Relative, normalized S/N penalty as mass increases, against a
#' user-supplied baseline mass (their own reference protein/instrument
#' setup). NOT an absolute intensity or detectability prediction -- absolute
#' signal depends on sample loading/instrument sensitivity this tool can't
#' know. Denaturing conditions are modeled with faster S/N decay than
#' native (signal divides across more charge states and a wider isotope
#' envelope as mass increases; the effect is stronger when denatured).
#'
#' @param mass neutral mass of the proteoform being scored, Da
#' @param baseline_mass mass of the user's reference protein (S/N == 1 there)
#' @param mode "denatured" or "native"
#' @return relative penalty; 1 at baseline_mass, below 1 as mass grows past it
sn_relative_penalty <- function(mass, baseline_mass, mode = c("denatured", "native")) {
  mode <- match.arg(mode)
  decay_exponent <- if (mode == "denatured") 1.0 else 0.5
  (baseline_mass / mass)^decay_exponent
}
