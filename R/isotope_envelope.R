# Isotope envelope width, via the Averagine model.
#
# pyteomics.mass.isotopologues() enumerates isotopic states combinatorially,
# which is fine for small molecules but computationally intractable at
# intact-protein scale (empirically: ~7s for a ~1 kDa test composition;
# unusable well before 16 kDa). Real deconvolution tools don't brute-force
# this for large molecules either -- they use a closed-form/polynomial
# approach. This implements the standard closed-form approximation: each
# element's isotope count is treated as an independent
# Binomial(n_atoms, natural abundance) of "heavy" substitutions, summed via
# the Central Limit Theorem into a Gaussian approximation of the whole
# envelope (mean shift from monoisotopic mass, and a width/FWHM).
#
# Uses the Averagine model (classic Senko/Rockwood/Marshall idealized
# "average amino acid residue" composition, scaled to the protein's own
# mass) rather than the real sequence's exact elemental composition, per
# the design spec's guidance that Averagine is the field-standard
# approximation above a few kDa -- this keeps results comparable to other
# TDP tools that make the same assumption.
#
# Isotope masses/natural abundances below are the dominant non-monoisotopic
# isotope per element, from pyteomics.mass.nist_mass (NIST atomic weights
# and isotopic compositions); minor secondary isotopes (e.g. 17-O) are
# ignored as a simplification.

AVERAGINE_UNIT_MASS <- 111.1254
AVERAGINE_UNIT_COMPOSITION <- c(C = 4.9384, H = 7.7583, N = 1.3577, O = 1.4773, S = 0.0417)

ELEMENT_ISOTOPE_SHIFTS <- data.frame(
  element = c("C", "H", "N", "O", "S"),
  mass_diff = c(1.0033548378, 1.0062767457, 0.9970348934, 2.0042463804, 1.9957959),
  abundance = c(0.0107, 0.000115, 0.00364, 0.00205, 0.0425)
)

#' Averagine elemental composition scaled to a target mass.
#'
#' @param mass neutral mass, Da
#' @return named numeric vector of atom counts (C, H, N, O, S)
averagine_composition <- function(mass) {
  n_units <- mass / AVERAGINE_UNIT_MASS
  round(AVERAGINE_UNIT_COMPOSITION * n_units)
}

#' Averagine-based isotope envelope statistics: how far the most abundant
#' isotope peak sits above the monoisotopic mass, and how wide the overall
#' envelope is.
#'
#' @param mass neutral monoisotopic mass, Da
#' @return list(mean_shift, sigma, fwhm), all in Da
averagine_isotope_envelope <- function(mass) {
  comp <- averagine_composition(mass)
  shifts <- ELEMENT_ISOTOPE_SHIFTS
  n <- comp[shifts$element]
  mean_shift <- sum(n * shifts$abundance * shifts$mass_diff)
  variance <- sum(n * shifts$abundance * (1 - shifts$abundance) * shifts$mass_diff^2)
  sigma <- sqrt(variance)
  list(mean_shift = mean_shift, sigma = sigma, fwhm = 2.3548 * sigma)
}
