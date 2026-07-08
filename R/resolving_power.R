# Orbitrap resolving-power model.
# R(m/z) = R_ref * sqrt(m/z_ref / m/z); resolving power is quoted at one
# reference m/z (conventionally 200) and degrades as m/z increases.
# Instrument-family specific: this scaling law is Orbitrap-only. Do not
# reuse for FT-ICR or other platforms without a different relationship.

PROTON_MASS <- 1.007276

#' Orbitrap resolving power at a given m/z.
#'
#' @param mz observed m/z
#' @param r_ref resolving power quoted at mz_ref (instrument setting)
#' @param mz_ref reference m/z the resolving power is quoted at (conventionally 200)
resolving_power <- function(mz, r_ref = 120000, mz_ref = 200) {
  r_ref * sqrt(mz_ref / mz)
}

#' m/z of a proteoform at a given charge state (singly protonated additions).
#'
#' @param mass neutral monoisotopic (or average) mass, Da
#' @param z charge state (positive integer)
mz_for_charge <- function(mass, z) {
  (mass + z * PROTON_MASS) / z
}

#' Mass-domain FWHM peak width at a given m/z: DeltaM_FWHM (Da) ~= M / R(m/z).
#'
#' @param mass neutral mass, Da
#' @param mz m/z the peak is observed at (i.e. the charge state's m/z)
#' @param r_ref,mz_ref passed to resolving_power()
fwhm_mass <- function(mass, mz, r_ref = 120000, mz_ref = 200) {
  mass / resolving_power(mz, r_ref, mz_ref)
}

#' m/z-domain FWHM peak width at a given m/z: mz / R(m/z). Used for
#' envelope-crowding checks, which compare peak positions directly in the
#' m/z domain (where charge-state interference between different masses
#' actually shows up), rather than the mass domain.
#'
#' @param mz m/z the peak is observed at
#' @param r_ref,mz_ref passed to resolving_power()
fwhm_mz <- function(mz, r_ref = 120000, mz_ref = 200) {
  mz / resolving_power(mz, r_ref, mz_ref)
}

#' FWHM (mass domain) at every charge state in a proteoform's charge-state
#' envelope, and which charge state gives the best (smallest) FWHM.
#'
#' R(m/z) increases monotonically as m/z decreases, so in practice the
#' highest available charge state always wins -- but every charge state is
#' evaluated explicitly (not just assumed) so this stays correct if the
#' envelope is later filtered by an intensity/S-N threshold.
#'
#' @param mass neutral mass, Da
#' @param charge_states integer vector of populated charge states
#' @param r_ref,mz_ref passed to resolving_power()
#' @return data.frame(z, mz, fwhm_mass), sorted by charge ascending, plus
#'   attr(., "best") giving the row index of the best (smallest FWHM) entry
fwhm_by_charge_state <- function(mass, charge_states, r_ref = 120000, mz_ref = 200) {
  charge_states <- sort(unique(as.integer(charge_states)))
  mz <- mz_for_charge(mass, charge_states)
  fwhm <- fwhm_mass(mass, mz, r_ref, mz_ref)
  result <- data.frame(z = charge_states, mz = mz, fwhm_mass = fwhm)
  attr(result, "best") <- which.min(fwhm)
  result
}
