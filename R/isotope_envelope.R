# Isotope envelope: summary statistics AND the full discrete peak pattern,
# both derived from the real, exact-composition isotope distribution
# (sequence_isotope_pattern() / proteoform_isotope_pattern(), in
# R/mass_calculation.R) rather than the Averagine model this file used to
# implement.
#
# Averagine (Senko/Rockwood/Marshall) approximates a protein's elemental
# composition from its MASS ALONE, by scaling an idealized "average amino
# acid residue" composition to the target mass. It exists for the case
# where the real sequence is unknown -- e.g. deconvoluting a raw intact-mass
# measurement with no sequence attached yet. Every proteoform in this app
# already has a known sequence, so approximating its composition from mass
# alone throws away information that's already sitting right there: summing
# the real residues' own elemental formulas (what
# sequence_isotope_pattern() does, via pyteomics.mass.Composition + IsoSpecPy
# -- see that function's doc comment) gives the protein's EXACT composition,
# at the same computational cost (milliseconds, even at intact-protein
# scale). There is no longer a reason to use Averagine anywhere in this app.

#' Summary statistics (mean shift from the pattern's own lowest-mass peak,
#' and width) of an isotope pattern -- e.g. sequence_isotope_pattern()'s or
#' proteoform_isotope_pattern()'s output. Probability-weighted mean and
#' variance over the binned peaks (the discrete-distribution equivalent of
#' Averagine's binomial-moment shortcut, but computed from the REAL pattern
#' instead of an approximated composition).
#'
#' Takes an already-computed pattern (not a sequence) so callers that also
#' need the full discrete peak list (e.g. building an MS1 chart) can compute
#' the pattern once via the Python bridge and derive both the full list and
#' these summary stats from that single result, rather than paying for two
#' separate reticulate round-trips.
#'
#' @param pattern data.frame(mass, prob), ascending by mass (e.g. from
#'   sequence_isotope_pattern())
#' @return list(mean_shift, sigma, fwhm), all in Da; mean_shift/sigma are
#'   relative to the pattern's own first (lowest-mass, i.e. monoisotopic) row
isotope_envelope_stats <- function(pattern) {
  mono <- pattern$mass[1]
  mean_shift <- sum(pattern$prob * (pattern$mass - mono))
  variance <- sum(pattern$prob * (pattern$mass - mono - mean_shift)^2)
  list(mean_shift = mean_shift, sigma = sqrt(variance), fwhm = 2.3548 * sqrt(variance))
}

#' Decide, for one charge state, whether its isotope peaks would actually be
#' individually resolved by the instrument -- and return exactly the points
#' needed to draw whichever case applies.
#'
#' The physics: a real protein isn't one exact mass -- natural isotopes
#' (mostly the ~1.1% of carbon atoms that are 13C) mean a population of
#' molecules spans a "comb" of peaks roughly 1 Da apart (in mass; 1/z apart
#' in m/z at charge z). Whether the instrument can actually SEE that comb as
#' separate peaks -- rather than one smooth blurred hump -- depends on
#' whether the comb's own tooth-spacing is bigger or smaller than the
#' instrument's resolving power at that m/z (fwhm_mz(), R/resolving_power.R
#' -- the SAME resolving-power model already used for the pairwise
#' resolvability verdict, reused here rather than a separate mass-cutoff
#' rule of thumb). Above roughly 15-20 kDa at typical Orbitrap settings this
#' usually tips into "unresolved," but it's charge-state and
#' resolving-power-setting dependent, not a fixed mass threshold -- which is
#' exactly why this is computed per charge state rather than hardcoded.
#'
#' @param pattern data.frame(mass, prob) neutral-mass isotope pattern (e.g.
#'   proteoform_isotope_pattern())
#' @param z charge state
#' @param r_ref,mz_ref passed to fwhm_mz() (R/resolving_power.R)
#' @return list(z, resolved, points) where points is a data.frame(mz, rel)
#'   -- one row per real isotope peak (rel = that peak's own probability,
#'   normalized to a max of 1) if resolved=TRUE, or a small set of samples
#'   along the smooth Gaussian envelope shape (mean/width from
#'   isotope_envelope_stats() on this same pattern) if resolved=FALSE
ms1_isotope_peaks_for_charge <- function(pattern, z, r_ref = 120000, mz_ref = 200) {
  mz <- mz_for_charge(pattern$mass, z)
  center_mz <- mz_for_charge(sum(pattern$mass * pattern$prob), z)
  # Real (probability-weighted) average spacing between adjacent isotope
  # peaks at this charge state -- not assumed to be exactly 1/z, since real
  # isotope mass differences (13C-12C ~= 1.00336 Da, etc.) aren't exactly
  # 1.0 Da either.
  spacing_mz <- if (length(mz) > 1) mean(diff(mz)) else 1 / z
  resolved <- fwhm_mz(center_mz, r_ref, mz_ref) < spacing_mz

  if (resolved) {
    points <- data.frame(mz = mz, rel = pattern$prob / max(pattern$prob))
  } else {
    stats <- isotope_envelope_stats(pattern)
    mono_mz <- mz_for_charge(pattern$mass[1], z)
    mean_mz <- mono_mz + stats$mean_shift / z
    sigma_mz <- stats$sigma / z
    # A handful of samples across +/-3 sigma is enough to draw a smooth
    # bump shape -- this is explicitly NOT claiming these 25 points are
    # real, individually-observable peaks (they're the opposite case: peaks
    # too close together to resolve), just enough samples of the envelope's
    # own shape for the renderer to draw a curve instead of discrete sticks.
    offsets <- seq(-3, 3, length.out = 25)
    sample_mz <- mean_mz + offsets * sigma_mz
    rel <- exp(-0.5 * offsets^2) # unit-height Gaussian, sigma_mz already applied above
    points <- data.frame(mz = sample_mz, rel = rel)
  }
  list(z = z, resolved = resolved, points = points)
}

#' Full MS1 peak prediction for a proteoform: which charge states are
#' populated and their relative abundance (predict_charge_envelope(),
#' R/charge_envelope.R) combined with, for each of those charge states, its
#' own isotope peaks -- either the real discrete comb if the instrument
#' would resolve it, or a smooth sampled envelope curve if not
#' (ms1_isotope_peaks_for_charge(), above). This is the single entry point
#' the app's MS1 chart(s) should call.
#'
#' @param p a proteoform object (see proteoform_schema.R)
#' @param mode "denatured" or "native"
#' @param average use average mass instead of monoisotopic (must match the
#'   caller's own convention elsewhere in the same analysis)
#' @param r_ref,mz_ref Orbitrap resolving-power settings
#' @return list(mass, charge_states = list(list(z=, resolved=, points=
#'   data.frame(mz, rel)), ...)) -- every point's `rel` is already scaled by
#'   its own charge state's relative abundance, so every point across every
#'   charge state sits on one shared 0-1 intensity axis, directly plottable
#'   without the caller needing to know which charge state it came from
predict_ms1_peaks <- function(p, mode = c("denatured", "native"), average = FALSE,
                               r_ref = 120000, mz_ref = 200) {
  mode <- match.arg(mode)
  mass <- proteoform_mass(p, average = average)$mass
  envelope <- predict_charge_envelope(p$sequence, mass, mode = mode)
  pattern <- proteoform_isotope_pattern(p, average = average)
  charge_states <- Map(function(z, rel) {
    cs <- ms1_isotope_peaks_for_charge(pattern, z, r_ref = r_ref, mz_ref = mz_ref)
    cs$points$rel <- cs$points$rel * rel
    cs
  }, envelope$z, envelope$relative_intensity)
  list(mass = mass, charge_states = unname(charge_states))
}

#' Build a lightweight, synthetic proteoform representing one b/y fragment
#' ion, so predict_ms1_peaks() can be reused unchanged for MS2 isotope
#' peaks: the fragment's own bare subsequence, plus a single synthetic PTM
#' entry (site = 1, never surfaced to the user) carrying TWO things bundled
#' into one flat delta, since predict_ms1_peaks()/proteoform_mass() only
#' ever need the total delta for mass/isotope-pattern purposes, never
#' per-site tracking:
#'   1. the SUM of whatever real PTM deltas fall within this fragment's span
#'   2. the fixed ion-type mass correction: proteoform_mass()/sequence_mass()
#'      compute a free peptide's mass (both termini intact), but a b-ion
#'      neutral fragment is missing the C-terminal OH a real free peptide
#'      would have -- exactly one water lighter, confirmed directly against
#'      sequence_ion_masses_batch()'s own (already-tested) ion-type-aware
#'      masses: for every subsequence tried, free-peptide-mass minus
#'      b-ion-mass was EXACTLY WATER_MASS_MONO (18.0105646863, R/digestion.R),
#'      sequence-independent, and free-peptide-mass minus y-ion-mass was
#'      exactly 0 (a y-ion, having gained a proton on its new N-terminus at
#'      the cleavage site, is chemically identical in mass to a standalone
#'      free peptide of that same subsequence) -- so only b-ions need this
#'      correction. Without it, every b-ion fragment's isotope pattern would
#'      be systematically ~18 Da too heavy.
#'
#' @param full_pf the intact proteoform the fragment ladder was built from
#' @param cleavage_position 1-based interior backbone cleavage position
#' @param ion_type "b" (N-terminal) or "y" (C-terminal)
#' @return a proteoform object for the fragment
fragment_pseudo_proteoform <- function(full_pf, cleavage_position, ion_type = c("b", "y")) {
  ion_type <- match.arg(ion_type)
  n <- nchar(full_pf$sequence)
  if (ion_type == "b") {
    frag_seq <- substr(full_pf$sequence, 1, cleavage_position)
    ptm_delta <- .ptm_delta_for_b_ion(full_pf$ptms, cleavage_position, "mass_delta_mono")
    ion_correction <- -WATER_MASS_MONO
  } else {
    frag_seq <- substr(full_pf$sequence, cleavage_position + 1, n)
    ptm_delta <- .ptm_delta_for_y_ion(full_pf$ptms, cleavage_position, "mass_delta_mono")
    ion_correction <- 0
  }
  delta <- ptm_delta + ion_correction
  ptms <- if (delta != 0) {
    list(ptm(site = 1, mass_delta_mono = delta, mass_delta_avg = delta, name = "cumulative"))
  } else {
    list()
  }
  proteoform(id = paste0(ion_type, cleavage_position), sequence = frag_seq, ptms = ptms, provenance = "manual")
}

