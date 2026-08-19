test_that("isotope_envelope_stats gives plausible mean shift and width for a ~16 kDa protein", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")

  seq <- paste(rep("AGCKRKWERTYAGCKRKWERTY", 6), collapse = "") # ~16.5 kDa
  pattern <- sequence_isotope_pattern(seq)
  stats <- isotope_envelope_stats(pattern)
  # apex of the isotope envelope for a ~16 kDa protein is typically several
  # Da above monoisotopic -- a sanity range, not a precise literature value
  expect_gt(stats$mean_shift, 5)
  expect_lt(stats$mean_shift, 15)
  expect_gt(stats$fwhm, 0)
  expect_equal(stats$fwhm, 2.3548 * stats$sigma)
})

test_that("isotope envelope width grows with sequence length (wider envelopes for bigger proteins)", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")

  small <- isotope_envelope_stats(sequence_isotope_pattern(paste(rep("AGCKRKWERTYAGCKRKWERTY", 12), collapse = "")))
  large <- isotope_envelope_stats(sequence_isotope_pattern(paste(rep("AGCKRKWERTYAGCKRKWERTY", 60), collapse = "")))
  expect_gt(large$fwhm, small$fwhm)
  expect_gt(large$mean_shift, small$mean_shift)
})

test_that("sequence_isotope_pattern matches sequence_mass at its own monoisotopic (lowest-mass) peak", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")

  seq <- paste(rep("AGCKRKWERTYAGCKRKWERTY", 20), collapse = "")
  pattern <- sequence_isotope_pattern(seq)
  expect_equal(pattern$mass[1], sequence_mass(seq, average = FALSE), tolerance = 1)
  expect_equal(sum(pattern$prob), 1, tolerance = 1e-6)
  expect_true(all(diff(pattern$mass) > 0)) # strictly ascending, one row per nominal-mass bin
})

test_that("proteoform_isotope_pattern shifts every peak by the proteoform's total PTM delta", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")

  seq <- paste(rep("AGCKRKWERTYAGCKRKWERTY", 10), collapse = "")
  bare <- proteoform(id = "bare", sequence = seq, provenance = "manual")
  modified <- proteoform(
    id = "modified", sequence = seq,
    ptms = list(ptm(site = 1, mass_delta_mono = 79.9663, name = "Phospho")),
    provenance = "manual"
  )
  bare_pattern <- proteoform_isotope_pattern(bare)
  modified_pattern <- proteoform_isotope_pattern(modified)
  expect_equal(modified_pattern$mass, bare_pattern$mass + 79.9663, tolerance = 1e-6)
  expect_equal(modified_pattern$prob, bare_pattern$prob) # shape unchanged, only position shifts
})

test_that("ms1_isotope_peaks_for_charge resolves individual isotope peaks at low mass/high resolving power", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")

  seq <- paste(rep("AGCKRKWERTYAGCKRKWERTY", 5), collapse = "") # small, ~2.7 kDa
  pattern <- sequence_isotope_pattern(seq)
  result <- ms1_isotope_peaks_for_charge(pattern, z = 1, r_ref = 120000, mz_ref = 200)
  expect_true(result$resolved)
  expect_equal(nrow(result$points), nrow(pattern)) # one plotted point per real isotope peak
  expect_equal(max(result$points$rel), 1) # normalized to the most-abundant peak
})

test_that("ms1_isotope_peaks_for_charge falls back to a smooth envelope at low resolving power", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")

  seq <- paste(rep("AGCKRKWERTYAGCKRKWERTY", 5), collapse = "")
  pattern <- sequence_isotope_pattern(seq)
  result <- ms1_isotope_peaks_for_charge(pattern, z = 1, r_ref = 1000, mz_ref = 200) # deliberately low R
  expect_false(result$resolved)
  expect_gt(nrow(result$points), 1) # still multiple points (a sampled curve), just not real discrete peaks
})

test_that("predict_ms1_peaks returns one charge-state entry per populated charge state, each on a 0-1 intensity axis", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")

  seq <- paste(rep("AGCKRKWERTYAGCKRKWERTY", 5), collapse = "")
  p <- proteoform(id = "test", sequence = seq, provenance = "manual")
  envelope <- predict_charge_envelope(p$sequence, proteoform_mass(p)$mass, mode = "denatured")
  result <- predict_ms1_peaks(p, mode = "denatured")

  expect_equal(length(result$charge_states), length(envelope$z))
  # every point's rel is already scaled by its own charge state's relative
  # abundance -- the MOST intense point across the whole chart should equal
  # the envelope's own max relative intensity, not 1.0 for every charge
  # state independently (that would defeat the point of putting every
  # charge state on one shared, comparable intensity axis).
  all_rel <- unlist(lapply(result$charge_states, function(cs) cs$points$rel))
  expect_equal(max(all_rel), max(envelope$relative_intensity), tolerance = 1e-6)
})

test_that("fragment_pseudo_proteoform's mass exactly matches generate_fragment_ladder()'s own b/y masses", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")

  seq <- paste(rep("AGCKRKWERTYAGCKRKWERTY", 10), collapse = "") # ~220 aa
  bare <- proteoform(id = "bare", sequence = seq, provenance = "manual")
  modified <- proteoform(
    id = "modified", sequence = seq,
    ptms = list(ptm(site = 5, mass_delta_mono = 14.0157, name = "Methyl")),
    provenance = "manual"
  )
  ladder_bare <- generate_fragment_ladder(bare)
  ladder_mod <- generate_fragment_ladder(modified)

  for (p in c(10, 100, nchar(seq) - 1)) {
    for (ion in c("b", "y")) {
      correct_bare <- if (ion == "b") ladder_bare$b_mass[p] else ladder_bare$y_mass[p]
      correct_mod <- if (ion == "b") ladder_mod$b_mass[p] else ladder_mod$y_mass[p]
      got_bare <- proteoform_mass(fragment_pseudo_proteoform(bare, p, ion))$mass
      got_mod <- proteoform_mass(fragment_pseudo_proteoform(modified, p, ion))$mass
      expect_equal(got_bare, correct_bare, tolerance = 1e-6)
      expect_equal(got_mod, correct_mod, tolerance = 1e-6)
    }
  }
})

