test_proteoform_sequence <- function(n_repeats = 20) {
  paste(rep("AGCKRKWERTYAGCKRKWERTY", n_repeats), collapse = "")
}

test_that("ms1_resolvability classifies resolvable/marginal/not-resolvable relative to best-case FWHM", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")

  seq <- test_proteoform_sequence()
  target <- proteoform(id = "target", sequence = seq, provenance = "manual")
  target_mass <- proteoform_mass(target)$mass

  envelope <- predict_charge_envelope(seq, target_mass, mode = "denatured")
  fwhm_tbl <- fwhm_by_charge_state(target_mass, envelope$z)
  best_fwhm <- min(fwhm_tbl$fwhm_mass)

  make_candidate <- function(id, delta) {
    proteoform(
      id = id, sequence = seq,
      ptms = list(ptm(site = 1, mass_delta_mono = delta, name = "test_delta")),
      provenance = "manual"
    )
  }

  result_resolvable <- ms1_resolvability(target, make_candidate("cand_resolvable", best_fwhm * 5))
  expect_equal(result_resolvable$verdict, "resolvable")

  result_marginal <- ms1_resolvability(target, make_candidate("cand_marginal", best_fwhm * 1.2))
  expect_equal(result_marginal$verdict, "marginal")

  result_not <- ms1_resolvability(target, make_candidate("cand_not", best_fwhm * 0.1))
  expect_equal(result_not$verdict, "not-resolvable")
})

test_that("ms1_resolvability flags envelope_interleave_risk above 25 kDa", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")

  small_seq <- test_proteoform_sequence(2) # small enough to stay under 25 kDa
  large_seq <- test_proteoform_sequence(40) # comfortably over 25 kDa

  small_target <- proteoform(id = "small", sequence = small_seq, provenance = "manual")
  small_candidate <- proteoform(
    id = "small_cand", sequence = small_seq,
    ptms = list(ptm(site = 1, mass_delta_mono = 5, name = "d")), provenance = "manual"
  )
  large_target <- proteoform(id = "large", sequence = large_seq, provenance = "manual")
  large_candidate <- proteoform(
    id = "large_cand", sequence = large_seq,
    ptms = list(ptm(site = 1, mass_delta_mono = 5, name = "d")), provenance = "manual"
  )

  expect_false(ms1_resolvability(small_target, small_candidate)$envelope_interleave_risk)
  expect_true(ms1_resolvability(large_target, large_candidate)$envelope_interleave_risk)
})

test_that("envelope_crowding_check flags identical-mass peaks from different proteoforms", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")

  seq <- test_proteoform_sequence()
  p1 <- proteoform(id = "p1", sequence = seq, provenance = "manual")
  p2 <- proteoform(id = "p2", sequence = seq, provenance = "manual") # identical mass/sequence

  result <- envelope_crowding_check(list(p1, p2), mode = "denatured")
  expect_gt(nrow(result$flags), 0)
  expect_true(all(result$flags$distance == 0))
})

test_that("envelope_crowding_check reports no flags when masses are far apart", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")

  small_seq <- test_proteoform_sequence(2)
  large_seq <- test_proteoform_sequence(40)
  p1 <- proteoform(id = "small", sequence = small_seq, provenance = "manual")
  p2 <- proteoform(id = "large", sequence = large_seq, provenance = "manual")

  result <- envelope_crowding_check(list(p1, p2), mode = "denatured")
  expect_equal(nrow(result$flags), 0)
})
