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

test_that("ms1_resolvability flags envelope_interleave_risk when Delta-mass is smaller than the isotope envelope width", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")

  seq <- test_proteoform_sequence()
  target <- proteoform(id = "target", sequence = seq, provenance = "manual")
  target_mass <- proteoform_mass(target)$mass
  isotope_fwhm <- averagine_isotope_envelope(target_mass)$fwhm

  make_candidate <- function(id, delta) {
    proteoform(
      id = id, sequence = seq,
      ptms = list(ptm(site = 1, mass_delta_mono = delta, name = "d")),
      provenance = "manual"
    )
  }

  # Delta-mass well inside the isotope envelope's own width -> interleaving risk
  small_delta_result <- ms1_resolvability(target, make_candidate("small_delta", isotope_fwhm * 0.1))
  expect_true(small_delta_result$envelope_interleave_risk)

  # Delta-mass well beyond the isotope envelope's width -> no interleaving risk
  large_delta_result <- ms1_resolvability(target, make_candidate("large_delta", isotope_fwhm * 5))
  expect_false(large_delta_result$envelope_interleave_risk)
})

test_that("isotope envelope interleaving risk can trip well below 25-30 kDa for small Delta-mass pairs", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")

  seq <- test_proteoform_sequence(6) # a modest-size proteoform, well under 25 kDa
  target <- proteoform(id = "small_target", sequence = seq, provenance = "manual")
  candidate <- proteoform(
    id = "small_candidate", sequence = seq,
    ptms = list(ptm(site = 1, mass_delta_mono = 0.5, name = "d")), provenance = "manual"
  )
  result <- ms1_resolvability(target, candidate)
  expect_lt(result$target_mass, 25000)
  expect_true(result$envelope_interleave_risk)
})

test_that("search_confounding_proteins sizes its window from FWHM x safety_margin and excludes the target", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")

  seq <- test_proteoform_sequence()
  target <- proteoform(id = "target", sequence = seq, provenance = "manual")
  target_mass <- proteoform_mass(target)$mass
  best_fwhm <- min(fwhm_by_charge_state(target_mass, predict_charge_envelope(seq, target_mass, "denatured")$z)$fwhm_mass)
  safety_margin <- DEFAULT_SAFETY_MARGIN
  window <- safety_margin * best_fwhm

  mass_index <- data.frame(
    id = c("target", "just_inside", "just_outside", "far_away"),
    sequence = seq,
    length = nchar(seq),
    mass = c(target_mass, target_mass + window * 0.5, target_mass + window * 3, target_mass + 5000),
    stringsAsFactors = FALSE
  )
  mass_index <- mass_index[order(mass_index$mass), ]

  result <- search_confounding_proteins(target, mass_index, mode = "denatured", safety_margin = safety_margin)

  expect_equal(result$window_da, window, tolerance = 1e-9)
  expect_false("target" %in% result$candidates$id)
  expect_true("just_inside" %in% result$candidates$id)
  expect_false("just_outside" %in% result$candidates$id)
  expect_false("far_away" %in% result$candidates$id)
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
