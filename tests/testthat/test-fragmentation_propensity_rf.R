test_that("fragmentation_propensity_rf requires a proteoform object", {
  skip_if_not(rf_model_available, "RF model not built/available in this environment")
  expect_error(fragmentation_propensity_rf(list(sequence = "AAAA"), rf_model = propensity_rf_model),
               "requires a proteoform")
})

test_that("fragmentation_propensity_rf errors clearly when no model is available", {
  p <- proteoform(id = "test", sequence = "AAAAADPAAAPAAAA", provenance = "manual")
  expect_error(fragmentation_propensity_rf(p, rf_model = NULL), "RF ranking mode unavailable")
})

test_that("fragmentation_propensity_rf rejects sequences shorter than 2 residues", {
  skip_if_not(rf_model_available, "RF model not built/available in this environment")
  p <- proteoform(id = "short", sequence = "A", provenance = "manual")
  expect_error(fragmentation_propensity_rf(p, rf_model = propensity_rf_model), "at least 2 residues")
})

test_that("fragmentation_propensity_rf returns one row per bond with scores in [0,1]", {
  skip_if_not(rf_model_available, "RF model not built/available in this environment")
  seq <- "AAAAADPAAAPAAAA"
  p <- proteoform(id = "test", sequence = seq, provenance = "manual")
  result <- fragmentation_propensity_rf(p, rf_model = propensity_rf_model)

  expect_equal(nrow(result), nchar(seq) - 1)
  expect_true(all(result$propensity_score >= 0 & result$propensity_score <= 1))
  expect_setequal(names(result), c("cleavage_position", "residue_before", "residue_after",
                                    "terminal_distance", "basic_density", "near_phospho", "propensity_score"))
})

test_that("fragmentation_propensity_rf's raw features match the GLM path's underlying calculations", {
  skip_if_not(rf_model_available, "RF model not built/available in this environment")
  seq <- "AAAAAKAAAADPAAAA"
  p <- proteoform(id = "test", sequence = seq, provenance = "manual")
  result <- fragmentation_propensity_rf(p, rf_model = propensity_rf_model)

  residues <- strsplit(seq, "")[[1]]
  bond6 <- result[result$cleavage_position == 6, ]
  expect_equal(bond6$basic_density, local_basic_density(residues, 6))
  expect_equal(bond6$terminal_distance, min(6, nchar(seq) - 6))
})

test_that("fragmentation_propensity_rf is not affected by proteoform length the way the GLM is", {
  # the whole point of this mode: a bond's score shouldn't collapse just
  # because the surrounding protein is long, unlike length_propensity()
  skip_if_not(rf_model_available, "RF model not built/available in this environment")
  short_seq <- paste0("AAAAADPAAAA", strrep("A", 20))
  long_seq <- paste0("AAAAADPAAAA", strrep("A", 600))
  p_short <- proteoform(id = "short", sequence = short_seq, provenance = "manual")
  p_long <- proteoform(id = "long", sequence = long_seq, provenance = "manual")

  score_short <- fragmentation_propensity_rf(p_short, rf_model = propensity_rf_model)
  score_long <- fragmentation_propensity_rf(p_long, rf_model = propensity_rf_model)

  # same D|P bond near the N-terminus in both -- RF score should be
  # identical or very close regardless of overall protein length, since
  # terminal_distance/basic_density/near_phospho near this bond don't
  # depend on total sequence length once the window is fully in-bounds
  bond6_short <- score_short[score_short$cleavage_position == 6, ]
  bond6_long <- score_long[score_long$cleavage_position == 6, ]
  expect_equal(bond6_short$propensity_score, bond6_long$propensity_score)
})
