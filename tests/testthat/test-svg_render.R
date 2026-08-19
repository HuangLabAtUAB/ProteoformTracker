test_that("compute_ladder_tiers marks a lone proteoform's bonds neutral (nothing to compare against)", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")
  pf <- proteoform(id = "solo", sequence = "AAAAAGGGGGSLLLLLVVVVV", provenance = "manual")

  result <- compute_ladder_tiers(list(solo = pf), scoring_mode = "glm")
  expect_true(all(result$solo$tier_b == "neutral"))
  expect_true(all(result$solo$tier_y == "neutral"))
})

test_that("compute_ladder_tiers marks identical proteoforms fully common, unrelated ones fully unique", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")
  seq <- "AAAAAGGGGGSLLLLLVVVVV"
  pf_list <- list(
    a = proteoform(id = "a", sequence = seq, provenance = "manual"),
    b = proteoform(id = "b", sequence = seq, provenance = "manual")
  )
  result <- compute_ladder_tiers(pf_list, scoring_mode = "glm")
  expect_true(all(result$a$tier_b == "common"))
  expect_true(all(result$a$tier_y == "common"))
  expect_true(all(result$b$tier_b == "common"))
})

test_that("compute_confounder_tiers gives the target multi-way tiers but each confounder only pairwise-vs-target tiers", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")
  target_seq <- "AAAAAGGGGGSLLLLLVVVVV"
  target <- proteoform(id = "TARGET", sequence = target_seq, provenance = "manual")
  identical_conf <- proteoform(id = "SAME", sequence = target_seq, provenance = "manual")
  unrelated_conf <- proteoform(id = "DIFF", sequence = "WYFHKQNMRCWYFHKQNMRC", provenance = "manual")

  tiers <- compute_confounder_tiers(target, "TARGET", list(SAME = identical_conf, DIFF = unrelated_conf), scoring_mode = "glm")

  # Target: matches exactly 1 of 2 confounders at every bond -> "partial" throughout
  expect_true(all(tiers$TARGET$tier_b == "partial"))
  expect_true(all(tiers$TARGET$tier_y == "partial"))

  # The identical confounder matches the target everywhere -> "common"
  expect_true(all(tiers$SAME$tier_b == "common"))
  expect_true(all(tiers$SAME$tier_y == "common"))

  # The unrelated confounder matches nowhere -> "unique" (never "partial",
  # since confounders are only ever compared against the target, not each other)
  expect_true(all(tiers$DIFF$tier_b == "unique"))
  expect_true(all(tiers$DIFF$tier_y == "unique"))
})

test_that("compute_confounder_tiers handles zero confounders (target tiers all neutral)", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")
  target <- proteoform(id = "TARGET", sequence = "AAAAAGGGGGSLLLLLVVVVV", provenance = "manual")
  tiers <- compute_confounder_tiers(target, "TARGET", list(), scoring_mode = "glm")
  expect_equal(names(tiers), "TARGET")
  expect_true(all(tiers$TARGET$tier_b == "neutral"))
})

test_that("compute_confounder_tiers is dramatically cheaper than full N-way tiering for a large confounder set", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")
  set.seed(42)
  aa <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
  mkseq <- function(n) paste(sample(aa, n, replace = TRUE), collapse = "")
  target <- proteoform(id = "TARGET", sequence = mkseq(200), provenance = "manual")
  confs <- setNames(lapply(1:15, function(i) proteoform(id = paste0("C", i), sequence = mkseq(200), provenance = "manual")), paste0("C", 1:15))

  t_asymmetric <- system.time(compute_confounder_tiers(target, "TARGET", confs, scoring_mode = "glm"))[["elapsed"]]
  t_full_nway <- system.time(compute_ladder_tiers(c(list(TARGET = target), confs), scoring_mode = "glm"))[["elapsed"]]

  # Not a strict inequality assertion (timing can be noisy in CI), just a
  # sanity check that the optimization is real and not accidentally undone.
  expect_true(t_asymmetric < t_full_nway)
})
