test_that("basic_residue_count counts R/K plus N-terminal amine", {
  # 2x K, 1x R -> 3 basic residues + 1 N-term amine = 4
  expect_equal(basic_residue_count("MAGCKRK"), 4)
  expect_equal(basic_residue_count("MAGCKRK", include_n_terminal_amine = FALSE), 3)
})

test_that("basic_residue_count floors at 1 for sequences with no basic residues", {
  expect_equal(basic_residue_count("AGCST", include_n_terminal_amine = FALSE), 1)
})

test_that("rayleigh_charge_ceiling scales with sqrt(mass)", {
  z_small <- rayleigh_charge_ceiling(16000)
  z_large <- rayleigh_charge_ceiling(64000) # 4x mass -> 2x ceiling
  expect_equal(z_large, 2 * z_small, tolerance = 1)
})

test_that("charge_state_ceiling uses basic-residue count for small denatured proteoforms", {
  seq <- "MAGCKRKAGCKRK" # several basic residues
  ceiling_denatured <- charge_state_ceiling(seq, mass = 16000, mode = "denatured")
  expect_equal(ceiling_denatured, basic_residue_count(seq))
})

test_that("charge_state_ceiling uses Rayleigh limit above 40 kDa even if denatured", {
  seq <- "MAGCKRKAGCKRK"
  ceiling_large <- charge_state_ceiling(seq, mass = 80000, mode = "denatured")
  expect_equal(ceiling_large, rayleigh_charge_ceiling(80000))
})

test_that("charge_state_ceiling uses Rayleigh limit for native mode regardless of size", {
  seq <- "MAGCKRKAGCKRK"
  ceiling_native <- charge_state_ceiling(seq, mass = 16000, mode = "native")
  expect_equal(ceiling_native, rayleigh_charge_ceiling(16000))
})

test_that("predict_charge_envelope returns a peak-normalized intensity envelope", {
  seq <- paste(rep("MAGCKRK", 20), collapse = "")
  mass <- 16000
  envelope <- predict_charge_envelope(seq, mass, mode = "denatured")
  expect_true(all(c("z", "mz", "relative_intensity") %in% names(envelope)))
  expect_equal(max(envelope$relative_intensity), 1)
  expect_true(all(envelope$z >= 1))
  expect_true(all(diff(envelope$z) == 1)) # contiguous charge states
})

test_that("native envelopes are narrower than denatured envelopes for the same proteoform", {
  seq <- paste(rep("MAGCKRK", 20), collapse = "")
  mass <- 16000
  denatured <- predict_charge_envelope(seq, mass, mode = "denatured")
  native <- predict_charge_envelope(seq, mass, mode = "native")
  expect_lt(nrow(native), nrow(denatured))
})

test_that("sn_relative_penalty equals 1 at the baseline mass and decreases as mass grows", {
  expect_equal(sn_relative_penalty(10000, baseline_mass = 10000, mode = "denatured"), 1)
  expect_lt(sn_relative_penalty(20000, baseline_mass = 10000, mode = "denatured"), 1)
})

test_that("sn_relative_penalty decays faster under denatured than native conditions", {
  denatured_penalty <- sn_relative_penalty(40000, baseline_mass = 10000, mode = "denatured")
  native_penalty <- sn_relative_penalty(40000, baseline_mass = 10000, mode = "native")
  expect_lt(denatured_penalty, native_penalty)
})
