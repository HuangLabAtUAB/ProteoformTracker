test_that("residue_pair_propensity is baseline when neither Pro nor Asp effect applies", {
  expect_equal(residue_pair_propensity("A", "A", mode = "denatured"), 1.0)
})

test_that("residue_pair_propensity enhances cleavage N-terminal to proline", {
  score <- residue_pair_propensity("A", "P", mode = "denatured")
  expect_equal(score, PROLINE_ENHANCEMENT$denatured)
})

test_that("residue_pair_propensity enhances cleavage C-terminal to aspartate", {
  score <- residue_pair_propensity("D", "A", mode = "denatured")
  expect_equal(score, ASPARTATE_ENHANCEMENT$denatured)
})

test_that("residue_pair_propensity multiplies both enhancements when a bond is D|P", {
  score <- residue_pair_propensity("D", "P", mode = "denatured")
  expect_equal(score, ASPARTATE_ENHANCEMENT$denatured * PROLINE_ENHANCEMENT$denatured)
})

test_that("residue_pair_propensity enhancement is stronger under native than denatured", {
  denatured_score <- residue_pair_propensity("D", "A", mode = "denatured")
  native_score <- residue_pair_propensity("D", "A", mode = "native")
  expect_gt(native_score, denatured_score)
})

test_that("positional_propensity matches the calibrated lookup table by terminal distance", {
  n <- 100
  expect_equal(positional_propensity(1, n), TERMINAL_PROXIMITY_WEIGHTS[1])
  expect_equal(positional_propensity(2, n), TERMINAL_PROXIMITY_WEIGHTS[2])
  expect_equal(positional_propensity(3, n), TERMINAL_PROXIMITY_WEIGHTS[3])
  expect_equal(positional_propensity(4, n), TERMINAL_PROXIMITY_WEIGHTS[4])
  # beyond the calibrated table's reach, clamps to the last (6+) entry
  expect_equal(positional_propensity(50, n), TERMINAL_PROXIMITY_WEIGHTS[6])
  expect_equal(positional_propensity(6, n), TERMINAL_PROXIMITY_WEIGHTS[6])
  # symmetric from the C-terminal side too
  expect_equal(positional_propensity(n - 1, n), TERMINAL_PROXIMITY_WEIGHTS[1])
})

test_that("positional_propensity peaks a few residues in from the terminus, not right at it", {
  n <- 100
  # distance-1 (touching the terminus) is the calibrated WEAKEST position
  expect_lt(positional_propensity(1, n), positional_propensity(2, n))
  expect_lt(positional_propensity(2, n), positional_propensity(3, n))
  expect_lt(positional_propensity(3, n), positional_propensity(4, n))
  # ... and relaxes back down past the peak toward the interior plateau
  expect_gt(positional_propensity(5, n), positional_propensity(50, n))
})

test_that("charge_density_propensity multiplier decreases with local K/R density", {
  residues <- strsplit("AAAAAAAAAAA", "")[[1]] # no basic residues nearby
  expect_equal(charge_density_propensity(residues, 6), CHARGE_DENSITY_WEIGHTS$none)

  residues_1k <- strsplit("AAAAAKAAAAA", "")[[1]] # exactly one K in the window
  expect_equal(charge_density_propensity(residues_1k, 6), CHARGE_DENSITY_WEIGHTS$low)

  residues_2k <- strsplit("AAAAKAKAAAA", "")[[1]] # two basic residues nearby
  expect_equal(charge_density_propensity(residues_2k, 6), CHARGE_DENSITY_WEIGHTS$high)
})

test_that("length_propensity favors shorter proteoforms and clamps at both ends", {
  expect_equal(length_propensity(LENGTH_REFERENCE), 1.0)
  expect_gt(length_propensity(50), length_propensity(200)) # shorter -> higher multiplier
  # clamped below LENGTH_CLAMP_MIN and above LENGTH_CLAMP_MAX
  expect_equal(length_propensity(5), length_propensity(LENGTH_CLAMP_MIN))
  expect_equal(length_propensity(10000), length_propensity(LENGTH_CLAMP_MAX))
})

test_that("phospho_proximity_propensity suppresses bonds near a phosphosite and is otherwise a no-op", {
  expect_equal(phospho_proximity_propensity(list(), 10), BASELINE_WEIGHT)

  far <- list(ptm(site = 50, mass_delta_mono = 79.966331, name = "Phospho", unimod_id = "UNIMOD:21"))
  expect_equal(phospho_proximity_propensity(far, 10), BASELINE_WEIGHT)

  near <- list(ptm(site = 12, mass_delta_mono = 79.966331, name = "Phospho", unimod_id = "UNIMOD:21"))
  expect_equal(phospho_proximity_propensity(near, 10), PHOSPHO_SUPPRESSION)

  # a non-phospho PTM nearby doesn't trigger suppression
  other <- list(ptm(site = 10, mass_delta_mono = 42.010565, name = "Acetyl", unimod_id = "UNIMOD:1"))
  expect_equal(phospho_proximity_propensity(other, 10), BASELINE_WEIGHT)

  # terminal-site PTMs have no residue position and are safely ignored
  terminal <- list(ptm(site = "N-term", mass_delta_mono = 79.966331, name = "Phospho", unimod_id = "UNIMOD:21"))
  expect_equal(phospho_proximity_propensity(terminal, 10), BASELINE_WEIGHT)
})

test_that("fragmentation_propensity extracts residue pairs correctly and combines effects multiplicatively", {
  # no K/R anywhere in this sequence, so charge_density_score is a 1.0
  # no-op throughout and residue-pair/positional effects can be checked
  # in isolation: i=6 D|P (both chemistry effects), i=10 A|P (Pro only)
  seq <- "AAAAADPAAAPAAAA"
  p <- proteoform(id = "test", sequence = seq, provenance = "manual")
  result <- fragmentation_propensity(p, mode = "denatured", method = "HCD")

  expect_equal(nrow(result), nchar(seq) - 1)
  expect_true(all(result$charge_density_score == CHARGE_DENSITY_WEIGHTS$none))
  expect_true(all(result$length_score == length_propensity(nchar(seq))))
  expect_true(all(result$phospho_score == BASELINE_WEIGHT)) # no ptms on this proteoform

  bond6 <- result[result$cleavage_position == 6, ]
  expect_equal(bond6$residue_before, "D")
  expect_equal(bond6$residue_after, "P")
  expect_equal(bond6$positional_score, positional_propensity(6, nchar(seq)))
  expect_equal(bond6$residue_pair_score, ASPARTATE_ENHANCEMENT$denatured * PROLINE_ENHANCEMENT$denatured)
  expect_equal(bond6$propensity_score,
    bond6$residue_pair_score * bond6$positional_score * bond6$charge_density_score *
      bond6$length_score * bond6$phospho_score * bond6$accessibility_score)

  bond10 <- result[result$cleavage_position == 10, ]
  expect_equal(bond10$residue_before, "A")
  expect_equal(bond10$residue_after, "P")
  expect_equal(bond10$residue_pair_score, PROLINE_ENHANCEMENT$denatured)

  bond8 <- result[result$cleavage_position == 8, ]
  expect_equal(bond8$residue_pair_score, 1.0)
  expect_equal(bond8$propensity_score, positional_propensity(8, nchar(seq)) * length_propensity(nchar(seq)))
})

test_that("fragmentation_propensity applies phospho suppression only near an actual phosphosite", {
  seq <- "AAAAADPAAAPAAAA" # 15 residues
  ptms <- list(ptm(site = 8, mass_delta_mono = 79.966331, name = "Phospho", unimod_id = "UNIMOD:21"))
  p <- proteoform(id = "test", sequence = seq, ptms = ptms, provenance = "manual")
  result <- fragmentation_propensity(p, mode = "denatured", method = "HCD")

  near <- result[result$cleavage_position == 8, ]
  expect_equal(near$phospho_score, PHOSPHO_SUPPRESSION)

  far <- result[result$cleavage_position == 1, ]
  expect_equal(far$phospho_score, BASELINE_WEIGHT)
})

test_that("fragmentation_propensity's accessibility term is always a 1.0 no-op for HCD/CID", {
  p <- proteoform(id = "test", sequence = "AAAAADPAAAPAAAA", provenance = "manual")
  result <- fragmentation_propensity(p, method = "HCD")
  expect_true(all(result$accessibility_score == 1.0))
})

test_that("fragmentation_propensity rejects unimplemented dissociation methods", {
  p <- proteoform(id = "test", sequence = "AAAAADPAAAPAAAA", provenance = "manual")
  expect_error(fragmentation_propensity(p, method = "ETD"), "HCD/CID")
  expect_error(fragmentation_propensity(p, method = "UVPD"), "HCD/CID")
})

test_that("fragmentation_propensity rejects sequences shorter than 2 residues", {
  p <- proteoform(id = "short", sequence = "A", provenance = "manual")
  expect_error(fragmentation_propensity(p), "at least 2 residues")
})

test_that("fragmentation_propensity requires a proteoform object", {
  expect_error(fragmentation_propensity(list(sequence = "AAAA")), "requires a proteoform")
})
