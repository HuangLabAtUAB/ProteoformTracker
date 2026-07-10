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

test_that("positional_propensity is elevated within the terminal zone and flat beyond it", {
  n <- 100
  expect_equal(positional_propensity(1, n), TERMINAL_PROXIMITY_WEIGHTS[1])
  expect_equal(positional_propensity(2, n), TERMINAL_PROXIMITY_WEIGHTS[2])
  expect_equal(positional_propensity(3, n), TERMINAL_PROXIMITY_WEIGHTS[3])
  expect_equal(positional_propensity(4, n), 1.0)
  expect_equal(positional_propensity(50, n), 1.0)
  # symmetric from the C-terminal side too
  expect_equal(positional_propensity(n - 1, n), TERMINAL_PROXIMITY_WEIGHTS[1])
})

test_that("positional_propensity weights decay monotonically within the terminal zone", {
  n <- 100
  expect_gt(positional_propensity(1, n), positional_propensity(2, n))
  expect_gt(positional_propensity(2, n), positional_propensity(3, n))
  expect_gt(positional_propensity(3, n), positional_propensity(4, n))
})

test_that("fragmentation_propensity extracts residue pairs correctly and combines effects multiplicatively", {
  # bonds of interest all placed >3 residues from either terminus, so
  # positional_score is baseline (1.0) and residue-pair effects can be
  # checked in isolation: i=6 D|P (both effects), i=10 A|P (Pro only)
  seq <- "AAAAADPAAAPAAAA"
  p <- proteoform(id = "test", sequence = seq, provenance = "manual")
  result <- fragmentation_propensity(p, mode = "denatured", method = "HCD")

  expect_equal(nrow(result), nchar(seq) - 1)

  bond6 <- result[result$cleavage_position == 6, ]
  expect_equal(bond6$residue_before, "D")
  expect_equal(bond6$residue_after, "P")
  expect_equal(bond6$positional_score, 1.0)
  expect_equal(bond6$residue_pair_score, ASPARTATE_ENHANCEMENT$denatured * PROLINE_ENHANCEMENT$denatured)
  expect_equal(bond6$propensity_score, bond6$residue_pair_score * bond6$positional_score * bond6$accessibility_score)

  bond10 <- result[result$cleavage_position == 10, ]
  expect_equal(bond10$residue_before, "A")
  expect_equal(bond10$residue_after, "P")
  expect_equal(bond10$residue_pair_score, PROLINE_ENHANCEMENT$denatured)

  bond8 <- result[result$cleavage_position == 8, ]
  expect_equal(bond8$residue_pair_score, 1.0)
  expect_equal(bond8$propensity_score, 1.0)
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
