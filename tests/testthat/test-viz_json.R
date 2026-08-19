test_that("tally_ms2_tiers counts each tier across both ion series", {
  tier_b <- c("unique", "common", "partial", "unique")
  tier_y <- c("common", "common", "unique", "neutral")
  result <- tally_ms2_tiers(tier_b, tier_y)

  expect_equal(result$total, 8)
  expect_equal(result$unique, 3)
  expect_equal(result$common, 3)
  expect_equal(result$partial, 1)
  expect_equal(result$neutral, 1)
})

test_that("tally_ms2_tiers handles the single-checked-proteoform (all-neutral) case", {
  result <- tally_ms2_tiers(rep("neutral", 5), rep("neutral", 5))
  expect_equal(result$total, 10)
  expect_equal(result$neutral, 10)
  expect_equal(result$unique, 0)
})
