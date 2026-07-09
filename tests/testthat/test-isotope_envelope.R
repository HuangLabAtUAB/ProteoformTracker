test_that("averagine_composition scales roughly linearly with mass", {
  comp_small <- averagine_composition(16000)
  comp_large <- averagine_composition(32000) # 2x mass
  expect_equal(unname(comp_large["C"]), unname(comp_small["C"]) * 2, tolerance = 2)
})

test_that("averagine_isotope_envelope gives plausible mean shift and width for a 16 kDa protein", {
  env <- averagine_isotope_envelope(16000)
  # apex of the isotope envelope for a ~16 kDa protein is typically several Da
  # above monoisotopic -- a sanity range, not a precise literature value
  expect_gt(env$mean_shift, 5)
  expect_lt(env$mean_shift, 15)
  expect_gt(env$fwhm, 0)
  expect_equal(env$fwhm, 2.3548 * env$sigma)
})

test_that("isotope envelope width grows with mass (wider envelopes for bigger proteins)", {
  small <- averagine_isotope_envelope(10000)
  large <- averagine_isotope_envelope(50000)
  expect_gt(large$fwhm, small$fwhm)
  expect_gt(large$mean_shift, small$mean_shift)
})
