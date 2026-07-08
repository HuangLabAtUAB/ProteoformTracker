test_that("resolving_power matches the design spec worked example", {
  # 16 kDa proteoform, z=20 -> m/z ~= 801; R(801) ~= 60,000 at R_ref=120,000 @ mz_ref=200
  mz <- mz_for_charge(16000, 20)
  expect_equal(mz, 801.0, tolerance = 0.1)

  r <- resolving_power(mz, r_ref = 120000, mz_ref = 200)
  expect_equal(r, 60000, tolerance = 1000)
})

test_that("fwhm_mass matches the design spec worked example (~0.27 Da)", {
  mz <- mz_for_charge(16000, 20)
  fwhm <- fwhm_mass(16000, mz, r_ref = 120000, mz_ref = 200)
  expect_equal(fwhm, 0.27, tolerance = 0.01)
})

test_that("resolving_power degrades (decreases) as m/z increases", {
  expect_gt(resolving_power(200), resolving_power(800))
  expect_gt(resolving_power(800), resolving_power(2000))
})

test_that("resolving_power at the reference m/z equals r_ref", {
  expect_equal(resolving_power(200, r_ref = 120000, mz_ref = 200), 120000)
})

test_that("higher charge state gives smaller (better) mass-domain FWHM", {
  mass <- 16000
  fwhm_low_z <- fwhm_mass(mass, mz_for_charge(mass, 8))
  fwhm_high_z <- fwhm_mass(mass, mz_for_charge(mass, 20))
  expect_lt(fwhm_high_z, fwhm_low_z)
})

test_that("fwhm_by_charge_state identifies the highest charge state as best", {
  tbl <- fwhm_by_charge_state(16000, charge_states = c(8, 12, 16, 20))
  best <- attr(tbl, "best")
  expect_equal(tbl$z[best], 20)
  expect_equal(tbl$fwhm_mass[best], min(tbl$fwhm_mass))
})
