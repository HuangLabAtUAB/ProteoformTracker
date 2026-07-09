test_that("unimod_ptm looks up a known modification by name", {
  m <- unimod_ptm("Phospho", site = 5)
  expect_s3_class(m, "ptm")
  expect_equal(m$mass_delta_mono, 79.966331)
  expect_equal(m$mass_delta_avg, 79.9799)
  expect_equal(m$name, "Phospho")
  expect_equal(m$unimod_id, "UNIMOD:21")
  expect_equal(m$site, 5)
})

test_that("unimod_ptm matches case-insensitively", {
  m <- unimod_ptm("phospho", site = 1)
  expect_equal(m$unimod_id, "UNIMOD:21")
  m2 <- unimod_ptm("OXIDATION", site = 1)
  expect_equal(m2$unimod_id, "UNIMOD:35")
})

test_that("unimod_ptm supports terminal sites", {
  m <- unimod_ptm("Acetyl", site = "N-term")
  expect_equal(m$site, "N-term")
})

test_that("unimod_ptm errors with available names for an unknown modification", {
  expect_error(unimod_ptm("NotARealMod", site = 1), "Unknown modification.*NotARealMod")
  expect_error(unimod_ptm("NotARealMod", site = 1), "Phospho")
})

test_that("list_unimod_mods returns the full curated table", {
  tbl <- list_unimod_mods()
  expect_true("Phospho" %in% tbl$name)
  expect_true(all(c("name", "unimod_id", "mass_delta_mono", "mass_delta_avg", "description") %in% names(tbl)))
})

test_that("unimod_ptm() output plugs directly into proteoform_mass()", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")
  p_bare <- proteoform(id = "bare", sequence = "AGCKST", provenance = "manual")
  p_phospho <- proteoform(
    id = "phospho", sequence = "AGCKST",
    ptms = list(unimod_ptm("Phospho", site = 6)),
    provenance = "manual"
  )
  expect_equal(
    proteoform_mass(p_phospho)$mass,
    proteoform_mass(p_bare)$mass + 79.966331,
    tolerance = 1e-6
  )
})
