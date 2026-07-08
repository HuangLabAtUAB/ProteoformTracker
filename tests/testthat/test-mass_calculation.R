test_that("sequence_mass matches a known reference value (Glycine monoisotopic)", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")
  # Glycine (G), monoisotopic mass with water added ~= 75.0320
  expect_equal(sequence_mass("G", average = FALSE), 75.0320, tolerance = 1e-3)
})

test_that("proteoform_mass sums bare-sequence mass and PTM deltas", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")
  p_no_ptm <- proteoform(id = "a", sequence = "AGCK", provenance = "manual")
  p_ptm <- proteoform(
    id = "b", sequence = "AGCK",
    ptms = list(ptm(site = 1, mass_delta_mono = 79.9663, name = "Phospho")),
    provenance = "manual"
  )

  m_no_ptm <- proteoform_mass(p_no_ptm)
  m_ptm <- proteoform_mass(p_ptm)

  expect_equal(m_no_ptm$ptm_delta, 0)
  expect_equal(m_ptm$ptm_delta, 79.9663)
  expect_equal(m_ptm$mass, m_no_ptm$mass + 79.9663, tolerance = 1e-6)
  expect_equal(m_ptm$sequence_mass, m_no_ptm$sequence_mass)
})

test_that("average mass is greater than monoisotopic mass for the same sequence", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")
  mono <- sequence_mass("MAGCKWERTY", average = FALSE)
  avg <- sequence_mass("MAGCKWERTY", average = TRUE)
  expect_gt(avg, mono)
})

test_that("sequence_masses_batch matches per-sequence sequence_mass", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")
  seqs <- c("G", "AGCK", "MAGCKWERTY")
  batch <- sequence_masses_batch(seqs, script_path = "../../python/ptracker_mass.py")
  singles <- vapply(seqs, sequence_mass, numeric(1))
  expect_equal(unname(batch), unname(singles), tolerance = 1e-6)
})
