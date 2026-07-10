test_that("generate_fragment_ladder produces one row per interior cleavage position", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")
  p <- proteoform(id = "test", sequence = "PEPTIDE", provenance = "manual")
  ladder <- generate_fragment_ladder(p, script_path = TEST_PY_SCRIPT)
  expect_equal(nrow(ladder), 6) # n=7 -> 6 interior cleavage positions
  expect_equal(ladder$cleavage_position, 1:6)
  expect_equal(ladder$n_term_length, 1:6)
  expect_equal(ladder$c_term_length, 6:1)
})

test_that("complementary bare b/y neutral masses sum to the intact bare-sequence mass", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")
  p <- proteoform(id = "test", sequence = "PEPTIDE", provenance = "manual")
  ladder <- generate_fragment_ladder(p, script_path = TEST_PY_SCRIPT)
  intact_mass <- sequence_mass(p$sequence)

  # b_i (neutral) + y_(n-i) (neutral) == M (neutral), for every complementary pair
  for (i in ladder$cleavage_position) {
    b_i <- ladder$b_mass[ladder$cleavage_position == i]
    y_complement <- ladder$y_mass[ladder$cleavage_position == i] # y at the SAME row is y_(n-i)
    expect_equal(b_i + y_complement, intact_mass, tolerance = 1e-6)
  }
})

test_that("generate_fragment_ladder matches direct pyteomics ion-mass calls for known subsequences", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")
  p <- proteoform(id = "test", sequence = "PEPTIDE", provenance = "manual")
  ladder <- generate_fragment_ladder(p, script_path = TEST_PY_SCRIPT)

  b3_expected <- sequence_ion_masses_batch("PEP", ion_type = "b", script_path = TEST_PY_SCRIPT)
  y4_expected <- sequence_ion_masses_batch("TIDE", ion_type = "y", script_path = TEST_PY_SCRIPT)
  expect_equal(ladder$b_mass[ladder$cleavage_position == 3], b3_expected, tolerance = 1e-6)
  expect_equal(ladder$y_mass[ladder$cleavage_position == 3], y4_expected, tolerance = 1e-6)
})

test_that("an internal PTM only shows up in fragments that include its site", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")
  bare <- proteoform(id = "bare", sequence = "PEPTIDE", provenance = "manual")
  modified <- proteoform(
    id = "modified", sequence = "PEPTIDE",
    ptms = list(ptm(site = 4, mass_delta_mono = 79.966331, name = "Phospho")), # T at position 4
    provenance = "manual"
  )

  bare_ladder <- generate_fragment_ladder(bare, script_path = TEST_PY_SCRIPT)
  mod_ladder <- generate_fragment_ladder(modified, script_path = TEST_PY_SCRIPT)

  # b-ions covering position 4 (cleavage_position >= 4) gain the delta; earlier ones don't
  expect_equal(mod_ladder$b_mass[1:3], bare_ladder$b_mass[1:3], tolerance = 1e-9)
  expect_equal(mod_ladder$b_mass[4:6], bare_ladder$b_mass[4:6] + 79.966331, tolerance = 1e-6)

  # y-ions covering position 4 (cleavage_position < 4, i.e. fragment starts at 4 or earlier
  # within the C-terminal side) gain the delta; later ones (fragment starts after site 4) don't
  expect_equal(mod_ladder$y_mass[1:3], bare_ladder$y_mass[1:3] + 79.966331, tolerance = 1e-6)
  expect_equal(mod_ladder$y_mass[4:6], bare_ladder$y_mass[4:6], tolerance = 1e-9)
})

test_that("an N-term PTM shows up in every b-ion and no y-ion", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")
  bare <- proteoform(id = "bare", sequence = "PEPTIDE", provenance = "manual")
  modified <- proteoform(
    id = "modified", sequence = "PEPTIDE",
    ptms = list(ptm(site = "N-term", mass_delta_mono = 42.0106, name = "Acetyl")),
    provenance = "manual"
  )

  bare_ladder <- generate_fragment_ladder(bare, script_path = TEST_PY_SCRIPT)
  mod_ladder <- generate_fragment_ladder(modified, script_path = TEST_PY_SCRIPT)

  expect_equal(mod_ladder$b_mass, bare_ladder$b_mass + 42.0106, tolerance = 1e-6)
  expect_equal(mod_ladder$y_mass, bare_ladder$y_mass, tolerance = 1e-9)
})

test_that("a C-term PTM shows up in every y-ion and no b-ion", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")
  bare <- proteoform(id = "bare", sequence = "PEPTIDE", provenance = "manual")
  modified <- proteoform(
    id = "modified", sequence = "PEPTIDE",
    ptms = list(ptm(site = "C-term", mass_delta_mono = -0.984016, name = "Amidated")),
    provenance = "manual"
  )

  bare_ladder <- generate_fragment_ladder(bare, script_path = TEST_PY_SCRIPT)
  mod_ladder <- generate_fragment_ladder(modified, script_path = TEST_PY_SCRIPT)

  expect_equal(mod_ladder$b_mass, bare_ladder$b_mass, tolerance = 1e-9)
  expect_equal(mod_ladder$y_mass, bare_ladder$y_mass - 0.984016, tolerance = 1e-6)
})

test_that("generate_fragment_ladder rejects sequences shorter than 2 residues", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")
  p <- proteoform(id = "short", sequence = "P", provenance = "manual")
  expect_error(generate_fragment_ladder(p, script_path = TEST_PY_SCRIPT), "at least 2 residues")
})

test_that("generate_fragment_ladder requires a proteoform object", {
  expect_error(generate_fragment_ladder(list(sequence = "PEPTIDE")), "requires a proteoform")
})
