test_that("proteoform() constructs a valid object", {
  p <- proteoform(id = "test1", sequence = "MAGCK", provenance = "manual")
  expect_s3_class(p, "proteoform")
  expect_equal(p$sequence, "MAGCK")
  expect_equal(length(p$ptms), 0)
})

test_that("proteoform() rejects non-standard residues", {
  expect_error(
    proteoform(id = "bad", sequence = "MAGCKZ1", provenance = "manual"),
    "non-standard amino acid"
  )
})

test_that("proteoform() rejects invalid provenance", {
  expect_error(
    proteoform(id = "bad", sequence = "MAGCK", provenance = "not_a_module"),
    "provenance"
  )
})

test_that("proteoform() rejects empty id or sequence", {
  expect_error(proteoform(id = "", sequence = "MAGCK", provenance = "manual"))
  expect_error(proteoform(id = "x", sequence = "", provenance = "manual"))
})

test_that("ptm() validates site and carries mass deltas", {
  m <- ptm(site = 5, mass_delta_mono = 79.9663, name = "Phospho")
  expect_s3_class(m, "ptm")
  expect_equal(m$mass_delta_avg, 79.9663)

  m2 <- ptm(site = "N-term", mass_delta_mono = 42.0106, name = "Acetyl")
  expect_equal(m2$site, "N-term")

  expect_error(ptm(site = -1, mass_delta_mono = 10))
  expect_error(ptm(site = 1.5, mass_delta_mono = 10))
})

test_that("proteoform() carries a list of ptms", {
  mods <- list(
    ptm(site = 5, mass_delta_mono = 79.9663, name = "Phospho"),
    ptm(site = "N-term", mass_delta_mono = 42.0106, name = "Acetyl")
  )
  p <- proteoform(id = "test2", sequence = "MAGCK", ptms = mods, provenance = "manual")
  expect_equal(length(p$ptms), 2)
})

test_that("predict_nterminal_met_excision applies the small-penultimate-residue rule", {
  # Ala at position 2 -> small side chain -> Met excised
  r1 <- predict_nterminal_met_excision("MAGCK")
  expect_true(r1$met_excised)
  expect_equal(r1$mature_sequence, "AGCK")

  # Trp at position 2 -> bulky side chain -> Met retained
  r2 <- predict_nterminal_met_excision("MWGCK")
  expect_false(r2$met_excised)
  expect_equal(r2$mature_sequence, "MWGCK")

  # No leading Met -> unchanged
  r3 <- predict_nterminal_met_excision("AGCK")
  expect_false(r3$met_excised)
})

test_that("build_proteoform applies Met excision and records provenance/notes", {
  p <- build_proteoform(id = "iso1", raw_sequence = "MAGCK", provenance = "module2_longread_orf")
  expect_equal(p$sequence, "AGCK")
  expect_equal(p$provenance, "module2_longread_orf")
  expect_true(any(grepl("Met excision", p$metadata$processing_notes)))
})

test_that("build_proteoform applies a user-supplied mature_start (signal peptide cleavage)", {
  # raw: signal peptide "MKLV" (positions 1-4), mature chain starts at position 5 = "SAGCK"
  p <- build_proteoform(
    id = "iso2", raw_sequence = "MKLVSAGCK", mature_start = 5,
    provenance = "manual"
  )
  expect_equal(p$sequence, "SAGCK")
  expect_true(any(grepl("signal peptide/propeptide removed", p$metadata$processing_notes)))
})
