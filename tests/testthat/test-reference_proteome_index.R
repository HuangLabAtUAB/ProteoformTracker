test_that("read_fasta parses ids and concatenated sequences", {
  seqs <- read_fasta("fixtures/test_proteome.fasta")
  expect_equal(names(seqs), c("PROT1", "PROT2", "PROT3", "PROT4"))
  expect_equal(unname(seqs["PROT1"]), "MAGCKWERTYAGCKWERTY")
})

test_that("parse_uniprot_accession extracts the accession from sp|/tr| headers", {
  expect_equal(parse_uniprot_accession("sp|P12345|FOO_HUMAN"), "P12345")
  expect_equal(parse_uniprot_accession("tr|Q9XYZ1|BAR_HUMAN"), "Q9XYZ1")
  expect_equal(parse_uniprot_accession("PROT1"), "PROT1")
})

test_that("build_reference_mass_index skips non-standard residues and sorts by mass", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")
  out_path <- tempfile(fileext = ".rds")
  index <- build_reference_mass_index(
    "fixtures/test_proteome.fasta", out_path,
    script_path = "../../python/ptracker_mass.py"
  )

  expect_true(file.exists(out_path))
  expect_false("PROT4" %in% index$id)
  expect_equal(index$mass, sort(index$mass))
  expect_equal(nrow(index), 3)
})

test_that("query_confounding_proteins returns entries within the mass window, excludes self", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")
  out_path <- tempfile(fileext = ".rds")
  index <- build_reference_mass_index(
    "fixtures/test_proteome.fasta", out_path,
    script_path = "../../python/ptracker_mass.py"
  )

  prot1_mass <- index$mass[index$id == "PROT1"]
  hits <- query_confounding_proteins(index, target_mass = prot1_mass, window_da = 5, exclude_id = "PROT1")

  expect_false("PROT1" %in% hits$id)
  expect_true(all(abs(hits$mass - prot1_mass) <= 5))

  # far-larger PROT3 should not appear in a narrow window around PROT1's mass
  expect_false("PROT3" %in% hits$id)
})

test_that("query_confounding_proteins returns zero rows when nothing falls in window", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")
  out_path <- tempfile(fileext = ".rds")
  index <- build_reference_mass_index(
    "fixtures/test_proteome.fasta", out_path,
    script_path = "../../python/ptracker_mass.py"
  )

  hits <- query_confounding_proteins(index, target_mass = -1000, window_da = 1)
  expect_equal(nrow(hits), 0)
})
