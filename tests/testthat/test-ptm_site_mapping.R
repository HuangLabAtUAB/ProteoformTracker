test_that("align_sequences gives an identity map for identical sequences", {
  map <- align_sequences("AGLVE", "AGLVE")
  expect_equal(map, 1:5)
})

test_that("align_sequences maps an insertion to NA and shifts downstream positions", {
  map <- align_sequences("AGLKVE", "AGLVE")
  expect_equal(map[1:3], 1:3)
  expect_true(is.na(map[4]))
  expect_equal(map[5:6], 4:5)
})

test_that("map_ptm_site maps terminal sites trivially", {
  mod <- ptm(site = "N-term", mass_delta_mono = 42.0106, name = "Acetyl")
  result <- map_ptm_site(mod, "AGCK", "MAGCK")
  expect_equal(result$site, "N-term")
  expect_true(result$applicable)
})

test_that("map_ptm_site maps an integer site across an upstream insertion", {
  target_seq <- "AAAAAGGGGGSLLLLLVVVVV"
  isoform_seq <- "AAAAAKKKGGGGGSLLLLLVVVVV"
  mod <- ptm(site = 11, mass_delta_mono = 79.966331, name = "Phospho")

  result <- map_ptm_site(mod, target_seq, isoform_seq)
  expect_true(result$applicable)
  expect_equal(result$site, 14)
})

test_that("map_ptm_site flags a site removed by a deletion as not applicable", {
  target_seq <- "AAAAAGGGGGSLLLLLVVVVV"
  isoform_seq <- "AAAAAGGGGGVVVVV" # the S/L block is spliced out
  mod <- ptm(site = 11, mass_delta_mono = 79.966331, name = "Phospho")

  result <- map_ptm_site(mod, target_seq, isoform_seq)
  expect_false(result$applicable)
  expect_match(result$reason, "splicing")
})

test_that("map_ptm_site flags a residue substitution at the mapped position", {
  target_seq <- "AAAAAGGGGGSLLLLLVVVVV"
  isoform_seq <- "AAAAAGGGGGQLLLLLVVVVV" # S -> Q right at the site
  mod <- ptm(site = 11, mass_delta_mono = 79.966331, name = "Phospho")

  result <- map_ptm_site(mod, target_seq, isoform_seq)
  expect_false(result$applicable)
  expect_match(result$reason, "mismatch")
})

test_that("propagate_ptms_to_relevant_set applies, shifts, and skips PTMs correctly across a mixed isoform set", {
  target_seq <- "AAAAAGGGGGSLLLLLVVVVV"
  target <- proteoform(
    id = "target", sequence = target_seq,
    ptms = list(ptm(site = 11, mass_delta_mono = 79.966331, name = "Phospho")),
    provenance = "manual"
  )

  iso_shifted <- proteoform(id = "iso_shifted", sequence = "AAAAAKKKGGGGGSLLLLLVVVVV", provenance = "manual")
  iso_deleted <- proteoform(id = "iso_deleted", sequence = "AAAAAGGGGGVVVVV", provenance = "manual")
  iso_identical <- proteoform(id = "iso_identical", sequence = target_seq, provenance = "manual")

  result <- propagate_ptms_to_relevant_set(target, list(iso_shifted, iso_deleted, iso_identical))

  by_id <- setNames(result$proteoforms, vapply(result$proteoforms, function(p) p$id, character(1)))

  expect_equal(length(by_id$iso_shifted$ptms), 1)
  expect_equal(by_id$iso_shifted$ptms[[1]]$site, 14)

  expect_equal(length(by_id$iso_deleted$ptms), 0)

  expect_equal(length(by_id$iso_identical$ptms), 1)
  expect_equal(by_id$iso_identical$ptms[[1]]$site, 11)

  expect_equal(nrow(result$report), 3)
  expect_true(result$report$applied[result$report$isoform_id == "iso_shifted"])
  expect_false(result$report$applied[result$report$isoform_id == "iso_deleted"])
  expect_true(result$report$applied[result$report$isoform_id == "iso_identical"])
})

test_that("propagate_ptms_to_relevant_set preserves an isoform's own pre-existing PTMs", {
  target_seq <- "AAAAAGGGGGSLLLLLVVVVV"
  target <- proteoform(
    id = "target", sequence = target_seq,
    ptms = list(ptm(site = 11, mass_delta_mono = 79.966331, name = "Phospho")),
    provenance = "manual"
  )
  iso_with_own_ptm <- proteoform(
    id = "iso_own_ptm", sequence = target_seq,
    ptms = list(ptm(site = "N-term", mass_delta_mono = 42.0106, name = "Acetyl")),
    provenance = "manual"
  )

  result <- propagate_ptms_to_relevant_set(target, list(iso_with_own_ptm))
  expect_equal(length(result$proteoforms[[1]]$ptms), 2)
})

test_that("propagate_ptms_to_relevant_set errors when there are no PTMs to propagate", {
  bare_target <- proteoform(id = "bare", sequence = "AAAAA", provenance = "manual")
  isoform <- proteoform(id = "iso", sequence = "AAAAA", provenance = "manual")
  expect_error(propagate_ptms_to_relevant_set(bare_target, list(isoform)), "no ptms to propagate")
})
