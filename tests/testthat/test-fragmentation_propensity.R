test_that("local_basic_density counts K/R fraction in the window around a bond", {
  residues <- strsplit("AAAAAAAAAAA", "")[[1]] # no basic residues nearby
  expect_equal(local_basic_density(residues, 6), 0)

  residues_1k <- strsplit("AAAAAKAAAAA", "")[[1]] # exactly one K in the +/-5 window
  expect_gt(local_basic_density(residues_1k, 6), 0)

  residues_2k <- strsplit("AAAAKAKAAAA", "")[[1]] # two basic residues nearby
  expect_gt(local_basic_density(residues_2k, 6), local_basic_density(residues_1k, 6))
})

test_that("is_near_phosphosite detects proximity to a real phospho PTM and is otherwise a no-op", {
  expect_false(is_near_phosphosite(list(), 10))

  far <- list(ptm(site = 50, mass_delta_mono = 79.966331, name = "Phospho", unimod_id = "UNIMOD:21"))
  expect_false(is_near_phosphosite(far, 10))

  near <- list(ptm(site = 12, mass_delta_mono = 79.966331, name = "Phospho", unimod_id = "UNIMOD:21"))
  expect_true(is_near_phosphosite(near, 10))

  # a non-phospho PTM nearby doesn't count
  other <- list(ptm(site = 10, mass_delta_mono = 42.010565, name = "Acetyl", unimod_id = "UNIMOD:1"))
  expect_false(is_near_phosphosite(other, 10))

  # terminal-site PTMs have no residue position and are safely ignored
  terminal <- list(ptm(site = "N-term", mass_delta_mono = 79.966331, name = "Phospho", unimod_id = "UNIMOD:21"))
  expect_false(is_near_phosphosite(terminal, 10))
})

test_that("fragmentation_propensity requires a fitted GLM model", {
  p <- proteoform(id = "test", sequence = "AAAAADPAAAPAAAA", provenance = "manual")
  expect_error(fragmentation_propensity(p, glm_model = NULL), "unavailable")
})

test_that("fragmentation_propensity returns a genuine probability (0-1) per bond, one row per cleavage position", {
  skip_if_not(glm_model_available, "GLM model not built/available in this environment")
  seq <- "AAAAADPAAAPAAAA"
  p <- proteoform(id = "test", sequence = seq, provenance = "manual")
  result <- fragmentation_propensity(p, mode = "denatured", method = "HCD")

  expect_equal(nrow(result), nchar(seq) - 1)
  expect_setequal(colnames(result), c("cleavage_position", "residue_before", "residue_after", "propensity_score"))
  expect_true(all(result$propensity_score >= 0 & result$propensity_score <= 1))

  bond6 <- result[result$cleavage_position == 6, ]
  expect_equal(bond6$residue_before, "D")
  expect_equal(bond6$residue_after, "P")
})

test_that("fragmentation_propensity ranks a Pro/Asp-adjacent bond above a plain-residue bond, all else equal", {
  skip_if_not(glm_model_available, "GLM model not built/available in this environment")
  # bond 6 is D|P (both chemistry effects); bond 3 is a plain A|A bond, same
  # distance-from-terminus tier (>=6, so both hit the same td_bucket) and no
  # K/R anywhere in the sequence (same charge-density bucket for both)
  seq <- "AAAAADPAAAAAAAAAAAAAA"
  p <- proteoform(id = "test", sequence = seq, provenance = "manual")
  result <- fragmentation_propensity(p, mode = "denatured", method = "HCD")
  bond6 <- result$propensity_score[result$cleavage_position == 6]
  bond3 <- result$propensity_score[result$cleavage_position == 3]
  expect_gt(bond6, bond3)
})

test_that("fragmentation_propensity is deterministic (same input, same output) despite averaging over nuisance covariates", {
  skip_if_not(glm_model_available, "GLM model not built/available in this environment")
  p <- proteoform(id = "test", sequence = "AAAAADPAAAPAAAA", provenance = "manual")
  r1 <- fragmentation_propensity(p, mode = "denatured", method = "HCD")
  r2 <- fragmentation_propensity(p, mode = "denatured", method = "HCD")
  expect_equal(r1$propensity_score, r2$propensity_score)
})

test_that("fragmentation_propensity applies phospho suppression only near an actual phosphosite", {
  skip_if_not(glm_model_available, "GLM model not built/available in this environment")
  seq <- "AAAAADPAAAPAAAAAAAAAAAAAAAAAAAAAA"
  ptms <- list(ptm(site = 8, mass_delta_mono = 79.966331, name = "Phospho", unimod_id = "UNIMOD:21"))
  p_with <- proteoform(id = "test", sequence = seq, ptms = ptms, provenance = "manual")
  p_without <- proteoform(id = "test", sequence = seq, provenance = "manual")
  with_phospho <- fragmentation_propensity(p_with, mode = "denatured", method = "HCD")
  without_phospho <- fragmentation_propensity(p_without, mode = "denatured", method = "HCD")

  near <- with_phospho$propensity_score[with_phospho$cleavage_position == 8]
  near_baseline <- without_phospho$propensity_score[without_phospho$cleavage_position == 8]
  expect_lt(near, near_baseline) # suppressed near the phosphosite

  far <- with_phospho$propensity_score[with_phospho$cleavage_position == 1]
  far_baseline <- without_phospho$propensity_score[without_phospho$cleavage_position == 1]
  expect_equal(far, far_baseline) # unaffected far from the phosphosite
})

test_that("fragmentation_propensity's length effect: shorter proteoforms score higher on average, all else equal", {
  skip_if_not(glm_model_available, "GLM model not built/available in this environment")
  set.seed(1)
  aa <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
  mkseq <- function(n) paste(sample(aa, n, replace = TRUE), collapse = "")
  short_pf <- proteoform(id = "short", sequence = mkseq(60), provenance = "manual")
  long_pf <- proteoform(id = "long", sequence = mkseq(300), provenance = "manual")
  short_scores <- fragmentation_propensity(short_pf, mode = "denatured", method = "HCD")$propensity_score
  long_scores <- fragmentation_propensity(long_pf, mode = "denatured", method = "HCD")$propensity_score
  expect_gt(median(short_scores), median(long_scores))
})

test_that("fragmentation_propensity rejects unimplemented dissociation methods", {
  skip_if_not(glm_model_available, "GLM model not built/available in this environment")
  p <- proteoform(id = "test", sequence = "AAAAADPAAAPAAAA", provenance = "manual")
  expect_error(fragmentation_propensity(p, method = "ETD"), "HCD/CID")
  expect_error(fragmentation_propensity(p, method = "UVPD"), "HCD/CID")
})

test_that("fragmentation_propensity rejects sequences shorter than 2 residues", {
  skip_if_not(glm_model_available, "GLM model not built/available in this environment")
  p <- proteoform(id = "short", sequence = "A", provenance = "manual")
  expect_error(fragmentation_propensity(p), "at least 2 residues")
})

test_that("fragmentation_propensity requires a proteoform object", {
  skip_if_not(glm_model_available, "GLM model not built/available in this environment")
  expect_error(fragmentation_propensity(list(sequence = "AAAA")), "requires a proteoform")
})
