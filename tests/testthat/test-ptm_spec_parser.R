test_that("parse_ptm_spec_text parses a single valid PTM", {
  seq <- "MDKFWWHAAWGLCLVPLSLAQIDLNITCRFAGVFHVEKNGRYSISRTEAADLCKAFNSTLPT"
  # residue 62 is T (Thr)
  expect_equal(substr(seq, 62, 62), "T")
  res <- parse_ptm_spec_text(seq, "62_Thr_Phospho", "test")
  expect_length(res$warnings, 0)
  expect_length(res$groups, 1)
  expect_length(res$groups[[1]], 1)
  expect_equal(res$groups[[1]][[1]]$site, 62L)
  expect_equal(res$groups[[1]][[1]]$name, "Phospho")
})

test_that("semicolons separate proteoforms, commas co-occur on one", {
  seq <- "MDKFWWHAAWGLCLVPLSLAQIDLNITCRFAGVFHVEKNGRYSISRTEAADLCKAFNSTLPT"
  res <- parse_ptm_spec_text(seq, "62_Thr_Phospho; 1_Met_Acetyl,7_His_Methyl", "test")
  expect_length(res$groups, 2)
  expect_length(res$groups[[1]], 1)
  expect_length(res$groups[[2]], 2)
})

test_that("residue/AA mismatch is rejected with a clear warning, group dropped", {
  seq <- "MDKFWWHAAWGLCLVPLSLAQIDLNITCRFAGVFHVEKNGRYSISRTEAADLCKAFNSTLPT"
  res <- parse_ptm_spec_text(seq, "62_Leu_Phospho", "test")
  expect_length(res$groups, 0)
  expect_length(res$warnings, 1)
  expect_match(res$warnings[1], "residue 62 is T.*not Leu")
})

test_that("unknown PTM name is rejected", {
  seq <- "MDKFWWHAAWGLCLVPLSLAQIDLNITCRFAGVFHVEKNGRYSISRTEAADLCKAFNSTLPT"
  res <- parse_ptm_spec_text(seq, "62_Thr_NotARealPTM", "test")
  expect_length(res$groups, 0)
  expect_match(res$warnings[1], 'unknown PTM')
})

test_that("out-of-range residue is rejected", {
  seq <- "MDKFWWHAAWGLCLVPLSLAQIDLNITCRFAGVFHVEKNGRYSISRTEAADLCKAFNSTLPT"
  res <- parse_ptm_spec_text(seq, "9999_Thr_Phospho", "test")
  expect_length(res$groups, 0)
  expect_match(res$warnings[1], "out of range")
})

test_that("more than the max PTMs per proteoform is rejected", {
  seq <- "MDKFWWHAAWGLCLVPLSLAQIDLNITCRFAGVFHVEKNGRYSISRTEAADLCKAFNSTLPT"
  spec <- paste(sprintf("%d_%s_Phospho", 1:6,
    c("Met", "Asp", "Lys", "Phe", "Trp", "Trp")), collapse = ",")
  res <- parse_ptm_spec_text(seq, spec, "test")
  expect_length(res$groups, 0)
  expect_match(res$warnings[1], "exceeds the limit")
})

test_that("empty spec text produces no groups and no warnings", {
  seq <- "MDKFWWHAAWGLCLVPLSLAQIDLNITCRFAGVFHVEKNGRYSISRTEAADLCKAFNSTLPT"
  res <- parse_ptm_spec_text(seq, "", "test")
  expect_length(res$groups, 0)
  expect_length(res$warnings, 0)
})

test_that("parsed PTMs apply the correct mass shift via the real proteoform/ladder pipeline", {
  skip_if_not(mass_engine_available, "Python mass engine not available")
  seq <- "MDKFWWHAAWGLCLVPLSLAQIDLNITCRFAGVFHVEKNGRYSISRTEAADLCKAFNSTLPTMAQMEKALSIGFETCRYGFIEGHVVIPRIHPNSICAANNTGVYILTSNTSQYDTYCFNASAPPEEDCTSVTDLPNAFDGPITITIVNRDGTRYVQKGEYRTNPEDIYPSNPTDDDVSSGSSSERSSTSGGYIFYTFSTVHPIPDEDSPWITDSTDRIPAT"
  res <- parse_ptm_spec_text(seq, "133_Thr_Phospho", "test")
  expect_length(res$groups, 1)

  pf_bare <- proteoform(id = "bare", sequence = seq, provenance = "manual")
  pf_mod <- proteoform(id = "mod", sequence = seq, ptms = res$groups[[1]], provenance = "manual")

  mass_delta <- proteoform_mass(pf_mod)$mass - proteoform_mass(pf_bare)$mass
  expect_equal(mass_delta, 79.966331, tolerance = 1e-4)

  ladder_bare <- generate_fragment_ladder(pf_bare)
  ladder_mod <- generate_fragment_ladder(pf_mod)
  expect_equal(ladder_mod$b_mass[132] - ladder_bare$b_mass[132], 0, tolerance = 1e-6)
  expect_equal(ladder_mod$b_mass[133] - ladder_bare$b_mass[133], 79.966331, tolerance = 1e-4)
})
