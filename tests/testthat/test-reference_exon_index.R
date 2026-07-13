# Fixture: real CDS rows for CD44 canonical (ENST00000428726, 742 aa,
# 18 exons) and isoform 11 (ENST00000434472, 429 aa, 10 exons), trimmed
# from a full Ensembl GTF. Expected residue boundaries were independently
# cross-validated against direct sequence comparison (UniProt P16070-1 vs
# P16070-11) before this test was written -- see the fragment_collision
# divergence-point fix in commit history.

test_that("build_reference_exon_index reproduces validated CD44 residue boundaries (canonical)", {
  skip_if_not(requireNamespace("rtracklayer", quietly = TRUE), "rtracklayer not available")
  index <- build_reference_exon_index("fixtures/cd44_cds.gtf", prefiltered = TRUE)

  canonical <- get_exon_structure(index, "CD44", canonical_only = TRUE)
  expect_equal(unique(canonical$transcript_id), "ENST00000428726")
  expect_equal(unique(canonical$protein_length), 742)
  expect_equal(nrow(canonical), 18)

  # exon 5 ends exactly where the sequence-level divergence point was found
  expect_equal(canonical$residue_end[canonical$exon_number == 5], 222)
  # exon 6 (start of the isoform-skipped variable-exon block)
  expect_equal(canonical$residue_start[canonical$exon_number == 6], 223)
  expect_equal(canonical$residue_end[canonical$exon_number == 18], 742)
})

test_that("build_reference_exon_index reproduces validated CD44 residue boundaries (isoform 11)", {
  skip_if_not(requireNamespace("rtracklayer", quietly = TRUE), "rtracklayer not available")
  index <- build_reference_exon_index("fixtures/cd44_cds.gtf", prefiltered = TRUE)

  iso11 <- get_exon_structure(index, transcript_id = "ENST00000434472")
  expect_equal(unique(iso11$protein_length), 429)
  expect_equal(nrow(iso11), 10)

  # exons 1-5 identical to canonical (same genomic coordinates -> same residues)
  expect_equal(iso11$residue_end[iso11$exon_number == 5], 222)
  # exon 6 here is canonical's exon 14 (skips the variable-exon block)
  expect_equal(iso11$residue_start[iso11$exon_number == 6], 223)
  expect_equal(iso11$residue_end[iso11$exon_number == 10], 429)
  expect_false(unique(iso11$is_canonical))
})

test_that("get_exon_structure errors without gene_symbol or transcript_id", {
  skip_if_not(requireNamespace("rtracklayer", quietly = TRUE), "rtracklayer not available")
  index <- build_reference_exon_index("fixtures/cd44_cds.gtf", prefiltered = TRUE)
  expect_error(get_exon_structure(index), "requires gene_symbol or transcript_id")
})

test_that("get_exon_structure returns all transcripts for a gene symbol when not canonical_only", {
  skip_if_not(requireNamespace("rtracklayer", quietly = TRUE), "rtracklayer not available")
  index <- build_reference_exon_index("fixtures/cd44_cds.gtf", prefiltered = TRUE)
  hits <- get_exon_structure(index, "CD44")
  expect_equal(sort(unique(hits$transcript_id)), sort(c("ENST00000428726", "ENST00000434472")))
})
