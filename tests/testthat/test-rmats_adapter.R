# RI/A5SS/A3SS parsing and matching -- SE/MXE (and build_rmats_arm_isoform()/
# build_rmats_full_alignment(), which need a live Ensembl REST fetch per
# candidate transcript) have never had automated tests either, so these
# focus on what's fully testable without network access: the parsers
# (pure, deterministic) and the matching logic against the real, local
# reference_exon_index (gated on exon_index_available, same pattern as the
# RF/GLM model tests -- see helper.R). tests/RI_test.txt, tests/A5SS_test.txt,
# and tests/A3SS_test.txt are real rMATS output rows, independently
# verified (via IsoPepTracker) to represent solid isoform comparisons.

test_that("parse_rmats_ri reads real RI output with correct 0-based/1-based coordinate conversion", {
  events <- parse_rmats_ri("../RI_test.txt")
  expect_equal(nrow(events), 1)
  e <- events[1, ]
  expect_equal(e$gene_id, "ENSG00000114416")
  expect_equal(e$gene_symbol, "FXR1")
  expect_equal(e$chr, "3")
  expect_equal(e$strand, "+")
  # riExonStart_0base=180962882 (0-based) -> 180962883 (1-based); riExonEnd
  # already 1-based inclusive, unchanged
  expect_equal(e$ri_start, 180962883)
  expect_equal(e$ri_end, 180963090)
  expect_equal(e$flank_lo_start, 180962883)
  expect_equal(e$flank_lo_end, 180962940)
  expect_equal(e$flank_hi_start, 180963028)
  expect_equal(e$flank_hi_end, 180963090)
  # the retained-intron exon is exactly the union of the two flanks
  expect_equal(e$ri_start, e$flank_lo_start)
  expect_equal(e$ri_end, e$flank_hi_end)
  expect_equal(e$inc_level_difference, -0.115)
})

test_that("parse_rmats_ri errors clearly on a non-RI file", {
  expect_error(parse_rmats_ri("../A5SS_test.txt"), "Not a valid rMATS RI file")
})

test_that("parse_rmats_a5ss reads real A5SS output with correct coordinates", {
  events <- parse_rmats_a5ss("../A5SS_test.txt")
  expect_equal(nrow(events), 1)
  e <- events[1, ]
  expect_equal(e$gene_symbol, "PARP2")
  expect_equal(e$strand, "+")
  expect_equal(e$long_start, 20344932)
  expect_equal(e$long_end, 20345126)
  expect_equal(e$short_start, 20344932)
  expect_equal(e$short_end, 20345087)
  expect_equal(e$flank_start, 20345394)
  expect_equal(e$flank_end, 20345464)
  # long and short forms share their start and differ only at the end (the
  # alternative 5' splice donor site) -- and the long form is strictly
  # longer, extending further toward the flank
  expect_equal(e$long_start, e$short_start)
  expect_gt(e$long_end, e$short_end)
  # on this + strand gene, the flank sits beyond BOTH forms' shared end
  expect_gt(e$flank_start, e$long_end)
})

test_that("parse_rmats_a3ss reads real A3SS output with correct coordinates", {
  events <- parse_rmats_a3ss("../A3SS_test.txt")
  expect_equal(nrow(events), 1)
  e <- events[1, ]
  expect_equal(e$gene_symbol, "PDLIM5")
  expect_equal(e$strand, "+")
  expect_equal(e$long_start, 94575616)
  expect_equal(e$long_end, 94576034)
  expect_equal(e$short_start, 94575943)
  expect_equal(e$short_end, 94576034)
  expect_equal(e$flank_start, 94573351)
  expect_equal(e$flank_end, 94573393)
  # long and short forms share their end and differ only at the start (the
  # alternative 3' splice acceptor site)
  expect_equal(e$long_end, e$short_end)
  expect_lt(e$long_start, e$short_start)
  # on this + strand gene, the flank sits before BOTH forms' shared start --
  # the opposite side from A5SS's flank, confirming the two event types
  # place their single reported flank on opposite sides of the alt exon
  expect_lt(e$flank_end, e$long_start)
})

test_that("parse_rmats_a3ss errors clearly on a non-A3SS-shaped file", {
  expect_error(parse_rmats_a3ss("../RI_test.txt"), "Not a valid rMATS A3SS file")
})

test_that("match_rmats_ri_retained_transcripts finds only exact whole-exon matches, no adjacency required", {
  gene_tbl <- data.frame(
    transcript_id = c("T1", "T1", "T2", "T2", "T2"),
    exon_number = c(1, 2, 1, 2, 3),
    start = c(100, 500, 100, 500, 900),
    end = c(200, 700, 200, 700, 950),
    stringsAsFactors = FALSE
  )
  # T1's own 2nd exon exactly matches a "retained" span of 500-700
  hits <- match_rmats_ri_retained_transcripts(gene_tbl, list(start = 500, end = 700))
  expect_setequal(hits$transcript_id, c("T1", "T2"))
  expect_true(all(hits$anchor == "both"))
  # no exon anywhere matches this span
  none <- match_rmats_ri_retained_transcripts(gene_tbl, list(start = 1, end = 2))
  expect_equal(nrow(none), 0)
})

test_that("match_rmats_arm_transcripts degrades a NULL-side flank to a one-flank-only requirement", {
  # a cassette exon (300-400) adjacent to a real flank (100-200) on its
  # five_prime side, with NOTHING known/required on the three_prime side
  # (the A5SS/A3SS case: only one flank is ever reported)
  gene_tbl <- data.frame(
    transcript_id = c("MATCHES", "MATCHES", "NO_MATCH", "NO_MATCH"),
    exon_number = c(1, 2, 1, 2),
    start = c(100, 300, 999, 300),
    end = c(200, 400, 1099, 400),
    stringsAsFactors = FALSE
  )
  flanks <- list(five_prime = list(start = 100, end = 200), three_prime = NULL)
  cassette <- list(start = 300, end = 400)
  result <- match_rmats_arm_transcripts(gene_tbl, flanks, cassette)
  expect_equal(result$transcript_id, "MATCHES")
  expect_equal(result$anchor, "five_prime")
})

test_that("match_rmats_ri_transcripts finds real candidates for FXR1's retained and spliced forms", {
  skip_if_not(exon_index_available, "reference_exon_index not available in this environment")
  events <- parse_rmats_ri("../RI_test.txt")
  m <- match_rmats_ri_transcripts(events[1, ], reference_exon_index)

  expect_setequal(names(m$arms), c("retained", "spliced"))
  expect_gt(nrow(m$arms$retained$candidates), 0)
  expect_gt(nrow(m$arms$spliced$candidates), 0)
  # retained-arm matches are always "both" (a full exact-coordinate match,
  # not a partial one-flank match)
  expect_true(all(m$arms$retained$candidates$anchor == "both"))
  # the canonical transcript splices the intron out normally, not a
  # retained-intron transcript, for this real event
  expect_true(m$canonical_transcript_id %in% m$arms$spliced$candidates$transcript_id)
  # per-arm flanks: retained needs none, spliced needs the real pair
  expect_null(m$flanks$retained$five_prime)
  expect_null(m$flanks$retained$three_prime)
  expect_false(is.null(m$flanks$spliced$five_prime))
  expect_false(is.null(m$flanks$spliced$three_prime))
  expect_length(m$highlight_regions, 3)
})

test_that("match_rmats_a5ss_transcripts finds real candidates anchored on the three_prime flank only", {
  skip_if_not(exon_index_available, "reference_exon_index not available in this environment")
  events <- parse_rmats_a5ss("../A5SS_test.txt")
  m <- match_rmats_a5ss_transcripts(events[1, ], reference_exon_index)

  expect_setequal(names(m$arms), c("long", "short"))
  expect_gt(nrow(m$arms$long$candidates), 0)
  expect_gt(nrow(m$arms$short$candidates), 0)
  # only the three_prime side was ever reported, so every real match must
  # be anchored there specifically, never five_prime or both
  expect_true(all(m$arms$long$candidates$anchor == "three_prime"))
  expect_true(all(m$arms$short$candidates$anchor == "three_prime"))
  expect_true(m$canonical_transcript_id %in% m$arms$short$candidates$transcript_id)
})

test_that("match_rmats_a3ss_transcripts finds real candidates anchored on the five_prime flank only", {
  skip_if_not(exon_index_available, "reference_exon_index not available in this environment")
  events <- parse_rmats_a3ss("../A3SS_test.txt")
  m <- match_rmats_a3ss_transcripts(events[1, ], reference_exon_index)

  expect_setequal(names(m$arms), c("long", "short"))
  expect_gt(nrow(m$arms$long$candidates), 0)
  expect_gt(nrow(m$arms$short$candidates), 0)
  expect_true(all(m$arms$long$candidates$anchor == "five_prime"))
  expect_true(all(m$arms$short$candidates$anchor == "five_prime"))
  expect_true(m$canonical_transcript_id %in% m$arms$long$candidates$transcript_id)
})

test_that("RMATS_PARSERS/RMATS_MATCHERS dispatch tables cover every RMATS_EVENT_TYPES entry", {
  expect_setequal(names(RMATS_PARSERS), RMATS_EVENT_TYPES)
  expect_setequal(names(RMATS_MATCHERS), RMATS_EVENT_TYPES)
  expect_true(all(vapply(RMATS_PARSERS, is.function, logical(1))))
  expect_true(all(vapply(RMATS_MATCHERS, is.function, logical(1))))
})

test_that("rmats_event_region_text handles all five event types", {
  ri <- parse_rmats_ri("../RI_test.txt")[1, ]
  expect_equal(rmats_event_region_text(ri, "RI"), "3:180962883-180963090")

  a5 <- parse_rmats_a5ss("../A5SS_test.txt")[1, ]
  expect_equal(rmats_event_region_text(a5, "A5SS"), "14:20344932-20345126 (long) / 14:20344932-20345087 (short)")

  a3 <- parse_rmats_a3ss("../A3SS_test.txt")[1, ]
  expect_equal(rmats_event_region_text(a3, "A3SS"), "4:94575616-94576034 (long) / 4:94575943-94576034 (short)")
})
