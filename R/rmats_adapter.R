# Option 3 pipeline: rMATS alternative-splicing event -> real annotated
# transcript(s) -> full-length proteoform.
#
# rMATS only reports the differential exon(s) and their immediate flanking
# exons -- not the rest of the transcript, so there is no way to compute a
# full protein mass/sequence from the event alone. But rMATS itself is run
# against a real annotation, so the differential exon's genomic coordinates
# almost always already belong to one or more REAL, already-annotated
# transcripts. This adapter's whole strategy is: don't try to invent a
# transcript from a handful of exons' worth of information -- instead, use
# the precomputed genome-wide exon index (R/reference_exon_index.R, the
# same one Option 1 uses) to find which real transcripts of the gene
# structurally match each "arm" of the event (e.g. exon-inclusion vs
# exon-skipping for SE; 1st-exon vs 2nd-exon for MXE), and hand the user
# those real transcripts directly.
#
# Every event type is represented the same way once matched -- a named list
# of "arms", each an (id, is_canonical, pairing_score) candidate table, plus
# a list of genomic highlight_regions (the exon(s) rMATS itself flagged as
# differential) -- so the server/UI code and the exon-alignment renderer
# don't need to know or care which event type produced them. SE, MXE, RI,
# A5SS, and A3SS are all implemented, using the same
# match_rmats_arm_transcripts()/build_rmats_arms_result() core (RI's
# intron-retained arm is the one exception -- see
# match_rmats_ri_retained_transcripts()'s doc comment for why).
#
# RMATS_PARSERS/RMATS_MATCHERS (bottom of this file) are the dispatch
# tables server.R uses to go from the user's chosen event type to the
# right parse_rmats_*()/match_rmats_*_transcripts() pair.

#' Parse an rMATS SE (skipped-exon) results file (tab-delimited, one row per
#' event). Coordinate convention (confirmed directly against a real event --
#' MYOM1, a minus-strand gene -- by cross-referencing against the
#' precomputed exon index): the *ES columns are 0-based starts (need +1 for
#' our 1-based-inclusive convention), the *EE columns are already 1-based
#' inclusive ends (no adjustment).
#'
#' @param path path to an rMATS SE .txt/.MATS.JC.txt file
#' @return data.frame, one row per event: event_id, gene_id (Ensembl,
#'   version-stripped), gene_symbol, chr (no "chr" prefix, matching
#'   reference_exon_index$seqname), strand, target_start, target_end
#'   (1-based inclusive, the differential/skipped exon), flank_lo_start,
#'   flank_lo_end (the genomically-LOWER-coordinate flanking exon -- rMATS'
#'   own "upstream" columns), flank_hi_start, flank_hi_end (the
#'   genomically-HIGHER-coordinate flanking exon -- rMATS' own "downstream"
#'   columns), inc_level_difference
parse_rmats_se <- function(path) {
  df <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("GeneID", "geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd",
                "upstreamES", "upstreamEE", "downstreamES", "downstreamEE")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(sprintf("Not a valid rMATS SE file -- missing column(s): %s", paste(missing, collapse = ", ")))
  }
  strip_quotes <- function(x) gsub('^"|"$', "", x)
  data.frame(
    event_id = df[[1]],
    gene_id = sub("\\.[0-9]+$", "", strip_quotes(df$GeneID)),
    gene_symbol = strip_quotes(df$geneSymbol),
    chr = sub("^chr", "", df$chr),
    strand = df$strand,
    target_start = as.integer(df$exonStart_0base) + 1L,
    target_end = as.integer(df$exonEnd),
    flank_lo_start = as.integer(df$upstreamES) + 1L,
    flank_lo_end = as.integer(df$upstreamEE),
    flank_hi_start = as.integer(df$downstreamES) + 1L,
    flank_hi_end = as.integer(df$downstreamEE),
    inc_level_difference = suppressWarnings(as.numeric(df$IncLevelDifference)),
    stringsAsFactors = FALSE
  )
}

#' Parse an rMATS MXE (mutually-exclusive-exons) results file. Same
#' coordinate convention as parse_rmats_se(): *ES columns 0-based (+1 here),
#' *EE columns already 1-based inclusive.
#'
#' @param path path to an rMATS MXE .txt/.MATS.JC.txt file
#' @return data.frame, one row per event: event_id, gene_id, gene_symbol,
#'   chr, strand, exon1_start, exon1_end (rMATS' "1st exon"), exon2_start,
#'   exon2_end (rMATS' "2nd exon"), flank_lo_start/end, flank_hi_start/end,
#'   inc_level_difference -- same flank convention as parse_rmats_se()
parse_rmats_mxe <- function(path) {
  df <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("GeneID", "geneSymbol", "chr", "strand",
                "1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd",
                "upstreamES", "upstreamEE", "downstreamES", "downstreamEE")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(sprintf("Not a valid rMATS MXE file -- missing column(s): %s", paste(missing, collapse = ", ")))
  }
  strip_quotes <- function(x) gsub('^"|"$', "", x)
  data.frame(
    event_id = df[[1]],
    gene_id = sub("\\.[0-9]+$", "", strip_quotes(df$GeneID)),
    gene_symbol = strip_quotes(df$geneSymbol),
    chr = sub("^chr", "", df$chr),
    strand = df$strand,
    exon1_start = as.integer(df[["1stExonStart_0base"]]) + 1L,
    exon1_end = as.integer(df[["1stExonEnd"]]),
    exon2_start = as.integer(df[["2ndExonStart_0base"]]) + 1L,
    exon2_end = as.integer(df[["2ndExonEnd"]]),
    flank_lo_start = as.integer(df$upstreamES) + 1L,
    flank_lo_end = as.integer(df$upstreamEE),
    flank_hi_start = as.integer(df$downstreamES) + 1L,
    flank_hi_end = as.integer(df$downstreamEE),
    inc_level_difference = suppressWarnings(as.numeric(df$IncLevelDifference)),
    stringsAsFactors = FALSE
  )
}

#' Parse an rMATS RI (retained-intron) results file. Same coordinate
#' convention as parse_rmats_se(): *ES columns 0-based (+1 here), *EE
#' columns already 1-based inclusive. The "retained" form is a single exon
#' spanning riExonStart-riExonEnd (upstream exon + intron + downstream exon
#' merged into one); the "spliced" form has the intron properly removed,
#' leaving upstreamES-upstreamEE and downstreamES-downstreamEE as two
#' separate, adjacent exons. Confirmed directly against tests/RI_test.txt:
#' the retained exon's own start/end exactly equal the upstream flank's
#' start and the downstream flank's end -- it IS the union of the two.
#'
#' @param path path to an rMATS RI .txt/.MATS.JC.txt file
#' @return data.frame, one row per event: event_id, gene_id, gene_symbol,
#'   chr, strand, ri_start, ri_end (the retained-intron exon, spanning both
#'   flanks), flank_lo_start, flank_lo_end, flank_hi_start, flank_hi_end
#'   (same genomic-lower/higher convention as parse_rmats_se()),
#'   inc_level_difference
parse_rmats_ri <- function(path) {
  df <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("GeneID", "geneSymbol", "chr", "strand", "riExonStart_0base", "riExonEnd",
                "upstreamES", "upstreamEE", "downstreamES", "downstreamEE")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(sprintf("Not a valid rMATS RI file -- missing column(s): %s", paste(missing, collapse = ", ")))
  }
  strip_quotes <- function(x) gsub('^"|"$', "", x)
  data.frame(
    event_id = df[[1]],
    gene_id = sub("\\.[0-9]+$", "", strip_quotes(df$GeneID)),
    gene_symbol = strip_quotes(df$geneSymbol),
    chr = sub("^chr", "", df$chr),
    strand = df$strand,
    ri_start = as.integer(df$riExonStart_0base) + 1L,
    ri_end = as.integer(df$riExonEnd),
    flank_lo_start = as.integer(df$upstreamES) + 1L,
    flank_lo_end = as.integer(df$upstreamEE),
    flank_hi_start = as.integer(df$downstreamES) + 1L,
    flank_hi_end = as.integer(df$downstreamEE),
    inc_level_difference = suppressWarnings(as.numeric(df$IncLevelDifference)),
    stringsAsFactors = FALSE
  )
}

#' Parse an rMATS A5SS (alternative 5' splice site) results file. Same
#' coordinate convention as parse_rmats_se(). The "long" and "short" forms
#' are two alternative boundaries of ONE exon (sharing one edge, differing
#' on the other -- the varying edge is the alternative 5' splice donor site
#' itself); unlike SE/MXE, only ONE flanking exon is reported (not an
#' upstream+downstream pair) -- see match_rmats_altss_transcripts()'s doc
#' comment for which side of the alt exon that flank sits on.
#'
#' @param path path to an rMATS A5SS .txt/.MATS.JC.txt file
#' @return data.frame, one row per event: event_id, gene_id, gene_symbol,
#'   chr, strand, long_start, long_end, short_start, short_end,
#'   flank_start, flank_end, inc_level_difference
parse_rmats_a5ss <- function(path) {
  df <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("GeneID", "geneSymbol", "chr", "strand", "longExonStart_0base", "longExonEnd",
                "shortES", "shortEE", "flankingES", "flankingEE")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(sprintf("Not a valid rMATS A5SS file -- missing column(s): %s", paste(missing, collapse = ", ")))
  }
  strip_quotes <- function(x) gsub('^"|"$', "", x)
  data.frame(
    event_id = df[[1]],
    gene_id = sub("\\.[0-9]+$", "", strip_quotes(df$GeneID)),
    gene_symbol = strip_quotes(df$geneSymbol),
    chr = sub("^chr", "", df$chr),
    strand = df$strand,
    long_start = as.integer(df$longExonStart_0base) + 1L,
    long_end = as.integer(df$longExonEnd),
    short_start = as.integer(df$shortES) + 1L,
    short_end = as.integer(df$shortEE),
    flank_start = as.integer(df$flankingES) + 1L,
    flank_end = as.integer(df$flankingEE),
    inc_level_difference = suppressWarnings(as.numeric(df$IncLevelDifference)),
    stringsAsFactors = FALSE
  )
}

#' Parse an rMATS A3SS (alternative 3' splice site) results file. Same file
#' shape and column names as parse_rmats_a5ss() (rMATS reuses
#' longExonStart_0base/longExonEnd/shortES/shortEE/flankingES/flankingEE
#' for both event types) -- the only difference between A5SS and A3SS is
#' which side of the alternative exon the single flank represents, handled
#' entirely in match_rmats_altss_transcripts(), not here.
#'
#' @param path path to an rMATS A3SS .txt/.MATS.JC.txt file
#' @return same columns as parse_rmats_a5ss()
parse_rmats_a3ss <- function(path) {
  df <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("GeneID", "geneSymbol", "chr", "strand", "longExonStart_0base", "longExonEnd",
                "shortES", "shortEE", "flankingES", "flankingEE")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(sprintf("Not a valid rMATS A3SS file -- missing column(s): %s", paste(missing, collapse = ", ")))
  }
  strip_quotes <- function(x) gsub('^"|"$', "", x)
  data.frame(
    event_id = df[[1]],
    gene_id = sub("\\.[0-9]+$", "", strip_quotes(df$GeneID)),
    gene_symbol = strip_quotes(df$geneSymbol),
    chr = sub("^chr", "", df$chr),
    strand = df$strand,
    long_start = as.integer(df$longExonStart_0base) + 1L,
    long_end = as.integer(df$longExonEnd),
    short_start = as.integer(df$shortES) + 1L,
    short_end = as.integer(df$shortEE),
    flank_start = as.integer(df$flankingES) + 1L,
    flank_end = as.integer(df$flankingEE),
    inc_level_difference = suppressWarnings(as.numeric(df$IncLevelDifference)),
    stringsAsFactors = FALSE
  )
}

#' rMATS' own "upstream"/"downstream" column names are by GENOMIC coordinate
#' (upstream = lower coordinate, downstream = higher), NOT by transcription
#' direction -- confirmed directly: for MYOM1 (minus strand), the real
#' annotated exon at the "downstream" (higher-coordinate) position has a
#' LOWER exon_number (i.e. comes FIRST in the transcript, 5' side) than the
#' target exon, which in turn has a lower exon_number than the "upstream"
#' (lower-coordinate) flank. On a "+" strand gene, transcript order matches
#' rMATS' naming directly; on "-" it's reversed.
#'
#' @param strand "+" or "-"
#' @param flank_lo list(start=,end=) -- the genomically-lower-coordinate flank
#' @param flank_hi list(start=,end=) -- the genomically-higher-coordinate flank
#' @return list(five_prime = list(start=,end=), three_prime = list(start=,end=))
genomic_flanks_in_transcript_order <- function(strand, flank_lo, flank_hi) {
  if (identical(strand, "-")) list(five_prime = flank_hi, three_prime = flank_lo)
  else list(five_prime = flank_lo, three_prime = flank_hi)
}

#' Looks up an event's gene in the exon index, by gene_id first (falling
#' back to gene_symbol if that finds nothing -- e.g. a stale/mismatched
#' Ensembl release).
lookup_rmats_gene_tbl <- function(event, exon_index) {
  if (is.null(exon_index)) return(NULL)
  gene_tbl <- exon_index[!is.na(exon_index$gene_id) & exon_index$gene_id == event$gene_id, ]
  if (nrow(gene_tbl) == 0 && !is.na(event$gene_symbol) && nzchar(event$gene_symbol)) {
    gene_tbl <- exon_index[!is.na(exon_index$gene_name) & exon_index$gene_name == event$gene_symbol, ]
  }
  if (nrow(gene_tbl) == 0) NULL else gene_tbl
}

#' Finds every annotated transcript in gene_tbl whose exon chain (ordered
#' 5'->3' by exon_number) contains `flanks$five_prime` and
#' `flanks$three_prime` as either:
#'   - immediately ADJACENT exons (a "skip" match, when middle_exon is
#'     NULL) -- used for SE's exon-skipping arm, or
#'   - separated by exactly one exon matching `middle_exon`'s coordinates
#'     (a "cassette exon included" match) -- used for SE's exon-inclusion
#'     arm and each of MXE's two mutually-exclusive-exon arms.
#' Checking flank ADJACENCY (not just "this exon happens to be present/
#' absent somewhere") is what distinguishes a transcript that genuinely
#' represents this arm of the event from an unrelated/truncated transcript
#' that simply doesn't extend into this locus.
#'
#' For a cassette-exon arm (middle_exon non-NULL: SE-inclusion, or either of
#' MXE's two arms), a transcript only needs ONE of the two flanks to match
#' AND be immediately adjacent to the cassette exon, not both. rMATS'
#' reported flank/cassette coordinates don't always exactly match a real
#' annotated exon boundary on BOTH sides -- confirmed directly on a real
#' MXE event (IL32): the reported "1st exon" and its upstream flank had no
#' exact annotated match anywhere in the gene, while the "2nd exon" and its
#' downstream flank matched many real transcripts (including canonical)
#' exactly. Requiring both flanks meant this real, useful match was missed
#' entirely, even though the downstream side alone is strong evidence for
#' which transcripts include the 2nd exon. For the no-cassette-exon
#' "skip"/exclusion arm (middle_exon NULL, SE only), both flanks are still
#' required to be immediately adjacent to EACH OTHER -- that adjacency IS
#' the entire evidence for "this transcript skips directly from one flank
#' to the other", so it can't be relaxed to just one flank the same way.
#'
#' A5SS/A3SS events only ever report ONE real flank (not an upstream+
#' downstream pair) -- see match_rmats_altss_transcripts(). Passing NULL
#' for whichever of five_prime/three_prime rMATS didn't report degrades
#' the "at least one flank matches" rule down to "the one real flank must
#' match", exactly what's needed when there's only one to check.
#'
#' @param gene_tbl the gene's own subset of reference_exon_index
#' @param flanks list(five_prime=, three_prime=), each either NULL or
#'   list(start=,end=) -- NULL means "no flank reported on this side, don't
#'   require anything here" (used by A5SS/A3SS, which only report one side)
#' @param middle_exon NULL, or list(start=,end=) for the cassette exon
#'   expected adjacent to at least one of the two flanks
#' @return data.frame(transcript_id, anchor) -- anchor is "both",
#'   "five_prime", or "three_prime" (which flank(s) confirmed the match);
#'   always "both" when middle_exon is NULL. 0 rows if nothing matches.
match_rmats_arm_transcripts <- function(gene_tbl, flanks, middle_exon = NULL) {
  fp <- flanks$five_prime; tp <- flanks$three_prime
  matches_coord <- function(row, coord) {
    !is.null(coord) && !is.na(coord$start) && !is.na(row$start) &&
      row$start == coord$start && row$end == coord$end
  }

  ids <- character(0)
  anchors <- character(0)
  tx_ids <- unique(gene_tbl$transcript_id)
  for (tid in tx_ids) {
    tx <- gene_tbl[gene_tbl$transcript_id == tid, ]
    tx <- tx[order(tx$exon_number), ]
    n <- nrow(tx)
    if (is.null(middle_exon)) {
      if (n < 2) next
      for (i in seq_len(n - 1)) {
        if (matches_coord(tx[i, ], fp) && matches_coord(tx[i + 1, ], tp)) {
          ids <- c(ids, tid); anchors <- c(anchors, "both"); break
        }
      }
    } else {
      if (n < 2) next
      for (i in seq_len(n)) {
        if (is.na(tx$start[i]) || tx$start[i] != middle_exon$start || tx$end[i] != middle_exon$end) next
        five_ok <- i > 1 && matches_coord(tx[i - 1, ], fp)
        three_ok <- i < n && matches_coord(tx[i + 1, ], tp)
        if (five_ok || three_ok) {
          anchor <- if (five_ok && three_ok) "both" else if (five_ok) "five_prime" else "three_prime"
          ids <- c(ids, tid); anchors <- c(anchors, anchor); break
        }
      }
    }
  }
  data.frame(transcript_id = ids, anchor = anchors, stringsAsFactors = FALSE)
}

#' Human-readable version of match_rmats_arm_transcripts()'s "anchor" value
#' -- translated into rMATS' own "upstream"/"downstream" terminology (the
#' columns the user's own file uses) rather than "five_prime"/"three_prime",
#' since upstream/downstream flips relative to 5'/3' on a minus-strand gene
#' (see genomic_flanks_in_transcript_order()'s own doc comment).
#'
#' @param anchor "both", "five_prime", or "three_prime"
#' @param strand "+" or "-"
#' @return "" for "both" (nothing extra worth noting); otherwise e.g.
#'   "upstream-flank match only" / "downstream-flank match only"
rmats_anchor_note <- function(anchor, strand) {
  if (identical(anchor, "both")) return("")
  is_upstream <- if (identical(strand, "-")) identical(anchor, "three_prime") else identical(anchor, "five_prime")
  sprintf("%s-flank match only", if (is_upstream) "upstream" else "downstream")
}

#' Builds the final "arms" result for an arbitrary number of named arms of
#' one rMATS event (2 for SE/MXE; kept general for any future event type).
#' Each candidate's pairing_score is the best exon-set Jaccard similarity
#' against any candidate in a DIFFERENT arm -- close to 1 means "the same
#' underlying transcript, +/- this event". This is what lets the UI suggest
#' a sensible default pair instead of picking among several
#' structurally-valid candidates arbitrarily.
#'
#' @param gene_tbl the gene's own subset of reference_exon_index (NULL if
#'   the gene wasn't found -- every arm then comes back with 0 candidates)
#' @param arm_matches named list, arm key -> data.frame(transcript_id,
#'   anchor) as returned by match_rmats_arm_transcripts()
#' @param arm_labels named character vector, arm key -> human-readable label
#' @return list(arms = list(key = list(key=, label=, candidates=data.frame(
#'   transcript_id, is_canonical, anchor, pairing_score)), ...)) --
#'   candidates is a 0-row data.frame (never NULL) when an arm has no match
build_rmats_arms_result <- function(gene_tbl, arm_matches, arm_labels) {
  empty <- data.frame(transcript_id = character(0), is_canonical = logical(0),
                       anchor = character(0), pairing_score = numeric(0), stringsAsFactors = FALSE)
  arm_keys <- names(arm_matches)
  if (is.null(gene_tbl)) {
    arms <- lapply(arm_keys, function(k) list(key = k, label = arm_labels[[k]], candidates = empty))
    names(arms) <- arm_keys
    return(list(arms = arms))
  }

  exon_key_set <- function(tid) {
    tx <- gene_tbl[gene_tbl$transcript_id == tid, ]
    paste(tx$start, tx$end, sep = "-")
  }
  jaccard <- function(a, b) {
    u <- length(union(a, b))
    if (u == 0) return(0)
    length(intersect(a, b)) / u
  }
  is_canonical_of <- function(tid) any(gene_tbl$transcript_id == tid & gene_tbl$is_canonical)
  arm_ids <- lapply(arm_matches, function(df) if (is.null(df)) character(0) else unique(df$transcript_id))

  arms <- lapply(arm_keys, function(k) {
    df <- arm_matches[[k]]
    ids <- arm_ids[[k]]
    other_ids <- unique(unlist(arm_ids[setdiff(arm_keys, k)]))
    anchor_of <- if (nrow(df) == 0) character(0) else setNames(df$anchor, df$transcript_id)
    candidates <- if (length(ids) == 0) {
      empty
    } else if (length(other_ids) == 0) {
      data.frame(transcript_id = ids, is_canonical = vapply(ids, is_canonical_of, logical(1)),
                 anchor = unname(anchor_of[ids]), pairing_score = 0, stringsAsFactors = FALSE)
    } else {
      other_keys <- lapply(other_ids, exon_key_set)
      scores <- vapply(ids, function(tid) {
        own <- exon_key_set(tid)
        max(vapply(other_keys, function(ok) jaccard(own, ok), numeric(1)))
      }, numeric(1))
      d <- data.frame(transcript_id = ids, is_canonical = vapply(ids, is_canonical_of, logical(1)),
                       anchor = unname(anchor_of[ids]), pairing_score = round(unname(scores), 3), stringsAsFactors = FALSE)
      d[order(-d$pairing_score), ]
    }
    list(key = k, label = arm_labels[[k]], candidates = candidates)
  })
  names(arms) <- arm_keys
  list(arms = arms)
}

#' For one parsed SE event, find every annotated transcript of its gene that
#' structurally represents the exon-INCLUSION form or the exon-EXCLUSION/
#' skipping form -- see match_rmats_arm_transcripts()'s doc comment for the
#' matching rule.
#'
#' @param event one row of parse_rmats_se()'s output
#' @param exon_index the precomputed reference_exon_index
#' @return list(arms = list(inclusion=, exclusion=), highlight_regions =
#'   list(list(start=,end=,label=))) -- see build_rmats_arms_result() and
#'   build_rmats_candidate_alignment() for how these are used
match_rmats_se_transcripts <- function(event, exon_index) {
  gene_tbl <- lookup_rmats_gene_tbl(event, exon_index)
  arm_labels <- c(inclusion = "Exon-inclusion form", exclusion = "Exon-skipping form")
  empty_match <- data.frame(transcript_id = character(0), anchor = character(0), stringsAsFactors = FALSE)
  flanks <- genomic_flanks_in_transcript_order(
    event$strand,
    list(start = event$flank_lo_start, end = event$flank_lo_end),
    list(start = event$flank_hi_start, end = event$flank_hi_end)
  )
  cassette <- list(inclusion = list(start = event$target_start, end = event$target_end), exclusion = NULL)
  highlight_regions <- list(
    list(start = event$flank_lo_start, end = event$flank_lo_end, label = "Upstream flank exon (rMATS)"),
    list(start = event$target_start, end = event$target_end, label = "Differential exon (rMATS)"),
    list(start = event$flank_hi_start, end = event$flank_hi_end, label = "Downstream flank exon (rMATS)")
  )
  # flanks is the SAME for both arms here (unlike RI, where the retained
  # arm needs no flanks at all) -- stored per-arm anyway so
  # build_rmats_full_alignment()/build_rmats_arm_isoform() callers can
  # always index matches$flanks[[arm_key]] regardless of event type.
  if (is.null(gene_tbl)) {
    result <- build_rmats_arms_result(NULL, list(inclusion = empty_match, exclusion = empty_match), arm_labels)
    result$highlight_regions <- highlight_regions
    result$flanks <- list(inclusion = flanks, exclusion = flanks)
    result$cassette <- cassette
    result$gene_transcript_ids <- character(0)
    result$canonical_transcript_id <- NA_character_
    return(result)
  }
  inclusion_matches <- match_rmats_arm_transcripts(gene_tbl, flanks, cassette$inclusion)
  exclusion_matches <- match_rmats_arm_transcripts(gene_tbl, flanks, NULL)
  result <- build_rmats_arms_result(gene_tbl, list(inclusion = inclusion_matches, exclusion = exclusion_matches), arm_labels)
  result$highlight_regions <- highlight_regions
  result$flanks <- list(inclusion = flanks, exclusion = flanks)
  result$cassette <- cassette
  result$gene_transcript_ids <- unique(gene_tbl$transcript_id)
  canon <- unique(gene_tbl$transcript_id[gene_tbl$is_canonical])
  result$canonical_transcript_id <- if (length(canon) > 0) canon[1] else NA_character_
  result
}

#' For one parsed MXE event, find every annotated transcript of its gene
#' that structurally includes the "1st exon" (with the "2nd exon" absent)
#' or the "2nd exon" (with the "1st" absent) -- see
#' match_rmats_arm_transcripts()'s doc comment for the matching rule.
#'
#' @param event one row of parse_rmats_mxe()'s output
#' @param exon_index the precomputed reference_exon_index
#' @return same shape as match_rmats_se_transcripts(), with arms exon1/exon2
#'   and two highlight_regions (one per mutually exclusive exon)
match_rmats_mxe_transcripts <- function(event, exon_index) {
  gene_tbl <- lookup_rmats_gene_tbl(event, exon_index)
  arm_labels <- c(exon1 = "1st exon form", exon2 = "2nd exon form")
  empty_match <- data.frame(transcript_id = character(0), anchor = character(0), stringsAsFactors = FALSE)
  flanks <- genomic_flanks_in_transcript_order(
    event$strand,
    list(start = event$flank_lo_start, end = event$flank_lo_end),
    list(start = event$flank_hi_start, end = event$flank_hi_end)
  )
  cassette <- list(exon1 = list(start = event$exon1_start, end = event$exon1_end),
                    exon2 = list(start = event$exon2_start, end = event$exon2_end))
  highlight_regions <- list(
    list(start = event$flank_lo_start, end = event$flank_lo_end, label = "Upstream flank exon (rMATS)"),
    list(start = event$exon1_start, end = event$exon1_end, label = "1st mutually exclusive exon (rMATS)"),
    list(start = event$exon2_start, end = event$exon2_end, label = "2nd mutually exclusive exon (rMATS)"),
    list(start = event$flank_hi_start, end = event$flank_hi_end, label = "Downstream flank exon (rMATS)")
  )
  if (is.null(gene_tbl)) {
    result <- build_rmats_arms_result(NULL, list(exon1 = empty_match, exon2 = empty_match), arm_labels)
    result$highlight_regions <- highlight_regions
    result$flanks <- list(exon1 = flanks, exon2 = flanks)
    result$cassette <- cassette
    result$gene_transcript_ids <- character(0)
    result$canonical_transcript_id <- NA_character_
    return(result)
  }
  exon1_matches <- match_rmats_arm_transcripts(gene_tbl, flanks, cassette$exon1)
  exon2_matches <- match_rmats_arm_transcripts(gene_tbl, flanks, cassette$exon2)
  result <- build_rmats_arms_result(gene_tbl, list(exon1 = exon1_matches, exon2 = exon2_matches), arm_labels)
  result$highlight_regions <- highlight_regions
  result$flanks <- list(exon1 = flanks, exon2 = flanks)
  result$cassette <- cassette
  result$gene_transcript_ids <- unique(gene_tbl$transcript_id)
  canon <- unique(gene_tbl$transcript_id[gene_tbl$is_canonical])
  result$canonical_transcript_id <- if (length(canon) > 0) canon[1] else NA_character_
  result
}

#' Finds every annotated transcript in gene_tbl with an exon EXACTLY
#' spanning ri_exon's coordinates. Unlike every other arm-matching function
#' in this file, this is NOT an adjacency check: the retained-intron form's
#' whole point is a single continuous exon that already merges the upstream
#' exon, the intron, and the downstream exon into one, so there's nothing
#' left for it to be adjacent to -- a transcript either has an exon at
#' exactly this position or it doesn't.
#'
#' @param gene_tbl the gene's own subset of reference_exon_index
#' @param ri_exon list(start=,end=) the full retained-intron exon span
#' @return data.frame(transcript_id, anchor) -- anchor is always "both",
#'   matching the schema match_rmats_arm_transcripts() returns (rmats_
#'   anchor_note() treats "both" as "nothing partial to flag", which is
#'   right here too: this is a full exact-coordinate match, not a
#'   one-flank-only one)
match_rmats_ri_retained_transcripts <- function(gene_tbl, ri_exon) {
  hit <- !is.na(gene_tbl$start) & gene_tbl$start == ri_exon$start & gene_tbl$end == ri_exon$end
  ids <- unique(gene_tbl$transcript_id[hit])
  data.frame(transcript_id = ids, anchor = rep("both", length(ids)), stringsAsFactors = FALSE)
}

#' For one parsed RI event, find every annotated transcript of its gene
#' that structurally represents the intron-RETAINED form (a single exon
#' exactly spanning ri_start-ri_end, see
#' match_rmats_ri_retained_transcripts()) or the intron-SPLICED form
#' (separate, immediately-adjacent upstream/downstream exons -- the same
#' "both flanks directly adjacent to each other" rule
#' match_rmats_arm_transcripts() uses for SE's exclusion arm).
#'
#' Unlike SE/MXE, the two arms need DIFFERENT flanks for constructing their
#' synthetic isoforms (build_rmats_arm_isoform()): the retained arm needs
#' NONE (its one merged exon already covers the whole locus, so anything
#' the backbone contributes there would just create duplicate/overlapping
#' exons), while the spliced arm needs both upstream and downstream. So
#' (unlike match_rmats_se_transcripts()/match_rmats_mxe_transcripts(),
#' where both arms share one flanks object) `result$flanks` here is
#' genuinely per-arm.
#'
#' @param event one row of parse_rmats_ri()'s output
#' @param exon_index the precomputed reference_exon_index
#' @return same shape as match_rmats_se_transcripts(), with arms
#'   retained/spliced and highlight_regions for the upstream flank, the
#'   full retained-intron exon, and the downstream flank
match_rmats_ri_transcripts <- function(event, exon_index) {
  gene_tbl <- lookup_rmats_gene_tbl(event, exon_index)
  arm_labels <- c(retained = "Intron-retained form", spliced = "Intron-spliced form")
  empty_match <- data.frame(transcript_id = character(0), anchor = character(0), stringsAsFactors = FALSE)
  flanks <- genomic_flanks_in_transcript_order(
    event$strand,
    list(start = event$flank_lo_start, end = event$flank_lo_end),
    list(start = event$flank_hi_start, end = event$flank_hi_end)
  )
  ri_exon <- list(start = event$ri_start, end = event$ri_end)
  no_flanks <- list(five_prime = NULL, three_prime = NULL)
  highlight_regions <- list(
    list(start = event$flank_lo_start, end = event$flank_lo_end, label = "Upstream flank exon (rMATS)"),
    list(start = event$ri_start, end = event$ri_end, label = "Retained-intron exon (rMATS)"),
    list(start = event$flank_hi_start, end = event$flank_hi_end, label = "Downstream flank exon (rMATS)")
  )
  if (is.null(gene_tbl)) {
    result <- build_rmats_arms_result(NULL, list(retained = empty_match, spliced = empty_match), arm_labels)
    result$highlight_regions <- highlight_regions
    result$flanks <- list(retained = no_flanks, spliced = flanks)
    result$cassette <- list(retained = ri_exon, spliced = NULL)
    result$gene_transcript_ids <- character(0)
    result$canonical_transcript_id <- NA_character_
    return(result)
  }
  retained_matches <- match_rmats_ri_retained_transcripts(gene_tbl, ri_exon)
  spliced_matches <- match_rmats_arm_transcripts(gene_tbl, flanks, NULL)
  result <- build_rmats_arms_result(gene_tbl, list(retained = retained_matches, spliced = spliced_matches), arm_labels)
  result$highlight_regions <- highlight_regions
  result$flanks <- list(retained = no_flanks, spliced = flanks)
  result$cassette <- list(retained = ri_exon, spliced = NULL)
  result$gene_transcript_ids <- unique(gene_tbl$transcript_id)
  canon <- unique(gene_tbl$transcript_id[gene_tbl$is_canonical])
  result$canonical_transcript_id <- if (length(canon) > 0) canon[1] else NA_character_
  result
}

#' Shared implementation for match_rmats_a5ss_transcripts()/
#' match_rmats_a3ss_transcripts() -- A5SS and A3SS are structurally
#' identical (a "long" vs "short" form of one exon, adjacent to a SINGLE
#' reported flanking exon) and differ only in which side of the
#' alternative exon that flank sits on. That side is a fixed property of
#' the event type, not the strand: A5SS varies the 5' splice DONOR site,
#' which by definition is the boundary between the alternative exon and
#' the intron immediately FOLLOWING it in transcript direction -- so the
#' flank is always on the alt exon's transcript-3' side. A3SS varies the
#' 3' splice ACCEPTOR site, the boundary with the intron immediately
#' PRECEDING the alt exon -- so the flank is always on the transcript-5'
#' side. Confirmed against real + strand data (tests/A5SS_test.txt,
#' tests/A3SS_test.txt): A5SS's flank sits at a HIGHER genomic coordinate
#' than the alt exon, A3SS's at a LOWER one -- consistent with "transcript-
#' 3'/5' side" on a + strand gene, where transcript order matches genomic
#' order. NOT independently re-verified against a minus-strand example (none
#' was available in the provided test files) -- the minus-strand case relies
#' on the same genomic<->transcript-direction remapping already validated
#' for SE/MXE's upstream/downstream columns (see
#' genomic_flanks_in_transcript_order()), applied here to a fixed
#' five_prime/three_prime side instead of two genomically-named columns.
#' No strand-dependent remapping of WHICH raw coordinate (start vs end) is
#' long/short's "shared" vs "varying" boundary is needed anywhere in this
#' code, since matching always uses both forms' complete (start,end) pairs
#' directly from the file, never decomposing them into shared/varying parts.
#'
#' @param event one row of parse_rmats_a5ss()/parse_rmats_a3ss()'s output
#'   (both share the same long_start/long_end/short_start/short_end/
#'   flank_start/flank_end columns)
#' @param exon_index the precomputed reference_exon_index
#' @param flank_side "five_prime" (A3SS) or "three_prime" (A5SS) -- which
#'   side of the alternative exon the single reported flank represents
#' @param arm_labels named c(long=, short=) human-readable arm labels
#' @return same shape as match_rmats_se_transcripts(), with arms long/short
match_rmats_altss_transcripts <- function(event, exon_index, flank_side, arm_labels) {
  gene_tbl <- lookup_rmats_gene_tbl(event, exon_index)
  empty_match <- data.frame(transcript_id = character(0), anchor = character(0), stringsAsFactors = FALSE)
  flank <- list(start = event$flank_start, end = event$flank_end)
  flanks <- if (identical(flank_side, "five_prime")) {
    list(five_prime = flank, three_prime = NULL)
  } else {
    list(five_prime = NULL, three_prime = flank)
  }
  long_exon <- list(start = event$long_start, end = event$long_end)
  short_exon <- list(start = event$short_start, end = event$short_end)
  highlight_regions <- list(
    list(start = event$flank_start, end = event$flank_end, label = "Flanking exon (rMATS)"),
    list(start = event$long_start, end = event$long_end, label = "Long-exon form (rMATS)"),
    list(start = event$short_start, end = event$short_end, label = "Short-exon form (rMATS)")
  )
  if (is.null(gene_tbl)) {
    result <- build_rmats_arms_result(NULL, list(long = empty_match, short = empty_match), arm_labels)
    result$highlight_regions <- highlight_regions
    result$flanks <- list(long = flanks, short = flanks)
    result$cassette <- list(long = long_exon, short = short_exon)
    result$gene_transcript_ids <- character(0)
    result$canonical_transcript_id <- NA_character_
    return(result)
  }
  long_matches <- match_rmats_arm_transcripts(gene_tbl, flanks, long_exon)
  short_matches <- match_rmats_arm_transcripts(gene_tbl, flanks, short_exon)
  result <- build_rmats_arms_result(gene_tbl, list(long = long_matches, short = short_matches), arm_labels)
  result$highlight_regions <- highlight_regions
  result$flanks <- list(long = flanks, short = flanks)
  result$cassette <- list(long = long_exon, short = short_exon)
  result$gene_transcript_ids <- unique(gene_tbl$transcript_id)
  canon <- unique(gene_tbl$transcript_id[gene_tbl$is_canonical])
  result$canonical_transcript_id <- if (length(canon) > 0) canon[1] else NA_character_
  result
}

#' For one parsed A5SS event, find every annotated transcript of its gene
#' that structurally represents the long-exon form or the short-exon form
#' -- see match_rmats_altss_transcripts()'s doc comment.
#'
#' @param event one row of parse_rmats_a5ss()'s output
#' @param exon_index the precomputed reference_exon_index
#' @return see match_rmats_altss_transcripts()
match_rmats_a5ss_transcripts <- function(event, exon_index) {
  match_rmats_altss_transcripts(event, exon_index, "three_prime",
                                 c(long = "Long-exon form", short = "Short-exon form"))
}

#' For one parsed A3SS event, find every annotated transcript of its gene
#' that structurally represents the long-exon form or the short-exon form
#' -- see match_rmats_altss_transcripts()'s doc comment.
#'
#' @param event one row of parse_rmats_a3ss()'s output
#' @param exon_index the precomputed reference_exon_index
#' @return see match_rmats_altss_transcripts()
match_rmats_a3ss_transcripts <- function(event, exon_index) {
  match_rmats_altss_transcripts(event, exon_index, "five_prime",
                                 c(long = "Long-exon form", short = "Short-exon form"))
}

#' Picks a sensible default "backbone" transcript for constructing an arm's
#' full synthetic isoform (build_rmats_arm_isoform()): prefers a real
#' annotated match for THIS arm itself (already structurally verified by
#' match_rmats_arm_transcripts()), falling back to the best-scoring real
#' match from ANOTHER arm of the SAME event (same gene, most likely sharing
#' the rest of the transcript structure -- exactly the MXE case where one
#' arm has real matches and the other doesn't), falling back to the gene's
#' canonical transcript, falling back to just the first transcript of the
#' gene. The user can always override this via the backbone dropdown in the
#' UI -- this is only ever the INITIAL suggestion.
#'
#' @param matches result of match_rmats_se_transcripts()/match_rmats_mxe_transcripts()
#' @param arm_key which arm to pick a backbone for
#' @return transcript_id string, or NA_character_ if nothing at all is available
default_backbone_for_arm <- function(matches, arm_key) {
  arm <- matches$arms[[arm_key]]
  if (nrow(arm$candidates) > 0) return(arm$candidates$transcript_id[1])
  other_keys <- setdiff(names(matches$arms), arm_key)
  other <- do.call(rbind, lapply(other_keys, function(k) matches$arms[[k]]$candidates))
  if (!is.null(other) && nrow(other) > 0) return(other$transcript_id[which.max(other$pairing_score)])
  if (!is.na(matches$canonical_transcript_id) && nzchar(matches$canonical_transcript_id)) return(matches$canonical_transcript_id)
  if (length(matches$gene_transcript_ids) > 0) return(matches$gene_transcript_ids[1])
  NA_character_
}

#' Constructs a full synthetic exon list for one arm of an rMATS event: the
#' backbone transcript's own exons OUTSIDE the local alternative-splicing
#' region (by genomic position, not exact coordinate matching -- see below
#' for why), spliced together with the event's OWN reported local exons
#' (upstream flank, this arm's cassette exon if any, downstream flank)
#' taken verbatim. This is the "construct the whole isoform" piece: rather
#' than only showing real ENST matches (which may not exist for one or
#' both arms -- confirmed directly, IL32/GPSM1), this ALWAYS produces a
#' labeled "constructed" track per arm using a user-selectable (defaulted)
#' backbone transcript for everything rMATS itself didn't report.
#'
#' Position-thresholding (not exact-coordinate matching against the
#' backbone) is deliberate: a backbone transcript's own exon at this locus
#' can genuinely extend further than rMATS' own reported flank boundary
#' (e.g. a partially-coding exon whose 5'UTR portion the backbone's real
#' annotation includes but rMATS' own flank coordinate does not -- confirmed
#' directly on IL32's own "1st exon"/upstream flank). Any backbone exon that
#' OVERLAPS the local region at all is dropped in favor of rMATS' own
#' reported boundary there, rather than requiring an exact match.
#'
#' @param backbone_tid Ensembl transcript id to use for everything outside
#'   the local AS region
#' @param flanks list(five_prime=, three_prime=) genomic coords
#'   (matches$flanks[[arm_key]]) -- EITHER side may be NULL, meaning rMATS
#'   didn't report a flank there (A5SS/A3SS, which only ever report one
#'   real flank) or none applies (RI's retained arm, whose single merged
#'   exon needs no flank override on either side -- see
#'   match_rmats_ri_transcripts()'s doc comment). A NULL side is simply
#'   left out of the local region entirely: the backbone's own exon
#'   structure there is kept as-is, un-clipped, since there's no rMATS-
#'   reported boundary to override it with.
#' @param cassette_exon NULL (skip/exclusion arm, or RI's spliced arm), or
#'   list(start=,end=) for this arm's own middle exon (SE-inclusion, either
#'   MXE arm, RI's retained arm, or either A5SS/A3SS arm)
#' @param exon_index the precomputed reference_exon_index (CDS-only), used
#'   to look up which sub-ranges of the resulting exon list are coding
#' @return NULL if the backbone transcript's exon structure can't be
#'   fetched, else list(exons = data.frame(start,end) full exon list in
#'   genomic order, cds = data.frame(start,end) coding sub-ranges, seqname,
#'   strand)
build_rmats_arm_isoform <- function(backbone_tid, flanks, cassette_exon, exon_index) {
  bb <- tryCatch(fetch_transcript_exons(backbone_tid), error = function(e) NULL)
  if (is.null(bb) || is.null(bb$exons) || nrow(bb$exons) == 0) return(NULL)
  full <- bb$exons

  fp <- flanks$five_prime; tp <- flanks$three_prime
  local_parts <- Filter(Negate(is.null), list(fp, cassette_exon, tp))
  local_bounds <- unlist(lapply(local_parts, function(x) c(x$start, x$end)))
  region_lo <- min(local_bounds); region_hi <- max(local_bounds)

  before <- full[full$end < region_lo, , drop = FALSE]
  after <- full[full$start > region_hi, , drop = FALSE]
  local_exons <- do.call(rbind, lapply(local_parts, function(x) data.frame(start = x$start, end = x$end)))
  all_exons <- rbind(before, local_exons, after)
  all_exons <- all_exons[order(all_exons$start), ]

  # Coding sub-ranges: the backbone's OWN CDS rows (for the before/after,
  # backbone-derived portions) plus whatever CDS evidence exists ANYWHERE in
  # this gene's annotation overlapping the local (rMATS-reported) exons --
  # a local exon's own reported boundary can include real 5'/3'UTR fused
  # onto a CDS portion (confirmed directly, IL32's own "1st exon"), so its
  # coding sub-range has to come from gene-wide CDS evidence, not assumed to
  # be 100% coding just because it's a differential exon.
  backbone_cds <- unique(exon_index[exon_index$transcript_id == backbone_tid, c("start", "end")])
  gene_id <- exon_index$gene_id[exon_index$transcript_id == backbone_tid][1]
  gene_cds <- if (!is.null(gene_id) && !is.na(gene_id)) {
    exon_index[!is.na(exon_index$gene_id) & exon_index$gene_id == gene_id, c("start", "end")]
  } else {
    backbone_cds
  }
  # Collapsed to ONE union interval per local exon (not one row per
  # overlapping transcript) -- a locus this common to the gene can have
  # dozens of transcripts with their own near-identical CDS row there,
  # which would otherwise all be kept as separate "coding" sub-ranges and
  # massively over-count how many coding segments a single exon has.
  local_cds <- do.call(rbind, lapply(seq_len(nrow(local_exons)), function(i) {
    s <- local_exons$start[i]; e <- local_exons$end[i]
    ov <- gene_cds[gene_cds$start <= e & gene_cds$end >= s, ]
    if (nrow(ov) == 0) return(NULL)
    data.frame(start = min(pmax(ov$start, s)), end = max(pmin(ov$end, e)))
  }))
  cds <- unique(rbind(backbone_cds, local_cds))

  list(exons = all_exons, cds = cds, seqname = bb$seqname, strand = bb$strand)
}

RMATS_GENOME_FASTA <- "reference/genome/GRCh38.fa"

.STANDARD_CODON_TABLE <- local({
  bases <- c("T", "C", "A", "G")
  # First base slowest-varying, third base fastest -- standard codon-table
  # ordering (TTT, TTC, TTA, TTG, TCT, ..., GGG), matching `aas` below.
  codons <- paste0(rep(bases, each = 16), rep(bases, each = 4, times = 4), rep(bases, times = 16))
  aas <- "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
  setNames(strsplit(aas, "")[[1]], codons)
})

#' Reverse-complement a DNA sequence (uppercase in, uppercase out).
#'
#' @param seq character(1) DNA sequence (A/C/G/T/N)
reverse_complement_dna <- function(seq) {
  comp <- chartr("ACGTN", "TGCAN", toupper(seq))
  paste(rev(strsplit(comp, "")[[1]]), collapse = "")
}

#' Translate a nucleotide sequence (assumed already in the correct 5'->3'
#' reading frame) using the standard genetic code, stopping at the first
#' in-frame stop codon if one is present (the stop itself excluded from the
#' returned protein) -- i.e. ordinary CDS translation. If NO in-frame stop
#' is found, the entire sequence is translated: reference_exon_index's own
#' CDS ranges (built from GTF "CDS" feature rows -- confirmed directly
#' against ENST00000008180/IL32, whose CDS rows total exactly 168*3 = 504 nt
#' with no in-frame stop anywhere in them) follow GENCODE/Ensembl convention,
#' where the CDS feature EXCLUDES the terminal stop codon -- requiring one
#' to be found would wrongly discard every transcript that supplied exactly
#' its CDS-without-stop, which is the common case. Returns NULL (rather than
#' a protein string containing "X") if an ambiguous ("N"-containing) codon
#' falls within the portion that would actually be kept (before any stop),
#' or the sequence isn't a non-empty multiple of 3 -- both indicate the
#' caller's CDS coordinates are wrong rather than a real biological outcome.
#'
#' @param nt_seq character(1) nucleotide sequence, 5'->3', already
#'   strand-corrected (i.e. reverse-complemented if the source was "-" strand)
#' @return character(1) protein sequence (standard 20 AA letters only), or
#'   NULL if translation didn't produce a clean protein
translate_dna_standard <- function(nt_seq) {
  nt_seq <- toupper(nt_seq)
  n <- nchar(nt_seq)
  if (n < 3 || n %% 3 != 0) return(NULL)
  codons <- substring(nt_seq, seq(1, n, by = 3), seq(3, n, by = 3))
  aa <- .STANDARD_CODON_TABLE[codons]
  stop_at <- which(aa == "*")
  keep_n <- if (length(stop_at) == 0) length(aa) else stop_at[1] - 1L
  if (keep_n == 0) return(NULL)
  kept <- aa[seq_len(keep_n)]
  if (any(is.na(kept))) return(NULL)
  paste(kept, collapse = "")
}

#' Fetch a genomic DNA range (plus/forward strand, uppercase, no header)
#' from the local reference genome FASTA via `samtools faidx` random access
#' (the .fai index makes this a fast seek, not a full-file scan even though
#' the FASTA itself is ~3GB).
#'
#' @param seqname chromosome/contig name (no "chr" prefix, matching
#'   reference_exon_index$seqname)
#' @param start,end 1-based inclusive genomic range
#' @param genome_fasta path to the indexed genome FASTA
#' @return character(1) uppercase DNA sequence, or NULL if the samtools call
#'   failed (missing binary/index, out-of-range coordinates, etc.)
fetch_genome_dna <- function(seqname, start, end, genome_fasta = RMATS_GENOME_FASTA) {
  region <- sprintf("%s:%d-%d", seqname, as.integer(start), as.integer(end))
  out <- tryCatch(
    system2("samtools", args = c("faidx", genome_fasta, region), stdout = TRUE, stderr = FALSE),
    error = function(e) NULL
  )
  if (is.null(out) || length(out) < 2) return(NULL)
  toupper(paste(out[-1], collapse = ""))
}

#' Translate a set of genomic CDS sub-ranges (already correctly identified --
#' e.g. build_rmats_arm_isoform()'s own `cds` output) into a protein
#' sequence: extracts each range's plus-strand DNA from the local genome
#' FASTA, concatenates in ascending genomic order, and reverse-complements
#' the WHOLE combined sequence if the transcript is "-" strand -- equivalent
#' to (and simpler than) reverse-complementing each range individually in
#' descending order, since revcomp(A+B+C) == revcomp(C)+revcomp(B)+revcomp(A).
#'
#' @param seqname chromosome/contig name
#' @param strand "+" or "-"
#' @param cds_df data.frame(start, end) genomic CDS sub-ranges, any order
#' @param genome_fasta path to the indexed genome FASTA
#' @return character(1) protein sequence, or NULL if any range's DNA
#'   couldn't be fetched or translation didn't produce a clean, complete CDS
#'   (see translate_dna_standard())
translate_cds_ranges <- function(seqname, strand, cds_df, genome_fasta = RMATS_GENOME_FASTA) {
  if (is.null(cds_df) || nrow(cds_df) == 0) return(NULL)
  cds_df <- cds_df[order(cds_df$start), ]
  fragments <- vapply(seq_len(nrow(cds_df)), function(i) {
    frag <- fetch_genome_dna(seqname, cds_df$start[i], cds_df$end[i], genome_fasta)
    if (is.null(frag)) NA_character_ else frag
  }, character(1))
  if (any(is.na(fragments))) return(NULL)
  combined <- paste(fragments, collapse = "")
  if (identical(strand, "-")) combined <- reverse_complement_dna(combined)
  translate_dna_standard(combined)
}

#' Builds a reference_exon_index-schema exon table (transcript_id, gene_id,
#' gene_name, is_canonical, exon_number, seqname, strand, start, end,
#' residue_start, residue_end, protein_length) for a CONSTRUCTED (synthetic)
#' isoform's CDS sub-ranges, so it row-binds unchanged alongside real
#' transcripts' rows into catalog()$exon_table for the shared MS1/MS2 exon
#' axis -- same approach as R/fasta_pipeline.R's build_novel_exon_table()
#' for Option 2's novel ORFs, adapted here to work directly from already-
#' known genomic CDS ranges rather than projecting a query-sequence ORF
#' range through a CIGAR alignment.
#'
#' @param cds_df data.frame(start, end) genomic CDS sub-ranges, any order
#' @param strand "+" or "-"
#' @param seqname chromosome/contig name
#' @param transcript_id id to assign this constructed transcript
#' @param gene_name matched gene name
#' @return data.frame, same columns/order as reference_exon_index
build_constructed_isoform_exon_table <- function(cds_df, strand, seqname, transcript_id, gene_name) {
  df <- cds_df[order(cds_df$start), , drop = FALSE]
  if (identical(strand, "-")) df <- df[rev(seq_len(nrow(df))), , drop = FALSE]
  nt_len <- df$end - df$start + 1L
  cum_nt <- cumsum(nt_len)
  prev_cum_nt <- cum_nt - nt_len
  df$residue_start <- prev_cum_nt %/% 3 + 1
  df$residue_end <- cum_nt %/% 3
  df$protein_length <- max(df$residue_end)
  df$exon_number <- seq_len(nrow(df))
  df$transcript_id <- transcript_id
  df$gene_id <- NA_character_
  df$gene_name <- gene_name
  df$is_canonical <- FALSE
  df$seqname <- seqname
  df$strand <- strand
  df[, c("transcript_id", "gene_id", "gene_name", "is_canonical", "exon_number",
         "seqname", "strand", "start", "end", "residue_start", "residue_end", "protein_length")]
}

#' Translates one arm's constructed isoform (build_rmats_arm_isoform()'s
#' output) all the way to a real protein sequence + exon table, so it can be
#' offered as a selectable candidate for "Add to comparison" alongside real
#' matched transcripts, not just shown in the exon-alignment visualization.
#'
#' @param synth result of build_rmats_arm_isoform()
#' @param transcript_id id to assign this constructed isoform (e.g.
#'   "CONSTRUCTED_exon1_ENST00000008180")
#' @param gene_name gene symbol
#' @return list(protein_sequence, exon_table) or NULL if translation failed
#'   (e.g. the spliced CDS wasn't a clean multiple of 3 with exactly one
#'   trailing stop -- can happen if the backbone's own CDS annotation and
#'   rMATS' local exon boundaries don't actually agree on frame)
translate_rmats_constructed_isoform <- function(synth, transcript_id, gene_name) {
  protein <- translate_cds_ranges(synth$seqname, synth$strand, synth$cds)
  if (is.null(protein)) return(NULL)
  exon_table <- build_constructed_isoform_exon_table(synth$cds, synth$strand, synth$seqname, transcript_id, gene_name)
  list(protein_sequence = protein, exon_table = exon_table)
}

#' Human-readable genomic region text for an rMATS event, for status
#' messages -- SE and RI each have one differential exon, MXE has two,
#' A5SS/A3SS report both alternative forms of one exon.
#'
#' @param event one parsed event row
#' @param event_type "SE", "MXE", "RI", "A5SS", or "A3SS"
rmats_event_region_text <- function(event, event_type) {
  if (identical(event_type, "MXE")) {
    sprintf("%s:%s-%s and %s:%s-%s", event$chr, event$exon1_start, event$exon1_end,
            event$chr, event$exon2_start, event$exon2_end)
  } else if (identical(event_type, "RI")) {
    sprintf("%s:%s-%s", event$chr, event$ri_start, event$ri_end)
  } else if (event_type %in% c("A5SS", "A3SS")) {
    sprintf("%s:%s-%s (long) / %s:%s-%s (short)", event$chr, event$long_start, event$long_end,
            event$chr, event$short_start, event$short_end)
  } else {
    sprintf("%s:%s-%s", event$chr, event$target_start, event$target_end)
  }
}

#' Builds the exon-alignment preview payload (same shape
#' PT.renderExonAlignment() expects, via build_multi_track_exon_alignment())
#' for one rMATS event, combining:
#'   - every REAL annotated transcript matched to any arm -- full exon
#'     structure fetched live via fetch_transcript_exons() (so 5'/3' UTRs
#'     render correctly, unlike reference_exon_index's CDS-only rows), with
#'     exon_index's own CDS rows as the coding sub-range overlay, and
#'   - one CONSTRUCTED synthetic isoform per arm (build_rmats_arm_isoform()),
#'     using backbone_choices[[arm_key]] if given, else
#'     default_backbone_for_arm()'s preset -- always present regardless of
#'     whether any real transcript matched, clearly labeled apart from real
#'     ENST tracks.
#' matches$highlight_regions (both flanking exons AND the cassette exon(s),
#' not just the differential one(s)) are carried through to
#' `payload$highlights`.
#'
#' @param matches result of match_rmats_se_transcripts()/
#'   match_rmats_mxe_transcripts()
#' @param exon_index the precomputed reference_exon_index
#' @param backbone_choices named list, arm key -> transcript_id to use as
#'   backbone for that arm's constructed isoform (falls back to
#'   default_backbone_for_arm() for any arm not given)
#' @return same shape as build_exon_alignment_preview(), plus `highlights`,
#'   or NULL if nothing at all could be built (no real matches AND no
#'   backbone available for any arm)
build_rmats_full_alignment <- function(matches, exon_index, backbone_choices = list()) {
  arm_of <- character(0)
  canonical_ids <- character(0)
  for (arm in matches$arms) {
    if (nrow(arm$candidates) == 0) next
    arm_of[arm$candidates$transcript_id] <- arm$label
    canonical_ids <- c(canonical_ids, arm$candidates$transcript_id[arm$candidates$is_canonical])
  }
  real_ids <- unique(names(arm_of))

  real_track_dfs <- lapply(real_ids, function(tid) {
    full <- tryCatch(fetch_transcript_exons(tid), error = function(e) NULL)
    if (is.null(full) || is.null(full$exons) || nrow(full$exons) == 0) return(NULL)
    full$exons[order(full$exons$start), ]
  })
  ok <- !vapply(real_track_dfs, is.null, logical(1))
  real_ids <- real_ids[ok]
  real_track_dfs <- real_track_dfs[ok]
  real_track_cds <- lapply(real_ids, function(tid) exon_index[exon_index$transcript_id == tid, c("start", "end")])
  real_labels <- vapply(real_ids, function(tid) {
    tags <- c(unname(arm_of[tid]), if (tid %in% canonical_ids) "canonical" else NA)
    tags <- tags[!is.na(tags)]
    if (length(tags) == 0) tid else sprintf("%s [%s]", tid, paste(tags, collapse = ", "))
  }, character(1), USE.NAMES = FALSE)

  arm_keys <- names(matches$arms)
  synth_track_dfs <- list(); synth_track_cds <- list(); synth_labels <- character(0)
  synth_seqname <- NA_character_; synth_strand <- NA_character_
  for (k in arm_keys) {
    backbone <- backbone_choices[[k]]
    if (is.null(backbone) || !nzchar(backbone)) backbone <- default_backbone_for_arm(matches, k)
    if (is.na(backbone)) next
    synth <- build_rmats_arm_isoform(backbone, matches$flanks[[k]], matches$cassette[[k]], exon_index)
    if (is.null(synth)) next
    synth_track_dfs[[length(synth_track_dfs) + 1]] <- synth$exons
    synth_track_cds[[length(synth_track_cds) + 1]] <- synth$cds
    synth_labels <- c(synth_labels, sprintf("Constructed: %s (backbone %s)", matches$arms[[k]]$label, backbone))
    if (is.na(synth_seqname)) { synth_seqname <- as.character(synth$seqname); synth_strand <- as.character(synth$strand) }
  }

  track_dfs <- c(real_track_dfs, synth_track_dfs)
  track_cds <- c(real_track_cds, synth_track_cds)
  track_labels <- c(real_labels, synth_labels)
  if (length(track_dfs) == 0) return(NULL)

  if (length(real_ids) > 0) {
    first_row <- exon_index[exon_index$transcript_id == real_ids[1], ][1, ]
    seqname <- as.character(first_row$seqname); strand <- as.character(first_row$strand)
  } else {
    seqname <- synth_seqname; strand <- synth_strand
  }

  build_multi_track_exon_alignment(track_dfs, track_labels, seqname, strand, track_cds = track_cds,
                                    highlight_regions = matches$highlight_regions)
}

#' Dispatch tables server.R uses to go from the user's chosen rMATS event
#' type (the "rmats_event_type" selectInput) to the right parser/matcher
#' pair, rather than a chain of if/else on the event type string at each
#' of several call sites (file parsing, event-list labeling, "Find matching
#' transcripts"). Keeping this here (not server.R) means adding a future
#' event type only ever requires touching this file.
RMATS_EVENT_TYPES <- c("SE", "MXE", "RI", "A5SS", "A3SS")
RMATS_PARSERS <- list(
  SE = parse_rmats_se, MXE = parse_rmats_mxe, RI = parse_rmats_ri,
  A5SS = parse_rmats_a5ss, A3SS = parse_rmats_a3ss
)
RMATS_MATCHERS <- list(
  SE = match_rmats_se_transcripts, MXE = match_rmats_mxe_transcripts, RI = match_rmats_ri_transcripts,
  A5SS = match_rmats_a5ss_transcripts, A3SS = match_rmats_a3ss_transcripts
)
