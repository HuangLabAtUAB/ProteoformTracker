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
# don't need to know or care which event type produced them. Only SE and
# MXE are implemented so far; A3SS/A5SS/RI are a planned follow-up using
# the same match_rmats_arm_transcripts()/build_rmats_arms_result() core.

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
#' @param gene_tbl the gene's own subset of reference_exon_index
#' @param flanks list(five_prime=, three_prime=), each list(start=,end=)
#' @param middle_exon NULL, or list(start=,end=) for the cassette exon
#'   expected between the two flanks
#' @return character vector of matching transcript ids (possibly empty)
match_rmats_arm_transcripts <- function(gene_tbl, flanks, middle_exon = NULL) {
  fp <- flanks$five_prime; tp <- flanks$three_prime
  matches_coord <- function(row, coord) !is.na(row$start) && row$start == coord$start && row$end == coord$end

  ids <- character(0)
  tx_ids <- unique(gene_tbl$transcript_id)
  for (tid in tx_ids) {
    tx <- gene_tbl[gene_tbl$transcript_id == tid, ]
    tx <- tx[order(tx$exon_number), ]
    n <- nrow(tx)
    if (is.null(middle_exon)) {
      if (n < 2) next
      for (i in seq_len(n - 1)) {
        if (matches_coord(tx[i, ], fp) && matches_coord(tx[i + 1, ], tp)) { ids <- c(ids, tid); break }
      }
    } else {
      if (n < 3) next
      for (i in seq_len(n - 2)) {
        if (matches_coord(tx[i, ], fp) &&
            tx$start[i + 1] == middle_exon$start && tx$end[i + 1] == middle_exon$end &&
            matches_coord(tx[i + 2, ], tp)) { ids <- c(ids, tid); break }
      }
    }
  }
  unique(ids)
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
#' @param arm_ids named list, arm key -> character vector of matched
#'   transcript ids for that arm
#' @param arm_labels named character vector, arm key -> human-readable label
#' @return list(arms = list(key = list(key=, label=, candidates=data.frame(
#'   transcript_id, is_canonical, pairing_score)), ...)) -- candidates is a
#'   0-row data.frame (never NULL) when an arm has no match
build_rmats_arms_result <- function(gene_tbl, arm_ids, arm_labels) {
  empty <- data.frame(transcript_id = character(0), is_canonical = logical(0),
                       pairing_score = numeric(0), stringsAsFactors = FALSE)
  arm_keys <- names(arm_ids)
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

  arms <- lapply(arm_keys, function(k) {
    ids <- unique(arm_ids[[k]])
    other_ids <- unique(unlist(arm_ids[setdiff(arm_keys, k)]))
    candidates <- if (length(ids) == 0) {
      empty
    } else if (length(other_ids) == 0) {
      data.frame(transcript_id = ids, is_canonical = vapply(ids, is_canonical_of, logical(1)),
                 pairing_score = 0, stringsAsFactors = FALSE)
    } else {
      other_keys <- lapply(other_ids, exon_key_set)
      scores <- vapply(ids, function(tid) {
        own <- exon_key_set(tid)
        max(vapply(other_keys, function(ok) jaccard(own, ok), numeric(1)))
      }, numeric(1))
      d <- data.frame(transcript_id = ids, is_canonical = vapply(ids, is_canonical_of, logical(1)),
                       pairing_score = round(unname(scores), 3), stringsAsFactors = FALSE)
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
  if (is.null(gene_tbl)) {
    result <- build_rmats_arms_result(NULL, list(inclusion = character(0), exclusion = character(0)), arm_labels)
    result$highlight_regions <- list()
    return(result)
  }
  flanks <- genomic_flanks_in_transcript_order(
    event$strand,
    list(start = event$flank_lo_start, end = event$flank_lo_end),
    list(start = event$flank_hi_start, end = event$flank_hi_end)
  )
  inclusion_ids <- match_rmats_arm_transcripts(gene_tbl, flanks, list(start = event$target_start, end = event$target_end))
  exclusion_ids <- match_rmats_arm_transcripts(gene_tbl, flanks, NULL)
  result <- build_rmats_arms_result(gene_tbl, list(inclusion = inclusion_ids, exclusion = exclusion_ids), arm_labels)
  result$highlight_regions <- list(list(start = event$target_start, end = event$target_end, label = "Differential exon (rMATS)"))
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
  if (is.null(gene_tbl)) {
    result <- build_rmats_arms_result(NULL, list(exon1 = character(0), exon2 = character(0)), arm_labels)
    result$highlight_regions <- list()
    return(result)
  }
  flanks <- genomic_flanks_in_transcript_order(
    event$strand,
    list(start = event$flank_lo_start, end = event$flank_lo_end),
    list(start = event$flank_hi_start, end = event$flank_hi_end)
  )
  exon1_ids <- match_rmats_arm_transcripts(gene_tbl, flanks, list(start = event$exon1_start, end = event$exon1_end))
  exon2_ids <- match_rmats_arm_transcripts(gene_tbl, flanks, list(start = event$exon2_start, end = event$exon2_end))
  result <- build_rmats_arms_result(gene_tbl, list(exon1 = exon1_ids, exon2 = exon2_ids), arm_labels)
  result$highlight_regions <- list(
    list(start = event$exon1_start, end = event$exon1_end, label = "1st mutually exclusive exon (rMATS)"),
    list(start = event$exon2_start, end = event$exon2_end, label = "2nd mutually exclusive exon (rMATS)")
  )
  result
}

#' Human-readable genomic region text for an rMATS event, for status
#' messages -- SE has one differential exon, MXE has two.
#'
#' @param event one parsed event row
#' @param event_type "SE" or "MXE"
rmats_event_region_text <- function(event, event_type) {
  if (identical(event_type, "MXE")) {
    sprintf("%s:%s-%s and %s:%s-%s", event$chr, event$exon1_start, event$exon1_end,
            event$chr, event$exon2_start, event$exon2_end)
  } else {
    sprintf("%s:%s-%s", event$chr, event$target_start, event$target_end)
  }
}

#' Builds the exon-alignment preview payload (same shape
#' PT.renderExonAlignment() expects, via build_multi_track_exon_alignment())
#' for the matched candidates of an rMATS event. Every track here is a
#' REAL, fully-annotated transcript with known coordinates already in
#' exon_index -- no REST fetch needed, and since reference_exon_index only
#' ever holds CDS exons (see its own doc comment), every block IS coding
#' sequence, so `coding` is set to the block's own full range rather than
#' left unknown. The event's own differential exon(s) (matches$
#' highlight_regions) are carried through to `payload$highlights` so the
#' renderer can box them on top of the usual common/partial/unique tier
#' coloring.
#'
#' @param matches result of match_rmats_se_transcripts()/
#'   match_rmats_mxe_transcripts()
#' @param exon_index the precomputed reference_exon_index
#' @return same shape as build_exon_alignment_preview(), plus a `highlights`
#'   field, or NULL if no candidate has exon rows
build_rmats_candidate_alignment <- function(matches, exon_index) {
  arm_of <- character(0)
  canonical_ids <- character(0)
  for (arm in matches$arms) {
    if (nrow(arm$candidates) == 0) next
    arm_of[arm$candidates$transcript_id] <- arm$label
    canonical_ids <- c(canonical_ids, arm$candidates$transcript_id[arm$candidates$is_canonical])
  }
  candidate_ids <- unique(names(arm_of))
  if (length(candidate_ids) == 0) return(NULL)

  track_dfs <- lapply(candidate_ids, function(tid) {
    sub <- exon_index[exon_index$transcript_id == tid, c("start", "end")]
    if (nrow(sub) == 0) NULL else sub[order(sub$start), ]
  })
  ok <- !vapply(track_dfs, is.null, logical(1))
  if (!any(ok)) return(NULL)
  candidate_ids <- candidate_ids[ok]
  track_dfs <- track_dfs[ok]

  track_labels <- vapply(candidate_ids, function(tid) {
    tags <- c(unname(arm_of[tid]), if (tid %in% canonical_ids) "canonical" else NA)
    tags <- tags[!is.na(tags)]
    if (length(tags) == 0) tid else sprintf("%s [%s]", tid, paste(tags, collapse = ", "))
  }, character(1), USE.NAMES = FALSE)

  first_row <- exon_index[exon_index$transcript_id == candidate_ids[1], ][1, ]
  # seqname/strand come out of reference_exon_index as factor columns --
  # as.character() them explicitly rather than relying on downstream code
  # (string comparisons, sprintf, JSON serialization) to coerce correctly.
  build_multi_track_exon_alignment(track_dfs, track_labels, as.character(first_row$seqname),
                                    as.character(first_row$strand), track_cds = track_dfs,
                                    highlight_regions = matches$highlight_regions)
}
