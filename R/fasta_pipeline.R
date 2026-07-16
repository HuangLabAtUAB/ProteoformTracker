# Option 2 pipeline: user-pasted/uploaded FASTA (spliced mRNA/cDNA or CDS
# nucleotide sequence) -> TransDecoder ORF translation -> minimap2 spliced
# genome alignment -> real CDS exon structure -> overlap search against the
# precomputed reference_exon_index to identify which known gene/isoforms
# this sequence belongs to, so it can be compared exactly like Option 1's
# gene/isoform pathway.
#
# Deliberately does NOT use StringTie or TransDecoder's genome-lifting Perl
# utilities (gtf_to_alignment_gff3.pl / cdna_alignment_orf_to_genome_orf.pl,
# which IsoPepTracker calls with hardcoded /usr/lib/transdecoder/util paths
# that don't exist on this machine): those exist to reconstruct transcript
# structure from many short RNA-seq reads. Here the input is already one
# full-length transcript sequence, so its own minimap2 CIGAR string directly
# *is* the exon structure -- no assembly step needed.

.GENETIC_CODE_STANDARD <- "Universal"

#' Parse a SAM CIGAR string into (op, length) pairs.
#' @param cigar CIGAR string, e.g. "1364M824N107M"
#' @return data.frame(op, length)
parse_cigar <- function(cigar) {
  m <- gregexpr("[0-9]+[MIDNSHP=X]", cigar)[[1]]
  toks <- regmatches(cigar, gregexpr("[0-9]+[MIDNSHP=X]", cigar))[[1]]
  ops <- sub("^[0-9]+", "", toks)
  lens <- as.integer(sub("[MIDNSHP=X]$", "", toks))
  data.frame(op = ops, length = lens, stringsAsFactors = FALSE)
}

#' Walk a CIGAR string (given the SAM POS of its first aligned base) and
#' return one row per contiguous aligned segment (an intron-bounded exon in
#' genomic coordinates), each carrying its constituent M/=/X micro-blocks
#' (qstart, qend, gstart, gend) so ORF (query) coordinates can later be
#' projected onto genomic coordinates within that segment.
#'
#' @param cigar CIGAR string
#' @param pos 1-based leftmost genomic mapping position (SAM POS field)
#' @return list of segments, each list(gstart, gend, blocks = data.frame(qstart, qend, gstart, gend))
cigar_to_segments <- function(cigar, pos) {
  ops <- parse_cigar(cigar)
  qpos <- 1L
  gpos <- as.integer(pos)
  segments <- list()
  cur_blocks <- list()
  for (i in seq_len(nrow(ops))) {
    op <- ops$op[i]; len <- ops$length[i]
    if (op %in% c("M", "=", "X")) {
      cur_blocks[[length(cur_blocks) + 1]] <- data.frame(
        qstart = qpos, qend = qpos + len - 1L, gstart = gpos, gend = gpos + len - 1L
      )
      qpos <- qpos + len; gpos <- gpos + len
    } else if (op %in% c("I", "S")) {
      qpos <- qpos + len
    } else if (op == "D") {
      gpos <- gpos + len
    } else if (op == "N") {
      if (length(cur_blocks) > 0) {
        blocks <- do.call(rbind, cur_blocks)
        segments[[length(segments) + 1]] <- list(gstart = min(blocks$gstart), gend = max(blocks$gend), blocks = blocks)
        cur_blocks <- list()
      }
      gpos <- gpos + len
    }
    # H and P consume neither query nor reference position for our purposes
  }
  if (length(cur_blocks) > 0) {
    blocks <- do.call(rbind, cur_blocks)
    segments[[length(segments) + 1]] <- list(gstart = min(blocks$gstart), gend = max(blocks$gend), blocks = blocks)
  }
  segments
}

#' Run TransDecoder.LongOrfs + TransDecoder.Predict on a FASTA file already
#' written to `work_dir`, returning ALL scored candidate ORFs (not just one)
#' sorted by TransDecoder's own coding-likelihood score, highest first.
#'
#' Deliberately independent of minimap2/gene-matching: this is purely a
#' translation step over the sequence's own nucleotides. A single "best"
#' ORF is a translation heuristic (Markov/hexamer scoring) that can pick a
#' short or incomplete candidate over the true CDS for a noisy or
#' partially-evidenced sequence -- surfacing the top candidates lets the
#' user judge which one (if any) looks like the real protein, using
#' completeness (start/stop codon presence) and length as evidence.
#'
#' @param fasta_path path to the input FASTA (single sequence)
#' @param work_dir directory to run in (TransDecoder writes many sidecar files here)
#' @param min_protein_length minimum ORF length in amino acids
#' @param max_candidates cap on how many candidates to return (highest-scoring first)
#' @return data.frame(orf_id, type, quality, has_start, has_stop, aa_len,
#'   score, nt_start, nt_end, orf_strand, protein_sequence), 0 rows if none found
list_transdecoder_orfs <- function(fasta_path, work_dir, min_protein_length = 30, max_candidates = 5) {
  fasta_name <- basename(fasta_path)
  args1 <- c("-t", fasta_name, "-m", as.character(min_protein_length))
  res1 <- system2("TransDecoder.LongOrfs", args = args1, wait = TRUE,
                   stdout = file.path(work_dir, "longorfs.log"), stderr = file.path(work_dir, "longorfs.err"))
  if (res1 != 0) {
    stop("TransDecoder.LongOrfs failed (exit code ", res1, "). See longorfs.err for details.")
  }

  args2 <- c("-t", fasta_name, "--no_refine_starts")
  res2 <- system2("TransDecoder.Predict", args = args2, wait = TRUE,
                   stdout = file.path(work_dir, "predict.log"), stderr = file.path(work_dir, "predict.err"))
  if (res2 != 0) {
    stop("TransDecoder.Predict failed (exit code ", res2, "). See predict.err for details.")
  }

  pep_path <- file.path(work_dir, paste0(fasta_name, ".transdecoder.pep"))
  if (!file.exists(pep_path) || file.size(pep_path) == 0) {
    return(data.frame(
      orf_id = character(0), type = character(0), quality = character(0),
      has_start = logical(0), has_stop = logical(0), aa_len = integer(0),
      score = numeric(0), nt_start = integer(0), nt_end = integer(0),
      orf_strand = character(0), protein_sequence = character(0)
    ))
  }
  pep_lines <- readLines(pep_path)
  header_idx <- which(grepl("^>", pep_lines))

  quality_labels <- c(
    complete = "Complete (start + stop)",
    "5prime_partial" = "Missing start codon",
    "3prime_partial" = "Missing stop codon",
    internal = "Missing start and stop"
  )

  rows <- lapply(seq_along(header_idx), function(i) {
    header <- pep_lines[header_idx[i]]
    seq_end <- if (i < length(header_idx)) header_idx[i + 1] - 1 else length(pep_lines)
    protein_sequence <- paste(pep_lines[(header_idx[i] + 1):seq_end], collapse = "")
    protein_sequence <- sub("\\*$", "", protein_sequence)

    m <- regmatches(header, regexec(":([0-9]+)-([0-9]+)\\(([+-])\\)\\s*$", header))[[1]]
    score_m <- regmatches(header, regexec("score=([-0-9.]+)", header))[[1]]
    orf_id_m <- regmatches(header, regexec("^>(\\S+)", header))[[1]]
    type <- sub(".*type:(\\S+).*", "\\1", header)

    data.frame(
      orf_id = orf_id_m[2],
      type = type,
      quality = unname(quality_labels[type]),
      has_start = type %in% c("complete", "3prime_partial"),
      has_stop = type %in% c("complete", "5prime_partial"),
      aa_len = nchar(protein_sequence),
      score = as.numeric(score_m[2]),
      nt_start = as.integer(m[2]),
      nt_end = as.integer(m[2]) + nchar(protein_sequence) * 3L - 1L,
      orf_strand = m[4],
      protein_sequence = protein_sequence,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)
  df <- df[order(-df$score), ]
  rownames(df) <- NULL
  head(df, max_candidates)
}

#' Spliced-align the input FASTA against the GRCh38 minimap2 index and parse
#' the primary alignment's genomic exon segments.
#'
#' @param fasta_path path to the input FASTA
#' @param work_dir directory to run in
#' @param mmi_path path to the minimap2 genome index
#' @return list(seqname, strand, mapq, segments) from cigar_to_segments(), or NULL if unmapped
run_minimap2_alignment <- function(fasta_path, work_dir, mmi_path) {
  sam_path <- file.path(work_dir, "aligned.sam")
  res <- system2("minimap2", args = c("-ax", "splice", "--secondary=no", mmi_path, fasta_path),
                  stdout = sam_path, stderr = file.path(work_dir, "minimap2.err"), wait = TRUE)
  if (res != 0) stop("minimap2 failed (exit code ", res, "). See minimap2.err for details.")

  sam_lines <- readLines(sam_path)
  aln_lines <- sam_lines[!grepl("^@", sam_lines)]
  if (length(aln_lines) == 0) return(NULL)
  fields <- strsplit(aln_lines[1], "\t")[[1]]
  flag <- as.integer(fields[2])
  if (bitwAnd(flag, 4L) != 0) return(NULL)  # unmapped
  seqname <- fields[3]
  pos <- as.integer(fields[4])
  mapq <- as.integer(fields[5])
  cigar <- fields[6]
  strand <- if (bitwAnd(flag, 16L) != 0) "-" else "+"

  list(seqname = seqname, strand = strand, mapq = mapq, segments = cigar_to_segments(cigar, pos))
}

#' Convert an ORF's nucleotide range -- reported by TransDecoder relative to
#' the original input sequence exactly as given (5' to 3', regardless of
#' which strand of that input the ORF itself reads from) -- into the
#' coordinate frame of the SAM record's CIGAR/query positions. SAM always
#' lists CIGAR left-to-right in genomic-forward order; when the alignment's
#' FLAG marks it reverse (0x10), that means the input sequence had to be
#' reverse-complemented to align, so position 1 of the CIGAR's query frame
#' is actually the LAST base of the original input, not the first -- the
#' range must be flipped end-for-end before it can be overlapped against
#' cigar_to_segments()'s (qstart, qend) blocks.
#'
#' @param nt_start,nt_end 1-based inclusive range in original-input coordinates
#' @param query_length total length of the original input sequence
#' @param genomic_strand "+" or "-", from the SAM alignment FLAG
#' @return list(nt_start, nt_end) in SAM query-coordinate frame
to_sam_query_coords <- function(nt_start, nt_end, query_length, genomic_strand) {
  if (genomic_strand == "-") {
    list(nt_start = query_length - nt_end + 1L, nt_end = query_length - nt_start + 1L)
  } else {
    list(nt_start = nt_start, nt_end = nt_end)
  }
}

#' Project the ORF's nucleotide range onto the genomic exon segments to
#' build a CDS-only exon table in the same schema as reference_exon_index,
#' so it can be row-bound alongside known transcripts and fed into
#' build_shared_exon_axis()/build_section1_payload() unchanged.
#'
#' @param segments output of cigar_to_segments()
#' @param strand genomic strand ("+" or "-") from the alignment
#' @param seqname chromosome/contig name
#' @param nt_start,nt_end 1-based inclusive ORF nucleotide range (coding
#'   portion only, stop codon excluded) within the input query sequence
#' @param transcript_id id to assign this novel transcript
#' @param gene_name matched gene name, or NA if none found
#' @return data.frame with the same columns as reference_exon_index, or NULL if no overlap
build_novel_exon_table <- function(segments, strand, seqname, nt_start, nt_end,
                                    transcript_id = "NOVEL_1", gene_name = NA_character_) {
  exon_rows <- list()
  for (seg in segments) {
    blocks <- seg$blocks
    ov_qs <- pmax(blocks$qstart, nt_start)
    ov_qe <- pmin(blocks$qend, nt_end)
    keep <- ov_qs <= ov_qe
    if (!any(keep)) next
    blocks <- blocks[keep, , drop = FALSE]
    ov_qs <- ov_qs[keep]; ov_qe <- ov_qe[keep]
    ov_gs <- blocks$gstart + (ov_qs - blocks$qstart)
    ov_ge <- blocks$gstart + (ov_qe - blocks$qstart)
    exon_rows[[length(exon_rows) + 1]] <- data.frame(
      start = min(ov_gs), end = max(ov_ge), nt_len = sum(ov_qe - ov_qs + 1L)
    )
  }
  if (length(exon_rows) == 0) return(NULL)
  df <- do.call(rbind, exon_rows)
  df <- df[order(df$start), ]
  if (strand == "-") df <- df[rev(seq_len(nrow(df))), ]

  cum_nt <- cumsum(df$nt_len)
  prev_cum_nt <- cum_nt - df$nt_len
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

#' Search the precomputed exon index for known transcripts whose exons
#' overlap the novel transcript's genomic footprint, on the same
#' seqname/strand -- this is the "find relevant transcripts sharing exons or
#' coding sequence" step.
#'
#' Deliberately searches the FULL aligned exon footprint (every segment
#' minimap2 reported), not just the CDS-restricted subset that overlaps
#' wherever TransDecoder's single best-scoring ORF happened to land.
#' TransDecoder's ORF choice is a translation heuristic and can pick a
#' short/incomplete candidate over the true CDS (e.g. for a noisy or
#' partially-evidenced long-read sequence) -- gene identity is a property of
#' the whole aligned transcript's genomic location, not of wherever that one
#' heuristic-picked ORF sits, so it must not be gated on it.
#'
#' @param exon_index reference_exon_index (or NULL if not loaded)
#' @param seqname,strand genomic location of the novel transcript
#' @param segments output of cigar_to_segments() -- the full alignment, not
#'   the CDS-restricted exon table
#' @return list(gene_name, matched_transcript_ids, n_overlap_exons) or NULL if no match
find_matched_gene <- function(exon_index, seqname, strand, segments) {
  if (is.null(exon_index) || length(segments) == 0) return(NULL)
  cand <- exon_index[exon_index$seqname == seqname & exon_index$strand == strand, ]
  if (nrow(cand) == 0) return(NULL)

  overlaps <- logical(nrow(cand))
  for (seg in segments) {
    overlaps <- overlaps | (cand$start <= seg$gend & cand$end >= seg$gstart)
  }
  cand <- cand[overlaps, ]
  if (nrow(cand) == 0) return(NULL)

  gene_counts <- sort(table(cand$gene_name), decreasing = TRUE)
  best_gene <- names(gene_counts)[1]
  matched_transcript_ids <- unique(cand$transcript_id[cand$gene_name == best_gene])
  list(gene_name = best_gene, matched_transcript_ids = matched_transcript_ids, n_overlap_exons = gene_counts[[1]])
}

#' Clean raw FASTA text (pasted or uploaded) down to just its nucleotide
#' sequence, robust to multi-line input (see the header-stripping bug this
#' replaced: a whole-string regex silently failed on multi-line text and let
#' stray letters like the "N" in ">ENST..." leak into the sequence).
#'
#' @param fasta_text raw FASTA text, header optional
#' @return character(1) uppercased nucleotide-only sequence
clean_fasta_sequence <- function(fasta_text) {
  lines <- strsplit(trimws(fasta_text), "\n")[[1]]
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0) stop("No sequence provided.")
  seq_lines <- lines[!grepl("^>", lines)]
  seq_only <- toupper(gsub("[^ACGTUNacgtun]", "", paste(seq_lines, collapse = "")))
  if (nchar(seq_only) < 30) stop("Sequence is too short (need at least 30 nt of coding sequence).")
  if (nchar(seq_only) > 200000) stop("Sequence is too long (200,000 nt limit) -- provide a single spliced transcript/CDS sequence, not a whole chromosome.")
  seq_only
}

#' Step 2 of Option 2 (independent of run_fasta_alignment()/minimap2):
#' translate the sequence with TransDecoder and return its top candidate
#' ORFs. Manages its own temp working directory, mirroring
#' run_fasta_alignment()'s pattern.
#'
#' @param fasta_text raw FASTA text (header optional)
#' @param transcript_id id to assign the resulting novel proteoform
#' @param max_candidates cap on how many candidates to return
#' @return data.frame, see list_transdecoder_orfs()
run_fasta_translation <- function(fasta_text, transcript_id = "NOVEL_1", max_candidates = 5) {
  seq_only <- clean_fasta_sequence(fasta_text)

  work_dir <- file.path(tempdir(), paste0("ptracker_fasta_orf_", as.integer(Sys.time()), "_", sample.int(1e6, 1)))
  dir.create(work_dir, recursive = TRUE)
  fasta_path <- file.path(work_dir, "input.fa")
  writeLines(c(paste0(">", transcript_id), seq_only), fasta_path)

  original_wd <- getwd()
  on.exit(setwd(original_wd), add = TRUE)
  setwd(work_dir)

  list_transdecoder_orfs(basename(fasta_path), work_dir, max_candidates = max_candidates)
}

#' Step 1 of Option 2, run independently of TransDecoder: spliced-align the
#' sequence against GRCh38 with minimap2 and identify which known gene (and
#' its known transcripts) it belongs to, from the FULL aligned exon
#' footprint. Nothing here depends on translation.
#'
#' @param fasta_text raw FASTA text (header optional)
#' @param exon_index reference_exon_index (for the overlap search), or NULL
#' @param mmi_path path to the minimap2 GRCh38 index
#' @param transcript_id id to assign the resulting novel proteoform
#' @return list(seq_only, seqname, strand, mapq, segments, matched_gene, warnings)
run_fasta_alignment <- function(fasta_text, exon_index, mmi_path = "reference/genome/GRCh38.mmi",
                                 transcript_id = "NOVEL_1") {
  warnings_out <- character(0)
  seq_only <- clean_fasta_sequence(fasta_text)

  work_dir <- file.path(tempdir(), paste0("ptracker_fasta_aln_", as.integer(Sys.time()), "_", sample.int(1e6, 1)))
  dir.create(work_dir, recursive = TRUE)
  fasta_path <- file.path(work_dir, "input.fa")
  writeLines(c(paste0(">", transcript_id), seq_only), fasta_path)

  original_wd <- getwd()
  mmi_abs <- normalizePath(mmi_path, mustWork = TRUE)
  on.exit(setwd(original_wd), add = TRUE)
  setwd(work_dir)

  aln <- run_minimap2_alignment(basename(fasta_path), work_dir, mmi_abs)
  if (is.null(aln)) {
    return(list(seq_only = seq_only, seqname = NA_character_, strand = NA_character_, mapq = NA_integer_,
                segments = list(), matched_gene = NULL,
                warnings = "Could not align this sequence to the reference genome (GRCh38)."))
  }

  matched <- find_matched_gene(exon_index, aln$seqname, aln$strand, aln$segments)
  if (is.null(matched)) {
    span <- do.call(rbind, lapply(aln$segments, function(s) data.frame(gstart = s$gstart, gend = s$gend)))
    warnings_out <- c(warnings_out, sprintf(
      "No known transcript overlaps this sequence's alignment (%s:%d-%d, strand %s) -- it may be a novel locus, or not human/not in the reference build.",
      aln$seqname, min(span$gstart), max(span$gend), aln$strand
    ))
  }

  list(seq_only = seq_only, seqname = aln$seqname, strand = aln$strand, mapq = aln$mapq,
       segments = aln$segments, matched_gene = matched, warnings = warnings_out)
}

#' Step 2 (independent of step 1): once the user has picked one of the
#' TransDecoder candidates from list_transdecoder_orfs(), project just that
#' candidate's nucleotide range onto the alignment's genomic exon segments
#' to build its CDS exon table -- same schema as reference_exon_index, so it
#' row-binds with a known transcript's rows unchanged.
#'
#' @param orf_row one row of list_transdecoder_orfs()'s result (a data.frame
#'   with 1 row, e.g. candidates[i, ])
#' @param alignment result of run_fasta_alignment()
#' @param transcript_id id to assign this novel transcript
#' @param gene_name matched gene name, or NA
#' @return data.frame (see build_novel_exon_table()), or NULL if no overlap
#'   between the selected ORF and the aligned exons
build_selected_orf_exon_table <- function(orf_row, alignment, transcript_id = "NOVEL_1", gene_name = NA_character_) {
  if (length(alignment$segments) == 0) return(NULL)
  sam_coords <- to_sam_query_coords(orf_row$nt_start, orf_row$nt_end, nchar(alignment$seq_only), alignment$strand)
  build_novel_exon_table(alignment$segments, alignment$strand, alignment$seqname,
                          sam_coords$nt_start, sam_coords$nt_end, transcript_id, gene_name)
}

#' Build a multi-row exon-presence comparison between the novel sequence's
#' full minimap2 alignment footprint and one or more known transcripts' full
#' exon lists (UTR included, via fetch_transcript_exons() -- NOT
#' reference_exon_index, which is CDS-only and would show spurious "gaps"
#' for real UTR exons it simply never indexed). Independent of any
#' ORF/TransDecoder choice.
#'
#' Genomic coverage (any track's exon at all) is merged into "coverage
#' super-intervals" purely to decide WHERE to compress -- a real intron gap
#' between super-intervals collapses to a small fixed-width spacer, but
#' *inside* a super-interval, genomic position maps directly (1bp = 1 axis
#' unit, strand-aware) to axis position. Each track's own block is rendered
#' at its OWN precise genomic boundaries projected into that linear space --
#' never the full super-interval -- so two tracks' exons that overlap but
#' don't share identical boundaries (e.g. the novel sequence's own
#' minimap2 segment vs. a known transcript's finer real exon) land at
#' genuinely overlapping (not merged-away, not force-identical) axis
#' positions, the same way a real genome-browser multi-track view would
#' show them. Each row draws a thin connecting line across its FULL width
#' first (the "intron"), then draws its own exon blocks on top -- a
#' super-interval stretch a track doesn't cover simply shows the bare line
#' through it, which is the "placed gap" alignment.
#'
#' @param alignment result of run_fasta_alignment()
#' @param known_transcript_ids character vector of Ensembl transcript ids to
#'   compare against (one or more)
#' @param exon_index reference_exon_index (CDS-only genomic ranges), used to
#'   mark which part of each known transcript's exon is actually coding --
#'   that index is built from GTF CDS rows alone, so its rows already ARE
#'   the coding sub-ranges. Pass NULL to skip coding/non-coding marking for
#'   known transcripts.
#' @param novel_label display label for the novel-sequence row
#' @param known_labels display labels for each known transcript row
#'   (defaults to the transcript ids themselves)
#' @param novel_cds_range data.frame(start, end) (or NULL) -- the novel
#'   sequence's own coding genomic sub-ranges, i.e. build_selected_orf_exon_table()'s
#'   result for whichever TransDecoder candidate is currently selected. NULL
#'   if no candidate is selected yet (the novel row then renders with no
#'   coding/non-coding distinction, since it isn't known).
#' @return list(seqname, strand, axis_length, tracks = list(list(label,
#'   blocks), ...)) with tracks[[1]] always the novel sequence, or NULL if
#'   the novel sequence didn't align or no known transcript's exons could be
#'   fetched. Each block is list(start, end, coding = list(list(start,end),
#'   ...), shared_count) where `coding` gives the axis sub-range(s) of that
#'   block which are coding sequence (absent entirely if coding status is
#'   unknown for that track; an empty list means the whole block is
#'   non-coding UTR), and `shared_count` is how many compared tracks
#'   (including this one) have real genomic overlap at this exon
build_exon_alignment_preview <- function(alignment, known_transcript_ids, exon_index = NULL,
                                          novel_label = "Your sequence", known_labels = known_transcript_ids,
                                          novel_cds_range = NULL) {
  if (length(alignment$segments) == 0) return(NULL)
  if (length(known_labels) != length(known_transcript_ids)) known_labels <- known_transcript_ids

  novel_df <- do.call(rbind, lapply(alignment$segments, function(s) data.frame(start = s$gstart, end = s$gend)))
  known_dfs <- lapply(known_transcript_ids, function(tid) {
    k <- tryCatch(fetch_transcript_exons(tid), error = function(e) NULL)
    if (is.null(k)) NULL else k$exons
  })
  ok <- !vapply(known_dfs, is.null, logical(1))
  known_transcript_ids <- known_transcript_ids[ok]
  known_labels <- known_labels[ok]
  known_dfs <- known_dfs[ok]
  if (length(known_dfs) == 0) return(NULL)

  track_dfs <- c(list(novel_df), known_dfs)
  track_labels <- c(novel_label, known_labels)
  track_ids <- c(NA_character_, known_transcript_ids)
  n_tracks <- length(track_dfs)

  track_cds <- lapply(seq_along(track_ids), function(ti) {
    if (ti == 1) {
      if (is.null(novel_cds_range) || nrow(novel_cds_range) == 0) NULL else novel_cds_range[, c("start", "end")]
    } else if (!is.null(exon_index)) {
      sub <- exon_index[exon_index$transcript_id == track_ids[ti], c("start", "end")]
      if (nrow(sub) == 0) NULL else sub
    } else {
      NULL
    }
  })

  # Coverage super-intervals: merge any OVERLAPPING exon boundaries across
  # all tracks (regardless of exact match) so real intron gaps -- stretches
  # no track has any exon in -- are the only thing compressed.
  all_exons <- do.call(rbind, lapply(track_dfs, function(d) d[, c("start", "end")]))
  all_exons <- all_exons[order(all_exons$start), ]
  supers <- list()
  cur_start <- all_exons$start[1]; cur_end <- all_exons$end[1]
  if (nrow(all_exons) > 1) {
    for (i in 2:nrow(all_exons)) {
      if (all_exons$start[i] <= cur_end) {
        cur_end <- max(cur_end, all_exons$end[i])
      } else {
        supers[[length(supers) + 1]] <- data.frame(start = cur_start, end = cur_end)
        cur_start <- all_exons$start[i]; cur_end <- all_exons$end[i]
      }
    }
  }
  supers[[length(supers) + 1]] <- data.frame(start = cur_start, end = cur_end)
  super_df <- do.call(rbind, supers)
  super_df$length <- super_df$end - super_df$start + 1

  strand <- alignment$strand
  if (strand == "-") super_df <- super_df[rev(seq_len(nrow(super_df))), ]

  n_supers <- nrow(super_df)
  gap <- if (n_supers > 1) max(30, round(0.12 * sum(super_df$length) / (n_supers - 1))) else 0
  super_df$axis_end <- cumsum(super_df$length + gap) - gap
  super_df$axis_start <- super_df$axis_end - super_df$length + 1
  axis_length <- super_df$axis_end[n_supers]

  # Map one genomic position to axis coordinates within its super-interval.
  # Strand-aware: on "+" strand, low genomic = 5'/left, so position moves
  # the same direction as axis position. On "-" strand, transcription runs
  # HIGH genomic to LOW, so within a super-interval the HIGH genomic end is
  # 5'/left and the LOW end is 3'/right -- the offset must be measured from
  # the super-interval's end, not its start, or a coding sub-range near a
  # minus-strand exon's low-genomic (3') boundary would incorrectly land on
  # the LEFT (5') side of the block instead of the right.
  to_axis <- function(gpos, si) {
    if (strand == "-") super_df$axis_start[si] + (super_df$end[si] - gpos)
    else super_df$axis_start[si] + (gpos - super_df$start[si])
  }
  find_super <- function(s, e) which(super_df$start <= e & super_df$end >= s)[1]

  # Each track's own exons, individually, at their OWN precise genomic
  # boundaries -- never widened to the full super-interval.
  track_blocks <- lapply(seq_along(track_dfs), function(ti) {
    df <- track_dfs[[ti]]
    cds <- track_cds[[ti]]
    blocks <- list()
    for (r in seq_len(nrow(df))) {
      s <- df$start[r]; e <- df$end[r]
      si <- find_super(s, e)
      if (is.na(si)) next
      ax_a <- to_axis(s, si); ax_b <- to_axis(e, si)
      block <- list(start = min(ax_a, ax_b), end = max(ax_a, ax_b), gstart = s, gend = e)
      if (!is.null(cds)) {
        ov_s <- pmax(cds$start, s); ov_e <- pmin(cds$end, e)
        keep <- ov_s <= ov_e
        block$coding <- if (any(keep)) {
          lapply(which(keep), function(j) {
            ca <- to_axis(ov_s[j], si); cb <- to_axis(ov_e[j], si)
            list(start = min(ca, cb), end = max(ca, cb))
          })
        } else {
          list()
        }
      }
      blocks[[length(blocks) + 1]] <- block
    }
    blocks
  })

  # Tier ("shared") count: real genomic overlap with each OTHER track's own
  # exons -- not axis-position matching, since two tracks' exons here can
  # genuinely have different (but overlapping) boundaries.
  tracks <- lapply(seq_along(track_dfs), function(ti) {
    blocks <- lapply(track_blocks[[ti]], function(b) {
      shared_count <- 1L
      for (tj in seq_along(track_dfs)) {
        if (tj == ti) next
        if (any(track_dfs[[tj]]$start <= b$gend & track_dfs[[tj]]$end >= b$gstart)) shared_count <- shared_count + 1L
      }
      list(start = b$start, end = b$end, coding = b$coding, shared_count = shared_count)
    })
    list(label = track_labels[ti], blocks = blocks)
  })

  list(seqname = alignment$seqname, strand = strand, axis_length = axis_length, tracks = tracks, n_tracks = n_tracks)
}
