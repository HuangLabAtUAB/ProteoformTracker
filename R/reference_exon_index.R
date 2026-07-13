# Precomputed genome-wide exon structure index: every protein-coding
# transcript's CDS exons mapped to residue coordinates, built once from a
# bulk Ensembl GTF (not per-gene REST calls -- the same "precompute once,
# query fast" pattern as the mass and m/z indices), so any gene/transcript
# a user asks about can be looked up instantly instead of re-fetched.
#
# Covers ALL transcripts (not canonical-only) -- this index exists for the
# relevant-isoform axis (comparing a gene's isoforms against each other,
# e.g. CD44 canonical vs isoform 11), which canonical-only data can't serve.

#' Build the exon index from a bulk Ensembl GTF (already downloaded, e.g.
#' via scripts/build_exon_index.R). Only CDS features are used -- residue
#' coordinates come from cumulative CDS nucleotide length per transcript,
#' the same method validated by hand against CD44 (canonical vs isoform 11)
#' and cross-checked against direct sequence comparison.
#'
#' Canonical-transcript detection deliberately does NOT rely on
#' rtracklayer's parsed "tag" column: GTF allows a repeated `tag` attribute
#' per line, and rtracklayer's import() keeps only the last occurrence,
#' which silently drops "Ensembl_canonical" for transcripts where it isn't
#' the last tag listed. Canonical status is instead determined by a direct
#' text search over the raw GTF lines.
#'
#' @param gtf_path path to a (optionally CDS-prefiltered) GTF file
#' @param prefiltered if FALSE, prefilters to CDS rows of protein_coding
#'   genes into a temp file first (faster/lighter for rtracklayer to parse)
#' @return data.frame: transcript_id, gene_id, gene_name, is_canonical,
#'   exon_number, seqname, strand, genomic_start, genomic_end,
#'   residue_start, residue_end, protein_length (repeated per transcript)
build_reference_exon_index <- function(gtf_path, prefiltered = FALSE) {
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("package 'rtracklayer' is required to build the exon index")
  }

  cds_path <- gtf_path
  if (!prefiltered) {
    cds_path <- tempfile(fileext = ".gtf")
    status <- system(sprintf(
      "gzcat -f %s | awk -F'\t' '$3==\"CDS\" && /gene_biotype \"protein_coding\"/' > %s",
      shQuote(gtf_path), shQuote(cds_path)
    ))
    if (status != 0) stop("prefiltering the GTF to CDS rows failed")
  }

  message("Identifying canonical transcripts (direct text search, not rtracklayer's tag column)...")
  canonical_lines <- system(sprintf("grep Ensembl_canonical %s", shQuote(cds_path)), intern = TRUE)
  canonical_ids <- unique(sub('.*transcript_id "([^"]+)".*', "\\1", canonical_lines))

  message("Parsing CDS coordinates with rtracklayer...")
  gtf <- rtracklayer::import(cds_path, format = "gtf")
  df <- as.data.frame(gtf)
  df <- df[, c("seqnames", "start", "end", "strand", "gene_id", "gene_name", "transcript_id", "exon_number")]
  df$exon_number <- as.integer(df$exon_number)
  df$is_canonical <- df$transcript_id %in% canonical_ids

  message("Computing residue coordinates per transcript...")
  df <- df[order(df$transcript_id, df$exon_number), ]
  nt_len <- df$end - df$start + 1
  cum_nt <- ave(nt_len, df$transcript_id, FUN = cumsum)
  prev_cum_nt <- cum_nt - nt_len
  df$residue_start <- prev_cum_nt %/% 3 + 1
  df$residue_end <- cum_nt %/% 3

  protein_length <- ave(df$residue_end, df$transcript_id, FUN = max)
  df$protein_length <- protein_length

  names(df)[names(df) == "seqnames"] <- "seqname"
  rownames(df) <- NULL
  df[, c(
    "transcript_id", "gene_id", "gene_name", "is_canonical", "exon_number",
    "seqname", "strand", "start", "end", "residue_start", "residue_end", "protein_length"
  )]
}

#' Load a previously built exon index.
load_reference_exon_index <- function(path) {
  readRDS(path)
}

#' Look up exon structure by gene symbol or transcript id.
#'
#' @param index exon index (from build_/load_reference_exon_index())
#' @param gene_symbol gene symbol, e.g. "CD44" (returns all its transcripts)
#' @param transcript_id specific Ensembl transcript id (returns just that one)
#' @param canonical_only if TRUE (with gene_symbol), only the canonical transcript
#' @return data.frame subset of the index, one row per exon
get_exon_structure <- function(index, gene_symbol = NULL, transcript_id = NULL, canonical_only = FALSE) {
  if (is.null(gene_symbol) && is.null(transcript_id)) {
    stop("get_exon_structure() requires gene_symbol or transcript_id")
  }
  hits <- if (!is.null(transcript_id)) {
    index[index$transcript_id == transcript_id, ]
  } else {
    hits <- index[!is.na(index$gene_name) & index$gene_name == gene_symbol, ]
    if (canonical_only) hits <- hits[hits$is_canonical, ]
    hits
  }
  hits[order(hits$transcript_id, hits$exon_number), ]
}
