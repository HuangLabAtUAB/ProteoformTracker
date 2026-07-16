# Builds the real isoform catalog + proteoform objects for Option 1 (gene ->
# isoform -> proteoform selection). Unlike the earlier standalone prototype
# (which spliced isoform sequences out of one hand-fetched canonical
# sequence, and only worked for transcripts sharing canonical's exact splice
# sites), this version fetches each transcript's real, independently
# translated protein sequence from Ensembl REST (R/ensembl_protein_fetch.R),
# so it works for every protein-coding transcript of a gene, not just a
# "clean exon subset" of one reference sequence.

#' Real isoform catalog for a gene: every protein-coding transcript in the
#' precomputed genome-wide exon index, sorted longest-first.
#'
#' @param gene_symbol HGNC gene symbol, e.g. "CD44"
#' @param exon_index the loaded reference_exon_index (see R/reference_exon_index.R)
#' @param max_isoforms cap on how many transcripts to return (some genes have
#'   very large numbers of annotated transcripts; the UI needs a sane limit)
#' @return NULL if the gene isn't found, else list(gene, exon_table,
#'   transcript_ids, protein_lengths) with transcript_ids ordered longest-first
build_gene_isoform_catalog <- function(gene_symbol, exon_index, max_isoforms = 40) {
  tbl <- get_exon_structure(exon_index, gene_symbol = gene_symbol)
  if (is.null(tbl) || nrow(tbl) == 0) return(NULL)

  tx_ids <- unique(tbl$transcript_id)
  protein_lengths <- vapply(tx_ids, function(tid) {
    max(tbl$residue_end[tbl$transcript_id == tid])
  }, numeric(1))

  ord <- order(-protein_lengths)
  tx_ids <- tx_ids[ord]
  protein_lengths <- protein_lengths[ord]

  if (length(tx_ids) > max_isoforms) {
    tx_ids <- tx_ids[seq_len(max_isoforms)]
    protein_lengths <- protein_lengths[seq_len(max_isoforms)]
  }

  list(gene = gene_symbol, exon_table = tbl,
       transcript_ids = tx_ids, protein_lengths = protein_lengths)
}

#' Compress a sorted vector of exon numbers into a compact range string,
#' e.g. c(1,2,3,5,7,8) -> "1-3,5,7-8", for display in the isoform table.
#'
#' @param exon_numbers integer vector of 1-based exon numbers
compress_exon_ranges <- function(exon_numbers) {
  v <- sort(unique(exon_numbers))
  if (length(v) == 0) return("")
  breaks <- c(1, which(diff(v) != 1) + 1, length(v) + 1)
  runs <- Map(function(s, e) v[s:(e - 1)], head(breaks, -1), tail(breaks, -1))
  paste(vapply(runs, function(r) {
    if (length(r) == 1) as.character(r) else paste0(r[1], "-", r[length(r)])
  }, character(1)), collapse = ",")
}

#' Given one transcript's own exon subtable (from the catalog's exon_table),
#' return its 1-based exon numbers *within that transcript's own annotation*
#' (not remapped to any other transcript's numbering -- each transcript's
#' exon_number column already reflects its own structure).
transcript_exon_numbers <- function(exon_table, transcript_id) {
  sort(unique(exon_table$exon_number[exon_table$transcript_id == transcript_id]))
}

#' Real-app wrapper around search_confounding_proteins(): the reference mass
#' index is keyed by UniProt accession (e.g. "P16070") while proteoforms
#' built from the gene/isoform pathway are keyed by Ensembl transcript id
#' (e.g. "ENST00000428726#bare") -- exclude_id can never match across that
#' id-scheme gap, so the target's own reviewed-proteome entry would otherwise
#' come back as a "confounder" of itself. This filters any candidate whose
#' sequence is identical to the target's instead, which is robust regardless
#' of id scheme.
#'
#' @param target proteoform object
#' @param mass_index reference mass index (id, sequence, length, mass)
#' @param ... passed through to search_confounding_proteins()
#' @return same shape as search_confounding_proteins(), with true
#'   self-sequence-matches removed from $candidates
search_confounding_proteins_real <- function(target, mass_index, ...) {
  res <- search_confounding_proteins(target, mass_index, ..., exclude_id = target$id)
  if (nrow(res$candidates) > 0 && "sequence" %in% names(res$candidates)) {
    res$candidates <- res$candidates[res$candidates$sequence != target$sequence, ]
  }
  res
}

#' Build real proteoform objects for confounding-protein search hits. The
#' reference mass index carries each entry's own sequence, so (unlike the
#' earlier standalone prototype, which only had full fragment ladders for 2
#' hand-picked confounders) every confounder the live search turns up gets a
#' real, correctly computed fragment ladder here -- no special-casing.
#'
#' @param candidate_ids character vector of UniProt accessions (mass_index$id)
#' @param mass_index the reference mass index (id, sequence, length, mass)
#' @return named list, id -> proteoform object
build_confounder_proteoforms <- function(candidate_ids, mass_index) {
  pfs <- lapply(candidate_ids, function(cid) {
    row <- mass_index[mass_index$id == cid, ]
    if (nrow(row) == 0) return(NULL)
    proteoform(id = cid, sequence = row$sequence[1], provenance = "manual")
  })
  names(pfs) <- candidate_ids
  pfs[!vapply(pfs, is.null, logical(1))]
}

#' Fetch real protein sequences for a set of transcripts and build
#' proteoform() objects, skipping (with a message) any transcript with no
#' annotated translation (e.g. a non-coding or NMD transcript slipped into
#' the index).
#'
#' @param transcript_ids character vector of Ensembl transcript ids
#' @return named list, transcript_id -> proteoform object (only for
#'   transcripts that returned a real sequence)
build_proteoforms_for_transcripts <- function(transcript_ids) {
  seqs <- fetch_transcript_proteins(transcript_ids)
  valid_ids <- transcript_ids[!is.na(seqs)]
  pf_list <- lapply(valid_ids, function(tid) {
    proteoform(id = tid, sequence = seqs[[tid]], provenance = "module1_isoform_selection")
  })
  names(pf_list) <- valid_ids
  pf_list
}
