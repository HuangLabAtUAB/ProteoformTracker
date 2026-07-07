# Offline precomputed reference-proteome mass index.
# Confounding proteins share similar intact mass regardless of sequence/exon
# relationship to the target -- candidates for cross-gene MS1/MS2 collisions.
# Computed once offline; queried as a fast sorted-mass range lookup, never
# recomputed per user query.

#' Minimal FASTA reader (id = text up to first whitespace on the header line).
#' Avoids a Biostrings dependency for this simple, size-unbounded use case.
read_fasta <- function(path) {
  lines <- readLines(path, warn = FALSE)
  header_idx <- grep("^>", lines)
  if (length(header_idx) == 0) stop("no FASTA headers found in ", path)
  ids <- sub("^>(\\S+).*$", "\\1", lines[header_idx])
  end_idx <- c(header_idx[-1] - 1, length(lines))
  sequences <- mapply(
    function(start, end) paste(lines[(start + 1):end], collapse = ""),
    header_idx, end_idx
  )
  names(sequences) <- ids
  sequences
}

#' Build the offline reference-proteome mass index.
#'
#' One-time (per proteome release) computation: translates every sequence's
#' theoretical intact mass (bare sequence, no PTMs -- PTMs are pair-specific
#' and applied downstream, not part of the background index), applies the
#' N-terminal Met excision heuristic, and stores a mass-sorted table for fast
#' range queries.
#'
#' @param fasta_path path to reference proteome FASTA (e.g. UniProt reference
#'   proteome)
#' @param output_path where to save the index (.rds)
#' @param average if TRUE, index average mass; otherwise monoisotopic
#' @return the index data.frame (invisibly), also written to output_path
build_reference_mass_index <- function(fasta_path, output_path, average = FALSE) {
  init_mass_calculation_engine()
  sequences <- read_fasta(fasta_path)

  valid <- grepl(paste0("^[", paste(STANDARD_AA, collapse = ""), "]+$"), toupper(sequences))
  if (any(!valid)) {
    message(sum(!valid), " of ", length(sequences), " entries skipped (non-standard residues)")
  }
  sequences <- sequences[valid]

  mature_sequences <- vapply(sequences, function(s) {
    predict_nterminal_met_excision(s)$mature_sequence
  }, character(1))

  masses <- vapply(mature_sequences, function(s) sequence_mass(s, average = average), numeric(1))

  index <- data.frame(
    id = names(sequences),
    sequence = unname(mature_sequences),
    length = nchar(mature_sequences),
    mass = unname(masses),
    stringsAsFactors = FALSE
  )
  index <- index[order(index$mass), ]
  rownames(index) <- NULL

  saveRDS(index, output_path)
  invisible(index)
}

#' Load a previously built reference-proteome mass index.
load_reference_mass_index <- function(path) {
  readRDS(path)
}

#' Fast range lookup: reference-proteome entries whose theoretical intact
#' mass falls within +/- window_da of target_mass. O(log n) via binary
#' search against the mass-sorted index (findInterval), not a per-query
#' full-table scan.
#'
#' @param index a reference mass index (from build_/load_reference_mass_index)
#' @param target_mass mass to search around (Da)
#' @param window_da half-width of the search window (Da)
#' @param exclude_id optional id to exclude from results (the target's own
#'   entry, if it is itself part of the reference proteome)
query_confounding_proteins <- function(index, target_mass, window_da, exclude_id = NULL) {
  lo <- target_mass - window_da
  hi <- target_mass + window_da
  lo_idx <- findInterval(lo, index$mass) + 1L
  hi_idx <- findInterval(hi, index$mass)
  if (lo_idx > hi_idx) {
    return(index[0, ])
  }
  hits <- index[lo_idx:hi_idx, ]
  if (!is.null(exclude_id)) {
    hits <- hits[hits$id != exclude_id, ]
  }
  hits
}
