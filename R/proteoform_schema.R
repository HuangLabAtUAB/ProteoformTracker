# Proteoform object: the shared internal schema all three input modules
# (isoform selection, long-read FASTA->ORF, rMATS event parsing) must resolve
# into before any downstream code (mass calc, charge envelope, resolvability
# scoring, fragment ladder) runs.

VALID_PROVENANCE <- c(
  "module1_isoform_selection",
  "module2_longread_orf",
  "module3_rmats_event",
  "module_middledown_digest",
  "manual"
)

STANDARD_AA <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]

#' Construct a single PTM entry (Unimod-style mass delta)
#'
#' @param site 1-based position in the mature sequence, or "N-term"/"C-term"
#' @param mass_delta_mono monoisotopic mass delta in Da
#' @param mass_delta_avg average mass delta in Da (defaults to mono if omitted)
#' @param name human-readable modification name (e.g. "Phospho")
#' @param unimod_id Unimod accession, e.g. "UNIMOD:21" (optional)
ptm <- function(site, mass_delta_mono, mass_delta_avg = mass_delta_mono,
                 name = NA_character_, unimod_id = NA_character_) {
  is_terminal_site <- identical(site, "N-term") || identical(site, "C-term")
  if (!is_terminal_site) {
    if (!is.numeric(site) || length(site) != 1 || site < 1 || site != as.integer(site)) {
      stop("ptm() site must be a positive integer position, or 'N-term'/'C-term'")
    }
    site <- as.integer(site)
  }
  if (!is.numeric(mass_delta_mono) || length(mass_delta_mono) != 1) {
    stop("ptm() mass_delta_mono must be a single numeric value")
  }
  structure(
    list(
      site = site,
      mass_delta_mono = mass_delta_mono,
      mass_delta_avg = mass_delta_avg,
      name = name,
      unimod_id = unimod_id
    ),
    class = "ptm"
  )
}

#' Construct a Proteoform object
#'
#' @param id identifier for this proteoform (e.g. gene_transcript_isoform tag)
#' @param sequence mature amino acid sequence: after N-terminal Met excision
#'   and signal peptide/propeptide removal (uppercase, standard 20 AA letters)
#' @param ptms list of ptm() objects; empty by default
#' @param provenance which input module produced this object; one of
#'   VALID_PROVENANCE
#' @param metadata optional free-form list (gene id, transcript id, notes on
#'   what processing was applied to reach `sequence`, etc.)
proteoform <- function(id, sequence, ptms = list(), provenance, metadata = list()) {
  obj <- structure(
    list(
      id = id,
      sequence = sequence,
      ptms = ptms,
      provenance = provenance,
      metadata = metadata
    ),
    class = "proteoform"
  )
  validate_proteoform(obj)
  obj
}

validate_proteoform <- function(x) {
  if (!inherits(x, "proteoform")) stop("not a proteoform object")
  if (!is.character(x$id) || length(x$id) != 1 || is.na(x$id) || nchar(x$id) == 0) {
    stop("proteoform id must be a single non-empty string")
  }
  if (!is.character(x$sequence) || length(x$sequence) != 1 || is.na(x$sequence) || nchar(x$sequence) == 0) {
    stop("proteoform sequence must be a single non-empty string")
  }
  seq_chars <- strsplit(toupper(x$sequence), "")[[1]]
  bad_chars <- setdiff(unique(seq_chars), STANDARD_AA)
  if (length(bad_chars) > 0) {
    stop(
      "proteoform sequence contains non-standard amino acid letters: ",
      paste(bad_chars, collapse = ", ")
    )
  }
  if (!is.list(x$ptms) || !all(vapply(x$ptms, inherits, logical(1), "ptm"))) {
    stop("proteoform ptms must be a list of ptm() objects")
  }
  if (!is.character(x$provenance) || length(x$provenance) != 1 || !(x$provenance %in% VALID_PROVENANCE)) {
    stop(
      "proteoform provenance must be one of: ",
      paste(VALID_PROVENANCE, collapse = ", ")
    )
  }
  invisible(TRUE)
}

print.proteoform <- function(x, ...) {
  cat(sprintf(
    "<proteoform %s> %d aa, %d PTM(s), provenance=%s\n",
    x$id, nchar(x$sequence), length(x$ptms), x$provenance
  ))
  invisible(x)
}

#' Predict whether N-terminal Met is excised, from the penultimate residue
#' (position 2 of the *unprocessed* ORF translation).
#'
#' Heuristic (methionine aminopeptidase specificity / "N-end rule"):
#' Met is excised when the penultimate residue has a small side chain
#' (A, C, G, P, S, T, V); retained otherwise. Context-dependent, not
#' guaranteed correct for every sequence -- flag as a heuristic in the UI.
#'
#' @param raw_sequence unprocessed ORF translation, Met-initiated
#' @return list(mature_sequence, met_excised = logical)
predict_nterminal_met_excision <- function(raw_sequence) {
  chars <- strsplit(toupper(raw_sequence), "")[[1]]
  if (length(chars) < 2 || chars[1] != "M") {
    return(list(mature_sequence = raw_sequence, met_excised = FALSE))
  }
  small_penultimate <- c("A", "C", "G", "P", "S", "T", "V")
  if (chars[2] %in% small_penultimate) {
    list(mature_sequence = paste(chars[-1], collapse = ""), met_excised = TRUE)
  } else {
    list(mature_sequence = raw_sequence, met_excised = FALSE)
  }
}

#' Build a Proteoform from a raw (unprocessed) ORF translation, applying
#' N-terminal Met excision and, optionally, a user-supplied mature-sequence
#' start position (signal peptide/propeptide cleavage site). No signal
#' peptide *predictor* is implemented -- the cleavage site must be supplied
#' by the caller when known/suspected (same policy as PTMs).
#'
#' @param id proteoform identifier
#' @param raw_sequence unprocessed ORF translation (Met-initiated)
#' @param mature_start 1-based position (in raw_sequence) where the mature
#'   chain begins, e.g. after signal peptide/propeptide removal. If NULL,
#'   only N-terminal Met excision is applied.
#' @param ptms list of ptm() objects
#' @param provenance which input module produced this object
#' @param metadata optional free-form list
build_proteoform <- function(id, raw_sequence, mature_start = NULL,
                              ptms = list(), provenance, metadata = list()) {
  processing_notes <- character(0)

  if (!is.null(mature_start)) {
    if (!is.numeric(mature_start) || mature_start < 1 || mature_start > nchar(raw_sequence)) {
      stop("mature_start must be a valid 1-based position within raw_sequence")
    }
    mature_sequence <- substring(raw_sequence, mature_start)
    processing_notes <- c(processing_notes, sprintf(
      "signal peptide/propeptide removed: mature chain starts at raw position %d (user-supplied)",
      mature_start
    ))
    met_result <- predict_nterminal_met_excision(mature_sequence)
    mature_sequence <- met_result$mature_sequence
  } else {
    met_result <- predict_nterminal_met_excision(raw_sequence)
    mature_sequence <- met_result$mature_sequence
  }

  processing_notes <- c(processing_notes, sprintf(
    "N-terminal Met excision: %s (heuristic, penultimate-residue rule)",
    if (met_result$met_excised) "applied" else "not applied"
  ))

  metadata$processing_notes <- processing_notes
  metadata$raw_sequence <- raw_sequence

  proteoform(
    id = id,
    sequence = mature_sequence,
    ptms = ptms,
    provenance = provenance,
    metadata = metadata
  )
}
