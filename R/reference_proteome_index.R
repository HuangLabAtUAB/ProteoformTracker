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

#' Extract the bare UniProt accession from a "sp|ACCESSION|ENTRY_NAME" or
#' "tr|ACCESSION|ENTRY_NAME" style FASTA header token. Falls back to the
#' input unchanged if it doesn't match that pattern (e.g. non-UniProt FASTA).
parse_uniprot_accession <- function(id) {
  ifelse(grepl("^(sp|tr)\\|", id), sub("^(sp|tr)\\|([^|]+)\\|.*$", "\\2", id), id)
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
#' @param script_path path to python/ptracker_mass.py, relative to the
#'   current working directory (passed through to sequence_masses_batch())
#' @return the index data.frame (invisibly), also written to output_path
build_reference_mass_index <- function(fasta_path, output_path, average = FALSE,
                                        script_path = "python/ptracker_mass.py") {
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

  masses <- sequence_masses_batch(mature_sequences, average = average, script_path = script_path)

  index <- data.frame(
    id = parse_uniprot_accession(names(sequences)),
    sequence = unname(mature_sequences),
    length = nchar(mature_sequences),
    mass = masses,
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

# Fixed precomputed mass range for the middle-down confounder pool below --
# independent of whatever mass window the user currently has the app's own
# "Min/max peptide mass" sliders set to. The per-target search window
# query_confounding_proteins() actually uses is tiny (a fraction of a Da to a
# few Da, from the resolving-power model), so the pool just needs broad
# enough coverage that any mass window a user would plausibly pick for
# middle-down sits inside it -- it does not need to track the live sliders.
DIGESTED_POOL_MASS_MIN_DA <- 1000
DIGESTED_POOL_MASS_MAX_DA <- 15000

#' Build the enzyme-digested confounder pool for one protease: every
#' proteome-scale digest fragment (any missed-cleavage count, any parent
#' protein) whose mass falls in DIGESTED_POOL_MASS_MIN/MAX_DA, in the same
#' (id, sequence, length, mass) shape as build_reference_mass_index()'s
#' output -- so query_confounding_proteins()/search_confounding_proteins()
#' work against it completely unchanged.
#'
#' Real confounders for a middle-down target peptide are OTHER proteins'
#' digest peptides, not other intact proteins -- an intact 4 kDa protein and
#' a 4 kDa peptide cut out of the middle of a 60 kDa protein are both real
#' possible confounders, but only the digest-pool captures the second case.
#'
#' This is a genuinely heavy, proteome-scale computation (~20k proteins x up
#' to a few hundred candidate windows each): pure-R cumulative-sum math
#' (digest_sequence_for_pool(), R/digestion.R) rather than a per-candidate
#' pyteomics call is what keeps this in the tens-of-seconds range rather
#' than many minutes -- verified to match pyteomics' own per-sequence
#' calculate_mass() to ~1e-9 Da. Meant to be run once per enzyme and cached
#' to disk (see get_digested_reference_pool()), not recomputed per request.
#'
#' @param fasta_path path to reference proteome FASTA
#' @param enzyme one of names(ENZYME_SPECS)
#' @param output_path where to save the index (.rds)
#' @param mass_min_da,mass_max_da override DIGESTED_POOL_MASS_MIN/MAX_DA
#' @return the index data.frame (invisibly), also written to output_path
build_digested_reference_pool <- function(fasta_path, enzyme, output_path,
                                           mass_min_da = DIGESTED_POOL_MASS_MIN_DA,
                                           mass_max_da = DIGESTED_POOL_MASS_MAX_DA) {
  sequences <- read_fasta(fasta_path)
  valid <- grepl(paste0("^[", paste(STANDARD_AA, collapse = ""), "]+$"), toupper(sequences))
  sequences <- sequences[valid]
  ids <- parse_uniprot_accession(names(sequences))

  rows <- vector("list", length(sequences))
  for (i in seq_along(sequences)) {
    rows[[i]] <- digest_sequence_for_pool(ids[i], sequences[[i]], enzyme, mass_min_da, mass_max_da)
  }
  index <- do.call(rbind, rows)
  if (is.null(index) || nrow(index) == 0) {
    stop("no digest fragments survived the mass window for enzyme ", enzyme)
  }
  index <- index[order(index$mass), ]
  rownames(index) <- NULL

  saveRDS(index, output_path)
  invisible(index)
}

# In-memory single-slot cache: holds at most one enzyme's digested pool at a
# time (each is large -- tens of millions of Da-sorted rows including full
# peptide sequences, order ~200+ MB), so keeping all enzymes loaded
# simultaneously would be wasteful; switching enzymes just reloads from the
# on-disk cache (a plain readRDS(), not a rebuild) instead.
.digested_pool_cache <- local({
  cached_enzyme <- NULL
  cached_index <- NULL
  function(enzyme, fasta_path, cache_dir, on_build_start = NULL) {
    if (!is.null(cached_enzyme) && identical(cached_enzyme, enzyme)) {
      return(cached_index)
    }
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
    cache_path <- file.path(cache_dir, sprintf("reference_digest_pool_%s.rds", gsub("[^A-Za-z0-9]", "", enzyme)))
    index <- if (file.exists(cache_path)) {
      readRDS(cache_path)
    } else {
      if (!is.null(on_build_start)) on_build_start()
      build_digested_reference_pool(fasta_path, enzyme, cache_path)
    }
    cached_enzyme <<- enzyme
    cached_index <<- index
    index
  }
})

#' Lazy-loaded, disk-cached enzyme-digested confounder pool: builds it once
#' per enzyme (a real, ~20-40 second computation over the whole reference
#' proteome) and caches the result to `cache_dir` so every later app restart
#' or enzyme reselection just reads the cached file back rather than
#' rebuilding. `on_build_start` lets callers (e.g. a Shiny withProgress()
#' block) surface that first-time cost to the user instead of it looking
#' like a hang.
#'
#' @param enzyme one of names(ENZYME_SPECS)
#' @param fasta_path reference proteome FASTA to digest if no cache exists yet
#' @param cache_dir directory for the per-enzyme .rds cache files
#' @param on_build_start optional zero-arg callback invoked right before a
#'   cache-miss build starts (e.g. to show a progress message)
#' @return data.frame(id, sequence, length, mass), same shape as
#'   load_reference_mass_index()'s output
get_digested_reference_pool <- function(enzyme, fasta_path = "data/UP000005640_9606_reviewed_canonical.fasta",
                                         cache_dir = "data", on_build_start = NULL) {
  .digested_pool_cache(enzyme, fasta_path, cache_dir, on_build_start)
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
