# Live per-transcript protein sequence fetch via the Ensembl REST API.
#
# The genome-wide exon index (R/reference_exon_index.R) gives us exon
# structure/coordinates for every transcript without a live call, but not
# translated protein sequence -- building that would mean shipping (or
# indexing into) a ~3GB genome FASTA. For arbitrary-gene support in the app,
# a small per-transcript REST call is simpler and faster than requiring the
# full genome locally: Ensembl's /sequence/id/:id?type=protein endpoint
# already applies the correct CDS/frame/splicing, so no local translation
# logic is needed. Results are cached to disk (keyed by transcript id) since
# a given transcript's sequence never changes within an Ensembl release.

ENSEMBL_REST_BASE <- "https://rest.ensembl.org"
ENSEMBL_PROTEIN_CACHE_DIR <- "data/ensembl_protein_cache"

#' GET a URL with retries for transient failures. Ensembl's public REST API
#' intermittently returns 5xx or times out under load (observed directly:
#' the same request failed once with a 500 and succeeded moments later with
#' a 200 and a real sequence) -- retrying a few times with backoff turns
#' those transient hiccups into a successful fetch instead of a permanent
#' "unavailable" result. A 400 is NOT retried: Ensembl uses it to mean "this
#' id genuinely has no translation/entry", a real negative answer, not a
#' server hiccup.
#'
#' @param url request URL
#' @param content_type value for the Content-Type header
#' @param max_attempts total attempts including the first
#' @return httr response object (its status may still be non-2xx if every
#'   attempt failed or the failure was a genuine 400)
.ensembl_get_with_retry <- function(url, content_type, max_attempts = 4) {
  resp <- NULL
  for (attempt in seq_len(max_attempts)) {
    resp <- tryCatch(
      httr::GET(url, httr::add_headers(`Content-Type` = content_type), httr::timeout(15)),
      error = function(e) e
    )
    transient_failure <- inherits(resp, "error") ||
      (httr::status_code(resp) >= 500 && httr::status_code(resp) < 600)
    if (!transient_failure) break
    if (attempt < max_attempts) Sys.sleep(0.5 * 2^(attempt - 1))
  }
  if (inherits(resp, "error")) stop(conditionMessage(resp))
  resp
}

#' Fetch a transcript's translated protein sequence from Ensembl REST,
#' caching the result to disk so repeat lookups (including across app
#' restarts) don't re-hit the network.
#'
#' @param transcript_id Ensembl transcript id, e.g. "ENST00000428726"
#' @param cache_dir directory for on-disk cache (one .txt file per transcript)
#' @return character(1) protein sequence (no header), or NA_character_ with
#'   a warning if the transcript genuinely has no annotated translation
#'   (e.g. a non-coding transcript) or the fetch failed after retries
fetch_transcript_protein <- function(transcript_id, cache_dir = ENSEMBL_PROTEIN_CACHE_DIR) {
  if (!is.character(transcript_id) || length(transcript_id) != 1 || !nzchar(transcript_id)) {
    stop("fetch_transcript_protein() requires a single transcript id string")
  }
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
  cache_path <- file.path(cache_dir, paste0(transcript_id, ".txt"))
  if (file.exists(cache_path)) {
    cached <- trimws(readLines(cache_path, warn = FALSE))
    if (length(cached) == 1 && nzchar(cached)) return(cached)
  }

  url <- sprintf("%s/sequence/id/%s?type=protein", ENSEMBL_REST_BASE, transcript_id)
  resp <- .ensembl_get_with_retry(url, "text/x-fasta")
  if (httr::status_code(resp) == 400) {
    warning("Ensembl REST: no protein translation for ", transcript_id,
            " (likely a non-coding transcript)")
    return(NA_character_)
  }
  httr::stop_for_status(resp, task = paste("fetch protein sequence for", transcript_id))

  body <- httr::content(resp, as = "text", encoding = "UTF-8")
  lines <- strsplit(body, "\n")[[1]]
  seq_lines <- lines[!startsWith(lines, ">")]
  protein_seq <- paste(trimws(seq_lines), collapse = "")
  if (!nzchar(protein_seq)) {
    stop("Ensembl REST returned an empty sequence for ", transcript_id)
  }

  writeLines(protein_seq, cache_path)
  protein_seq
}

ENSEMBL_EXON_CACHE_DIR <- "data/ensembl_exon_cache"

#' Fetch a transcript's FULL exon structure (every exon, including UTR-only
#' ones) from Ensembl REST, caching to disk like fetch_transcript_protein().
#'
#' Deliberately a live per-transcript call rather than reference_exon_index:
#' that index is built from CDS-only GTF rows (see R/reference_exon_index.R),
#' so it omits non-coding exons (e.g. TP53's own exon 1). Comparing a novel
#' sequence's full minimap2 alignment footprint against a CDS-only exon list
#' would show spurious "gaps" for real UTR exons that were simply never in
#' the index -- this fetches the true, complete exon list instead so the
#' comparison is apples-to-apples.
#'
#' @param transcript_id Ensembl transcript id, e.g. "ENST00000714408"
#' @param cache_dir directory for on-disk cache (one .rds file per transcript)
#' @return list(seqname, strand, exons = data.frame(start, end)) ordered by
#'   genomic start, or NULL if the transcript isn't found
fetch_transcript_exons <- function(transcript_id, cache_dir = ENSEMBL_EXON_CACHE_DIR) {
  if (!is.character(transcript_id) || length(transcript_id) != 1 || !nzchar(transcript_id)) {
    stop("fetch_transcript_exons() requires a single transcript id string")
  }
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
  cache_path <- file.path(cache_dir, paste0(transcript_id, ".rds"))
  if (file.exists(cache_path)) return(readRDS(cache_path))

  url <- sprintf("%s/lookup/id/%s?expand=1", ENSEMBL_REST_BASE, transcript_id)
  resp <- .ensembl_get_with_retry(url, "application/json")
  if (httr::status_code(resp) == 400) {
    warning("Ensembl REST: transcript not found: ", transcript_id)
    return(NULL)
  }
  httr::stop_for_status(resp, task = paste("fetch exon structure for", transcript_id))

  parsed <- jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"), simplifyVector = FALSE)
  exon_list <- parsed$Exon
  if (is.null(exon_list) || length(exon_list) == 0) {
    warning("Ensembl REST returned no exons for ", transcript_id)
    return(NULL)
  }
  exons <- data.frame(
    start = vapply(exon_list, function(e) e$start, numeric(1)),
    end = vapply(exon_list, function(e) e$end, numeric(1))
  )
  exons <- exons[order(exons$start), ]
  result <- list(
    seqname = as.character(parsed$seq_region_name),
    strand = if (identical(parsed$strand, 1L) || identical(parsed$strand, 1)) "+" else "-",
    exons = exons
  )
  saveRDS(result, cache_path)
  result
}

#' Batch convenience wrapper: fetch protein sequences for several transcript
#' ids at once, returning a named list (skips/NA entries for any that fail).
#'
#' @param transcript_ids character vector of Ensembl transcript ids
#' @param cache_dir passed through to fetch_transcript_protein()
#' @return named list, transcript_id -> protein sequence (or NA_character_)
fetch_transcript_proteins <- function(transcript_ids, cache_dir = ENSEMBL_PROTEIN_CACHE_DIR) {
  seqs <- lapply(transcript_ids, function(tid) {
    tryCatch(
      fetch_transcript_protein(tid, cache_dir = cache_dir),
      error = function(e) {
        warning("Failed to fetch protein sequence for ", tid, ": ", conditionMessage(e))
        NA_character_
      }
    )
  })
  names(seqs) <- transcript_ids
  seqs
}
