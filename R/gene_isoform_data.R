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

#' Combined confounder search: the mass-domain window search
#' (search_confounding_proteins_real(), candidates whose own INTACT MASS is
#' close to the target's) UNIONed with the m/z-domain collision search
#' (search_mz_collisions(), R/mz_collision_index.R -- candidates whose
#' charge-state peaks land on the target's peaks REGARDLESS of how far
#' their own intact mass is). These are two genuinely different real risks:
#' a candidate found only by mass could still turn out resolvable once you
#' check actual charge-state peaks; a candidate found only by m/z collision
#' can have a totally different mass yet still produce a real ambiguous
#' peak in the spectrum (e.g. a 12 kDa and a 16 kDa protein colliding at a
#' 4:3 charge-state ratio) -- a mass-window search alone would never surface
#' it. Both were already-implemented, already-tested pieces
#' (search_mz_collisions() just wasn't wired into the live confounder
#' search flow before now).
#'
#' @param target proteoform object
#' @param mass_index reference mass index (id, sequence, length, mass) --
#'   if this is a FILTERED subset (e.g. server.R's "usable top-down mass
#'   range"), an m/z collision from an id outside that filter is dropped
#'   rather than surfaced (no sequence available to build its proteoform
#'   from), even though the collision itself is real
#' @param mz_index reference m/z collision index (R/mz_collision_index.R),
#'   or NULL to skip the m/z-domain search (e.g. none built for this mode,
#'   or a digested-peptide confounder pool with no matching m/z index)
#' @param mode,r_ref,mz_ref,safety_margin passed to both underlying searches
#' @param max_candidates cap on how many candidates to return (default 30,
#'   same spirit as build_gene_isoform_catalog()'s max_isoforms). The m/z-
#'   collision search alone can surface hundreds of hits for some targets
#'   (every reference protein whose charge envelope happens to touch the
#'   target's, across its whole envelope) -- each one gets a real fragment
#'   ladder + eagerly-computed isotope pattern built downstream, so an
#'   unbounded list turns "Run analysis" from seconds into minutes. Kept
#'   candidates are the highest-priority ones: "both" (close in mass AND
#'   colliding) first, then "mass"-only, then "mz"-only sorted by
#'   n_colliding_peaks (worst offenders first) -- never an arbitrary/
#'   first-N-found subset.
#' @return list(window_da, best_fwhm_mass, best_charge_state,
#'   candidates = data.frame(id, sequence, length, mass, found_via,
#'   n_colliding_peaks) -- found_via is "mass"/"mz"/"both",
#'   n_colliding_peaks counts distinct TARGET peaks this candidate collides
#'   with (0 for mass-only hits), mz_collision_detail = raw per-peak-pair
#'   rows from search_mz_collisions(), for stats/detail display,
#'   n_candidates_before_cap = how many candidates existed before
#'   max_candidates was applied, so the UI can say "showing top 30 of 248")
search_confounding_proteins_combined <- function(target, mass_index, mz_index = NULL, mode = "denatured",
                                                  r_ref = 120000, mz_ref = 200,
                                                  safety_margin = DEFAULT_SAFETY_MARGIN,
                                                  max_candidates = 30) {
  mass_result <- search_confounding_proteins_real(target, mass_index, mode = mode,
                                                    r_ref = r_ref, mz_ref = mz_ref, safety_margin = safety_margin)
  mass_candidates <- mass_result$candidates
  mass_ids <- if (nrow(mass_candidates) > 0) mass_candidates$id else character(0)

  mz_detail <- NULL
  n_colliding_peaks <- integer(0)
  if (!is.null(mz_index)) {
    mz_detail <- tryCatch(
      search_mz_collisions(target, mz_index, mode = mode, r_ref = r_ref, mz_ref = mz_ref,
                            safety_margin = safety_margin, exclude_id = target$id),
      error = function(e) NULL
    )
    if (!is.null(mz_detail) && nrow(mz_detail) > 0) {
      # distinct TARGET peaks colliding, not raw row count -- a candidate
      # whose two different charge states both land on the SAME target peak
      # shouldn't inflate "how many of the target's own peaks are compromised"
      distinct_pairs <- unique(mz_detail[, c("id", "target_z")])
      tab <- table(distinct_pairs$id)
      n_colliding_peaks <- as.integer(tab)
      names(n_colliding_peaks) <- names(tab)
    }
  }
  mz_ids <- names(n_colliding_peaks)

  keep_cols <- c("id", "sequence", "length", "mass")
  # gene_symbol is only present on the real reference-proteome index (parsed
  # from each UniProt FASTA header's GN= field, R/reference_proteome_index.R)
  # -- optional so hand-built test fixtures without that column keep working.
  if ("gene_symbol" %in% names(mass_index)) keep_cols <- c(keep_cols, "gene_symbol")
  # mz_index may have been built from a LARGER (unfiltered) pool than
  # mass_index -- e.g. server.R restricts mass_index to the "usable
  # top-down mass range" (10-220 kDa by default), but the m/z-collision
  # search runs against the full reference proteome regardless, so it can
  # surface an id with no matching row here. Drop those defensively rather
  # than crash on match()'s resulting NA: without a sequence to build a
  # real proteoform from, that candidate can't get a fragment ladder here
  # anyway (a real, if narrow, limitation -- not silently wrong).
  extra_ids <- intersect(setdiff(mz_ids, mass_ids), mass_index$id)
  extra_rows <- if (length(extra_ids) > 0) {
    mass_index[match(extra_ids, mass_index$id), keep_cols]
  } else {
    mass_candidates[0, keep_cols]
  }
  candidates <- rbind(
    if (nrow(mass_candidates) > 0) mass_candidates[, keep_cols] else mass_candidates[0, keep_cols],
    extra_rows
  )
  if (nrow(candidates) > 0) {
    candidates$found_via <- ifelse(candidates$id %in% mass_ids & candidates$id %in% mz_ids, "both",
                              ifelse(candidates$id %in% mass_ids, "mass", "mz"))
    candidates$n_colliding_peaks <- ifelse(candidates$id %in% mz_ids,
                                            unname(n_colliding_peaks[candidates$id]), 0L)
  } else {
    candidates$found_via <- character(0)
    candidates$n_colliding_peaks <- integer(0)
  }

  n_before_cap <- nrow(candidates)
  if (n_before_cap > max_candidates) {
    via_rank <- c(both = 2L, mass = 1L, mz = 0L)[candidates$found_via]
    priority <- via_rank * 1e6 + candidates$n_colliding_peaks # via tier dominates, peaks break ties within it
    candidates <- candidates[order(-priority), ][seq_len(max_candidates), ]
    rownames(candidates) <- NULL
  }

  list(
    window_da = mass_result$window_da,
    best_fwhm_mass = mass_result$best_fwhm_mass,
    best_charge_state = mass_result$best_charge_state,
    candidates = candidates,
    mz_collision_detail = mz_detail,
    n_candidates_before_cap = n_before_cap
  )
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

#' Collapse proteoforms that share an IDENTICAL sequence down to one
#' representative each -- multiple Ensembl transcripts (e.g. splice variants
#' differing only in UTRs) commonly translate to the exact same protein, but
#' EACH still gets its own distinct protein accession (ENSP) from Ensembl;
#' confirmed directly against the live REST API (RBMS1's ENST00000269305,
#' ENST00000714359, and ENST00000905353 have three DIFFERENT ENSP ids --
#' ENSP00000269305/519626/575412 -- despite an identical, byte-for-byte
#' fetched protein sequence). So accession matching can't be used to spot
#' these groups; comparing the actual fetched sequences is what finds them.
#'
#' @param pf_list named list of proteoform objects, id -> proteoform (order
#'   determines representatives: the first id (in input order) in each
#'   sequence group wins)
#' @return list(
#'   pf_list = deduplicated named list, representative ids only,
#'   synonyms = named list, representative id -> character vector of every
#'     original id (including itself) that shared its sequence
#' )
dedupe_proteoforms_by_sequence <- function(pf_list) {
  if (length(pf_list) == 0) return(list(pf_list = pf_list, synonyms = list()))
  ids <- names(pf_list)
  seqs <- vapply(pf_list, function(p) p$sequence, character(1))
  groups <- split(ids, seqs)
  # split() orders groups by the sorted unique VALUES of `seqs`, not by
  # first-seen position -- re-order so representative/output order matches
  # the input list's order instead of alphabetical-by-sequence.
  first_pos <- vapply(groups, function(g) min(match(g, ids)), integer(1))
  groups <- groups[order(first_pos)]
  rep_ids <- vapply(groups, function(g) g[1], character(1))
  list(
    pf_list = pf_list[rep_ids],
    synonyms = setNames(groups, rep_ids)
  )
}
