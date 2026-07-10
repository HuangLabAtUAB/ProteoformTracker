# Propagating a target's PTM(s) onto its relevant-isoform set.
#
# A PTM's `site` is a position in one specific sequence. Reusing that same
# numeric index against a different isoform is only safe if the two
# sequences are perfectly colinear -- alternative splicing shifts numbering
# whenever an exon is inserted/removed upstream of the site, and can delete
# the modified residue entirely. This maps a site from the target's
# sequence to each relevant isoform's own coordinates via pairwise
# alignment, and skips (with a stated reason) any isoform where the site
# doesn't have a valid counterpart -- it never guesses.
#
# Deliberately NOT run against the confounding-protein set: PTM occupancy
# is sample/condition-specific biology, not a fixed proteome property, so
# there is no canonical "modified state" of an unrelated background protein
# to propagate onto. See known limitations in ptracker_vibe.md.

#' Global pairwise alignment (Needleman-Wunsch, linear gap penalty).
#' A minimal, dependency-free alignment sufficient for mapping positions
#' between closely related isoform sequences (differ mainly by indels from
#' alternative splicing, not distant homology) -- not a general-purpose
#' substitution-matrix aligner. Quadratic in sequence length; fine for
#' typical proteoform sizes, but expect it to be slow on very large
#' proteins (many thousands of residues).
#'
#' @param seq1,seq2 amino acid sequences
#' @param match_score,mismatch_score,gap_penalty simple linear scoring
#' @return integer vector, length nchar(seq1): position in seq2 that each
#'   position in seq1 aligns to, or NA where seq1's residue aligns to a gap
align_sequences <- function(seq1, seq2, match_score = 2, mismatch_score = -1, gap_penalty = -2) {
  s1 <- strsplit(seq1, "")[[1]]
  s2 <- strsplit(seq2, "")[[1]]
  n <- length(s1)
  m <- length(s2)

  score <- matrix(0, nrow = n + 1, ncol = m + 1)
  score[, 1] <- (0:n) * gap_penalty
  score[1, ] <- (0:m) * gap_penalty
  trace <- matrix("", nrow = n + 1, ncol = m + 1)
  trace[, 1] <- "U"
  trace[1, ] <- "L"
  trace[1, 1] <- ""

  for (i in seq_len(n)) {
    for (j in seq_len(m)) {
      diag_score <- score[i, j] + if (s1[i] == s2[j]) match_score else mismatch_score
      up_score <- score[i, j + 1] + gap_penalty
      left_score <- score[i + 1, j] + gap_penalty
      best <- max(diag_score, up_score, left_score)
      score[i + 1, j + 1] <- best
      trace[i + 1, j + 1] <- if (best == diag_score) "D" else if (best == up_score) "U" else "L"
    }
  }

  map <- rep(NA_integer_, n)
  i <- n + 1
  j <- m + 1
  while (i > 1 || j > 1) {
    dir <- trace[i, j]
    if (dir == "D") {
      map[i - 1] <- j - 1
      i <- i - 1
      j <- j - 1
    } else if (dir == "U") {
      i <- i - 1
    } else {
      j <- j - 1
    }
  }
  map
}

#' Map a single PTM's site from one sequence to another.
#'
#' Terminal sites ("N-term"/"C-term") always map trivially. Integer sites
#' are mapped via align_sequences(); the result is only marked applicable
#' if the site has a counterpart in `to_sequence` AND the residue there
#' matches the modified residue in `from_sequence` -- a position that
#' exists but has a different residue (e.g. a substitution right at an exon
#' boundary) is not treated as carrying "the same" modification.
#'
#' @param mod a ptm() object, with site relative to from_sequence
#' @param from_sequence the sequence the ptm's site is currently expressed in
#' @param to_sequence the sequence to map the site onto
#' @return list(site, applicable, reason)
map_ptm_site <- function(mod, from_sequence, to_sequence) {
  if (identical(mod$site, "N-term") || identical(mod$site, "C-term")) {
    return(list(site = mod$site, applicable = TRUE, reason = "terminal site"))
  }

  if (mod$site > nchar(from_sequence)) {
    stop("ptm site ", mod$site, " is out of range for from_sequence (length ", nchar(from_sequence), ")")
  }

  map <- align_sequences(from_sequence, to_sequence)
  mapped_site <- map[mod$site]

  if (is.na(mapped_site)) {
    return(list(
      site = NA_integer_, applicable = FALSE,
      reason = "site absent in target sequence (likely removed by alternative splicing)"
    ))
  }

  from_residue <- substr(from_sequence, mod$site, mod$site)
  to_residue <- substr(to_sequence, mapped_site, mapped_site)
  if (from_residue != to_residue) {
    return(list(
      site = mapped_site, applicable = FALSE,
      reason = sprintf("residue mismatch at mapped position (%s in source, %s here)", from_residue, to_residue)
    ))
  }

  list(site = mapped_site, applicable = TRUE, reason = "mapped")
}

#' Propagate a set of PTMs (default: the target's own) onto a list of
#' relevant-isoform proteoforms, skipping (with a reason, not silently) any
#' isoform where a site isn't a valid counterpart.
#'
#' @param target proteoform object the PTMs are currently expressed on
#' @param relevant_isoforms named or unnamed list of proteoform objects to
#'   propagate the PTMs onto
#' @param ptms list of ptm() objects to propagate (default: target$ptms)
#' @return list(proteoforms = list of new proteoform objects, one per input
#'   isoform, carrying whichever ptms mapped successfully; report =
#'   data.frame(isoform_id, ptm_name, mapped_site, applied, reason))
propagate_ptms_to_relevant_set <- function(target, relevant_isoforms, ptms = target$ptms) {
  if (!inherits(target, "proteoform")) {
    stop("propagate_ptms_to_relevant_set() requires a proteoform target")
  }
  if (length(ptms) == 0) {
    stop("no ptms to propagate -- pass `ptms` explicitly or use a target that has some")
  }

  report_rows <- list()
  new_proteoforms <- vector("list", length(relevant_isoforms))

  for (k in seq_along(relevant_isoforms)) {
    isoform <- relevant_isoforms[[k]]
    if (!inherits(isoform, "proteoform")) {
      stop("propagate_ptms_to_relevant_set() requires a list of proteoform objects")
    }

    mapped_ptms <- list()
    for (mod in ptms) {
      result <- map_ptm_site(mod, target$sequence, isoform$sequence)
      report_rows[[length(report_rows) + 1]] <- data.frame(
        isoform_id = isoform$id,
        ptm_name = if (is.na(mod$name)) "unnamed" else mod$name,
        mapped_site = if (is.na(result$site)) NA_character_ else as.character(result$site),
        applied = result$applicable,
        reason = result$reason,
        stringsAsFactors = FALSE
      )
      if (result$applicable) {
        mapped_ptms[[length(mapped_ptms) + 1]] <- ptm(
          site = result$site,
          mass_delta_mono = mod$mass_delta_mono,
          mass_delta_avg = mod$mass_delta_avg,
          name = mod$name,
          unimod_id = mod$unimod_id
        )
      }
    }

    new_proteoforms[[k]] <- proteoform(
      id = isoform$id,
      sequence = isoform$sequence,
      ptms = c(isoform$ptms, mapped_ptms),
      provenance = isoform$provenance,
      metadata = isoform$metadata
    )
  }

  list(proteoforms = new_proteoforms, report = do.call(rbind, report_rows))
}
