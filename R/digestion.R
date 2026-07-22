# Middle-down in silico digestion: simulate the large, partially-cleaved
# peptides a limited-time protease digestion produces (as opposed to a full
# bottom-up digest), by enumerating every fragment between any two cleavage
# sites (i.e. every missed-cleavage count at once) and keeping only the ones
# whose mass lands in the user's target window. No cleavage-kinetics model is
# attempted -- real per-site cleavage rates aren't data this app has, and the
# mass window is what actually stands in for "limited digestion" here: a
# fragment requiring many missed cleavages only survives if it happens to
# land in-window anyway, so implausible sky-high-missed-cleavage fragments
# are excluded by the same mass filter as everything else, not by a separate
# arbitrary cutoff.

#' Sanitize a candidate/proteoform id for use as (part of) an HTML element
#' id -- e.g. "ENST00000714408#bare_pep321_351" contains "#", which breaks
#' jQuery/Shiny's id-based lookup for dynamically created inputs (checkboxes
#' named "pep_chk_<that id>") since "#" is the CSS/jQuery id-selector prefix
#' character. Server-side R code should keep using the real (unsanitized)
#' candidate id for all data lookups; only the actual widget id needs this.
sanitize_html_id <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)

# Each cleavage-site rule is expressed in the same "cleavage_position" 1..n-1
# convention already used throughout fragment_ladder.R / fragmentation_propensity.R:
# position p means the bond BETWEEN residue p and residue p+1.
#   kind = "after"  cuts C-terminal to a residue in `residues` (site = residue's own position)
#   kind = "before" cuts N-terminal to a residue in `residues` (site = residue's position - 1)
#   kind = "dibasic" cuts between two consecutive basic residues (OmpT)
# `block_next`/`block_prev` encode known inhibited contexts (e.g. trypsin/Lys-C
# do not cleave Lys-Pro; Glu-C is likewise blocked before Pro).
ENZYME_SPECS <- list(
  "OmpT"  = list(kind = "dibasic"),
  "Lys-C" = list(kind = "after", residues = c("K"), block_next = c("P")),
  "Lys-N" = list(kind = "before", residues = c("K")),
  "Glu-C" = list(kind = "after", residues = c("E"), block_next = c("P")),
  "Asp-N" = list(kind = "before", residues = c("D"))
)

#' Enumerate every cleavage-site position (1..n-1) a given enzyme would cut
#' at in a sequence, per its ENZYME_SPECS rule.
#'
#' @param sequence amino acid sequence
#' @param enzyme one of names(ENZYME_SPECS)
#' @return sorted integer vector of cleavage positions (possibly empty)
find_cleavage_sites <- function(sequence, enzyme) {
  spec <- ENZYME_SPECS[[enzyme]]
  if (is.null(spec)) stop("unknown enzyme: ", enzyme, " -- must be one of: ", paste(names(ENZYME_SPECS), collapse = ", "))
  chars <- strsplit(toupper(sequence), "")[[1]]
  n <- length(chars)
  if (n < 2) return(integer(0))

  if (identical(spec$kind, "dibasic")) {
    basic <- c("K", "R")
    return(which(chars[1:(n - 1)] %in% basic & chars[2:n] %in% basic))
  }
  if (identical(spec$kind, "after")) {
    pos <- which(chars[1:(n - 1)] %in% spec$residues)
    if (!is.null(spec$block_next) && length(pos) > 0) {
      pos <- pos[!(chars[pos + 1] %in% spec$block_next)]
    }
    return(sort(unique(pos)))
  }
  if (identical(spec$kind, "before")) {
    idx <- which(chars %in% spec$residues)
    idx <- idx[idx > 1] # nothing precedes position 1, so no bond to cut there
    return(sort(unique(idx - 1L)))
  }
  stop("unrecognized ENZYME_SPECS kind for ", enzyme)
}

# Monoisotopic residue masses (pyteomics.mass.calculate_mass convention: bare
# peptide mass = sum of residue masses + one water), hardcoded so bulk
# reference-proteome digestion (see build_digested_reference_pool() in
# R/reference_proteome_index.R) can compute fragment masses via a single
# cumulative sum per protein instead of one reticulate round-trip per
# candidate fragment -- the whole point being that a proteome-scale digest
# (tens of thousands of proteins x hundreds of candidate fragments each) is
# only tractable to build if mass computation itself is pure-R arithmetic.
# Values spot-checked directly against pyteomics.mass.calculate_mass() and
# match to ~1e-9 Da; this table has NO PTM support (background pool masses
# are always bare, same as build_reference_mass_index()'s existing choice).
RESIDUE_MASS_MONO <- c(
  A = 71.0371137847, C = 103.0091847847, D = 115.0269430238, E = 129.0425930880,
  F = 147.0684139130, G = 57.0214637206, H = 137.0589118585, I = 113.0840639771,
  K = 128.0949630140, L = 113.0840639771, M = 131.0404849130, N = 114.0429274411,
  P = 97.0527638489, Q = 128.0585775053, R = 156.1011110236, S = 87.0320284043,
  T = 101.0476784684, V = 99.0684139130, W = 186.0793129499, Y = 163.0633285325
)
WATER_MASS_MONO <- 18.0105646863

#' Digest a single bare sequence into every mass-in-window fragment (any
#' missed-cleavage count), for building the proteome-scale confounder pool.
#' Pure-R cumulative-sum arithmetic (no PTMs, no MS1/MS2 scoring) -- this is
#' the bulk-background-pool counterpart to digest_proteoform(), which does
#' the fuller per-candidate scoring but is too slow (one reticulate call
#' worth of overhead per candidate) to run at proteome scale.
#'
#' @param id sequence identifier (e.g. UniProt accession)
#' @param sequence amino acid sequence (standard 20 AA letters only; returns
#'   NULL if it contains anything else, e.g. an ambiguous code)
#' @param enzyme one of names(ENZYME_SPECS)
#' @param mass_min_da,mass_max_da target mass window, Da
#' @return data.frame(id, sequence, length, mass) -- one row per surviving
#'   candidate fragment, same column shape as the reference_mass_index this
#'   is meant to stand in for (see R/reference_proteome_index.R) -- or NULL
#'   if nothing survives (no cleavage sites, or nothing lands in-window)
digest_sequence_for_pool <- function(id, sequence, enzyme, mass_min_da, mass_max_da) {
  n <- nchar(sequence)
  chars <- strsplit(sequence, "")[[1]]
  resm <- RESIDUE_MASS_MONO[chars]
  if (anyNA(resm)) return(NULL) # non-standard residue somewhere in this sequence
  sites <- find_cleavage_sites(sequence, enzyme)
  boundaries <- c(0L, sites, n)
  k <- length(boundaries)
  if (k < 2) return(NULL)

  cum <- c(0, cumsum(resm)) # cum[i+1] = sum of residue masses for the first i residues
  idx_i <- rep(seq_len(k - 1), times = rev(seq_len(k - 1)))
  idx_j <- unlist(lapply(seq_len(k - 1), function(i) seq(i + 1, k)))
  starts <- boundaries[idx_i] + 1L
  ends <- boundaries[idx_j]
  mass <- (cum[ends + 1] - cum[starts]) + WATER_MASS_MONO

  keep <- mass >= mass_min_da & mass <= mass_max_da
  if (!any(keep)) return(NULL)
  starts <- starts[keep]; ends <- ends[keep]; mass <- mass[keep]

  data.frame(
    id = sprintf("%s_pep%d_%d", id, starts, ends),
    sequence = substring(sequence, starts, ends),
    length = ends - starts + 1L,
    mass = mass,
    stringsAsFactors = FALSE
  )
}

#' Builds the single "intact protein" fallback candidate for a proteoform
#' whose digest has no real fragment (at any missed-cleavage count) landing
#' in the mass window -- see digest_proteoform()'s doc comment for why this
#' exists. Bypasses the mass window entirely (that's the point: guarantee
#' this proteoform is representable one way or another) and reuses the
#' parent's own PTMs unchanged (no site remapping needed since start=1).
.intact_fallback_candidate <- function(pf, parent_id, mode, r_ref, mz_ref, average) {
  n <- nchar(pf$sequence)
  pep_id <- sprintf("%s_intact", parent_id)
  pep_pf <- proteoform(
    id = pep_id, sequence = pf$sequence, ptms = pf$ptms,
    provenance = "module_middledown_digest",
    metadata = list(parent_id = parent_id, start = 1L, end = n, missed_cleavages = NA_integer_, is_intact = TRUE)
  )
  mass <- proteoform_mass(pf, average = average)$mass
  fwhm_info <- .best_case_fwhm(pep_pf, mode, average, r_ref, mz_ref)
  ms2_prop <- if (n >= 2) mean(fragmentation_propensity(pep_pf, mode = mode)$propensity_score) else NA_real_
  candidates <- data.frame(
    id = pep_id, parent_id = parent_id, start = 1L, end = n, length = n,
    missed_cleavages = NA_integer_, mass = mass,
    ms1_fwhm_da = fwhm_info$best_fwhm, ms1_best_z = fwhm_info$best_z,
    ms2_avg_propensity = ms2_prop, ptm_sites_covered = length(pf$ptms),
    is_intact = TRUE, stringsAsFactors = FALSE
  )
  list(candidates = candidates, peptides = setNames(list(pep_pf), pep_id))
}

#' Simulate limited-digestion middle-down peptides from a full-length
#' proteoform: every fragment between any two cleavage sites (any number of
#' missed cleavages), filtered to a target mass window (applied to the
#' PTM-inclusive mass, since middle-down's point is usually to preserve PTM
#' co-occurrence). For each surviving candidate, also computes the two-tier
#' selection metrics discussed with the user:
#'   Tier 1 (feasibility) -- missed_cleavages (fewer = more likely generated),
#'     ms1_fwhm_da/ms1_best_z (this peptide's own best-case resolvability,
#'     reusing the same charge-envelope/resolving-power model used for
#'     full proteoforms), ms2_avg_propensity (mean HCD/CID fragmentation
#'     propensity across its own bonds, reusing the existing propensity model)
#'   Tier 2 (informativeness) -- ptm_sites_covered (how many of the parent's
#'     annotated PTM sites this peptide fully contains)
#' Non-redundant coverage (the other Tier 2 point) is inherently relative to
#' whichever candidates the user has *already* selected, so it isn't a
#' per-row property computed here -- see digestion_coverage_summary().
#'
#' @param pf parent proteoform object (full-length, PTMs already applied)
#' @param enzyme one of names(ENZYME_SPECS)
#' @param mass_min_da,mass_max_da target mass window, Da
#' @param mode "denatured" or "native" (feeds both MS1 envelope + MS2 propensity)
#' @param r_ref,mz_ref Orbitrap resolving-power settings
#' @param average use average mass instead of monoisotopic (should match
#'   whatever convention the rest of the app's masses use -- default mono)
#' @param parent_id id to key candidates/metadata off of -- defaults to the
#'   proteoform's own $id, but callers that key their proteoform lists by a
#'   DIFFERENT id (e.g. derived()$rows' "<transcript>#bare"/"<transcript>#N"
#'   keys, which differ from the bare proteoform object's own $id field)
#'   should pass that key explicitly so candidate ids/parent_id stay
#'   consistent with however the caller looks its own proteoforms up.
#' @return list(candidates = data.frame or NULL, peptides = named list of
#'   proteoform objects, one per surviving candidate id). candidates$is_intact
#'   is TRUE for the intact-protein fallback row (see below), FALSE for a
#'   real digest fragment.
digest_proteoform <- function(pf, enzyme, mass_min_da, mass_max_da,
                               mode = c("denatured", "native"),
                               r_ref = 120000, mz_ref = 200, average = FALSE,
                               parent_id = pf$id) {
  mode <- match.arg(mode)
  if (!inherits(pf, "proteoform")) stop("digest_proteoform() requires a proteoform object")
  sequence <- pf$sequence
  n <- nchar(sequence)
  sites <- find_cleavage_sites(sequence, enzyme)
  boundaries <- c(0L, sites, n)
  k <- length(boundaries)
  # A sparse-cleavage-site enzyme (e.g. OmpT's rare dibasic sites) on a
  # protein where none of the resulting fragments happen to land in the
  # mass window previously meant this protein contributed NOTHING to
  # "Run analysis" even when explicitly checked in the Result proteoform
  # table -- confirmed as a real point of user confusion, since a checked-
  # but-silently-excluded proteoform looks identical to a bug. Falling back
  # to the intact protein as its own single candidate guarantees every
  # checked proteoform is representable in the comparison one way or
  # another; it's clearly flagged (is_intact) so it reads as "no digest
  # fragment qualified" rather than a real middle-down peptide.
  if (k < 2) return(.intact_fallback_candidate(pf, parent_id, mode, r_ref, mz_ref, average))

  # Every (i, j) pair of boundary indices, i < j, is one candidate fragment;
  # missed_cleavages counts how many real cleavage sites fall strictly
  # between them. Realistic protein lengths/enzyme site densities keep this
  # comfortably small (a few thousand pairs at most), so brute-force
  # enumeration is simpler and less bug-prone than trying to bound it cleverly.
  idx_i <- rep(seq_len(k - 1), times = rev(seq_len(k - 1)))
  idx_j <- unlist(lapply(seq_len(k - 1), function(i) seq(i + 1, k)))
  starts <- boundaries[idx_i] + 1L
  ends <- boundaries[idx_j]
  missed <- idx_j - idx_i - 1L

  substrings <- substring(sequence, starts, ends)
  ptm_field <- if (average) "mass_delta_avg" else "mass_delta_mono"
  bare_masses <- if (length(substrings) > 0) sequence_masses_batch(substrings, average = average) else numeric(0)

  ptm_delta <- vapply(seq_along(starts), function(idx) {
    s <- starts[idx]; e <- ends[idx]
    if (length(pf$ptms) == 0) return(0)
    sum(vapply(pf$ptms, function(m) {
      included <- if (identical(m$site, "N-term")) {
        s == 1
      } else if (identical(m$site, "C-term")) {
        e == n
      } else {
        m$site >= s && m$site <= e
      }
      if (included) m[[ptm_field]] else 0
    }, numeric(1)))
  }, numeric(1))

  mass <- bare_masses + ptm_delta
  keep <- mass >= mass_min_da & mass <= mass_max_da
  if (!any(keep)) return(.intact_fallback_candidate(pf, parent_id, mode, r_ref, mz_ref, average))

  starts <- starts[keep]; ends <- ends[keep]; missed <- missed[keep]
  substrings <- substrings[keep]; mass <- mass[keep]

  peptides <- list()
  rows <- vector("list", length(starts))
  for (idx in seq_along(starts)) {
    s <- starts[idx]; e <- ends[idx]
    local_ptms <- Filter(Negate(is.null), lapply(pf$ptms, function(m) {
      if (identical(m$site, "N-term")) {
        if (s == 1) return(ptm(1L, m$mass_delta_mono, m$mass_delta_avg, m$name, m$unimod_id))
        return(NULL)
      }
      if (identical(m$site, "C-term")) {
        if (e == n) return(ptm(e - s + 1L, m$mass_delta_mono, m$mass_delta_avg, m$name, m$unimod_id))
        return(NULL)
      }
      if (m$site >= s && m$site <= e) {
        return(ptm(as.integer(m$site - s + 1L), m$mass_delta_mono, m$mass_delta_avg, m$name, m$unimod_id))
      }
      NULL
    }))

    pep_id <- sprintf("%s_pep%d_%d", parent_id, s, e)
    pep_pf <- proteoform(
      id = pep_id, sequence = substrings[idx], ptms = local_ptms,
      provenance = "module_middledown_digest",
      metadata = list(parent_id = parent_id, start = s, end = e, missed_cleavages = missed[idx], enzyme = enzyme)
    )
    peptides[[pep_id]] <- pep_pf

    fwhm_info <- .best_case_fwhm(pep_pf, mode, average, r_ref, mz_ref)
    ms2_prop <- if (nchar(substrings[idx]) >= 2) {
      mean(fragmentation_propensity(pep_pf, mode = mode)$propensity_score)
    } else {
      NA_real_
    }

    rows[[idx]] <- data.frame(
      id = pep_id, parent_id = parent_id, start = s, end = e, length = e - s + 1L,
      missed_cleavages = missed[idx], mass = mass[idx],
      ms1_fwhm_da = fwhm_info$best_fwhm, ms1_best_z = fwhm_info$best_z,
      ms2_avg_propensity = ms2_prop, ptm_sites_covered = length(local_ptms),
      is_intact = FALSE, stringsAsFactors = FALSE
    )
  }
  candidates <- do.call(rbind, rows)
  rownames(candidates) <- NULL
  # Default ordering follows the two-tier framework: Tier 1 feasibility
  # (fewer missed cleavages, then tighter/better-resolved MS1 peak) first,
  # Tier 2 informativeness (more PTM sites covered) only as a tiebreaker.
  candidates <- candidates[order(candidates$missed_cleavages, candidates$ms1_fwhm_da, -candidates$ptm_sites_covered), ]
  list(candidates = candidates, peptides = peptides)
}

#' Digest every proteoform in a named list, tagging each candidate row with
#' its parent's own label/iso_key so a combined table across several checked
#' proteoforms can still be traced back to where each candidate came from.
#'
#' @param pf_list named list of proteoform objects (e.g. derived()$rows' $pf)
#' @param labels named character vector, same names as pf_list, for display
#' @param iso_keys named character vector, same names as pf_list, mapping
#'   each parent to the transcript id its exon axis should be built from
#' @param enzyme,mass_min_da,mass_max_da,mode,r_ref,mz_ref,average passed to digest_proteoform()
#' @return list(candidates = combined data.frame (NULL if nothing survived),
#'   peptides = combined named list of proteoform objects,
#'   parent_label = named vector id -> parent's display label,
#'   parent_iso_key = named vector id -> parent's iso_key,
#'   parent_counts = data.frame(parent_id, label, n_candidates, is_intact_only)
#'     -- ALWAYS one row per proteoform in pf_list. digest_proteoform() now
#'     always returns at least one row (falling back to the intact protein
#'     when no real digest fragment lands in the mass window -- see its doc
#'     comment), so n_candidates is never really 0 here; is_intact_only flags
#'     the fallback case so callers can still surface "no digest fragment in
#'     window for this protease" instead of it silently looking like a
#'     normal single-candidate result)
digest_proteoform_set <- function(pf_list, labels, iso_keys, enzyme, mass_min_da, mass_max_da,
                                   mode = c("denatured", "native"), r_ref = 120000, mz_ref = 200, average = FALSE) {
  mode <- match.arg(mode)
  all_candidates <- list()
  all_peptides <- list()
  parent_label <- character(0)
  parent_iso_key <- character(0)
  parent_counts <- vector("list", length(pf_list))
  names(parent_counts) <- names(pf_list)

  for (id in names(pf_list)) {
    # parent_id = id (the CALLER's own list key, e.g. derived()$rows'
    # "<transcript>#bare"), not the proteoform object's own $id field --
    # those two differ (the bare proteoform's $id is just "<transcript>",
    # no "#bare" suffix), and candidate ids/parent_id need to match `id`
    # here since that's what labels/iso_keys (and downstream lookups by
    # derived()$rows name) are keyed by.
    d <- digest_proteoform(pf_list[[id]], enzyme, mass_min_da, mass_max_da,
                            mode = mode, r_ref = r_ref, mz_ref = mz_ref, average = average,
                            parent_id = id)
    n_here <- if (is.null(d$candidates)) 0L else nrow(d$candidates)
    is_intact_only <- !is.null(d$candidates) && all(d$candidates$is_intact)
    parent_counts[[id]] <- data.frame(parent_id = id, label = labels[[id]] %||% id, n_candidates = n_here,
                                       is_intact_only = is_intact_only, stringsAsFactors = FALSE)
    if (is.null(d$candidates)) next
    all_candidates[[id]] <- d$candidates
    all_peptides <- c(all_peptides, d$peptides)
    new_ids <- d$candidates$id
    parent_label[new_ids] <- labels[[id]] %||% id
    parent_iso_key[new_ids] <- iso_keys[[id]] %||% id
  }

  candidates <- if (length(all_candidates) > 0) do.call(rbind, all_candidates) else NULL
  if (!is.null(candidates)) rownames(candidates) <- NULL
  list(candidates = candidates, peptides = all_peptides, parent_label = parent_label,
       parent_iso_key = parent_iso_key, parent_counts = do.call(rbind, parent_counts))
}

#' Re-base a parent's exon/residue table onto a peptide's own local numbering
#' (1..peptide_length), clipping to the peptide's covered range. Used for
#' Section 2 (confounder search), whose exon_blocks are plain local residue
#' ranges with no separate shared-axis step (unlike Section 1, which instead
#' offsets the LOCAL position before looking it up on the parent's shared
#' axis -- see build_section1_payload()'s residue_offset argument).
#'
#' @param exon_table data.frame with residue_start/residue_end columns, in
#'   the PARENT's own full-length numbering
#' @param start,end the peptide's own start/end, in that same parent numbering
#' @return exon_table subset, clipped and re-based so residue_start/residue_end
#'   are 1-based positions in the peptide's own sequence
shift_exon_table_for_peptide <- function(exon_table, start, end) {
  if (is.null(exon_table) || nrow(exon_table) == 0) return(exon_table)
  te <- exon_table[exon_table$residue_end >= start & exon_table$residue_start <= end, ]
  if (nrow(te) == 0) return(te)
  te$residue_start <- pmax(te$residue_start, start) - start + 1L
  te$residue_end <- pmin(te$residue_end, end) - start + 1L
  te
}

#' How much of each parent's own sequence is covered by the currently
#' selected (checked) middle-down peptide candidates -- the "non-redundant
#' coverage" Tier 2 point, which only makes sense relative to what's already
#' selected, so it's a live summary rather than a per-row column.
#'
#' @param candidates combined candidate data.frame (from digest_proteoform_set())
#' @param selected_ids character vector of candidate ids currently checked
#' @param parent_lengths named integer vector, parent_id -> full sequence length
#' @return data.frame(parent_id, covered_residues, total_residues, pct)
digestion_coverage_summary <- function(candidates, selected_ids, parent_lengths) {
  if (is.null(candidates) || length(selected_ids) == 0) {
    return(data.frame(parent_id = character(0), covered_residues = integer(0),
                       total_residues = integer(0), pct = numeric(0)))
  }
  sel <- candidates[candidates$id %in% selected_ids, ]
  parents <- unique(sel$parent_id)
  rows <- lapply(parents, function(pid) {
    rows_p <- sel[sel$parent_id == pid, ]
    covered <- rep(FALSE, parent_lengths[[pid]] %||% 0L)
    for (r in seq_len(nrow(rows_p))) {
      covered[rows_p$start[r]:rows_p$end[r]] <- TRUE
    }
    total <- length(covered)
    data.frame(parent_id = pid, covered_residues = sum(covered), total_residues = total,
               pct = if (total > 0) round(100 * sum(covered) / total, 1) else NA_real_)
  })
  do.call(rbind, rows)
}
