# Cross-proteoform MS2 fragment collision check + structural divergence-
# point finder (design spec section 6).
#
# Two distinct reasons a confounding protein matters at MS2, flagged
# separately per spec: database/scoring ambiguity (a fragment mass shared
# with another proteoform is weaker confirmatory evidence) vs
# co-isolation/chimeric spectra risk (an MS1-level precursor-window
# question, already handled by search_confounding_proteins()/
# search_mz_collisions() in R/ms1_scoring.R and R/mz_collision_index.R).
# This file only addresses the first: does a given fragment's mass collide
# with another proteoform's fragment ladder at all.
#
# Because cross-gene proteins essentially never share extended sequence
# identity, most fragment collisions should clear as resolvable -- this
# module's value is catching exceptions, and (for relevant isoforms)
# locating exactly where two ladders structurally diverge.

#' Check which of a target's b/y fragments have a mass-matching counterpart
#' in each candidate's own fragment ladder. Works uniformly whether
#' `candidates` are relevant isoforms or confounding proteins -- per the
#' design spec, that's exactly what should differ between the two:
#' cross-gene proteins rarely share fragments even at similar intact mass,
#' while isoforms typically do, below their sequence-divergence point.
#'
#' Match tolerance is derived from the resolving-power model (same
#' architecture as the MS1 searches), not a fixed Da value: two fragment
#' masses are called "matched" if they fall within safety_margin x FWHM at
#' the given fragment charge state.
#'
#' @param target proteoform object
#' @param candidates named or unnamed list of proteoform objects
#' @param average use average mass instead of monoisotopic
#' @param charge fragment charge state used to size the mass tolerance
#' @param r_ref,mz_ref Orbitrap resolving-power settings
#' @param safety_margin multiplier applied to FWHM to size the match tolerance
#' @param script_path path to python/ptracker_mass.py
#' @return list(target_ladder, per_candidate = list of data.frames --
#'   ion_type, cleavage_position, mass, matched -- one per candidate, in the
#'   same order as `candidates`)
fragment_mass_collision_check <- function(target, candidates, average = FALSE, charge = 1,
                                           r_ref = 120000, mz_ref = 200,
                                           safety_margin = DEFAULT_SAFETY_MARGIN,
                                           script_path = "python/ptracker_mass.py") {
  if (!inherits(target, "proteoform")) {
    stop("fragment_mass_collision_check() requires a proteoform target")
  }
  if (!is.list(candidates) || !all(vapply(candidates, inherits, logical(1), "proteoform"))) {
    stop("fragment_mass_collision_check() requires a list of proteoform candidates")
  }

  target_ladder <- generate_fragment_ladder(target, average = average, script_path = script_path)
  target_long <- rbind(
    data.frame(
      ion_type = "b", cleavage_position = target_ladder$cleavage_position,
      mass = target_ladder$b_mass, stringsAsFactors = FALSE
    ),
    data.frame(
      ion_type = "y", cleavage_position = target_ladder$cleavage_position,
      mass = target_ladder$y_mass, stringsAsFactors = FALSE
    )
  )
  target_long$tolerance <- safety_margin * fwhm_mass(
    target_long$mass, mz_for_charge(target_long$mass, charge), r_ref, mz_ref
  )

  per_candidate <- lapply(candidates, function(cand) {
    cand_ladder <- generate_fragment_ladder(cand, average = average, script_path = script_path)
    cand_masses <- c(cand_ladder$b_mass, cand_ladder$y_mass)

    result <- target_long
    result$matched <- vapply(seq_len(nrow(result)), function(i) {
      any(abs(cand_masses - result$mass[i]) <= result$tolerance[i])
    }, logical(1))
    result
  })
  if (!is.null(names(candidates))) {
    names(per_candidate) <- names(candidates)
  }

  list(target_ladder = target_ladder, per_candidate = per_candidate)
}

#' Where a target's and a candidate's fragment ladders structurally
#' diverge. The b-ion ladder is identical up to the N-terminal colinear
#' prefix length; the y-ion ladder up to the C-terminal colinear suffix
#' length. Per the design spec's key structural property: two isoforms'
#' ladders are identical up to the sequence divergence point (typically a
#' differential exon) and uniformly offset beyond it.
#'
#' "Last identical residue from a terminus" is exactly the longest common
#' prefix/suffix -- a direct O(n) character comparison from each end, not
#' a general alignment problem. An earlier version routed this through
#' align_sequences() (Needleman-Wunsch, R/ptm_site_mapping.R), which is
#' the right tool for mapping an *internal* PTM site across an indel, but
#' is the wrong tool here: with a plain linear gap penalty, a long indel
#' can have more than one equal-scoring placement (no cost difference
#' between one contiguous gap and several fragmented ones), so the
#' aligner can return *a* valid alignment whose gap start doesn't match
#' the true longest common prefix. Caught validating against real CD44
#' isoform data -- the true divergence (a 313-residue exon-skip) is at
#' residue 223, but the aligner placed its gap at residue 203, undercounting
#' by 20 residues. Direct prefix/suffix comparison has no such ambiguity.
#'
#' @param target,candidate proteoform objects
#' @return list(b_shared_length, y_shared_length) -- number of residues
#'   from the N-/C-terminus, respectively, where the two sequences remain
#'   identical (same residue at the same relative position)
find_fragment_divergence_point <- function(target, candidate) {
  if (!inherits(target, "proteoform") || !inherits(candidate, "proteoform")) {
    stop("find_fragment_divergence_point() requires proteoform objects")
  }
  target_chars <- strsplit(target$sequence, "")[[1]]
  candidate_chars <- strsplit(candidate$sequence, "")[[1]]
  n_target <- length(target_chars)
  n_candidate <- length(candidate_chars)

  max_prefix <- min(n_target, n_candidate)
  b_shared_length <- 0
  for (i in seq_len(max_prefix)) {
    if (target_chars[i] != candidate_chars[i]) break
    b_shared_length <- i
  }

  y_shared_length <- 0
  for (k in seq_len(max_prefix)) {
    if (target_chars[n_target - k + 1] != candidate_chars[n_candidate - k + 1]) break
    y_shared_length <- k
  }

  list(b_shared_length = b_shared_length, y_shared_length = y_shared_length)
}
