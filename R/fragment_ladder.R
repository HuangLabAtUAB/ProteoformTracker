# Fragment ladder generation (deterministic arithmetic, per design spec).
# Generates cumulative N-terminal (b-ion) and C-terminal (y-ion) fragment
# masses for every backbone cleavage position -- HCD/CID ion series for v1
# (ETD/ECD's c/z ions and UVPD are deferred).
#
# Because terminal fragment ladders are cumulative sums along the sequence,
# two proteoforms' ladders are identical up to their sequence divergence
# point and uniformly offset by Delta-mass beyond it (see
# R/fragment_comparison.R for the comparison that operationalizes this).

#' Sum the PTM mass deltas that fall within an N-terminal (b-ion) fragment
#' covering positions [1, cleavage_position]. A b-ion fragment always
#' retains the free N-terminus and never reaches the true C-terminus (the
#' cleavage position is always interior, 1..n-1).
.ptm_delta_for_b_ion <- function(ptms, cleavage_position, field) {
  if (length(ptms) == 0) {
    return(0)
  }
  sum(vapply(ptms, function(mod) {
    included <- if (identical(mod$site, "N-term")) {
      TRUE
    } else if (identical(mod$site, "C-term")) {
      FALSE
    } else {
      mod$site <= cleavage_position
    }
    if (included) mod[[field]] else 0
  }, numeric(1)))
}

#' Sum the PTM mass deltas that fall within a C-terminal (y-ion) fragment
#' covering positions [cleavage_position + 1, n]. A y-ion fragment always
#' retains the true C-terminus and never reaches the free N-terminus.
.ptm_delta_for_y_ion <- function(ptms, cleavage_position, field) {
  if (length(ptms) == 0) {
    return(0)
  }
  sum(vapply(ptms, function(mod) {
    included <- if (identical(mod$site, "N-term")) {
      FALSE
    } else if (identical(mod$site, "C-term")) {
      TRUE
    } else {
      mod$site > cleavage_position
    }
    if (included) mod[[field]] else 0
  }, numeric(1)))
}

#' Generate the full b/y fragment ladder for a proteoform: neutral fragment
#' mass at every backbone cleavage position, PTM deltas included wherever
#' the modified site falls inside that fragment.
#'
#' @param proteoform a proteoform object (see proteoform_schema.R)
#' @param average use average mass instead of monoisotopic
#' @param script_path path to python/ptracker_mass.py, relative to the
#'   current working directory
#' @return data.frame(cleavage_position, n_term_length, c_term_length,
#'   b_mass, y_mass) -- one row per interior backbone position (1..n-1)
generate_fragment_ladder <- function(proteoform, average = FALSE,
                                      script_path = "python/ptracker_mass.py") {
  if (!inherits(proteoform, "proteoform")) {
    stop("generate_fragment_ladder() requires a proteoform object")
  }
  sequence <- proteoform$sequence
  n <- nchar(sequence)
  if (n < 2) {
    stop("sequence must have at least 2 residues to generate fragment ions")
  }

  cleavage_positions <- seq_len(n - 1)
  n_term_prefixes <- substring(sequence, 1, cleavage_positions)
  c_term_suffixes <- substring(sequence, cleavage_positions + 1, n)

  b_mass_bare <- sequence_ion_masses_batch(n_term_prefixes, ion_type = "b", average = average, script_path = script_path)
  y_mass_bare <- sequence_ion_masses_batch(c_term_suffixes, ion_type = "y", average = average, script_path = script_path)

  ptm_field <- if (average) "mass_delta_avg" else "mass_delta_mono"
  b_ptm_delta <- vapply(cleavage_positions, function(i) {
    .ptm_delta_for_b_ion(proteoform$ptms, i, ptm_field)
  }, numeric(1))
  y_ptm_delta <- vapply(cleavage_positions, function(i) {
    .ptm_delta_for_y_ion(proteoform$ptms, i, ptm_field)
  }, numeric(1))

  data.frame(
    cleavage_position = cleavage_positions,
    n_term_length = cleavage_positions,
    c_term_length = n - cleavage_positions,
    b_mass = b_mass_bare + b_ptm_delta,
    y_mass = y_mass_bare + y_ptm_delta
  )
}
