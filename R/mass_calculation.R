# Proteoform mass calculation.
# Wraps pyteomics (via reticulate) for residue-mass tables rather than
# reimplementing them. PTM contributions are applied as flat mass deltas
# from the Proteoform object's $ptms list (Unimod-style).

PTRACKER_PY_ENV <- "proteoformtracker-py"

#' Point reticulate at the project's dedicated virtualenv and confirm
#' pyteomics is importable. Call once at app/session startup.
init_mass_calculation_engine <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("package 'reticulate' is required for mass calculation")
  }
  if (reticulate::virtualenv_exists(PTRACKER_PY_ENV)) {
    reticulate::use_virtualenv(PTRACKER_PY_ENV, required = TRUE)
  }
  if (!reticulate::py_module_available("pyteomics")) {
    stop(
      "pyteomics is not available in the active Python environment. ",
      "Run: Rscript scripts/setup_python_env.R"
    )
  }
  invisible(TRUE)
}

.pyteomics_mass <- function() {
  reticulate::import("pyteomics.mass", convert = TRUE)
}

#' Lazily source python/ptracker_mass.py once and cache the resulting
#' environment (batch_calculate_mass, batch_ion_mass). `script_path` is
#' relative to the caller's working directory (project root for the
#' app/scripts; tests pass their own relative path).
.ptracker_mass_env <- local({
  env <- NULL
  function(script_path = "python/ptracker_mass.py") {
    if (is.null(env)) {
      # Only cache on success: assigning `env` before source_python() could
      # throw would permanently poison the cache with an empty environment
      # after any transient failure (wrong cwd, missing file), breaking
      # every subsequent call for the rest of the R session.
      new_env <- new.env()
      reticulate::source_python(script_path, envir = new_env)
      env <<- new_env
    }
    env
  }
})

.batch_mass_fn <- function(script_path = "python/ptracker_mass.py") {
  .ptracker_mass_env(script_path)$batch_calculate_mass
}

.batch_ion_mass_fn <- function(script_path = "python/ptracker_mass.py") {
  .ptracker_mass_env(script_path)$batch_ion_mass
}

.isotope_pattern_fn <- function(script_path = "python/ptracker_mass.py") {
  .ptracker_mass_env(script_path)$isotope_pattern
}

#' Exact-composition isotope pattern for a bare amino acid sequence, binned
#' to nominal-mass peaks (one peak per integer Da offset from the
#' monoisotopic mass -- see python/ptracker_mass.py's isotope_pattern() doc
#' comment for why binning, not the raw fine-structure isotopologue list, is
#' the right thing to compare against a real spectrum).
#'
#' Deliberately NOT the Averagine model: Averagine approximates elemental
#' composition from mass alone for when the sequence is unknown. Every
#' caller here already has the real sequence, so the real composition
#' (summed real residue formulas) is used -- strictly more accurate at the
#' same computational cost (milliseconds, even for an intact protein).
#'
#' @param sequence bare amino acid sequence (standard 20 AA letters). PTM
#'   mass deltas are NOT applied here -- add them as a flat shift to the
#'   returned masses afterwards (see proteoform_isotope_pattern()), the same
#'   flat-delta treatment PTMs get everywhere else in this codebase.
#' @param prob_to_cover total probability mass to enumerate (default 99.99%
#'   -- the excluded tail is genuinely negligible)
#' @param script_path path to python/ptracker_mass.py
#' @return data.frame(mass, prob), ascending by mass, prob summing to 1
sequence_isotope_pattern <- function(sequence, prob_to_cover = 0.9999,
                                      script_path = "python/ptracker_mass.py") {
  fn <- .isotope_pattern_fn(script_path)
  result <- fn(sequence, prob_to_cover)
  data.frame(mass = unlist(result$masses), prob = unlist(result$probs))
}

#' proteoform_mass()'s isotope-pattern equivalent: the bare sequence's exact
#' isotope pattern (sequence_isotope_pattern()) with every peak's mass
#' shifted by the proteoform's total PTM mass delta -- a PTM changes total
#' mass but this codebase only ever tracks that as a flat delta (no
#' elemental-formula representation for modifications), so the pattern's
#' SHAPE (relative peak spacing/intensities) is assumed unchanged by the
#' PTM and only its position shifts. Reasonable for a single small PTM;
#' would understate the shape change for something isotopically unusual
#' (e.g. a heavy-isotope metabolic label) -- not a case this app handles.
#'
#' @param p a proteoform object (see proteoform_schema.R)
#' @param average if TRUE, shift by the average-mass PTM delta; otherwise
#'   monoisotopic (must match how the rest of the caller's pipeline computes
#'   mass, same convention as proteoform_mass())
#' @param prob_to_cover passed to sequence_isotope_pattern()
#' @return data.frame(mass, prob), ascending by mass, prob summing to 1
proteoform_isotope_pattern <- function(p, average = FALSE, prob_to_cover = 0.9999) {
  if (!inherits(p, "proteoform")) stop("proteoform_isotope_pattern() requires a proteoform object")
  pattern <- sequence_isotope_pattern(p$sequence, prob_to_cover = prob_to_cover)
  ptm_delta <- sum(vapply(
    p$ptms,
    function(m) if (average) m$mass_delta_avg else m$mass_delta_mono,
    numeric(1)
  ))
  pattern$mass <- pattern$mass + ptm_delta
  pattern
}

#' Monoisotopic or average mass of a bare amino acid sequence (no PTMs),
#' via pyteomics.mass.calculate_mass.
#'
#' @param sequence amino acid sequence (standard 20 AA letters)
#' @param average if TRUE, return average mass; otherwise monoisotopic
sequence_mass <- function(sequence, average = FALSE) {
  pm <- .pyteomics_mass()
  pm$calculate_mass(sequence = sequence, average = average)
}

#' Monoisotopic or average masses for many sequences in a single Python
#' round-trip (a list comprehension over pyteomics.mass.calculate_mass),
#' rather than one reticulate call per sequence. Used when building the
#' reference-proteome mass index (tens of thousands of sequences).
#'
#' @param sequences character vector of amino acid sequences
#' @param average if TRUE, return average mass; otherwise monoisotopic
#' @param script_path path to python/ptracker_mass.py, relative to the
#'   current working directory
#' @return numeric vector of masses, same length/order as sequences
sequence_masses_batch <- function(sequences, average = FALSE,
                                   script_path = "python/ptracker_mass.py") {
  batch_fn <- .batch_mass_fn(script_path)
  # unname(): a named list would convert to a Python dict, and iterating a
  # dict yields keys (the names) instead of the sequences themselves.
  unlist(batch_fn(as.list(unname(sequences)), average))
}

#' Neutral fragment-ion masses (b/y/c/z, etc. -- see pyteomics.mass.std_ion_comp)
#' for many subsequences in a single Python round-trip. Neutral, not
#' charged: add PROTON_MASS * z and divide by z (mz_for_charge(), see
#' R/resolving_power.R) to get an m/z at a given fragment charge state,
#' the same convention used for intact-proteoform mass throughout this
#' project.
#'
#' @param sequences character vector of amino acid subsequences (e.g.
#'   N-terminal prefixes for b-ions, C-terminal suffixes for y-ions)
#' @param ion_type ion series, e.g. "b" or "y"
#' @param average if TRUE, return average mass; otherwise monoisotopic
#' @param script_path path to python/ptracker_mass.py, relative to the
#'   current working directory
#' @return numeric vector of neutral masses, same length/order as sequences
sequence_ion_masses_batch <- function(sequences, ion_type, average = FALSE,
                                       script_path = "python/ptracker_mass.py") {
  batch_fn <- .batch_ion_mass_fn(script_path)
  unlist(batch_fn(as.list(unname(sequences)), ion_type, average))
}

#' Total mass of a Proteoform: bare-sequence mass + sum of PTM mass deltas.
#'
#' @param p a proteoform object (see proteoform_schema.R)
#' @param average if TRUE, use average mass; otherwise monoisotopic
#' @return list(mass, sequence_mass, ptm_delta, average)
proteoform_mass <- function(p, average = FALSE) {
  if (!inherits(p, "proteoform")) stop("proteoform_mass() requires a proteoform object")
  base_mass <- sequence_mass(p$sequence, average = average)
  ptm_delta <- sum(vapply(
    p$ptms,
    function(m) if (average) m$mass_delta_avg else m$mass_delta_mono,
    numeric(1)
  ))
  list(
    mass = base_mass + ptm_delta,
    sequence_mass = base_mass,
    ptm_delta = ptm_delta,
    average = average
  )
}
