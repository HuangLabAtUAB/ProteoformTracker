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

#' Monoisotopic or average mass of a bare amino acid sequence (no PTMs),
#' via pyteomics.mass.calculate_mass.
#'
#' @param sequence amino acid sequence (standard 20 AA letters)
#' @param average if TRUE, return average mass; otherwise monoisotopic
sequence_mass <- function(sequence, average = FALSE) {
  pm <- .pyteomics_mass()
  pm$calculate_mass(sequence = sequence, average = average)
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
