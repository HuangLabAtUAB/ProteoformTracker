# SETUP AND DATA LOADING
#===============================================================================

library(shiny)
library(reticulate)

# SOURCE CORE MODULES (Phase 1: data model + mass engine)
source("R/proteoform_schema.R")
source("R/unimod_table.R")
source("R/ptm_site_mapping.R")
source("R/mass_calculation.R")
source("R/reference_proteome_index.R")
source("R/resolving_power.R")
source("R/charge_envelope.R")
source("R/isotope_envelope.R")
source("R/ms1_scoring.R")
source("R/mz_collision_index.R")
source("R/fragment_ladder.R")
source("R/fragmentation_propensity.R")
source("R/fragment_collision.R")
source("R/reference_exon_index.R")

# Point reticulate at the project's Python venv and confirm pyteomics loads.
tryCatch(
  init_mass_calculation_engine(),
  error = function(e) {
    message("WARNING: mass calculation engine not ready: ", conditionMessage(e))
  }
)

reference_mass_index_path <- "data/reference_mass_index.rds"
reference_mass_index <- if (file.exists(reference_mass_index_path)) {
  load_reference_mass_index(reference_mass_index_path)
} else {
  message("No reference-proteome mass index found at ", reference_mass_index_path,
          " -- confounding-protein search unavailable until build_reference_mass_index() is run.")
  NULL
}

reference_mz_index_paths <- list(
  denatured = "data/reference_mz_index_denatured.rds",
  native = "data/reference_mz_index_native.rds"
)
reference_mz_index <- lapply(reference_mz_index_paths, function(p) {
  if (file.exists(p)) load_reference_mz_index(p) else NULL
})
if (all(vapply(reference_mz_index, is.null, logical(1)))) {
  message("No m/z collision index found -- run scripts/build_mz_collision_index.R to enable m/z-domain confounding-protein search.")
}

reference_exon_index_path <- "data/reference_exon_index.rds"
reference_exon_index <- if (file.exists(reference_exon_index_path)) {
  load_reference_exon_index(reference_exon_index_path)
} else {
  message("No exon index found -- run scripts/build_exon_index.R to enable exon-track lookups by gene/transcript.")
  NULL
}
