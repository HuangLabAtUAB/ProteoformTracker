# SETUP AND DATA LOADING
#===============================================================================

library(shiny)
library(reticulate)

# SOURCE CORE MODULES (Phase 1: data model + mass engine)
source("R/proteoform_schema.R")
source("R/mass_calculation.R")
source("R/reference_proteome_index.R")

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
