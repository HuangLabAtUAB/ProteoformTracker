library(testthat)

source("../../R/proteoform_schema.R")
source("../../R/unimod_table.R")
source("../../R/ptm_site_mapping.R")
source("../../R/mass_calculation.R")
source("../../R/reference_proteome_index.R")
source("../../R/resolving_power.R")
source("../../R/charge_envelope.R")
source("../../R/isotope_envelope.R")
source("../../R/ms1_scoring.R")
source("../../R/mz_collision_index.R")
source("../../R/fragment_ladder.R")

TEST_PY_SCRIPT <- "../../python/ptracker_mass.py"

mass_engine_available <- tryCatch({
  init_mass_calculation_engine()
  TRUE
}, error = function(e) FALSE)
