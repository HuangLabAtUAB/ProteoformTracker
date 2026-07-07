library(testthat)

source("../../R/proteoform_schema.R")
source("../../R/mass_calculation.R")
source("../../R/reference_proteome_index.R")

mass_engine_available <- tryCatch({
  init_mass_calculation_engine()
  TRUE
}, error = function(e) FALSE)
