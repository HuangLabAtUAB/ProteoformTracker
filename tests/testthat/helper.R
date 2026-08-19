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
source("../../R/gene_isoform_data.R")
source("../../R/svg_render.R")
source("../../R/viz_json.R")
source("../../R/fragment_ladder.R")
source("../../R/fragmentation_propensity.R")
source("../../R/fragmentation_propensity_rf.R")
source("../../R/fragment_collision.R")
source("../../R/digestion.R")
source("../../R/reference_exon_index.R")
source("../../R/ptm_spec_parser.R")
source("../../R/ensembl_protein_fetch.R")

TEST_PY_SCRIPT <- "../../python/ptracker_mass.py"

mass_engine_available <- tryCatch({
  init_mass_calculation_engine()
  TRUE
}, error = function(e) FALSE)

TEST_RF_MODEL_PATH <- "../../data/fragmentation_propensity_rf.rds"
propensity_rf_model <- tryCatch({
  if (!requireNamespace("ranger", quietly = TRUE) || !file.exists(TEST_RF_MODEL_PATH)) stop("unavailable")
  readRDS(TEST_RF_MODEL_PATH)
}, error = function(e) NULL)
rf_model_available <- !is.null(propensity_rf_model)
