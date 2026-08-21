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
# assign() into .GlobalEnv explicitly: testthat::test_dir() sources this file
# into a private helper environment, not .GlobalEnv, but fragmentation_propensity_rf()
# is defined (via the R/*.R source() calls above, which default to local = FALSE)
# with .GlobalEnv as its closure -- so its `rf_model = propensity_rf_model` default
# argument can only resolve if the object actually lives in .GlobalEnv. This mirrors
# how the deployed app works too: Shiny sources global.R into .GlobalEnv.
propensity_rf_model <- tryCatch({
  if (!requireNamespace("ranger", quietly = TRUE) || !file.exists(TEST_RF_MODEL_PATH)) stop("unavailable")
  readRDS(TEST_RF_MODEL_PATH)
}, error = function(e) NULL)
assign("propensity_rf_model", propensity_rf_model, envir = .GlobalEnv)
rf_model_available <- !is.null(propensity_rf_model)

TEST_GLM_MODEL_PATH <- "../../data/fragmentation_propensity_glm.rds"
propensity_glm_model <- tryCatch({
  if (!file.exists(TEST_GLM_MODEL_PATH)) stop("unavailable")
  readRDS(TEST_GLM_MODEL_PATH)
}, error = function(e) NULL)
assign("propensity_glm_model", propensity_glm_model, envir = .GlobalEnv)
glm_model_available <- !is.null(propensity_glm_model)
