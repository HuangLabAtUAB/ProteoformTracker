# Builds the m/z-domain collision index (R/mz_collision_index.R) from the
# already-built reference mass index (scripts/build_reference_proteome.R
# must have been run first). Built once per acquisition mode, since the
# charge-state envelope depends on it.

source("R/proteoform_schema.R")
source("R/resolving_power.R")
source("R/charge_envelope.R")
source("R/mz_collision_index.R")

mass_index_path <- "data/reference_mass_index.rds"
if (!file.exists(mass_index_path)) {
  stop("No reference mass index at ", mass_index_path, " -- run scripts/build_reference_proteome.R first")
}
mass_index <- readRDS(mass_index_path)

for (mode in c("denatured", "native")) {
  out_path <- paste0("data/reference_mz_index_", mode, ".rds")
  message("Building m/z collision index (", mode, ")...")
  mz_index <- build_reference_mz_index(mass_index, mode = mode)
  saveRDS(mz_index, out_path)
  message("  ", nrow(mz_index), " peaks across ", length(unique(mz_index$id)), " proteins -> ", out_path)
}
