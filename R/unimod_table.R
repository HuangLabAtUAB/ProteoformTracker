# Curated Unimod lookup table: common PTMs relevant to top-down proteomics,
# so a caller can attach a modification by name (e.g. "Phospho") instead of
# having to already know its mass delta. Not the full Unimod database --
# unimod_ptm() errors with the available names if a modification isn't in
# this table; use ptm() directly with a manual mass delta for anything else.
# Mass deltas are standard Unimod values (monoisotopic and average, Da).

UNIMOD_PTM_TABLE <- data.frame(
  name = c(
    "Phospho", "Acetyl", "Methyl", "Dimethyl", "Trimethyl", "Oxidation",
    "GG", "Deamidated", "Gln->pyro-Glu", "Glu->pyro-Glu", "Palmitoyl",
    "Sulfo", "Carbamidomethyl", "Formyl", "Nitro", "Amidated"
  ),
  unimod_id = c(
    "UNIMOD:21", "UNIMOD:1", "UNIMOD:34", "UNIMOD:36", "UNIMOD:37", "UNIMOD:35",
    "UNIMOD:121", "UNIMOD:7", "UNIMOD:28", "UNIMOD:27", "UNIMOD:47",
    "UNIMOD:40", "UNIMOD:4", "UNIMOD:122", "UNIMOD:354", "UNIMOD:2"
  ),
  mass_delta_mono = c(
    79.966331, 42.010565, 14.015650, 28.031300, 42.046950, 15.994915,
    114.042927, 0.984016, -17.026549, -18.010565, 238.229666,
    79.956815, 57.021464, 27.994915, 44.985078, -0.984016
  ),
  mass_delta_avg = c(
    79.9799, 42.0367, 14.0266, 28.0532, 42.0797, 15.9994,
    114.1026, 0.9848, -17.0305, -18.0153, 238.4136,
    80.0642, 57.0513, 28.0101, 44.9976, -0.9848
  ),
  description = c(
    "Phosphorylation (Ser/Thr/Tyr)",
    "Acetylation (protein N-term or Lys)",
    "Methylation (Lys/Arg)",
    "Dimethylation (Lys/Arg)",
    "Trimethylation (Lys)",
    "Oxidation (Met/Trp/Pro/Cys)",
    "Ubiquitin/SUMO remnant (di-glycine, Lys)",
    "Deamidation (Asn/Gln)",
    "Pyroglutamate formation from N-terminal Gln",
    "Pyroglutamate formation from N-terminal Glu",
    "Palmitoylation (Cys)",
    "Sulfation (Tyr)",
    "Carbamidomethylation (Cys, common alkylation artifact)",
    "Formylation (protein N-term or Lys)",
    "Nitration (Tyr)",
    "C-terminal amidation"
  ),
  stringsAsFactors = FALSE
)

#' Build a ptm() entry by looking up its mass deltas in UNIMOD_PTM_TABLE by
#' name (case-insensitive), instead of supplying mass_delta_mono directly.
#'
#' @param name modification name, e.g. "Phospho" (see list_unimod_mods())
#' @param site 1-based position in the mature sequence, or "N-term"/"C-term"
unimod_ptm <- function(name, site) {
  idx <- match(tolower(name), tolower(UNIMOD_PTM_TABLE$name))
  if (is.na(idx)) {
    stop(
      "Unknown modification '", name, "'. Available: ",
      paste(UNIMOD_PTM_TABLE$name, collapse = ", "),
      ". For anything else, use ptm() directly with a manual mass delta."
    )
  }
  row <- UNIMOD_PTM_TABLE[idx, ]
  ptm(
    site = site,
    mass_delta_mono = row$mass_delta_mono,
    mass_delta_avg = row$mass_delta_avg,
    name = row$name,
    unimod_id = row$unimod_id
  )
}

#' List the modifications available to unimod_ptm().
list_unimod_mods <- function() {
  UNIMOD_PTM_TABLE[, c("name", "unimod_id", "mass_delta_mono", "mass_delta_avg", "description")]
}
