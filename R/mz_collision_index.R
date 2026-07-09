# Offline precomputed m/z-domain collision index for the reference
# proteome: every reference protein's predicted charge-state peaks, sorted
# by m/z, for fast range queries.
#
# This is a different axis of confounding-protein selection than
# search_confounding_proteins() (R/ms1_scoring.R), which only finds
# candidates whose *intact mass* is near the target's. A protein whose mass
# is far from the target can still produce a charge-state peak that lands
# on the target's own m/z (e.g. a 12 kDa and a 16 kDa protein whose masses
# happen to sit in a ~4:3 ratio will collide at several charge-state
# pairs) -- a mass-window search would never surface it. Built per
# acquisition mode, since the charge-state envelope depends on it.

#' Build the m/z collision index: every entry in a reference mass index
#' (from build_/load_reference_mass_index()) gets its predicted charge-state
#' envelope computed and flattened into one big (id, z, mz, fwhm_mz) table,
#' sorted by m/z. Pure R (no reticulate calls) -- mass is read from the
#' already-computed mass_index, not recomputed.
#'
#' @param mass_index reference-proteome mass index (id, sequence, mass, ...)
#' @param mode "denatured" or "native"
#' @param r_ref,mz_ref Orbitrap resolving-power settings
#' @return data.frame(id, z, mz, relative_intensity, fwhm_mz), sorted by mz
build_reference_mz_index <- function(mass_index, mode = c("denatured", "native"),
                                      r_ref = 120000, mz_ref = 200) {
  mode <- match.arg(mode)

  rows <- vector("list", nrow(mass_index))
  for (i in seq_len(nrow(mass_index))) {
    envelope <- predict_charge_envelope(mass_index$sequence[i], mass_index$mass[i], mode = mode)
    envelope$id <- mass_index$id[i]
    envelope$fwhm_mz <- fwhm_mz(envelope$mz, r_ref = r_ref, mz_ref = mz_ref)
    rows[[i]] <- envelope
  }

  index <- do.call(rbind, rows)
  index <- index[order(index$mz), c("id", "z", "mz", "relative_intensity", "fwhm_mz")]
  rownames(index) <- NULL
  index
}

#' Load a previously built m/z collision index.
load_reference_mz_index <- function(path) {
  readRDS(path)
}

#' Query the m/z collision index for entries colliding with one target peak.
#' Two-stage: a generous range lookup (findInterval binary search, since
#' candidate FWHM at nearby m/z is close in magnitude to the target peak's
#' own FWHM but not known in advance), then an exact distance/window check
#' using each candidate's actual FWHM.
#'
#' @param mz_index m/z collision index (from build_/load_reference_mz_index())
#' @param peak_mz,peak_fwhm_mz the target's own peak position and FWHM
#' @param safety_margin multiplier applied to max(target FWHM, candidate FWHM)
#' @param exclude_id id to exclude from results
#' @return data.frame subset of mz_index (plus distance column) that collide
#'   with this one peak
query_mz_collisions_for_peak <- function(mz_index, peak_mz, peak_fwhm_mz,
                                          safety_margin = 1.0, exclude_id = NULL) {
  search_window <- 2 * safety_margin * peak_fwhm_mz # generous net; refined exactly below
  lo_idx <- findInterval(peak_mz - search_window, mz_index$mz) + 1L
  hi_idx <- findInterval(peak_mz + search_window, mz_index$mz)
  if (lo_idx > hi_idx) {
    return(mz_index[0, ])
  }

  candidates <- mz_index[lo_idx:hi_idx, ]
  if (!is.null(exclude_id)) {
    candidates <- candidates[candidates$id != exclude_id, ]
  }
  if (nrow(candidates) == 0) {
    return(candidates)
  }

  distance <- abs(candidates$mz - peak_mz)
  window <- safety_margin * pmax(candidates$fwhm_mz, peak_fwhm_mz)
  hits <- candidates[distance <= window, ]
  if (nrow(hits) > 0) {
    hits$distance <- distance[distance <= window]
  }
  hits
}

#' Search the m/z collision index for every reference-proteome peak that
#' collides with any of a target proteoform's own predicted charge-state
#' peaks.
#'
#' @param target proteoform object
#' @param mz_index m/z collision index built in the same mode as `mode`
#' @param mode "denatured" or "native"
#' @param average use average mass instead of monoisotopic
#' @param r_ref,mz_ref Orbitrap resolving-power settings
#' @param safety_margin multiplier applied to the FWHM window
#' @param exclude_id id to exclude from results (defaults to the target's own id)
#' @return data.frame: one row per colliding (target peak, reference peak) pair
search_mz_collisions <- function(target, mz_index, mode = c("denatured", "native"),
                                  average = FALSE, r_ref = 120000, mz_ref = 200,
                                  safety_margin = 1.0, exclude_id = target$id) {
  mode <- match.arg(mode)
  if (!inherits(target, "proteoform")) {
    stop("search_mz_collisions() requires a proteoform object")
  }

  mass <- proteoform_mass(target, average = average)$mass
  target_peaks <- predict_charge_envelope(target$sequence, mass, mode = mode)
  target_peaks$fwhm_mz <- fwhm_mz(target_peaks$mz, r_ref = r_ref, mz_ref = mz_ref)

  results <- lapply(seq_len(nrow(target_peaks)), function(i) {
    hits <- query_mz_collisions_for_peak(
      mz_index, target_peaks$mz[i], target_peaks$fwhm_mz[i],
      safety_margin = safety_margin, exclude_id = exclude_id
    )
    if (nrow(hits) == 0) {
      return(NULL)
    }
    hits$target_z <- target_peaks$z[i]
    hits$target_mz <- target_peaks$mz[i]
    hits
  })
  results <- results[!vapply(results, is.null, logical(1))]

  if (length(results) == 0) {
    return(data.frame(
      id = character(0), z = integer(0), mz = numeric(0),
      relative_intensity = numeric(0), fwhm_mz = numeric(0), distance = numeric(0),
      target_z = integer(0), target_mz = numeric(0)
    ))
  }
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}
