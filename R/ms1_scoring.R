# MS1 resolvability and envelope-crowding scoring.
# Ties together the resolving-power model (R/resolving_power.R) and charge-
# envelope prediction (R/charge_envelope.R) into the pairwise verdict and
# multi-proteoform crowding check described in the design spec, section 3.

DEFAULT_SAFETY_MARGIN <- 1.75 # spec recommends 1.5-2x over DeltaM_FWHM

#' MS1 resolvability verdict for a candidate pair (target vs. a relevant
#' isoform, or target vs. a confounding protein).
#'
#' Predicts the target's charge-state envelope, takes the best (smallest)
#' mass-domain FWHM across that envelope, and compares
#' best-case ΔM_FWHM x safety_margin against the pair's actual Δmass.
#'
#' Separately (per the design spec) flags `envelope_interleave_risk`: above
#' ~25-30 kDa, isotope envelopes themselves widen enough that even a
#' perfectly resolving instrument may not cleanly separate close
#' proteoforms. This is computed, not a hardcoded mass cutoff -- via
#' averagine_isotope_envelope() (R/isotope_envelope.R), risk is flagged
#' whenever the pair's actual Δmass is smaller than the width of the
#' isotope envelope itself: even though each individual isotope peak may be
#' instrumentally well-resolved, the two proteoforms' isotope combs overlap
#' peak-for-peak across most of their span, so a given peak can't be
#' unambiguously assigned to one proteoform vs the other. This can trip at
#' masses below 25-30 kDa for small-Δmass pairs, and not trip above it for
#' large-Δmass pairs -- consistent with the spec's caveat but not reducible
#' to a single mass threshold.
#'
#' @param target,candidate proteoform objects (see proteoform_schema.R)
#' @param mode "denatured" or "native"
#' @param average use average mass instead of monoisotopic
#' @param r_ref,mz_ref Orbitrap resolving-power settings
#' @param safety_margin multiplier applied to best-case ΔM_FWHM before
#'   comparing to actual Δmass (spec: 1.5-2x)
#' @return list with masses, delta_mass, best charge state/FWHM, verdict,
#'   isotope envelope stats, envelope_interleave_risk flag, and the full
#'   per-charge-state FWHM table
ms1_resolvability <- function(target, candidate, mode = c("denatured", "native"),
                               average = FALSE, r_ref = 120000, mz_ref = 200,
                               safety_margin = DEFAULT_SAFETY_MARGIN) {
  mode <- match.arg(mode)
  if (!inherits(target, "proteoform") || !inherits(candidate, "proteoform")) {
    stop("ms1_resolvability() requires proteoform objects")
  }

  target_mass <- proteoform_mass(target, average = average)$mass
  candidate_mass <- proteoform_mass(candidate, average = average)$mass
  delta_mass <- abs(target_mass - candidate_mass)

  envelope <- predict_charge_envelope(target$sequence, target_mass, mode = mode)
  fwhm_tbl <- fwhm_by_charge_state(target_mass, envelope$z, r_ref = r_ref, mz_ref = mz_ref)
  best_idx <- attr(fwhm_tbl, "best")
  best_fwhm <- fwhm_tbl$fwhm_mass[best_idx]
  best_z <- fwhm_tbl$z[best_idx]

  verdict <- if (delta_mass >= safety_margin * best_fwhm) {
    "resolvable"
  } else if (delta_mass >= best_fwhm) {
    "marginal"
  } else {
    "not-resolvable"
  }

  isotope_envelope <- averagine_isotope_envelope(target_mass)

  list(
    target_id = target$id,
    candidate_id = candidate$id,
    target_mass = target_mass,
    candidate_mass = candidate_mass,
    delta_mass = delta_mass,
    best_charge_state = best_z,
    best_fwhm_mass = best_fwhm,
    safety_margin = safety_margin,
    verdict = verdict,
    isotope_envelope_mean_shift = isotope_envelope$mean_shift,
    isotope_envelope_fwhm = isotope_envelope$fwhm,
    envelope_interleave_risk = delta_mass < isotope_envelope$fwhm,
    fwhm_table = fwhm_tbl,
    mode = mode
  )
}

#' Envelope-crowding check: predicted (m/z, charge) peaks from multiple
#' proteoforms, flagging any pair of peaks from *different* proteoforms
#' whose m/z falls within one resolvable window of each other. Run against
#' both the relevant-isoform set and the confounding-protein set.
#'
#' Unlike ms1_resolvability() (mass domain, one target vs. one candidate),
#' this compares directly in the m/z domain across every populated charge
#' state of every proteoform supplied -- the scenario a pairwise Δmass check
#' alone cannot see (different charge states of different-mass proteins
#' landing on the same m/z).
#'
#' @param proteoforms named list of proteoform objects (names used as ids
#'   in the output if the proteoform's own $id is not unique)
#' @param mode "denatured" or "native"
#' @param average use average mass instead of monoisotopic
#' @param r_ref,mz_ref Orbitrap resolving-power settings
#' @param safety_margin multiplier applied to the FWHM window before
#'   flagging a pair as crowded
#' @return list(peaks = data.frame of every predicted peak, flags =
#'   data.frame of flagged cross-proteoform peak pairs)
envelope_crowding_check <- function(proteoforms, mode = c("denatured", "native"),
                                     average = FALSE, r_ref = 120000, mz_ref = 200,
                                     safety_margin = 1.0) {
  mode <- match.arg(mode)
  if (!is.list(proteoforms) || !all(vapply(proteoforms, inherits, logical(1), "proteoform"))) {
    stop("envelope_crowding_check() requires a list of proteoform objects")
  }

  peaks <- do.call(rbind, lapply(proteoforms, function(p) {
    mass <- proteoform_mass(p, average = average)$mass
    envelope <- predict_charge_envelope(p$sequence, mass, mode = mode)
    envelope$id <- p$id
    envelope$mass <- mass
    envelope$fwhm_mz <- fwhm_mz(envelope$mz, r_ref = r_ref, mz_ref = mz_ref)
    envelope
  }))
  rownames(peaks) <- NULL

  n <- nrow(peaks)
  flags <- list()
  if (n >= 2) {
    for (i in seq_len(n - 1)) {
      for (j in seq((i + 1), n)) {
        if (peaks$id[i] == peaks$id[j]) next # only cross-proteoform crowding
        window <- safety_margin * max(peaks$fwhm_mz[i], peaks$fwhm_mz[j])
        distance <- abs(peaks$mz[i] - peaks$mz[j])
        if (distance <= window) {
          flags[[length(flags) + 1]] <- data.frame(
            id1 = peaks$id[i], z1 = peaks$z[i], mz1 = peaks$mz[i],
            id2 = peaks$id[j], z2 = peaks$z[j], mz2 = peaks$mz[j],
            distance = distance, window = window
          )
        }
      }
    }
  }

  flags_df <- if (length(flags) > 0) do.call(rbind, flags) else {
    data.frame(
      id1 = character(0), z1 = integer(0), mz1 = numeric(0),
      id2 = character(0), z2 = integer(0), mz2 = numeric(0),
      distance = numeric(0), window = numeric(0)
    )
  }
  rownames(flags_df) <- NULL

  list(peaks = peaks, flags = flags_df)
}
