# Parses the user-facing PTM specification text box into proteoform() PTM
# lists, reusing unimod_ptm() (R/unimod_table.R) for mass deltas so a parsed
# spec produces exactly the same proteoform object generate_fragment_ladder()
# and fragmentation_propensity() already know how to handle -- no separate
# "apply PTM to ladder" logic needed here.
#
# Grammar (one isoform's PTM text box):
#   "133_Ser_Phospho; 210_Thr_Phospho,215_Tyr_Sulfo"
#   semicolon = separate proteoforms derived from this isoform
#   comma     = multiple PTMs co-occurring on one proteoform (capped, see
#               PTM_SPEC_MAX_PER_PROTEOFORM)
#   each PTM  = "<1-based residue>_<3-letter AA code>_<Unimod name>"

PTM_SPEC_MAX_PER_PROTEOFORM <- 5

`%||%` <- function(a, b) if (is.null(a)) b else a

.AA3_TO_1 <- c(
  Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C", Gln = "Q", Glu = "E",
  Gly = "G", His = "H", Ile = "I", Leu = "L", Lys = "K", Met = "M", Phe = "F",
  Pro = "P", Ser = "S", Thr = "T", Trp = "W", Tyr = "Y", Val = "V"
)

#' Parse one isoform's PTM spec text box into a list of PTM groups, each a
#' list of ptm() objects ready to pass as proteoform(ptms = ...).
#'
#' Invalid groups (bad format, out-of-range residue, AA/sequence mismatch,
#' unknown PTM name, or more than PTM_SPEC_MAX_PER_PROTEOFORM entries) are
#' dropped with a human-readable warning rather than erroring, so one bad
#' group doesn't block the other valid ones in the same text box.
#'
#' @param sequence the isoform's amino acid sequence, for AA validation
#' @param spec_text raw text box contents (possibly empty/blank)
#' @param label display label for this isoform, used only in warning text
#' @return list(groups = list of ptm-object-lists, warnings = character vector)
parse_ptm_spec_text <- function(sequence, spec_text, label = "isoform") {
  groups_out <- list()
  warnings_out <- character(0)
  if (is.null(spec_text) || !nzchar(trimws(spec_text %||% ""))) {
    return(list(groups = groups_out, warnings = warnings_out))
  }

  groups_raw <- trimws(strsplit(spec_text, ";")[[1]])
  groups_raw <- groups_raw[nzchar(groups_raw)]

  for (gi in seq_along(groups_raw)) {
    specs <- trimws(strsplit(groups_raw[gi], ",")[[1]])
    specs <- specs[nzchar(specs)]
    if (length(specs) > PTM_SPEC_MAX_PER_PROTEOFORM) {
      warnings_out <- c(warnings_out, sprintf(
        "%s group %d: %d PTMs exceeds the limit of %d -- group skipped.",
        label, gi, length(specs), PTM_SPEC_MAX_PER_PROTEOFORM
      ))
      next
    }

    ptms <- list()
    ok <- TRUE
    for (spec in specs) {
      parts <- strsplit(spec, "_")[[1]]
      if (length(parts) != 3) {
        warnings_out <- c(warnings_out, sprintf(
          '%s group %d: "%s" does not match residue_AA_PTM format -- group skipped.',
          label, gi, spec
        ))
        ok <- FALSE
        break
      }
      pos <- suppressWarnings(as.integer(parts[1]))
      aa3 <- parts[2]
      ptm_name <- parts[3]

      if (is.na(pos) || pos < 1 || pos > nchar(sequence)) {
        warnings_out <- c(warnings_out, sprintf(
          "%s group %d: residue %s out of range (1-%d) -- group skipped.",
          label, gi, parts[1], nchar(sequence)
        ))
        ok <- FALSE
        break
      }
      aa3_norm <- paste0(toupper(substr(aa3, 1, 1)), tolower(substr(aa3, 2, nchar(aa3))))
      aa1 <- .AA3_TO_1[[aa3_norm]]
      if (is.null(aa1)) {
        warnings_out <- c(warnings_out, sprintf(
          '%s group %d: unknown residue code "%s" -- group skipped.', label, gi, aa3
        ))
        ok <- FALSE
        break
      }
      actual <- substr(sequence, pos, pos)
      if (actual != aa1) {
        warnings_out <- c(warnings_out, sprintf(
          "%s group %d: residue %d is %s in this isoform, not %s -- group skipped.",
          label, gi, pos, actual, aa3
        ))
        ok <- FALSE
        break
      }
      ptm_obj <- tryCatch(
        unimod_ptm(ptm_name, site = pos),
        error = function(e) NULL
      )
      if (is.null(ptm_obj)) {
        warnings_out <- c(warnings_out, sprintf(
          '%s group %d: unknown PTM "%s" -- group skipped.', label, gi, ptm_name
        ))
        ok <- FALSE
        break
      }
      ptms[[length(ptms) + 1]] <- ptm_obj
    }
    if (ok && length(ptms) > 0) groups_out[[length(groups_out) + 1]] <- ptms
  }

  list(groups = groups_out, warnings = warnings_out)
}
