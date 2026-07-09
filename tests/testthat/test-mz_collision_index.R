# These tests build a mass_index by hand (id, sequence, mass) rather than
# via build_reference_mass_index(), so they don't need the Python mass
# engine: build_reference_mz_index() reads mass from the index directly and
# only touches the sequence for basic-residue-count charge ceilings.

# 1 basic residue (K) per 19-char unit -> ceiling = n_repeats + 1 (N-term
# amine), giving a modest charge-state envelope comparable to a real ~16-20
# kDa protein (e.g. 15 repeats -> ceiling 16, envelope spans z=6..16).
basic_rich_sequence <- function(n_repeats = 15) {
  paste(rep("AGSTVLINQFDEHWYCPMK", n_repeats), collapse = "")
}

test_that("build_reference_mz_index produces a sorted, well-formed peak table", {
  seq <- basic_rich_sequence()
  mass_index <- data.frame(
    id = c("p1", "p2"), sequence = c(seq, seq),
    mass = c(16000, 20000), stringsAsFactors = FALSE
  )
  mz_index <- build_reference_mz_index(mass_index, mode = "denatured")

  expect_true(all(c("id", "z", "mz", "relative_intensity", "fwhm_mz") %in% names(mz_index)))
  expect_equal(mz_index$mz, sort(mz_index$mz))
  expect_true(all(c("p1", "p2") %in% mz_index$id))
})

test_that("query_mz_collisions_for_peak finds a colliding entry and excludes by id", {
  mz_index <- data.frame(
    id = c("a", "b", "c"), z = c(10L, 8L, 12L),
    mz = c(1000.00, 1000.02, 1500.00),
    relative_intensity = c(1, 1, 1),
    fwhm_mz = c(0.02, 0.02, 0.02),
    stringsAsFactors = FALSE
  )

  hits <- query_mz_collisions_for_peak(mz_index, peak_mz = 1000.00, peak_fwhm_mz = 0.02, safety_margin = 1.0)
  expect_true("b" %in% hits$id)
  expect_false("c" %in% hits$id)

  hits_excl <- query_mz_collisions_for_peak(
    mz_index,
    peak_mz = 1000.00, peak_fwhm_mz = 0.02, safety_margin = 1.0, exclude_id = "b"
  )
  expect_false("b" %in% hits_excl$id)
})

test_that("search_mz_collisions finds a real cross-mass charge-state collision that a mass-only search would miss", {
  # M_A/z_A = M_B/z_B exactly cancels the proton-mass term, so scaling a
  # target's mass by z_B/z_A creates a genuine, exact m/z collision at that
  # charge-state pair -- the same phenomenon found with the real O60384
  # example (masses in a ~4:3 ratio colliding at z=8/6, 12/9, 16/12).
  seq <- basic_rich_sequence()
  target_mass <- 20000
  colliding_mass <- 20000 * 6 / 8 # collides with target's z=8 peak at z=6

  # Bypass proteoform_mass()/pyteomics entirely: predict peaks directly from
  # a known mass, and build the mz_index the same way.
  target_peaks <- predict_charge_envelope(seq, target_mass, mode = "denatured")
  target_peaks$fwhm_mz <- fwhm_mz(target_peaks$mz)

  mass_index <- data.frame(
    id = c("far_but_colliding", "unrelated"),
    sequence = c(seq, seq),
    mass = c(colliding_mass, 45000),
    stringsAsFactors = FALSE
  )
  mz_index <- build_reference_mz_index(mass_index, mode = "denatured")

  peak_at_z8 <- target_peaks[target_peaks$z == 8, ]
  hits <- query_mz_collisions_for_peak(mz_index, peak_at_z8$mz, peak_at_z8$fwhm_mz, safety_margin = 1.0)

  expect_true("far_but_colliding" %in% hits$id)
  # the mass difference here (5000 Da) is far larger than any resolvable window --
  # a mass-domain-only search would never have surfaced this candidate
  expect_gt(abs(colliding_mass - target_mass), 100)
})
