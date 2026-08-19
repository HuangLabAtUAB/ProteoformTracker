# search_confounding_proteins_combined() needs a real target proteoform's
# actual computed mass (to construct realistic mass-window and m/z-collision
# candidates relative to it), so these are gated on the Python mass engine
# unlike test-mz_collision_index.R's lower-level tests, which bypass
# proteoform_mass() entirely by predicting peaks from a hand-picked mass.
basic_rich_sequence <- function(n_repeats = 15) {
  paste(rep("AGSTVLINQFDEHWYCPMK", n_repeats), collapse = "")
}

test_that("search_confounding_proteins_combined finds mass-window and m/z-collision candidates, and merges a candidate found by both", {
  skip_if_not(mass_engine_available, "Python mass engine not available")

  seq <- basic_rich_sequence()
  target <- proteoform(id = "target", sequence = seq, provenance = "manual")
  target_mass <- proteoform_mass(target)$mass

  target_peaks <- predict_charge_envelope(seq, target_mass, mode = "denatured")
  # a genuinely different sequence so it doesn't ALSO turn up as a mass-window
  # hit by coincidence of identical composition
  other_seq <- basic_rich_sequence(14)

  mass_window_mass <- target_mass + 0.01 # well within any realistic FWHM window
  far_colliding_mass <- target_mass * 6 / 8 # collides with target's z=8 peak at z=6
  unrelated_mass <- target_mass + 50000 # far in both mass and m/z

  mass_index <- data.frame(
    id = c("close_in_mass", "far_but_colliding", "unrelated"),
    sequence = c(other_seq, other_seq, other_seq),
    length = nchar(other_seq),
    mass = c(mass_window_mass, far_colliding_mass, unrelated_mass),
    stringsAsFactors = FALSE
  )
  mass_index <- mass_index[order(mass_index$mass), ] # query_confounding_proteins() requires sorted mass
  mz_index <- build_reference_mz_index(mass_index, mode = "denatured")

  result <- search_confounding_proteins_combined(target, mass_index, mz_index = mz_index, mode = "denatured")
  cands <- result$candidates

  # "close_in_mass" is ~0.01 Da from the target, so its own charge-state
  # envelope naturally collides at essentially every charge state too --
  # "both" is the physically correct outcome here, not "mass" alone
  expect_true("close_in_mass" %in% cands$id)
  expect_true(cands$found_via[cands$id == "close_in_mass"] %in% c("mass", "both"))

  # the real point of this test: "far_but_colliding" is 5000+ Da from the
  # target -- far outside any resolvable mass window -- yet still turns up,
  # found ONLY via its m/z collision, which a mass-domain-only search
  # (the app's behavior before this change) would never have surfaced
  expect_true("far_but_colliding" %in% cands$id)
  expect_equal(cands$found_via[cands$id == "far_but_colliding"], "mz")
  expect_gt(cands$n_colliding_peaks[cands$id == "far_but_colliding"], 0)

  expect_false("unrelated" %in% cands$id)
})

test_that("search_confounding_proteins_combined drops an m/z-colliding id absent from a FILTERED mass_index, instead of crashing", {
  # regression test: server.R's top-down mass filter passes a SUBSET of the
  # full reference proteome as mass_index, while mz_index is built from the
  # unfiltered full proteome -- a real collision from an id outside that
  # filtered subset must be dropped gracefully (no sequence to build a
  # proteoform from), not produce an NA id that crashes proteoform() downstream
  skip_if_not(mass_engine_available, "Python mass engine not available")

  seq <- basic_rich_sequence()
  target <- proteoform(id = "target", sequence = seq, provenance = "manual")
  target_mass <- proteoform_mass(target)$mass
  other_seq <- basic_rich_sequence(14)
  far_colliding_mass <- target_mass * 6 / 8

  full_mass_index <- data.frame(
    id = c("in_filtered_pool", "outside_filtered_pool"),
    sequence = c(other_seq, other_seq),
    length = nchar(other_seq),
    mass = c(target_mass + 0.01, far_colliding_mass),
    stringsAsFactors = FALSE
  )
  full_mass_index <- full_mass_index[order(full_mass_index$mass), ]
  mz_index <- build_reference_mz_index(full_mass_index, mode = "denatured")

  # the FILTERED index server.R would actually pass in -- excludes
  # "outside_filtered_pool" entirely, even though mz_index still knows about it
  filtered_mass_index <- full_mass_index[full_mass_index$id == "in_filtered_pool", ]

  result <- expect_no_error(
    search_confounding_proteins_combined(target, filtered_mass_index, mz_index = mz_index, mode = "denatured")
  )
  expect_false("outside_filtered_pool" %in% result$candidates$id)
  expect_false(any(is.na(result$candidates$id)))
})

test_that("search_confounding_proteins_combined caps candidates and prioritizes 'both' over mass-only over mz-only", {
  skip_if_not(mass_engine_available, "Python mass engine not available")

  seq <- basic_rich_sequence()
  target <- proteoform(id = "target", sequence = seq, provenance = "manual")
  target_mass <- proteoform_mass(target)$mass
  other_seq <- basic_rich_sequence(14)

  # 5 candidates, all within a trivially-close mass window (real
  # mass-domain hits), so found_via == "mass" for all of them without
  # needing to construct real m/z collisions
  mass_index <- data.frame(
    id = paste0("cand", 1:5),
    sequence = other_seq,
    length = nchar(other_seq),
    mass = target_mass + seq(0.001, 0.005, length.out = 5),
    stringsAsFactors = FALSE
  )
  mass_index <- mass_index[order(mass_index$mass), ]

  full_result <- search_confounding_proteins_combined(target, mass_index, mz_index = NULL, mode = "denatured", max_candidates = 1000)
  expect_gte(nrow(full_result$candidates), 5) # sanity: the window really does catch all 5

  capped <- search_confounding_proteins_combined(target, mass_index, mz_index = NULL, mode = "denatured", max_candidates = 2)
  expect_equal(nrow(capped$candidates), 2)
  expect_equal(capped$n_candidates_before_cap, nrow(full_result$candidates))
})

test_that("search_confounding_proteins_combined falls back to mass-only when mz_index is NULL", {
  skip_if_not(mass_engine_available, "Python mass engine not available")

  seq <- basic_rich_sequence()
  target <- proteoform(id = "target", sequence = seq, provenance = "manual")
  target_mass <- proteoform_mass(target)$mass
  other_seq <- basic_rich_sequence(14)

  mass_index <- data.frame(
    id = "close_in_mass", sequence = other_seq, length = nchar(other_seq),
    mass = target_mass + 0.01, stringsAsFactors = FALSE
  )

  result <- search_confounding_proteins_combined(target, mass_index, mz_index = NULL, mode = "denatured")
  expect_equal(result$candidates$found_via, "mass")
  expect_null(result$mz_collision_detail)
})
