test_that("tally_ms2_tiers counts each tier across both ion series", {
  tier_b <- c("unique", "common", "partial", "unique")
  tier_y <- c("common", "common", "unique", "neutral")
  result <- tally_ms2_tiers(tier_b, tier_y)

  expect_equal(result$total, 8)
  expect_equal(result$unique, 3)
  expect_equal(result$common, 3)
  expect_equal(result$partial, 1)
  expect_equal(result$neutral, 1)
})

test_that("tally_ms2_tiers handles the single-checked-proteoform (all-neutral) case", {
  result <- tally_ms2_tiers(rep("neutral", 5), rep("neutral", 5))
  expect_equal(result$total, 10)
  expect_equal(result$neutral, 10)
  expect_equal(result$unique, 0)
})

test_that("propensity_group applies the GLM thresholds (elevated: >, high: >=), no Very high tier", {
  scores <- c(0.5, 1.2, 1.200001, 4, 9.999, 10)
  expect_equal(
    propensity_group(scores, "glm"),
    c("Baseline", "Baseline", "Elevated", "High", "High", "High")
  )
})

test_that("propensity_group applies the RF thresholds, numerically distinct from GLM, no Very high tier", {
  scores <- c(0.01, 0.08, 0.081, 0.15, 0.29, 0.30)
  expect_equal(
    propensity_group(scores, "rf"),
    c("Baseline", "Baseline", "Elevated", "High", "High", "High")
  )
})

# ---- Fixtures: hand-built tier lists in the exact shape compute_ladder_tiers()
# produces, so build_section1/2_payload() can be tested without needing the
# mass-calculation engine (reticulate) or the RF model loaded. ----
.fake_pf <- function(id, sequence) proteoform(id = id, sequence = sequence, provenance = "manual")
.fake_tier <- function(n_bonds, propensity = rep(2, n_bonds), tier_b = rep("unique", n_bonds), tier_y = rep("common", n_bonds)) {
  list(
    ladder = list(b_mass = seq_len(n_bonds) * 100 + 0.5, y_mass = seq_len(n_bonds) * 90 + 0.3),
    propensity = propensity, tier_b = tier_b, tier_y = tier_y
  )
}

test_that("build_section1_payload carries transcript_id (iso_key) and a top-level gene_symbol", {
  pf_list <- list(PF1 = .fake_pf("PF1", "ACDEFG"), PF2 = .fake_pf("PF2", "GHIKLMN"))
  iso_key_of <- c(PF1 = "ENST00001", PF2 = "ENST00002")
  masses <- c(PF1 = 1000.1, PF2 = 1200.2)
  tiers <- list(PF1 = .fake_tier(5), PF2 = .fake_tier(6))

  payload <- build_section1_payload(pf_list, iso_key_of, masses, tiers, exon_table = NULL,
                                     scoring_mode = "glm", gene_symbol = "CD44")

  expect_equal(payload$gene_symbol, "CD44")
  expect_equal(payload$proteoforms[[1]]$transcript_id, "ENST00001")
  expect_equal(payload$proteoforms[[2]]$transcript_id, "ENST00002")
})

test_that("build_section2_payload keys tiers by target_id and gives each confounder its own ladder/tiers/gene", {
  target_pf <- .fake_pf("TARGET", "ACDEFGHIK")
  conf_pfs <- list(CONFA = .fake_pf("CONFA", "LMNPQRST"), CONFB = .fake_pf("CONFB", "VWYACDEFGH"))
  tiers <- list(TARGET = .fake_tier(8), CONFA = .fake_tier(7, tier_b = rep("partial", 7)), CONFB = .fake_tier(9))

  cands <- data.frame(
    id = c("CONFA", "CONFB"), sequence = c("LMNPQRST", "VWYACDEFGH"), length = c(8L, 10L),
    mass = c(5000, 6000), gene_symbol = c("GENEA", NA), found_via = c("mass", "both"),
    n_colliding_peaks = c(0L, 2L), stringsAsFactors = FALSE
  )
  confounder_search <- list(candidates = cands, window_da = 0.5, best_charge_state = 10, n_candidates_before_cap = 2)
  confounder_envs <- list(CONFA = NULL, CONFB = NULL)

  payload <- build_section2_payload(target_pf, 9000, tiers, "TARGET", target_exon_table = NULL,
                                     confounder_search = confounder_search, confounder_envs = confounder_envs,
                                     scoring_mode = "rf", gene_symbol = "CD44")

  expect_equal(payload$target$id, "TARGET")
  expect_equal(payload$target$gene_symbol, "CD44")
  expect_length(payload$target$tier_b, 8)

  expect_equal(payload$confounders[[1]]$id, "CONFA")
  expect_equal(payload$confounders[[1]]$gene_symbol, "GENEA")
  expect_equal(payload$confounders[[1]]$tier_b, rep("partial", 7))
  expect_true(is.na(payload$confounders[[2]]$gene_symbol))
  expect_length(payload$confounders[[2]]$tier_b, 9)
})

test_that("flatten_section1_peaks produces one row per MS1 point and per MS2 bond/ion, with propensity_group filled in", {
  payload <- list(
    gene_symbol = "CD44", scoring_mode = "glm",
    proteoforms = list(list(
      id = "PF1", label = "PF1", transcript_id = "ENST00001", mass = 1000.1, len = 6,
      env = list(list(z = 10L, resolved = TRUE, points = list(list(mz = 500.1, rel = 1.0), list(mz = 500.2, rel = 0.8)))),
      b_mass = c(100.1, 200.2), y_mass = c(90.1, 190.2), propensity = c(0.5, 2.0),
      tier_b = c("unique", "common"), tier_y = c("common", "unique")
    ))
  )
  df <- flatten_section1_peaks(payload)

  expect_equal(nrow(df), 2 + 2 * 2) # 2 MS1 points + 2 bonds x 2 ion types
  expect_equal(sum(df$peak_type == "MS1"), 2)
  expect_equal(sum(df$peak_type == "MS2"), 4)
  expect_true(all(df$gene_id == "CD44"))
  expect_true(all(df$transcript_id == "ENST00001"))

  ms2 <- df[df$peak_type == "MS2", ]
  expect_setequal(ms2$ion_type, c("b", "y"))
  # position-1 bond had propensity 0.5 -> Baseline in GLM mode; position-2 had 2.0 -> Elevated
  expect_equal(ms2$propensity_group[ms2$position == 1][1], "Baseline")
  expect_equal(ms2$propensity_group[ms2$position == 2][1], "Elevated")
})

test_that("flatten_section2_peaks labels target vs confounder rows and carries found_via/n_colliding_peaks only for confounders", {
  payload <- list(
    scoring_mode = "rf",
    target = list(
      id = "TARGET", label = NULL, gene_symbol = "CD44", mass = 9000, len = 3,
      env = list(), b_mass = c(1.0, 2.0), y_mass = c(1.5, 2.5), propensity = c(0.01, 0.5),
      tier_b = c("unique", "unique"), tier_y = c("unique", "unique")
    ),
    confounders = list(list(
      id = "CONFA", gene_symbol = "GENEA", mass = 5000, len = 2, found_via = "mass", n_colliding_peaks = 0L,
      env = list(), b_mass = c(3.0), y_mass = c(3.5), propensity = c(0.4), tier_b = c("partial"), tier_y = c("partial")
    ))
  )
  df <- flatten_section2_peaks(payload)

  expect_setequal(df$role, c("target", "confounder"))
  target_rows <- df[df$role == "target", ]
  conf_rows <- df[df$role == "confounder", ]
  expect_true(all(is.na(target_rows$found_via)))
  expect_true(all(is.na(target_rows$n_colliding_peaks)))
  expect_true(all(conf_rows$found_via == "mass"))
  expect_true(all(conf_rows$gene_id == "GENEA"))
})
