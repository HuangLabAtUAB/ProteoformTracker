test_that("find_fragment_divergence_point locates the correct shared prefix/suffix across an insertion", {
  target_seq <- "AAAAAGGGGGSLLLLLVVVVV"
  iso_shifted_seq <- "AAAAAKKKGGGGGSLLLLLVVVVV" # 3-residue insertion after position 5

  target <- proteoform(id = "target", sequence = target_seq, provenance = "manual")
  iso_shifted <- proteoform(id = "iso_shifted", sequence = iso_shifted_seq, provenance = "manual")

  result <- find_fragment_divergence_point(target, iso_shifted)
  expect_equal(result$b_shared_length, 5)
  expect_equal(result$y_shared_length, 16)
  # a single contiguous insertion splits the sequence cleanly: shared
  # prefix + shared suffix should sum to exactly the target's full length
  expect_equal(result$b_shared_length + result$y_shared_length, nchar(target_seq))
})

test_that("find_fragment_divergence_point detects a same-length point substitution (no indel)", {
  # Real bug found validating against HBG1/HBG2 (P69891/P69892): they're the
  # same length, differing only by point substitutions, so the alignment
  # never introduces a gap (map[i] == i everywhere) -- b/y_shared_length
  # must still stop at the substitution, not just check alignment position.
  target_seq <- "AAAAAGGGGGSLLLLLVVVVV"
  substituted_seq <- "AAAAAGGGGGQLLLLLVVVVV" # S -> Q at position 11, same length

  target <- proteoform(id = "target", sequence = target_seq, provenance = "manual")
  substituted <- proteoform(id = "substituted", sequence = substituted_seq, provenance = "manual")

  result <- find_fragment_divergence_point(target, substituted)
  expect_equal(result$b_shared_length, 10)
  expect_equal(result$y_shared_length, 10)
})

test_that("find_fragment_divergence_point locates the correct shared prefix/suffix across a deletion", {
  target_seq <- "AAAAAGGGGGSLLLLLVVVVV"
  iso_deleted_seq <- "AAAAAGGGGGVVVVV" # the S/L block is spliced out

  target <- proteoform(id = "target", sequence = target_seq, provenance = "manual")
  iso_deleted <- proteoform(id = "iso_deleted", sequence = iso_deleted_seq, provenance = "manual")

  result <- find_fragment_divergence_point(target, iso_deleted)
  expect_equal(result$b_shared_length, 10)
  expect_equal(result$y_shared_length, 5)
})

test_that("find_fragment_divergence_point reports full-length sharing for identical sequences", {
  target_seq <- "AAAAAGGGGGSLLLLLVVVVV"
  target <- proteoform(id = "target", sequence = target_seq, provenance = "manual")
  identical_copy <- proteoform(id = "copy", sequence = target_seq, provenance = "manual")

  result <- find_fragment_divergence_point(target, identical_copy)
  expect_equal(result$b_shared_length, nchar(target_seq))
  expect_equal(result$y_shared_length, nchar(target_seq))
})

test_that("find_fragment_divergence_point requires proteoform objects", {
  p <- proteoform(id = "target", sequence = "AAAA", provenance = "manual")
  expect_error(find_fragment_divergence_point(p, list(sequence = "AAAA")), "requires proteoform")
})

test_that("fragment_mass_collision_check flags identical-ladder fragments as matched", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")
  seq <- "AAAAAGGGGGSLLLLLVVVVV"
  target <- proteoform(id = "target", sequence = seq, provenance = "manual")
  identical_copy <- proteoform(id = "copy", sequence = seq, provenance = "manual")

  result <- fragment_mass_collision_check(target, list(copy = identical_copy), script_path = TEST_PY_SCRIPT)
  expect_true(all(result$per_candidate$copy$matched))
})

test_that("fragment_mass_collision_check finds no matches against an unrelated sequence", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")
  target <- proteoform(id = "target", sequence = "AAAAAGGGGGSLLLLLVVVVV", provenance = "manual")
  unrelated <- proteoform(id = "unrelated", sequence = "WYFHKQNMRC", provenance = "manual")

  result <- fragment_mass_collision_check(target, list(unrelated = unrelated), script_path = TEST_PY_SCRIPT)
  expect_false(any(result$per_candidate$unrelated$matched))
})

test_that("fragment_mass_collision_check shows partial matches for an isoform that shares a prefix", {
  skip_if_not(mass_engine_available, "pyteomics/reticulate not available")
  target_seq <- "AAAAAGGGGGSLLLLLVVVVV"
  iso_deleted_seq <- "AAAAAGGGGGVVVVV"
  target <- proteoform(id = "target", sequence = target_seq, provenance = "manual")
  iso_deleted <- proteoform(id = "iso_deleted", sequence = iso_deleted_seq, provenance = "manual")

  result <- fragment_mass_collision_check(target, list(iso_deleted = iso_deleted), script_path = TEST_PY_SCRIPT)
  cand_result <- result$per_candidate$iso_deleted

  # early b-ions (within the shared AAAAAGGGGG prefix) should match
  expect_true(all(cand_result$matched[cand_result$ion_type == "b" & cand_result$cleavage_position <= 9]))
  # some fragment beyond the divergence point should NOT match
  expect_true(any(!cand_result$matched[cand_result$ion_type == "b" & cand_result$cleavage_position > 10]))
})

test_that("fragment_mass_collision_check requires a proteoform target and a list of proteoform candidates", {
  p <- proteoform(id = "target", sequence = "AAAA", provenance = "manual")
  expect_error(fragment_mass_collision_check(list(sequence = "AAAA"), list(p)), "requires a proteoform target")
  expect_error(fragment_mass_collision_check(p, list(list(sequence = "AAAA"))), "requires a list of proteoform candidates")
})
