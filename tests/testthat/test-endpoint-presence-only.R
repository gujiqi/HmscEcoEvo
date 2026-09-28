test_that("observation-design audit refuses to call occurrence-grid zeros absences", {
  y <- rbind(c(1, 0), c(0, 1), c(1, 0))
  out <- hee_audit_observation_design(y)

  expect_true(out$all_rows_have_focal_record)
  expect_false(out$bernoulli_absence_likelihood_eligible)
  expect_equal(out$observation_design,
               "target_occurrence_grid_non_detections_not_absences")

  surveyed <- hee_audit_observation_design(y, survey_complete = rep(TRUE, 3))
  expect_true(surveyed$bernoulli_absence_likelihood_eligible)
})

test_that("presence-only score rewards support at known record locations", {
  occurrence <- rbind(c(2, 0), c(0, 3), c(0, 0))
  good <- rbind(c(.9, .1), c(.1, .9), c(.1, .1))
  bad <- rbind(c(.1, .9), c(.9, .1), c(.1, .1))
  colnames(occurrence) <- colnames(good) <- colnames(bad) <- c("sp1", "sp2")
  rownames(occurrence) <- rownames(good) <- rownames(bad) <- c("a", "b", "c")

  good_score <- hee_endpoint_presence_only_score(
    good, occurrence, effort = c(1, 1, 1), data_role = "spatial_holdout"
  )
  bad_score <- hee_endpoint_presence_only_score(
    bad, occurrence, effort = c(1, 1, 1), data_role = "spatial_holdout"
  )
  expect_equal(good_score$summary$endpoint_type,
               "presence_only_conditional_spatial_support")
  expect_gt(good_score$summary$mean_log_support_lift_over_uniform,
            bad_score$summary$mean_log_support_lift_over_uniform)
  expect_equal(sum(good_score$taxa$n_occurrence_records_raw), 5)
  expect_equal(sum(good_score$taxa$n_occurrence_units), 2)
})

test_that("formal presence-only scoring requires effort and defaults to binary grids", {
  occurrence <- rbind(c(3, 0), c(0, 4), c(0, 0))
  probability <- rbind(c(.9, .1), c(.1, .9), c(.1, .1))

  expect_error(
    hee_endpoint_presence_only_score(
      probability, occurrence, data_role = "independent_validation"
    ),
    "requires a declared effort"
  )
  out <- hee_endpoint_presence_only_score(
    probability, occurrence, effort = c(1, 2, 1),
    data_role = "independent_validation"
  )
  expect_equal(out$summary$occurrence_representation, "binary_grid")
  expect_equal(out$summary$n_occurrence_records_raw, 7)
  expect_equal(out$summary$n_occurrence_units, 2)
})
