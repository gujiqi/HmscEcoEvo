test_that("presence-only endpoint particles are normalized without absences", {
  occurrence <- rbind(c(1, 0), c(0, 1), c(0, 0))
  rownames(occurrence) <- c("a", "b", "c")
  colnames(occurrence) <- c("x", "y")
  focused <- rbind(c(0.9, 0.1), c(0.1, 0.9), c(0.1, 0.1))
  rownames(focused) <- rownames(occurrence)
  colnames(focused) <- colnames(occurrence)
  diffuse <- matrix(0.5, 3, 2, dimnames = dimnames(occurrence))

  weights <- hee_endpoint_presence_only_particle_weights(
    list(focused = focused, diffuse = diffuse), occurrence,
    effort = c(1, 1, 1), data_role = "spatial_holdout",
    particle_representation = "probability_map"
  )
  expect_equal(sum(weights$weight), 1)
  expect_gt(weights$weight[weights$particle_id == "focused"], 0.5)
  expect_true(all(weights$scientific_role ==
                    "scenario_probability_map_reweighting_not_particle_posterior"))
})

test_that("formal particle reweighting requires binary terminal occupancy states", {
  occurrence <- matrix(c(1, 0), ncol = 1,
                       dimnames = list(c("a", "b"), "x"))
  binary_particles <- list(
    present = matrix(c(1, 0), ncol = 1, dimnames = dimnames(occurrence)),
    broad = matrix(c(1, 1), ncol = 1, dimnames = dimnames(occurrence))
  )
  out <- hee_endpoint_presence_only_particle_weights(
    binary_particles, occurrence, effort = c(1, 1),
    data_role = "spatial_holdout"
  )
  expect_true(all(out$scientific_role ==
                    "presence_only_endpoint_conditioned_binary_particle_weight"))
  expect_error(
    hee_endpoint_presence_only_particle_weights(
      list(not_a_state = matrix(c(.8, .2), ncol = 1)), occurrence,
      effort = c(1, 1), data_role = "spatial_holdout"
    ),
    "not a binary terminal occupancy state"
  )
})

test_that("training role remains an explicitly non-posterior diagnostic", {
  occurrence <- matrix(c(1, 0), ncol = 1,
                       dimnames = list(c("a", "b"), "x"))
  particles <- list(one = matrix(c(0.8, 0.2), ncol = 1,
                                 dimnames = dimnames(occurrence)))
  out <- hee_endpoint_presence_only_particle_weights(
    particles, occurrence, data_role = "training_diagnostic",
    particle_representation = "probability_map"
  )
  expect_match(out$scientific_role,
               "diagnostic_only_not_an_independent_posterior")
})
