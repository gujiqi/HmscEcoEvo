test_that("reference thresholds are calculated from occupied response values", {
  y <- matrix(c(1, 0, 1, 0, 0, 1), nrow = 3,
              dimnames = list(NULL, c("a", "b")))
  eta <- matrix(c(1, -2, 3, -1, 0, 2), nrow = 3,
                dimnames = list(NULL, c("a", "b")))
  out <- hee_environmental_filtering_reference_threshold(
    y, eta, presence_quantile = 0
  )
  expect_equal(out$eta_reference[out$lineage == "a"], 1)
  expect_equal(out$eta_reference[out$lineage == "b"], 2)
  expect_equal(out$n_observed_presence, c(2, 1))
})

test_that("relative support is centred on the declared response reference", {
  eta <- matrix(c(-1, 0, 1, 2), nrow = 2,
                dimnames = list(NULL, c("a", "b")))
  out <- hee_environmental_filtering_relative_support(
    eta, c(a = 0, b = 1), scale = c(a = 1, b = 2)
  )
  expect_equal(unname(out$relative_support[1, "a"]), stats::plogis(-1))
  expect_equal(unname(out$relative_support[2, "a"]), 0.5)
  expect_equal(unname(out$relative_support[1, "b"]), stats::plogis(0))
  expect_equal(unname(out$relative_support[2, "b"]), stats::plogis(0.5))
  expect_equal(unname(out$binary_environmental_support[2, "a"]), 1L)
})

test_that("prevalence-matched references reproduce the training prevalence rank", {
  y <- matrix(c(1, 0, 0, 0, 1, 0, 1, 0), ncol = 1,
              dimnames = list(NULL, "a"))
  eta <- matrix(c(1, -2, 3, 0, -1, 4, 2, -3), ncol = 1,
                dimnames = list(NULL, "a"))
  out <- hee_environmental_filtering_reference_threshold(
    y, eta, method = "prevalence_matched"
  )
  expect_equal(out$training_prevalence, 3 / 8)
  expect_equal(out$reference_method, "prevalence_matched")
  expect_equal(out$reference_source,
               "modern_training_prevalence_matched_X_beta_quantile")
  expect_true(is.finite(out$eta_reference))
})

test_that("relative process gradients do not treat prevalence references as survival penalties", {
  eta <- cbind(a = c(-2, 0, 2, NA), b = c(4, 4, 4, 4))
  got <- hee_environmental_filtering_relative_gradient(
    eta, habitat_state = c(1, 1, 1, 0), clip = 3
  )
  expect_equal(got$process_eta[1:3, "a"], c(-1.349, 0, 1.349), tolerance = 1e-10)
  expect_equal(got$process_eta[1:3, "b"], c(0, 0, 0))
  expect_true(is.na(got$process_eta[4, "a"]))
  expect_equal(got$audit$n_finite_available_habitat, c(3L, 3L))
})
