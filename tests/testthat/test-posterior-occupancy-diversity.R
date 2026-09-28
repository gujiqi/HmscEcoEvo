test_that("posterior occupancy diversity is calculated from binary draw states", {
  draws <- list(
    d1 = rbind(a = c(1, 0), b = c(1, 1)),
    d2 = rbind(a = c(1, 1), b = c(0, 1))
  )
  out <- hee_posterior_occupancy_diversity(draws)

  expect_equal(out$expected_sampled_lineage_richness_mean, c(1.5, 1.5))
  expect_equal(out$posterior_incidence_shannon_entropy_mean,
               c(log(2) / 2, log(2) / 2))
  expect_equal(out$posterior_incidence_gini_simpson_mean, c(0.25, 0.25))
  expect_equal(out$n_posterior_draws, c(2, 2))
})

test_that("posterior occupancy diversity refuses a mean probability field", {
  expect_error(
    hee_posterior_occupancy_diversity(list(matrix(c(.2, .8), ncol = 1))),
    "binary"
  )
})
