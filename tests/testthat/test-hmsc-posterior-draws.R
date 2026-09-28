test_that("HMSC Beta posterior extractor samples real retained chain states", {
  model <- list(
    covNames = c("(Intercept)", "temp", "precip"),
    spNames = c("sp1", "sp2"),
    postList = list(
      list(
        list(Beta = matrix(c(1, 2, 3, 11, 12, 13), nrow = 3)),
        list(Beta = matrix(c(4, 5, 6, 14, 15, 16), nrow = 3))
      ),
      list(
        list(Beta = matrix(c(7, 8, 9, 17, 18, 19), nrow = 3)),
        list(Beta = matrix(c(10, 20, 30, 40, 50, 60), nrow = 3))
      )
    )
  )
  out <- hee_hmsc_beta_posterior_draws(model, n_draws = 3, seed = 9)

  expect_equal(nrow(out), 6L)
  expect_setequal(unique(out$lineage), c("sp1", "sp2"))
  expect_equal(length(unique(out$response_draw)), 3L)
  expect_true(all(out$chain %in% c(1L, 2L)))
  expect_true(all(out$intercept == 0))
  expect_true(all(is.finite(out$hmsc_intercept_original)))
  expect_setequal(names(out), c("lineage", "response_draw", "chain", "iteration",
                                "response_weight", "intercept",
                                "hmsc_intercept_original", "temp", "precip"))
  expect_equal(sum(unique(out[, c("response_draw", "response_weight")])$response_weight), 1)
  expect_equal(attr(out, "n_available_posterior_samples"), 4L)
  expect_equal(attr(out, "n_selected_posterior_samples"), 3L)
})

test_that("HMSC Beta posterior extractor validates requested axes and draw count", {
  model <- list(
    covNames = c("(Intercept)", "temp"), spNames = "sp1",
    postList = list(list(list(Beta = matrix(c(1, 2), nrow = 2))) )
  )
  expect_error(hee_hmsc_beta_posterior_draws(model, n_draws = 2),
               "exceeds the number")
  expect_error(hee_hmsc_beta_posterior_draws(model, axes = "missing"),
               "missing from model")
})
