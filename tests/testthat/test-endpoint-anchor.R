test_that("endpoint nowcast anchor respects endpoints and dimensions", {
  forward <- matrix(c(0.1, 0.8, 0.4, 0.2), nrow = 2,
                    dimnames = list(c("c1", "c2"), c("sp1", "sp2")))
  nowcast <- matrix(c(0.3, 0.6, 0.5, 0.7), nrow = 2,
                    dimnames = dimnames(forward))
  expect_equal(hee_endpoint_nowcast_anchor(forward, nowcast, 0), forward)
  expect_equal(hee_endpoint_nowcast_anchor(forward, nowcast, 1), nowcast)
  expect_equal(
    hee_endpoint_nowcast_anchor(forward, nowcast, 0.25),
    0.75 * forward + 0.25 * nowcast
  )
  expect_error(hee_endpoint_nowcast_anchor(forward, nowcast[, 1], 0.5),
               "identical dimensions")
  expect_error(hee_endpoint_nowcast_anchor(forward, nowcast, 1.1), "\\[0, 1\\]")
})

test_that("prevalence-calibrated endpoint anchor uses the smallest valid weight", {
  forward <- matrix(c(0.8, 0.8, 0.8, 0.8), ncol = 1)
  nowcast <- matrix(c(0.2, 0.2, 0.2, 0.2), ncol = 1)
  out <- hee_endpoint_anchor_weight(forward, nowcast, c(0.9, 1.1))
  expect_equal(out$forward_to_endpoint_ratio, 4)
  expect_equal(out$selected_anchor_weight, 2.9 / 3, tolerance = 1e-12)
  expect_equal(out$anchored_to_endpoint_ratio, 1.1, tolerance = 1e-12)
  expect_true(out$tolerance_met)
  expect_equal(
    hee_endpoint_anchor_weight(nowcast, nowcast, c(0.9, 1.1))$selected_anchor_weight,
    0
  )
})
