test_that("global grid dynamic occupancy uses arrival, establishment, persistence and hard land state", {
  suitability <- data.frame(
    lineage = "L1",
    cell_id = rep(c("c1", "c2", "c3"), 2),
    time_ma = rep(c(10, 0), each = 3),
    eta = c(2, 2, 2, 2, 2, 2),
    suitability = 0.9,
    H_state = c(1, 0, 1, 1, 1, 1),
    stringsAsFactors = FALSE
  )
  movement <- data.frame(
    lineage = "L1",
    from_cell_id = c("c1", "c1", "c2", "c2", "c3", "c3"),
    to_cell_id = c("c1", "c2", "c2", "c3", "c3", "c1"),
    time_ma = 0,
    K_movement = c(1, 0.5, 1, 0.5, 1, 0),
    stringsAsFactors = FALSE
  )
  root <- data.frame(lineage = "L1", cell_id = "c1", q_root = 1)

  out <- hee_global_grid_dynamic_occupancy(
    suitability,
    movement_kernel = movement,
    root_state = root,
    establishment_intercept = 0,
    establishment_slope = 0,
    persistence_intercept = 0,
    persistence_slope = 0,
    internal_dt = 10,
    summarise_draws = FALSE
  )

  draws <- out$draws
  expect_true(all(draws$probability_draw >= 0 & draws$probability_draw <= 1))
  expect_equal(draws$probability_draw[draws$cell_id == "c2" &
                                        draws$time_ma == 10], 0)
  expect_equal(draws$colonisation_probability[draws$cell_id == "c3" &
                                                draws$time_ma == 0], 0)
  expect_gt(draws$arrival_pressure[draws$cell_id == "c2" &
                                     draws$time_ma == 0], 0)
  expect_gt(draws$colonisation_probability[draws$cell_id == "c2" &
                                             draws$time_ma == 0], 0)
  expect_gt(draws$probability_draw[draws$cell_id == "c2" &
                                     draws$time_ma == 0], 0)
  expect_equal(draws$transition_case[draws$cell_id == "c2" &
                                       draws$time_ma == 0],
               "global_grid_ctmc_update")
})

test_that("global grid dynamic occupancy can summarize response draws with weights", {
  suitability <- data.frame(
    lineage = "L1",
    response_draw = c("d1", "d2", "d1", "d2"),
    response_weight = c(0.25, 0.75, 0.25, 0.75),
    cell_id = "c1",
    time_ma = c(10, 10, 0, 0),
    eta = 1,
    suitability = 0.8,
    H_state = 1,
    stringsAsFactors = FALSE
  )
  root <- data.frame(lineage = "L1", cell_id = "c1",
                     response_draw = c("d1", "d2"),
                     q_root = c(0.2, 0.6))

  out <- hee_global_grid_dynamic_occupancy(
    suitability,
    movement_kernel = NULL,
    root_state = root,
    persistence_intercept = 10,
    persistence_slope = 0,
    internal_dt = 10,
    summarise_draws = TRUE
  )

  expect_true(all(c("probability_mean", "probability_q025",
                    "probability_q975", "n_response_draws") %in%
                    names(out$summary)))
  p10 <- out$summary$probability_mean[out$summary$time_ma == 10]
  expect_equal(p10, 0.25 * 0.2 + 0.75 * 0.6, tolerance = 1e-8)
  expect_equal(out$formula_catalog$formula_id[1], "global_grid_ctmc_update")
})
