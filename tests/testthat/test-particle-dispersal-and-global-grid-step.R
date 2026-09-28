test_that("particle dispersal preserves the total source kernel and seed", {
  edges <- data.frame(
    lineage = "lin1",
    time_ma = 20,
    from_cell_id = c("a", "a", "a", "b", "b"),
    to_cell_id = c("b", "c", "d", "a", "c"),
    K_movement = c(0.5, 0.25, 0.25, 0.8, 0.2)
  )
  x <- hee_dispersal_particle_kernel(edges, particles_per_source = 200L,
                                     seed = 17)
  y <- hee_dispersal_particle_kernel(edges, particles_per_source = 200L,
                                     seed = 17)
  expect_equal(x, y)
  expect_true(all(x$particle_count >= 0L))
  expect_equal(sum(x$particle_count[x$from_cell_id == "a"]), 200L)
  expect_equal(sum(x$particle_count[x$from_cell_id == "b"]), 200L)
  expect_equal(sum(x$particle_kernel[x$from_cell_id == "a"]),
               sum(edges$K_movement[edges$from_cell_id == "a"]))
  expect_equal(sum(x$particle_kernel[x$from_cell_id == "b"]),
               sum(edges$K_movement[edges$from_cell_id == "b"]))
})

test_that("source-conserving rate kernel keeps total source movement bounded", {
  edges <- data.frame(
    lineage = "lin1", time_ma = 20,
    from_cell_id = c("a", "a", "b"), to_cell_id = c("b", "c", "a"),
    movement_weight = c(0.8, 0.2, 0.5),
    K_movement = c(0.8, 0.2, 0.5)
  )
  rates <- hee_dispersal_rate_kernel(edges, base_rate = 0.6)
  expect_equal(sum(rates$movement_rate[rates$from_cell_id == "a"]), 0.3)
  expect_equal(sum(rates$movement_rate[rates$from_cell_id == "b"]), 0.3)
  expect_true(all(rates$movement_rate >= 0))
  expect_true(all(rates$source_total_rate <= 0.6))
})

test_that("global grid step enforces arrival, habitat, and deterministic seeds", {
  q <- matrix(c(1, 0), ncol = 1,
              dimnames = list(c("source", "target"), "lin1"))
  movement <- matrix(c(0, 0.3, 0, 0), nrow = 2,
                     dimnames = list(rownames(q), rownames(q)))
  eta <- matrix(0, 2, 1, dimnames = dimnames(q))
  step <- hee_projection_global_grid_step(q, movement, eta, delta_t = 1)
  expect_equal(step$arrival_hazard["source", "lin1"], 0)
  expect_equal(step$arrival_hazard["target", "lin1"], 0.3)
  expect_equal(step$colonisation_hazard["source", "lin1"], 0)
  expect_true(all(step$occupancy_next >= 0 & step$occupancy_next <= 1))
  expect_equal(step$arrival_probability_delta_t["target", "lin1"],
               1 - exp(-0.3))
  quarter <- hee_projection_global_grid_step(
    q, movement, eta, delta_t = 0.25, reporting_interval_myr = 1
  )
  expect_equal(quarter$arrival_probability_delta_t["target", "lin1"],
               1 - exp(-0.3 * 0.25))
  expect_equal(quarter$arrival_probability_reporting_interval["target", "lin1"],
               1 - exp(-0.3))

  blocked <- hee_projection_global_grid_step(
    q, movement, eta, delta_t = 1, habitat_state = c(1, 0),
    colonisation_biotic = c(0, 0.2)
  )
  expect_equal(blocked$occupancy_next["target", "lin1"], 0)
  stochastic_a <- hee_projection_global_grid_step(
    q, movement, eta, delta_t = 1, stochastic = TRUE, seed = 99
  )
  stochastic_b <- hee_projection_global_grid_step(
    q, movement, eta, delta_t = 1, stochastic = TRUE, seed = 99
  )
  expect_equal(stochastic_a$occupancy_next, stochastic_b$occupancy_next)
  expect_true(all(stochastic_a$occupancy_next %in% c(0, 1)))

  half_step <- hee_projection_global_grid_step(
    q, movement, eta, delta_t = 0.5, persistence_reference_interval = 1
  )
  full_step <- hee_projection_global_grid_step(
    q, movement, eta, delta_t = 1, persistence_reference_interval = 1
  )
  expect_equal(half_step$local_loss_hazard, full_step$local_loss_hazard)
  expect_equal(half_step$persistence_reference_interval, 1)
})

test_that("global grid step accepts lineage-specific process baselines", {
  q <- matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2,
              dimnames = list(c("c1", "c2"), c("slow", "fast")))
  eta <- matrix(0, nrow = 2, ncol = 2, dimnames = dimnames(q))
  movement <- Matrix::Diagonal(2, x = 0)
  step <- hee_projection_global_grid_step(
    occupancy = q, movement_matrix = movement, eta = eta, delta_t = 1,
    establishment_intercept = c(slow = -3, fast = 3),
    persistence_intercept = c(slow = -3, fast = 3),
    establishment_slope = 0, persistence_slope = 0
  )

  expect_true(all(step$establishment_probability[, "fast"] >
                  step$establishment_probability[, "slow"]))
  expect_true(all(step$persistence_probability[, "fast"] >
                  step$persistence_probability[, "slow"]))
})

test_that("global grid step accepts a serialized sparse movement cache", {
  q <- matrix(c(1, 0), ncol = 1,
              dimnames = list(c("source", "target"), "lin1"))
  eta <- matrix(0, 2, 1, dimnames = dimnames(q))
  movement <- Matrix::sparseMatrix(i = 2, j = 1, x = 0.3,
                                   dims = c(2, 2),
                                   dimnames = list(rownames(q), rownames(q)))
  cache <- tempfile(fileext = ".rds")
  saveRDS(movement, cache)
  cached_movement <- readRDS(cache)
  step <- hee_projection_global_grid_step(q, cached_movement, eta,
                                          delta_t = 1)
  expect_equal(step$arrival_hazard["target", "lin1"], 0.3)
})

test_that("global grid source movement multipliers act before arrival", {
  q <- matrix(c(1, 1), ncol = 2,
              dimnames = list("source", c("blocked", "mobile")))
  movement <- matrix(0.5, nrow = 1, ncol = 1,
                     dimnames = list("source", "source"))
  eta <- matrix(0, nrow = 1, ncol = 2, dimnames = dimnames(q))
  step <- hee_projection_global_grid_step(
    q, movement, eta, delta_t = 1,
    movement_multiplier = c(blocked = 0, mobile = 2)
  )
  expect_equal(step$arrival_hazard["source", "blocked"], 0)
  expect_equal(step$arrival_hazard["source", "mobile"], 1)
})
