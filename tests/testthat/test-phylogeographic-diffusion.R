test_that("location transitions retain probability mass without demographic rates", {
  edges <- data.frame(
    from_cell_id = c("a", "a", "b"),
    to_cell_id = c("b", "c", "a"),
    directional_weight = c(0.75, 0.25, 1),
    movement_rate_per_myr = c(0.3, 0.1, 0.2)
  )
  transition <- hee_dispersal_location_transition(
    edges, cell_order = c("a", "b", "c"), delta_t_myr = 1,
    sparse = FALSE
  )
  expect_equal(unname(base::colSums(transition$transition)), rep(1, 3), tolerance = 1e-12)
  expect_equal(
    transition$transition["a", "a"], exp(-0.4), tolerance = 1e-12
  )
  step <- hee_dispersal_location_distribution_step(
    matrix(c(a = 1, b = 0, c = 0), ncol = 1,
           dimnames = list(c("a", "b", "c"), "lineage_1")),
    transition
  )
  expect_equal(sum(step$location_mass), 1, tolerance = 1e-12)
  expect_equal(step$location_mass["b", "lineage_1"],
               (1 - exp(-0.4)) * 0.75, tolerance = 1e-12)
  expect_match(step$interpretation, "not colonisation")
})

test_that("terrestrial sink removes marine probability without renormalising it", {
  states <- c("land_a", "sea", "land_b")
  P <- matrix(c(0.5, 0.5, 0,
                0, 1, 0,
                0, 0.25, 0.75), 3, 3,
              dimnames = list(states, states))
  z <- hee_dispersal_terrestrial_sink_transition(
    P, source_land = c(TRUE, FALSE, TRUE),
    target_land = c(TRUE, FALSE, TRUE))
  expect_equal(unname(Matrix::colSums(z$transition)), rep(1, 4))
  expect_equal(as.numeric(z$transition["land_a", "land_a"]), 0.5)
  expect_equal(as.numeric(z$transition[z$sink_state, "land_a"]), 0.5)
  expect_equal(as.numeric(z$transition[z$sink_state, "sea"]), 1)
  expect_equal(as.numeric(z$transition["land_b", "land_b"]), 0.75)
  expect_error(hee_dispersal_terrestrial_sink_transition(
    P, c(TRUE, FALSE), c(TRUE, FALSE, TRUE)), "source_land")
})

test_that("rare jumps reach land hubs but never marine cells", {
  states <- c("land_a", "marine", "land_b")
  P <- diag(3); dimnames(P) <- list(states, states)
  out <- hee_dispersal_rare_land_jump_transition(
    P, c(TRUE, FALSE, TRUE), "land_b", jump_rate_per_myr = 0.1,
    delta_t_myr = 1)
  expect_equal(unname(Matrix::colSums(out)), rep(1, 3))
  expect_equal(as.numeric(out["land_b", "land_a"]), 1 - exp(-0.1))
  expect_equal(as.numeric(out["marine", "land_a"]), 0)
  expect_equal(as.numeric(out["marine", "marine"]), 1)
})

test_that("unavailable geography enters an explicit sink rather than persistence", {
  transition <- hee_dispersal_location_transition(
    data.frame(
      from_cell_id = "a", to_cell_id = "b", directional_weight = 1,
      movement_rate_per_myr = 1
    ),
    cell_order = c("a", "b"), delta_t_myr = 1,
    available_cells = c(a = 0, b = 1), sparse = FALSE
  )
  sink <- transition$sink_state
  expect_equal(transition$transition[sink, "a"], 1)
  step <- hee_dispersal_location_distribution_step(
    matrix(c(a = 1, b = 0, `__outside_available_geography__` = 0), ncol = 1,
           dimnames = list(rownames(transition$transition), "lineage_1")),
    transition
  )
  expect_equal(step$sink_mass[1, 1], 1)
  expect_error(
    hee_dispersal_location_distribution_step(
      matrix(c(a = 1, b = 0, `__outside_available_geography__` = 0), ncol = 1,
             dimnames = list(rownames(transition$transition), "lineage_1")),
      transition, condition_on_non_sink = TRUE
    ),
    "all mass entered"
  )
})

test_that("a submerged source can export before its retained mass enters the sink", {
  transition <- hee_dispersal_location_transition(
    data.frame(
      from_cell_id = "a", to_cell_id = "b", directional_weight = 1,
      movement_rate_per_myr = 1
    ),
    cell_order = c("a", "b"), delta_t_myr = 1,
    source_available_cells = c(a = 1, b = 1),
    target_available_cells = c(a = 0, b = 1), sparse = FALSE
  )
  sink <- transition$sink_state
  expect_equal(transition$transition["b", "a"], 1 - exp(-1), tolerance = 1e-12)
  expect_equal(transition$transition[sink, "a"], exp(-1), tolerance = 1e-12)
})

test_that("a sparse movement-rate matrix is accepted without rebuilding edges", {
  rate <- Matrix::sparseMatrix(
    i = c(2, 3, 1), j = c(1, 1, 2), x = c(0.3, 0.1, 0.2),
    dims = c(3, 3), dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
  )
  transition <- hee_dispersal_location_transition(rate, delta_t_myr = 1)
  expect_equal(unname(Matrix::colSums(transition$transition)), rep(1, 3), tolerance = 1e-12)
  expect_equal(transition$source_emigration_rate_per_myr[["a"]], 0.4)
})

test_that("time-ordered pruning matches a two-tip symmetric reconstruction", {
  tree <- structure(list(
    edge = matrix(c(3, 1, 3, 2), ncol = 2, byrow = TRUE),
    tip.label = c("a", "b"), Nnode = 1L
  ), class = "phylo")
  K <- matrix(0.5, 2, 2, dimnames = list(c("x", "y"), c("x", "y")))
  captured <- list()
  result <- hee_dispersal_time_ordered_pruning(
    tree,
    tip_likelihood = list(a = c(x = 1, y = 0), b = c(x = 0, y = 1)),
    branch_steps = list(
      a = list(list(transition = K, time_ma = 0)),
      b = list(list(transition = K, time_ma = 0))
    ),
    state_callback = function(branch_key, child_node, time_ma, stage_index,
                              posterior) {
      captured[[branch_key]] <<- posterior
    }
  )
  expect_true(is.finite(result$log_likelihood))
  expect_equal(unname(result$root_posterior), c(0.5, 0.5), tolerance = 1e-12)
  expect_equal(length(captured), 2L)
  expect_true(all(vapply(captured, function(x) abs(sum(x) - 1) < 1e-12,
                         logical(1))))
  expect_match(result$interpretation, "colonisation and persistence are deliberately absent")
})

test_that("plate carriage composes after horizontal diffusion", {
  D <- matrix(c(0.8, 0.2, 0.2, 0.8), 2, 2,
              dimnames = list(c("a", "b"), c("a", "b")))
  A <- matrix(c(0, 1, 1, 0), 2, 2,
              dimnames = list(c("a", "b"), c("a", "b")))
  result <- hee_dispersal_compose_plate_transport(D, A)
  expect_s3_class(result, "hee_dispersal_location_transition")
  expect_equal(as.numeric(result$transition), as.numeric(A %*% D))
  expect_equal(unname(colSums(result$transition)), c(1, 1), tolerance = 1e-12)
})
