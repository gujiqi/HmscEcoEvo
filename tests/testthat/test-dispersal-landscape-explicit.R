test_that("movement landscape separates hard geography from permeability", {
  cells <- data.frame(
    cell_id = c("land_a", "land_b", "ocean"),
    land_exists = c(1, 1, 0),
    permeability = c(1, 0.4, 1),
    cell_area_km2 = c(100, 200, 100)
  )
  out <- hee_dispersal_landscape_weights(
    cells, hard_mask_col = "land_exists",
    movement_weight_col = "permeability", cell_area_col = "cell_area_km2"
  )
  expect_equal(out$cells$movement_allowed, c(1L, 1L, 0L))
  expect_equal(out$cells$movement_landscape_weight, c(1, 0.4, 0))
  expect_equal(out$cells$destination_measure[3], 0)
  expect_equal(out$diagnostics$n_traversable_cells, 2)
})

test_that("spherical kernel uses great-circle geometry, area, and landscape weights", {
  cells <- data.frame(
    cell_id = c("source", "across_date_line", "far", "ocean"),
    lon = c(179, -179, 160, 175), lat = c(0, 0, 0, 0), time_ma = 10,
    land_exists = c(1, 1, 1, 0), permeability = c(1, 1, 0.5, 1),
    cell_area_km2 = c(1, 1, 1, 1)
  )
  landscape <- hee_dispersal_landscape_weights(
    cells, hard_mask_col = "land_exists",
    movement_weight_col = "permeability", cell_area_col = "cell_area_km2"
  )
  kernel <- hee_dispersal_spherical_kernel(
    landscape, delta_t_myr = 1,
    diffusion_variance_rad2_per_myr = 0.002,
    source_emigration_rate_per_myr = 0.4, max_distance_km = 3000
  )
  from_source <- kernel$edges[kernel$edges$from_cell_id == "source", ]
  expect_false(any(kernel$edges$to_cell_id == "ocean"))
  expect_equal(sum(from_source$directional_weight), 1, tolerance = 1e-12)
  expect_equal(sum(from_source$movement_rate_per_myr), 0.4, tolerance = 1e-12)
  expect_lt(
    from_source$distance_km[from_source$to_cell_id == "across_date_line"],
    from_source$distance_km[from_source$to_cell_id == "far"]
  )
  expect_gt(
    from_source$directional_weight[from_source$to_cell_id == "across_date_line"],
    from_source$directional_weight[from_source$to_cell_id == "far"]
  )
})

test_that("within-region spherical movement cannot cross an imposed BGB boundary", {
  cells <- data.frame(
    cell_id = c("a", "b", "c"), lon = c(0, 1, 2), lat = 0, time_ma = 5,
    region = c("R1", "R1", "R2"), movement_allowed = 1,
    movement_landscape_weight = 1, cell_measure_weight = 1
  )
  kernel <- hee_dispersal_spherical_kernel(
    cells, delta_t_myr = 1, diffusion_precision_myr_per_rad2 = 100,
    max_distance_km = 500, spatial_domain = "within_region"
  )
  expect_false(any(kernel$edges$from_cell_id == "a" & kernel$edges$to_cell_id == "c"))
  expect_true(any(kernel$edges$from_cell_id == "a" & kernel$edges$to_cell_id == "b"))
})

test_that("a documented topographic edge cost attenuates active movement once", {
  cells <- data.frame(
    cell_id = c("source", "easy", "steep"), lon = c(0, 1, 1), lat = 0,
    time_ma = 5, movement_allowed = 1, movement_landscape_weight = 1,
    cell_measure_weight = 1
  )
  edges <- data.frame(
    from_cell_id = c("source", "source"),
    to_cell_id = c("easy", "steep"),
    distance_km = c(100, 100),
    effective_cost_km = c(100, 500)
  )
  kernel <- hee_dispersal_spherical_kernel(
    cells, edges = edges, delta_t_myr = 1,
    diffusion_variance_rad2_per_myr = 0.002
  )
  source_edges <- kernel$edges[kernel$edges$from_cell_id == "source", ]
  expect_gt(
    source_edges$directional_weight[source_edges$to_cell_id == "easy"],
    source_edges$directional_weight[source_edges$to_cell_id == "steep"]
  )
  expect_equal(
    source_edges$effective_movement_distance_km[source_edges$to_cell_id == "steep"],
    500
  )
  expect_equal(
    kernel$metadata$edge_cost_col,
    "effective_cost_km"
  )
})

test_that("a time-indexed sparse graph is filtered to the active slice", {
  cells <- data.frame(
    cell_id = c("a", "b", "c"), lon = c(0, 1, 2), lat = 0,
    time_ma = 5, movement_allowed = 1
  )
  edges <- data.frame(
    from_cell_id = c("a", "a"), to_cell_id = c("b", "c"),
    time_ma = c(5, 10), distance_km = c(100, 200)
  )
  kernel <- hee_dispersal_spherical_kernel(
    cells, edges = edges, delta_t_myr = 1,
    diffusion_variance_rad2_per_myr = 0.002
  )
  expect_equal(nrow(kernel$edges), 1L)
  expect_equal(kernel$edges$to_cell_id, "b")
})

test_that("spherical edge kernel converts to a target-by-source movement matrix", {
  kernel <- structure(list(edges = data.frame(
    from_cell_id = c("a", "a", "b"), to_cell_id = c("b", "c", "a"),
    movement_rate_per_myr = c(0.3, 0.1, 0.4)
  )), class = "hee_dispersal_spherical_kernel")
  movement <- hee_dispersal_transition_matrix(kernel, c("a", "b", "c"), sparse = FALSE)
  expect_equal(movement["b", "a"], 0.3)
  expect_equal(movement["c", "a"], 0.1)
  expect_equal(sum(movement[, "a"]), 0.4)
})

test_that("a disconnected time slice returns a labelled zero movement matrix", {
  cells <- data.frame(
    cell_id = c("a", "b"), lon = c(0, 1), lat = c(0, 0), time_ma = 5,
    land_exists = 0
  )
  landscape <- hee_dispersal_landscape_weights(cells, hard_mask_col = "land_exists")
  kernel <- hee_dispersal_spherical_kernel(
    landscape, delta_t_myr = 1, diffusion_variance_rad2_per_myr = 0.002,
    max_distance_km = 500
  )
  movement <- hee_dispersal_transition_matrix(kernel, c("a", "b"), sparse = FALSE)
  expect_equal(dim(movement), c(2L, 2L))
  expect_equal(sum(movement), 0)
})

test_that("spacetime graph keeps plate carriage separate and path search is directed", {
  horizontal <- rbind(
    data.frame(from_cell_id = "a", to_cell_id = "b", time_ma = 10,
               effective_cost_km = 1),
    data.frame(from_cell_id = "b2", to_cell_id = "c", time_ma = 5,
               effective_cost_km = 1)
  )
  transport <- data.frame(
    source_cell_id = "b", target_cell_id = "b2", target_weight = 1,
    time_from_ma = 10, time_to_ma = 5
  )
  graph <- hee_dispersal_spacetime_graph(horizontal, transport)
  expect_equal(graph$diagnostics$n_plate_carriage_edges, 1)
  expect_false(graph$edges$active_dispersal[graph$edges$edge_kind == "plate_carriage"])
  path <- hee_dispersal_path_diagnostic(graph, "a", "c", 10, 5)
  expect_true(path$found)
  expect_equal(path$path_cost, 2)
  expect_equal(path$n_plate_carriage_edges, 1)
  no_carry <- hee_dispersal_path_diagnostic(
    graph, "a", "c", 10, 5, include_plate_carriage = FALSE
  )
  expect_false(no_carry$found)
})

test_that("explicit transition pruning produces a finite root posterior", {
  tree <- structure(list(
    edge = matrix(c(3, 1, 3, 2), ncol = 2, byrow = TRUE),
    tip.label = c("a", "b"), Nnode = 1L
  ), class = "phylo")
  K <- matrix(0.5, 2, 2, dimnames = list(c("x", "y"), c("x", "y")))
  result <- hee_dispersal_pruning_likelihood(
    tree,
    tip_likelihood = list(a = c(x = 1, y = 0), b = c(x = 0, y = 1)),
    branch_transition = list(a = K, b = K),
    root_prior = c(x = 0.5, y = 0.5)
  )
  expect_true(is.finite(result$log_likelihood))
  expect_equal(sum(result$root_posterior), 1, tolerance = 1e-12)
  expect_equal(unname(result$root_posterior), c(0.5, 0.5), tolerance = 1e-12)
})

test_that("canonical propagule sampler is a stable alias for the legacy name", {
  edges <- data.frame(
    from_cell_id = c("a", "a", "b"), to_cell_id = c("b", "c", "a"),
    K_movement = c(0.8, 0.2, 1)
  )
  expect_equal(
    hee_dispersal_propagule_particles(edges, particles_per_source = 20L, seed = 7),
    hee_dispersal_particle_kernel(edges, particles_per_source = 20L, seed = 7)
  )
})
