test_that("elevation-derived resistance uses local relief rather than absolute elevation", {
  landscape <- expand.grid(
    lon = c(0, 1, 2), lat = c(0, 1, 2), time_ma = 10,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  landscape$cell_id <- paste0("c", seq_len(nrow(landscape)))
  landscape$land_mask_dem <- 1
  landscape$elevation_m <- 1000
  out <- hee_dispersal_resistance_surface(
    landscape, resolution_deg = 1, relief_ref_m = 500
  )
  expect_true(all(out$cells$movement_allowed == 1L))
  expect_true(all(out$cells$local_relief_m == 0))
  expect_true(all(out$cells$topographic_resistance == 0))
  expect_true(all(out$cells$n_terrain_neighbours > 0))
})

test_that("ocean has no local movement and elevation steps raise edge cost", {
  landscape <- data.frame(
    cell_id = c("a", "b", "c", "o"), time_ma = 10,
    lon = c(0, 1, 2, 1), lat = c(0, 0, 0, 1),
    elevation_m = c(100, 150, 2000, NA),
    land_mask_dem = c(1, 1, 1, 0)
  )
  resistance <- hee_dispersal_resistance_surface(
    landscape, resolution_deg = 1, relief_ref_m = 500
  )
  expect_equal(resistance$cells$movement_allowed[4], 0L)
  expect_true(is.na(resistance$cells$topographic_resistance[4]))
  graph <- hee_dispersal_connectivity_graph(
    resistance, resolution_deg = 1, step_ref_m = 500,
    weight_relief = 0.5, weight_step = 0.5, omega_topo = 1
  )
  expect_false(any(graph$edges$from_cell_id == "o" | graph$edges$to_cell_id == "o"))
  ab <- graph$edges[graph$edges$from_cell_id == "a" &
                      graph$edges$to_cell_id == "b", ]
  bc <- graph$edges[graph$edges$from_cell_id == "b" &
                      graph$edges$to_cell_id == "c", ]
  expect_equal(nrow(ab), 1)
  expect_equal(nrow(bc), 1)
  expect_gt(bc$effective_cost_km, ab$effective_cost_km)
})

test_that("zero topographic strength reduces effective cost to geographic distance", {
  landscape <- data.frame(
    cell_id = c("a", "b"), time_ma = 10, lon = c(0, 1), lat = 0,
    elevation_m = c(0, 2000), land_mask_dem = 1
  )
  resistance <- hee_dispersal_resistance_surface(landscape, resolution_deg = 1)
  graph <- hee_dispersal_connectivity_graph(
    resistance, resolution_deg = 1, omega_topo = 0
  )
  expect_equal(graph$edges$effective_cost_km, graph$edges$distance_km)
})

test_that("kernel and arrival distinguish movement allocation from arrival pressure", {
  edges <- data.frame(
    from_cell_id = c("a", "a"), to_cell_id = c("b", "c"), time_ma = 0,
    effective_cost_km = c(100, 500)
  )
  kernel <- hee_dispersal_edge_kernel(edges, dispersal_scale_km = 200,
                                      normalise_by_source = TRUE)
  expect_equal(sum(kernel$K_movement), 1)
  expect_gt(kernel$K_movement[kernel$to_cell_id == "b"],
            kernel$K_movement[kernel$to_cell_id == "c"])
  arrival <- hee_dispersal_arrival_pressure(
    data.frame(lineage = "global", cell_id = "a", q = 0.8), kernel
  )
  expect_true(all(arrival$raw_arrival > 0))
  expect_true(all(arrival$arrival_pressure > 0 & arrival$arrival_pressure < 1))
})
