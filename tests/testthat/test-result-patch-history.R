test_that("patches use terrestrial adjacency and preserve evidence types", {
  cells <- data.frame(cell_id = letters[1:5],
    support = c(.9, .8, .7, .9, .95),
    land_area_km2 = c(10, 10, 0, 10, 10),
    location_density = c(.3, .2, 0, .1, .4))
  edges <- data.frame(from_cell_id = c("a", "b", "c", "d"),
                      to_cell_id = c("b", "c", "d", "e"))
  out <- hee_result_patches(cells, edges, threshold = .6,
                            lineage_id = "L", time_ma = 10)
  expect_equal(nrow(out$patches), 2L)
  expect_equal(sort(out$patches$area_km2), c(20, 20))
  expect_equal(sum(out$patches$location_mass), 1)
  expect_false("c" %in% out$membership$cell_id)
  expect_error(hee_result_patches(cells, edges, threshold = 1.1,
                                   lineage_id = "L", time_ma = 10))
})

test_that("carriage links preserve splits and exclude unsupported cells", {
  old <- data.frame(patch_id = "old", cell_id = c("a", "b"),
    land_area_km2 = c(10, 10), time_ma = 10)
  young <- data.frame(patch_id = c("left", "right"),
    cell_id = c("x", "y"), time_ma = 5)
  tr <- data.frame(source_cell_id = c("a", "b"),
    target_cell_id = c("x", "y"), carriage_weight = 1)
  links <- hee_result_patch_links(old, young, tr)
  expect_equal(nrow(links), 2L)
  expect_equal(links$source_share, c(.5, .5))
  partial <- hee_result_patch_links(old, young[1L, , drop = FALSE],
    transform(tr[1L, , drop = FALSE], carriage_weight = .5))
  expect_equal(partial$source_share, .25)
  expect_error(hee_result_patch_links(old, young,
    rbind(tr, tr[1L, , drop = FALSE])), "at most one")
  expect_error(hee_result_patch_links(young, old, tr))
})

test_that("stability and contraction refugia do not imply occupancy", {
  p <- data.frame(patch_id = c("a", "b", "c"), time_ma = c(10, 5, 5),
    area_km2 = c(100, 30, 20), mean_support = c(.8, .9, .7),
    location_mass = c(.6, .2, .01))
  l <- data.frame(from_patch_id = c("a", "a"), to_patch_id = c("b", "c"),
    time_from_ma = 10, time_to_ma = 5, carried_area_km2 = c(30, 20),
    source_share = c(.3, .2))
  stable <- hee_result_patch_stability(p, l)
  expect_equal(nrow(stable$networks), 1L)
  expect_equal(stable$networks$n_patch_nodes, 3L)
  r <- hee_result_refugia_candidates(p, l, contraction_fraction = .25)
  expect_equal(nrow(r$episodes), 1L)
  expect_equal(nrow(r$candidates), 2L)
  expect_true(all(grepl("candidate", r$candidates$interpretation)))
  expect_setequal(r$candidates$geographic_evidence,
                  c("location_mass_above_threshold", "low_location_mass"))
})

test_that("corridors require joint edges and ignore plate carriage", {
  a <- data.frame(cell_id = c("a", "b"), patch_id = c("P", "Q"))
  b <- data.frame(cell_id = c("x", "y"), patch_id = c("R", "S"))
  f <- data.frame(from_cell_id = c("a", "a", "b"),
    to_cell_id = c("y", "x", "x"), time_from_ma = 10,
    time_to_ma = 5, posterior_edge_mass = c(.2, .7, .1),
    edge_kind = c("active_dispersal", "plate_carriage", "active_dispersal"))
  out <- hee_result_corridor_flow(f, a, b)
  expect_equal(sum(out$corridors$posterior_movement_mass), .3)
  expect_equal(out$unassigned_edge_mass, 0)
  same_id <- b; same_id$patch_id[[2L]] <- "P"
  same <- hee_result_corridor_flow(f, a, same_id)
  expect_equal(sum(same$corridors$posterior_movement_mass), .3)
  roles <- hee_result_donor_recipient(out$corridors, minimum_flow = .05)
  expect_true("dispersal_exporter" %in% roles$role)
  expect_true("dispersal_recipient" %in% roles$role)
  expect_error(hee_result_corridor_flow(f[, -5L], a, b),
               "marginal maps are insufficient")
})
