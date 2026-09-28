test_that("plate carrier grid separates process state from output interpolation", {
  skip_if_not_installed("RANN")
  carrier <- data.frame(track_id = c("a", "b"), lon = c(0, 2), lat = c(0, 0))
  grid <- data.frame(cell_id = c("x", "y", "ocean"), lon = c(0, 2, 30),
                     lat = c(0, 0, 0), H_state = c(1, 1, 0), climate = c(3, 4, 9))
  matched <- hee_plate_carrier_match_grid(carrier, grid, time_ma = 10,
                                          max_match_distance_km = 250)
  expect_true(all(matched$active_land))
  q <- matrix(c(.2, .8), ncol = 1, dimnames = list(c("a", "b"), "lin"))
  rendered <- hee_plate_carrier_to_grid(q, matched, grid, k = 1,
                                         max_output_distance_km = 250)
  expect_equal(unname(rendered$values[c("x", "y"), "lin"]), c(.2, .8))
  expect_true(all(rendered$coverage$output_covered))
  expect_true(is.na(rendered$values["ocean", "lin"]))
})

test_that("carrier interpolation supports additive map summaries", {
  skip_if_not_installed("RANN")
  carrier <- data.frame(track_id = c("a", "b"), paleo_lon = c(0, 1),
                        paleo_lat = c(0, 0), active_land = TRUE)
  grid <- data.frame(cell_id = c("x", "y"), lon = c(0, 1), lat = c(0, 0),
                     H_state = 1)
  values <- matrix(c(2, 4), ncol = 1,
                   dimnames = list(c("a", "b"), "expected_richness"))
  out <- hee_plate_carrier_to_grid(
    values, carrier, grid, max_output_distance_km = 200,
    require_probability = FALSE
  )
  expect_true(all(is.finite(out$values)))
  expect_true(all(out$values >= 2 & out$values <= 4))
})

test_that("cached carrier output weights reproduce direct output interpolation", {
  skip_if_not_installed("RANN")
  skip_if_not_installed("Matrix")
  carrier <- data.frame(track_id = c("a", "b"), paleo_lon = c(0, 1),
                        paleo_lat = c(0, 0), active_land = TRUE)
  grid <- data.frame(cell_id = c("x", "y", "ocean"), lon = c(0, 1, 4),
                     lat = c(0, 0, 0), H_state = c(1, 1, 0))
  values <- matrix(c(.25, .75), ncol = 1,
                   dimnames = list(c("a", "b"), "q"))
  cached <- hee_plate_carrier_output_weights(
    carrier, grid, k = 1, max_output_distance_km = 200
  )
  direct <- hee_plate_carrier_to_grid(
    values, carrier, grid, k = 1, max_output_distance_km = 200
  )
  via_cache <- hee_plate_carrier_apply_output_weights(values, cached)
  expect_equal(via_cache, direct$values)
  expect_true(is.na(via_cache["ocean", "q"]))
})
