test_that("plate-grid carriage remaps occupancy before active dispersal", {
  transport <- data.frame(
    source_cell_id = c("old_a", "old_b", "old_c"),
    target_cell_id = c("young_x", "young_x", "young_y"),
    target_weight = c(0.25, 0.75, 1),
    stringsAsFactors = FALSE
  )
  tr <- hee_plate_grid_transport(
    transport,
    source_cells = c("old_a", "old_b", "old_c"),
    target_cells = c("young_x", "young_y", "new_land")
  )
  expect_s3_class(tr, "hee_plate_grid_transport")
  expect_equal(tr$diagnostics$n_unmapped_target_cells, 1)
  q <- matrix(c(1, 0, 0.5), ncol = 1,
              dimnames = list(c("old_a", "old_b", "old_c"), "lin"))
  carried <- hee_plate_carry_occupancy(
    q, tr, target_cells = c("young_x", "young_y", "new_land")
  )
  expect_equal(carried$occupancy["young_x", "lin"], 0.25)
  expect_equal(carried$occupancy["young_y", "lin"], 0.5)
  expect_equal(carried$occupancy["new_land", "lin"], 0)
})

test_that("plate-grid carriage rejects non-normalised target mixtures", {
  bad <- data.frame(
    source_cell_id = c("a", "b"), target_cell_id = c("x", "x"),
    target_weight = c(0.7, 0.7)
  )
  expect_error(hee_plate_grid_transport(bad), "target weights exceed one")
})
