test_that("forward plate transport conserves occupied area", {
  transport <- data.frame(
    source_cell_id = c("old_a", "old_b"),
    target_cell_id = c("young_x", "young_y"),
    source_weight = c(1, 1),
    source_area_km2 = c(2, 3),
    target_area_km2 = c(2, 3),
    stringsAsFactors = FALSE
  )
  tr <- hee_plate_forward_transport(
    transport,
    source_cells = c("old_a", "old_b"),
    target_cells = c("young_x", "young_y", "new_land")
  )
  expect_s3_class(tr, "hee_plate_forward_transport")
  q <- matrix(c(0.5, 1), ncol = 1,
              dimnames = list(c("old_a", "old_b"), "lin"))
  carried <- hee_plate_carry_forward_occupancy(
    q, tr, target_cells = c("young_x", "young_y", "new_land")
  )
  expect_equal(unname(carried$occupancy[, "lin"]), c(0.5, 1, 0))
  expect_equal(carried$source_expected_area_km2[["lin"]],
               carried$target_expected_area_km2[["lin"]])
  expect_true(carried$new_target_land[[3L]])
})

test_that("forward plate transport rejects undeclared source loss and overfill", {
  x <- data.frame(
    source_cell_id = "old_a", target_cell_id = "young_x",
    source_weight = 1, source_area_km2 = 1, target_area_km2 = 1
  )
  expect_error(
    hee_plate_forward_transport(x, source_cells = c("old_a", "old_b")),
    "mapped or explicitly declared"
  )
  expect_error(
    hee_plate_forward_transport(
      rbind(x, transform(x, source_cell_id = "old_b", source_area_km2 = 2)),
      source_cells = c("old_a", "old_b")
    ),
    "overfills"
  )
})
