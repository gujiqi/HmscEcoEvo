test_that("time axis reads the complete cube axis", {
  cube <- hee_load_timecube(object = make_test_timecube_540(times = c(0, 252, 540)))
  expect_equal(hee_time_axis_from_env_cube(cube), c(540, 252, 0))
  expect_equal(hee_time_axis_from_env_cube(cube, decreasing = FALSE), c(0, 252, 540))
})

test_that("full mode rejects a two-slice time axis", {
  expect_error(hee_assert_full_time_axis(c(540, 0), quick = FALSE),
               "Full mode requires")
  expect_silent(hee_assert_full_time_axis(c(540, 0), quick = TRUE))
})
