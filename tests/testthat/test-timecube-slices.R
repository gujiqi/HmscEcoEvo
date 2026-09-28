test_that("timecube slice writer preserves complete per-time grids", {
  cube <- data.frame(
    cell_id = rep(c("c1", "c2"), 2),
    time_ma = rep(c(10, 0), each = 2),
    lon = c(0.5, 1.5, 0.5, 1.5),
    temperature = c(20, 21, 18, 19),
    stringsAsFactors = FALSE
  )
  out_dir <- tempfile("hee_timecube_slices_")
  index <- hee_write_timecube_slices(
    cube, out_dir, keep_cols = names(cube), key_cols = "cell_id",
    prefix = "earth", compression = FALSE
  )
  expect_equal(nrow(index), 2L)
  expect_true(all(index$n_rows == 2L))
  expect_true(all(file.exists(index$file)))
  expect_equal(readRDS(index$file[index$time_ma == 10])$temperature, c(20, 21))
  expect_true(file.exists(file.path(out_dir, "earth_index.csv")))
})

test_that("timecube slice writer rejects duplicate cell-time keys", {
  cube <- data.frame(cell_id = c("c1", "c1"), time_ma = c(0, 0), x = 1:2)
  expect_error(hee_write_timecube_slices(cube, tempfile("hee_duplicate_")),
               "Duplicate key")
})
