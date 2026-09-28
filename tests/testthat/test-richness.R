test_that("richness sums occurrence probabilities by cell and time", {
  proj <- data.frame(cell_id = c(1, 1, 2), time_ma = c(0, 0, 0),
                     species = c("sp1", "sp2", "sp1"),
                     probability = c(0.2, 0.7, 0.4))
  rich <- hee_richness(proj)
  expect_equal(rich$expected_richness[rich$cell_id == 1], 0.9)
  bin <- hee_richness(proj, threshold = 0.5)
  expect_equal(bin$binary_richness[bin$cell_id == 1], 1)
})

test_that("richness rejects inconsistent coordinates within one output key", {
  proj <- data.frame(cell_id = c(1, 1), time_ma = c(0, 0),
                     species = c("sp1", "sp2"),
                     probability = c(0.2, 0.7),
                     lon = c(10, 11), lat = c(0, 0))
  expect_error(hee_richness(proj), "Duplicate key")
})
