test_that("target-group HMSC design retains zero semantics and aggregates cells", {
  sites <- data.frame(site_id = c("a", "b", "c"),
                      lon = c(1, 3, 7), lat = c(1, 1, 1))
  y <- matrix(c(1, 0,
                0, 1,
                0, 1), nrow = 3, byrow = TRUE,
              dimnames = list(c("a", "b", "c"), c("sp1", "sp2")))
  grid <- data.frame(site_id = c("a", "b", "c", "d"),
                     lon = c(1, 3, 7, 11), lat = c(1, 1, 1, 1))
  records <- data.frame(selected_species_id = c("sp1", "sp2", "sp2"),
                        lon = c(1, 3, 7), lat = c(1, 1, 1))
  out <- hee_hmsc_target_group_design(
    sites, y, grid, records, grid_deg = 4, min_target_group_records = 1
  )
  expect_equal(nrow(out$Y), 2L)
  expect_equal(colnames(out$Y), c("sp1", "sp2"))
  expect_true(all(out$sites$target_group_record_count >= 1))
  expect_match(out$metadata$value[out$metadata$field == "zero_semantics"],
               "not_confirmed_absence")
})
