test_that("master time axis retains exact tree and fossil events", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "((a:5,b:5):5,c:10);")
  out <- hee_master_time_axis(c(10, 5, 0), tree = tree,
                              fossil_times = 7.5, additional_times = 2.5)
  expect_true(all(c(10, 7.5, 5, 2.5, 0) %in% out$time_ma))
  expect_true(out$is_tree_node[match(5, out$time_ma)])
  row <- out[match(7.5, out$time_ma), ]
  expect_equal(row$earth_slice_older_ma, 10)
  expect_equal(row$earth_slice_younger_ma, 5)
  expect_false(row$earth_slice_exact)
})
