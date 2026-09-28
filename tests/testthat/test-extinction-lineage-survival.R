test_that("extant-lineage survival conditioning is the correct conditional probability", {
  q <- cbind(a = c(0.1, 0.2), b = c(0.5, 0))
  got <- hee_extinction_condition_extant_lineage_survival(q)
  p_a <- 1 - (1 - 0.1) * (1 - 0.2)
  p_b <- 0.5
  expect_equal(got$occupancy[, "a"], q[, "a"] / p_a)
  expect_equal(got$occupancy[, "b"], q[, "b"] / p_b)
  expect_equal(got$audit$survival_probability_before, c(p_a, p_b))
  expect_true(all(colSums(got$occupancy) >= 1))
})

test_that("conditioning can target known extant-tree branches only", {
  q <- cbind(a = c(0.1, 0.2), b = c(0.3, 0.4))
  got <- hee_extinction_condition_extant_lineage_survival(q, "a")
  expect_equal(got$occupancy[, "b"], q[, "b"])
  expect_false(got$audit$conditioned[got$audit$lineage == "b"])
})

test_that("an impossible all-zero lineage is reported rather than repaired", {
  q <- cbind(a = c(0, 0), b = c(0.3, 0.4))
  expect_error(
    hee_extinction_condition_extant_lineage_survival(q),
    "incompatible with the required extant-tree branch survival"
  )
})
