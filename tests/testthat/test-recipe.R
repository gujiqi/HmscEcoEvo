test_that("recipe locks modern scaling and rejects missing variables", {
  X <- data.frame(bio1 = c(10, 12, 14), bio12 = c(800, 900, 1000),
                  row.names = paste0("site", 1:3))
  recipe <- hee_lock_recipe(X, ~ bio1 + bio12)
  Xp <- hee_apply_recipe(recipe, X)
  expect_equal(rownames(Xp), rownames(X))
  expect_equal(round(colMeans(Xp), 8), c(bio1 = 0, bio12 = 0))
  expect_error(hee_apply_recipe(recipe, data.frame(bio1 = 1)), "bio12")
})

test_that("recipe rejects unseen factor levels", {
  X <- data.frame(region = factor(c("A", "B", "A")))
  recipe <- hee_lock_recipe(X, ~ region)
  expect_error(hee_apply_recipe(recipe, data.frame(region = "C")), "unseen factor")
})

test_that("recipe rejects non-finite or all-missing numeric training variables", {
  expect_error(
    hee_lock_recipe(data.frame(bio1 = c(1, Inf, 3)), ~ bio1),
    "non-finite"
  )
  expect_error(
    hee_lock_recipe(data.frame(bio1 = c(NA_real_, NA_real_)), ~ bio1),
    "no finite values"
  )
  recipe <- hee_lock_recipe(data.frame(bio1 = c(2, NA_real_, 2)), ~ bio1)
  expect_equal(hee_apply_recipe(recipe, data.frame(bio1 = 3))$bio1, 1)
})
