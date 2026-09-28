test_that("extrapolation risk supports full env-cube interface", {
  cube <- hee_load_timecube(object = make_test_timecube_540(times = c(540, 0)))
  env_now <- data.frame(bio1 = c(10, 20), bio12 = c(850, 950), elev = c(0, 100))
  risk <- hee_extrapolation_risk(env_now, cube,
                                 variables = c("bio1", "bio12", "elev"),
                                 times = c(540, 0),
                                 method = "range")
  expect_true(all(c("extrapolation_score", "extrapolation_flag",
                    "no_analog_fraction", "variable_out_of_range") %in% names(risk)))
  expect_true(nrow(risk) > 0)
  expect_true(any(grepl("^out_of_range_", names(risk))))
})

test_that("extrapolation score applies threshold once when producing weights", {
  cube <- hee_load_timecube(object = make_test_timecube_540(times = c(540, 0)))
  env_now <- data.frame(bio1 = c(14.5, 15.5))
  risk <- hee_extrapolation_risk(env_now, cube, variables = "bio1",
                                 times = 540, threshold = 0.5,
                                 method = "range")
  expect_true(any(risk$extrapolation_flag))
  expect_true(any(risk$extrapolation_weight < 1))
})

test_that("extrapolation risk gives clear errors for missing or unusable variables", {
  cube <- hee_load_timecube(object = make_test_timecube_540(times = c(540, 0)))
  expect_error(
    hee_extrapolation_risk(data.frame(bio1 = 1:2), cube,
                           variables = c("bio1", "bio12"), times = 0),
    "Modern reference data are missing"
  )
  expect_error(
    hee_extrapolation_risk(data.frame(bio1 = c(NA_real_, NA_real_)), cube,
                           variables = "bio1", times = 0),
    "no finite values"
  )
})

test_that("extrapolation risk treats missing paleo rows as high-risk rather than reliable", {
  env_now <- data.frame(bio1 = c(10, 20), bio12 = c(800, 900))
  past <- data.frame(cell_id = c("c1", "c2"), time_ma = c(100, 100),
                     bio1 = c(15, NA_real_), bio12 = c(850, NA_real_))
  recipe <- hee_lock_recipe(env_now, formula = ~ bio1 + bio12)
  risk_recipe <- hee_extrapolation_risk(recipe, past, threshold = 2)
  expect_false(risk_recipe$extrapolation_flag[1])
  expect_true(risk_recipe$extrapolation_flag[2])
  expect_lt(risk_recipe$extrapolation_weight[2], risk_recipe$extrapolation_weight[1])

  partial <- past
  partial$bio12[2] <- 850
  risk_partial <- hee_extrapolation_risk(recipe, partial, threshold = 2)
  expect_true(risk_partial$extrapolation_flag[2])
  expect_lt(risk_partial$extrapolation_weight[2], risk_partial$extrapolation_weight[1])
})
