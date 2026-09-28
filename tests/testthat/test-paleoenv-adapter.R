test_that("paleoenvironment role map prefers DEM land and primary climate layers", {
  vars <- c(
    "LANDFRAC_li_fraction", "land_mask_dem", "land_area_km2",
    "PHIS_li_m", "elevation_m", "MAT_gridded_model_mean_C",
    "MAT_pohl_C", "MAP_gridded_model_mean_mm_yr", "MAP_pohl_mm_yr",
    "moisture_availability_index_z", "bryophyte_moisture_score_z"
  )
  mp <- hee_paleoenv_default_variable_map(variables = vars)
  expect_equal(mp$variable[mp$role == "land_mask"], "land_mask_dem")
  expect_equal(mp$variable[mp$role == "elevation"], "elevation_m")
  expect_equal(mp$variable[mp$role == "temperature"], "MAT_pohl_C")
  expect_equal(mp$variable[mp$role == "precipitation"], "MAP_pohl_mm_yr")
  expect_equal(mp$variable[mp$role == "moisture"], "moisture_availability_index_z")
})

test_that("paleoenvironment adapter creates arena and opportunity inputs from a cube", {
  cube <- hee_load_timecube(object = make_test_timecube_540(times = c(10, 0)))
  out <- hee_prepare_paleoenv_geoprocess_inputs(
    cube,
    times = c(10, 0),
    hmsc_variables = c("bio1", "bio12"),
    include_opportunity = TRUE
  )
  expect_s3_class(out, "hee_paleoenv_geoprocess_inputs")
  expect_true(all(c("paleo_grid", "arena", "landscape_state",
                    "land_age", "landscape_events",
                    "ecological_opportunity") %in% names(out)))
  expect_equal(nrow(out$arena), nrow(out$paleo_grid))
  expect_true(all(out$arena$L_land >= 0 & out$arena$L_land <= 1))
  expect_true(all(out$arena$G_arena >= 0 & out$arena$G_arena <= 1))
  expect_true(any(out$arena$L_land == 0))
  expect_true(all(out$ecological_opportunity$ecological_opportunity >= 0 &
                    out$ecological_opportunity$ecological_opportunity <= 1))
})

test_that("habitat strategies are explicit and do not silently invent habitat", {
  cube <- hee_load_timecube(object = make_test_timecube_540(times = c(0)))
  expect_error(
    hee_prepare_paleoenv_geoprocess_inputs(
      cube,
      times = 0,
      habitat_strategy = "bryophyte"
    ),
    "matching palaeoenvironment layer is not available"
  )
  out <- hee_prepare_paleoenv_geoprocess_inputs(
    cube,
    times = 0,
    habitat_strategy = "land_only",
    include_opportunity = FALSE
  )
  expect_equal(out$arena$G_arena, out$arena$L_land)
})
