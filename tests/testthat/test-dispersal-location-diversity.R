test_that("terminal range-support calibration recovers the declared support fraction", {
  density <- rbind(
    a = c(sp1 = 0.7, sp2 = 0.0),
    b = c(sp1 = 0.3, sp2 = 0.4),
    c = c(sp1 = 0.0, sp2 = 0.6)
  )
  calibration <- hee_dispersal_range_support_calibrate(
    density, c(sp1 = 2, sp2 = 2), target_fraction = 0.9
  )
  expect_equal(calibration$achieved_terminal_support_cells,
               calibration$target_terminal_support_cells, tolerance = 1e-6)
  expect_true(all(calibration$support_scale > 0))
})

test_that("location-support diversity is bounded and never claims occupancy", {
  density <- rbind(
    a = c(sp1 = 0.5, sp2 = 0.0),
    b = c(sp1 = 0.5, sp2 = 1.0)
  )
  result <- hee_dispersal_range_support_diversity(
    density, c(sp1 = 1, sp2 = 1), return_support = TRUE
  )
  expect_equal(nrow(result$summary), 2)
  expect_true(all(result$summary$support_mixture_shannon >= 0))
  expect_true(all(result$summary$support_mixture_effective_lineages >= 1))
  expect_true(all(result$support >= 0 & result$support <= 1))
  expect_match(result$interpretation, "not posterior occupancy", fixed = TRUE)
})

test_that("tiny positive support still gives a normalized lineage mixture", {
  density <- rbind(
    tiny = c(sp1 = 8e-17, sp2 = 8e-17),
    main = c(sp1 = 1 - 8e-17, sp2 = 1 - 8e-17)
  )
  result <- hee_dispersal_range_support_diversity(
    density, c(sp1 = 1, sp2 = 1)
  )$summary
  expect_equal(result$support_mixture_effective_lineages, c(2, 2),
               tolerance = 1e-12)
  expect_lte(max(result$support_mixture_effective_lineages), 2)
})
