test_that("core process catalog uses six explicit process names", {
  catalog <- hee_core_process_catalog()
  expected <- c(
    "environmental_filtering",
    "dispersal",
    "biotic_filtering",
    "evolution",
    "speciation",
    "extinction"
  )

  core <- catalog[catalog$layer == "core_process", ]
  expect_equal(core$process_slug, expected)
  expect_true(all(c("process_name_zh", "process_name_en",
                    "function_prefix") %in% names(core)))
  expect_false(any(grepl("^P[0-9]+$", core$process_slug)))
  expect_true(any(grepl("BioGeoBEARS", catalog$interpretation_boundary,
                        fixed = TRUE)))
})

test_that("function catalog maps canonical names to layers and legacy names", {
  catalog <- hee_function_catalog()
  required <- c(
    "hee_environmental_filtering_suitability",
    "hee_dispersal_landscape_weights",
    "hee_dispersal_spherical_kernel",
    "hee_dispersal_transition_matrix",
    "hee_dispersal_spacetime_graph",
    "hee_dispersal_path_diagnostic",
    "hee_dispersal_pruning_likelihood",
    "hee_dispersal_propagule_particles",
    "hee_dispersal_kernel",
    "hee_biotic_filtering_effect",
    "hee_evolution_ancestral_response",
    "hee_speciation_events",
    "hee_extinction_probability",
    "hee_global_grid_dynamic_occupancy"
  )

  expect_true(all(required %in% catalog$canonical_name))
  legacy <- catalog[catalog$canonical_name %in% c(
    "hee_colonisation_arrival_pressure",
    "hee_colonisation_establishment_probability",
    "hee_colonisation_probability_from_arrival",
    "hee_colonisation_entry_seed",
    "hee_persistence_probability"
  ), ]
  expect_equal(nrow(legacy), 5L)
  expect_true(all(legacy$layer == "legacy_compatibility"))
  expect_true(all(legacy$process_slug == "legacy_occupancy_transition"))
  expect_equal(
    catalog$legacy_name[catalog$canonical_name ==
                          "hee_environmental_filtering_suitability"],
    "hee_lineage_suitability"
  )
  expect_equal(
    catalog$process_slug[catalog$canonical_name == "hee_bgb_region_history"],
    "historical_constraint"
  )
  expect_true(any(catalog$status == "legacy_compatibility"))
  expect_equal(
    catalog$status[catalog$canonical_name == "hee_dispersal_kernel"],
    "legacy_compatibility"
  )
  expect_false(any(grepl("^P[0-9]+$", catalog$process_slug)))
})

test_that("six-process aliases and legacy occupancy helpers preserve calculations", {
  earth <- data.frame(cell_id = "c1", time_ma = 0, region = "R1",
                      land = 1, bio1 = 1)
  E <- hee_paleo_earth_state(earth, env_cols = "bio1")
  response <- data.frame(lineage = "sp1", intercept = 0, bio1 = 1)

  expect_equal(
    hee_environmental_filtering_suitability(E, response, basis_cols = "bio1"),
    hee_lineage_suitability(E, response, basis_cols = "bio1")
  )

  paths <- data.frame(lineage = "sp1", from_cell_id = "c1", to_cell_id = "c1",
                      from_region = "R1", to_region = "R1", time_ma = 0,
                      least_cost_distance_km = 0)
  expect_equal(
    hee_dispersal_kernel(paths, movement_intercept = 0,
                         distance_decay = 0, delta_t = 1),
    hee_within_region_movement_kernel(paths, movement_intercept = 0,
                                      distance_decay = 0, delta_t = 1)
  )

  previous <- data.frame(lineage = "sp1", cell_id = "c1", q = 0.5)
  kernel <- hee_dispersal_kernel(paths, movement_intercept = 0,
                                 distance_decay = 0, delta_t = 1)
  expect_equal(
    hee_colonisation_arrival_pressure(previous, kernel),
    hee_arrival_pressure(previous, kernel)
  )

  x <- data.frame(eta = 0, arrival_pressure = 0.4)
  est <- hee_colonisation_establishment_probability(x)
  expect_equal(est$establishment_probability,
               hee_establishment_probability(x)$establishment_probability)
  expect_equal(
    hee_colonisation_probability_from_arrival(est),
    hee_colonisation_from_arrival(est)
  )
})

test_that("biotic filtering is explicit and neutral without external evidence", {
  x <- data.frame(lineage = "sp1", cell_id = "c1", time_ma = 0)
  out <- hee_biotic_filtering_effect(x)

  expect_equal(out$biotic_filtering_effect, 0)
  expect_match(out$biotic_filtering_source,
               "not_parameterised", fixed = TRUE)

  y <- hee_biotic_filtering_effect(data.frame(raw_interaction = c(-2, NA, 3)),
                                   biotic_col = "raw_interaction")
  expect_equal(y$biotic_filtering_effect, c(-2, 0, 3))
  expect_match(y$biotic_filtering_source[1], "external", fixed = TRUE)
})
