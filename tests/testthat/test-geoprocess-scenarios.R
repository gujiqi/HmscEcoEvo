expect_geoprob01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  expect_false(anyNA(x))
  expect_true(all(is.finite(x)))
  expect_true(all(x >= 0 & x <= 1))
}

test_that("geoprocess scenario suite covers normal and stress scenarios", {
  suite <- hee_geoprocess_scenario_suite(species = c("sp1", "sp2", "sp3"),
                                         seed = 20260727)
  expect_s3_class(suite, "hee_geoprocess_scenario_suite")
  expect_equal(
    suite$scenario_summary$status[
      match("duplicate_key", suite$scenario_summary$scenario)
    ],
    "EXPECTED_ERROR"
  )
  ok <- suite$scenario_summary$scenario != "duplicate_key"
  expect_true(all(suite$scenario_summary$status[ok] == "PASS"))
  expect_true(all(suite$scenario_summary$n_validation_failures == 0))

  required_inputs <- c("landscape_state", "projection", "traits",
                       "phylo_mask", "species_origin")
  required_results <- c("landscape_events", "ecological_opportunity",
                        "region_state", "source_pressure",
                        "rescue_effect", "colonisation_probability",
                        "extinction_probability", "speciation_opportunity",
                        "dynamic_occupancy", "species_pool",
                        "mechanism_attribution", "reliability",
                        "richness", "turnover", "refugia",
                        "network_node_metrics", "network_summary",
                        "corridor_persistence", "local_extinction",
                        "regional_extinction", "lineage_extinction",
                        "cradle_museum_grave", "dynamic_model_family",
                        "mechanism_hypotheses")
  expect_true(all(required_inputs %in% names(suite$input_tables)))
  expect_true(all(required_results %in% names(suite$result_tables)))
  expect_gt(nrow(suite$formula_catalog), 0)
  expect_gt(nrow(suite$interpretation_table), 0)
  expect_true(all(paste0("M", 0:8) %in%
                    suite$result_tables$dynamic_model_family$model_id))
  expect_true(all(paste0("H", 1:8) %in%
                    suite$result_tables$mechanism_hypotheses$hypothesis_id))
  expect_true(all(c(
    "core_environmental_filtering", "core_dispersal",
    "core_biotic_filtering", "core_evolution",
    "core_speciation", "core_extinction"
  ) %in% suite$formula_catalog$process))
  expect_false(any(c("core_colonisation", "core_persistence") %in%
                   suite$formula_catalog$process))
  expect_true("core_process_mapping" %in% names(suite$formula_catalog))
  expect_true(any(grepl("not a biological process",
                        suite$formula_catalog$core_process_mapping,
                        fixed = TRUE)))
})

test_that("geoprocess scenario suite preserves scientific boundary conditions", {
  suite <- hee_geoprocess_scenario_suite(species = c("sp1", "sp2"),
                                         seed = 11)

  dyn <- suite$result_tables$dynamic_occupancy
  for (nm in c("colonisation_probability", "extinction_probability",
               "occupancy_probability", "geographic_existence")) {
    expect_geoprob01(dyn[[nm]])
  }
  absent <- dyn$scenario == "land_appearance_loss" &
    dyn$geographic_existence <= 0
  expect_true(any(absent))
  expect_true(all(dyn$occupancy_probability[absent] == 0))

  src_missing <- suite$result_tables$source_pressure[
    suite$result_tables$source_pressure$scenario == "missing_optional",
    "source_pressure"
  ]
  expect_true(length(src_missing) > 0)
  expect_true(all(src_missing == 0))

  events <- suite$result_tables$landscape_events
  loss_events <- events$event_type[events$scenario == "land_appearance_loss"]
  expect_true(all(c("emergence", "submergence") %in% loss_events))

  intervals <- suite$result_tables$region_state$interval_myr[
    suite$result_tables$region_state$scenario == "uneven_time"
  ]
  expect_gt(length(unique(intervals)), 1)

  net <- suite$result_tables$network_summary
  expect_true(all(c("fragmentation_index", "largest_component_fraction",
                    "network_components") %in% names(net)))
  expect_geoprob01(net$fragmentation_index)
  expect_geoprob01(net$largest_component_fraction)

  lin <- suite$result_tables$lineage_extinction
  expect_true(all(c("lineage_extinction_proxy",
                    "lineage_persistence_proxy") %in% names(lin)))
  expect_geoprob01(lin$lineage_extinction_proxy)
  expect_geoprob01(lin$lineage_persistence_proxy)

  cmg <- suite$result_tables$cradle_museum_grave
  expect_true(all(c("cradle_score", "museum_score", "grave_score",
                    "dominant_region_status") %in% names(cmg)))
  expect_true(all(cmg$dominant_region_status %in%
                    c("cradle", "museum", "grave", "source", "sink",
                      "undetermined")))
})

test_that("geoprocess scenario suite records duplicate-key failures without repairing them", {
  suite <- hee_geoprocess_scenario_suite(scenarios = "duplicate_key",
                                         species = c("sp1", "sp2"))
  expect_equal(suite$scenario_summary$status, "EXPECTED_ERROR")
  expect_match(suite$scenario_summary$error_message, "Duplicate key")
  expect_equal(suite$validation_table$status, "PASS")
  expect_match(suite$validation_table$check, "expected_clear_error")
})

test_that("ecological opportunity tolerates env column names that collide with internal summaries", {
  land <- expand.grid(cell_id = c("c1", "c2"),
                      time_ma = c(10, 0),
                      stringsAsFactors = FALSE)
  land$region <- "R1"
  land$cell_area_km2 <- c(1, 2, 3, 4)
  land$habitat_heterogeneity <- c(0.2, 0.8, 0.4, 0.7)
  out <- hee_ecological_opportunity(
    land,
    env_cols = "habitat_heterogeneity",
    area_col = "cell_area_km2"
  )
  expect_equal(nrow(out), nrow(land))
  expect_true("heterogeneity_index" %in% names(out))
  expect_false(".hee_habitat_heterogeneity" %in% names(out))
  expect_geoprob01(out$ecological_opportunity)
})

test_that("geoprocess scenario plots return printable ggplot objects", {
  skip_if_not_installed("ggplot2")
  suite <- hee_geoprocess_scenario_suite(
    scenarios = c("normal", "land_appearance_loss"),
    species = c("sp1", "sp2")
  )
  plots <- plot_geoprocess_scenario_suite(suite)
  expect_s3_class(plots, "hee_plot_list")
  expect_true(all(c("process_heatmap", "validation_status",
                    "landscape_events", "dynamic_occupancy",
                    "network_fragmentation", "extinction_layers",
                    "cradle_museum_grave") %in%
                    names(plots)))
  expect_true(all(vapply(plots, inherits, logical(1), what = "ggplot")))
})

test_that("single-time geoprocess scenario plots points without line warnings", {
  skip_if_not_installed("ggplot2")
  suite <- hee_geoprocess_scenario_suite(
    scenarios = "single_time",
    species = c("sp1", "sp2")
  )
  plots <- plot_geoprocess_scenario_suite(suite)
  expect_s3_class(plots$dynamic_occupancy, "ggplot")
  expect_silent(print(plots$dynamic_occupancy))
  expect_silent(print(plots$network_fragmentation))
  expect_silent(print(plots$extinction_layers))
})
