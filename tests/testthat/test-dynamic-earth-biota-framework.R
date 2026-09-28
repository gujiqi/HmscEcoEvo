test_that("HmscEE catalog no longer exposes obsolete multiplicative core", {
  p <- hee_hmscee_process_catalog()
  expect_equal(
    p$process_slug[p$layer == "core_process"],
    c("environmental_filtering", "dispersal",
      "biotic_filtering", "evolution", "speciation",
      "extinction")
  )
  expect_false(any(p$process_name_en[p$layer == "core_process"] %in%
                     c("BioGeoBEARS", "Fossils", "Extrapolation")))
  expect_true(any(grepl("BioGeoBEARS", p$interpretation_boundary,
                        fixed = TRUE)))

  f <- hee_hmscee_formula_catalog()
  expect_true(all(c("regional_history", "movement_kernel",
                    "speciation_identity", "final_probability") %in%
                    f$formula_id))
  expect_false(any(c("colonisation", "persistence", "establishment") %in%
                   f$formula_id))
  expect_true(all(c("framework_layer", "core_process") %in% names(f)))
  expect_equal(
    f$core_process[match(c("ancestral_environmental_response",
                           "regional_history", "movement_kernel",
                           "extrapolation_uncertainty"), f$formula_id)],
    c("evolution", "dispersal/speciation", "dispersal", "not_core_process")
  )
  expect_false(any(grepl("A_hist.*D_dynamic.*Q_extrap", f$formula)))
  expect_match(
    f$interpretation[f$formula_id == "final_probability"],
    "conditional lineage-location mass",
    fixed = TRUE
  )
})

test_that("old product-core entry points fail loudly instead of multiplying", {
  expect_error(
    hee_dynamic_earth_biota_probability(),
    "obsolete multiplicative HmscEE core",
    fixed = TRUE
  )
  expect_error(
    hee_static_snapshot_approximation(data.frame()),
    "obsolete multiplicative HmscEE core",
    fixed = TRUE
  )
  expect_error(
    hee_regional_species_pool_transition(data.frame()),
    "obsolete multiplicative HmscEE core",
    fixed = TRUE
  )
})

test_that("palaeo-Earth state treats land and habitat as hard states", {
  earth <- data.frame(
    cell_id = c("c1", "c2", "c3"),
    time_ma = 0,
    region = "R1",
    land = c(1, 1, NA),
    habitat = c(1, 0, 1),
    bio1 = c(0, 0, 0)
  )
  out <- hee_paleo_earth_state(earth, habitat_col = "habitat",
                               env_cols = "bio1")
  expect_equal(out$H_state, c(1L, 0L, 0L))
  expect_warning(
    alias <- hee_geographic_stage(earth, habitat_col = "habitat"),
    "hard-state `H_state`",
    fixed = TRUE
  )
  expect_equal(alias$H_state, out$H_state)
})

test_that("BioGeoBEARS region history is binary and uniquely keyed", {
  region_history <- data.frame(
    lineage = "sp1", region = "R1", time_ma = c(10, 0),
    history_draw = "h1", accessibility = c(0.3, 0)
  )
  out <- hee_bsm_region_history(region_history)
  expect_equal(out$R_region, c(1L, 0L))
  dup <- rbind(region_history, region_history[1, ])
  expect_error(hee_bsm_region_history(dup), "duplicate", ignore.case = TRUE)
})

test_that("BioGeoBEARS histories are conditional on geography scenario, not BGB e", {
  earth <- expand.grid(
    geography_scenario = c("g_landbridge", "g_no_bridge"),
    climate_scenario = "c_model1",
    cell_id = "c1",
    time_ma = 0,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  earth$region <- "R1"
  earth$land <- 1
  earth$bio1 <- 0
  E <- hee_paleo_earth_state(earth, env_cols = "bio1")
  S <- hee_lineage_suitability(
    E,
    data.frame(lineage = "sp1", intercept = 1, bio1 = 0),
    basis_cols = "bio1"
  )
  region_history <- data.frame(
    geography_scenario = c("g_landbridge", "g_landbridge", "g_no_bridge"),
    history_draw = c("h1", "h2", "h1"),
    lineage = "sp1",
    region = "R1",
    time_ma = 0,
    R_region = c(1, 0, 0),
    geography_weight = c(0.5, 0.5, 0.5),
    history_weight_given_geography = c(0.8, 0.2, 1)
  )
  fit <- hee_nested_region_cell_occupancy(
    S, region_history, E, initial_rho = 0.5, min_region_q = 0
  )
  expect_equal(fit$summary$probability_mean, 0.2, tolerance = 1e-8)
  expect_true(all(c("geography_scenario", "climate_scenario",
                    "geography_weight", "history_weight_given_geography",
                    "climate_weight_given_geography") %in% names(fit$draws)))
  expect_equal(unique(fit$draws$geography_weight[
    fit$draws$geography_scenario == "g_landbridge"
  ]),
               0.5)
  expect_equal(sum(unique(fit$draws$history_weight_given_geography[
    fit$draws$geography_scenario == "g_landbridge"
  ])), 1)
})

test_that("palaeoclimate scenarios are weighted conditional on geography", {
  earth <- expand.grid(
    geography_scenario = "g1",
    climate_scenario = c("c_cool", "c_warm"),
    cell_id = c("c1", "c2"),
    time_ma = 0,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  earth$region <- "R1"
  earth$land <- 1
  earth$bio1 <- ifelse(earth$cell_id == "c1", 2, 0)
  earth$climate_weight_given_geography <- ifelse(
    earth$climate_scenario == "c_cool", 0.25, 0.75
  )
  E <- hee_paleo_earth_state(
    earth, env_cols = "bio1",
    geography_scenario_col = "geography_scenario",
    climate_scenario_col = "climate_scenario",
    climate_weight_col = "climate_weight_given_geography"
  )
  S <- hee_lineage_suitability(
    E,
    data.frame(lineage = "sp1", response_draw = "s1",
               response_weight = 1, intercept = 0, bio1 = 1),
    basis_cols = "bio1"
  )
  R <- data.frame(geography_scenario = "g1", history_draw = "h1",
                  lineage = "sp1", region = "R1", time_ma = 0,
                  R_region = 1, geography_weight = 1,
                  history_weight_given_geography = 1)
  fit <- hee_nested_region_cell_occupancy(S, R, E, initial_rho = 1,
                                          min_region_q = 0)
  c1 <- fit$summary$probability_mean[fit$summary$cell_id == "c1"]
  c2 <- fit$summary$probability_mean[fit$summary$cell_id == "c2"]
  expect_gt(c1, c2)
  expect_equal(unique(fit$draws$climate_weight_given_geography[
    fit$draws$climate_scenario == "c_cool"
  ]), 0.25)
  expect_equal(unique(fit$summary$n_climate_scenarios), 2)
})

test_that("trait-mediated ancestral response separates traits, beta and geography", {
  ancestor_traits <- data.frame(
    lineage = "anc1", time_ma = 10, response_draw = "d1",
    leaf = 2, wood = 1, region = "not_a_trait"
  )
  Gamma <- matrix(
    c(0.5, -0.1, 0.2, 0.3),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("leaf", "wood"), c("bio1", "bio12"))
  )
  u <- data.frame(lineage = "anc1", time_ma = 10, response_draw = "d1",
                  bio1 = 0.3, bio12 = 0)
  beta <- hee_trait_mediated_ancestral_response(
    ancestor_traits, Gamma, residual_response = u
  )
  expect_equal(beta$bio1, 1.5)
  expect_equal(beta$bio12, 0.1)
  expect_match(beta$response_source, "Gamma_T_plus_u")
  expect_error(
    hee_trait_mediated_ancestral_response(
      ancestor_traits, Gamma, trait_cols = c("leaf", "region")
    ),
    "BioGeoBEARS range columns cannot be used as ancestral niche traits",
    fixed = TRUE
  )

  earth <- data.frame(
    cell_id = c("c1", "c1"), time_ma = c(10, 0), region = "R1",
    land = 1, bio1 = c(1, 1), bio12 = c(0, 0)
  )
  E <- hee_paleo_earth_state(earth, env_cols = c("bio1", "bio12"))
  S <- hee_lineage_suitability(E, beta, basis_cols = c("bio1", "bio12"))
  expect_equal(unique(S$time_ma), 10)
})

test_that("within-region movement kernel blocks ordinary cross-region dispersal", {
  paths <- data.frame(
    lineage = "sp1",
    from_cell_id = c("c1", "c1"),
    to_cell_id = c("c2", "c3"),
    from_region = "R1",
    to_region = c("R1", "R2"),
    time_ma = 0,
    least_cost_distance_km = c(10, 10)
  )
  out <- hee_within_region_movement_kernel(paths, delta_t = 1,
                                           movement_intercept = 0,
                                           distance_decay = 0)
  expect_gt(out$K_movement[out$to_region == "R1"], 0)
  expect_equal(out$K_movement[out$to_region == "R2"], 0)
})

test_that("arrival and colonisation are arrival times establishment only", {
  previous <- data.frame(lineage = "sp1", cell_id = "c1", q = 0.5)
  kernel <- data.frame(
    lineage = "sp1", from_cell_id = "c1", to_cell_id = "c2",
    time_ma = 0, K_movement = 0.4
  )
  target <- data.frame(lineage = "sp1", cell_id = "c2", time_ma = 0,
                       eta = 0)
  arr <- hee_arrival_pressure(previous, kernel, target)
  expect_equal(arr$arrival_pressure, 0.2)
  est <- hee_establishment_probability(arr, intercept = 0, slope = 0)
  gamma <- hee_colonisation_from_arrival(est)
  expect_equal(gamma$colonisation_probability, 0.1)
})

test_that("nested occupancy honours BGB region state and hard habitat state", {
  earth <- expand.grid(
    cell_id = c("c1", "c2", "c3"),
    time_ma = c(10, 0),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  earth$region <- c("R1", "R1", "R2", "R1", "R1", "R2")
  earth$land <- c(1, 0, 1, 1, 1, 1)
  earth$bio1 <- c(0, 0, 0, 0, 0, 0)
  E <- hee_paleo_earth_state(earth, env_cols = "bio1")
  responses <- data.frame(lineage = "sp1", intercept = 1, bio1 = 0)
  S <- hee_lineage_suitability(E, responses, basis_cols = "bio1")
  region_history <- data.frame(
    lineage = "sp1",
    region = rep(c("R1", "R2"), each = 2),
    time_ma = c(10, 0, 10, 0),
    R_region = c(1, 1, 0, 1)
  )
  paths <- data.frame(
    lineage = "sp1",
    from_cell_id = c("c1", "c1", "c2", "c2", "c3"),
    to_cell_id = c("c1", "c2", "c1", "c2", "c3"),
    from_region = c("R1", "R1", "R1", "R1", "R2"),
    to_region = c("R1", "R1", "R1", "R1", "R2"),
    time_ma = 0,
    least_cost_distance_km = c(0, 1, 1, 0, 0)
  )
  K <- hee_within_region_movement_kernel(paths, delta_t = 10,
                                         movement_intercept = 0,
                                         distance_decay = 0)
  fit <- hee_nested_region_cell_occupancy(S, region_history, E,
                                          movement_kernel = K,
                                          initial_rho = 0.5)
  draws <- fit$draws
  expect_true(all(draws$probability_draw >= 0 & draws$probability_draw <= 1))
  expect_equal(draws$probability_draw[draws$cell_id == "c2" &
                                        draws$time_ma == 10], 0)
  expect_equal(draws$probability_draw[draws$region == "R2" &
                                        draws$time_ma == 10], 0)
  entered <- draws[draws$region == "R2" & draws$time_ma == 0, ]
  expect_true(any(entered$transition_case == "BioGeoBEARS_region_entry"))
  expect_true(any(entered$probability_draw > 0))
})

test_that("region pool is derived from cell occupancy", {
  x <- data.frame(
    lineage = "sp1", region = "R1", time_ma = 0,
    cell_id = c("c1", "c2"), probability_mean = c(0.2, 0.5)
  )
  out <- hee_region_pool_from_cell_occupancy(x)
  expect_equal(out$regional_probability, 1 - (1 - 0.2) * (1 - 0.5))
})

test_that("extrapolation report does not downweight biological probability", {
  p <- data.frame(cell_id = "c1", time_ma = 0, probability_mean = 0.8)
  q <- data.frame(cell_id = "c1", time_ma = 0, extrapolation_score = 1)
  out <- hee_extrapolation_uncertainty_report(p, q)
  expect_equal(out$probability_mean, 0.8)
  expect_equal(out$extrapolation_score, 1)
  expect_match(out$probability_interpretation,
               "not_downweighted",
               fixed = TRUE)
})
