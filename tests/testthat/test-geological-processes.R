test_that("geoprocess model family is explicit", {
  models <- hee_define_geoprocess_models()
  expect_equal(models$model_id, paste0("M", 0:8))
  expect_true(any(grepl("dynamic", models$model_name)))
  expect_true(any(grepl("speciation", models$mechanisms)))
})

test_that("landscape events detect arena creation and transformation", {
  landscape <- data.frame(
    region = "A",
    time_ma = c(10, 5, 0),
    land = c(0, 1, 1),
    habitat = c(0, 1, 0),
    area_km2 = c(0, 100, 40),
    elev_m = c(0, 30, 250)
  )
  ev <- hee_landscape_events(
    landscape,
    habitat_col = "habitat",
    area_col = "area_km2",
    elevation_col = "elev_m",
    elevation_threshold = 100
  )
  expect_true(all(c("emergence", "habitat_gain", "area_contraction", "uplift") %in%
                    ev$event_type))
  expect_true(all(ev$time_start_ma >= ev$time_end_ma))
})

test_that("land age and ecological opportunity follow arena dynamics", {
  landscape <- data.frame(
    cell_id = "c1",
    region = "A",
    time_ma = c(10, 5, 0),
    land = c(0, 1, 1),
    cell_area_km2 = c(0, 10, 20),
    bio1 = c(0, 2, 4)
  )
  aged <- hee_land_age(landscape)
  expect_true(is.na(aged$land_age_myr[aged$time_ma == 10]))
  expect_equal(aged$land_age_myr[aged$time_ma == 5], 0)
  expect_equal(aged$land_age_myr[aged$time_ma == 0], 5)
  opp <- hee_ecological_opportunity(aged, env_cols = "bio1")
  expect_true(all(opp$ecological_opportunity >= 0 & opp$ecological_opportunity <= 1))
  expect_true(any(opp$young_land_index > 0, na.rm = TRUE))
})

test_that("connectivity cube and isolation history are bounded and species-specific", {
  distances <- data.frame(
    from_region = rep("A", 3),
    to_region = rep("B", 3),
    time_ma = c(10, 5, 0),
    paleodistance = c(1, 5, 5),
    barrier_strength = c(0, 0.7, 0.7)
  )
  traits <- data.frame(species = c("good", "poor"),
                       dispersal_distance = c(10, 1))
  cc <- hee_build_connectivity_cube(distances, traits = traits,
                                    distance_scale = 5)
  expect_true(all(cc$connectivity >= 0 & cc$connectivity <= 1))
  expect_gt(cc$connectivity[cc$species == "good" & cc$time_ma == 10],
            cc$connectivity[cc$species == "poor" & cc$time_ma == 10])

  iso <- hee_isolation_history(cc, threshold = 0.2)
  poor <- iso[iso$species == "poor", ]
  poor <- poor[order(poor$time_ma, decreasing = TRUE), ]
  expect_equal(poor$isolation_duration_ma[1], 0)
  expect_gte(tail(poor$isolation_duration_ma, 1), 5)
})

test_that("structural, functional, and climatic connectivity are separated", {
  pairs <- data.frame(from_region = "A", to_region = "B", time_ma = 0,
                      paleodistance = 10, barrier_strength = 0.5,
                      environment_distance = 2)
  s <- hee_structural_connectivity(pairs, alpha = 0.01)
  expect_true("structural_connectivity" %in% names(s))
  f <- hee_functional_connectivity(s, traits = data.frame(species = c("a", "b"),
                                                          dispersal = c(0, 10)),
                                   expand_species = TRUE)
  expect_true(all(f$functional_connectivity >= 0 & f$functional_connectivity <= 1))
  c <- hee_climatic_connectivity(f, phi = 0.5)
  expect_true(all(c$climate_connectivity >= 0 & c$climate_connectivity <= 1))
  expect_true(all(c$connectivity >= 0 & c$connectivity <= 1))
})

test_that("explicit missing climate corridor values are conservative, while omitted climate is neutral", {
  explicit <- hee_climatic_connectivity(
    data.frame(from_region = "A", to_region = "B", time_ma = 0,
               climate_connectivity = NA_real_)
  )
  expect_equal(explicit$climate_connectivity, 0)
  expect_equal(explicit$connectivity, 0)

  corridor <- hee_climatic_connectivity(
    data.frame(from_region = "A", to_region = "B", time_ma = 0,
               corridor_suitability = NA_real_)
  )
  expect_equal(corridor$climate_connectivity, 0)
  expect_equal(corridor$connectivity, 0)

  omitted <- hee_climatic_connectivity(
    data.frame(from_region = "A", to_region = "B", time_ma = 0)
  )
  expect_equal(omitted$climate_connectivity, 1)
  expect_equal(omitted$connectivity, 1)

  cube <- hee_build_connectivity_cube(
    data.frame(from_region = "A", to_region = "B", time_ma = 0,
               paleodistance = 1),
    climate = data.frame(from_region = "A", to_region = "B", time_ma = 0,
                         climate_connectivity = NA_real_)
  )
  expect_equal(cube$climate_connectivity, 0)
  expect_equal(cube$connectivity, 0)

  barrier_cube <- hee_build_connectivity_cube(
    data.frame(from_region = "A", to_region = c("B", "C"), time_ma = 0,
               paleodistance = 1),
    barrier = data.frame(from_region = "A", to_region = "B", time_ma = 0,
                         barrier_passability = NA_real_)
  )
  expect_equal(barrier_cube$barrier_passability, c(0, 0))
  expect_equal(barrier_cube$structural_connectivity, c(0, 0))
})

test_that("barrier strength without passability is a partial barrier, not a closed route", {
  cube <- hee_build_connectivity_cube(
    data.frame(from_region = "A", to_region = "B", time_ma = 0,
               paleodistance = 0),
    barrier = data.frame(from_region = "A", to_region = "B", time_ma = 0,
                         barrier_strength = 0.2)
  )
  expect_equal(cube$barrier_passability, 0.8)
  expect_equal(cube$structural_connectivity, 0.8)
  expect_equal(cube$connectivity, 0.8)
})

test_that("negative palaeodistance and climate distance are rejected, not rescaled to high connectivity", {
  expect_error(
    hee_structural_connectivity(
      data.frame(from_region = "A", to_region = "B", time_ma = 0,
                 paleodistance = -10)
    ),
    "non-negative"
  )
  expect_error(
    hee_functional_connectivity(
      data.frame(from_region = "A", to_region = "B", time_ma = 0,
                 paleodistance = -10),
      traits = data.frame(species = "sp1", dispersal = 5),
      expand_species = TRUE
    ),
    "non-negative"
  )
  expect_error(
    hee_climatic_connectivity(
      data.frame(from_region = "A", to_region = "B", time_ma = 0,
                 environment_distance = -1)
    ),
    "non-negative"
  )
  expect_error(
    hee_build_connectivity_cube(
      data.frame(from_region = "A", to_region = "B", time_ma = 0,
                 paleodistance = -1)
    ),
    "non-negative"
  )
})

test_that("connectivity cube rejects duplicate process keys before joining", {
  distances <- data.frame(from_region = "A", to_region = "B", time_ma = 0,
                          paleodistance = 10)
  dup_dist <- rbind(distances, distances)
  expect_error(hee_build_connectivity_cube(dup_dist), "Duplicate key")

  dup_barrier <- data.frame(from_region = "A", to_region = "B", time_ma = 0,
                            barrier_strength = c(0.1, 0.9))
  expect_error(hee_build_connectivity_cube(distances, barrier = dup_barrier),
               "Duplicate key")

  dup_climate <- data.frame(from_region = "A", to_region = "B", time_ma = 0,
                            climate_connectivity = c(0.1, 0.9))
  expect_error(hee_build_connectivity_cube(distances, climate = dup_climate),
               "Duplicate key")
})

test_that("functional connectivity does not create trait-only rows", {
  conn <- data.frame(
    from_region = c("A", "A"),
    to_region = c("B", "C"),
    time_ma = 0,
    species = c("sp1", "sp2"),
    paleodistance = c(10, 20),
    barrier_strength = c(0.2, 0.3),
    structural_connectivity = c(0.5, 0.4)
  )
  traits <- data.frame(species = c("sp1", "sp2", "sp3"),
                       dispersal = c(1, 10, 5))
  out <- hee_functional_connectivity(conn, traits = traits)
  expect_equal(nrow(out), nrow(conn))
  expect_false(any(is.na(out$from_region)))
  expect_false("sp3" %in% out$species)
})

test_that("functional connectivity rejects ambiguous species expansion and duplicated traits", {
  conn_no_species <- data.frame(from_region = "A", to_region = "B", time_ma = 0,
                                paleodistance = 10)
  traits <- data.frame(species = c("sp1", "sp2"), dispersal = c(1, 2))
  expect_error(
    hee_functional_connectivity(conn_no_species, traits = traits,
                                expand_species = FALSE),
    "expand_species"
  )

  conn_species <- data.frame(from_region = "A", to_region = "B", time_ma = 0,
                             species = "sp1", paleodistance = 10)
  dup_traits <- data.frame(species = c("sp1", "sp1"), dispersal = c(1, 5))
  expect_error(hee_functional_connectivity(conn_species, traits = dup_traits),
               "Duplicate key")
})

test_that("dispersal traits are finite non-negative and zero dispersal is not replaced by a median", {
  conn <- data.frame(from_region = "A", to_region = "B", time_ma = 0,
                     paleodistance = 10)
  expect_error(
    hee_functional_connectivity(conn,
                                traits = data.frame(species = "sp1",
                                                    dispersal = -1),
                                expand_species = TRUE),
    "finite non-negative"
  )
  expect_error(
    hee_build_connectivity_cube(conn,
                                traits = data.frame(species = "sp1",
                                                    dispersal_distance = NA_real_)),
    "finite non-negative"
  )
  expect_error(
    hee_build_connectivity_cube(conn,
                                traits = data.frame(species = "sp1",
                                                    dispersal_distance = -1)),
    "finite non-negative"
  )
  z <- hee_build_connectivity_cube(
    conn,
    traits = data.frame(species = c("zero", "good"),
                        dispersal_distance = c(0, 20)),
    distance_scale = 5
  )
  expect_equal(z$functional_connectivity[z$species == "zero"], 0)
  expect_gt(z$functional_connectivity[z$species == "good"],
            z$functional_connectivity[z$species == "zero"])
})

test_that("colonisation pressure handles matched and unmatched networks", {
  prev <- data.frame(species = "sp1", region = "A", time_ma = 5,
                     occupancy_probability = 0.8)
  cc <- data.frame(from_region = "A", to_region = "B", time_ma = 5,
                   paleodistance = 2, connectivity = 0.5)
  cp <- hee_colonisation_pressure(prev, cc, dispersal_scale = c(sp1 = 5))
  expect_true(all(c("source_pressure", "n_source_regions") %in% names(cp)))
  expect_gt(cp$source_pressure[cp$region == "B"], 0)

  none <- hee_colonisation_pressure(prev, cc[FALSE, ], target = data.frame(
    species = "sp1", region = "C", time_ma = 5
  ))
  expect_equal(none$source_pressure, 0)
})

test_that("source pressure and rescue do not count self-loops by default", {
  prev <- data.frame(species = "sp1", region = "A", time_ma = 0,
                     occupancy_probability = 1)
  cc <- data.frame(from_region = c("A", "A"),
                   to_region = c("A", "B"),
                   time_ma = 0,
                   paleodistance = c(0, 1),
                   connectivity = c(1, 0.5))
  target_a <- data.frame(species = "sp1", region = "A", time_ma = 0)
  target_b <- data.frame(species = "sp1", region = "B", time_ma = 0)

  self_default <- hee_colonisation_pressure(prev, cc, target = target_a)
  self_allowed <- hee_colonisation_pressure(prev, cc, target = target_a,
                                            include_self = TRUE)
  rescue_default <- hee_rescue_effect(prev, cc, target = target_a)
  between_region <- hee_colonisation_pressure(prev, cc, target = target_b)

  expect_equal(self_default$source_pressure, 0)
  expect_equal(rescue_default$rescue_effect, 0)
  expect_equal(self_allowed$source_pressure, 1)
  expect_gt(between_region$source_pressure, 0)
})

test_that("colonisation and extinction probabilities use delta_t and bounds", {
  suit <- data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                     suitability = 0.8)
  no_source <- hee_colonisation_probability(suit, delta_t = 5)
  expect_equal(no_source$source_pressure, 0)
  expect_equal(no_source$colonisation_probability, 0)

  gamma <- hee_colonisation_probability(
    suit,
    source_pressure = data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                                 source_pressure = 0.5),
    delta_t = 5
  )
  expect_true(gamma$colonisation_probability > 0)
  expect_true(gamma$colonisation_probability <= 1)
  expect_lt(gamma$colonisation_probability, 0.5)

  expect_error(
    hee_colonisation_probability(
      suit,
      source_pressure = data.frame(species = c("sp1", "sp1"), cell_id = "c1",
                                   time_ma = 0, source_pressure = c(0.2, 0.8)),
      delta_t = 1
    ),
    "Duplicate key"
  )

  expect_error(
    hee_colonisation_probability(
      suit,
      source_pressure = data.frame(species = "sp1", cell_id = "c1",
                                   time_ma = 0, bad_weight = 0.5),
      delta_t = 1
    ),
    "source_pressure must contain"
  )

  eps <- hee_extinction_probability(data.frame(suitability = 0.1, land_loss = 1),
                                    delta_t = 5)
  expect_equal(eps$extinction_probability, 1)
  expect_true(eps$forced_geographic_extinction)
})

test_that("colonisation and extinction reject invalid intervals and honour geography", {
  suit <- data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                     suitability = 0.8, bad_dt = -5)
  expect_error(hee_colonisation_probability(suit, delta_t = "bad_dt"),
               "delta_t")
  expect_error(hee_extinction_probability(suit, delta_t = -1),
               "delta_t")

  eps <- hee_extinction_probability(
    data.frame(suitability = 0.9, geographic_existence = 0),
    delta_t = 1
  )
  expect_equal(eps$extinction_probability, 1)
  expect_true(eps$forced_geographic_extinction)

  eps_na <- hee_extinction_probability(
    data.frame(suitability = 0.9, geographic_existence = NA_real_),
    delta_t = 1
  )
  risk_na <- hee_extinction_risk(
    data.frame(suitability = 0.9, geographic_existence = NA_real_)
  )
  expect_equal(eps_na$extinction_probability, 1)
  expect_true(eps_na$forced_geographic_extinction)
  expect_equal(risk_na$extinction_probability, 1)
  expect_true(risk_na$forced_geographic_extinction)
})

test_that("source pressure rejects duplicated previous source keys", {
  prev <- data.frame(species = c("sp1", "sp1"), region = c("A", "A"),
                     time_ma = c(5, 5), occupancy_probability = c(0.5, 0.5))
  cc <- data.frame(from_region = "A", to_region = "B", time_ma = 5,
                   connectivity = 0.5)
  expect_error(hee_colonisation_pressure(prev, cc), "Duplicate key")
})

test_that("extinction risk increases with land loss and low suitability", {
  x <- data.frame(
    suitability = c(0.9, 0.1),
    land_loss = c(0, 1),
    isolation = c(0.1, 0.9),
    area_km2 = c(100, 1),
    rescue_effect = c(0.8, 0)
  )
  risk <- hee_extinction_risk(x)
  expect_lt(risk$extinction_probability[1], risk$extinction_probability[2])
  expect_true(all(risk$extinction_probability >= 0 & risk$extinction_probability <= 1))
})

test_that("speciation opportunity indices are bounded and directional", {
  x <- data.frame(
    isolation_duration_ma = c(0, 10),
    niche_divergence = c(0.1, 0.8),
    area_km2 = c(10, 100),
    habitat_heterogeneity = c(0.2, 0.8),
    population_persistence = c(0.2, 0.9),
    source_pressure = c(0.1, 0.5),
    connectivity = c(0.9, 0.2),
    ecological_opportunity = c(0.2, 0.8),
    land_age_ma = c(1, 20)
  )
  sp <- hee_speciation_opportunity(x)
  expect_true(all(sp$speciation_opportunity >= 0 & sp$speciation_opportunity <= 1))
  expect_lt(sp$speciation_opportunity[1], sp$speciation_opportunity[2])
  expect_true(all(is.finite(sp$speciation_rate_lambda)))
})

test_that("individual speciation opportunity helpers are bounded", {
  x <- data.frame(range_fragmentation = 1, isolation_duration_ma = 10,
                  population_persistence = 0.8, rare_dispersal = 0.2,
                  isolation_after_arrival = 0.7, ecological_opportunity = 0.9,
                  habitat_heterogeneity = 0.6, area_km2 = 100,
                  connectivity = 0.2, colonisation_input = 0.5,
                  richness = 3, persistence = 0.9)
  expect_true(hee_allopatric_speciation_opportunity(x)$allopatric_opportunity <= 1)
  expect_true(hee_founder_speciation_opportunity(x)$founder_event_opportunity <= 1)
  expect_true(hee_insitu_speciation_opportunity(x)$in_situ_speciation_index <= 1)
  expect_true(hee_radiation_opportunity(x)$radiation_opportunity <= 1)
  expect_gt(hee_allopatric_speciation_opportunity(x)$allopatric_opportunity, 0)
  expect_gt(hee_insitu_speciation_opportunity(x)$in_situ_speciation_index, 0)
})

test_that("missing speciation process evidence is conservative", {
  allopatric_missing <- hee_allopatric_speciation_opportunity(
    data.frame(isolation_duration_ma = 10)
  )
  insitu_missing_connectivity <- hee_insitu_speciation_opportunity(
    data.frame(ecological_opportunity = 1, habitat_heterogeneity = 1,
               area_km2 = 100)
  )
  combined_missing <- hee_speciation_opportunity(
    data.frame(isolation_duration_ma = 10, niche_divergence = 1,
               area_km2 = 100, habitat_heterogeneity = 1)
  )

  expect_equal(allopatric_missing$allopatric_opportunity, 0)
  expect_equal(insitu_missing_connectivity$in_situ_speciation_index, 0)
  expect_equal(combined_missing$founder_event_opportunity, 0)
  expect_equal(combined_missing$in_situ_speciation_index, 0)
  expect_equal(combined_missing$radiation_opportunity_index, 0)
  expect_equal(combined_missing$speciation_opportunity_lambda,
               combined_missing$speciation_rate_lambda)
})

test_that("species pool updates from regional events", {
  init <- data.frame(region = c("A", "B"), species_pool_size = c(2, 1))
  events <- data.frame(
    time_ma = c(5, 5, 0),
    event_type = c("speciation", "colonisation", "local_extinction"),
    region = c("A", NA, "A"),
    from_region = c(NA, "A", NA),
    to_region = c(NA, "B", NA),
    probability = c(1, 1, 0.5)
  )
  pool <- hee_update_species_pool(init, events, times = c(10, 5, 0))
  expect_true(all(pool$species_pool_size >= 0))
  expect_gt(pool$species_pool_size[pool$region == "A" & pool$time_ma == 5],
            pool$species_pool_size[pool$region == "A" & pool$time_ma == 10])
  expect_gt(pool$species_pool_size[pool$region == "B" & pool$time_ma == 5],
            pool$species_pool_size[pool$region == "B" & pool$time_ma == 10])
})

test_that("dynamic assembly enforces geographic and lineage masks", {
  suit <- data.frame(species = "sp1", cell_id = "c1", time_ma = c(10, 5, 0),
                     suitability = c(0.8, 0.7, 0.6))
  prev <- data.frame(species = "sp1", cell_id = "c1", time_ma = 10,
                     occupancy_probability = 0.6)
  dyn <- hee_dynamic_assembly(
    previous_state = prev,
    suitability = suit,
    colonisation = data.frame(species = "sp1", cell_id = "c1",
                              time_ma = c(5, 0),
                              source_pressure = c(0.5, 0.5)),
    geo_mask = data.frame(cell_id = "c1", time_ma = 5, land = 0),
    phylo_mask = data.frame(species = "sp1", time_ma = 0, E_phylo = 0)
  )
  expect_equal(dyn$occupancy_probability[dyn$time_ma == 5], 0)
  expect_equal(dyn$occupancy_probability[dyn$time_ma == 0], 0)
  expect_true(all(dyn$occupancy_probability >= 0 & dyn$occupancy_probability <= 1))
})

test_that("supplied dynamic geo and phylo masks are conservative for unmatched rows", {
  suit <- data.frame(species = "sp1", cell_id = "c1", time_ma = c(10, 0),
                     suitability = c(1, 1))
  prev <- data.frame(species = "sp1", cell_id = "c1", time_ma = 10,
                     occupancy_probability = 1)
  dyn_geo <- hee_dynamic_assembly(
    previous_state = prev,
    suitability = suit,
    geo_mask = data.frame(cell_id = "c1", time_ma = 10,
                          geographic_existence = 1)
  )
  expect_equal(dyn_geo$geographic_existence[dyn_geo$time_ma == 0], 0)
  expect_equal(dyn_geo$occupancy_probability[dyn_geo$time_ma == 0], 0)

  dyn_phy <- hee_dynamic_assembly(
    previous_state = prev,
    suitability = suit,
    phylo_mask = data.frame(species = "sp1", time_ma = 10, E_phylo = 1)
  )
  expect_equal(dyn_phy$lineage_exists[dyn_phy$time_ma == 0], 0)
  expect_equal(dyn_phy$occupancy_probability[dyn_phy$time_ma == 0], 0)

  expect_error(
    hee_dynamic_assembly(previous_state = prev, suitability = suit,
                         geo_mask = data.frame(cell_id = "c1", time_ma = 0,
                                               bad = 1)),
    "geo_mask must contain"
  )
  expect_error(
    hee_dynamic_assembly(previous_state = prev, suitability = suit,
                         phylo_mask = data.frame(species = "sp1", time_ma = 0,
                                                 bad = 1)),
    "phylo_mask must contain"
  )
})

test_that("dynamic assembly is sequential, not repeated from the initial state", {
  suit <- data.frame(species = "sp1", cell_id = "c1", time_ma = c(10, 5, 0),
                     suitability = c(1, 1, 1))
  prev <- data.frame(species = "sp1", cell_id = "c1", time_ma = 10,
                     occupancy_probability = 1)
  ext <- data.frame(species = "sp1", cell_id = "c1", time_ma = c(5, 0),
                    extinction_probability = c(0.5, 0))
  col <- data.frame(species = "sp1", cell_id = "c1", time_ma = c(5, 0),
                    colonisation_probability = c(0, 0))
  dyn <- hee_dynamic_assembly(prev, suit, colonisation = col, extinction = ext)
  expect_equal(dyn$occupancy_probability[dyn$time_ma == 5], 0.5)
  expect_equal(dyn$occupancy_probability[dyn$time_ma == 0], 0.5)
})

test_that("dynamic assembly converts source pressure with the colonisation formula and time interval", {
  suit <- data.frame(species = "sp1", cell_id = "c1", time_ma = c(10, 0),
                     suitability = c(0.8, 0.8))
  prev <- data.frame(species = "sp1", cell_id = "c1", time_ma = 10,
                     occupancy_probability = 0)
  src <- data.frame(species = "sp1", cell_id = "c1", time_ma = c(10, 0),
                    source_pressure = c(1, 1))
  dyn <- hee_dynamic_assembly(
    previous_state = prev,
    suitability = suit,
    colonisation = src,
    extinction = data.frame(species = "sp1", cell_id = "c1",
                            time_ma = c(10, 0),
                            extinction_probability = 0)
  )
  expected <- hee_colonisation_probability(
    data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
               suitability = 0.8),
    source_pressure = 1,
    delta_t = 10
  )$colonisation_probability
  expect_equal(dyn$interval_myr[dyn$time_ma == 0], 10)
  expect_equal(dyn$colonisation_probability[dyn$time_ma == 0], expected)
  expect_equal(dyn$occupancy_probability[dyn$time_ma == 0], expected)
})

test_that("missing colonisation process is neutral in dynamic assembly", {
  suit <- data.frame(species = "sp1", cell_id = "c1", time_ma = c(10, 0),
                     suitability = c(1, 1))
  prev <- data.frame(species = "sp1", cell_id = "c1", time_ma = 10,
                     occupancy_probability = 0)
  ext <- data.frame(species = "sp1", cell_id = "c1", time_ma = c(10, 0),
                    extinction_probability = c(0, 0))
  dyn <- hee_dynamic_assembly(previous_state = prev, suitability = suit,
                              extinction = ext)
  expect_equal(dyn$occupancy_probability[dyn$time_ma == 0], 0)
})

test_that("missing extinction process is neutral in dynamic assembly", {
  suit <- data.frame(species = "sp1", cell_id = "c1", time_ma = c(10, 0),
                     suitability = c(0.2, 0.2))
  prev <- data.frame(species = "sp1", cell_id = "c1", time_ma = 10,
                     occupancy_probability = 1)
  col <- data.frame(species = "sp1", cell_id = "c1", time_ma = c(10, 0),
                    colonisation_probability = c(0, 0))
  dyn <- hee_dynamic_assembly(previous_state = prev, suitability = suit,
                              colonisation = col)
  expect_equal(dyn$occupancy_probability[dyn$time_ma == 0], 1)
})

test_that("dynamic assembly rejects duplicate times and duplicate optional joins", {
  dup_time <- data.frame(species = "sp1", cell_id = "c1",
                         time_ma = c(5, 5), suitability = c(0.5, 0.6))
  expect_error(hee_dynamic_assembly(suitability = dup_time), "Duplicate time")

  suit <- data.frame(species = "sp1", cell_id = "c1", time_ma = c(5, 0),
                     suitability = 0.8)
  dup_col <- data.frame(species = "sp1", cell_id = "c1", time_ma = 5,
                        colonisation_probability = c(0.2, 0.8))
  expect_error(hee_dynamic_assembly(suitability = suit, colonisation = dup_col),
               "Duplicate key")
})

test_that("mechanism attribution follows scientific priority order", {
  x <- data.frame(
    geographic_existence = c(0, 1, 1, 1),
    lineage_exists = c(1, 0, 1, 1),
    suitability = c(0.8, 0.8, 0.1, 0.8),
    source_pressure = c(0.5, 0.5, 0.5, 0.01),
    connectivity = c(0.8, 0.8, 0.8, 0.8),
    extinction_probability = c(0.1, 0.1, 0.1, 0.1)
  )
  z <- hee_mechanism_attribution(x)
  expect_equal(z$dominant_mechanism[1], "arena_absent_or_lost")
  expect_equal(z$dominant_mechanism[2], "lineage_absent")
  expect_equal(z$dominant_mechanism[3], "niche_limited")
  expect_equal(z$dominant_mechanism[4], "source_limited")
})

test_that("accessibility, static combine, limitation maps, and reliability are callable", {
  A <- data.frame(species = "sp1", region = "R1", time_ma = 0,
                  accessibility = 0.7)
  cells <- data.frame(cell_id = "c1", region = "R1", time_ma = 0,
                      geographic_existence = 1)
  Ac <- hee_accessibility_to_cells(A, cells)
  expect_equal(Ac$cell_id, "c1")

  S <- data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                  suitability = 0.8)
  P <- hee_combine_static(
    S,
    accessibility = Ac,
    geo_mask = cells,
    static_filter = data.frame(cell_id = "c1", time_ma = 0, D_static = 0.5)
  )
  expect_equal(P$probability, 0.8 * 0.7 * 0.5)

  P_na_geo <- hee_combine_static(
    S,
    geo_mask = data.frame(cell_id = "c1", time_ma = 0,
                          geographic_existence = NA_real_)
  )
  expect_equal(P_na_geo$land_weight, 0)
  expect_equal(P_na_geo$probability, 0)
  expect_error(
    hee_combine_static(
      S,
      geo_mask = data.frame(cell_id = "c1", time_ma = 0,
                            unrelated_geo = 1)
    ),
    "geo_mask must contain"
  )

  P2 <- hee_combine_static(
    S,
    dynamic_filter = data.frame(cell_id = "c1", time_ma = 0, D_dynamic = 0.25)
  )
  expect_equal(P2$probability, 0.8 * 0.25)

  att <- hee_mechanism_attribution(data.frame(cell_id = "c1", time_ma = 0,
                                              geographic_existence = 1,
                                              lineage_exists = 1,
                                              suitability = 0.1))
  lim <- hee_limitation_map(att)
  expect_equal(lim$dominant_community_mechanism, "niche_limited")

  rel <- hee_prediction_reliability(
    extrapolation = data.frame(cell_id = "c1", time_ma = 0,
                               extrapolation_score = 0.2),
    model_agreement = data.frame(cell_id = "c1", time_ma = 0,
                                 model_agreement = 0.5)
  )
  expect_equal(rel$prediction_reliability, 0.8 * 0.5)

  rel_vec <- hee_prediction_reliability(model_agreement = c(0.5, 0.25))
  expect_equal(rel_vec$prediction_reliability, c(0.5, 0.25))
})

test_that("prediction reliability handles unbounded risk, NA components, and weights", {
  rel <- hee_prediction_reliability(
    extrapolation = data.frame(cell_id = "c1", time_ma = 0,
                               extrapolation_score = 3),
    model_agreement = data.frame(cell_id = "c1", time_ma = 0,
                                 model_agreement = NA_real_)
  )
  expect_gt(rel$prediction_reliability, 0)
  expect_lt(rel$prediction_reliability, 1)
  expect_equal(rel$model_agreement, 1)

  rel2 <- hee_prediction_reliability(
    extrapolation = c(0.5, 0.5),
    model_agreement = c(0.1, 0.9),
    weights = c(extrapolation = 1, model_agreement = 0)
  )
  expect_equal(rel2$prediction_reliability, c(0.5, 0.5))

  expect_error(
    hee_prediction_reliability(
      model_agreement = data.frame(cell_id = "c1", time_ma = 0,
                                   bad_weight = 0.5)
    ),
    "model_agreement must contain"
  )
})

test_that("accessibility to cells rejects duplicated cell-time-region keys", {
  A <- data.frame(species = "sp1", region = "R1", time_ma = 0,
                  accessibility = 0.7)
  cells <- data.frame(cell_id = "c1", region = "R1", time_ma = 0,
                      lon = c(1, 2), lat = 0)
  expect_error(hee_accessibility_to_cells(A, cells), "Duplicate key")
})

test_that("cell-level accessibility validates keys and bounds before returning", {
  cells <- data.frame(cell_id = "c1", region = "R1", time_ma = 0)
  dup <- data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                    accessibility = c(2, 2))
  expect_error(hee_accessibility_to_cells(dup, cells), "Duplicate key")

  A <- data.frame(species = "sp1", cell_id = c("c1", "c2"), time_ma = 0,
                  accessibility = c(2, NA_real_))
  out <- hee_accessibility_to_cells(A, cells)
  expect_equal(out$accessibility, c(1, 0))
})

test_that("geoprocess diagnostics rejects duplicated region lookup keys", {
  land <- data.frame(cell_id = c("c1", "c1"), region = c("R1", "R2"),
                     time_ma = c(0, 0), land_mask_dem = 1,
                     land_area_km2 = 1, bio1 = c(1, 2))
  proj <- data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                     probability = 0.4, suitability = 0.5)
  expect_error(hee_geoprocess_diagnostics(land, proj, env_cols = "bio1"),
               "Duplicate key")
})

test_that("geoprocess diagnostics respects custom region and time columns in reliability joins", {
  land <- data.frame(cell_id = "c1", realm = "R1", age = 0,
                     land_mask_dem = 1, land_area_km2 = 1, bio1 = 1)
  proj <- data.frame(species = "sp1", cell_id = "c1", realm = "R1", age = 0,
                     probability = 0.4, suitability = 0.5)
  ex <- data.frame(cell_id = "c1", age = 0, extrapolation_score = 0.4)
  out <- hee_geoprocess_diagnostics(
    landscape_state = land,
    projection = proj,
    extrapolation = ex,
    env_cols = "bio1",
    region_col = "realm",
    time_col = "age"
  )
  expect_equal(out$reliability$model_agreement, 0.9)
  expect_true("realm" %in% names(out$reliability))
  expect_true("age" %in% names(out$reliability))
})

test_that("complete geological process diagnostics returns expected layers", {
  land <- expand.grid(cell_id = c("c1", "c2"), time_ma = c(10, 5, 0),
                      stringsAsFactors = FALSE)
  land$region <- ifelse(land$cell_id == "c1", "R1", "R2")
  land$land_mask_dem <- 1
  land$land_area_km2 <- ifelse(land$region == "R1",
                               c(8, 10, 12, 8, 10, 12)[seq_len(nrow(land))],
                               c(5, 7, 8, 5, 7, 8)[seq_len(nrow(land))])
  land$bio1 <- seq_len(nrow(land))
  proj <- expand.grid(species = c("sp1", "sp2"),
                      cell_id = c("c1", "c2"),
                      time_ma = c(10, 5, 0),
                      stringsAsFactors = FALSE)
  proj <- merge(proj, land[, c("cell_id", "region", "time_ma")],
                by = c("cell_id", "time_ma"), all.x = TRUE)
  proj$suitability <- 0.7
  proj$probability <- 0.3
  proj$accessibility <- 0.8
  proj$phylo_existence <- 1
  conn <- expand.grid(from_region = c("R1", "R2"),
                      to_region = c("R1", "R2"),
                      time_ma = c(10, 5, 0),
                      stringsAsFactors = FALSE)
  conn$connectivity <- ifelse(conn$from_region == conn$to_region, 1, 0.35)
  out <- hee_geoprocess_diagnostics(
    landscape_state = land,
    projection = proj,
    connectivity_cube = conn,
    phylo_mask = matrix(1, nrow = 2, ncol = 3,
                        dimnames = list(c("sp1", "sp2"), paste0(c(10, 5, 0), "Ma"))),
    tip_ranges = matrix(1, nrow = 2, ncol = 2,
                        dimnames = list(c("sp1", "sp2"), c("R1", "R2"))),
    times = c(10, 5, 0),
    env_cols = "bio1"
  )
  expect_true(all(c("region_state", "source_pressure",
                    "colonisation_probability", "extinction_probability",
                    "speciation_opportunity", "dynamic_occupancy",
                    "mechanism_attribution", "species_pool",
                    "process_summary") %in% names(out)))
  expect_true(all(out$colonisation_probability$colonisation_probability >= 0 &
                    out$colonisation_probability$colonisation_probability <= 1))
  expect_true(all(out$extinction_probability$extinction_probability >= 0 &
                    out$extinction_probability$extinction_probability <= 1))
  expect_gt(nrow(out$process_summary), 0)
})

test_that("geoprocess diagnostics are conservative for unmatched accessibility and phylo masks", {
  land <- data.frame(cell_id = c("c1", "c2"), region = c("R1", "R2"),
                     time_ma = 0, land_mask_dem = 1,
                     land_area_km2 = c(10, 10), bio1 = c(1, 2))
  proj <- data.frame(species = "sp1", cell_id = c("c1", "c2"),
                     region = c("R1", "R2"), time_ma = 0,
                     suitability = 0.7, probability = 0.3)
  acc <- data.frame(species = "sp1", region = "R1", time_ma = 0,
                    accessibility = 0.5)
  phy <- matrix(1, nrow = 1, ncol = 1, dimnames = list("sp1", "5Ma"))
  out <- hee_geoprocess_diagnostics(
    landscape_state = land,
    projection = proj,
    accessibility = acc,
    phylo_mask = phy,
    times = 0,
    env_cols = "bio1"
  )
  r2 <- out$region_state[out$region_state$region == "R2", ]
  expect_equal(r2$accessibility, 0)
  expect_equal(out$region_state$phylo_existence, c(0, 0))
})
