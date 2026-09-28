expect_bounded01 <- function(x, allow_na = FALSE) {
  x <- suppressWarnings(as.numeric(x))
  if (!allow_na) {
    expect_false(anyNA(x))
  }
  ok <- if (allow_na) !is.na(x) else rep(TRUE, length(x))
  expect_true(all(is.finite(x[ok])))
  expect_true(all(x[ok] >= 0 & x[ok] <= 1))
}

expect_unique_keys <- function(x, keys) {
  expect_false(any(duplicated(as.data.frame(x)[, keys, drop = FALSE])))
}

test_that("geological process inputs handle empty, zero-row, one-row, mixed types, and extra columns", {
  expect_error(hee_land_age(data.frame()), "missing")
  expect_error(hee_colonisation_probability(data.frame()), "suitability")

  zero_land <- data.frame(cell_id = character(), time_ma = numeric(),
                          land = numeric(), unused = character())
  zero_age <- hee_land_age(zero_land)
  expect_equal(nrow(zero_age), 0)

  zero_suit <- data.frame(species = character(), cell_id = character(),
                          time_ma = numeric(), suitability = numeric())
  expect_equal(nrow(hee_colonisation_probability(zero_suit)), 0)
  expect_equal(nrow(hee_extinction_probability(zero_suit)), 0)

  one <- data.frame(cell_id = factor("c1"), time_ma = "540",
                    land = 1L, extra = "ignored")
  age <- hee_land_age(one)
  expect_equal(age$geographic_existence, 1)
  expect_equal(age$land_age_myr, 0)

  mixed <- data.frame(region = factor("R1"), cell_id = 1L,
                      time_ma = as.integer(0), cell_area_km2 = 5,
                      novelty = "3", extra = "ignored")
  opp <- hee_ecological_opportunity(mixed, novelty_col = "novelty",
                                    weights = c(area = 0, heterogeneity = 0,
                                                novelty = 1, young_land = 0))
  expect_gt(opp$ecological_opportunity, 0)
  expect_bounded01(opp$ecological_opportunity)
})

test_that("geological process inputs reject missing columns, invalid numeric values, and duplicate keys", {
  expect_error(hee_landscape_events(data.frame(region = "R1", land = 1)),
               "missing")
  expect_error(hee_land_age(data.frame(cell_id = "c1", time_ma = -1, land = 1)),
               "negative")
  expect_error(hee_land_age(data.frame(cell_id = "c1", time_ma = Inf, land = 1)),
               "finite")
  expect_error(hee_land_age(data.frame(cell_id = "c1", time_ma = c(0, 0),
                                       land = c(1, 1))),
               "Duplicate key")
  expect_error(hee_ecological_opportunity(data.frame(cell_id = "c1",
                                                     time_ma = c(0, 0),
                                                     cell_area_km2 = c(1, 2))),
               "Duplicate key")
  expect_error(hee_structural_connectivity(data.frame(from_region = "A",
                                                      to_region = "B",
                                                      time_ma = c(0, 0),
                                                      paleodistance = c(1, 2))),
               "Duplicate key")
})

test_that("time ordering and unequal intervals are explicit in land age and dynamic assembly", {
  land_old_to_new <- data.frame(cell_id = "c1", time_ma = c(540, 500, 0),
                                land = c(1, 1, 1))
  land_new_to_old <- land_old_to_new[c(3, 2, 1), ]
  a <- hee_land_age(land_old_to_new)
  b <- hee_land_age(land_new_to_old)
  a <- a[order(a$time_ma, decreasing = TRUE), ]
  b <- b[order(b$time_ma, decreasing = TRUE), ]
  expect_equal(a$land_age_myr, c(0, 40, 540))
  expect_equal(a$land_age_myr, b$land_age_myr)

  suit <- data.frame(species = "sp1", cell_id = "c1",
                     time_ma = c(0, 540, 500), suitability = 1)
  prev <- data.frame(species = "sp1", cell_id = "c1", time_ma = 540,
                     occupancy_probability = 1)
  ext <- data.frame(species = "sp1", cell_id = "c1",
                    time_ma = c(540, 500, 0),
                    extinction_probability = c(0, 0.5, 0))
  col <- data.frame(species = "sp1", cell_id = "c1",
                    time_ma = c(540, 500, 0),
                    colonisation_probability = 0)
  dyn <- hee_dynamic_assembly(prev, suit, colonisation = col, extinction = ext)
  dyn <- dyn[order(dyn$time_ma, decreasing = TRUE), ]
  expect_equal(dyn$occupancy_probability[1], 1)
  expect_equal(dyn$occupancy_probability[2], 0.5)
  expect_equal(dyn$occupancy_probability[3], 0.5)
  expect_bounded01(dyn$occupancy_probability)

  gamma_dt <- hee_colonisation_probability(
    data.frame(suitability = 0.8, delta = c(1, 1000)),
    source_pressure = c(0.8, 0.8),
    accessibility = c(1, 1),
    delta_t = "delta"
  )
  expect_lt(gamma_dt$colonisation_probability[1],
            gamma_dt$colonisation_probability[2])
  expect_error(hee_extinction_probability(data.frame(suitability = 1),
                                          delta_t = NaN),
               "delta_t")
})

test_that("probability-like geological process outputs stay bounded under extreme inputs", {
  x <- data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                  suitability = c(0, 1, NA, Inf, -Inf, NaN))
  gamma <- hee_colonisation_probability(
    x,
    source_pressure = c(0, 1, NA, Inf, -Inf, NaN),
    accessibility = c(1, 1, NA, Inf, -Inf, NaN),
    ecological_opportunity = c(0, 1, NA, Inf, -Inf, NaN),
    delta_t = 1
  )
  expect_bounded01(gamma$colonisation_probability)

  eps <- hee_extinction_probability(
    data.frame(suitability = c(0, 1, NA, Inf, -Inf, NaN),
               land_loss = c(0, 1, NA, Inf, -Inf, NaN),
               isolation = c(0, 1, NA, Inf, -Inf, NaN),
               rescue_effect = c(0, 1, NA, Inf, -Inf, NaN),
               area_km2 = c(0, 1, NA, Inf, -Inf, NaN)),
    delta_t = 1
  )
  expect_bounded01(eps$extinction_probability)

  s <- hee_structural_connectivity(
    data.frame(from_region = "A", to_region = paste0("B", 1:6), time_ma = 0,
               paleodistance = c(0, 1e9, NA, Inf, -Inf, NaN),
               barrier_strength = c(0, 1, NA, Inf, -Inf, NaN),
               route_open = c(1, 0, NA, Inf, -Inf, NaN))
  )
  expect_bounded01(s$structural_connectivity)
  s_unknown <- hee_structural_connectivity(
    data.frame(from_region = "A", to_region = c("B", "C"), time_ma = 0,
               paleodistance = 0,
               barrier_strength = c(NA_real_, 0),
               route_open = c(1, NA_real_))
  )
  expect_lt(s_unknown$structural_connectivity[1], 1)
  expect_equal(s_unknown$structural_connectivity[2], 0)

  f <- hee_functional_connectivity(
    data.frame(from_region = "A", to_region = paste0("B", 1:6), time_ma = 0,
               species = "sp1",
               paleodistance = c(0, 1e9, NA, Inf, -Inf, NaN),
               barrier_strength = c(0, 1, NA, Inf, -Inf, NaN),
               structural_connectivity = 1),
    traits = data.frame(species = "sp1", dispersal = 5)
  )
  expect_bounded01(f$functional_connectivity)
  expect_bounded01(f$connectivity)
  f_unknown <- hee_functional_connectivity(
    data.frame(from_region = "A", to_region = c("B", "C"), time_ma = 0,
               species = "sp1", paleodistance = 0,
               barrier_strength = c(NA_real_, 0),
               structural_connectivity = c(1, NA_real_)),
    traits = data.frame(species = "sp1", dispersal = 5)
  )
  expect_lte(f_unknown$functional_connectivity[1],
             f_unknown$functional_connectivity[2])
  expect_equal(f_unknown$connectivity[2], 0)

  clim <- hee_climatic_connectivity(
    data.frame(from_region = "A", to_region = "B", time_ma = 0,
               environment_distance = c(0, 1e9, NA, Inf, -Inf, NaN),
               structural_connectivity = 1,
               functional_connectivity = 1)
  )
  expect_bounded01(clim$climate_connectivity)
  expect_bounded01(clim$connectivity)

  sp <- hee_speciation_opportunity(
    data.frame(isolation_duration_ma = c(0, 1, 1e9, NA, Inf, -Inf, NaN),
               niche_divergence = c(0, 1, 0.5, NA, Inf, -Inf, NaN),
               area_km2 = c(0, 1, 1e9, NA, Inf, -Inf, NaN),
               habitat_heterogeneity = c(0, 1, 0.5, NA, Inf, -Inf, NaN),
               population_persistence = c(0, 1, 0.5, NA, Inf, -Inf, NaN),
               source_pressure = c(0, 1, 0.5, NA, Inf, -Inf, NaN),
               connectivity = c(0, 1, 0.5, NA, Inf, -Inf, NaN),
               ecological_opportunity = c(0, 1, 0.5, NA, Inf, -Inf, NaN))
  )
  expect_bounded01(sp$speciation_opportunity)

  rel <- hee_prediction_reliability(extrapolation = c(0, 1, 1e9, NA, Inf, -Inf, NaN),
                                    model_agreement = c(0, 1, NA, Inf, -Inf, NaN, 0.5))
  expect_bounded01(rel$prediction_reliability)
})

test_that("scientific directionality is monotone for colonisation, connectivity, extinction, and speciation", {
  base <- data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                     suitability = c(0.2, 0.8))
  by_suit <- hee_colonisation_probability(base, source_pressure = 1,
                                          accessibility = 1, delta_t = 5)
  expect_lte(by_suit$colonisation_probability[1],
             by_suit$colonisation_probability[2])

  by_src <- hee_colonisation_probability(
    data.frame(suitability = c(0.8, 0.8)),
    source_pressure = c(0.1, 0.9),
    accessibility = 1,
    delta_t = 5
  )
  expect_lte(by_src$colonisation_probability[1],
             by_src$colonisation_probability[2])

  by_acc <- hee_colonisation_probability(
    data.frame(suitability = c(0.8, 0.8)),
    source_pressure = 0.8,
    accessibility = c(0.2, 0.9),
    delta_t = 5
  )
  expect_lte(by_acc$colonisation_probability[1],
             by_acc$colonisation_probability[2])

  dist <- hee_structural_connectivity(
    data.frame(from_region = "A", to_region = c("B1", "B2"), time_ma = 0,
               paleodistance = c(1, 10), barrier_strength = 0)
  )
  expect_gte(dist$structural_connectivity[1], dist$structural_connectivity[2])

  barr <- hee_structural_connectivity(
    data.frame(from_region = "A", to_region = c("B1", "B2"), time_ma = 0,
               paleodistance = 1, barrier_strength = c(0, 0.9))
  )
  expect_gte(barr$structural_connectivity[1], barr$structural_connectivity[2])

  func <- hee_functional_connectivity(
    data.frame(from_region = "A", to_region = "B", time_ma = 0,
               species = c("poor", "good"), paleodistance = 5,
               barrier_strength = 0),
    traits = data.frame(species = c("poor", "good"), dispersal = c(1, 10))
  )
  expect_lte(func$functional_connectivity[func$species == "poor"],
             func$functional_connectivity[func$species == "good"])

  climate <- hee_climatic_connectivity(
    data.frame(from_region = "A", to_region = "B", time_ma = 0,
               environment_distance = c(1, 10))
  )
  expect_gte(climate$climate_connectivity[1], climate$climate_connectivity[2])

  rescue <- hee_extinction_probability(
    data.frame(suitability = 0.2, isolation = 0.8, rescue_effect = c(0, 1)),
    delta_t = 5
  )
  expect_gte(rescue$extinction_probability[1], rescue$extinction_probability[2])

  allo <- hee_allopatric_speciation_opportunity(
    data.frame(range_fragmentation = 1,
               isolation_duration_ma = c(0, 10),
               population_persistence = 1)
  )
  expect_lte(allo$allopatric_opportunity[1], allo$allopatric_opportunity[2])
})

test_that("functional connectivity is neutral without traits and errors on missing species traits", {
  structural <- data.frame(from_region = "A", to_region = "B", time_ma = 0,
                           species = "sp1", paleodistance = 0,
                           barrier_strength = 0,
                           structural_connectivity = 0.8)
  no_traits <- hee_functional_connectivity(structural, traits = NULL)
  expect_equal(no_traits$functional_connectivity, 1)
  expect_equal(no_traits$connectivity, structural$structural_connectivity)

  expect_error(
    hee_functional_connectivity(
      structural,
      traits = data.frame(species = "other_sp", dispersal = 10)
    ),
    "missing dispersal values"
  )
})

test_that("zero area increases extinction vulnerability and extreme speciation rates stay finite", {
  eps <- hee_extinction_probability(
    data.frame(suitability = 1, area_km2 = c(0, 1000),
               land_loss = 0, isolation = 0, rescue_effect = 0),
    delta_t = 1
  )
  expect_equal(eps$extinction_probability[1], 1)
  expect_true(eps$forced_geographic_extinction[1])
  expect_gt(eps$extinction_probability[1], eps$extinction_probability[2])

  sp <- hee_speciation_opportunity(
    data.frame(isolation_duration_ma = Inf, niche_divergence = 1,
               area_km2 = 1e12, habitat_heterogeneity = 1,
               population_persistence = 1, ecological_opportunity = 1)
  )
  expect_true(is.finite(sp$speciation_rate_lambda))
})

test_that("species-pool updates use the oldest initial row independent of row order", {
  init <- data.frame(
    region = c("A", "A", "B", "B"),
    time_ma = c(0, 10, 0, 10),
    species_pool_size = c(1, 5, 2, 6)
  )
  events <- data.frame(time_ma = 0, event_type = "extinction",
                       region = "A", probability = 1)
  out1 <- hee_update_species_pool(init, events = events, times = c(10, 0))
  out2 <- hee_update_species_pool(init[c(3, 1, 4, 2), ],
                                  events = events, times = c(10, 0))
  expect_equal(out1, out2)
  expect_equal(out1$species_pool_size[out1$region == "A" & out1$time_ma == 10], 5)
  expect_equal(out1$species_pool_size[out1$region == "A" & out1$time_ma == 0], 4)
})

test_that("occupancy update obeys geographic, colonisation, and extinction constraints exactly", {
  expect_equal(hee_update_occupancy(p_prev = 0.9, gamma = 1, epsilon = 0,
                                    phylo_mask = 1, geo_mask = 0), 0)
  expect_equal(hee_update_occupancy(p_prev = 0, gamma = 0, epsilon = 0), 0)
  expect_equal(hee_update_occupancy(p_prev = 1, gamma = 0, epsilon = 0), 1)
  expect_equal(hee_update_occupancy(p_prev = 0, gamma = 1, epsilon = 0), 1)
  expect_equal(hee_update_occupancy(p_prev = 1, gamma = 0, epsilon = 1), 0)

  tab <- hee_update_occupancy(
    data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
               previous_probability = NA_real_),
    gamma = data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                       colonisation_probability = 1),
    epsilon = data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                         extinction_probability = 0)
  )
  expect_equal(tab$occupancy_probability, 1)
})

test_that("species, connectivity, and accessibility joins preserve intended keys and row counts", {
  conn <- data.frame(from_region = c("A", "A", "B"),
                     to_region = c("B", "C", "A"),
                     time_ma = 0,
                     species = c("sp2", "sp1", "sp3"),
                     paleodistance = c(2, 3, 4),
                     structural_connectivity = c(0.5, 0.4, 0.3))
  traits <- data.frame(species = c("sp1", "sp2", "extra"),
                       dispersal = c(5, 10, 100))
  expect_error(hee_functional_connectivity(conn, traits = traits),
               "missing dispersal values")
  traits <- rbind(traits, data.frame(species = "sp3", dispersal = 3))
  out <- hee_functional_connectivity(conn, traits = traits)
  expect_equal(nrow(out), nrow(conn))
  expect_equal(out$species, conn$species)
  expect_false("extra" %in% out$species)
  expect_false(any(is.na(out$from_region)))

  conn_no_species <- conn[, setdiff(names(conn), "species")]
  expanded <- hee_functional_connectivity(conn_no_species, traits = traits,
                                          expand_species = TRUE)
  expect_equal(nrow(expanded), nrow(conn_no_species) * nrow(traits))
  expect_error(hee_functional_connectivity(conn_no_species, traits = traits,
                                           expand_species = FALSE),
               "expand_species")

  expect_error(
    hee_functional_connectivity(conn,
                                traits = data.frame(species = c("sp1", "sp1"),
                                                    dispersal = c(1, 2))),
    "Duplicate key"
  )

  A <- expand.grid(species = c("sp1", "sp2"), region = c("R1", "R2"),
                   time_ma = 0, KEEP.OUT.ATTRS = FALSE)
  A$accessibility <- c(1, 0.5, 0.25, 0.75)
  cells <- data.frame(cell_id = c("c1", "c2", "c3"),
                      region = c("R1", "R2", "R2"),
                      time_ma = 0)
  Ac <- hee_accessibility_to_cells(A, cells)
  expect_equal(nrow(Ac), 6)
  expect_unique_keys(Ac, c("species", "cell_id", "time_ma"))
  expect_false(any(is.na(Ac$cell_id)))
})

test_that("minimal two-species three-cell workflow preserves keys, masks, and probability ranges", {
  times <- c(540, 500, 0)
  cells <- data.frame(cell_id = c("c1", "c2", "c3"),
                      region = c("A", "B", "B"),
                      lon = c(0, 10, 20),
                      lat = 0,
                      stringsAsFactors = FALSE)
  suit <- expand.grid(species = c("sp1", "sp2"), cell_id = cells$cell_id,
                      time_ma = times, KEEP.OUT.ATTRS = FALSE)
  suit <- merge(suit, cells, by = "cell_id", all.x = TRUE, sort = FALSE)
  suit$suitability <- ifelse(suit$species == "sp1", 0.8, 0.6)

  pairs <- expand.grid(from_region = c("A", "B"),
                       to_region = c("A", "B"),
                       time_ma = times, KEEP.OUT.ATTRS = FALSE)
  pairs$paleodistance <- ifelse(pairs$from_region == pairs$to_region, 0, 5)
  pairs$barrier_strength <- ifelse(pairs$time_ma == 500 &
                                     pairs$to_region == "B", 0.5, 0)
  structural <- hee_structural_connectivity(pairs)
  traits <- data.frame(species = c("sp1", "sp2"), dispersal = c(10, 2))
  functional <- hee_functional_connectivity(structural, traits = traits,
                                            expand_species = TRUE)
  climatic <- hee_climatic_connectivity(functional,
                                        corridor_suitability_col = NULL)
  expect_true(all(c("structural_connectivity", "functional_connectivity",
                    "connectivity") %in% names(climatic)))

  prev_region <- expand.grid(species = c("sp1", "sp2"), region = c("A", "B"),
                             time_ma = 540, KEEP.OUT.ATTRS = FALSE)
  prev_region$occupancy_probability <- ifelse(prev_region$region == "A", 1, 0)
  target_region <- expand.grid(species = c("sp1", "sp2"), region = c("A", "B"),
                               time_ma = times, KEEP.OUT.ATTRS = FALSE)
  src <- hee_colonisation_pressure(prev_region, climatic,
                                   target = target_region)
  gamma <- hee_colonisation_probability(
    data.frame(species = src$species, region = src$region,
               time_ma = src$time_ma, suitability = 0.8),
    source_pressure = src,
    delta_t = 40,
    keys = c("species", "region", "time_ma")
  )
  rescue <- hee_rescue_effect(prev_region, climatic, target = target_region)
  eps <- hee_extinction_probability(
    data.frame(species = rescue$species, region = rescue$region,
               time_ma = rescue$time_ma, suitability = 0.8,
               rescue_effect = rescue$rescue_effect),
    delta_t = 40
  )
  expect_bounded01(src$source_pressure)
  expect_bounded01(rescue$rescue_effect)
  expect_bounded01(gamma$colonisation_probability)
  expect_bounded01(eps$extinction_probability)

  geo <- data.frame(region = c("A", "B", "A", "B", "A", "B"),
                    time_ma = rep(times, each = 2),
                    geographic_existence = c(1, 1, 1, 0, 1, 1))
  dyn_suit <- unique(suit[, c("species", "region", "time_ma", "suitability")])
  dyn <- hee_dynamic_assembly(
    suitability = dyn_suit,
    previous_state = prev_region,
    colonisation = gamma[, c("species", "region", "time_ma",
                             "colonisation_probability")],
    extinction = eps[, c("species", "region", "time_ma",
                         "extinction_probability")],
    geo_mask = geo,
    species_col = "species",
    cell_col = "region",
    time_col = "time_ma"
  )
  expect_equal(dyn$occupancy_probability[dyn$region == "B" &
                                           dyn$time_ma == 500], c(0, 0))
  expect_unique_keys(dyn, c("species", "region", "time_ma"))
  expect_bounded01(dyn$occupancy_probability)

  A <- data.frame(species = rep(c("sp1", "sp2"), each = 2),
                  region = rep(c("A", "B"), times = 2),
                  time_ma = 0,
                  accessibility = c(1, 0.5, 1, 0.5))
  cell_time0 <- cells[, c("cell_id", "region")]
  cell_time0$time_ma <- 0
  Ac <- hee_accessibility_to_cells(A, cell_time0)
  static <- data.frame(cell_id = cells$cell_id, time_ma = 0,
                       D_static = c(1, 0.5, 0.5),
                       weight = c(1, 0.5, 0.5))
  P <- hee_combine_static(suit[suit$time_ma == 0, ],
                          accessibility = Ac,
                          geo_mask = data.frame(cell_id = cells$cell_id,
                                                time_ma = 0,
                                                geographic_existence = 1),
                          static_filter = static)
  expect_unique_keys(P, c("species", "cell_id", "time_ma"))
  expect_bounded01(P$probability)
  expect_true(all(c("D_static", "land_weight", "accessibility") %in% names(P)))

  rel <- hee_prediction_reliability(
    extrapolation = data.frame(cell_id = cells$cell_id, time_ma = 0,
                               extrapolation_score = c(0, 0.2, 2)),
    model_agreement = data.frame(cell_id = cells$cell_id, time_ma = 0,
                                 model_agreement = c(1, 0.8, 0.6)),
    phylo_validity = data.frame(cell_id = cells$cell_id, time_ma = 0,
                                phylo_validity = 1)
  )
  expect_unique_keys(rel, c("cell_id", "time_ma"))
  expect_bounded01(rel$prediction_reliability)
})

test_that("regressions keep neutral defaults, explicit static weights, scalar reliability, and order invariance", {
  no_access <- hee_colonisation_probability(data.frame(suitability = 0.8),
                                            source_pressure = 0.5,
                                            delta_t = 5)
  access_one <- hee_colonisation_probability(data.frame(suitability = 0.8),
                                             source_pressure = 0.5,
                                             accessibility = 1,
                                             delta_t = 5)
  expect_equal(no_access$colonisation_probability,
               access_one$colonisation_probability)

  zero_source <- hee_colonisation_probability(data.frame(suitability = 1),
                                              source_pressure = 0,
                                              delta_t = 5)
  expect_lt(zero_source$colonisation_probability, 1e-6)

  unmatched_access <- hee_colonisation_probability(
    data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
               suitability = 0.9),
    source_pressure = data.frame(species = "sp1", cell_id = "c1",
                                 time_ma = 0, source_pressure = 1),
    accessibility = data.frame(species = "sp1", cell_id = "other",
                               time_ma = 0, accessibility = 1)
  )
  expect_equal(unmatched_access$accessibility, 0)
  expect_equal(unmatched_access$colonisation_probability, 0)

  explicit_zero_access <- hee_colonisation_probability(
    data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
               suitability = 1),
    source_pressure = data.frame(species = "sp1", cell_id = "c1",
                                 time_ma = 0, source_pressure = 1),
    accessibility = data.frame(species = "sp1", cell_id = "c1",
                               time_ma = 0, accessibility = 0),
    delta_t = 5
  )
  expect_equal(explicit_zero_access$colonisation_probability, 0)

  explicit_zero_suitability <- hee_colonisation_probability(
    data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
               suitability = 0),
    source_pressure = data.frame(species = "sp1", cell_id = "c1",
                                 time_ma = 0, source_pressure = 1),
    accessibility = data.frame(species = "sp1", cell_id = "c1",
                               time_ma = 0, accessibility = 1),
    delta_t = 5
  )
  expect_equal(explicit_zero_suitability$colonisation_probability, 0)

  occ_tab <- hee_update_occupancy(
    data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
               previous_probability = 1),
    gamma = data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                       colonisation_probability = 0),
    epsilon = data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                         extinction_probability = 0),
    geo_mask = data.frame(cell_id = "other", time_ma = 0,
                          geographic_existence = 1)
  )
  expect_equal(occ_tab$occupancy_probability, 0)
  expect_equal(hee_update_occupancy(1, 0, 0, geo_mask = NA_real_), 0)

  expect_error(
    hee_colonisation_probability(
      data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                 suitability = 0.8),
      source_pressure = data.frame(species = c("sp1", "sp1"),
                                   cell_id = "c1", time_ma = 0,
                                   source_pressure = c(0.2, 0.8))
    ),
    "Duplicate key"
  )

  one_opp <- hee_ecological_opportunity(
    data.frame(cell_id = "c1", time_ma = 0, novelty = 5),
    novelty_col = "novelty",
    weights = c(area = 0, heterogeneity = 0, novelty = 1, young_land = 0)
  )
  expect_gt(one_opp$ecological_opportunity, 0)

  absent_opp <- hee_ecological_opportunity(
    data.frame(cell_id = "c1", time_ma = 0, novelty = 5,
               geographic_existence = 0, habitat_availability = 1),
    novelty_col = "novelty",
    weights = c(area = 0, heterogeneity = 0, novelty = 1, young_land = 0)
  )
  unknown_opp <- hee_ecological_opportunity(
    data.frame(cell_id = "c1", time_ma = 0, novelty = 5,
               geographic_existence = NA_real_, habitat_availability = 1),
    novelty_col = "novelty",
    weights = c(area = 0, heterogeneity = 0, novelty = 1, young_land = 0)
  )
  expect_equal(absent_opp$ecological_opportunity, 0)
  expect_equal(unknown_opp$ecological_opportunity, 0)

  missing_distance_pressure <- hee_colonisation_pressure(
    previous_state = data.frame(species = "sp1", region = "A",
                                time_ma = 0, occupancy_probability = 1),
    connectivity_cube = data.frame(from_region = "A", to_region = "B",
                                   time_ma = 0, connectivity = 1,
                                   paleodistance = NA_real_),
    target = data.frame(species = "sp1", region = "B", time_ma = 0),
    dispersal_scale = c(sp1 = 10)
  )
  expect_equal(missing_distance_pressure$source_pressure, 0)

  base_opp <- hee_ecological_opportunity(
    data.frame(cell_id = "c1", time_ma = 0, land_age_ma = 1),
    weights = c(area = 0, heterogeneity = 0, novelty = 0, young_land = 1)
  )
  expect_warning(
    bad_weight_opp <- hee_ecological_opportunity(
      data.frame(cell_id = "c1", time_ma = 0, land_age_ma = 1),
      weights = c(area = 0, heterogeneity = 0, novelty = 0,
                  young_land = 1, unused_component = 100)
    ),
    "Ignoring unknown ecological-opportunity weight"
  )
  expect_equal(bad_weight_opp$ecological_opportunity,
               base_opp$ecological_opportunity)

  S <- data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                  suitability = 0.8)
  p_d <- hee_combine_static(S, static_filter = data.frame(cell_id = "c1",
                                                          time_ma = 0,
                                                          D_static = 0.25))
  p_w <- hee_combine_static(S, static_filter = data.frame(cell_id = "c1",
                                                          time_ma = 0,
                                                          weight = 0.25))
  expect_equal(p_d$probability, 0.8 * 0.25)
  expect_equal(p_w$probability, 0.8 * 0.25)

  rel_scalar <- hee_prediction_reliability(extrapolation = 0.2,
                                           model_agreement = 0.5)
  rel_vector <- hee_prediction_reliability(extrapolation = c(0.2, 0.4),
                                           model_agreement = c(0.5, 0.5))
  expect_equal(rel_scalar$prediction_reliability, 0.8 * 0.5)
  expect_equal(rel_vector$prediction_reliability, c(0.8, 0.6) * 0.5)

  suit <- data.frame(species = "sp1", cell_id = "c1",
                     time_ma = c(540, 500, 0), suitability = 1)
  prev <- data.frame(species = "sp1", cell_id = "c1", time_ma = 540,
                     occupancy_probability = 1)
  ext <- data.frame(species = "sp1", cell_id = "c1",
                    time_ma = c(540, 500, 0),
                    extinction_probability = c(0, 0.5, 0.25))
  col <- data.frame(species = "sp1", cell_id = "c1",
                    time_ma = c(540, 500, 0),
                    colonisation_probability = 0)
  dyn1 <- hee_dynamic_assembly(prev, suit, colonisation = col, extinction = ext)
  dyn2 <- hee_dynamic_assembly(prev, suit[c(3, 1, 2), ],
                               colonisation = col[c(2, 3, 1), ],
                               extinction = ext[c(3, 1, 2), ])
  dyn1 <- dyn1[order(dyn1$time_ma, decreasing = TRUE), ]
  dyn2 <- dyn2[order(dyn2$time_ma, decreasing = TRUE), ]
  expect_equal(dyn1$occupancy_probability, dyn2$occupancy_probability)
})
