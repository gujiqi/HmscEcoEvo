test_that("combine multiplies HMSC suitability by masks and filters", {
  suit <- data.frame(species = c("sp1", "sp1"), cell_id = c(1, 1),
                     time_ma = c(50, 0), suitability = c(0.8, 0.8))
  phylo <- hee_phylo_time_mask(times = c(50, 0), species_origin = c(sp1 = 10))
  land <- data.frame(cell_id = 1, time_ma = c(50, 0), land_weight = 1)
  dyn <- data.frame(cell_id = 1, time_ma = c(50, 0), dynamic_weight = c(0.5, 1))
  out <- hee_combine(suit, phylo_mask = phylo, landmask = land, dynamic_filter = dyn)
  expect_equal(out$probability[out$time_ma == 50], 0)
  expect_equal(out$probability[out$time_ma == 0], 0.8)
})

test_that("combine sets probability to zero on non-land cells", {
  suit <- data.frame(species = "sp1", cell_id = 1, time_ma = 0, suitability = 0.9)
  land <- data.frame(cell_id = 1, time_ma = 0, land_weight = 0)
  out <- hee_combine(suit, landmask = land)
  expect_equal(out$probability, 0)
})

test_that("combine hard-zero masks override missing suitability", {
  suit <- data.frame(species = "sp1", cell_id = "c1", time_ma = 50,
                     suitability = NA_real_)
  phylo <- data.frame(species = "sp1", time_ma = 50, phylo_existence = 0)
  out_phylo <- hee_combine(suit, phylo_mask = phylo)
  expect_equal(out_phylo$phylo_existence, 0)
  expect_equal(out_phylo$probability, 0)

  land <- data.frame(cell_id = "c1", time_ma = 50, geographic_existence = 0)
  out_land <- hee_combine(suit, landmask = land)
  expect_equal(out_land$land_weight, 0)
  expect_equal(out_land$probability, 0)
})

test_that("combine uses explicit land and static process columns, not arbitrary numeric columns", {
  suit <- data.frame(species = "sp1", cell_id = "c1", time_ma = 0, suitability = 0.9)
  land <- data.frame(cell_id = "c1", time_ma = 0, temp = 12, land_mask_dem = 0)
  out <- hee_combine(suit, landmask = land)
  expect_equal(out$land_weight, 0)
  expect_equal(out$probability, 0)

  static <- data.frame(
    species = "sp1",
    cell_id = "c1",
    time_ma = 0,
    structural_connectivity = 0.5,
    functional_connectivity = 0.2,
    climatic_connectivity = 0.5
  )
  out2 <- hee_combine(suit, static_filter = static)
  expect_equal(out2$D_static, 0.5 * 0.2 * 0.5)
  expect_equal(out2$probability, 0.9 * 0.5 * 0.2 * 0.5)
  expect_true(all(c("structural_connectivity", "functional_connectivity",
                    "climatic_connectivity") %in% names(out2)))
})

test_that("combine rejects duplicated filter keys instead of silently expanding rows", {
  suit <- data.frame(species = "sp1", cell_id = "c1", time_ma = 0, suitability = 0.8)
  land <- data.frame(cell_id = "c1", time_ma = 0, land_weight = c(1, 0))
  expect_error(hee_combine(suit, landmask = land), "Duplicate key")
})

test_that("supplied combine components must expose recognised value columns", {
  suit <- data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                     suitability = 0.8)
  bad_cell_time <- data.frame(cell_id = "c1", time_ma = 0,
                              bad_weight = 0.5)
  bad_species_cell_time <- data.frame(species = "sp1", cell_id = "c1",
                                      time_ma = 0, bad_weight = 0.5)
  expect_error(hee_combine(suit, accessibility = bad_species_cell_time),
               "accessibility was supplied")
  expect_error(hee_combine(suit, landmask = bad_cell_time),
               "landmask was supplied")
  expect_error(hee_combine(suit, static_filter = bad_species_cell_time),
               "static_filter was supplied")
  expect_error(hee_combine(suit, dynamic_filter = bad_species_cell_time),
               "dynamic_filter was supplied")
  expect_error(hee_combine(suit, extrapolation_weight = bad_cell_time),
               "extrapolation_weight was supplied")
})

test_that("combine applies named matrix accessibility by species and time", {
  suit <- data.frame(
    species = c("sp1", "sp2"),
    cell_id = c("c1", "c1"),
    time_ma = c(0, 0),
    suitability = c(1, 1)
  )
  acc <- matrix(c(0.2, 0.8), nrow = 2,
                dimnames = list(c("sp1", "sp2"), "0Ma"))
  out <- hee_combine(suit, accessibility = acc)
  expect_equal(out$accessibility, c(0.2, 0.8))
  expect_equal(out$probability, c(0.2, 0.8))
})

test_that("species-level accessibility and phylo masks are not silently broadcast", {
  suit <- data.frame(
    species = c("sp1", "sp2"),
    cell_id = c("c1", "c2"),
    region = "R1",
    time_ma = 0,
    suitability = 1
  )
  acc <- data.frame(region = "R1", time_ma = 0, accessibility = 0.2)
  expect_error(hee_combine(suit, accessibility = acc),
               "without a `species` column")
  shared <- hee_combine(suit, accessibility = acc,
                        allow_species_broadcast = TRUE)
  expect_equal(shared$probability, c(0.2, 0.2))

  phylo <- data.frame(time_ma = 0, phylo_existence = 1)
  expect_error(hee_combine(suit, phylo_mask = phylo),
               "without a `species` column")
})

test_that("combine rejects unnamed matrix weights instead of treating them as neutral", {
  suit <- data.frame(species = "sp1", cell_id = "c1", time_ma = 0, suitability = 1)
  acc <- matrix(0.2, nrow = 1)
  expect_error(hee_combine(suit, accessibility = acc), "must have row and column names")
})

test_that("supplied filters with missing keys do not default to full accessibility", {
  suit <- data.frame(
    species = "sp1",
    cell_id = c("c1", "c2"),
    time_ma = 0,
    suitability = 1
  )
  acc <- data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                    accessibility = 0.5)
  land <- data.frame(cell_id = "c1", time_ma = 0, land_weight = 1)
  static <- data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                       D_static = 0.8)
  dyn <- data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                    D_dynamic = 0.7)
  out <- hee_combine(suit, accessibility = acc, landmask = land,
                     static_filter = static, dynamic_filter = dyn)
  expect_equal(out$accessibility[out$cell_id == "c2"], 0)
  expect_equal(out$land_weight[out$cell_id == "c2"], 0)
  expect_equal(out$D_static[out$cell_id == "c2"], 0)
  expect_equal(out$D_dynamic[out$cell_id == "c2"], 0)
  expect_equal(out$probability[out$cell_id == "c2"], 0)
})

test_that("supplied extrapolation weights with missing keys are conservative", {
  suit <- data.frame(
    species = "sp1",
    cell_id = c("c1", "c2"),
    time_ma = 0,
    suitability = 1
  )
  q <- data.frame(cell_id = "c1", time_ma = 0, Q_extrap = 0)
  out <- hee_combine(suit, extrapolation_weight = q)
  expect_equal(out$Q_extrap[out$cell_id == "c1"], 0)
  expect_equal(out$Q_extrap[out$cell_id == "c2"], 0)
  expect_equal(out$probability, c(0, 0))
})

test_that("explicit NA land or static weights are conservative, not full access", {
  suit <- data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                     suitability = 1)
  land <- data.frame(cell_id = "c1", time_ma = 0, land_weight = NA_real_)
  static <- data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                       D_static = NA_real_)
  out_land <- hee_combine(suit, landmask = land)
  out_static <- hee_combine(suit, static_filter = static)
  expect_equal(out_land$land_weight, 0)
  expect_equal(out_land$probability, 0)
  expect_equal(out_static$D_static, 0)
  expect_equal(out_static$probability, 0)

  local_land <- suit
  local_land$geographic_existence <- NA_real_
  out_local <- hee_combine(local_land)
  expect_equal(out_local$land_weight, 0)
  expect_equal(out_local$probability, 0)

  static_components <- data.frame(species = "sp1", cell_id = "c1",
                                  time_ma = 0,
                                  structural_connectivity = NA_real_,
                                  functional_connectivity = 1,
                                  climatic_connectivity = 1)
  out_components <- hee_combine(suit, static_filter = static_components)
  expect_equal(out_components$structural_connectivity, 0)
  expect_equal(out_components$D_static, 0)
  expect_equal(out_components$probability, 0)

  static_omitted_component <- data.frame(species = "sp1", cell_id = "c1",
                                         time_ma = 0,
                                         functional_connectivity = 0.5)
  out_omitted <- hee_combine(suit, static_filter = static_omitted_component)
  expect_equal(out_omitted$structural_connectivity, 1)
  expect_equal(out_omitted$D_static, 0.5)
})

test_that("supplied phylogenetic masks treat unmatched species or times as absent", {
  suit <- data.frame(
    species = c("sp1", "sp2", "sp1"),
    cell_id = "c1",
    time_ma = c(0, 0, 10),
    suitability = 1
  )
  mask <- matrix(1, nrow = 1, dimnames = list("sp1", "0Ma"))
  out <- hee_combine(suit, phylo_mask = mask)
  expect_equal(out$phylo_existence[out$species == "sp2"], 0)
  expect_equal(out$phylo_existence[out$time_ma == 10], 0)
  expect_equal(out$probability[out$species == "sp2"], 0)
  expect_equal(out$probability[out$time_ma == 10], 0)
})

test_that("combine supports unique reconstructed tracks without a cell_id", {
  suitability <- data.frame(
    species = "sp1", track_id = "track_a", time_ma = c(10, 0),
    suitability = c(0.4, 0.8), stringsAsFactors = FALSE
  )
  out <- hee_combine(suitability)
  expect_equal(out$track_id, c("track_a", "track_a"))
  expect_equal(out$probability, c(0.4, 0.8))
  duplicated <- rbind(suitability[1, ], suitability[1, ])
  expect_error(hee_combine(duplicated), "Duplicate key")
})

test_that("combine accepts Dynamic Earth-Biota framework aliases", {
  suit <- data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                     S_HMSC = 0.9)
  acc <- data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                    A_hist = 0.5)
  phylo <- data.frame(species = "sp1", time_ma = 0, E_phylo = 1)
  land <- data.frame(cell_id = "c1", time_ma = 0, G_arena = 0.8)
  dyn <- data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                    D_dynamic = 0.25)
  q <- data.frame(cell_id = "c1", time_ma = 0, Q_extrap = 0.5)
  out <- hee_combine(suit, accessibility = acc, phylo_mask = phylo,
                     landmask = land, dynamic_filter = dyn,
                     extrapolation_weight = q)
  expect_equal(out$probability, 0.9 * 0.5 * 1 * 0.8 * 0.25 * 0.5)
  expect_equal(out$S_HMSC, 0.9)
  expect_equal(out$A_hist, 0.5)
  expect_equal(out$E_phylo, 1)
  expect_equal(out$L_land, 0.8)
})

test_that("missing HMSC suitability is conservative and auditable", {
  suit <- data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                     S_HMSC = NA_real_)
  out <- hee_combine(suit)
  expect_true(out$missing_suitability)
  expect_equal(out$probability, 0)
})
