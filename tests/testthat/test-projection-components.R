test_that("projection formulas and numeric example are transparent", {
  defs <- hee_projection_formula(extrapolation_mode = "flag")
  expect_equal(nrow(defs), 5)
  expect_true(any(grepl("D_dynamic", defs$formula)))

  ex <- hee_projection_component_example(
    S_HMSC = 0.8, A_BGB = 0.25, E_phylo = 1,
    L_land = 1, D_static = 0.6, D_dynamic = 0.1
  )
  expect_equal(ex$probability[ex$model_id == "M1_env_only"], 0.8)
  expect_equal(ex$probability[ex$model_id == "M5_env_phylo_bgb_dynamic_dispersal"], 0.02)
})

test_that("M1-M5 wrapper refuses complete labels when required components are missing", {
  suit <- data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                     suitability = 0.8)
  expect_error(
    hee_combine_projection_models(suit),
    "Strict M1-M5 projection requires all named process components"
  )

  land <- data.frame(cell_id = "c1", time_ma = 0, land_weight = 1)
  phylo <- data.frame(species = "sp1", time_ma = 0, phylo_existence = 1)
  acc <- data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                    accessibility = 0.5)
  static <- data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                       D_static = 0.5)
  dynamic <- data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
                        D_dynamic = 0.25)
  out <- hee_combine_projection_models(
    suit,
    accessibility = acc,
    phylo_mask = phylo,
    landmask = land,
    static_filter = static,
    dynamic_filter = dynamic
  )
  expect_equal(sort(unique(out$model_id)), sort(hee_projection_formula()$model_id))
  expect_equal(out$probability[out$model_id == "M5_env_phylo_bgb_dynamic_dispersal"],
               0.8 * 0.5 * 1 * 1 * 0.25)

  neutral_demo <- hee_combine_projection_models(suit, strict = FALSE)
  expect_equal(length(unique(neutral_demo$model_id)), 5)
})

test_that("extrapolation is a flag by default and a downweight only on request", {
  x <- data.frame(cell_id = 1:2, time_ma = 0, extrapolation_score = c(0, 2))
  flag <- hee_extrapolation_weight(x, mode = "flag")
  down <- hee_extrapolation_weight(x, mode = "downweight", lambda = 1)
  expect_equal(flag$Q_extrap, c(1, 1))
  expect_lt(down$Q_extrap[2], 1)

  missing_score <- hee_extrapolation_weight(
    data.frame(cell_id = 1, time_ma = 0, extrapolation_score = NA_real_),
    mode = "downweight"
  )
  expect_equal(missing_score$Q_extrap, 1)

  infinite_score <- hee_extrapolation_weight(
    data.frame(cell_id = 1, time_ma = 0, extrapolation_score = Inf),
    mode = "downweight"
  )
  expect_equal(infinite_score$Q_extrap, 0)

  negative_score <- hee_extrapolation_weight(
    data.frame(cell_id = 1, time_ma = 0, extrapolation_score = -Inf),
    mode = "downweight"
  )
  expect_equal(negative_score$Q_extrap, 1)
})

test_that("static dispersal decays with distance from accessible source", {
  cells <- data.frame(cell_id = 1:3, time_ma = 0, lon = c(0, 1, 4),
                      lat = 0, region = c("A", "A", "B"))
  acc <- data.frame(species = "sp1", region = c("A", "B"), time_ma = 0,
                    accessibility = c(1, 0))
  d <- hee_static_dispersal(cells, acc, dispersal_scale = c(sp1 = 1))
  expect_true(all(c("static_weight", "D_static") %in% names(d)))
  expect_gte(d$D_static[d$cell_id == 1], d$D_static[d$cell_id == 3])
})

test_that("coordinate dispersal distances are great-circle kilometres, not degrees", {
  cells <- data.frame(species = "sp1", cell_id = c("near", "far"),
                      time_ma = 0, lon = c(1, 2), lat = 0)
  src <- data.frame(species = "sp1", cell_id = "source", time_ma = 0,
                    lon = 0, lat = 0)
  out <- hee_static_dispersal(cells, source_cells = src,
                              dispersal_scale = 111.195)
  expect_equal(out$distance_to_accessible_source[out$cell_id == "near"],
               hee_great_circle_distance_km(1, 0, 0, 0),
               tolerance = 0.01)
  expect_equal(out$D_static[out$cell_id == "near"], exp(-1),
               tolerance = 0.01)
})

test_that("static dispersal rejects duplicate accessibility and barrier keys", {
  cells <- data.frame(cell_id = 1:2, time_ma = 0, lon = c(0, 1),
                      lat = 0, region = "A")
  dup_acc <- data.frame(species = "sp1", region = "A", time_ma = 0,
                        accessibility = c(1, 0.5))
  expect_error(hee_static_dispersal(cells, dup_acc), "Duplicate key")
  acc <- data.frame(species = "sp1", region = "A", time_ma = 0,
                    accessibility = 1)
  dup_barrier <- data.frame(cell_id = 1, time_ma = 0,
                            barrier_passability = c(1, 0.5))
  expect_error(hee_static_dispersal(cells, acc, barrier = dup_barrier),
               "Duplicate key")
  global_barrier <- data.frame(barrier_passability = 0.5)
  out <- hee_static_dispersal(cells, acc, barrier = global_barrier)
  expect_equal(out$barrier_passability, c(0.5, 0.5))

  sparse_barrier <- data.frame(cell_id = 1, time_ma = 0,
                               barrier_passability = NA_real_)
  sparse <- hee_static_dispersal(cells, acc, barrier = sparse_barrier)
  expect_equal(sparse$barrier_passability, c(0, 0))
  expect_equal(sparse$D_static, c(0, 0))
  expect_error(hee_static_dispersal(cells, acc, dispersal_scale = -1),
               "finite non-negative")
  expect_error(hee_make_dispersal_matrices(c("A", "B"), scale = -1),
               "finite non-negative")
  W0 <- hee_make_dispersal_matrices(
    c("A", "B"),
    distance = matrix(c(0, 10, NA, 0), nrow = 2,
                      dimnames = list(c("A", "B"), c("A", "B"))),
    scale = 0
  )
  expect_equal(unname(diag(W0)), c(1, 1))
  expect_equal(W0[1, 2], 0)
  expect_equal(W0[2, 1], 0)
  expect_error(
    hee_make_dispersal_matrices(
      c("A", "B"),
      distance = matrix(c(0, -1, 1, 0), nrow = 2),
      scale = 1
    ),
    "non-negative"
  )
})

test_that("dynamic filter uses previous time slice under forward direction", {
  prev <- data.frame(species = "sp1", cell_id = 1, time_ma = 10,
                     lon = 0, lat = 0, probability = 1)
  target <- data.frame(species = "sp1", cell_id = c(1, 2), time_ma = 0,
                       lon = c(0, 5), lat = 0)
  d <- hee_dynamic_filter(prev, cell_table = target, dispersal_scale = c(sp1 = 1),
                          time_direction = "forward", default = 0)
  expect_equal(d$D_dynamic[d$cell_id == 1], 1)
  expect_lt(d$D_dynamic[d$cell_id == 2], d$D_dynamic[d$cell_id == 1])
  d_barrier <- hee_dynamic_filter(
    prev,
    cell_table = target,
    dispersal_scale = c(sp1 = 1),
    time_direction = "forward",
    default = 0,
    barrier = data.frame(cell_id = 1, time_ma = 0,
                         barrier_passability = NA_real_)
  )
  expect_equal(d_barrier$barrier_passability, c(0, 0))
  expect_equal(d_barrier$D_dynamic, c(0, 0))
  expect_error(hee_dynamic_filter(prev, cell_table = target,
                                  dispersal_scale = -1),
               "finite non-negative")
})

test_that("dynamic filter reports minimum distance in kilometres", {
  prev <- data.frame(species = "sp1", cell_id = "src", time_ma = 10,
                     lon = 0, lat = 0, probability = 1)
  target <- data.frame(species = "sp1", cell_id = "dst", time_ma = 0,
                       lon = 1, lat = 0)
  d <- hee_dynamic_filter(prev, cell_table = target,
                          dispersal_scale = 111.195,
                          time_direction = "forward", default = 0)
  expect_equal(d$dynamic_min_distance,
               hee_great_circle_distance_km(1, 0, 0, 0),
               tolerance = 0.01)
  expect_equal(d$D_dynamic, exp(-1), tolerance = 0.01)
})

test_that("region-level dispersal filters preserve rows and use prior slices", {
  cells <- expand.grid(species = "sp1", time_ma = c(20, 10, 0),
                       region = c("A", "B"), cell_id = 1:2,
                       KEEP.OUT.ATTRS = FALSE)
  cells$lon <- ifelse(cells$region == "A", 0, 10)
  cells$lat <- 0
  acc <- data.frame(species = "sp1", region = "A", time_ma = c(20, 10, 0),
                    accessibility = 1)

  s <- hee_static_region_dispersal(cells, acc, dispersal_scale = 5)
  expect_equal(nrow(s), nrow(cells))
  expect_true(all(c("static_weight", "D_static") %in% names(s)))
  expect_equal(s$D_static[s$region == "A"][1], 1)
  expect_lt(s$D_static[s$region == "B"][1], 1)

  acc_zero <- acc
  acc_zero$accessibility <- 0
  s_zero <- hee_static_region_dispersal(cells, acc_zero, dispersal_scale = 5)
  expect_equal(s_zero$D_static, rep(0, nrow(cells)))

  cell_zero <- hee_static_dispersal(
    cells[c("species", "cell_id", "time_ma", "lon", "lat", "region")],
    acc_zero,
    dispersal_scale = 5
  )
  expect_equal(cell_zero$D_static, rep(0, nrow(cells)))

  expect_error(hee_static_region_dispersal(cells, acc,
                                           dispersal_scale = -5),
               "finite non-negative")

  prev <- cells
  prev$probability <- ifelse(prev$time_ma == 20 & prev$region == "A", 0.8, 0.1)
  d <- hee_dynamic_region_filter(prev, cell_table = cells, dispersal_scale = 5,
                                 time_direction = "forward", default = 0.25)
  expect_equal(nrow(d), nrow(cells))
  expect_equal(d$D_dynamic[d$time_ma == 20][1], 0.25)
  expect_gt(max(d$D_dynamic[d$time_ma == 10 & d$region == "A"]), 0.7)
  expect_error(hee_dynamic_region_filter(prev, cell_table = cells,
                                         dispersal_scale = -5),
               "finite non-negative")
})

test_that("dynamic dispersal filters do not use first-slice default for missing source evidence", {
  prev <- data.frame(species = "sp1", region = "A", time_ma = 20,
                     probability = 1, lon = 0, lat = 0)
  cells <- expand.grid(species = c("sp1", "sp2"), time_ma = c(20, 10),
                       region = c("A", "B"), KEEP.OUT.ATTRS = FALSE)
  cells$lon <- ifelse(cells$region == "A", 0, 10)
  cells$lat <- 0
  region_dyn <- hee_dynamic_region_filter(prev, cell_table = cells,
                                          dispersal_scale = 5,
                                          time_direction = "forward",
                                          default = 0.75)
  expect_equal(unique(region_dyn$D_dynamic[region_dyn$species == "sp2" &
                                             region_dyn$time_ma == 20]), 0.75)
  expect_equal(unique(region_dyn$D_dynamic[region_dyn$species == "sp2" &
                                             region_dyn$time_ma == 10]), 0)

  prev_point <- data.frame(species = "sp1", cell_id = "c0", region = "A",
                           time_ma = 20, probability = 1, lon = 0, lat = 0)
  point_dyn <- hee_dynamic_filter(prev_point, cell_table = data.frame(
    species = c("sp1", "sp2"), cell_id = c("c1", "c2"),
    time_ma = 10, lon = c(0, 10), lat = 0
  ), dispersal_scale = 5, time_direction = "forward", default = 0.75)
  expect_gt(point_dyn$D_dynamic[point_dyn$species == "sp1"], 0)
  expect_equal(point_dyn$D_dynamic[point_dyn$species == "sp2"], 0)

  expect_error(
    hee_dynamic_filter(prev_point[, c("species", "cell_id", "time_ma", "probability")],
                       cell_table = data.frame(species = "sp1", cell_id = "c1",
                                               time_ma = 10, lon = 0, lat = 0),
                       time_direction = "forward"),
    "lon/lat"
  )
})

test_that("region static dispersal rejects duplicate accessibility keys", {
  cells <- expand.grid(species = "sp1", time_ma = 0, region = c("A", "B"),
                       cell_id = 1:2, KEEP.OUT.ATTRS = FALSE)
  cells$lon <- ifelse(cells$region == "A", 0, 10)
  cells$lat <- 0
  dup_acc <- data.frame(species = "sp1", region = "A", time_ma = 0,
                        accessibility = c(1, 0.8))
  expect_error(hee_static_region_dispersal(cells, dup_acc), "Duplicate key")
})

test_that("hee_combine exposes static, dynamic, and Q components", {
  suit <- data.frame(species = "sp1", cell_id = 1, time_ma = 0,
                     suitability = 0.8)
  static <- data.frame(species = "sp1", cell_id = 1, time_ma = 0,
                       static_weight = 0.5)
  dyn <- data.frame(species = "sp1", cell_id = 1, time_ma = 0,
                    dynamic_weight = 0.25)
  q <- data.frame(cell_id = 1, time_ma = 0, Q_extrap = 0.5)
  p <- hee_combine(suit, static_filter = static, dynamic_filter = dyn,
                   extrapolation_weight = q)
  expect_equal(p$probability, 0.8 * 0.5 * 0.25 * 0.5)
  expect_true(all(c("D_static", "D_dynamic", "Q_extrap") %in% names(p)))
})
