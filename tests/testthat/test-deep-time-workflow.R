make_mock_timecube <- function() {
  e <- new.env(parent = emptyenv())
  e$age_ma <- c(10, 0)
  e$lon <- c(-2, 2)
  e$lat <- c(-1, 1)
  e$variable_names <- c("temp", "precip", "land_mask_dem", "land_area_km2")
  e$product <- list(grid = "mock 4 degree")
  grid <- expand.grid(lon = e$lon, lat = e$lat)
  e$get_slice <- function(age, variables = NULL, as_data_table = TRUE) {
    if (is.null(variables)) variables <- e$variable_names
    d <- data.frame(time_ma = age, lon = grid$lon, lat = grid$lat)
    if ("temp" %in% variables) d$temp <- 10 + grid$lat + age / 10
    if ("precip" %in% variables) d$precip <- 100 + grid$lon * 2 - age
    if ("land_mask_dem" %in% variables) d$land_mask_dem <- c(1, 1, 0, 1)
    if ("land_area_km2" %in% variables) d$land_area_km2 <- 10
    d
  }
  e$get_metadata <- function(variable = NULL) {
    data.frame(variable = e$variable_names, unit = "mock", stringsAsFactors = FALSE)
  }
  class(e) <- "mock_compact_cube"
  e
}

test_that("deep-time timecube, recipe, and projection workflow runs", {
  cube <- hee_load_timecube(object = make_mock_timecube())
  idx <- hee_index_timecube(cube, variables = c("temp", "precip"), times = c(10, 0))
  expect_equal(nrow(idx), 4)
  chk <- hee_check_timecube(cube, variables = c("temp", "precip"), times = c(10, 0))
  expect_true(chk$ok)

  modern <- hee_make_paleo_grid(cube, times = 0, variables = c("temp", "precip"))
  recipe <- hee_lock_recipe(modern, formula = ~ temp + precip)
  beta <- matrix(c(-0.2, 0.5, -0.3,
                   -0.1, -0.4, 0.2),
                 nrow = 2, byrow = TRUE,
                 dimnames = list(c("sp1", "sp2"), c("(Intercept)", "temp", "precip")))
  origin <- c(sp1 = 12, sp2 = 5)
  res <- hee_run_deep_time_pipeline(beta = beta, env_cube = cube, recipe = recipe,
                                    times = c(10, 0), variables = c("temp", "precip"),
                                    species_origin = origin)
  expect_s3_class(res, "hee_deep_time_result")
  expect_true(all(c("suitability", "projection", "richness", "turnover", "refugia") %in% names(res)))
  expect_true(all(res$projection$probability >= 0 & res$projection$probability <= 1))
  expect_true(any(res$projection$species == "sp2" & res$projection$time_ma == 10 &
                    res$projection$phylo_existence == 0))
})

test_that("deep-time track extraction and table projection run", {
  cube <- hee_load_timecube(object = make_mock_timecube())
  modern <- hee_make_paleo_grid(cube, times = 0, variables = c("temp", "precip"))
  recipe <- hee_lock_recipe(modern, formula = ~ temp + precip)
  pts <- data.frame(track_id = "a", time_ma = c(10, 0),
                    paleo_lon = c(-2, 2), paleo_lat = c(-1, 1))
  tr <- hee_extract_paleoenv_track(pts, cube, variables = c("temp", "precip"))
  expect_equal(nrow(tr), 2)
  beta <- matrix(c(0, 0.3, -0.2), nrow = 1,
                 dimnames = list("sp1", c("(Intercept)", "temp", "precip")))
  suit <- hee_project_hmsc_table(beta = beta, X_past = tr, recipe = recipe)
  proj <- hee_combine(suit, phylo_mask = hee_phylo_time_mask(times = c(10, 0),
                                                             species_origin = c(sp1 = 20)))
  rich <- hee_richness(proj, group_cols = c("track_id", "time_ma"))
  ref <- hee_refugia(rich, group_col = "track_id")
  expect_equal(nrow(rich), 2)
  expect_equal(nrow(ref), 1)
})

test_that("paleo grid and track extraction reject missing requested variables", {
  cube <- hee_load_timecube(object = make_mock_timecube())
  expect_error(
    hee_make_paleo_grid(cube, times = 0, variables = c("temp", "missing_var")),
    "missing requested paleo variable"
  )
  pts <- data.frame(track_id = "a", time_ma = 0, paleo_lon = 0, paleo_lat = 0)
  expect_error(
    hee_extract_paleoenv_track(pts, cube, variables = c("temp", "missing_var")),
    "missing requested paleo track variable"
  )
  expect_error(
    hee_make_paleo_grid(cube, times = 0, variables = "temp",
                        landmask = "missing_land", land_only = TRUE),
    "land_only = TRUE requires landmask"
  )
})

test_that("HMSC projection rejects missing Beta axes instead of dropping them", {
  cube <- hee_load_timecube(object = make_mock_timecube())
  beta <- matrix(c(0, 0.3, -0.2), nrow = 1,
                 dimnames = list("sp1", c("(Intercept)", "temp", "missing_axis")))
  expect_error(
    hee_project_hmsc_timecube(beta = beta, env_cube = cube,
                              times = 0, variables = "temp"),
    "missing Beta axis"
  )
  expect_error(
    hee_project_hmsc_table(beta = beta,
                           X_past = data.frame(time_ma = 0, temp = 1)),
    "missing Beta axis"
  )
})

test_that("refugia scores are finite for single-row and constant richness inputs", {
  one <- data.frame(cell_id = "c1", time_ma = 0, expected_richness = 2)
  ref_one <- hee_refugia(one)
  expect_true(is.finite(ref_one$refugia_score))

  constant <- data.frame(cell_id = c("c1", "c2"), time_ma = c(0, 0),
                         expected_richness = c(2, 2))
  ref_constant <- hee_refugia(constant)
  expect_true(all(is.finite(ref_constant$refugia_score)))
})

test_that("deep-time utility checks catch leakage and summarize BGB tables", {
  tax <- hee_make_taxon_table(c("Species one", "Species/two"))
  expect_true(all(c("taxon_id", "file_safe_name") %in% names(tax)))
  roles <- hee_tag_predictors(c("temp", "colonization_age"),
                              c(temp = "exogenous",
                                colonization_age = "composition_derived"))
  expect_error(hee_detect_leakage(~ temp + colonization_age, roles),
               "Information leakage")
  cmp <- hee_bgb_compare(data.frame(model = c("DEC", "DEC+J"),
                                    logLik = c(-10, -8),
                                    n_parameters = c(2, 3)), n = 20)
  expect_true("AICc_weight" %in% names(cmp))
  lrt <- hee_bgb_lrt(cmp[cmp$model == "DEC", ], cmp[cmp$model == "DEC+J", ])
  expect_true(all(c("LR", "p_value") %in% names(lrt)))
})

test_that("all deep-time Word workflow functions are exported", {
  notes <- system.file("extdata", "deep_time_workflow_notes", package = "HmscEcoEvo")
  skip_if(!dir.exists(notes), "deep-time Word notes are not installed")
  files <- list.files(notes, pattern = "\\.txt$", full.names = TRUE)
  txt <- paste(unlist(lapply(files, readLines, warn = FALSE)), collapse = "\n")
  hits <- gregexpr("\\bhee_[A-Za-z][A-Za-z0-9_]*\\b", txt, perl = TRUE)
  funcs <- unique(regmatches(txt, hits)[[1]])
  funcs <- funcs[vapply(funcs, function(z) length(gregexpr("hee_", z, fixed = TRUE)[[1]]) == 1,
                        logical(1))]
  missing <- setdiff(funcs, getNamespaceExports("HmscEcoEvo"))
  expect_equal(missing, character())
  expect_gte(length(funcs), 100)
})

test_that("deep-time extension helpers compute lightweight outputs", {
  sim <- hee_simulate_example_data(n_species = 4, n_site = 8,
                                   variables = c("temp", "precip"))
  dat <- hee_prepare_data(sim$comm, sim$env_now)
  audit <- hee_audit_data(dat)
  expect_true(all(c("component", "ok", "n_warning") %in% names(audit)))

  hmsc_input <- hee_as_hmsc_input(dat, XFormula = ~ temp + precip)
  expect_equal(dim(hmsc_input$Y), c(8, 4))
  expect_equal(names(hmsc_input$XData), c("temp", "precip"))

  expect_equal(hee_cv_random(10, k = 3, seed = 1), hee_cv_random(10, k = 3, seed = 1))
  expect_equal(length(hee_cv_spatial_block(data.frame(lon = seq(-10, 10, length.out = 8),
                                                     lat = 0), k = 4)), 8)
  expect_equal(length(hee_cv_environment_block(sim$env_now, variable = "temp", k = 4)), 8)

  filtered <- hee_filter_collinearity(data.frame(a = 1:10, b = 1:10, c = c(1:5, 5:1)))
  expect_true("b" %in% filtered$drop)
  imputed <- hee_impute_traits(data.frame(trait1 = c(1, NA, 3),
                                          trait2 = c("a", NA, "a")))
  expect_false(anyNA(imputed$traits))
  expect_true(any(imputed$missing_flags$trait1_missing))

  metrics <- hee_evaluate_hmsc(c(0, 0, 1, 1), c(0.1, 0.2, 0.8, 0.9))
  expect_true(all(c("AUC", "Tjur_R2", "RMSE", "calibration_slope") %in% names(metrics)))
  expect_equal(metrics$AUC, 1)
})

test_that("deep-time extension helpers handle external interfaces honestly", {
  expect_error(hee_bgb_run_models(), "not run internally")
  expect_error(hee_reconstruct_points(data.frame(lon = 0, lat = 0)),
               "Plate reconstruction")
  expect_error(
    hee_import_plate_points(data.frame(lon = 0, lat = 0, time_ma = 100)),
    "require plate-reconstructed"
  )
  modern_pt <- hee_import_plate_points(data.frame(lon = 0, lat = 0, time_ma = 0))
  expect_equal(modern_pt$paleo_lon, modern_pt$lon)
  expect_equal(modern_pt$paleo_lat, modern_pt$lat)
  expect_equal(modern_pt$coord_status, "modern_coordinates_0Ma")
  expect_error(
    hee_compare_plate_corrected(
      data.frame(point_id = c("p1", "p1"), lon = c(0, 1), lat = 0),
      data.frame(point_id = "p1", paleo_lon = 1, paleo_lat = 1)
    ),
    "Duplicate key"
  )

  acc <- hee_bgb_accessibility(list(
    DEC = data.frame(species = "sp1", region = "A", time_ma = 0, accessibility = 0.2),
    DECJ = data.frame(species = "sp1", region = "A", time_ma = 0, accessibility = 0.8)
  ), weights = c(DEC = 0.25, DECJ = 0.75))
  expect_equal(acc$accessibility, 0.65)

  ev <- hee_bgb_event_summary(data.frame(
    event_type = "dispersal", from_region = "A", to_region = "B",
    region = "B", time_ma = 10, event_prob = 0.5
  ))
  expect_equal(ev$dispersal$dispersal_rate, 0.5)

  project <- tempfile("hee_project_")
  expect_false(dir.exists(project))
  hee_init_project(project)
  expect_true(dir.exists(file.path(project, "config")))
})
