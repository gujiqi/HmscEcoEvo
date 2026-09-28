#!/usr/bin/env Rscript

args <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)][1])
script_dir <- if (!is.na(file_arg) && nzchar(file_arg)) dirname(normalizePath(file_arg)) else getwd()
pkg_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = TRUE)

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(pkg_root, quiet = TRUE)
} else {
  library(HmscEcoEvo)
}
if (!requireNamespace("ape", quietly = TRUE)) {
  stop("Package 'ape' is required to generate dated toy phylogenies.", call. = FALSE)
}

set.seed(20260608)

root <- script_dir
source_rds <- "C:/Users/Google/Documents/HMSC-HIST/Phanerozoic_Environment_540_0Ma_5Myr_4deg_landharmonized_v2/08_rds_export/Phanerozoic_Environment_540_0Ma_5Myr_Level2_4deg_landharmonized_v2_compact_float32.rds"
if (!file.exists(source_rds)) {
  stop("Cannot find 4-degree Phanerozoic RDS: ", source_rds, call. = FALSE)
}

old_case_dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
old_case_dirs <- old_case_dirs[grepl("^case0[12]_", basename(old_case_dirs))]
if (length(old_case_dirs) > 0) unlink(old_case_dirs, recursive = TRUE, force = TRUE)
old_outputs <- file.path(root, "outputs")
if (dir.exists(old_outputs)) unlink(old_outputs, recursive = TRUE, force = TRUE)

cube <- hee_load_timecube(source_rds)
times <- sort(unique(cube$times), decreasing = TRUE)
if (!all(seq(540, 0, by = -5) %in% times)) {
  stop("The selected RDS does not contain the full 540-0 Ma / 5 Myr time axis.",
       call. = FALSE)
}
times <- seq(540, 0, by = -5)

variables <- c(
  "MAT_pohl_C",
  "MAP_pohl_mm_yr",
  "P_seasonality_pohl_sd",
  "elevation_m",
  "distance_to_coast_km",
  "moisture_availability_index_z",
  "wetland_potential_index_z",
  "aridity_index_P_over_E"
)
missing_vars <- setdiff(variables, cube$variables)
if (length(missing_vars) > 0) {
  stop("The source RDS is missing variable(s): ", paste(missing_vars, collapse = ", "),
       call. = FALSE)
}

species <- sprintf("sp%02d", 1:20)
clade <- setNames(rep(LETTERS[1:5], each = 4), species)
region_names <- paste0("Region_", LETTERS[1:6])

make_time_tree <- function(case_id = 1) {
  root_age <- 540
  crowns <- if (case_id == 1) {
    c(A = 520, B = 420, C = 315, D = 210, E = 115)
  } else {
    c(A = 500, B = 385, C = 285, D = 185, E = 95)
  }
  split_a <- if (case_id == 1) {
    c(A = 390, B = 285, C = 205, D = 125, E = 45)
  } else {
    c(A = 430, B = 255, C = 165, D = 95, E = 35)
  }
  split_b <- if (case_id == 1) {
    c(A = 470, B = 350, C = 250, D = 160, E = 75)
  } else {
    c(A = 455, B = 315, C = 225, D = 135, E = 60)
  }
  clade_tips <- split(species, clade)
  subtrees <- vapply(names(clade_tips), function(g) {
    sp <- clade_tips[[g]]
    crown <- crowns[g]
    s1 <- min(split_a[g], crown - 5)
    s2 <- min(split_b[g], crown - 5)
    sprintf("((%s:%.3f,%s:%.3f):%.3f,(%s:%.3f,%s:%.3f):%.3f):%.3f",
            sp[1], s1, sp[2], s1, crown - s1,
            sp[3], s2, sp[4], s2, crown - s2,
            root_age - crown)
  }, character(1))
  tr <- ape::read.tree(text = paste0("(", paste(subtrees, collapse = ","), ");"))
  tr$root.time <- root_age
  tr
}

make_traits <- function(case_id = 1) {
  idx <- seq_along(species)
  cl <- as.numeric(factor(clade[species]))
  x <- data.frame(
    body_size = scale(log(idx + 3) + 0.35 * cl + sin(idx / 3)),
    dispersal = scale(runif(length(species), 0.2, 1.4) + 0.10 * cl),
    thermal_tolerance = scale(cos(idx / 2.6) - 0.25 * cl),
    moisture_affinity = scale(sin(idx / 4.2) + 0.22 * cl),
    coastal_affinity = scale(rep(c(1.2, 0.6, -0.3, -0.8), 5) + rnorm(length(species), 0, 0.12)),
    elevation_tolerance = scale(rep(c(-0.8, -0.2, 0.5, 1.0), 5) + 0.08 * cl),
    row.names = species
  )
  as.data.frame(x) + case_id * 0.015
}

make_beta <- function(traits, case_id = 1) {
  coef_mat <- matrix(c(
     0.42, -0.18,  0.10, -0.20,  0.05, -0.08,  0.14, -0.16,
    -0.14,  0.28, -0.08,  0.02,  0.22,  0.18, -0.10,  0.20,
     0.30, -0.10,  0.16, -0.06, -0.08,  0.22,  0.04, -0.20,
    -0.08,  0.22,  0.18, -0.18,  0.04,  0.32,  0.28, -0.26,
     0.10, -0.06, -0.02,  0.08,  0.36, -0.10,  0.24, -0.08,
    -0.18,  0.02,  0.10,  0.40, -0.12, -0.16, -0.06,  0.18
  ), nrow = ncol(traits), byrow = TRUE,
  dimnames = list(colnames(traits), variables))
  clade_shift <- matrix(c(
     0.25, -0.20,  0.10, -0.08,  0.10, -0.10,  0.06, -0.12,
    -0.10,  0.28, -0.06,  0.06, -0.12,  0.18,  0.08, -0.04,
     0.12, -0.06,  0.26,  0.10, -0.02,  0.04,  0.18, -0.10,
    -0.18,  0.08, -0.10,  0.30,  0.06, -0.16, -0.08,  0.22,
     0.08, -0.14,  0.04, -0.10,  0.24,  0.12,  0.26, -0.18
  ), nrow = 5, byrow = TRUE, dimnames = list(LETTERS[1:5], variables))
  B <- as.matrix(traits) %*% coef_mat + clade_shift[clade[species], , drop = FALSE]
  B <- B / max(abs(B), na.rm = TRUE) * (0.85 + 0.03 * case_id)
  intercept <- -0.35 + 0.10 * as.numeric(traits$body_size) -
    0.08 * as.numeric(traits$dispersal) + rnorm(length(species), 0, 0.04)
  cbind("(Intercept)" = intercept, B)
}

make_tip_ranges <- function(case_id = 1) {
  mat <- matrix(0, nrow = length(species), ncol = length(region_names),
                dimnames = list(species, region_names))
  for (i in seq_along(species)) {
    primary <- ((i + case_id - 1) %% length(region_names)) + 1
    secondary <- ((primary + as.numeric(factor(clade[species[i]])) - 1) %% length(region_names)) + 1
    mat[i, unique(c(primary, secondary))] <- 1
  }
  mat
}

make_bgb_objects <- function(tip_ranges, species_origin) {
  models <- data.frame(
    model = c("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J"),
    logLik = c(-124.4, -119.8, -127.0, -121.2, -132.1, -126.4),
    n_parameters = c(2, 3, 2, 3, 2, 3)
  )
  compare <- hee_bgb_compare(models, n = nrow(tip_ranges) * ncol(tip_ranges))
  event_times <- seq(540, 0, by = -20)
  events <- data.frame(
    event_type = rep(c("dispersal", "extinction", "range expansion"), length.out = length(event_times)),
    from_region = rep(region_names, length.out = length(event_times)),
    to_region = rep(rev(region_names), length.out = length(event_times)),
    region = rep(region_names, length.out = length(event_times)),
    time_ma = event_times,
    event_prob = pmin(0.95, pmax(0.05, 0.45 + 0.35 * sin(event_times / 55))),
    stringsAsFactors = FALSE
  )
  accessibility_models <- lapply(compare$model, function(model_id) {
    weight_shift <- compare$AICc_weight[compare$model == model_id]
    rows <- expand.grid(species = rownames(tip_ranges), region = region_names,
                        time_ma = seq(540, 0, by = -5),
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    origin <- species_origin[rows$species]
    present <- as.numeric(rows$time_ma <= origin)
    region_pref <- as.numeric(tip_ranges[cbind(rows$species, rows$region)])
    rows$accessibility <- present * pmin(1, 0.20 + 0.55 * region_pref +
                                           0.20 * exp(-rows$time_ma / 220) +
                                           0.10 * weight_shift)
    rows
  })
  names(accessibility_models) <- compare$model
  accessibility <- hee_bgb_accessibility(accessibility_models,
                                         weights = stats::setNames(compare$AICc_weight, compare$model))
  list(
    prepared = hee_bgb_prepare(tip_ranges, regions = region_names),
    models = models,
    compare = compare,
    lrt_DEC_J = hee_bgb_lrt(compare[compare$model == "DEC", ],
                            compare[compare$model == "DEC+J", ]),
    events = events,
    event_summary = hee_bgb_event_summary(events),
    accessibility = accessibility
  )
}

make_modern_validation <- function(beta, recipe, modern_grid, tree, traits, species_origin) {
  modern_sites <- modern_grid[sample(seq_len(nrow(modern_grid)), min(900, nrow(modern_grid))), , drop = FALSE]
  modern_sites$site_id <- sprintf("site%04d", seq_len(nrow(modern_sites)))
  rownames(modern_sites) <- modern_sites$site_id
  modern_suit <- hee_project_hmsc_table(
    beta = beta, X_past = modern_sites, recipe = recipe, link = "probit",
    id_cols = c("site_id", "time_ma", "lon", "lat", "cell_id")
  )
  modern_proj <- hee_combine(
    modern_suit,
    phylo_mask = hee_phylo_time_mask(times = 0, species_origin = species_origin)
  )
  wide_prob <- reshape(modern_proj[, c("site_id", "species", "probability")],
                       idvar = "site_id", timevar = "species", direction = "wide")
  prob <- as.matrix(wide_prob[, paste0("probability.", species)])
  colnames(prob) <- species
  rownames(prob) <- wide_prob$site_id
  comm <- matrix(rbinom(length(prob), 1, pmin(pmax(prob, 0), 1)),
                 nrow = nrow(prob), dimnames = dimnames(prob))
  env_now <- modern_sites[rownames(comm), variables, drop = FALSE]
  dat <- hee_prepare_data(comm = comm, env_now = env_now, tree = tree,
                          traits = traits, response_type = "presence")
  input <- hee_as_hmsc_input(dat, XFormula = stats::as.formula(
    paste("~", paste(variables, collapse = " + "))
  ), TrFormula = ~ body_size + dispersal + thermal_tolerance +
    moisture_affinity + coastal_affinity + elevation_tolerance)
  eval <- hee_evaluate_hmsc(as.numeric(comm), as.numeric(prob))
  threshold <- hee_threshold_train(as.numeric(comm), as.numeric(prob))
  list(
    sites = modern_sites,
    suitability = modern_suit,
    projection = modern_proj,
    comm = comm,
    data = dat,
    audit = hee_audit_data(dat),
    hmsc_input = input,
    evaluation = eval,
    threshold = threshold,
    folds = list(
      random = hee_cv_random(nrow(modern_sites), k = 5, seed = 20260608),
      spatial = hee_cv_spatial_block(modern_sites[, c("lon", "lat")], k = 5),
      environment = hee_cv_environment_block(modern_sites[, variables], variable = "MAT_pohl_C", k = 5)
    )
  )
}

make_global_case <- function() {
  case_dir <- file.path(root, "case01_global540_20")
  dir.create(case_dir, recursive = TRUE, showWarnings = FALSE)
  tree <- make_time_tree(1)
  species_origin <- hee_estimate_species_origin(tree)
  traits <- make_traits(1)
  beta <- make_beta(traits, 1)
  tip_ranges <- make_tip_ranges(1)
  bgb <- make_bgb_objects(tip_ranges, species_origin)
  formula <- stats::as.formula(paste("~", paste(variables, collapse = " + ")))
  modern <- hee_make_paleo_grid(cube, times = 0, variables = variables, land_only = TRUE)
  recipe <- hee_lock_recipe(modern, formula = formula)
  roles <- hee_tag_predictors(variables, stats::setNames(rep("deep_time_exogenous", length(variables)), variables))
  leakage <- hee_detect_leakage(formula, roles)
  paleo_grid <- hee_make_paleo_grid(cube, times = times, variables = variables, land_only = TRUE)
  extrapolation <- hee_extrapolation_risk(recipe, paleo_grid, threshold = 3)
  extrapolation_weight <- extrapolation[, c("cell_id", "time_ma", "extrapolation_weight"), drop = FALSE]
  suitability <- hee_project_hmsc_timecube(beta = beta, env_cube = cube, recipe = recipe,
                                           times = times, variables = variables, link = "probit")
  phylo_mask <- hee_phylo_time_mask(times = times, species_origin = species_origin)
  projection <- hee_combine(suitability, phylo_mask = phylo_mask,
                            extrapolation_weight = extrapolation_weight)
  richness <- hee_richness(projection)
  turnover <- hee_turnover(projection)
  refugia <- hee_refugia(richness)
  validation <- make_modern_validation(beta, recipe, modern, tree, traits, species_origin)

  fossil_candidates <- projection[projection$probability > 0.45 & projection$phylo_existence > 0, ]
  fossil_candidates <- fossil_candidates[order(fossil_candidates$time_ma, -fossil_candidates$probability), ]
  fossil_candidates <- fossil_candidates[!duplicated(fossil_candidates$time_ma), ]
  fossils <- data.frame(
    species = fossil_candidates$species,
    lon = fossil_candidates$lon,
    lat = fossil_candidates$lat,
    age_min = pmax(fossil_candidates$time_ma - 2.5, 0),
    age_max = pmin(fossil_candidates$time_ma + 2.5, 540),
    evidence_type = "simulated_full_axis_fossil",
    certainty = 0.85,
    stringsAsFactors = FALSE
  )
  fossil_validation <- hee_validate_fossils(fossils, projection, threshold = validation$threshold$threshold)

  output <- list(
    case = "case01_global540_20",
    source_rds = source_rds,
    source_level = "4-degree Phanerozoic land-harmonized compact RDS",
    variables = variables,
    times = times,
    species = species,
    clade = clade,
    tree = tree,
    species_origin = species_origin,
    traits = traits,
    beta = beta,
    tip_ranges = tip_ranges,
    bgb = bgb,
    recipe = recipe,
    predictor_roles = roles,
    leakage = leakage,
    timecube_index = hee_index_timecube(cube, variables = variables, times = times),
    timecube_check = hee_check_timecube(cube, variables = variables, times = times),
    paleo_grid = paleo_grid,
    extrapolation = extrapolation,
    suitability = suitability,
    projection = projection,
    richness = richness,
    turnover = turnover,
    refugia = refugia,
    modern_validation = validation,
    fossils = fossils,
    fossil_validation = fossil_validation,
    model_set = hee_define_hmsc_models(formula, traits = TRUE, phylogeny = TRUE,
                                       spatial = TRUE, historical_exogenous = variables)
  )
  saveRDS(output, file.path(case_dir, "case01_global540_20.rds"), compress = "xz")
  ape::write.tree(tree, file.path(case_dir, "dated_phylogeny_540Ma.tre"))
  utils::write.csv(data.frame(species = species, clade = clade[species],
                              origin_ma = species_origin[species]),
                   file.path(case_dir, "species_origin.csv"), row.names = FALSE)
  utils::write.csv(beta, file.path(case_dir, "beta.csv"))
  utils::write.csv(traits, file.path(case_dir, "traits.csv"))
  utils::write.csv(roles, file.path(case_dir, "predictor_roles.csv"), row.names = FALSE)
  utils::write.csv(output$timecube_index, file.path(case_dir, "timecube_index.csv"), row.names = FALSE)
  utils::write.csv(bgb$compare, file.path(case_dir, "bgb_model_compare.csv"), row.names = FALSE)
  invisible(output)
}

make_tracks <- function(times) {
  base <- data.frame(
    region = region_names,
    track_id = paste0("track", seq_along(region_names)),
    lon0 = c(-120, -65, 15, 70, 120, -20),
    lat0 = c(-10, 25, -35, 10, 45, 0)
  )
  tracks <- merge(expand.grid(region = region_names, time_ma = times,
                              KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE),
                  base, by = "region", sort = FALSE)
  phase <- as.numeric(factor(tracks$region))
  tracks$paleo_lon <- tracks$lon0 +
    65 * sin((540 - tracks$time_ma) / 95 + phase * 0.7) -
    0.10 * tracks$time_ma
  tracks$paleo_lat <- tracks$lat0 +
    32 * cos((540 - tracks$time_ma) / 120 + phase * 0.5) +
    0.03 * tracks$time_ma
  tracks$paleo_lon <- ((tracks$paleo_lon + 180) %% 360) - 180
  tracks$paleo_lat <- pmax(pmin(tracks$paleo_lat, 84), -84)
  tracks
}

make_track_dynamic_filter <- function(track_projection) {
  base <- unique(track_projection[, intersect(c("species", "region", "time_ma", "probability"),
                                              names(track_projection)), drop = FALSE])
  base$dynamic_weight <- pmin(1, pmax(0.15, 0.35 + 0.65 * base$probability))
  base[, c("species", "region", "time_ma", "dynamic_weight")]
}

make_track_case <- function() {
  case_dir <- file.path(root, "case02_tracks540_20")
  dir.create(case_dir, recursive = TRUE, showWarnings = FALSE)
  tree <- make_time_tree(2)
  species_origin <- hee_estimate_species_origin(tree)
  traits <- make_traits(2)
  beta <- make_beta(traits, 2)
  tip_ranges <- make_tip_ranges(2)
  bgb <- make_bgb_objects(tip_ranges, species_origin)
  formula <- stats::as.formula(paste("~", paste(variables, collapse = " + ")))
  modern <- hee_make_paleo_grid(cube, times = 0, variables = variables, land_only = TRUE)
  recipe <- hee_lock_recipe(modern, formula = formula)
  tracks <- make_tracks(times)
  plate_points <- hee_reconstruct_points(tracks, precomputed = hee_import_plate_points(tracks))
  track_env <- hee_extract_paleoenv_track(plate_points, cube, variables = variables,
                                          land_only = TRUE)
  track_env$nearest_cell_id <- track_env$cell_id
  track_env$cell_id <- track_env$track_id
  track_extrapolation <- hee_extrapolation_risk(recipe, track_env, threshold = 3)
  track_suitability <- hee_project_hmsc_table(
    beta = beta, X_past = track_env, recipe = recipe, link = "probit",
    id_cols = c("track_id", "region", "time_ma", "matched_time_ma", "lon", "lat",
                "paleo_lon", "paleo_lat", "cell_id", "nearest_cell_id")
  )
  phylo_mask <- hee_phylo_time_mask(times = times, species_origin = species_origin)
  dynamic_filter <- NULL
  first_projection <- hee_combine(track_suitability, accessibility = bgb$accessibility,
                                  phylo_mask = phylo_mask)
  dynamic_filter <- make_track_dynamic_filter(first_projection)
  track_projection <- hee_combine(track_suitability, accessibility = bgb$accessibility,
                                  phylo_mask = phylo_mask,
                                  dynamic_filter = dynamic_filter,
                                  extrapolation_weight = track_extrapolation[, c("cell_id", "time_ma", "extrapolation_weight")])
  track_richness <- hee_richness(track_projection, group_cols = c("track_id", "region", "time_ma"))
  track_turnover <- hee_turnover(track_projection)
  track_refugia <- hee_refugia(track_richness, group_col = "track_id")
  validation <- make_modern_validation(beta, recipe, modern, tree, traits, species_origin)
  dispersal <- hee_make_dispersal_matrices(region_names, scale = 2)

  output <- list(
    case = "case02_tracks540_20",
    source_rds = source_rds,
    source_level = "4-degree Phanerozoic land-harmonized compact RDS",
    variables = variables,
    times = times,
    species = species,
    clade = clade,
    tree = tree,
    species_origin = species_origin,
    traits = traits,
    beta = beta,
    tip_ranges = tip_ranges,
    bgb = bgb,
    recipe = recipe,
    tracks = tracks,
    plate_points = plate_points,
    track_environment = track_env,
    track_extrapolation = track_extrapolation,
    dispersal_matrix = dispersal,
    dynamic_filter = dynamic_filter,
    suitability = track_suitability,
    projection = track_projection,
    richness = track_richness,
    turnover = track_turnover,
    refugia = track_refugia,
    modern_validation = validation,
    model_set = hee_define_hmsc_models(formula, traits = TRUE, phylogeny = TRUE,
                                       spatial = FALSE, historical_exogenous = variables)
  )
  saveRDS(output, file.path(case_dir, "case02_tracks540_20.rds"), compress = "xz")
  ape::write.tree(tree, file.path(case_dir, "dated_phylogeny_540Ma.tre"))
  utils::write.csv(data.frame(species = species, clade = clade[species],
                              origin_ma = species_origin[species]),
                   file.path(case_dir, "species_origin.csv"), row.names = FALSE)
  utils::write.csv(beta, file.path(case_dir, "beta.csv"))
  utils::write.csv(traits, file.path(case_dir, "traits.csv"))
  utils::write.csv(tracks, file.path(case_dir, "plate_corrected_tracks_540_0Ma.csv"), row.names = FALSE)
  utils::write.csv(track_env, file.path(case_dir, "track_environment_540_0Ma.csv"), row.names = FALSE)
  utils::write.csv(bgb$compare, file.path(case_dir, "bgb_model_compare.csv"), row.names = FALSE)
  invisible(output)
}

case01 <- make_global_case()
case02 <- make_track_case()

writeLines(c(
  "# Full Phanerozoic Deep-Time HmscEcoEvo Cases",
  "",
  "These two cases use the 4-degree Phanerozoic land-harmonized compact RDS product supplied with the project.",
  "",
  "Time axis: 540, 535, ..., 0 Ma (109 slices).",
  "Species: 20 simulated taxa with a dated 540 Ma phylogeny and unequal species-origin ages.",
  "",
  "case01_global540_20: global land-grid projection, extrapolation risk, phylogenetic time mask, richness, turnover, refugia, model validation, simulated fossil validation, and BioGeoBEARS-style precomputed accessibility/event objects.",
  "case02_tracks540_20: plate-corrected track projection for six regions across the full 540-0 Ma axis, model-averaged accessibility, dynamic dispersal filtering, richness, turnover, refugia, and validation objects.",
  "",
  paste("Source RDS:", source_rds)
), file.path(root, "README.md"))

cat("Generated full 540-0 Ma deep-time cases under ", root, "\n", sep = "")
cat("Case 01 rows: projection=", nrow(case01$projection),
    " richness=", nrow(case01$richness),
    " turnover=", nrow(case01$turnover), "\n", sep = "")
cat("Case 02 rows: projection=", nrow(case02$projection),
    " richness=", nrow(case02$richness),
    " turnover=", nrow(case02$turnover), "\n", sep = "")
