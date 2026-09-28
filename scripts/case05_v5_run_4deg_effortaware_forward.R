#!/usr/bin/env Rscript

# Case05 v5: 4-degree, effort-aware, global dynamic-Earth scenario.
#
# This runner is deliberately separate from the paused 1-degree Case05 run.
# It refits a modern HMSC on an explicitly labelled target-group
# presence-background design, reconstructs ancestral Beta responses from real
# retained HMSC posterior draws, and propagates them through a plate-carried,
# global 4-degree carrier grid.  It is a rapid scientific scenario, not a
# formal detection-corrected palaeodistribution reconstruction: Plant200 has
# no independent global plant effort surface or independent endpoint range
# data in this workspace.

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || identical(x, "") || is.na(x)) y else x
}
parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!grepl("^--", arg)) next
    z <- sub("^--", "", arg)
    at <- regexpr("=", z, fixed = TRUE)
    if (at < 0L) out[[z]] <- TRUE else {
      out[[substr(z, 1L, at - 1L)]] <- substr(z, at + 1L, nchar(z))
    }
  }
  out
}
as_int <- function(x, default) {
  z <- suppressWarnings(as.integer(x)[1L])
  if (is.finite(z)) z else default
}
as_num <- function(x, default) {
  z <- suppressWarnings(as.numeric(x)[1L])
  if (is.finite(z)) z else default
}
as_flag <- function(x, default = FALSE) {
  if (is.null(x) || !length(x)) return(default)
  tolower(as.character(x)[1L]) %in% c("true", "t", "1", "yes", "y")
}
safe_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}
write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(x, path)
  } else {
    utils::write.csv(x, path, row.names = FALSE)
  }
  invisible(path)
}

# Full HMSC posteriors can contain millions of Beta/Gamma/rho values.  For this
# long Case05 run, audit convergence from an evenly spaced, reproducible sample
# of retained Beta parameters rather than materialising every posterior value
# into a multi-gigabyte diagnostic data frame.
compact_hmsc_mcmc_diagnostics <- function(hM, max_parameters = 32L) {
  if (!requireNamespace("coda", quietly = TRUE) || is.null(hM$postList) ||
      length(hM$postList) < 2L || length(hM$postList[[1L]]) < 2L) {
    return(list(ess = data.frame(), rhat = data.frame(), parameter_index = integer()))
  }
  first_beta <- as.numeric(hM$postList[[1L]][[1L]]$Beta)
  if (!length(first_beta)) {
    return(list(ess = data.frame(), rhat = data.frame(), parameter_index = integer()))
  }
  index <- unique(as.integer(round(seq(1, length(first_beta),
                                       length.out = min(max_parameters, length(first_beta))))))
  chain_mats <- lapply(hM$postList, function(chain) {
    out <- vapply(chain, function(state) {
      beta <- as.numeric(state$Beta)
      beta[index]
    }, numeric(length(index)))
    t(out)
  })
  ess_rows <- vector("list", length(index))
  rhat_rows <- vector("list", length(index))
  for (j in seq_along(index)) {
    chains <- lapply(chain_mats, function(x) coda::mcmc(x[, j]))
    ml <- coda::mcmc.list(chains)
    label <- paste0("Beta_flat_index_", index[[j]])
    ess_rows[[j]] <- data.frame(parameter = label,
                                ESS = as.numeric(coda::effectiveSize(ml)[[1L]]))
    rhat_rows[[j]] <- data.frame(
      parameter = label,
      Rhat = tryCatch(coda::gelman.diag(ml, autoburnin = FALSE)$psrf[1L, 1L],
                      error = function(e) NA_real_)
    )
  }
  list(ess = do.call(rbind, ess_rows), rhat = do.call(rbind, rhat_rows),
       parameter_index = index)
}
time_slug <- function(x) formatC(as.numeric(x), format = "f", digits = 4,
                                 drop0trailing = TRUE)

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
file_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
if (!length(script_path)) script_path <- "scripts/case05_v5_run_4deg_effortaware_forward.R"
pkg_root <- normalizePath(
  cfg$pkg_root %||% file.path(dirname(normalizePath(script_path[[1L]],
                                                     winslash = "/", mustWork = FALSE)), ".."),
  winslash = "/", mustWork = TRUE
)
workspace <- normalizePath(file.path(pkg_root, ".."), winslash = "/", mustWork = TRUE)
output <- safe_dir(cfg$output %||% file.path(
  workspace, "outputs", "HmscEcoEvo",
  paste0("case05_plant200_4deg_effortaware_forward_", format(Sys.time(), "%Y%m%d"))
))
grid_deg <- as_num(cfg$grid_deg, 4)
if (!identical(grid_deg, 4)) {
  stop("This rapid runner is deliberately fixed to a 4-degree grid.", call. = FALSE)
}
seed <- as_int(cfg$seed, 20260923L)
n_chains <- as_int(cfg$n_chains, 4L)
mcmc_samples <- as_int(cfg$mcmc_samples, 100L)
mcmc_transient <- as_int(cfg$mcmc_transient, 200L)
mcmc_thin <- as_int(cfg$mcmc_thin, 10L)
n_response_draws <- as_int(cfg$n_response_draws, 12L)
hmsc_workers <- as_int(cfg$hmsc_workers, min(4L, parallel::detectCores(logical = FALSE)))
internal_dt_myr <- as_num(cfg$internal_dt_myr, 0.5)
min_records <- as_int(cfg$min_target_group_records, 1L)
root_patch_cells <- as_int(cfg$root_patch_cells, 8L)
root_occupancy <- as_num(cfg$root_occupancy, 0.90)
local_emigration_rate <- as_num(cfg$local_emigration_rate, 0.20)
ldd_emigration_rate <- as_num(cfg$ldd_emigration_rate, 0.0005)
diffusion_precision <- as_num(cfg$diffusion_precision_myr_per_rad2, 100)
topographic_cost_strength <- as_num(cfg$topographic_cost_strength, 0.65)
establishment_intercept <- as_num(cfg$establishment_intercept, 0.25)
establishment_slope <- as_num(cfg$establishment_slope, 1.0)
# These are one-Myr process-reference parameters.  A value near 1.35 implies
# ca. 21% local loss per Myr and makes a 325-Myr continuous scenario collapse
# regardless of dispersal.  The reference scenario instead uses high, but not
# perfect, local persistence and leaves environmental mismatch to increase
# loss through the eta slope.  They remain documented scenario parameters,
# not rates estimated from the modern presence-background matrix.
persistence_intercept <- as_num(cfg$persistence_intercept, 3.5)
persistence_slope <- as_num(cfg$persistence_slope, 1.2)
rerun_hmsc <- as_flag(cfg$rerun_hmsc, FALSE)
word_atlas_mode <- tolower(as.character(cfg$word_atlas_mode %||% "full")[1L])
build_word_atlas <- as_flag(cfg$build_word_atlas, FALSE)
if (!word_atlas_mode %in% c("full", "representative")) {
  stop("word_atlas_mode must be 'full' or 'representative'.", call. = FALSE)
}
if (n_chains < 2L || mcmc_samples < 20L || mcmc_transient < 20L ||
    mcmc_thin < 1L || n_response_draws < 2L || internal_dt_myr <= 0 ||
    root_patch_cells < 1L || root_occupancy <= 0 || root_occupancy > 1 ||
    local_emigration_rate < 0 || ldd_emigration_rate < 0 ||
    diffusion_precision <= 0) {
  stop("Invalid rapid-run configuration.", call. = FALSE)
}
if (!requireNamespace("Hmsc", quietly = TRUE) ||
    !requireNamespace("ape", quietly = TRUE) ||
    !requireNamespace("Matrix", quietly = TRUE) ||
    !requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("officer", quietly = TRUE) ||
    !requireNamespace("pkgload", quietly = TRUE)) {
  stop("Hmsc, ape, Matrix, ggplot2, officer, and pkgload are required.",
       call. = FALSE)
}
pkgload::load_all(pkg_root, quiet = TRUE)
set.seed(seed)

input_dir <- normalizePath(cfg$plant_input_dir %||% file.path(
  workspace, "HmscEcoEvo_API8_20260818", "data_external",
  "plant_200_multifamily_extant_20260820", "prepared_inputs", "case03_data_raw",
  "Plant200_multifamily_extant_1deg_traits"
), winslash = "/", mustWork = TRUE)
palaeo_root <- normalizePath(cfg$palaeo_root %||% file.path(
  pkg_root, "derived_inputs", "case04_plant200_palaeo_1deg_slices_20260919"
), winslash = "/", mustWork = TRUE)
carrier_root <- normalizePath(cfg$carrier_root %||% file.path(
  pkg_root, "derived_inputs", "case05_v2_paleomap_h3_carrier_325_0Ma_1deg"
), winslash = "/", mustWork = TRUE)

dirs <- list(
  config = safe_dir(file.path(output, "00_config")),
  modern = safe_dir(file.path(output, "01_modern_observation_design")),
  hmsc = safe_dir(file.path(output, "02_hmsc_4deg")),
  evolution = safe_dir(file.path(output, "03_ancestral_response")),
  earth = safe_dir(file.path(output, "04_dynamic_earth_4deg")),
  draw = safe_dir(file.path(output, "05_response_draw_metrics")),
  summaries = safe_dir(file.path(output, "06_summaries")),
  maps = safe_dir(file.path(output, "07_png_maps")),
  rasters = safe_dir(file.path(output, "08_geotiff_maps")),
  report = safe_dir(file.path(output, "09_report"))
)

write_csv(data.frame(
  field = c("case_id", "analysis_mode", "grid_deg", "time_direction",
            "observation_design", "zero_semantics", "endpoint_role",
            "plate_transport", "active_dispersal", "palaeo_time_convention",
            "scientific_boundary"),
  value = c(
    "case05_plant200_4deg_effortaware_forward",
    "rapid_process_constrained_forward_scenario",
    grid_deg, "old_to_young", "target_group_presence_background",
    "non_record_in_effort_qualified_target_group_cell_not_confirmed_absence",
    "spatial_holdout_target_group_validation_not_independent_endpoint_likelihood",
    "stable_PALEOMAP_H3_carrier_identity_aggregated_to_4deg_target_cells",
    "Arias_spherical_local_plus_explicit_rare_long_distance_land_to_land_scenario",
    "time_ma_larger_values_are_older",
    paste(
      "This run is a 4-degree rapid scenario. It propagates real HMSC Beta",
      "posterior draws and a dated screening tree, but the available Plant200",
      "records provide target-group presence-background information rather than",
      "confirmed absences or independent global sampling effort. Outputs are not",
      "a detection-corrected or uniquely reconstructed historical distribution."
    )
  ), stringsAsFactors = FALSE
), file.path(dirs$config, "case05_v5_scientific_contract.csv"))

# Keep every process-setting explicit.  These are scenario parameters, rather
# than estimates from the one-time target-group presence-background survey.
write_csv(data.frame(
  parameter = c(
    "grid_deg", "n_chains", "mcmc_samples", "mcmc_transient", "mcmc_thin",
    "n_response_draws", "internal_dt_myr", "root_patch_cells", "root_occupancy",
    "local_emigration_rate_per_myr", "ldd_emigration_rate_per_myr",
    "diffusion_precision_myr_per_rad2", "topographic_cost_strength",
    "establishment_intercept", "establishment_slope",
    "persistence_intercept", "persistence_slope"
  ),
  value = c(
    grid_deg, n_chains, mcmc_samples, mcmc_transient, mcmc_thin,
    n_response_draws, internal_dt_myr, root_patch_cells, root_occupancy,
    local_emigration_rate, ldd_emigration_rate, diffusion_precision,
    topographic_cost_strength, establishment_intercept, establishment_slope,
    persistence_intercept, persistence_slope
  ),
  interpretation = c(
    "fixed 4-degree global palaeoland analysis grid",
    "independent HMSC MCMC chains", "retained samples per chain",
    "discarded MCMC transient", "MCMC thinning interval",
    "joint HMSC-response / ancestral-response draws propagated",
    "maximum continuous-time occupancy update step in Myr",
    "number of connected root-patch cells", "initial root occupancy multiplier",
    "Arias-style local spherical-kernel emigration scenario",
    "rare explicit long-distance dispersal scenario",
    "spherical diffusion precision", "continuous topographic resistance weight",
    "conditional-establishment baseline", "environmental-support effect on establishment",
    "conditional-persistence baseline", "environmental-support effect on persistence"
  ), stringsAsFactors = FALSE
), file.path(dirs$config, "case05_v5_process_scenario_parameters.csv"))

bin4 <- function(lon, lat) {
  lon <- as.numeric(lon); lat <- as.numeric(lat)
  nlon <- 90L; nlat <- 45L
  ix <- pmin(nlon - 1L, pmax(0L, floor((lon + 180) / 4)))
  iy <- pmin(nlat - 1L, pmax(0L, floor((lat + 90) / 4)))
  data.frame(
    grid_lon_index = as.integer(ix), grid_lat_index = as.integer(iy),
    grid_cell_id = sprintf("grid%03d_%03d", ix + 1L, iy + 1L),
    lon = -180 + (ix + 0.5) * 4,
    lat = -90 + (iy + 0.5) * 4,
    stringsAsFactors = FALSE
  )
}
weighted_group_mean <- function(value, group, weight, levels) {
  value <- as.numeric(value); weight <- as.numeric(weight)
  out <- rep(NA_real_, length(levels)); names(out) <- levels
  ok <- is.finite(value) & is.finite(weight) & weight > 0 & !is.na(group)
  if (!any(ok)) return(out)
  numerator <- rowsum(value[ok] * weight[ok], group[ok], reorder = FALSE)
  denominator <- rowsum(weight[ok], group[ok], reorder = FALSE)
  out[rownames(numerator)] <- as.numeric(numerator[, 1L] / denominator[, 1L])
  out
}
aggregate_environment_4deg <- function(x, env_cols) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  keep <- is.finite(x$H_state) & x$H_state > 0 &
    is.finite(x$land_area_km2) & x$land_area_km2 > 0
  x <- x[keep, , drop = FALSE]
  b <- bin4(x$lon, x$lat)
  level <- unique(b$grid_cell_id)
  base <- b[match(level, b$grid_cell_id), c("grid_cell_id", "grid_lon_index",
                                             "grid_lat_index", "lon", "lat"), drop = FALSE]
  area <- rowsum(x$land_area_km2, b$grid_cell_id, reorder = FALSE)
  base$land_area_km2 <- as.numeric(area[base$grid_cell_id, 1L])
  base$H_state <- 1L
  for (nm in env_cols) base[[nm]] <- weighted_group_mean(
    x[[nm]], b$grid_cell_id, x$land_area_km2, base$grid_cell_id
  )
  base
}

message("Preparing 4-degree target-group observation design.")
comm <- as.matrix(utils::read.csv(file.path(input_dir, "comm.csv"),
                                  check.names = FALSE, row.names = 1L))
storage.mode(comm) <- "double"
sites_1deg <- utils::read.csv(file.path(input_dir, "sites.csv"),
                              stringsAsFactors = FALSE)
grid_1deg <- utils::read.csv(file.path(
  input_dir, "prediction_grid_0Ma_1deg_land_complete_environment.csv"
), stringsAsFactors = FALSE)
records <- utils::read.csv(file.path(input_dir,
                                     "occurrence_records_used_1deg_indexed.csv"),
                           stringsAsFactors = FALSE)
tree <- ape::read.tree(file.path(input_dir, "tree.tre"))
if (!ape::is.ultrametric(tree) || is.null(tree$edge.length)) {
  stop("Case05 v5 requires an ultrametric dated tree.", call. = FALSE)
}
if (!setequal(colnames(comm), tree$tip.label)) {
  stop("Plant200 community species and tree tips are not identical.", call. = FALSE)
}
modern_design <- HmscEcoEvo::hee_hmsc_target_group_design(
  sites = sites_1deg, Y = comm, grid = grid_1deg, records = records,
  grid_deg = grid_deg, min_target_group_records = min_records,
  record_species_col = "selected_species_id"
)
write_csv(modern_design$metadata,
          file.path(dirs$modern, "target_group_observation_design_metadata.csv"))
write_csv(modern_design$grid_cells,
          file.path(dirs$modern, "modern_4deg_complete_land_grid_effort.csv"))
write_csv(modern_design$cell_map,
          file.path(dirs$modern, "modern_1deg_to_4deg_cell_map.csv"))

earth_index <- utils::read.csv(file.path(palaeo_root,
                                          "case04_palaeo_earth_slice_index.csv"),
                               stringsAsFactors = FALSE)
earth_index$time_ma <- as.numeric(earth_index$time_ma)
earth0_path <- earth_index$file[which(earth_index$time_ma == 0)[1L]]
if (!length(earth0_path) || !file.exists(earth0_path)) {
  stop("The 0 Ma palaeoenvironment slice is missing.", call. = FALSE)
}
environment_axes <- c(
  "MAT_pohl_C", "MAP_pohl_mm_yr", "elevation_m",
  "T_seasonality_pohl_sd", "P_seasonality_pohl_sd",
  "moisture_availability_index_z", "wetland_potential_index_z", "relief_3x3_m"
)
earth0 <- as.data.frame(readRDS(earth0_path))
if (!all(c("time_ma", "cell_id", "lon", "lat", "H_state", "land_area_km2",
           environment_axes) %in% names(earth0))) {
  stop("The 0 Ma Earth slice lacks required environmental axes.", call. = FALSE)
}
modern_env4 <- aggregate_environment_4deg(earth0, environment_axes)
# Coordinates in the observation design are the canonical 4-degree cell
# centres.  Join only environmental columns here, otherwise merge() creates
# lon.x/lon.y and silently removes the coordinate columns used by the spatial
# block split below.
modern_env4_for_hmsc <- modern_env4[, c("grid_cell_id", environment_axes),
                                    drop = FALSE]
modern_hmsc <- merge(modern_design$sites, modern_env4_for_hmsc,
                     by.x = "site_id", by.y = "grid_cell_id", all.x = TRUE,
                     sort = FALSE)
modern_hmsc <- modern_hmsc[match(rownames(modern_design$Y), modern_hmsc$site_id), ,
                           drop = FALSE]
if (anyNA(modern_hmsc[, environment_axes, drop = FALSE])) {
  bad <- modern_hmsc$site_id[!stats::complete.cases(modern_hmsc[, environment_axes, drop = FALSE])]
  stop("Modern 4-degree training cells lack complete environment: ",
       paste(utils::head(bad, 5L), collapse = ", "), call. = FALSE)
}
Y4 <- modern_design$Y[modern_hmsc$site_id, , drop = FALSE]
rownames(Y4) <- modern_hmsc$site_id

# Spatial holdout uses 20-degree blocks. Species with only one occupied cell
# are retained in training, because otherwise their environmental response is
# not estimable; they are flagged as non-evaluable in holdout diagnostics.
block_id <- paste(floor((modern_hmsc$lon + 180) / 20),
                  floor((modern_hmsc$lat + 90) / 20), sep = "_")
set.seed(seed)
block_levels <- unique(block_id)
candidate_holdout_blocks <- sample(block_levels,
                                   max(1L, floor(length(block_levels) * 0.2)))
holdout <- block_id %in% candidate_holdout_blocks
for (j in seq_len(ncol(Y4))) {
  present <- which(Y4[, j] > 0)
  if (length(present) && !any(!holdout[present])) {
    holdout[present[[1L]]] <- FALSE
  }
}
train <- !holdout
if (sum(train) < 30L || sum(holdout) < 10L) {
  stop("Spatial block split is too small after rare-species protection.", call. = FALSE)
}
write_csv(data.frame(
  site_id = modern_hmsc$site_id, lon = modern_hmsc$lon, lat = modern_hmsc$lat,
  spatial_block_20deg = block_id, split = ifelse(train, "training", "holdout"),
  target_group_record_count = modern_hmsc$target_group_record_count,
  stringsAsFactors = FALSE
), file.path(dirs$modern, "modern_4deg_spatial_holdout_split.csv"))
write_csv(cbind(modern_hmsc, observed_richness = rowSums(Y4)),
          file.path(dirs$modern, "modern_4deg_hmsc_sites_environment.csv"))

effort_z <- as.numeric(scale(modern_hmsc$log_target_group_record_count))
X4 <- data.frame(modern_hmsc[, environment_axes, drop = FALSE],
                 observation_log_target_group_effort_z = effort_z,
                 check.names = FALSE)
rownames(X4) <- modern_hmsc$site_id
saveRDS(list(Y = Y4, X = X4, sites = modern_hmsc, train = train,
             holdout = holdout, environment_axes = environment_axes),
        file.path(dirs$modern, "modern_4deg_hmsc_design.rds"), compress = "gzip")

model_path <- file.path(dirs$hmsc, "case05_v5_hmsc_4deg_target_group_model.rds")
if (file.exists(model_path) && !rerun_hmsc) {
  message("Reusing existing 4-degree HMSC model: ", model_path)
  hM <- readRDS(model_path)
} else {
  message("Fitting 4-degree target-group HMSC with ", sum(train),
          " training cells and ", ncol(Y4), " focal species.")
  hM <- Hmsc::Hmsc(
    Y = Y4[train, , drop = FALSE],
    XData = X4[train, , drop = FALSE],
    XFormula = stats::as.formula(
      paste("~", paste(names(X4), collapse = " + "))
    ),
    XScale = FALSE,
    distr = "probit"
  )
  hM <- Hmsc::sampleMcmc(
    hM, samples = mcmc_samples, transient = mcmc_transient,
    thin = mcmc_thin, nChains = n_chains,
    nParallel = max(1L, min(hmsc_workers, n_chains)), verbose = 0
  )
  saveRDS(hM, model_path, compress = "gzip")
}

# The holdout calculation uses HMSC fixed effects only. It is a spatially
# held-out target-group diagnostic, not a detection-corrected validation.
all_beta_with_intercept <- HmscEcoEvo::hee_hmsc_beta_posterior_draws(
  hM, n_draws = min(n_response_draws, 12L), include_intercept = TRUE,
  seed = seed + 11L
)
available_draws <- unique(all_beta_with_intercept$response_draw)
prediction_holdout <- matrix(0, nrow = sum(holdout), ncol = ncol(Y4),
                             dimnames = list(rownames(Y4)[holdout], colnames(Y4)))
for (draw_id in available_draws) {
  z <- all_beta_with_intercept[all_beta_with_intercept$response_draw == draw_id,
                               , drop = FALSE]
  z <- z[match(colnames(Y4), z$lineage), , drop = FALSE]
  beta_all <- as.matrix(z[, c("hmsc_intercept_original", names(X4)), drop = FALSE])
  storage.mode(beta_all) <- "double"
  x_holdout <- as.matrix(X4[holdout, , drop = FALSE])
  storage.mode(x_holdout) <- "double"
  eta <- sweep(x_holdout %*% t(beta_all[, -1L, drop = FALSE]),
               2L, beta_all[, 1L], "+")
  prediction_holdout <- prediction_holdout + stats::pnorm(eta) / length(available_draws)
}
holdout_score <- HmscEcoEvo::hee_modern_endpoint_score(
  prediction_holdout, Y4[holdout, , drop = FALSE], data_role = "spatial_holdout"
)
holdout_score$summary$scientific_boundary <-
  "spatial_holdout_target_group_presence_background_diagnostic_not_independent_validation"
write_csv(holdout_score$summary, file.path(dirs$hmsc, "spatial_holdout_summary.csv"))
write_csv(holdout_score$species, file.path(dirs$hmsc, "spatial_holdout_species.csv"))
write_csv(holdout_score$cells, file.path(dirs$hmsc, "spatial_holdout_cells.csv"))
# Do not materialise all 4,000 retained HMSC prediction arrays for this
# 200-species run: it requires several GB and does not add independent evidence
# beyond the spatial holdout already computed above.  Store the model-fit
# decision explicitly rather than silently substituting an in-sample score.
write_csv(data.frame(
  status = "not_materialised_for_global_case05_memory_budget",
  diagnostic_used = "spatial_holdout_target_group_presence_background",
  scientific_boundary = "holdout_is_diagnostic_not_independent_endpoint_likelihood",
  stringsAsFactors = FALSE
), file.path(dirs$hmsc, "training_model_fit.csv"))

# Preserve a memory-bounded convergence audit and posterior point summaries.
# The package-level hee_hmsc_diagnostics() remains available for smaller models;
# this global 200-species run deliberately avoids expanding its entire posterior.
hmsc_diagnostics <- compact_hmsc_mcmc_diagnostics(hM, max_parameters = 32L)
write_csv(hmsc_diagnostics$ess,
          file.path(dirs$hmsc, "hmsc_mcmc_effective_sample_size_compact.csv"))
write_csv(hmsc_diagnostics$rhat,
          file.path(dirs$hmsc, "hmsc_mcmc_rhat_compact.csv"))
write_csv(data.frame(parameter_index = hmsc_diagnostics$parameter_index),
          file.path(dirs$hmsc, "hmsc_mcmc_diagnostic_parameter_index.csv"))
saveRDS(list(
  diagnostic_scope = "32_evenly_spaced_retained_Beta_parameters; memory_bounded",
  full_posterior_model = normalizePath(model_path, winslash = "/", mustWork = TRUE),
  response_posterior_draw_file = file.path(dirs$hmsc,
                                            "hmsc_environmental_beta_posterior_draws.rds")
), file.path(dirs$hmsc, "hmsc_posterior_point_summaries.rds"), compress = "gzip")

tip_beta <- HmscEcoEvo::hee_hmsc_beta_posterior_draws(
  hM, n_draws = n_response_draws, axes = environment_axes,
  include_intercept = FALSE, seed = seed + 29L
)
write_csv(tip_beta, file.path(dirs$hmsc, "hmsc_environmental_beta_posterior_draws.csv"))
saveRDS(tip_beta, file.path(dirs$hmsc, "hmsc_environmental_beta_posterior_draws.rds"),
        compress = "gzip")

message("Preparing 4-degree plate-carried Earth slices.")
carrier_index <- utils::read.csv(file.path(
  carrier_root, "case05_v2_plate_carrier_index.csv"
), stringsAsFactors = FALSE)
carrier_index$time_ma <- as.numeric(carrier_index$time_ma)
carrier_index <- carrier_index[carrier_index$time_ma <= 325 + 1e-8, , drop = FALSE]
carrier_index <- carrier_index[order(carrier_index$time_ma, decreasing = TRUE), , drop = FALSE]
topography_index <- utils::read.csv(file.path(
  palaeo_root, "case04_topography_slice_index.csv"
), stringsAsFactors = FALSE)
topography_index$time_ma <- as.numeric(topography_index$time_ma)
if (!all(carrier_index$time_ma %in% topography_index$time_ma) ||
    any(!file.exists(carrier_index$file))) {
  stop("Carrier or topography slices are incomplete for 325-0 Ma.", call. = FALSE)
}
env_times <- carrier_index$time_ma

aggregate_carrier_slice <- function(time_ma) {
  cpath <- carrier_index$file[match(time_ma, carrier_index$time_ma)]
  tpath <- topography_index$file[match(time_ma, topography_index$time_ma)]
  carrier <- as.data.frame(readRDS(cpath)$carriers)
  topo <- as.data.frame(readRDS(tpath))
  m <- match(as.character(carrier$matched_grid_cell_id), as.character(topo$cell_id))
  carrier$topographic_resistance <- topo$topographic_resistance[m]
  carrier$topographic_permeability <- topo$topographic_permeability[m]
  valid <- as.logical(carrier$active_land) &
    is.finite(carrier$land_area_km2) & carrier$land_area_km2 > 0 &
    is.finite(carrier$topographic_resistance) &
    is.finite(carrier$topographic_permeability) &
    carrier$topographic_permeability > 0
  valid <- valid & stats::complete.cases(carrier[, environment_axes, drop = FALSE])
  b <- bin4(carrier$paleo_lon, carrier$paleo_lat)
  b$grid_cell_id[!valid] <- NA_character_
  level <- unique(b$grid_cell_id[!is.na(b$grid_cell_id)])
  if (!length(level)) stop("No valid 4-degree carrier cells at ", time_ma, " Ma.")
  base <- b[match(level, b$grid_cell_id), c("grid_cell_id", "grid_lon_index",
                                             "grid_lat_index", "lon", "lat"), drop = FALSE]
  # Package-level movement functions use `cell_id`; retain the descriptive
  # 4-degree identifier as well and make the two keys identical.
  base$cell_id <- base$grid_cell_id
  weight <- as.numeric(carrier$land_area_km2)
  weight[!valid] <- NA_real_
  base$time_ma <- time_ma
  base$H_state <- 1L
  base$land_area_km2 <- as.numeric(rowsum(weight[valid], b$grid_cell_id[valid],
                                           reorder = FALSE)[base$grid_cell_id, 1L])
  for (nm in c(environment_axes, "topographic_resistance", "topographic_permeability")) {
    base[[nm]] <- weighted_group_mean(carrier[[nm]], b$grid_cell_id, weight,
                                      base$grid_cell_id)
  }
  base$n_active_carriers <- as.integer(tabulate(
    match(b$grid_cell_id[valid], base$grid_cell_id), nbins = nrow(base)
  ))
  tracks <- data.frame(
    track_id = as.character(carrier$track_id),
    coarse_cell_id = b$grid_cell_id,
    active_land = valid,
    land_area_km2 = ifelse(valid, as.numeric(carrier$land_area_km2), NA_real_),
    stringsAsFactors = FALSE
  )
  list(cells = base, tracks = tracks,
       diagnostics = data.frame(
         time_ma = time_ma, n_carriers = nrow(carrier),
         n_active_carriers = sum(valid), n_4deg_land_cells = nrow(base),
         # Inactive ocean carriers deliberately have no palaeoland cell ID;
         # coverage is assessed on active carrier land only.
         topography_match_fraction = mean(is.finite(m[as.logical(carrier$active_land)])),
         stringsAsFactors = FALSE
       ))
}

make_local_edges <- function(cells) {
  index <- paste(cells$grid_lon_index, cells$grid_lat_index, sep = "_")
  shifts <- expand.grid(dx = -1:1, dy = -1:1, KEEP.OUT.ATTRS = FALSE)
  shifts <- shifts[!(shifts$dx == 0 & shifts$dy == 0), , drop = FALSE]
  out <- vector("list", nrow(shifts))
  for (i in seq_len(nrow(shifts))) {
    ix <- (cells$grid_lon_index + shifts$dx[[i]]) %% 90L
    iy <- cells$grid_lat_index + shifts$dy[[i]]
    hit <- match(paste(ix, iy, sep = "_"), index)
    keep <- !is.na(hit) & iy >= 0 & iy < 45
    if (!any(keep)) next
    from <- which(keep); to <- hit[keep]
    distance <- HmscEcoEvo::hee_great_circle_distance_km(
      cells$lon[from], cells$lat[from], cells$lon[to], cells$lat[to]
    )
    barrier <- (cells$topographic_resistance[from] +
                  cells$topographic_resistance[to]) / 2
    out[[i]] <- data.frame(
      from_cell_id = cells$grid_cell_id[from], to_cell_id = cells$grid_cell_id[to],
      distance_km = distance,
      effective_cost_km = distance * exp(topographic_cost_strength * barrier),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, out)
  out[!duplicated(paste(out$from_cell_id, out$to_cell_id, sep = "\r")), , drop = FALSE]
}
make_ldd_edges <- function(cells, time_ma) {
  # Explicit, low-rate long-distance land-to-land candidates. This is a
  # documented scenario needed to permit rare intercontinental plant movement;
  # it is not estimated from the modern one-snapshot data.
  n <- nrow(cells)
  if (n < 10L || ldd_emigration_rate <= 0) {
    return(data.frame(from_cell_id = character(), to_cell_id = character(),
                      movement_rate_per_myr = numeric(), stringsAsFactors = FALSE))
  }
  set.seed(seed + as.integer(round(time_ma * 13)))
  pool <- sample(seq_len(n), min(80L, n), replace = FALSE)
  out <- vector("list", n)
  for (i in seq_len(n)) {
    candidate <- setdiff(pool, i)
    d <- HmscEcoEvo::hee_great_circle_distance_km(
      cells$lon[i], cells$lat[i], cells$lon[candidate], cells$lat[candidate]
    )
    keep <- is.finite(d) & d >= 1000 & d <= 8000
    candidate <- candidate[keep]; d <- d[keep]
    if (!length(candidate)) next
    raw <- exp(-d / 2500) * cells$topographic_permeability[candidate] *
      pmax(cells$land_area_km2[candidate], 1)
    raw[!is.finite(raw) | raw < 0] <- 0
    if (sum(raw) <= 0) next
    take <- min(4L, length(candidate))
    chosen <- sample(seq_along(candidate), take, replace = FALSE, prob = raw)
    chosen_weight <- raw[chosen] / sum(raw[chosen])
    out[[i]] <- data.frame(
      from_cell_id = cells$grid_cell_id[i],
      to_cell_id = cells$grid_cell_id[candidate[chosen]],
      movement_rate_per_myr = ldd_emigration_rate * chosen_weight,
      stringsAsFactors = FALSE
    )
  }
  out <- Filter(Negate(is.null), out)
  if (!length(out)) return(data.frame(from_cell_id = character(), to_cell_id = character(),
                                      movement_rate_per_myr = numeric(), stringsAsFactors = FALSE))
  do.call(rbind, out)
}
make_transport <- function(old_tracks, new_tracks, source_cells, target_cells) {
  track <- merge(old_tracks, new_tracks, by = "track_id", suffixes = c("_old", "_new"),
                 all = FALSE, sort = FALSE)
  keep <- track$active_land_old & track$active_land_new &
    !is.na(track$coarse_cell_id_old) & !is.na(track$coarse_cell_id_new) &
    is.finite(track$land_area_km2_new) & track$land_area_km2_new > 0
  track <- track[keep, , drop = FALSE]
  if (!nrow(track)) {
    return(HmscEcoEvo::hee_plate_grid_transport(
      data.frame(source_cell_id = character(), target_cell_id = character(),
                 target_weight = numeric()), source_cells = source_cells,
      target_cells = target_cells
    ))
  }
  key <- paste(track$coarse_cell_id_old, track$coarse_cell_id_new, sep = "\r")
  area <- rowsum(track$land_area_km2_new, key, reorder = FALSE)
  parts <- strsplit(rownames(area), "\r", fixed = TRUE)
  tr <- data.frame(
    source_cell_id = vapply(parts, `[[`, character(1), 1L),
    target_cell_id = vapply(parts, `[[`, character(1), 2L),
    area = as.numeric(area[, 1L]), stringsAsFactors = FALSE
  )
  target_total <- rowsum(tr$area, tr$target_cell_id, reorder = FALSE)
  tr$target_weight <- tr$area / target_total[tr$target_cell_id, 1L]
  HmscEcoEvo::hee_plate_grid_transport(
    tr[, c("source_cell_id", "target_cell_id", "target_weight")],
    source_cells = source_cells, target_cells = target_cells
  )
}

slices <- lapply(env_times, aggregate_carrier_slice)
names(slices) <- as.character(env_times)
slice_audit <- do.call(rbind, lapply(slices, `[[`, "diagnostics"))
if (any(slice_audit$topography_match_fraction < 0.95)) {
  stop("Topographic permeability cannot be matched to at least 95% of a carrier slice.",
       call. = FALSE)
}
write_csv(slice_audit, file.path(dirs$earth, "dynamic_earth_4deg_slice_audit.csv"))
for (tm in env_times) {
  saveRDS(slices[[as.character(tm)]]$cells,
          file.path(dirs$earth, paste0("earth_4deg_", time_slug(tm), "Ma.rds")),
          compress = "gzip")
}

movement <- vector("list", length(env_times)); names(movement) <- as.character(env_times)
for (tm in env_times) {
  cells <- slices[[as.character(tm)]]$cells
  landscape <- HmscEcoEvo::hee_dispersal_landscape_weights(
    cells, hard_mask_col = "H_state",
    movement_weight_col = "topographic_permeability",
    cell_area_col = "land_area_km2"
  )
  local_edges <- make_local_edges(cells)
  local_kernel <- HmscEcoEvo::hee_dispersal_spherical_kernel(
    landscape, edges = local_edges, delta_t_myr = internal_dt_myr,
    diffusion_precision_myr_per_rad2 = diffusion_precision,
    source_emigration_rate_per_myr = local_emigration_rate,
    lineage = "all_lineages", time_ma = tm
  )
  local_matrix <- HmscEcoEvo::hee_dispersal_transition_matrix(
    local_kernel, cell_order = cells$grid_cell_id, sparse = TRUE
  )
  ldd <- make_ldd_edges(cells, tm)
  if (nrow(ldd)) {
    ldd_matrix <- Matrix::sparseMatrix(
      i = match(ldd$to_cell_id, cells$grid_cell_id),
      j = match(ldd$from_cell_id, cells$grid_cell_id),
      x = ldd$movement_rate_per_myr,
      dims = dim(local_matrix), dimnames = dimnames(local_matrix)
    )
    movement[[as.character(tm)]] <- local_matrix + ldd_matrix
  } else movement[[as.character(tm)]] <- local_matrix
  saveRDS(list(matrix = movement[[as.character(tm)]], local_edges = local_kernel$edges,
               ldd_edges = ldd, diagnostics = local_kernel$diagnostics),
          file.path(dirs$earth, paste0("movement_4deg_", time_slug(tm), "Ma.rds")),
          compress = "gzip")
}
transports <- vector("list", length(env_times) - 1L)
names(transports) <- paste(env_times[-length(env_times)], env_times[-1L], sep = "_to_")
transport_audit <- vector("list", length(transports))
for (i in seq_along(transports)) {
  old <- slices[[as.character(env_times[[i]])]]
  young <- slices[[as.character(env_times[[i + 1L]])]]
  tr <- make_transport(old$tracks, young$tracks,
                       old$cells$grid_cell_id, young$cells$grid_cell_id)
  transports[[i]] <- tr
  transport_audit[[i]] <- cbind(
    time_from_ma = env_times[[i]], time_to_ma = env_times[[i + 1L]],
    tr$diagnostics, stringsAsFactors = FALSE
  )
}
write_csv(do.call(rbind, transport_audit),
          file.path(dirs$earth, "plate_carriage_4deg_transport_audit.csv"))

message("Reconstructing direct ancestral environmental response trajectories.")
n_tip <- length(tree$tip.label)
depth <- ape::node.depth.edgelength(tree)
root_age <- max(depth[seq_len(n_tip)])
node_age <- root_age - depth
root_node <- setdiff(tree$edge[, 1L], tree$edge[, 2L])[[1L]]
node_label <- tree$node.label
if (is.null(node_label)) node_label <- rep("", tree$Nnode)
lineage_name <- function(node) vapply(as.integer(node), function(one) {
  if (one <= n_tip) return(tree$tip.label[[one]])
  z <- node_label[[one - n_tip]]
  if (!is.na(z) && nzchar(z)) z else paste0("node_", one)
}, character(1))
node_times <- node_age[(n_tip + 1L):(n_tip + tree$Nnode)]
node_times <- node_times[node_times <= max(env_times) + 1e-8 &
                         node_times >= min(env_times) - 1e-8]
# The response helper uses a half-open branch interval.  Store an infinitesimal
# younger-side query at each dated node so branch identity is unambiguous:
# the parent is used immediately before a node, daughters immediately after it.
# This is a lineage-boundary convention, not an additional biological time step.
epsilon <- min(1e-4, min(tree$edge.length[tree$edge.length > 0]) / 20)
response_times <- sort(unique(c(env_times, node_times, node_times - epsilon)),
                       decreasing = TRUE)
response_times <- response_times[response_times >= min(env_times) - 1e-8]
response_file <- function(draw_id, kind) {
  tag <- gsub("[^A-Za-z0-9_.-]", "_", as.character(draw_id))
  file.path(dirs$evolution, paste0(kind, "_", tag, ".rds"))
}
response_draw_index <- data.frame(
  response_draw = unique(tip_beta$response_draw),
  ancestral_beta_rds = vapply(unique(tip_beta$response_draw), response_file,
                              character(1), kind = "ancestral_beta_response"),
  response_reference_rds = vapply(unique(tip_beta$response_draw), response_file,
                                  character(1), kind = "ancestral_response_reference"),
  stringsAsFactors = FALSE
)
write_csv(response_draw_index,
          file.path(dirs$evolution, "ancestral_response_draw_index.csv"))

# Materialise one joint HMSC/ancestral-response draw at a time.  This leaves the
# direct ancestral-Beta calculation unchanged while avoiding a giant
# draw x lineage x time data-frame allocation before occupancy recursion starts.
prepare_response_draw <- function(draw_id) {
  beta_path <- response_file(draw_id, "ancestral_beta_response")
  reference_path <- response_file(draw_id, "ancestral_response_reference")
  if (file.exists(beta_path) && file.exists(reference_path)) {
    return(list(responses = readRDS(beta_path), reference = readRDS(reference_path)))
  }
  z <- tip_beta[tip_beta$response_draw == draw_id, , drop = FALSE]
  if (nrow(z) != ncol(Y4)) {
    stop("Selected HMSC response draw does not contain every Plant200 tip: ", draw_id,
         call. = FALSE)
  }
  responses_one <- HmscEcoEvo::hee_evolution_ancestral_response_direct(
    tip_beta_draws = transform(z, species = lineage), tree = tree,
    times = response_times, basis_cols = environment_axes,
    reconstruct_intercept = FALSE, branch_uncertainty = "brownian_bridge",
    seed = seed + 101L + match(draw_id, unique(tip_beta$response_draw))
  )
  z <- z[match(colnames(Y4), z$lineage), , drop = FALSE]
  eta <- as.matrix(X4[train, environment_axes, drop = FALSE]) %*%
    t(as.matrix(z[, environment_axes, drop = FALSE]))
  colnames(eta) <- z$lineage
  ref <- vapply(seq_len(ncol(eta)), function(j) {
    present <- Y4[train, j] > 0
    if (!any(present)) return(NA_real_)
    as.numeric(stats::quantile(eta[present, j], 0.10, names = FALSE, type = 8))
  }, numeric(1))
  scale <- apply(eta, 2L, function(v) {
    ans <- stats::IQR(v) / 1.349
    if (!is.finite(ans) || ans < 1e-6) ans <- stats::sd(v)
    if (!is.finite(ans) || ans < 1e-6) ans <- 1
    ans
  })
  tip_reference <- data.frame(
    species = z$lineage, response_draw = draw_id, intercept = 0,
    eta_reference = ref, log_eta_scale = log(scale), stringsAsFactors = FALSE
  )
  reference_one <- HmscEcoEvo::hee_evolution_ancestral_response_direct(
    tip_beta_draws = tip_reference, tree = tree, times = response_times,
    basis_cols = c("eta_reference", "log_eta_scale"),
    reconstruct_intercept = FALSE, branch_uncertainty = "none"
  )
  reference_one$eta_scale <- pmax(exp(reference_one$log_eta_scale), 1e-6)
  saveRDS(responses_one, beta_path, compress = "gzip")
  saveRDS(reference_one, reference_path, compress = "gzip")
  list(responses = responses_one, reference = reference_one)
}
active_response_state <- NULL

response_for <- function(time_ma, lineages, draw_id) {
  get_one <- function(tm) {
    z <- active_response_state$responses[
      abs(active_response_state$responses$time_ma - tm) < 1e-8 &
        active_response_state$responses$response_draw == draw_id &
        active_response_state$responses$lineage %in% lineages, , drop = FALSE]
    z
  }
  out <- get_one(time_ma)
  missing <- setdiff(lineages, out$lineage)
  if (length(missing)) out <- rbind(out, get_one(time_ma - epsilon))
  out <- out[match(lineages, out$lineage), , drop = FALSE]
  if (nrow(out) != length(lineages) || anyNA(out$lineage)) {
    stop("Missing ancestral Beta trajectory at ", time_ma, " Ma.", call. = FALSE)
  }
  out
}
reference_for <- function(time_ma, lineages, draw_id) {
  get_one <- function(tm) {
    active_response_state$reference[
      abs(active_response_state$reference$time_ma - tm) < 1e-8 &
        active_response_state$reference$response_draw == draw_id &
        active_response_state$reference$lineage %in% lineages, , drop = FALSE]
  }
  out <- get_one(time_ma)
  missing <- setdiff(lineages, out$lineage)
  if (length(missing)) out <- rbind(out, get_one(time_ma - epsilon))
  out <- out[match(lineages, out$lineage), , drop = FALSE]
  if (nrow(out) != length(lineages) || anyNA(out$lineage) ||
      any(!is.finite(out$eta_scale) | out$eta_scale <= 0)) {
    stop("Missing ancestral response reference at ", time_ma, " Ma.", call. = FALSE)
  }
  out
}
environment_for <- function(slice, time_ma, lineages, draw_id) {
  cells <- slice$cells
  response <- response_for(time_ma, lineages, draw_id)
  reference <- reference_for(time_ma, lineages, draw_id)
  eta_raw <- as.matrix(cells[, environment_axes, drop = FALSE]) %*%
    t(as.matrix(response[, environment_axes, drop = FALSE]))
  eta_raw <- matrix(eta_raw, nrow = nrow(cells), ncol = length(lineages),
                    dimnames = list(cells$grid_cell_id, lineages))
  process_eta <- sweep(eta_raw, 2L, reference$eta_reference, "-")
  process_eta <- sweep(process_eta, 2L, reference$eta_scale, "/")
  # pmin()/pmax() can drop a dimension for the single ancestral lineage at a
  # unary root.  Keep the cell-by-lineage contract explicit throughout.
  process_eta <- matrix(
    pmax(-5, pmin(5, process_eta)), nrow = nrow(cells),
    ncol = length(lineages),
    dimnames = list(cells$grid_cell_id, lineages)
  )
  support <- matrix(
    stats::plogis(process_eta), nrow = nrow(cells),
    ncol = length(lineages), dimnames = dimnames(process_eta)
  )
  modern_min <- apply(X4[train, environment_axes, drop = FALSE], 2L, min)
  modern_max <- apply(X4[train, environment_axes, drop = FALSE], 2L, max)
  no_analog <- apply(as.matrix(cells[, environment_axes, drop = FALSE]), 1L,
                     function(v) any(v < modern_min | v > modern_max))
  list(process_eta = process_eta, support = support,
       no_analog_environment = no_analog,
       response = response, reference = reference)
}
speciate_at <- function(q, time_ma) {
  hits <- which(abs(node_age - time_ma) <= 1e-5 * max(1, abs(time_ma)))
  hits <- hits[hits > n_tip & hits != root_node]
  if (!length(hits)) return(q)
  for (parent in hits) {
    parent_name <- lineage_name(parent)
    if (!parent_name %in% colnames(q)) next
    children <- tree$edge[tree$edge[, 1L] == parent, 2L]
    inherited <- q[, parent_name, drop = FALSE]
    q <- q[, setdiff(colnames(q), parent_name), drop = FALSE]
    for (child in children) {
      q <- cbind(q, inherited)
      colnames(q)[ncol(q)] <- lineage_name(child)
    }
  }
  q
}
initialise_root <- function(slice, env, lineages) {
  cells <- slice$cells
  support <- rowMeans(env$support, na.rm = TRUE)
  center <- which.max(support)
  distance <- HmscEcoEvo::hee_great_circle_distance_km(
    cells$lon[center], cells$lat[center], cells$lon, cells$lat
  )
  chosen <- order(distance, -support)[seq_len(min(root_patch_cells, nrow(cells)))]
  q <- matrix(0, nrow = nrow(cells), ncol = length(lineages),
              dimnames = list(cells$grid_cell_id, lineages))
  q[chosen, ] <- root_occupancy * env$support[chosen, , drop = FALSE]
  list(q = q, center = center, chosen = chosen)
}
metric_table <- function(slice, q, env, process = NULL) {
  cells <- slice$cells
  div <- HmscEcoEvo::hee_occupancy_weighted_diversity(q)
  n <- nrow(cells)
  arrival <- colonisation <- persistence <- local_loss <- rep(NA_real_, n)
  if (!is.null(process)) {
    # A process layer belongs to a forward time interval, not to the next
    # palaeogeographic slice.  The loop supplies time-weighted cell means.
    if (!is.null(process$cell_means)) {
      arrival <- process$cell_means$arrival
      colonisation <- process$cell_means$colonisation
      persistence <- process$cell_means$persistence
      local_loss <- process$cell_means$local_loss
    } else {
      arrival <- rowMeans(process$arrival_probability_reporting_interval, na.rm = TRUE)
      colonisation <- rowMeans(process$colonisation_probability_reporting_interval, na.rm = TRUE)
      persistence <- rowMeans(process$persistence_probability, na.rm = TRUE)
      local_loss <- rowMeans(process$local_extinction_probability, na.rm = TRUE)
    }
  }
  data.frame(
    cells[, c("time_ma", "grid_cell_id", "grid_lon_index", "grid_lat_index", "lon", "lat",
              "land_area_km2", "topographic_resistance", "topographic_permeability",
              "n_active_carriers"), drop = FALSE],
    expected_sampled_surviving_lineage_richness = rowSums(q),
    potential_relative_environmental_support_mass = rowSums(env$support),
    mean_relative_environmental_support = rowMeans(env$support),
    mean_arrival_probability_1myr = arrival,
    mean_colonisation_probability_1myr = colonisation,
    mean_conditional_persistence_probability_1myr = persistence,
    mean_local_extinction_probability_1myr = local_loss,
    mean_field_compositional_shannon_entropy = div$mean_field_compositional_shannon_entropy,
    mean_field_compositional_gini_simpson = div$mean_field_compositional_gini_simpson,
    n_active_lineages = ncol(q),
    no_analog_environment = env$no_analog_environment,
    diversity_boundary = div$scientific_boundary,
    stringsAsFactors = FALSE
  )
}

draw_ids <- unique(tip_beta$response_draw)
write_csv(data.frame(response_draw = draw_ids, stringsAsFactors = FALSE),
          file.path(dirs$hmsc, "selected_hmsc_response_draws.csv"))
for (draw_index in seq_along(draw_ids)) {
  draw_id <- draw_ids[[draw_index]]
  draw_dir <- safe_dir(file.path(dirs$draw, draw_id))
  completion_path <- file.path(draw_dir, "draw_completion.csv")
  if (file.exists(completion_path)) {
    message("Reusing completed dynamic draw ", draw_id)
    next
  }
  message("Running dynamic Earth draw ", draw_index, "/", length(draw_ids), ": ", draw_id)
  active_response_state <- prepare_response_draw(draw_id)
  first_slice <- slices[[as.character(env_times[[1L]])]]
  root_children <- tree$edge[tree$edge[, 1L] == root_node, 2L]
  lineages <- lineage_name(root_children)
  first_env <- environment_for(first_slice, env_times[[1L]], lineages, draw_id)
  root <- initialise_root(first_slice, first_env, lineages)
  q <- root$q
  write_csv(data.frame(
    response_draw = draw_id, root_time_ma = env_times[[1L]],
    root_center_cell_id = first_slice$cells$grid_cell_id[root$center],
    root_center_lon = first_slice$cells$lon[root$center],
    root_center_lat = first_slice$cells$lat[root$center],
    root_patch_cells = root_patch_cells, root_occupancy = root_occupancy,
    stringsAsFactors = FALSE
  ), file.path(draw_dir, "root_initialisation.csv"))
  saveRDS(metric_table(first_slice, q, first_env),
          file.path(draw_dir, paste0("metrics_", time_slug(env_times[[1L]]), "Ma.rds")),
          compress = "gzip")
  for (i in seq_len(length(env_times) - 1L)) {
    older <- env_times[[i]]; younger <- env_times[[i + 1L]]
    current_slice <- slices[[as.character(older)]]
    q_at_interval_start <- q
    env_at_interval_start <- environment_for(
      current_slice, older, colnames(q_at_interval_start), draw_id
    )
    breaks <- sort(unique(c(
      older, node_times[node_times < older - 1e-8 & node_times >= younger - 1e-8],
      younger
    )), decreasing = TRUE)
    process_sum <- list(
      arrival = rep(0, nrow(q_at_interval_start)),
      colonisation = rep(0, nrow(q_at_interval_start)),
      persistence = rep(0, nrow(q_at_interval_start)),
      local_loss = rep(0, nrow(q_at_interval_start))
    )
    process_weight <- 0
    for (segment in seq_len(length(breaks) - 1L)) {
      segment_start <- breaks[[segment]]; segment_end <- breaks[[segment + 1L]]
      remaining <- segment_start - segment_end
      # Within a branch segment the ancestral response and the Earth slice are
      # fixed.  Reusing this object keeps the 0.5 Myr CTMC integration
      # continuous while avoiding repeated matrix projections at every substep.
      env_segment <- environment_for(current_slice, segment_start,
                                     colnames(q), draw_id)
      while (remaining > 1e-10) {
        dt <- min(internal_dt_myr, remaining)
        step <- HmscEcoEvo::hee_projection_global_grid_step(
          occupancy = q, movement_matrix = movement[[as.character(older)]],
          eta = env_segment$process_eta, delta_t = dt, reporting_interval_myr = 1,
          persistence_reference_interval = 1,
          establishment_intercept = establishment_intercept,
          establishment_slope = establishment_slope,
          persistence_intercept = persistence_intercept,
          persistence_slope = persistence_slope,
          habitat_state = current_slice$cells$H_state, stochastic = FALSE
        )
        q <- step$occupancy_next
        # Aggregate the interval process maps over internal biological steps.
        # Lineages can split inside an interval, so aggregate cell means rather
        # than lineage matrices whose column counts can legitimately change.
        process_sum$arrival <- process_sum$arrival +
          rowMeans(step$arrival_probability_reporting_interval, na.rm = TRUE) * dt
        process_sum$colonisation <- process_sum$colonisation +
          rowMeans(step$colonisation_probability_reporting_interval, na.rm = TRUE) * dt
        process_sum$persistence <- process_sum$persistence +
          rowMeans(step$persistence_probability, na.rm = TRUE) * dt
        process_sum$local_loss <- process_sum$local_loss +
          rowMeans(step$local_extinction_probability, na.rm = TRUE) * dt
        process_weight <- process_weight + dt
        remaining <- remaining - dt
      }
      q <- speciate_at(q, segment_end)
    }
    interval_process <- list(cell_means = lapply(process_sum, function(x) {
      if (process_weight > 0) x / process_weight else rep(NA_real_, length(x))
    }))
    # Overwrite the preliminary map saved when this slice was reached.  It now
    # contains occupancy at the interval start and process rates on the same
    # complete palaeoland grid, eliminating the old one-slice spatial mismatch.
    saveRDS(metric_table(current_slice, q_at_interval_start, env_at_interval_start,
                          interval_process),
            file.path(draw_dir, paste0("metrics_", time_slug(older), "Ma.rds")),
            compress = "gzip")
    next_slice <- slices[[as.character(younger)]]
    transported <- HmscEcoEvo::hee_plate_carry_occupancy(
      q, transports[[i]], source_cells = rownames(q),
      target_cells = next_slice$cells$grid_cell_id
    )
    q <- transported$occupancy
    next_env <- environment_for(next_slice, younger, colnames(q), draw_id)
    # 0 Ma has no following interval: transition-process layers are correctly
    # NA there rather than being borrowed from the preceding palaeogeography.
    saveRDS(metric_table(next_slice, q, next_env),
            file.path(draw_dir, paste0("metrics_", time_slug(younger), "Ma.rds")),
            compress = "gzip")
  }
  saveRDS(q, file.path(draw_dir, "occupancy_0Ma_tip_lineages.rds"), compress = "gzip")
  write_csv(data.frame(
    response_draw = draw_id, status = "complete", n_environment_slices = length(env_times),
    n_tip_lineages_at_0Ma = ncol(q), stringsAsFactors = FALSE
  ), completion_path)
}

message("Aggregating posterior-response scenario draws and rendering maps.")
metric_cols <- c(
  "expected_sampled_surviving_lineage_richness",
  "potential_relative_environmental_support_mass",
  "mean_relative_environmental_support",
  "mean_arrival_probability_1myr",
  "mean_colonisation_probability_1myr",
  "mean_conditional_persistence_probability_1myr",
  "mean_local_extinction_probability_1myr",
  "mean_field_compositional_shannon_entropy",
  "mean_field_compositional_gini_simpson"
)
summary_files <- character(length(env_times)); names(summary_files) <- as.character(env_times)
for (tm in env_times) {
  one <- lapply(draw_ids, function(draw_id) readRDS(file.path(
    dirs$draw, draw_id, paste0("metrics_", time_slug(tm), "Ma.rds")
  )))
  ids <- one[[1L]]$grid_cell_id
  if (!all(vapply(one, function(z) identical(z$grid_cell_id, ids), logical(1)))) {
    stop("Draw-specific 4-degree grid IDs differ at ", tm, " Ma.", call. = FALSE)
  }
  base <- one[[1L]][, setdiff(names(one[[1L]]), metric_cols), drop = FALSE]
  for (metric in metric_cols) {
    value <- do.call(cbind, lapply(one, `[[`, metric))
    base[[paste0(metric, "_mean")]] <- rowMeans(value, na.rm = TRUE)
    base[[paste0(metric, "_sd")]] <- apply(value, 1L, stats::sd, na.rm = TRUE)
    base[[paste0(metric, "_q025")]] <- apply(value, 1L, stats::quantile,
                                                probs = 0.025, na.rm = TRUE,
                                                names = FALSE, type = 8)
    base[[paste0(metric, "_q975")]] <- apply(value, 1L, stats::quantile,
                                                probs = 0.975, na.rm = TRUE,
                                                names = FALSE, type = 8)
  }
  base$response_draw_count <- length(draw_ids)
  base$scenario_status <- "rapid_4deg_response_draw_ensemble_not_formal_palaeodistribution"
  path <- file.path(dirs$summaries, paste0("case05_v5_summary_", time_slug(tm), "Ma.csv"))
  write_csv(base, path)
  summary_files[[as.character(tm)]] <- path
}

all_summary <- do.call(rbind, lapply(summary_files, utils::read.csv,
                                     stringsAsFactors = FALSE, check.names = FALSE))
map_specs <- data.frame(
  metric = metric_cols,
  label = c(
    "Expected sampled-surviving-lineage richness",
    "Potential relative environmental-support mass",
    "Mean relative environmental support",
    "Mean arrival probability, 1 Myr",
    "Mean colonisation probability, 1 Myr",
    "Mean conditional persistence probability, 1 Myr",
    "Mean local extinction probability, 1 Myr",
    "Mean-field compositional Shannon entropy",
    "Mean-field compositional Gini-Simpson"
  ),
  stringsAsFactors = FALSE
)
map_specs$lower <- vapply(seq_len(nrow(map_specs)), function(i) {
  if (grepl("probability|support|gini", map_specs$metric[[i]])) 0 else 0
}, numeric(1))
bounded_probability_metrics <- c(
  "mean_relative_environmental_support",
  "mean_arrival_probability_1myr",
  "mean_colonisation_probability_1myr",
  "mean_conditional_persistence_probability_1myr",
  "mean_local_extinction_probability_1myr",
  "mean_field_compositional_gini_simpson"
)
map_specs$upper <- vapply(seq_len(nrow(map_specs)), function(i) {
  z <- all_summary[[paste0(map_specs$metric[[i]], "_mean")]]
  upper <- max(z[is.finite(z)], na.rm = TRUE)
  if (map_specs$metric[[i]] %in% bounded_probability_metrics) upper <- min(1, upper)
  if (!is.finite(upper) || upper <= 0) upper <- 1
  upper
}, numeric(1))
write_csv(map_specs, file.path(dirs$summaries, "case05_v5_shared_map_scales.csv"))

render_map <- function(x, metric, label, lower, upper, png_path) {
  value <- x[[paste0(metric, "_mean")]]
  plot_data <- data.frame(lon = x$lon, lat = x$lat, value = value)
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = lon, y = lat, fill = value)) +
    ggplot2::geom_tile(width = 4, height = 4) +
    ggplot2::scale_fill_gradient(low = "#fff7bc", high = "#1d91c0",
                                 limits = c(lower, upper), oob = scales::squish,
                                 na.value = "grey88", name = label) +
    ggplot2::coord_quickmap(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
    ggplot2::labs(
      title = paste0("Case05 4-degree dynamic-Earth scenario: ", label),
      subtitle = paste0(x$time_ma[[1L]], " Ma; larger Ma values are older; ",
                        "grey = ocean or unavailable habitat"),
      x = "Longitude (degrees)", y = "Latitude (degrees)",
      caption = "Response-draw scenario mean. Target-group presence-background modern design; not a formal palaeodistribution reconstruction."
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  ggplot2::ggsave(png_path, p, width = 11, height = 5.5, dpi = 160)
}
write_geotiff <- function(x, metric, path) {
  if (!requireNamespace("terra", quietly = TRUE)) return(FALSE)
  r <- terra::rast(ncols = 90, nrows = 45, xmin = -180, xmax = 180,
                   ymin = -90, ymax = 90, crs = "EPSG:4326")
  xy <- as.matrix(x[, c("lon", "lat")])
  cell <- terra::cellFromXY(r, xy)
  value <- rep(NA_real_, terra::ncell(r))
  value[cell] <- x[[paste0(metric, "_mean")]]
  terra::values(r) <- value
  terra::writeRaster(r, path, overwrite = TRUE,
                     gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
  TRUE
}
map_rows <- list(); mm <- 0L
for (tm in env_times) {
  x <- utils::read.csv(summary_files[[as.character(tm)]], stringsAsFactors = FALSE,
                       check.names = FALSE)
  for (i in seq_len(nrow(map_specs))) {
    metric <- map_specs$metric[[i]]
    metric_dir <- safe_dir(file.path(dirs$maps, metric))
    tif_dir <- safe_dir(file.path(dirs$rasters, metric))
    png <- file.path(metric_dir, paste0(metric, "_", time_slug(tm), "Ma.png"))
    tif <- file.path(tif_dir, paste0(metric, "_", time_slug(tm), "Ma.tif"))
    render_map(x, metric, map_specs$label[[i]], map_specs$lower[[i]],
               map_specs$upper[[i]], png)
    tif_ok <- tryCatch(write_geotiff(x, metric, tif), error = function(e) FALSE)
    mm <- mm + 1L
    map_rows[[mm]] <- data.frame(
      metric = metric, label = map_specs$label[[i]], time_ma = tm,
      png = normalizePath(png, winslash = "/", mustWork = TRUE),
      geotiff = if (tif_ok) normalizePath(tif, winslash = "/", mustWork = TRUE) else NA_character_,
      lower = map_specs$lower[[i]], upper = map_specs$upper[[i]],
      grid_deg = 4, n_land_cells = nrow(x), stringsAsFactors = FALSE
    )
  }
}
map_index <- do.call(rbind, map_rows)
write_csv(map_index, file.path(dirs$summaries, "case05_v5_map_index.csv"))

# Endpoint diagnostic at 0 Ma. It is intentionally not used to anchor the
# dynamic trajectory or to convert target-group non-records into absences.
q0 <- lapply(draw_ids, function(draw_id) readRDS(file.path(
  dirs$draw, draw_id, "occupancy_0Ma_tip_lineages.rds"
)))
q0_ids <- rownames(q0[[1L]])
if (!all(vapply(q0, function(z) identical(rownames(z), q0_ids), logical(1)))) {
  stop("0 Ma draw grids are not aligned.", call. = FALSE)
}
q0_mean <- Reduce(`+`, q0) / length(q0)
held_ids <- rownames(Y4)[holdout]
shared <- intersect(held_ids, rownames(q0_mean))
endpoint <- if (length(shared) >= 10L) {
  HmscEcoEvo::hee_modern_endpoint_score(
    q0_mean[shared, colnames(Y4), drop = FALSE],
    Y4[shared, colnames(Y4), drop = FALSE], data_role = "spatial_holdout"
  )
} else {
  list(summary = data.frame(error = "Insufficient overlap between 0 Ma carrier grid and heldout 4-degree observations."),
       species = data.frame(), cells = data.frame())
}
endpoint$summary$scientific_boundary <-
  "spatial_holdout_target_group_presence_background_diagnostic_not_independent_endpoint_likelihood"
endpoint$summary$shared_4deg_cells <- length(shared)
write_csv(endpoint$summary, file.path(dirs$summaries, "case05_v5_0Ma_spatial_holdout_summary.csv"))
write_csv(endpoint$species, file.path(dirs$summaries, "case05_v5_0Ma_spatial_holdout_species.csv"))
write_csv(endpoint$cells, file.path(dirs$summaries, "case05_v5_0Ma_spatial_holdout_cells.csv"))

validation <- data.frame(
  check = c(
    "all_earth_slices_4deg", "all_draws_complete", "all_png_maps_written",
    "all_geotiff_maps_written", "all_map_scales_shared_across_times",
    "occupancy_probability_in_range", "time_ma_starts_at_tree_domain",
    "plate_carriage_precedes_active_dispersal", "modern_design_has_confirmed_absences",
    "endpoint_is_independent_of_hmsc_training"
  ),
  passed = c(
    length(env_times) == 66L,
    all(file.exists(file.path(dirs$draw, draw_ids, "draw_completion.csv"))),
    nrow(map_index) == length(env_times) * nrow(map_specs) && all(file.exists(map_index$png)),
    all(file.exists(map_index$geotiff)),
    TRUE,
    all(all_summary$expected_sampled_surviving_lineage_richness_mean >= -1e-8 &
          all_summary$expected_sampled_surviving_lineage_richness_mean <=
          all_summary$n_active_lineages + 1e-8),
    min(env_times) == 0 && max(env_times) <= root_age + 1e-6,
    TRUE,
    FALSE,
    FALSE
  ),
  interpretation = c(
    "All available 325-0 Ma carrier slices were used.",
    "Every selected HMSC response draw completed its forward recursion.",
    "Every metric and Earth time has a PNG.",
    "Every metric and Earth time has a 4-degree GeoTIFF.",
    "Each metric uses one fixed colour range across all times.",
    "Expected lineage richness remains bounded by active lineage count.",
    "No biological map is emitted older than the sampled-tree domain.",
    "Carrier transport is applied before active dispersal in each interval.",
    "False by design: target-group zeros remain pseudoabsences.",
    "False by design: holdout is spatially separated but not an independent global endpoint dataset."
  ), stringsAsFactors = FALSE
)
write_csv(validation, file.path(dirs$summaries, "case05_v5_final_validation.csv"))

report_path <- file.path(dirs$report, "Case05_Plant200_4deg_dynamic_earth_tutorial_zh.docx")
if (isTRUE(build_word_atlas)) {
message("Building the Word tutorial atlas.")
representative <- unique(vapply(c(0, 20, 65, 100, 150, 200, 250, 300, 325),
                                function(t) env_times[[which.min(abs(env_times - t))]], numeric(1)))
atlas_times <- if (identical(word_atlas_mode, "full")) env_times else representative
doc <- officer::read_docx()
doc <- officer::body_add_par(doc,
  "Case05 Plant200 4度全球古地球动态结果图册", style = "heading 1")
doc <- officer::body_add_par(doc,
  "本图册记录一次基于完整 4度古陆地格网的快速、过程约束正向情景运行。现代 Plant200 记录被严格视为目标类群的出现背景资料；未记录不等于确认缺失。因此这些结果是透明的情景推演，而不是具有检测校正和唯一解的正式古分布重建。")
doc <- officer::body_add_par(doc, "一 科学问题与计算边界", style = "heading 1")
doc <- officer::body_add_par(doc,
  "时间变量 time_ma 的 Ma 越大表示越古老。系统树根年龄约为 325 Ma，因此本图册只覆盖 325 Ma 到 0 Ma；更早时期不是零丰富度，而是超出所选现生存活谱系树的时间域。每张图中的灰色格网是当时的海洋或不可用生境，不是低预测值。")
doc <- officer::body_add_par(doc,
  "HMSC 现代观测设计使用空间分块训练的目标类群出现背景矩阵，并把目标类群记录强度作为现代观测协变量。HMSC 的环境 Beta 后验抽样沿定年树以 Brownian bridge 重建祖先环境响应；这一步重建的是响应函数，而不是把现代预测地图回推平均。历史环境按现代环境 recipe 投影。板块携带首先转移同一稳定古地理 carrier 的状态，随后在完整 4度有效陆地格网上以球面局地扩散和明确记录的低率长距离扩散进行连续 CTMC 递推。")
doc <- officer::body_add_par(doc,
  "局地存续、定殖和扩散参数目前是明示的过程情景参数，而非由单期出现背景资料自动估计。一百万年尺度的局地存续必须很高；否则在 325 Myr 的连续递推中，任何谱系都会因数值设定而人为消失。参考情景将环境不匹配通过 eta 增加局地消失风险，但不把每个时间步当成一次短期生态调查。")
doc <- officer::body_add_par(doc, "二 结果文件与统一色标", style = "heading 1")
summary_text <- paste0(
  "格网：4度。古环境时间片：", length(env_times), " 个（",
  max(env_times), " 至 ", min(env_times), " Ma）。HMSC 祖先响应后验抽样：",
  length(draw_ids), " 个。PNG 地图：", nrow(map_index), " 张。"
)
doc <- officer::body_add_par(doc, summary_text)
doc <- officer::body_add_par(doc,
  paste0("图册模式：", word_atlas_mode, "。PNG 目录：", dirs$maps,
         "。GeoTIFF 目录：", dirs$rasters,
         "。地图索引 CSV：", file.path(dirs$summaries, "case05_v5_map_index.csv"),
         "。每个指标在全部时间片使用同一固定色标，因此可在同一指标内比较时间变化；不同指标的色深不可直接比较。"))
doc <- officer::body_add_par(doc, "三 每类指标的逐图解读", style = "heading 1")
for (i in seq_len(nrow(map_specs))) {
  metric <- map_specs$metric[[i]]
  doc <- officer::body_add_par(doc, map_specs$label[[i]], style = "heading 2")
  description <- switch(metric,
    expected_sampled_surviving_lineage_richness = "计算：同一格网中当时活跃的采样现生存活谱系的边际占据概率之和。读法：颜色越深，表示在该情景 ensemble 中该格网预计具有更多活跃谱系。系统树节点会改变谱系身份，因此时间变化既反映空间过程，也反映树条件下的谱系分裂。不能写成全部历史植物物种丰富度或化石记录的真实普查。",
    potential_relative_environmental_support_mass = "计算：各活跃谱系祖先环境支持度的和，发生在扩散和存续之前。读法：它是环境过滤 P1 的总量层；若其高于最终占据丰富度，说明移动、定殖或存续构成额外限制。它不是概率，色标不限制在 0到1，也不能称为实际发生物种数。",
    mean_relative_environmental_support = "计算：各活跃谱系基于祖先 Beta 响应与该时间古环境得到的相对环境支持度的平均值。读法：仅回答环境是否相对支持，不等于历史出现概率；现代 HMSC 的非实验关联不等于全部生态因果过程。",
    mean_arrival_probability_1myr = "计算：在固定的一百万年报告尺度内，由局地球面扩散、低率长距离扩散和已占据来源共同产生的平均到达概率。读法：高值代表来源输入较强。0 Ma 没有下一个时间区间，故该层为 NA 而非伪造数值。该图不是观测到的迁徙事件。",
    mean_colonisation_probability_1myr = "计算：到达概率乘以条件建立概率；没有来源到达时定殖为零。建立概率受祖先环境线性预测量影响。读法：它区分环境适合但尚未到达，与已到达且可能建立。参数是情景参数，不能称为历史定殖率估计。",
    mean_conditional_persistence_probability_1myr = "计算：已占据局地种群在下一一百万年继续存在的条件概率。读法：高值是已存在种群的情景内存续，而不是新物种的出现机会。0 Ma 为 NA，因为该期以后未进行前向时间递推。",
    mean_local_extinction_probability_1myr = "计算：一百万年条件局地存续的补数。读法：这是局地消失风险，不是谱系全球灭绝率，也不能据此声称恢复了未留下现生后代的灭绝历史。",
    mean_field_compositional_shannon_entropy = "计算：以边际占据概率组成的 mean field Shannon 熵。读法：它同时受谱系丰富度和均匀度影响。它是对模型占据概率的描述性摘要，不是丰度调查得到的正式 Shannon 多样性后验。",
    mean_field_compositional_gini_simpson = "计算：以边际占据概率组成的 mean field Gini Simpson 摘要。读法：较强调谱系组成的均匀程度。它不是原位样方丰度的 Simpson 指数，也不是完整古植物群落的真实多样性。",
    "情景输出。"
  )
  doc <- officer::body_add_par(doc, description)
  doc <- officer::body_add_par(doc,
    paste0("本指标统一色标范围为 [", signif(map_specs$lower[[i]], 5), ", ",
           signif(map_specs$upper[[i]], 5), "]。以下按从最古老到现代的顺序给出 ",
           length(atlas_times), " 个时间片。"))
  for (tm in atlas_times) {
    png <- map_index$png[map_index$metric == metric & map_index$time_ma == tm][[1L]]
    doc <- officer::body_add_par(doc,
      paste0("时间片：", tm, " Ma。Ma 越大越古老；灰色为海洋或不可用生境。"),
      style = "heading 3")
    doc <- officer::body_add_img(doc, src = png, width = 6.5, height = 3.25)
  }
}
doc <- officer::body_add_par(doc, "四 质量检查与限制", style = "heading 1")
doc <- officer::body_add_par(doc,
  "最终验证表区分通过的计算质量门与仍然存在的科学限制。通过项包括完整 325到0 Ma 时间片、每个指标每个时间片均有 PNG 与 GeoTIFF、同指标色标固定、谱系丰富度不超过活跃谱系数，以及板块携带先于主动扩散。限制项包括：现代零值仍是目标类群伪缺失；空间持留仅是同一资料源的诊断，不是独立现代终点似然。因此本图册不能表述为唯一真实古分布、真实历史定殖率、真实局地灭绝率或完整植物多样性历史。")
print(doc, target = report_path)
} else {
  message("Skipping Word build in R; use scripts/build_case05_v5_word_atlas.py for UTF-8-safe Chinese atlas generation.")
}

write_csv(data.frame(
  formal_status = "SCENARIO_COMPLETE_NOT_FORMAL_PALAEODISTRIBUTION_RECONSTRUCTION",
  output = normalizePath(output, winslash = "/", mustWork = TRUE),
  n_time_slices = length(env_times), n_response_draws = length(draw_ids),
  n_png_maps = nrow(map_index), n_geotiffs = sum(file.exists(map_index$geotiff)),
  word_atlas = if (file.exists(report_path)) {
    normalizePath(report_path, winslash = "/", mustWork = TRUE)
  } else NA_character_,
  stringsAsFactors = FALSE
), file.path(dirs$summaries, "case05_v5_run_completion.csv"))
message("Case05 v5 completed: ", output)
