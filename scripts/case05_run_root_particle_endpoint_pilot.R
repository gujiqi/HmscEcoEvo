#!/usr/bin/env Rscript

# Run a complete-grid Case05 root-particle ensemble, reweight only complete
# histories at the modern endpoint, then merge the aligned cell-time outputs.
# This is a pilot for rate/root calibration before posterior-draw production.

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x)) y else x

parse_args <- function(x) {
  out <- list()
  for (arg in x) {
    if (!grepl("^--", arg)) next
    key_value <- sub("^--", "", arg)
    pos <- regexpr("=", key_value, fixed = TRUE)
    if (pos < 0L) out[[key_value]] <- TRUE else {
      out[[substr(key_value, 1L, pos - 1L)]] <-
        substr(key_value, pos + 1L, nchar(key_value))
    }
  }
  out
}
as_int <- function(x, default) {
  z <- suppressWarnings(as.integer(x %||% default))
  if (!is.finite(z)) default else z
}
as_num <- function(x, default) {
  z <- suppressWarnings(as.numeric(x %||% default))
  if (!is.finite(z)) default else z
}
safe_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}
write_csv <- function(x, path) {
  safe_dir(dirname(path))
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(x, path)
  } else utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
  invisible(path)
}
scalar_arg <- function(name, value) paste0("--", name, "=", value)
safe_id <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))

cfg <- parse_args(args)
pkg_root <- normalizePath(cfg$pkg_root %||%
  "C:/Users/Google/Documents/HMSC-HIST/HmscEcoEvo", winslash = "/", mustWork = TRUE)
rscript <- cfg$rscript %||% file.path(R.home("bin"), "Rscript.exe")
engine <- normalizePath(cfg$engine_script %||% file.path(
  pkg_root, "scripts", "case04_plant200_global_dynamic_no_bgb_final.R"
), winslash = "/", mustWork = TRUE)
merger <- normalizePath(cfg$merger_script %||% file.path(
  pkg_root, "scripts", "case04_merge_posterior_draw_shards.R"
), winslash = "/", mustWork = TRUE)
weighter <- normalizePath(cfg$weight_script %||% file.path(
  pkg_root, "scripts", "case05_weight_endpoint_shards.R"
), winslash = "/", mustWork = TRUE)
scorer <- normalizePath(cfg$endpoint_scorer %||% file.path(
  pkg_root, "scripts", "case05_score_modern_endpoint.R"
), winslash = "/", mustWork = TRUE)
response_history <- normalizePath(cfg$response_history %||% file.path(
  pkg_root, "derived_inputs", "case04_plant200_hmsc_beta_posterior_25draws_20260919",
  "case04_hmsc_beta_posterior_25draws.csv"
), winslash = "/", mustWork = TRUE)
process_parameter_history <- cfg$process_parameter_history %||% ""
if (nzchar(process_parameter_history)) {
  process_parameter_history <- normalizePath(
    process_parameter_history, winslash = "/", mustWork = TRUE
  )
}
slice_root <- file.path(pkg_root, "derived_inputs", "case04_plant200_palaeo_1deg_slices_20260919")
palaeo_slice_index <- normalizePath(cfg$palaeo_slice_index %||% file.path(
  slice_root, "case04_palaeo_earth_slice_index.csv"
), winslash = "/", mustWork = TRUE)
topographic_slice_index <- normalizePath(cfg$topographic_slice_index %||% file.path(
  slice_root, "case04_topography_slice_index.csv"
), winslash = "/", mustWork = TRUE)
output <- safe_dir(cfg$output %||% file.path(
  pkg_root, "outputs", "case05_root_particle_endpoint_pilot"
))
shard_root <- safe_dir(file.path(output, "02_root_particle_shards"))
logs <- safe_dir(file.path(output, "01_run_logs"))
root_particles <- as_int(cfg$root_particles, 12L)
n_time_slices <- as_int(cfg$n_time_slices, 20L)
max_species <- as_int(cfg$max_species, 200L)
base_seed <- as_int(cfg$seed, 20260924L)
response_draws <- strsplit(cfg$response_draw_ids %||% "chain01_iter0090",
                           ",", fixed = TRUE)[[1L]]
response_draws <- trimws(response_draws)
response_draws <- response_draws[nzchar(response_draws)]
if (root_particles < 2L || n_time_slices < 2L || max_species < 2L ||
    !length(response_draws)) {
  stop("root_particles, n_time_slices, max_species, and response_draw_ids are invalid.",
       call. = FALSE)
}

process_args <- list(
  movement_intercept = as_num(cfg$movement_intercept, -0.50),
  local_dispersal_scale_km = as_num(cfg$local_dispersal_scale_km, 200),
  establishment_intercept = as_num(cfg$establishment_intercept, -0.25),
  establishment_slope = as_num(cfg$establishment_slope, 0.85),
  persistence_intercept = as_num(cfg$persistence_intercept, 1.25),
  persistence_slope = as_num(cfg$persistence_slope, 0.75),
  root_rho = as_num(cfg$root_rho, 0.55)
)
transport_mode <- cfg$transport_mode %||% "paleomap_h3_targetcentred"
dispersal_scheme <- cfg$dispersal_scheme %||% "topographic_limited"
data_role <- cfg$data_role %||% "training_diagnostic"
likelihood_power <- as_num(cfg$endpoint_likelihood_power,
                           if (data_role == "training_diagnostic") NA_real_ else 1)

tasks <- expand.grid(response_draw = response_draws,
                     root_particle = seq_len(root_particles),
                     stringsAsFactors = FALSE)
tasks$root_seed <- base_seed + tasks$root_particle - 1L
tasks$shard <- file.path(
  shard_root,
  paste0("particle_", safe_id(tasks$response_draw), "_root",
         sprintf("%03d", tasks$root_particle))
)
tasks$log <- file.path(logs, paste0("particle_", safe_id(tasks$response_draw),
                                    "_root", sprintf("%03d", tasks$root_particle),
                                    ".log"))
tasks$status <- "pending"
write_csv(tasks, file.path(output, "00_config", "root_particle_task_table.csv"))
write_csv(data.frame(
  parameter = c("analysis_mode", "scientific_role", "grid", "transport_mode",
                "root_particle_count", "response_draw_count", "endpoint_data_role",
                "endpoint_likelihood_power", "process_parameter_history",
                names(process_args)),
  value = c(
    "empirical_core_endpoint_particle_pilot",
    "Root uncertainty is represented by complete forward-history particles. Endpoint reweighting does not modify individual past states.",
    "complete_1deg_palaeo_land_grid",
    transport_mode, root_particles, length(response_draws), data_role,
    if (is.finite(likelihood_power)) likelihood_power else "default_by_data_role",
    if (nzchar(process_parameter_history)) process_parameter_history else "not_supplied_shared_scalar_process_scenario",
    unlist(process_args)
  ),
  stringsAsFactors = FALSE
), file.path(output, "00_config", "root_particle_pilot_config.csv"))

for (i in seq_len(nrow(tasks))) {
  task <- tasks[i, , drop = FALSE]
  if (file.exists(file.path(task$shard, "20_modern_endpoint_validation",
                            paste0("modern_tip_endpoint_0Ma_", task$response_draw, ".rds")))) {
    tasks$status[[i]] <- "complete_existing"
    write_csv(tasks, file.path(output, "00_config", "root_particle_task_table.csv"))
    next
  }
  tasks$status[[i]] <- "running"
  write_csv(tasks, file.path(output, "00_config", "root_particle_task_table.csv"))
  engine_args <- c(
    engine, scalar_arg("pkg_root", pkg_root), "--quick=false",
    scalar_arg("response_history", response_history),
    scalar_arg("response_draw_ids", task$response_draw),
    scalar_arg("palaeo_slice_index", palaeo_slice_index),
    scalar_arg("topographic_slice_index", topographic_slice_index),
    scalar_arg("max_species", max_species), "--time_selection=evenly_spaced",
    scalar_arg("n_time_slices", n_time_slices),
    scalar_arg("transport_mode", transport_mode),
    scalar_arg("dispersal_scheme", dispersal_scheme), "--ldd_per_source=0",
    "--run_biotic_filtering=false", "--run_speciation_demo=false",
    "--run_lineage_extinction_demo=false", "--make_png=false", "--make_tif=false",
    "--write_attribution_state=false", "--write_modern_tip_occupancy=true",
    "--metric_storage=rds_per_time", "--checkpoint_every_draw=true",
    "--resume_checkpoint=true", "--cache_movement_matrices=true",
    "--cache_movement_in_memory=false", "--internal_dt=1",
    "--max_hazard_per_step=0.10", "--persistence_reference_myr=1",
    "--reporting_interval_myr=1", "--root_initialisation=compact_shared_root_patch",
    "--root_seed_count=1", "--root_max_cells=600", "--root_support_power=2",
    scalar_arg("movement_intercept", process_args$movement_intercept),
    scalar_arg("local_dispersal_scale_km", process_args$local_dispersal_scale_km),
    scalar_arg("establishment_intercept", process_args$establishment_intercept),
    scalar_arg("establishment_slope", process_args$establishment_slope),
    scalar_arg("persistence_intercept", process_args$persistence_intercept),
    scalar_arg("persistence_slope", process_args$persistence_slope),
    scalar_arg("root_rho", process_args$root_rho),
    scalar_arg("seed", task$root_seed), scalar_arg("particle_seed", task$root_seed),
    scalar_arg("output", task$shard)
  )
  if (nzchar(process_parameter_history)) {
    engine_args <- c(
      engine_args,
      scalar_arg("process_parameter_history", process_parameter_history)
    )
  }
  status <- system2(rscript, args = shQuote(engine_args),
                    stdout = task$log, stderr = task$log)
  tasks$status[[i]] <- if (identical(as.integer(status), 0L)) "complete" else "failed"
  write_csv(tasks, file.path(output, "00_config", "root_particle_task_table.csv"))
  if (!identical(as.integer(status), 0L)) {
    stop("Root particle failed. See: ", task$log, call. = FALSE)
  }
}

weight_args <- c(
  weighter, scalar_arg("pkg_root", pkg_root),
  scalar_arg("shard_root", shard_root), scalar_arg("output", output),
  scalar_arg("data_role", data_role)
)
if (is.finite(likelihood_power)) {
  weight_args <- c(weight_args, scalar_arg("likelihood_power", likelihood_power))
}
weight_status <- system2(rscript, args = shQuote(weight_args))
if (!identical(as.integer(weight_status), 0L)) {
  stop("Endpoint particle weighting failed.", call. = FALSE)
}
weight_file <- file.path(output, "20_modern_endpoint_validation",
                         "complete_history_endpoint_weights.csv")
merge_status <- system2(rscript, args = shQuote(c(
  merger, scalar_arg("shard_root", shard_root), scalar_arg("output", output),
  scalar_arg("shard_weights", weight_file)
)))
if (!identical(as.integer(merge_status), 0L)) {
  stop("Weighted root-particle merge failed.", call. = FALSE)
}
score_status <- system2(rscript, args = shQuote(c(
  scorer, scalar_arg("pkg_root", pkg_root), scalar_arg("scenario_output", output),
  scalar_arg("data_role", data_role)
)))
if (!identical(as.integer(score_status), 0L)) {
  stop("Merged endpoint score failed.", call. = FALSE)
}
cat("Case05 root-particle endpoint pilot completed: ", output, "\n", sep = "")
