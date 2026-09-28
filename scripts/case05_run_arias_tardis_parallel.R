#!/usr/bin/env Rscript

# Case05 standard/full executor: one genuine HMSC response draw per worker,
# complete 1 degree grid at every tree-domain time slice, three controlled P2
# schemes, then an exact cell-time posterior-mean merge.  Shards are separate
# output folders so each draw has its own checkpoint and a failed worker can be
# resumed without recomputing any finished draw.

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x)) y else x

parse_args <- function(x) {
  out <- list()
  for (arg in x) {
    if (!grepl("^--", arg)) next
    value <- sub("^--", "", arg)
    at <- regexpr("=", value, fixed = TRUE)
    if (at < 0L) out[[value]] <- TRUE else {
      out[[substr(value, 1L, at - 1L)]] <- substr(value, at + 1L, nchar(value))
    }
  }
  out
}

as_int <- function(x, default) {
  y <- suppressWarnings(as.integer(x)[1L])
  if (is.na(y) || y < 1L) default else y
}
as_num <- function(x, default) {
  y <- suppressWarnings(as.numeric(x)[1L])
  if (!is.finite(y)) default else y
}
as_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(as.character(x)[1L]) %in% c("true", "t", "1", "yes", "y")
}
safe_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}
write_csv <- function(x, path) {
  safe_dir(dirname(path))
  if (requireNamespace("data.table", quietly = TRUE)) data.table::fwrite(x, path) else {
    utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
  }
  invisible(path)
}
scalar_arg <- function(name, value) paste0("--", name, "=", value)
safe_id <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))

cfg <- parse_args(args)
pkg_root <- normalizePath(cfg$pkg_root %||%
                            "C:/Users/Google/Documents/HMSC-HIST/HmscEcoEvo",
                          winslash = "/", mustWork = TRUE)
rscript <- cfg$rscript %||% file.path(R.home("bin"), "Rscript.exe")
engine <- normalizePath(cfg$engine_script %||%
                          file.path(pkg_root,
                                    "scripts/case05_plant200_arias_tardis.R"),
                        winslash = "/", mustWork = TRUE)
merger <- normalizePath(cfg$merger_script %||%
                          file.path(pkg_root,
                                    "scripts/case04_merge_posterior_draw_shards.R"),
                        winslash = "/", mustWork = TRUE)
case05 <- normalizePath(cfg$case05_script %||%
                          file.path(pkg_root,
                                    "scripts/case05_render_arias_tardis_maps.R"),
                          winslash = "/", mustWork = TRUE)
endpoint_scorer <- normalizePath(cfg$endpoint_scorer %||%
  file.path(pkg_root, "scripts", "case05_score_modern_endpoint.R"),
  winslash = "/", mustWork = TRUE)
tutorial_builder <- normalizePath(cfg$tutorial_builder %||%
  file.path(pkg_root, "scripts", "build_case05_tutorial_docx.py"),
  winslash = "/", mustWork = FALSE)
tutorial_python <- cfg$tutorial_python %||%
  "C:/Users/Google/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe"
response_history <- normalizePath(cfg$response_history %||%
  file.path(pkg_root, "derived_inputs",
            "case04_plant200_hmsc_beta_posterior_25draws_20260919",
            "case04_hmsc_beta_posterior_25draws.csv"), winslash = "/", mustWork = TRUE)
slice_root <- file.path(pkg_root, "derived_inputs",
                        "case04_plant200_palaeo_1deg_slices_20260919")
palaeo_slice_index <- normalizePath(cfg$palaeo_slice_index %||%
  file.path(slice_root, "case04_palaeo_earth_slice_index.csv"),
  winslash = "/", mustWork = TRUE)
topographic_slice_index <- normalizePath(cfg$topographic_slice_index %||%
  file.path(slice_root, "case04_topography_slice_index.csv"),
  winslash = "/", mustWork = TRUE)
plate_transport_index <- normalizePath(cfg$plate_transport_index %||%
  file.path(pkg_root, "derived_inputs",
            "case05_paleomap_h3_grid_transport_325_0Ma_5Myr_1deg_v2_targetcentred",
            "paleomap_h3_transport_index.csv"),
  winslash = "/", mustWork = TRUE)

draw_tab <- utils::read.csv(response_history, check.names = FALSE,
                            stringsAsFactors = FALSE)
if (!"response_draw" %in% names(draw_tab)) stop("response_history lacks response_draw.")
draw_ids <- unique(as.character(draw_tab$response_draw))
requested <- cfg$response_draw_ids %||% ""
if (nzchar(requested)) {
  draw_ids <- trimws(strsplit(requested, ",", fixed = TRUE)[[1L]])
  absent <- setdiff(draw_ids, unique(as.character(draw_tab$response_draw)))
  if (length(absent)) stop("Unknown response_draw_ids: ", paste(absent, collapse = ", "))
}
if (!length(draw_ids)) stop("No posterior response draws selected.")

workers <- min(as_int(cfg$workers, 20L), length(draw_ids) * 3L,
               parallel::detectCores(logical = FALSE))
workers <- max(1L, workers)
max_species <- as_int(cfg$max_species, NA_integer_)
max_times <- as_int(cfg$max_times, NA_integer_)
time_selection <- tolower(cfg$time_selection %||%
                            if (is.finite(max_times)) "contiguous_oldest" else "all")
n_time_slices <- as_int(cfg$n_time_slices,
                        if (time_selection == "all") NA_integer_ else max_times)
if (!time_selection %in% c("all", "contiguous_oldest", "evenly_spaced")) {
  stop("time_selection must be all, contiguous_oldest or evenly_spaced.", call. = FALSE)
}
if (time_selection != "all" && (!is.finite(n_time_slices) || n_time_slices < 2L)) {
  stop("n_time_slices must be >= 2 for a non-all time selection.", call. = FALSE)
}
internal_dt <- suppressWarnings(as.numeric(cfg$internal_dt %||% 5))
movement_kernel_interval_myr <- as_num(cfg$movement_kernel_interval_myr, 1)
max_hazard_per_step <- suppressWarnings(as.numeric(cfg$max_hazard_per_step %||% 0.10))
persistence_reference_myr <- suppressWarnings(as.numeric(cfg$persistence_reference_myr %||% 5))
reporting_interval_myr <- as_num(cfg$reporting_interval_myr, 1)
map_n_time_slices <- as_int(cfg$map_n_time_slices, 20L)
transport_mode <- tolower(cfg$transport_mode %||% "paleomap_h3_targetcentred")
transport_min_target_coverage <- as_num(cfg$transport_min_target_coverage, 0.90)
root_initialisation <- tolower(cfg$root_initialisation %||% "compact_shared_root_patch")
root_seed_count <- as_int(cfg$root_seed_count, 1L)
root_support_power <- as_num(cfg$root_support_power, 2)
root_rho <- as_num(cfg$root_rho, 0.55)
root_max_cells <- as_int(cfg$root_max_cells, 600L)
diffusion_variance_rad2_per_myr <- as_num(
  cfg$diffusion_variance_rad2_per_myr, 0.00020
)
source_emigration_rate_per_myr <- as_num(
  cfg$source_emigration_rate_per_myr, 0.05
)
establishment_intercept <- as_num(cfg$establishment_intercept, -0.25)
establishment_slope <- as_num(cfg$establishment_slope, 0.85)
persistence_intercept <- as_num(cfg$persistence_intercept, 1.25)
persistence_slope <- as_num(cfg$persistence_slope, 0.75)
# The empirical-core runner requires a declared terminal consistency layer.
# It uses the smallest 0 Ma fixed-effect HMSC weight that meets the stated
# prevalence tolerance on modern 1-degree sampling cells. This remains a
# training diagnostic unless a spatially held-out or independent endpoint is
# provided; it does not edit the pre-0 Ma forward trajectory.
terminal_anchor_mode <- tolower(cfg$terminal_anchor_mode %||%
                                  "hmsc_fixed_effect_prevalence_calibrated")
terminal_anchor_weight <- as_num(cfg$terminal_anchor_weight, 0.8)
terminal_anchor_prevalence_tolerance <- cfg$terminal_anchor_prevalence_tolerance %||%
  "0.90,1.10"
tolerance_values <- suppressWarnings(as.numeric(
  strsplit(terminal_anchor_prevalence_tolerance, ",", fixed = TRUE)[[1L]]
))
if (!terminal_anchor_mode %in% c(
  "none", "hmsc_fixed_effect_convex",
  "hmsc_fixed_effect_prevalence_calibrated"
) ||
    !is.finite(terminal_anchor_weight) || terminal_anchor_weight < 0 ||
    terminal_anchor_weight > 1 || length(tolerance_values) != 2L ||
    any(!is.finite(tolerance_values)) || min(tolerance_values) <= 0 ||
    min(tolerance_values) > 1 || max(tolerance_values) < 1) {
  stop("terminal-anchor settings are invalid: use mode none, hmsc_fixed_effect_convex, or hmsc_fixed_effect_prevalence_calibrated; weight in [0,1]; and a two-value tolerance containing 1.",
       call. = FALSE)
}
seed <- as_int(cfg$seed, 20260919L)
if (!is.finite(internal_dt) || internal_dt <= 0 ||
    !is.finite(movement_kernel_interval_myr) || movement_kernel_interval_myr <= 0 ||
    !is.finite(max_hazard_per_step) || max_hazard_per_step <= 0 ||
    !is.finite(persistence_reference_myr) || persistence_reference_myr <= 0 ||
    !is.finite(reporting_interval_myr) || reporting_interval_myr <= 0 ||
    !is.finite(transport_min_target_coverage) || transport_min_target_coverage <= 0 ||
    transport_min_target_coverage > 1 ||
    !transport_mode %in% c("paleomap_h3_targetcentred") ||
    !identical(root_initialisation, "compact_shared_root_patch") ||
    root_seed_count != 1L || !is.finite(root_support_power) || root_support_power <= 0 ||
    !is.finite(diffusion_variance_rad2_per_myr) ||
    diffusion_variance_rad2_per_myr <= 0 ||
    !is.finite(source_emigration_rate_per_myr) ||
    source_emigration_rate_per_myr < 0 || !is.finite(establishment_intercept) ||
    !is.finite(establishment_slope) || !is.finite(persistence_intercept) ||
    !is.finite(persistence_slope)) {
  stop("The corrected Case05 runner requires valid numerical controls, PALEOMAP target-centred transport, and one compact root patch.", call. = FALSE)
}
maps_after_merge <- as_bool(cfg$maps_after_merge, TRUE)
map_cores <- min(as_int(cfg$map_cores, workers), parallel::detectCores(logical = FALSE))
cleanup_merged_shards <- as_bool(cfg$cleanup_merged_shards, TRUE)
# Run one dispersal scheme at a time, merge it, then remove its raw shards.
# This prevents the three 25-draw grids from accumulating simultaneously on
# an ordinary workstation.  Draws within a scheme still use the worker pool.
scenario_sequential <- as_bool(cfg$scenario_sequential, TRUE)
metric_storage <- tolower(cfg$metric_storage %||% "rds_per_time")
if (!metric_storage %in% c("rds_per_time", "csv")) {
  stop("metric_storage must be rds_per_time or csv.", call. = FALSE)
}
if (!scenario_sequential) {
  stop("Case05 now requires scenario_sequential=TRUE to keep complete-grid posterior shards within workstation disk limits.",
       call. = FALSE)
}
output <- safe_dir(cfg$output %||% file.path(
  "C:/Users/Google/Documents/HMSC-HIST/outputs/HmscEcoEvo",
  "case05_plant200_arias_tardis_full25_all66_20260922"
))
logs <- safe_dir(file.path(output, "01_run_logs"))
scenario_root <- safe_dir(file.path(output, "02_scenario_runs"))
shared_kernel_cache <- safe_dir(file.path(output, "00_config", "shared_movement_kernel_cache"))

schemes <- data.frame(
  scenario_id = c("01_arias_spherical_landscape",
                  "02_arias_tardis_topography"),
  dispersal_scheme = c("arias_spherical_landscape",
                       "arias_tardis_topography"), stringsAsFactors = FALSE
)
jobs <- do.call(rbind, lapply(seq_len(nrow(schemes)), function(i) {
  data.frame(scenario_id = schemes$scenario_id[[i]],
             dispersal_scheme = schemes$dispersal_scheme[[i]],
             response_draw = draw_ids, stringsAsFactors = FALSE)
}))
jobs$scenario_output <- file.path(scenario_root, jobs$scenario_id)
jobs$shard_output <- file.path(jobs$scenario_output, "posterior_shards",
                               paste0("draw_", safe_id(jobs$response_draw)))
jobs$log <- file.path(logs, paste0(jobs$scenario_id, "__",
                                   safe_id(jobs$response_draw), ".log"))

job_complete <- function(path, terminal_anchor_mode,
                         terminal_anchor_weight,
                         terminal_anchor_prevalence_tolerance) {
  config <- file.path(path, "00_config", "case04_final_config.csv")
  metrics <- file.path(path, "16_state_space_inference")
  validation <- file.path(path, "00_config", "case04_final_validation_checks.csv")
  endpoint <- list.files(file.path(path, "20_modern_endpoint_validation"),
                         pattern = "^modern_tip_endpoint_0Ma_.*\\.rds$",
                         full.names = TRUE)
  if (!file.exists(config) || !dir.exists(metrics) || !file.exists(validation) ||
      length(list.files(metrics, pattern = "^cell_metrics_.*\\.(csv|rds)$")) == 0L ||
      length(endpoint) != 1L) {
    return(FALSE)
  }
  prior_config <- tryCatch(utils::read.csv(config, stringsAsFactors = FALSE),
                           error = function(e) NULL)
  if (is.null(prior_config) || !all(c("parameter", "value") %in% names(prior_config))) {
    return(FALSE)
  }
  config_value <- function(parameter) {
    values <- prior_config$value[prior_config$parameter == parameter]
    if (length(values) != 1L) NA_character_ else as.character(values[[1L]])
  }
  requested_tolerance <- paste(sort(as.numeric(strsplit(
    terminal_anchor_prevalence_tolerance, ",", fixed = TRUE)[[1L]]
  )), collapse = ",")
  recorded_tolerance <- config_value("terminal_anchor_prevalence_tolerance")
  identical(config_value("terminal_anchor_mode"), terminal_anchor_mode) &&
    identical(config_value("terminal_anchor_weight"), as.character(terminal_anchor_weight)) &&
    identical(recorded_tolerance, requested_tolerance)
}
jobs$status <- ifelse(vapply(
  jobs$shard_output, job_complete, logical(1),
  terminal_anchor_mode = terminal_anchor_mode,
  terminal_anchor_weight = terminal_anchor_weight,
  terminal_anchor_prevalence_tolerance = terminal_anchor_prevalence_tolerance
),
                      "reused_complete", "pending")
write_csv(jobs, file.path(logs, "case05_posterior_draw_jobs_initial.csv"))

worker_fun <- function(job, engine, rscript, pkg_root, response_history,
                       palaeo_slice_index, topographic_slice_index,
                       max_species, max_times, time_selection, n_time_slices,
                       internal_dt, movement_kernel_interval_myr,
                       max_hazard_per_step, persistence_reference_myr,
                       reporting_interval_myr, map_n_time_slices,
                       transport_mode, plate_transport_index,
                       transport_min_target_coverage, root_initialisation,
                       root_seed_count, root_support_power, root_rho, root_max_cells,
                       diffusion_variance_rad2_per_myr,
                       source_emigration_rate_per_myr,
                       establishment_intercept, establishment_slope,
                       persistence_intercept, persistence_slope,
                       terminal_anchor_mode, terminal_anchor_weight,
                       terminal_anchor_prevalence_tolerance,
                       shared_kernel_cache, metric_storage, seed) {
  dir.create(job$shard_output, recursive = TRUE, showWarnings = FALSE)
  args <- c(
    engine, paste0("--pkg_root=", pkg_root), "--quick=false",
    paste0("--response_history=", response_history),
    paste0("--response_draw_ids=", job$response_draw),
    paste0("--palaeo_slice_index=", palaeo_slice_index),
    paste0("--topographic_slice_index=", topographic_slice_index),
    "--diagnostic_points_per_time=40", "--cache_movement_matrices=true",
    "--cache_movement_in_memory=false",
    paste0("--movement_kernel_cache_dir=", shared_kernel_cache),
    "--hazard_rate_quantile=0.999",
    paste0("--internal_dt=", internal_dt),
    paste0("--movement_kernel_interval_myr=", movement_kernel_interval_myr),
    paste0("--max_hazard_per_step=", max_hazard_per_step),
    paste0("--persistence_reference_myr=", persistence_reference_myr),
    paste0("--reporting_interval_myr=", reporting_interval_myr),
    paste0("--map_n_time_slices=", map_n_time_slices),
    paste0("--transport_mode=", transport_mode),
    paste0("--plate_transport_index=", plate_transport_index),
    paste0("--transport_min_target_coverage=", transport_min_target_coverage),
    paste0("--root_initialisation=", root_initialisation),
    paste0("--root_seed_count=", root_seed_count),
    paste0("--root_support_power=", root_support_power),
    paste0("--root_rho=", root_rho),
    paste0("--root_max_cells=", root_max_cells),
    paste0("--diffusion_variance_rad2_per_myr=", diffusion_variance_rad2_per_myr),
    paste0("--source_emigration_rate_per_myr=", source_emigration_rate_per_myr),
    paste0("--establishment_intercept=", establishment_intercept),
    paste0("--establishment_slope=", establishment_slope),
    paste0("--persistence_intercept=", persistence_intercept),
    paste0("--persistence_slope=", persistence_slope),
    paste0("--terminal_anchor_mode=", terminal_anchor_mode),
    paste0("--terminal_anchor_weight=", terminal_anchor_weight),
    paste0("--terminal_anchor_prevalence_tolerance=", terminal_anchor_prevalence_tolerance),
    paste0("--seed=", seed),
    paste0("--time_selection=", time_selection),
    paste0("--dispersal_scheme=", job$dispersal_scheme), "--ldd_per_source=0",
    "--run_biotic_filtering=false", "--run_speciation_demo=false",
    "--run_lineage_extinction_demo=false", "--make_png=false", "--make_tif=false",
    "--write_attribution_state=false", "--write_modern_tip_occupancy=true",
    paste0("--metric_storage=", metric_storage),
    "--checkpoint_every_draw=true", "--resume_checkpoint=true",
    paste0("--output=", job$shard_output)
  )
  if (is.finite(n_time_slices)) args <- c(args, paste0("--n_time_slices=", n_time_slices))
  if (is.finite(max_species)) args <- c(args, paste0("--max_species=", max_species))
  if (is.finite(max_times)) args <- c(args, paste0("--max_times=", max_times))
  status <- tryCatch(
    system2(rscript, args = shQuote(args), stdout = job$log, stderr = job$log),
    error = function(e) 1L
  )
  list(scenario_id = job$scenario_id, dispersal_scheme = job$dispersal_scheme,
       response_draw = job$response_draw, shard_output = job$shard_output,
       log = job$log, exit_status = as.integer(status),
       status = if (identical(as.integer(status), 0L)) "complete" else "failed")
}

scenario_is_merged <- function(path) {
  file.exists(file.path(path, "00_config", "case04_final_config.csv")) &&
    file.exists(file.path(path, "22_report", "case04_final_time_summary.csv")) &&
    file.exists(file.path(path, "20_modern_endpoint_validation",
                          "modern_tip_endpoint_0Ma_posterior_mean.rds")) &&
    length(list.files(file.path(path, "16_state_space_inference"),
                      pattern = "^cell_metrics_.*\\.csv$")) > 0L
}
all_results <- list()
already_merged <- vapply(schemes$scenario_id, function(id) {
  scenario_is_merged(file.path(scenario_root, id))
}, logical(1))
if (any(already_merged)) {
  jobs$status[jobs$scenario_id %in% schemes$scenario_id[already_merged]] <- "reused_aggregate"
}

run_jobs <- function(pending) {
  if (!nrow(pending)) return(data.frame())
  cluster <- parallel::makeCluster(
    min(workers, nrow(pending)),
    outfile = file.path(logs, "parallel_worker_startup.log")
  )
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  result <- parallel::parLapply(
    cluster, split(pending, seq_len(nrow(pending))), worker_fun,
    engine = engine, rscript = rscript, pkg_root = pkg_root,
    response_history = response_history, palaeo_slice_index = palaeo_slice_index,
    topographic_slice_index = topographic_slice_index,
    max_species = max_species, max_times = max_times,
    time_selection = time_selection, n_time_slices = n_time_slices,
    internal_dt = internal_dt,
    movement_kernel_interval_myr = movement_kernel_interval_myr,
    max_hazard_per_step = max_hazard_per_step,
    persistence_reference_myr = persistence_reference_myr,
    reporting_interval_myr = reporting_interval_myr,
    map_n_time_slices = map_n_time_slices,
    transport_mode = transport_mode, plate_transport_index = plate_transport_index,
    transport_min_target_coverage = transport_min_target_coverage,
    root_initialisation = root_initialisation, root_seed_count = root_seed_count,
    root_support_power = root_support_power, root_rho = root_rho,
    root_max_cells = root_max_cells,
    diffusion_variance_rad2_per_myr = diffusion_variance_rad2_per_myr,
    source_emigration_rate_per_myr = source_emigration_rate_per_myr,
    establishment_intercept = establishment_intercept,
    establishment_slope = establishment_slope,
    persistence_intercept = persistence_intercept,
    persistence_slope = persistence_slope,
    terminal_anchor_mode = terminal_anchor_mode,
    terminal_anchor_weight = terminal_anchor_weight,
    terminal_anchor_prevalence_tolerance = terminal_anchor_prevalence_tolerance,
    shared_kernel_cache = shared_kernel_cache, metric_storage = metric_storage,
    seed = seed
  )
  do.call(rbind, lapply(result, as.data.frame, stringsAsFactors = FALSE))
}

for (ss in seq_len(nrow(schemes))) {
  scenario <- schemes[ss, , drop = FALSE]
  scenario_output <- file.path(scenario_root, scenario$scenario_id)
  if (scenario_is_merged(scenario_output)) next
  pending <- jobs[jobs$scenario_id == scenario$scenario_id &
                    jobs$status == "pending", , drop = FALSE]
  if (nrow(pending)) {
    # A response draw is conditionally independent within one P2 scheme.  Run
    # those draws in parallel, then merge and remove raw RDS shards before the
    # next scheme starts.  This is intentional disk-bounded streaming, not a
    # change to the complete-grid calculation.
    results <- run_jobs(pending)
    all_results[[length(all_results) + 1L]] <- results
    for (ii in seq_len(nrow(results))) {
      hit <- jobs$scenario_id == results$scenario_id[[ii]] &
        jobs$response_draw == results$response_draw[[ii]]
      jobs$status[hit] <- results$status[[ii]]
    }
    write_csv(do.call(rbind, all_results),
              file.path(logs, "case05_posterior_draw_job_results.csv"))
    if (any(results$status != "complete")) {
      stop("At least one posterior-draw shard failed. See 01_run_logs and resume the same output directory.",
           call. = FALSE)
    }
  }
  shard_dir <- file.path(scenario_output, "posterior_shards")
  merge_log <- file.path(logs, paste0(scenario$scenario_id, "_merge.log"))
  merge_status <- system2(rscript, args = shQuote(c(
    merger, scalar_arg("shard_root", shard_dir), scalar_arg("output", scenario_output)
  )), stdout = merge_log, stderr = merge_log)
  if (!identical(as.integer(merge_status), 0L)) {
    stop("Shard merge failed for ", scenario$scenario_id, ". See ", merge_log,
         call. = FALSE)
  }
  if (cleanup_merged_shards) {
    unlink(shard_dir, recursive = TRUE, force = TRUE)
    writeLines("Raw per-draw grids were merged then removed to conserve disk space; rerun remains deterministic from the shard manifest, input hashes and response draw IDs.",
               file.path(scenario_output, "00_config", "posterior_shard_cleanup_note.txt"))
  }
}
if (length(all_results)) write_csv(do.call(rbind, all_results),
                                  file.path(logs, "case05_posterior_draw_job_results.csv"))
write_csv(jobs, file.path(logs, "case05_posterior_draw_jobs_final.csv"))

# A forward scenario is not promoted to a calibrated deep-time result until its
# 0 Ma dynamic endpoint is scored against the HMSC nowcast and declared modern
# observations. The default data role is deliberately only a training
# diagnostic; an independent or spatial holdout table can be supplied later to
# the same scoring script.
endpoint_status <- vector("list", nrow(schemes))
for (i in seq_len(nrow(schemes))) {
  scenario <- schemes[i, , drop = FALSE]
  endpoint_log <- file.path(logs, paste0(scenario$scenario_id, "_endpoint_score.log"))
  scenario_output <- file.path(scenario_root, scenario$scenario_id)
  status <- system2(rscript, args = shQuote(c(
    endpoint_scorer, scalar_arg("pkg_root", pkg_root),
    scalar_arg("scenario_output", scenario_output),
    "--data_role=training_diagnostic",
    paste0("--max_prevalence_ratio=", max(tolerance_values))
  )), stdout = endpoint_log, stderr = endpoint_log)
  endpoint_status[[i]] <- data.frame(
    scenario_id = scenario$scenario_id,
    exit_status = as.integer(status),
    log = endpoint_log,
    status = if (identical(as.integer(status), 0L)) "complete" else "failed",
    stringsAsFactors = FALSE
  )
}
endpoint_status <- do.call(rbind, endpoint_status)
write_csv(endpoint_status, file.path(logs, "case05_modern_endpoint_scoring_status.csv"))
if (any(endpoint_status$status != "complete")) {
  stop("At least one modern endpoint score failed. See 01_run_logs.", call. = FALSE)
}

write_csv(data.frame(
  response_history = response_history,
  n_real_hmsc_posterior_draws = length(draw_ids),
  workers = workers,
  time_selection = time_selection,
  n_time_slices = if (is.finite(n_time_slices)) n_time_slices else NA_integer_,
  max_species = if (is.finite(max_species)) max_species else NA_integer_,
  max_times = if (is.finite(max_times)) max_times else NA_integer_,
  internal_dt_myr = internal_dt,
  movement_kernel_interval_myr = movement_kernel_interval_myr,
  max_hazard_per_step = max_hazard_per_step,
  persistence_reference_myr = persistence_reference_myr,
  reporting_interval_myr = reporting_interval_myr,
  n_map_times = map_n_time_slices,
  transport_mode = transport_mode,
  plate_transport_index = plate_transport_index,
  transport_min_target_coverage = transport_min_target_coverage,
  root_initialisation = root_initialisation,
  root_seed_count = root_seed_count,
  root_support_power = root_support_power,
  root_rho = root_rho,
  root_max_cells = root_max_cells,
  diffusion_variance_rad2_per_myr = diffusion_variance_rad2_per_myr,
  source_emigration_rate_per_myr = source_emigration_rate_per_myr,
  establishment_intercept = establishment_intercept,
  establishment_slope = establishment_slope,
  persistence_intercept = persistence_intercept,
  persistence_slope = persistence_slope,
  terminal_anchor_mode = terminal_anchor_mode,
  terminal_anchor_weight = terminal_anchor_weight,
  terminal_anchor_prevalence_tolerance = terminal_anchor_prevalence_tolerance,
  biotic_filtering = "OFF_EMPIRICAL_CORE_NO_INDEPENDENT_INTERACTION_DATA",
  lineage_extinction = "NOT_ESTIMATED_EXTANT_TREE",
  seed = seed,
  shared_movement_kernel_cache = shared_kernel_cache,
  grid = "complete_1_degree_palaeo_land_at_every_time",
  diagnostic_points_per_time = 40L,
  posterior_aggregation = "equal_weight_draw_shards",
  scenario_execution = "sequential_scheme_parallel_draws_merge_then_cleanup",
  raw_shard_metric_storage = metric_storage,
  raw_shard_attribution_state = "disabled_to_prevent_duplicate_capped_diagnostic_tables",
  modern_endpoint_payload = "per_draw_tip_by_complete_0Ma_grid_matrices_merged_before_shard_cleanup",
  modern_endpoint_scoring = "required_training_diagnostic_before_map_and_report_phase",
  stringsAsFactors = FALSE
), file.path(output, "00_config", "case05_parallel_posterior_execution.csv"))

if (maps_after_merge) {
  status <- system2(rscript, args = shQuote(c(
    case05, scalar_arg("pkg_root", pkg_root), scalar_arg("output", output),
    "--reuse_complete_scenarios=true", "--render_maps=true",
    "--render_png=true", "--render_tif=true", scalar_arg("map_cores", map_cores),
    scalar_arg("time_selection", time_selection),
    scalar_arg("response_history", response_history),
    scalar_arg("palaeo_slice_index", palaeo_slice_index),
    scalar_arg("topographic_slice_index", topographic_slice_index),
    scalar_arg("n_time_slices", n_time_slices),
    scalar_arg("map_n_time_slices", map_n_time_slices),
    scalar_arg("transport_mode", transport_mode),
    scalar_arg("plate_transport_index", plate_transport_index),
    scalar_arg("transport_min_target_coverage", transport_min_target_coverage),
    scalar_arg("root_initialisation", root_initialisation),
    scalar_arg("root_seed_count", root_seed_count),
    scalar_arg("root_support_power", root_support_power),
    scalar_arg("root_rho", root_rho),
    scalar_arg("root_max_cells", root_max_cells),
    scalar_arg("reporting_interval_myr", reporting_interval_myr),
    "--diagnostic_points_per_time=40"
  )), stdout = file.path(logs, "case05_post_merge_maps.log"),
  stderr = file.path(logs, "case05_post_merge_maps.log"))
  if (!identical(as.integer(status), 0L)) stop("Case05 mapping/comparison phase failed.")
}

# The tutorial is a result artefact, rather than a static description.  Refresh
# it only after the final GeoTIFF/PNG map index exists, so it contains time-ordered
# atlases of the actual maps instead of placeholder graphics.
tutorial_log <- file.path(logs, "case05_post_merge_tutorial.log")
if (file.exists(tutorial_builder) && file.exists(tutorial_python)) {
  tutorial_status <- system2(
    tutorial_python,
    args = shQuote(c(tutorial_builder, "--output", output)),
    stdout = tutorial_log,
    stderr = tutorial_log
  )
  if (!identical(as.integer(tutorial_status), 0L)) {
    warning("Case05 numeric results completed, but the Word tutorial refresh failed. See ",
            tutorial_log, call. = FALSE)
  }
} else {
  writeLines(
    paste("Tutorial refresh skipped. Missing builder or Python runtime:",
          tutorial_builder, tutorial_python),
    tutorial_log
  )
}

cat("Case05 posterior draw shards completed and merged: ", output, "\n", sep = "")
