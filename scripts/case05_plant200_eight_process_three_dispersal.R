#!/usr/bin/env Rscript

# Legacy Case05: Plant200 global-grid eight-term occupancy sensitivity under
# three dispersal formulations. This is not the current six-process empirical
# Case05. This script deliberately delegates the dynamic
# calculations to Case04 and the HmscEcoEvo package. It changes only the
# dispersal kernel; all other inputs and process controls are held fixed.

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

as_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(as.character(x)[1L]) %in% c("true", "t", "1", "yes", "y")
}

as_int <- function(x, default = NA_integer_) {
  if (is.null(x) || identical(x, "")) return(default)
  y <- suppressWarnings(as.integer(x)[1L])
  if (is.na(y)) default else y
}

as_num <- function(x, default = NA_real_) {
  if (is.null(x) || identical(x, "")) return(default)
  y <- suppressWarnings(as.numeric(x)[1L])
  if (is.na(y)) default else y
}

safe_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

write_csv <- function(x, path) {
  safe_dir(dirname(path))
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(x, path)
  } else {
    utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
  }
  invisible(path)
}

read_csv <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    as.data.frame(data.table::fread(path, data.table = FALSE))
  } else utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
time_slug <- function(x) paste0(format(x, trim = TRUE, scientific = FALSE), "Ma")
scalar_arg <- function(name, value) paste0("--", name, "=", value)

cfg <- parse_args(args)
pkg_root <- normalizePath(cfg$pkg_root %||%
                            "C:/Users/Google/Documents/HMSC-HIST/HmscEcoEvo",
                          winslash = "/", mustWork = FALSE)
engine_script <- normalizePath(cfg$engine_script %||%
                                 file.path(pkg_root,
                                           "scripts/case04_plant200_global_dynamic_no_bgb_final.R"),
                               winslash = "/", mustWork = FALSE)
response_history <- cfg$response_history %||% NA_character_
palaeo_slice_index <- cfg$palaeo_slice_index %||% NA_character_
topographic_slice_index <- cfg$topographic_slice_index %||% NA_character_
rscript <- cfg$rscript %||% file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(engine_script)) stop("Missing Case04 engine: ", engine_script,
                                     call. = FALSE)
if (!file.exists(rscript)) stop("Missing Rscript executable: ", rscript,
                                call. = FALSE)

quick <- as_bool(cfg$quick, FALSE)
max_times <- as_int(cfg$max_times, if (quick) 5L else NA_integer_)
time_selection <- tolower(cfg$time_selection %||%
                            if (is.finite(max_times)) "contiguous_oldest" else "all")
n_time_slices <- as_int(cfg$n_time_slices,
                        if (time_selection == "all") NA_integer_ else max_times)
max_species <- as_int(cfg$max_species, if (quick) 20L else NA_integer_)
max_draws <- as_int(cfg$max_draws, if (quick) 1L else NA_integer_)
internal_dt <- as_num(cfg$internal_dt, if (quick) 10 else 5)
max_hazard_per_step <- as_num(cfg$max_hazard_per_step, 0.10)
persistence_reference_myr <- as_num(cfg$persistence_reference_myr, 5)
reporting_interval_myr <- as_num(cfg$reporting_interval_myr, 1)
map_n_time_slices <- as_int(cfg$map_n_time_slices, 20L)
transport_mode <- tolower(cfg$transport_mode %||% "paleomap_h3_targetcentred")
plate_transport_index <- cfg$plate_transport_index %||% file.path(
  pkg_root, "derived_inputs",
  "case05_paleomap_h3_grid_transport_325_0Ma_5Myr_1deg_v2_targetcentred",
  "paleomap_h3_transport_index.csv"
)
transport_min_target_coverage <- as_num(cfg$transport_min_target_coverage, 0.90)
root_initialisation <- tolower(cfg$root_initialisation %||% "compact_shared_root_patch")
root_seed_count <- as_int(cfg$root_seed_count, 1L)
root_support_power <- as_num(cfg$root_support_power, 2)
root_rho <- as_num(cfg$root_rho, 0.55)
root_max_cells <- as_int(cfg$root_max_cells, 600L)
diagnostic_points_per_time <- as_int(cfg$diagnostic_points_per_time, 40L)
cache_movement_matrices <- as_bool(cfg$cache_movement_matrices, TRUE)
cache_movement_in_memory <- as_bool(cfg$cache_movement_in_memory, quick)
particles_per_source <- as_int(cfg$particles_per_source, if (quick) 16L else 64L)
seed <- as_int(cfg$seed, 20260918L)
render_maps <- as_bool(cfg$render_maps, TRUE)
render_png <- as_bool(cfg$render_png, TRUE)
render_tif <- as_bool(cfg$render_tif, TRUE)
# Allow a map-only repair without recomputing the dynamic scenarios.  The two
# switches are deliberately separate because difference maps carry different
# semantic states from the ordinary continuous process maps.
render_regular_maps <- isTRUE(render_maps) && as_bool(cfg$render_regular_maps, TRUE)
render_difference_maps <- isTRUE(render_maps) && as_bool(cfg$render_difference_maps, TRUE)
map_dpi <- as_int(cfg$map_dpi, 150L)
map_cores <- as_int(cfg$map_cores, 1L)
map_cores <- max(1L, min(map_cores, parallel::detectCores(logical = FALSE)))
reuse_complete_scenarios <- as_bool(cfg$reuse_complete_scenarios, FALSE)

output <- safe_dir(cfg$output %||% file.path(
  pkg_root, "outputs", paste0("case05_plant200_eight_process_three_dispersal_", stamp())
))
dirs <- list(
  config = safe_dir(file.path(output, "00_config")),
  logs = safe_dir(file.path(output, "01_run_logs")),
  scenarios = safe_dir(file.path(output, "02_scenario_runs")),
  maps = safe_dir(file.path(output, "03_shared_scale_maps")),
  differences = safe_dir(file.path(output, "04_scheme_difference_maps")),
  summaries = safe_dir(file.path(output, "05_summaries")),
  report = safe_dir(file.path(output, "06_report"))
)

if (!requireNamespace("terra", quietly = TRUE) && render_tif) {
  stop("terra is required to create GIS GeoTIFF maps.", call. = FALSE)
}
if (!requireNamespace("ggplot2", quietly = TRUE) && render_png) {
  stop("ggplot2 is required to create PNG maps.", call. = FALSE)
}

scenarios <- data.frame(
  scenario_id = c("01_forward_geographic", "02_topographic_limited",
                  "03_particle_topographic"),
  dispersal_scheme = c("forward_geographic", "topographic_limited",
                       "particle_topographic"),
  title = c("Forward geographic dispersal",
            "Topographic-limited dispersal",
            "Particle topographic dispersal"),
  movement_cost_basis = c("great_circle_distance_km_no_topographic_penalty",
                          "topographic_effective_cost_km",
                          "topographic_effective_cost_km"),
  particle_interpretation = c(
    "not_applicable",
    "not_applicable",
    "finite_propagule_Monte_Carlo_approximation_not_posterior_particle_inference"
  ),
  stringsAsFactors = FALSE
)
scenarios$output <- vapply(scenarios$scenario_id, function(id) {
  safe_dir(file.path(dirs$scenarios, id))
}, character(1))

process_contract <- data.frame(
  process = c("Environmental filtering", "Dispersal", "Colonisation",
              "Biotic filtering", "Persistence", "Evolution", "Speciation",
              "Extinction"),
  implementation = c(
    "HMSC posterior Beta -> dated-tree ancestral response -> palaeoenvironmental support",
    "Only factor varied across the three scenarios; global land-only cell graph after PALEOMAP plate carriage",
    "Arrival hazard multiplied by conditional establishment within a CTMC update",
    "Disabled in empirical core because no independent interaction data are available",
    "Environment- and arrival-dependent local persistence; interaction term fixed to zero",
    "Direct ancestral environmental-response reconstruction from HMSC posterior draws",
    "Dated-tree copy-then-diverge lineage identity update; nodes are observed tree events",
    "Local and geography-forced loss; complete lineage extinction is not estimable from an extant-only tree"
  ),
  empirical_status = c("data_informed", "scenario_calibrated", "scenario_calibrated",
                       "off_no_independent_data", "scenario_calibrated", "data_informed",
                       "tree_conditioned", "local_and_geography_only"),
  stringsAsFactors = FALSE
)
write_csv(process_contract, file.path(dirs$config, "case05_eight_process_contract.csv"))

case05_config <- data.frame(
  parameter = c("case_id", "quick", "time_selection", "n_time_slices", "max_times", "max_species", "max_draws",
                "internal_dt_myr", "max_hazard_per_step", "persistence_reference_myr",
                "reporting_interval_myr", "map_n_time_slices",
                "transport_mode", "plate_transport_index", "transport_min_target_coverage",
                "root_initialisation", "root_seed_count", "root_support_power", "root_rho", "root_max_cells",
                "diagnostic_points_per_time",
                "palaeo_slice_index", "topographic_slice_index",
                "cache_movement_matrices", "cache_movement_in_memory",
                "particles_per_source", "seed", "render_maps",
                "render_png", "render_tif", "render_regular_maps",
                "render_difference_maps", "map_cores", "reuse_complete_scenarios",
                "response_history", "comparison_rule"),
  value = c("case05_plant200_eight_process_three_dispersal", quick, time_selection,
            n_time_slices, max_times,
            max_species, max_draws, internal_dt, max_hazard_per_step,
            persistence_reference_myr, reporting_interval_myr, map_n_time_slices,
            transport_mode, plate_transport_index, transport_min_target_coverage,
            root_initialisation, root_seed_count, root_support_power, root_rho, root_max_cells,
            diagnostic_points_per_time,
            palaeo_slice_index, topographic_slice_index,
            cache_movement_matrices, cache_movement_in_memory,
            particles_per_source, seed,
            render_maps, render_png, render_tif, render_regular_maps,
            render_difference_maps, map_cores, reuse_complete_scenarios,
            response_history,
            "same_inputs_same_tree_same_HMSC_draws_same_root_and_process_settings_change_P2_only"),
  stringsAsFactors = FALSE
)
write_csv(case05_config, file.path(dirs$config, "case05_config.csv"))
write_csv(scenarios, file.path(dirs$config, "case05_dispersal_scenarios.csv"))

run_engine <- function(scenario) {
  engine_args <- c(
    engine_script,
    scalar_arg("pkg_root", pkg_root),
    scalar_arg("quick", tolower(as.character(quick))),
    scalar_arg("time_selection", time_selection),
    scalar_arg("internal_dt", internal_dt),
    scalar_arg("max_hazard_per_step", max_hazard_per_step),
    scalar_arg("persistence_reference_myr", persistence_reference_myr),
    scalar_arg("reporting_interval_myr", reporting_interval_myr),
    scalar_arg("map_n_time_slices", map_n_time_slices),
    scalar_arg("transport_mode", transport_mode),
    scalar_arg("plate_transport_index", plate_transport_index),
    scalar_arg("transport_min_target_coverage", transport_min_target_coverage),
    scalar_arg("root_initialisation", root_initialisation),
    scalar_arg("root_seed_count", root_seed_count),
    scalar_arg("root_support_power", root_support_power),
    scalar_arg("root_rho", root_rho),
    scalar_arg("root_max_cells", root_max_cells),
    scalar_arg("diagnostic_points_per_time", diagnostic_points_per_time),
    scalar_arg("cache_movement_matrices", tolower(as.character(cache_movement_matrices))),
    scalar_arg("cache_movement_in_memory", tolower(as.character(cache_movement_in_memory))),
    scalar_arg("seed", seed),
    scalar_arg("particle_seed", seed),
    scalar_arg("particles_per_source", particles_per_source),
    scalar_arg("dispersal_scheme", scenario$dispersal_scheme),
    scalar_arg("ldd_per_source", 0L),
    scalar_arg("make_png", "false"),
    scalar_arg("make_tif", "false"),
    scalar_arg("run_biotic_filtering", "false"),
    scalar_arg("run_speciation_demo", "false"),
    scalar_arg("run_lineage_extinction_demo", "false"),
    scalar_arg("write_attribution_state", "true"),
    scalar_arg("output", scenario$output)
  )
  if (is.finite(n_time_slices)) engine_args <- c(engine_args, scalar_arg("n_time_slices", n_time_slices))
  if (is.finite(max_times)) engine_args <- c(engine_args, scalar_arg("max_times", max_times))
  if (is.finite(max_species)) engine_args <- c(engine_args, scalar_arg("max_species", max_species))
  if (is.finite(max_draws)) engine_args <- c(engine_args, scalar_arg("max_draws", max_draws))
  if (!is.na(response_history) && nzchar(response_history)) {
    engine_args <- c(engine_args, scalar_arg("response_history", response_history))
  }
  if (!is.na(palaeo_slice_index) && nzchar(palaeo_slice_index)) {
    engine_args <- c(engine_args, scalar_arg("palaeo_slice_index", palaeo_slice_index))
  }
  if (!is.na(topographic_slice_index) && nzchar(topographic_slice_index)) {
    engine_args <- c(engine_args, scalar_arg("topographic_slice_index", topographic_slice_index))
  }
  log_path <- file.path(dirs$logs, paste0(scenario$scenario_id, ".log"))
  status_path <- file.path(dirs$logs, paste0(scenario$scenario_id, "_status.csv"))
  started <- Sys.time()
  output_text <- tryCatch(
    system2(rscript, args = shQuote(engine_args), stdout = TRUE, stderr = TRUE),
    error = function(e) structure(conditionMessage(e), status = 1L)
  )
  exit_status <- attr(output_text, "status") %||% 0L
  writeLines(as.character(output_text), log_path, useBytes = TRUE)
  state <- data.frame(
    scenario_id = scenario$scenario_id,
    dispersal_scheme = scenario$dispersal_scheme,
    status = if (identical(as.integer(exit_status), 0L)) "complete" else "failed",
    exit_status = as.integer(exit_status),
    started_at = format(started, tz = "UTC"),
    finished_at = format(Sys.time(), tz = "UTC"),
    log = normalizePath(log_path, winslash = "/", mustWork = FALSE),
    output = scenario$output,
    stringsAsFactors = FALSE
  )
  write_csv(state, status_path)
  if (!identical(as.integer(exit_status), 0L)) {
    stop("Case05 scenario failed: ", scenario$scenario_id,
         ". See ", log_path, call. = FALSE)
  }
  state
}

scenario_is_reusable <- function(scenario) {
  config_path <- file.path(scenario$output, "00_config", "case04_final_config.csv")
  summary_path <- file.path(scenario$output, "22_report", "case04_final_time_summary.csv")
  metric_path <- file.path(scenario$output, "16_state_space_inference")
  endpoint_path <- file.path(scenario$output, "20_modern_endpoint_validation",
                             "modern_tip_endpoint_0Ma_posterior_mean.rds")
  if (!file.exists(config_path) || !file.exists(summary_path) || !dir.exists(metric_path) ||
      !file.exists(endpoint_path)) {
    return(FALSE)
  }
  config <- tryCatch(read_csv(config_path), error = function(e) NULL)
  if (is.null(config) || !all(c("parameter", "value") %in% names(config))) return(FALSE)
  scheme <- config$value[match("dispersal_scheme", config$parameter)]
  identical(as.character(scheme), as.character(scenario$dispersal_scheme)) &&
    length(list.files(metric_path, pattern = "^cell_metrics_.*\\.csv$")) > 0L
}

status_rows <- vector("list", nrow(scenarios))
for (i in seq_len(nrow(scenarios))) {
  scenario <- scenarios[i, , drop = FALSE]
  if (reuse_complete_scenarios && scenario_is_reusable(scenario)) {
    message("[Case05] Reusing verified output for ", scenario$scenario_id, ".")
    status_rows[[i]] <- data.frame(
      scenario_id = scenario$scenario_id,
      dispersal_scheme = scenario$dispersal_scheme,
      status = "reused_complete", exit_status = 0L,
      started_at = format(Sys.time(), tz = "UTC"),
      finished_at = format(Sys.time(), tz = "UTC"),
      log = NA_character_, output = scenario$output, stringsAsFactors = FALSE
    )
  } else {
    message("[Case05] Running ", scenario$scenario_id, ".")
    status_rows[[i]] <- run_engine(scenario)
  }
}
status_table <- do.call(rbind, status_rows)
write_csv(status_table, file.path(dirs$logs, "case05_scenario_run_status.csv"))

metric_dirs <- c(
  expected_lineage_richness = "expected_lineage_richness",
  relative_environmental_support_sum = "relative_environmental_support_sum",
  binary_lineage_richness = "binary_lineage_richness",
  mean_occupancy_probability = "mean_occupancy_probability",
  mean_environmental_support = "mean_environmental_support",
  mean_arrival_hazard_per_myr = "mean_arrival_hazard_per_myr",
  mean_arrival_probability_reporting_interval = "mean_arrival_probability_reporting_interval",
  mean_colonisation_hazard_per_myr = "mean_colonisation_hazard_per_myr",
  mean_colonisation_probability_reporting_interval = "mean_colonisation_probability_reporting_interval",
  mean_persistence_probability = "mean_persistence_probability",
  mean_local_extinction_probability = "mean_local_extinction_probability",
  weighted_endemism = "weighted_endemism",
  shannon_diversity = "occupancy_weighted_shannon_entropy",
  simpson_diversity = "occupancy_weighted_gini_simpson",
  occupancy_weighted_effective_lineages_shannon = "occupancy_weighted_effective_lineages_shannon",
  occupancy_weighted_effective_lineages_simpson = "occupancy_weighted_effective_lineages_simpson",
  occupancy_weighted_diversity_defined = "occupancy_weighted_diversity_defined",
  niche_response_dispersion = "niche_response_dispersion"
)
metric_file_id <- c(
  expected_lineage_richness = "lineage_richness",
  relative_environmental_support_sum = "relative_environmental_support_sum",
  binary_lineage_richness = "binary_richness",
  mean_occupancy_probability = "occupancy",
  mean_environmental_support = "environmental_support",
  mean_arrival_hazard_per_myr = "arrival_hazard_per_myr",
  mean_arrival_probability_reporting_interval = "arrival_probability_reporting_interval",
  mean_colonisation_hazard_per_myr = "colonisation_hazard_per_myr",
  mean_colonisation_probability_reporting_interval = "colonisation_probability_reporting_interval",
  mean_persistence_probability = "persistence",
  mean_local_extinction_probability = "local_extinction",
  weighted_endemism = "weighted_endemism",
  shannon_diversity = "occupancy_weighted_shannon_entropy",
  simpson_diversity = "occupancy_weighted_gini_simpson",
  occupancy_weighted_effective_lineages_shannon = "effective_lineages_shannon",
  occupancy_weighted_effective_lineages_simpson = "effective_lineages_simpson",
  occupancy_weighted_diversity_defined = "occupancy_weighted_diversity_defined",
  niche_response_dispersion = "niche_response_dispersion"
)
# Directory names and index labels remain fully descriptive. File stems use
# compact stable codes so a legitimate 1-degree GIS map cannot silently lose
# its extension at Windows' 260-character path boundary.
metric_map_code <- c(
  expected_lineage_richness = "linrich",
  relative_environmental_support_sum = "envsupportsum",
  binary_lineage_richness = "binrich",
  mean_occupancy_probability = "occupancy",
  mean_environmental_support = "envsupport",
  mean_arrival_hazard_per_myr = "arrivalhaz",
  mean_arrival_probability_reporting_interval = "arrivalprob",
  mean_colonisation_hazard_per_myr = "colonisehaz",
  mean_colonisation_probability_reporting_interval = "coloniseprob",
  mean_persistence_probability = "persist",
  mean_local_extinction_probability = "localext",
  weighted_endemism = "endemism",
  shannon_diversity = "owshannon",
  simpson_diversity = "owgini",
  occupancy_weighted_effective_lineages_shannon = "oweffshannon",
  occupancy_weighted_effective_lineages_simpson = "oweffsimpson",
  occupancy_weighted_diversity_defined = "owdefined",
  niche_response_dispersion = "nrd"
)
comparison_file_code <- c(
  topographic_minus_forward = "topo_minus_forward",
  particle_minus_topographic = "particle_minus_topo"
)
# Long, self-describing metric directory names are useful in a spreadsheet but
# exceed the Win32 path limit once combined with an output root, comparison and
# image filename.  Keep those names in the map index; use compact, stable
# folders on disk so that .png is never silently truncated to .p.
metric_map_folder <- stats::setNames(
  paste0("m_", unname(metric_map_code)), names(metric_map_code)
)
comparison_map_folder <- c(
  topographic_minus_forward = "d_topo_forward",
  particle_minus_topographic = "d_particle_topo"
)
endpoint_map_code <- c(
  hmsc_fixed_effect_nowcast_expected_richness = "hmscnowcast",
  forward_raw_expected_lineage_richness = "forwardraw",
  terminal_anchor_difference_expected_lineage_richness = "anchordelta",
  dynamic_endpoint_expected_lineage_richness = "dynamicendpoint"
)
map_codebook <- data.frame(
  metric = names(metric_map_folder),
  metric_label = unname(metric_dirs[names(metric_map_folder)]),
  map_folder = unname(metric_map_folder),
  map_file_code = unname(metric_map_code[names(metric_map_folder)]),
  stringsAsFactors = FALSE
)
map_codebook <- rbind(
  map_codebook,
  data.frame(
    metric = names(endpoint_map_code),
    metric_label = c(
      "0 Ma HMSC fixed-effect nowcast expected richness",
      "0 Ma forward-only expected sampled-surviving-lineage richness",
      "0 Ma terminal-anchor minus forward-only expected richness",
      "0 Ma endpoint-ensemble expected sampled-surviving-lineage richness"
    ),
    map_folder = paste0("m_", unname(endpoint_map_code)),
    map_file_code = unname(endpoint_map_code),
    stringsAsFactors = FALSE
  )
)
write_csv(map_codebook, file.path(dirs$config, "case05_map_folder_codebook.csv"))
write_csv(data.frame(
  comparison_id = names(comparison_map_folder),
  map_folder = unname(comparison_map_folder),
  file_code = unname(comparison_file_code[names(comparison_map_folder)]),
  stringsAsFactors = FALSE
), file.path(dirs$config, "case05_difference_folder_codebook.csv"))
probability_metrics <- c(
  "mean_occupancy_probability", "mean_environmental_support",
  "mean_arrival_probability_reporting_interval",
  "mean_colonisation_probability_reporting_interval",
  "mean_persistence_probability", "mean_local_extinction_probability",
  "simpson_diversity", "occupancy_weighted_diversity_defined"
)
richness_metrics <- c("expected_lineage_richness",
                       "relative_environmental_support_sum",
                       "binary_lineage_richness",
                       "occupancy_weighted_effective_lineages_shannon",
                       "occupancy_weighted_effective_lineages_simpson")

scenario_metric_files <- function(scenario_output) {
  files <- list.files(file.path(scenario_output, "16_state_space_inference"),
                      pattern = "^cell_metrics_.*\\.csv$", full.names = TRUE)
  if (!length(files)) stop("No per-time metric tables in ", scenario_output,
                           call. = FALSE)
  # A reused output directory can contain diagnostic tables from an earlier
  # smoke run with a different selected time window.  The latest engine summary
  # is the run manifest, so retain only its exact time axis rather than mixing
  # stale CSVs into scale calculations or validation.
  summary_path <- file.path(scenario_output, "22_report", "case04_final_time_summary.csv")
  if (!file.exists(summary_path)) {
    stop("Missing engine time manifest: ", summary_path, call. = FALSE)
  }
  target_times <- suppressWarnings(as.numeric(read_csv(summary_path)$time_ma))
  target_times <- sort(unique(target_times[is.finite(target_times)]), decreasing = TRUE)
  if (!length(target_times)) stop("Engine time manifest has no finite time values: ",
                                  summary_path, call. = FALSE)
  file_times <- suppressWarnings(as.numeric(sub(
    "^cell_metrics_(-?[0-9]+(?:\\.[0-9]+)?)Ma\\.csv$", "\\1", basename(files)
  )))
  keep <- is.finite(file_times) & file_times %in% target_times
  files <- files[keep]
  file_times <- file_times[keep]
  if (!setequal(file_times, target_times) || anyDuplicated(file_times)) {
    stop("Per-time metric files do not match the latest engine time manifest in ",
         scenario_output, call. = FALSE)
  }
  files[order(file_times, decreasing = TRUE)]
}
metric_files <- setNames(lapply(scenarios$output, scenario_metric_files),
                         scenarios$scenario_id)

read_time_metric <- function(path) {
  x <- read_csv(path)
  required <- c("cell_id", "lon", "lat", "time_ma", names(metric_dirs))
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Metric table is incomplete: ", path,
                            ". Missing ", paste(missing, collapse = ", "),
                            call. = FALSE)
  x
}

first_table <- read_time_metric(metric_files[[1L]][[1L]])
# Keep the mapper backwards-compatible with already completed legacy runs, but
# do not silently draw absent process layers. New endpoint-calibrated runs add
# the occupancy-weighted diversity fields below.
metric_dirs <- metric_dirs[names(metric_dirs) %in% names(first_table)]
metric_file_id <- metric_file_id[names(metric_dirs)]
metric_map_code <- metric_map_code[names(metric_dirs)]
probability_metrics <- probability_metrics[probability_metrics %in% names(metric_dirs)]
richness_metrics <- richness_metrics[richness_metrics %in% names(metric_dirs)]
n_species <- max(first_table$expected_lineage_richness, na.rm = TRUE)
engine_config <- read_csv(file.path(scenarios$output[[1L]], "00_config",
                                    "case04_final_config.csv"))
n_species_input <- suppressWarnings(as.numeric(
  engine_config$value[match("n_species", engine_config$parameter)]
))
if (is.finite(n_species_input)) n_species <- max(n_species, n_species_input)
if (!is.finite(n_species) || n_species <= 0) n_species <- 200

global_max <- setNames(rep(0, length(metric_dirs)), names(metric_dirs))
time_values <- numeric()
scenario_time_files <- vector("list", length(metric_files))
names(scenario_time_files) <- names(metric_files)
for (scenario_id in names(metric_files)) {
  table_times <- numeric(length(metric_files[[scenario_id]]))
  for (path in metric_files[[scenario_id]]) {
    x <- read_time_metric(path)
    unique_time <- unique(x$time_ma)
    if (length(unique_time) != 1L || !is.finite(unique_time)) {
      stop("Metric file must contain one finite time slice: ", path, call. = FALSE)
    }
    table_times[[which(metric_files[[scenario_id]] == path)[1L]]] <- unique_time
    time_values <- c(time_values, unique(x$time_ma))
    for (metric in names(metric_dirs)) {
      value <- suppressWarnings(as.numeric(x[[metric]]))
      candidate <- suppressWarnings(max(value[is.finite(value)], na.rm = TRUE))
      if (is.finite(candidate)) global_max[[metric]] <- max(global_max[[metric]], candidate)
    }
  }
  if (anyDuplicated(table_times)) {
    stop("Scenario has duplicate per-time metric tables: ", scenario_id, call. = FALSE)
  }
  scenario_time_files[[scenario_id]] <- stats::setNames(metric_files[[scenario_id]],
                                                         as.character(table_times))
}
times <- sort(unique(time_values), decreasing = TRUE)
if (!length(times)) stop("No projected time values were recovered.", call. = FALSE)
map_times <- times
if (is.finite(map_n_time_slices) && map_n_time_slices >= 2L &&
    length(times) > map_n_time_slices) {
  map_times <- times[unique(round(seq.int(
    1L, length(times), length.out = map_n_time_slices
  )))]
}
write_csv(data.frame(
  time_ma = times,
  selected_for_map_output = times %in% map_times,
  n_total_dynamic_time_slices = length(times),
  n_rendered_time_slices = length(map_times),
  stringsAsFactors = FALSE
), file.path(dirs$config, "case05_time_axis_and_map_selection.csv"))

# Audit the one directional causal constraint before drawing maps: zero P2
# arrival must force zero colonisation. The converse is not required because
# a finite arrival can still have effectively zero establishment in a strongly
# unsuitable cell (or underflow after summary output).
colonisation_zero_audit <- list()
for (scenario_id in names(metric_files)) {
  for (path in metric_files[[scenario_id]]) {
    x <- read_time_metric(path)
    arrival <- suppressWarnings(as.numeric(
      x$mean_arrival_probability_reporting_interval
    ))
    colonisation <- suppressWarnings(as.numeric(
      x$mean_colonisation_probability_reporting_interval
    ))
    finite <- is.finite(arrival) & is.finite(colonisation)
    zero_arrival <- finite & abs(arrival) <= 1e-15
    zero_colonisation <- finite & abs(colonisation) <= 1e-15
    colonisation_zero_audit[[length(colonisation_zero_audit) + 1L]] <- data.frame(
      scenario_id = scenario_id,
      time_ma = unique(x$time_ma)[1L],
      n_available_cells = sum(finite),
      zero_arrival_cells = sum(zero_arrival),
      zero_colonisation_cells = sum(zero_colonisation),
      zero_arrival_zero_colonisation_cells = sum(zero_arrival & zero_colonisation),
      nonzero_colonisation_when_arrival_zero = sum(zero_arrival & !zero_colonisation),
      zero_colonisation_when_arrival_nonzero = sum(!zero_arrival & zero_colonisation),
      interpretation = "arrival_zero_must_imply_colonisation_zero; positive_arrival_with_zero_colonisation_can_be_establishment_limited",
      stringsAsFactors = FALSE
    )
  }
}
colonisation_zero_audit <- do.call(rbind, colonisation_zero_audit)
colonisation_zero_audit$zero_colonisation_pct <- ifelse(
  colonisation_zero_audit$n_available_cells > 0,
  100 * colonisation_zero_audit$zero_colonisation_cells /
    colonisation_zero_audit$n_available_cells,
  NA_real_
)
write_csv(colonisation_zero_audit, file.path(dirs$summaries,
                                             "case05_colonisation_zero_arrival_audit.csv"))
zero_report_times <- sort(unique(c(max(times), 40, 0)), decreasing = TRUE)
zero_report_times <- zero_report_times[zero_report_times %in% unique(colonisation_zero_audit$time_ma)]
zero_report_rows <- colonisation_zero_audit[
  colonisation_zero_audit$time_ma %in% zero_report_times,
  c("scenario_id", "time_ma", "n_available_cells", "zero_arrival_cells",
    "zero_colonisation_cells", "zero_colonisation_pct"), drop = FALSE
]
zero_report_lines <- c(
  "# Case05 colonisation-zero interpretation",
  "",
  "A zero colonisation value is not a missing map cell. In this model, zero P2 arrival forces zero colonisation. The converse is not guaranteed: a cell with non-zero arrival can still have a zero or numerically negligible colonisation value when conditional establishment is extremely low. Therefore, source limitation must be read from the arrival map and this audit, not inferred from a colonisation map alone.",
  "",
  "The corrected controlled scenarios use complete 1-degree palaeoland grids, external PALEOMAP target-centred plate carriage before active movement, `ldd_per_source = 0`, and a single contiguous environmental root patch. Valid land cells with zero process intensity remain explicit zeroes rather than omitted map cells. Because the run is forward and not endpoint/fossil conditioned, a zero arrival value is a scenario outcome, not observed palaeofloral absence.",
  "",
  "On revised difference maps: medium grey = both scenarios have a valid zero value; white = equal non-zero; red/blue = numerical difference; dark grey = ocean or unavailable habitat. For colonisation, consult the paired arrival map/audit to distinguish source-limited from establishment-limited zeros.",
  "",
  "## Representative audit rows",
  "",
  "```text",
  capture.output(print(zero_report_rows, row.names = FALSE)),
  "```",
  "",
  "See `case05_colonisation_zero_arrival_audit.csv` for every time slice and scenario, and `case05_difference_state_audit.csv` for every difference map."
)
writeLines(zero_report_lines, file.path(dirs$report,
                                        "case05_colonisation_zero_interpretation.md"),
           useBytes = TRUE)

metric_limits <- data.frame(
  metric = names(metric_dirs),
  lower = 0,
  upper = unname(global_max[names(metric_dirs)]),
  scale_scope = "all_three_dispersal_scenarios_all_projected_times",
  stringsAsFactors = FALSE
)
metric_limits$upper[metric_limits$metric %in% probability_metrics] <- 1
metric_limits$upper[metric_limits$metric %in% richness_metrics] <- n_species
metric_limits$upper[metric_limits$metric == "shannon_diversity"] <- log(max(2, n_species))
metric_limits$upper[!is.finite(metric_limits$upper) | metric_limits$upper <= 0] <- 1
write_csv(metric_limits, file.path(dirs$config, "case05_shared_metric_limits.csv"))

global_raster <- function(df, value_col) {
  r <- terra::rast(ncols = 360, nrows = 180, xmin = -180, xmax = 180,
                   ymin = -90, ymax = 90, crs = "EPSG:4326")
  ids <- terra::cellFromXY(r, as.matrix(df[, c("lon", "lat"), drop = FALSE]))
  ok <- !is.na(ids) & is.finite(df[[value_col]])
  r[ids[ok]] <- as.numeric(df[[value_col]][ok])
  r
}

clip_to_limits <- function(x, limits) pmin(pmax(x, limits[[1L]]), limits[[2L]])
oob_squish <- function(x, range, ...) pmin(pmax(x, range[[1L]]), range[[2L]])

write_map <- function(df, metric, label, limits, out_dir, scenario_id,
                      title_prefix = "Case05 Plant200 predicted") {
  safe_dir(out_dir)
  tm <- unique(df$time_ma)
  if (length(tm) != 1L) stop("A map table must contain exactly one time slice.",
                             call. = FALSE)
  stem <- paste0(metric_map_code[[metric]], "_", time_slug(tm), "_", scenario_id)
  tif <- file.path(out_dir, paste0(stem, ".tif"))
  png <- file.path(out_dir, paste0(stem, ".png"))
  values <- clip_to_limits(suppressWarnings(as.numeric(df[[metric]])), limits)
  map_df <- df[, c("lon", "lat"), drop = FALSE]
  map_df$value <- values
  if (render_tif) {
    raster <- global_raster(transform(map_df, value = values), "value")
    terra::writeRaster(raster, tif, overwrite = TRUE,
                       gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
  }
  if (render_png) {
    plot <- ggplot2::ggplot(map_df, ggplot2::aes(lon, lat)) +
      ggplot2::geom_tile(ggplot2::aes(fill = value), width = 1, height = 1) +
      ggplot2::coord_equal(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
      ggplot2::scale_fill_gradientn(
        colours = c("#fff7bc", "#9ecae1", "#3182bd", "#08519c", "#08306b"),
        limits = limits, oob = oob_squish, na.value = "#6c757d", name = label
      ) +
      ggplot2::labs(
        title = paste(title_prefix, label),
        subtitle = paste0("", tm, " Ma; fixed legend across all three dispersal scenarios and all time slices"),
        x = "Longitude (degrees)", y = "Latitude (degrees)",
        caption = "All valid palaeo-land cells are mapped; pale yellow is a valid zero and grey is ocean or unavailable habitat. Values are model outputs, not observed palaeodistributions."
      ) +
      ggplot2::theme_minimal(base_size = 9) +
      ggplot2::theme(panel.grid = ggplot2::element_line(colour = "grey88", linewidth = 0.2),
                     panel.background = ggplot2::element_rect(fill = "#6c757d", colour = NA))
    ggplot2::ggsave(png, plot, width = 10, height = 5.5, dpi = map_dpi)
    if (!file.exists(png)) stop("PNG map was not written: ", png, call. = FALSE)
  }
  data.frame(
    scenario_id = scenario_id, metric = metric, label = label, time_ma = tm,
    lower = limits[[1L]], upper = limits[[2L]],
    tif = normalizePath(tif, winslash = "/", mustWork = FALSE),
    png = normalizePath(png, winslash = "/", mustWork = FALSE),
    crs = "EPSG:4326", resolution_deg = 1,
    source = "Case05 HmscEE eight-process scenario output",
    interpretation = "fixed_colour_scale_across_all_schemes_and_times",
    stringsAsFactors = FALSE
  )
}

run_regular_map_job <- function(job) {
  x <- read_time_metric(job$path)
  write_map(x, job$metric, job$label, job$limits, job$out_dir, job$scenario_id)
}

run_map_jobs <- function(jobs, worker) {
  if (!length(jobs)) return(list())
  if (map_cores <= 1L) return(lapply(jobs, worker))
  cluster <- parallel::makeCluster(map_cores)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  export_names <- c("safe_dir", "read_csv", "read_time_metric", "global_raster",
                    "clip_to_limits", "oob_squish", "time_slug", "metric_file_id",
                    "metric_map_code", "comparison_file_code",
                    "write_map", "write_difference_map", "run_regular_map_job",
                    "run_difference_map_job", "render_tif", "render_png", "map_dpi",
                    "metric_dirs")
  export_names <- export_names[vapply(export_names, exists, logical(1),
                                      envir = .GlobalEnv, inherits = FALSE)]
  parallel::clusterExport(cluster, varlist = export_names, envir = .GlobalEnv)
  parallel::parLapply(cluster, jobs, worker)
}

map_index_path <- file.path(dirs$summaries, "case05_map_index.csv")
existing_map_index <- if (file.exists(map_index_path)) read_csv(map_index_path) else data.frame()
drop_existing_ids <- c(
  if (render_regular_maps) scenarios$scenario_id else character(),
  if (render_difference_maps) c("topographic_minus_forward",
                                "particle_minus_topographic") else character()
)
preserved_map_index <- if (nrow(existing_map_index) && length(drop_existing_ids)) {
  existing_map_index[!existing_map_index$scenario_id %in% drop_existing_ids, , drop = FALSE]
} else existing_map_index
map_index <- if (nrow(preserved_map_index)) list(preserved_map_index) else list()

if (render_regular_maps) {
  regular_jobs <- list()
  for (i in seq_len(nrow(scenarios))) {
    scenario_id <- scenarios$scenario_id[[i]]
    for (tm in map_times) {
      path <- unname(scenario_time_files[[scenario_id]][[as.character(tm)]])
      for (metric in names(metric_dirs)) {
        limits <- unlist(metric_limits[metric_limits$metric == metric,
                                       c("lower", "upper"), drop = FALSE])
        regular_jobs[[length(regular_jobs) + 1L]] <- list(
          path = path, metric = metric, label = metric_dirs[[metric]], limits = limits,
          out_dir = file.path(dirs$maps, scenario_id, metric_map_folder[[metric]]),
          scenario_id = scenario_id
        )
      }
    }
  }
  map_index <- c(map_index, run_map_jobs(regular_jobs, run_regular_map_job))
}

# At 0 Ma, the appropriate comparison for a dynamic endpoint is a fixed-effect
# HMSC nowcast that restores the modern tip intercept. The intercept-free
# ancestral response score is intentionally not used for this endpoint check.
endpoint_map_rows <- list()
endpoint_summary_rows <- list()
for (i in seq_len(nrow(scenarios))) {
  scenario_id <- scenarios$scenario_id[[i]]
  endpoint_path <- file.path(scenarios$output[[i]],
                             "20_modern_endpoint_validation",
                             "modern_tip_endpoint_0Ma_posterior_mean.rds")
  if (!file.exists(endpoint_path)) next
  endpoint <- readRDS(endpoint_path)
  if (!all(c("cell_id", "lon", "lat", "species",
             "dynamic_occupancy_probability",
             "hmsc_fixed_effect_nowcast_probability") %in% names(endpoint))) {
    stop("Endpoint payload has an incomplete schema: ", endpoint_path, call. = FALSE)
  }
  endpoint_df <- data.frame(
    cell_id = endpoint$cell_id,
    lon = endpoint$lon,
    lat = endpoint$lat,
    time_ma = 0,
    hmsc_fixed_effect_nowcast_expected_richness = rowSums(
      endpoint$hmsc_fixed_effect_nowcast_probability, na.rm = TRUE
    ),
    dynamic_endpoint_expected_lineage_richness = rowSums(
      endpoint$dynamic_occupancy_probability, na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
  if ("forward_raw_occupancy_probability" %in% names(endpoint) &&
      !is.null(endpoint$forward_raw_occupancy_probability)) {
    endpoint_df$forward_raw_expected_lineage_richness <- rowSums(
      endpoint$forward_raw_occupancy_probability, na.rm = TRUE
    )
    endpoint_df$terminal_anchor_difference_expected_lineage_richness <-
      endpoint_df$dynamic_endpoint_expected_lineage_richness -
      endpoint_df$forward_raw_expected_lineage_richness
  }
  endpoint_summary_rows[[scenario_id]] <- transform(
    endpoint_df[, c("cell_id", "lon", "lat",
                    intersect(c(
                      "hmsc_fixed_effect_nowcast_expected_richness",
                      "forward_raw_expected_lineage_richness",
                      "terminal_anchor_difference_expected_lineage_richness",
                      "dynamic_endpoint_expected_lineage_richness"
                    ), names(endpoint_df)))],
    scenario_id = scenario_id,
    dynamic_to_nowcast_ratio = ifelse(
      hmsc_fixed_effect_nowcast_expected_richness > 0,
      dynamic_endpoint_expected_lineage_richness /
        hmsc_fixed_effect_nowcast_expected_richness,
      NA_real_
    )
  )
  for (metric in intersect(names(endpoint_map_code), names(endpoint_df))) {
    endpoint_map_rows[[paste(scenario_id, metric, sep = "\r")]] <- list(
      scenario_id = scenario_id, metric = metric, values = endpoint_df,
      terminal_anchor_mode = endpoint$terminal_anchor_mode %||% "none",
      terminal_anchor_weight = endpoint$terminal_anchor_weight %||% 0
    )
  }
}
if (length(endpoint_summary_rows)) {
  endpoint_summary <- do.call(rbind, endpoint_summary_rows)
  write_csv(endpoint_summary, file.path(dirs$summaries,
                                        "case05_0Ma_dynamic_vs_hmsc_nowcast_cell_summary.csv"))
}
if (render_regular_maps && length(endpoint_map_rows)) {
  endpoint_upper <- n_species
  for (job in endpoint_map_rows) {
    scenario_id <- job$scenario_id
    metric <- job$metric
    endpoint_df <- job$values
    # Temporarily add the endpoint map code to the standard writer's compact
    # filename registry; all schemes use the same [0, n_species] scale.
    metric_map_code[[metric]] <- endpoint_map_code[[metric]]
    out_dir <- file.path(dirs$maps, scenario_id,
                         paste0("m_", endpoint_map_code[[metric]]))
    label <- switch(
      metric,
      hmsc_fixed_effect_nowcast_expected_richness =
        "0 Ma HMSC fixed-effect nowcast expected richness",
      forward_raw_expected_lineage_richness =
        "0 Ma forward-only expected sampled-surviving-lineage richness",
      terminal_anchor_difference_expected_lineage_richness =
        "0 Ma terminal-anchor minus forward-only expected richness",
      dynamic_endpoint_expected_lineage_richness = paste0(
        "0 Ma endpoint-ensemble expected sampled-surviving-lineage richness; ",
        job$terminal_anchor_mode, " weight=", job$terminal_anchor_weight
      )
    )
    endpoint_limits <- if (identical(metric,
                                     "terminal_anchor_difference_expected_lineage_richness")) {
      c(0, max(endpoint_df[[metric]], na.rm = TRUE))
    } else c(0, endpoint_upper)
    endpoint_index <- write_map(endpoint_df, metric, label,
                                endpoint_limits, out_dir, scenario_id,
                                title_prefix = "Case05 Plant200 endpoint diagnostic")
    endpoint_index$source <- "0 Ma merged HMSC posterior endpoint payload"
    endpoint_index$interpretation <-
      paste(
        "endpoint_diagnostic; raw-forward, endpoint ensemble and HMSC nowcast",
        "are retained separately; not independent validation unless a holdout",
        "or independent data_role is supplied"
      )
    map_index <- c(map_index, list(endpoint_index))
  }
}

comparison_metrics <- c("expected_lineage_richness", "mean_occupancy_probability",
                          "mean_arrival_probability_reporting_interval",
                          "mean_colonisation_probability_reporting_interval",
                         "mean_local_extinction_probability", "weighted_endemism",
                         "shannon_diversity", "simpson_diversity",
                          "niche_response_dispersion")
comparison_metrics <- comparison_metrics[comparison_metrics %in% names(metric_dirs)]
comparison_specs <- data.frame(
  comparison_id = c("topographic_minus_forward", "particle_minus_topographic"),
  numerator = c("02_topographic_limited", "03_particle_topographic"),
  denominator = c("01_forward_geographic", "02_topographic_limited"),
  stringsAsFactors = FALSE
)

time_file <- function(paths, tm) {
  path <- unname(paths[[as.character(tm)]])
  if (is.null(path) || !length(path) || !file.exists(path)) {
    stop("No time table at ", tm, " Ma.", call. = FALSE)
  }
  read_time_metric(path)
}

diff_limits <- list()
for (cc in seq_len(nrow(comparison_specs))) {
  comparison <- comparison_specs$comparison_id[[cc]]
  max_abs <- stats::setNames(rep(0, length(comparison_metrics)), comparison_metrics)
  for (tm in times) {
    a <- time_file(scenario_time_files[[comparison_specs$numerator[[cc]]]], tm)
    b <- time_file(scenario_time_files[[comparison_specs$denominator[[cc]]]], tm)
    for (metric in comparison_metrics) {
      z <- merge(a[, c("cell_id", metric)], b[, c("cell_id", metric)], by = "cell_id",
                 suffixes = c("_a", "_b"), all = FALSE, sort = FALSE)
      delta <- z[[paste0(metric, "_a")]] - z[[paste0(metric, "_b")]]
      candidate <- suppressWarnings(max(abs(delta[is.finite(delta)]), na.rm = TRUE))
      if (is.finite(candidate)) max_abs[[metric]] <- max(max_abs[[metric]], candidate)
    }
  }
  for (metric in comparison_metrics) {
    diff_limits[[paste(comparison, metric, sep = "\r")]] <-
      max(max_abs[[metric]], .Machine$double.eps)
  }
}

write_difference_map <- function(df, metric, comparison_id, limits, out_dir) {
  safe_dir(out_dir)
  tm <- unique(df$time_ma)
  stem <- paste0(metric_map_code[[metric]], "_difference_", comparison_file_code[[comparison_id]],
                 "_", time_slug(tm))
  tif <- file.path(out_dir, paste0(stem, ".tif"))
  png <- file.path(out_dir, paste0(stem, ".png"))
  state_tif <- file.path(out_dir, paste0(stem, "_comparison_state.tif"))
  if (!all(c("comparison_state", "comparison_state_code") %in% names(df))) {
    stop("Difference-map data must include comparison_state and comparison_state_code.",
         call. = FALSE)
  }
  if (render_tif) {
    raster <- global_raster(df, "difference")
    terra::writeRaster(raster, tif, overwrite = TRUE,
                       gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
    # This companion raster makes the visual distinction machine readable:
    # 1 = both scenarios have a structural zero for this metric; 2 = equal,
    # non-zero values; 3 = a numerical scenario difference.  Cells not in the
    # raster remain NA and are ocean or unavailable palaeo-Earth cells.
    state_raster <- global_raster(df, "comparison_state_code")
    terra::writeRaster(state_raster, state_tif, overwrite = TRUE,
                       gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
  }
  if (render_png) {
    both_zero <- df[df$comparison_state == "both_zero", , drop = FALSE]
    equal_nonzero <- df[df$comparison_state == "equal_nonzero", , drop = FALSE]
    different <- df[df$comparison_state == "different", , drop = FALSE]
    # White is reserved for equality among non-zero values.  A medium-grey
    # land tile explicitly means both scenarios have zero process intensity,
    # avoiding the false impression that valid palaeo-grid cells were cut out.
    plot <- ggplot2::ggplot() +
      ggplot2::geom_tile(data = both_zero, ggplot2::aes(lon, lat),
                          fill = "grey68", width = 1, height = 1) +
      ggplot2::geom_tile(data = equal_nonzero, ggplot2::aes(lon, lat),
                          fill = "white", width = 1, height = 1) +
      ggplot2::geom_tile(data = different,
                          ggplot2::aes(lon, lat, fill = difference),
                          width = 1, height = 1) +
      ggplot2::coord_equal(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
      ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                                    midpoint = 0, limits = c(-limits, limits),
                                    oob = oob_squish, na.value = "grey55",
                                    name = "Difference") +
      ggplot2::labs(
        title = paste("Case05 dispersal-scheme difference:", metric),
        subtitle = paste0(comparison_id, "; ", tm,
                          " Ma; symmetric fixed scale across all comparison times"),
        x = "Longitude (degrees)", y = "Latitude (degrees)",
        caption = paste(
          "Red/blue: numerical scenario difference; white: equal non-zero values;",
          "medium grey: both scenarios are zero; dark grey: ocean or unavailable habitat."
        )
      ) +
      ggplot2::theme_minimal(base_size = 9) +
      ggplot2::theme(panel.grid = ggplot2::element_line(colour = "grey88", linewidth = 0.2),
                     panel.background = ggplot2::element_rect(fill = "grey55", colour = NA))
    ggplot2::ggsave(png, plot, width = 10, height = 5.5, dpi = map_dpi)
    if (!file.exists(png)) stop("Difference PNG map was not written: ", png,
                                 call. = FALSE)
  }
  data.frame(
    scenario_id = comparison_id, metric = paste0(metric, "_difference"),
    label = paste0(metric, " difference"), time_ma = tm, lower = -limits,
    upper = limits, tif = normalizePath(tif, winslash = "/", mustWork = FALSE),
    png = normalizePath(png, winslash = "/", mustWork = FALSE), crs = "EPSG:4326",
    resolution_deg = 1, source = "Case05 same-input dispersal comparison",
    interpretation = "symmetric_zero_centred_colour_scale_across_all_times",
    comparison_state_tif = normalizePath(state_tif, winslash = "/", mustWork = FALSE),
    both_zero_cells = nrow(df[df$comparison_state == "both_zero", , drop = FALSE]),
    equal_nonzero_cells = nrow(df[df$comparison_state == "equal_nonzero", , drop = FALSE]),
    different_cells = nrow(df[df$comparison_state == "different", , drop = FALSE]),
    stringsAsFactors = FALSE
  )
}

run_difference_map_job <- function(job) {
  a <- read_time_metric(job$path_a)
  b <- read_time_metric(job$path_b)
  merged <- merge(a[, c("cell_id", "lon", "lat", job$metric)],
                  b[, c("cell_id", job$metric)], by = "cell_id",
                  suffixes = c("_a", "_b"), all = FALSE, sort = FALSE)
  merged$time_ma <- job$time_ma
  a_value <- suppressWarnings(as.numeric(merged[[paste0(job$metric, "_a")]]))
  b_value <- suppressWarnings(as.numeric(merged[[paste0(job$metric, "_b")]]))
  merged$difference <- a_value - b_value
  tolerance <- 1e-15
  both_zero <- is.finite(a_value) & is.finite(b_value) &
    abs(a_value) <= tolerance & abs(b_value) <= tolerance
  equal_nonzero <- is.finite(merged$difference) &
    abs(merged$difference) <= tolerance & !both_zero
  merged$comparison_state <- ifelse(both_zero, "both_zero",
                                    ifelse(equal_nonzero, "equal_nonzero", "different"))
  merged$comparison_state_code <- ifelse(both_zero, 1L,
                                         ifelse(equal_nonzero, 2L, 3L))
  write_difference_map(
    merged[, c("lon", "lat", "time_ma", "difference", "comparison_state",
               "comparison_state_code"), drop = FALSE],
    job$metric, job$comparison_id, job$limit, job$out_dir
  )
}

if (render_difference_maps) {
  difference_jobs <- list()
  for (cc in seq_len(nrow(comparison_specs))) {
    comparison <- comparison_specs$comparison_id[[cc]]
    for (tm in map_times) {
      for (metric in comparison_metrics) {
        limit <- diff_limits[[paste(comparison, metric, sep = "\r")]]
        difference_jobs[[length(difference_jobs) + 1L]] <- list(
          path_a = unname(scenario_time_files[[comparison_specs$numerator[[cc]]]][[as.character(tm)]]),
          path_b = unname(scenario_time_files[[comparison_specs$denominator[[cc]]]][[as.character(tm)]]),
          metric = metric, time_ma = tm, comparison_id = comparison, limit = limit,
          out_dir = file.path(dirs$differences, comparison_map_folder[[comparison]],
                              metric_map_folder[[metric]])
        )
      }
    }
  }
  map_index <- c(map_index, run_map_jobs(difference_jobs, run_difference_map_job))
}

bind_fill <- function(x) {
  cols <- unique(unlist(lapply(x, names), use.names = FALSE))
  out <- lapply(x, function(z) {
    missing <- setdiff(cols, names(z))
    for (nm in missing) z[[nm]] <- NA
    z[, cols, drop = FALSE]
  })
  do.call(rbind, out)
}

map_index_df <- if (length(map_index)) {
  bind_fill(map_index)
} else {
  data.frame(
    scenario_id = character(), metric = character(), label = character(),
    time_ma = numeric(), lower = numeric(), upper = numeric(), tif = character(),
    png = character(), crs = character(), resolution_deg = numeric(),
    source = character(), interpretation = character(), stringsAsFactors = FALSE
  )
}
write_csv(map_index_df, map_index_path)

expected_regular_maps <- nrow(scenarios) * length(map_times) * length(metric_dirs)
expected_difference_maps <- nrow(comparison_specs) * length(map_times) *
  length(comparison_metrics)
regular_index <- map_index_df[map_index_df$scenario_id %in% scenarios$scenario_id,
                              , drop = FALSE]
regular_core_index <- regular_index[regular_index$metric %in% names(metric_dirs),
                                    , drop = FALSE]
difference_index <- map_index_df[map_index_df$scenario_id %in%
                                    comparison_specs$comparison_id, , drop = FALSE]
common_scale_ok <- if (nrow(regular_index)) {
  by_metric <- split(regular_index, regular_index$metric)
  all(vapply(by_metric, function(x) {
    length(unique(paste(x$lower, x$upper, sep = "\r"))) == 1L
  }, logical(1)))
} else NA
map_file_ok <- if (nrow(map_index_df)) {
  all(file.exists(map_index_df$tif)) && all(file.exists(map_index_df$png))
} else NA
state_file_ok <- if (nrow(difference_index) &&
                     "comparison_state_tif" %in% names(difference_index)) {
  all(!is.na(difference_index$comparison_state_tif) &
        file.exists(difference_index$comparison_state_tif))
} else NA
if (nrow(difference_index) &&
    all(c("both_zero_cells", "equal_nonzero_cells", "different_cells") %in%
        names(difference_index))) {
  difference_audit <- difference_index[, c(
    "scenario_id", "metric", "time_ma", "both_zero_cells",
    "equal_nonzero_cells", "different_cells", "comparison_state_tif"
  ), drop = FALSE]
  difference_audit$total_available_cells <- rowSums(
    difference_audit[, c("both_zero_cells", "equal_nonzero_cells", "different_cells"),
                     drop = FALSE], na.rm = TRUE
  )
  difference_audit$both_zero_pct <- ifelse(
    difference_audit$total_available_cells > 0,
    100 * difference_audit$both_zero_cells / difference_audit$total_available_cells,
    NA_real_
  )
  write_csv(difference_audit, file.path(dirs$summaries,
                                        "case05_difference_state_audit.csv"))
}
write_csv(data.frame(
  comparison_state_code = c(1L, 2L, 3L),
  comparison_state = c("both_zero", "equal_nonzero", "different"),
  interpretation = c(
    "Both scenarios have zero value in a valid palaeo-Earth cell; for colonisation this means no arrival-derived colonisation under either scenario.",
    "Both scenarios have the same finite non-zero value.",
    "The numerator and denominator scenarios differ numerically."
  ),
  stringsAsFactors = FALSE
), file.path(dirs$summaries, "case05_difference_state_codebook.csv"))
scenario_validation <- lapply(seq_len(nrow(scenarios)), function(i) {
  cfg_i <- read_csv(file.path(scenarios$output[[i]], "00_config",
                              "case04_final_config.csv"))
  val_i <- read_csv(file.path(scenarios$output[[i]], "00_config",
                              "case04_final_validation_checks.csv"))
  endpoint_score_path <- file.path(scenarios$output[[i]],
                                   "20_modern_endpoint_validation",
                                   "modern_endpoint_score_summary.csv")
  endpoint_gate_path <- file.path(scenarios$output[[i]],
                                  "20_modern_endpoint_validation",
                                  "modern_endpoint_calibration_gate.csv")
  endpoint_score <- if (file.exists(endpoint_score_path)) read_csv(endpoint_score_path) else data.frame()
  endpoint_gate <- if (file.exists(endpoint_gate_path)) read_csv(endpoint_gate_path) else data.frame()
  lookup <- function(x) val_i$status[match(x, val_i$check)]
  data.frame(
    scenario_id = scenarios$scenario_id[[i]],
    n_time_tables = length(metric_files[[scenarios$scenario_id[[i]]]]),
    n_projected_times = as.integer(cfg_i$value[match("n_projected_times", cfg_i$parameter)]),
    probability_range = lookup("probability_range"),
    global_grid_complete = lookup("global_grid_complete"),
    no_biogeobears = lookup("no_biogeobears"),
    no_discrete_regions = lookup("no_discrete_regions"),
    biotic_status = lookup("biotic_filtering_process"),
    speciation_status = lookup("speciation_event_table"),
    extinction_boundary = lookup("lineage_extinction_process_boundary"),
    endpoint_score_available = nrow(endpoint_score) == 2L,
    endpoint_gate_available = nrow(endpoint_gate) >= 1L,
    endpoint_data_role = if (nrow(endpoint_score)) endpoint_score$data_role[[1L]] else NA_character_,
    dynamic_endpoint_gate_pass = if (nrow(endpoint_gate)) {
      all(endpoint_gate$status[endpoint_gate$check != "endpoint_matrix_complete"] == "PASS")
    } else FALSE,
    stringsAsFactors = FALSE
  )
})
scenario_validation <- do.call(rbind, scenario_validation)
write_csv(scenario_validation, file.path(dirs$summaries,
                                         "case05_scenario_validation_detail.csv"))
required_process_columns <- c(
  "mean_environmental_support", "mean_arrival_hazard_per_myr",
  "mean_arrival_probability_reporting_interval",
  "mean_colonisation_hazard_per_myr",
  "mean_colonisation_probability_reporting_interval",
  "mean_persistence_probability", "expected_lineage_richness",
  "occupancy_weighted_diversity_defined"
)
validation <- data.frame(
  check = c("three_scenarios_complete", "same_time_axis_all_scenarios",
            "all_raw_time_tables_present", "probabilities_in_range",
             "complete_global_land_grid", "no_biogeobears_or_regions",
             "active_core_process_output_columns", "disabled_process_layers_not_mapped",
             "modern_endpoint_payload_and_score", "modern_endpoint_gate_training_diagnostic",
             "particle_kernel_recorded",
            "shared_legend_contract", "regular_map_count", "difference_map_count",
            "all_indexed_map_files_exist", "difference_zero_state_encoded",
            "arrival_zero_forces_colonisation_zero"),
  status = c(
    if (all(status_table$status %in% c("complete", "reused_complete"))) "PASS" else "FAIL",
    if (all(vapply(scenario_time_files, function(x) identical(sort(as.numeric(names(x)),
                                                                    decreasing = TRUE), times), logical(1)))) "PASS" else "FAIL",
    if (all(scenario_validation$n_time_tables == scenario_validation$n_projected_times)) "PASS" else "FAIL",
    if (all(scenario_validation$probability_range == "PASS")) "PASS" else "FAIL",
    if (all(scenario_validation$global_grid_complete == "PASS")) "PASS" else "FAIL",
    if (all(scenario_validation$no_biogeobears == "PASS") &&
        all(scenario_validation$no_discrete_regions == "PASS")) "PASS" else "FAIL",
     if (all(required_process_columns %in% names(first_table))) "PASS" else "FAIL",
     if (!any(map_index_df$metric %in% c("biotic_competition_pressure",
                                         "biotic_facilitation_pressure",
                                         "speciation_inheritance_footprint",
                                         "scenario_lineage_extinction_pressure"))) "PASS" else "FAIL",
     if (all(scenario_validation$endpoint_score_available) &&
         all(scenario_validation$endpoint_gate_available)) "PASS" else "FAIL",
     if (all(scenario_validation$endpoint_data_role == "training_diagnostic")) {
       "BOUNDARY_TRAINING_DIAGNOSTIC_ONLY"
     } else if (all(scenario_validation$dynamic_endpoint_gate_pass)) "PASS" else "FAIL",
     if (file.exists(file.path(scenarios$output[[3L]], "11_dispersal",
                              "movement_matrix_summary.csv"))) "PASS" else "FAIL",
    if (is.na(common_scale_ok)) "NOT_RUN" else if (common_scale_ok) "PASS" else "FAIL",
     if (!render_regular_maps) "NOT_RUN" else if (nrow(regular_core_index) == expected_regular_maps) "PASS" else "FAIL",
    if (!render_difference_maps) "NOT_RUN" else if (nrow(difference_index) == expected_difference_maps) "PASS" else "FAIL",
    if (is.na(map_file_ok)) "NOT_RUN" else if (map_file_ok) "PASS" else "FAIL",
    if (is.na(state_file_ok)) "NOT_RUN" else if (state_file_ok) "PASS" else "FAIL",
    if (all(colonisation_zero_audit$nonzero_colonisation_when_arrival_zero == 0L)) "PASS" else "FAIL"
  ),
  detail = c(
    "All three controlled P2 scenarios completed or were verified for reuse.",
    "Each scenario uses exactly the same tree-domain palaeoenvironment time values.",
    "A cell_metrics table exists for every projected time slice in every scenario.",
    "Engine validation confirms occupancy probabilities are in [0,1].",
    "Engine validation confirms every valid land/habitat cell is included in its time slice.",
    "The no-BioGeoBEARS no-discrete-region model boundary is enforced by the engine.",
     "P1 environmental support, P2 arrival/colonisation, persistence and probability-weighted diversity are explicit output columns. The intercept-free environmental-support sum is not labelled as modern expected richness.",
     "P3 biotic filtering, demonstration speciation footprint and demonstration lineage-extinction pressure are OFF in empirical core and therefore omitted from the formal map atlas.",
     "Every scenario preserves a merged full-grid 0 Ma dynamic endpoint and HMSC fixed-effect nowcast, with an auditable modern score table.",
     "The current supplied observations are labelled training_diagnostic; a spatial holdout or independent distribution source is required before the gate is evidence of out-of-sample calibration.",
     "Particle-topographic movement summary records the finite-propagule kernel and seed.",
    "The shared per-metric scale file defines one legend range across all regular maps.",
    paste0("Expected and indexed regular maps: ", expected_regular_maps, "."),
    paste0("Expected and indexed difference maps: ", expected_difference_maps, "."),
    "Every row in the map index resolves to both a GeoTIFF and a PNG.",
    "Every difference map has a categorical companion GeoTIFF distinguishing both-zero, equal-nonzero, and different cells.",
    "Across every scenario and time slice, zero P2 arrival always forces zero colonisation. A zero colonisation value with finite arrival is retained as an establishment-limited or numerical-negligibility state, not mislabelled as source limitation."
  ),
  stringsAsFactors = FALSE
)
write_csv(validation, file.path(dirs$summaries, "case05_final_validation.csv"))

time_summaries <- lapply(seq_len(nrow(scenarios)), function(i) {
  x <- read_csv(file.path(scenarios$output[[i]], "22_report", "case04_final_time_summary.csv"))
  x$scenario_id <- scenarios$scenario_id[[i]]
  x$dispersal_scheme <- scenarios$dispersal_scheme[[i]]
  x
})
time_summaries <- do.call(rbind, time_summaries)
write_csv(time_summaries, file.path(dirs$summaries, "case05_global_time_scheme_summary.csv"))

metadata <- data.frame(
  output_type = c("input", "HMSC_posterior", "Evolution", "Environmental_filtering",
                  "Dispersal", "Colonisation", "Biotic_filtering", "Persistence",
                  "Speciation", "Extinction", "model_boundaries", "uncertainty", "maps"),
  source_or_rule = c(
    "Plant200 comm/sites/traits/dated tree and 1 degree palaeo-Earth grid",
    "Existing modern HMSC posterior Beta response draws; no posterior maps are averaged backward",
    "Dated-tree direct ancestral response reconstruction",
    "Palaeoclimate plus ancestral beta gives environmental support, not final occupancy",
    "Three named dispersal scenarios; same local scale and movement intercept",
    "Arrival hazard multiplied by conditional establishment in the global-grid CTMC",
    "Scenario-based explicit competition/facilitation, not an HMSC residual-association interpretation",
    "Local persistence within the CTMC update",
    "Dated-tree lineage identity and copy-then-diverge inheritance",
    "Local/geography-forced loss plus a clearly labelled lineage-extinction pressure scenario",
    "Not used: no BioGeoBEARS, BSM, region gate, A_hist, or discrete regional accessibility",
    "No Q_extrap multiplicative absence penalty; extrapolation remains an uncertainty context",
    "All three schemes and all time slices share per-metric legend limits; differences use symmetric zero-centred limits"
  ),
  stringsAsFactors = FALSE
)
write_csv(metadata, file.path(dirs$summaries, "case05_interpretation_metadata.csv"))

report <- c(
  "# Case05 Plant200 empirical-core, three-dispersal comparison",
  "",
  "This is a global 1-degree, no-BioGeoBEARS case. It carries the sampled extant Plant200 tree backward through dated ancestral branches and forward through the palaeo-Earth grid.",
  "",
  "## Controlled comparison",
  "",
  "All three runs use identical HMSC response draws, dated tree, ancestral-response reconstruction, complete palaeo-land grid, PALEOMAP target-centred plate carriage, compact root prior, colonisation, persistence, seed, and internal time step. Only P2 dispersal differs. All actual time slices are integrated; only the configured representative time slices are rendered as maps.",
  "",
  "- `forward_geographic`: land-only 8-neighbour movement with great-circle edge distance, no relief or elevation-step penalty.",
  "- `topographic_limited`: same graph, with effective edge cost equal to geographic distance multiplied by continuous local-relief/elevation-step resistance.",
  "- `particle_topographic`: the same topographic kernel after finite-propagule Monte Carlo routing. It is not a posterior particle filter and not a gene-flow estimate.",
  "",
  "## Scientific boundaries",
  "",
  "Environmental support is not historical occupancy. Biotic filtering is OFF because no independent interaction data are available, and complete lineage extinction is not estimated from the extant-only tree. The PALEOMAP step is target-centred nearest-track plate carriage, not exact polygon overlap. Outputs are sampled-surviving-lineage forward scenarios, not the total global flora or a unique reconstructed history.",
  "",
  "## Reading maps",
  "",
  "Every mapped metric has one fixed colour range across all dynamic time slices and all three dispersal scenarios. Difference maps are numerator minus denominator, with red positive and blue negative values on a symmetric zero-centred scale. Pale yellow on regular maps is a valid zero; dark grey is ocean or unavailable habitat.",
  "On difference maps, medium grey means both scenarios have a valid grid cell but a zero value. For colonisation, consult the paired arrival map and zero audit before calling that state source-limited: low establishment can also yield a zero/numerically negligible colonisation value. Dark grey is ocean or unavailable habitat, and white denotes equal non-zero values. The categorical companion GeoTIFF uses codes 1=both_zero, 2=equal_nonzero, 3=different.",
  "",
  "## Core files",
  "",
  "- `00_config/case05_eight_process_contract.csv`: process roles and empirical/scenario status.",
  "- `00_config/case05_shared_metric_limits.csv`: the common legend contract.",
  "- `02_scenario_runs/`: raw, reproducible Case04 engine output for each scheme.",
  "- `03_shared_scale_maps/`: GIS GeoTIFF and PNG maps with matched legend scales.",
  "- `04_scheme_difference_maps/`: fixed zero-centred difference maps.",
  "- `05_summaries/case05_colonisation_zero_arrival_audit.csv`: verifies that zero arrival always forces zero colonisation and separates finite-arrival, establishment-limited zeros.",
  "- `05_summaries/case05_difference_state_audit.csv`: quantifies both-zero, equal-nonzero and different cells for every comparison map.",
  "- `05_summaries/case05_global_time_scheme_summary.csv`: global time trajectories.",
  ""
)
writeLines(report, file.path(dirs$report, "case05_readme.md"), useBytes = TRUE)

message("Case05 complete: ", output)
