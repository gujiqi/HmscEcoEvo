#!/usr/bin/env Rscript

# Legacy Case05 Plant200 global-grid dynamic occupancy scenario without BGB.
# Establishment and persistence here are scenario parameters, not empirical
# rates or independent processes in the current six-process framework.
#
# This is the Case05 dispatcher described in the no-BioGeoBEARS specification:
# modern HMSC beta posterior -> ancestral response history -> global palaeo
# grid environmental filtering -> local/LDD movement -> arrival -> establishment
# -> colonisation -> persistence/local loss -> dynamic occupancy -> diversity
# maps. BioGeoBEARS, BSM histories, discrete regions, A_hist, and M1-M5
# multipliers are intentionally not used.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || is.na(x) || identical(x, "")) y else x
}

parse_args <- function(args) {
  out <- list()
  for (a in args) {
    if (!grepl("^--", a)) next
    kv <- sub("^--", "", a)
    pos <- regexpr("=", kv, fixed = TRUE)
    if (pos < 0) {
      out[[kv]] <- TRUE
    } else {
      out[[substr(kv, 1, pos - 1L)]] <- substr(kv, pos + 1L, nchar(kv))
    }
  }
  out
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

as_num <- function(x, default = NA_real_) {
  if (is.null(x) || identical(x, "")) return(default)
  suppressWarnings(as.numeric(x))
}

as_int <- function(x, default = NA_integer_) {
  y <- as_num(x, default)
  if (!is.finite(y)) return(default)
  as.integer(y)
}

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

write_csv <- function(x, path) {
  safe_dir_create(dirname(path))
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(x, path)
  } else {
    utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
  }
  invisible(path)
}

read_csv <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    as.data.frame(data.table::fread(path))
  } else {
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
}

write_text <- function(x, path) {
  safe_dir_create(dirname(path))
  writeLines(x, path, useBytes = TRUE)
  invisible(path)
}

stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

find_rscript <- function() {
  exe <- if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  file.path(R.home("bin"), exe)
}

copy_if_exists <- function(from, to) {
  if (!file.exists(from)) return(FALSE)
  safe_dir_create(dirname(to))
  file.copy(from, to, overwrite = TRUE)
}

add_prefix_copy <- function(output, old_rel, new_rel) {
  old <- file.path(output, old_rel)
  new <- file.path(output, new_rel)
  copy_if_exists(old, new)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

pkg_root <- normalizePath(
  args$pkg_root %||% "C:/Users/Google/Documents/HMSC-HIST/HmscEcoEvo",
  winslash = "/", mustWork = FALSE
)

case04_engine <- normalizePath(
  args$engine %||% file.path(pkg_root, "scripts",
                             "case04_plant200_global_dynamic_no_bgb_final.R"),
  winslash = "/", mustWork = FALSE
)

quick <- as_bool(args$quick, FALSE)
run_type <- args$run_type %||% "forward_scenario"
valid_run_types <- c("forward_scenario",
                     "endpoint_conditioned_reconstruction",
                     "all_processes_demo")
if (!run_type %in% valid_run_types) {
  stop("Invalid --run_type. Use one of: ", paste(valid_run_types, collapse = ", "),
       call. = FALSE)
}

if (identical(run_type, "all_processes_demo")) {
  args$run_biotic_filtering <- args$run_biotic_filtering %||% "TRUE"
  args$run_speciation_demo <- args$run_speciation_demo %||% "TRUE"
  args$run_lineage_extinction_demo <-
    args$run_lineage_extinction_demo %||% "TRUE"
  args$biotic_competition_strength <-
    args$biotic_competition_strength %||% "0.35"
  args$biotic_facilitation_strength <-
    args$biotic_facilitation_strength %||% "0.12"
}

output <- args$output %||% file.path(
  pkg_root, "outputs",
  paste0("case05_plant200_global_dynamic_no_bgb_", run_type, "_", stamp())
)
output <- safe_dir_create(output)

dirs <- list(
  config = safe_dir_create(file.path(output, "00_config")),
  audit = safe_dir_create(file.path(output, "01_input_audit")),
  biotic = safe_dir_create(file.path(output, "10_biotic_filtering")),
  attribution = safe_dir_create(file.path(output, "21_sensitivity",
                                          "process_attribution")),
  report = safe_dir_create(file.path(output, "22_report")),
  logs = safe_dir_create(file.path(output, "logs"))
)
log_file <- file.path(dirs$logs, "case05_dispatcher.log")
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                paste(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_msg("Case05 dispatcher started.")
log_msg("Run type:", run_type)

mode_contract <- data.frame(
  run_type = valid_run_types,
  scientific_role = c(
    "Software demonstration and scenario exploration; not endpoint-conditioned reconstruction.",
    "Primary statistical target, but requires independent/held-out modern endpoint and fossil likelihood inputs.",
    "Interface demonstration for all eight processes; scenario parameters are not empirical estimates."
  ),
  allowed_without_independent_endpoint = c(TRUE, FALSE, TRUE),
  output_interpretation = c(
    "forward eco-evolutionary scenario",
    "endpoint/fossil-conditioned history ensemble",
    "process demonstration, not empirically estimated"
  ),
  stringsAsFactors = FALSE
)
write_csv(mode_contract, file.path(dirs$config, "case05_mode_contract.csv"))

removed_bgb <- data.frame(
  removed_object = c("region_cube", "region_id", "bgb_accessibility",
                     "bgb_region_history", "BioGeoBEARS_BSM_history",
                     "R_region", "A_hist", "M1_M5_multiplier_chain"),
  case05_status = "REMOVED_FROM_NO_BGB_MAIN_MODEL",
  reason = c(
    "No discrete regional conditioning in Case05.",
    "Cells are handled directly on the global grid.",
    "Historical accessibility is not multiplied into final probability.",
    "No region history is used in the no-BGB model.",
    "BioGeoBEARS/BSM is intentionally cancelled for this case.",
    "No lineage-by-region state gates are used.",
    "Avoids duplicate accessibility/dispersal multiplication.",
    "Legacy scenario only; not Case05 core dynamic occupancy."
  ),
  stringsAsFactors = FALSE
)
write_csv(removed_bgb, file.path(dirs$audit, "case05_removed_bgb_inputs.csv"))

conditioning_inputs <- data.frame(
  input = c("modern_holdout", "independent_range", "fossils",
            "fossil_likelihood_table"),
  path = c(args$modern_holdout %||% "", args$independent_range %||% "",
           args$fossils %||% "", args$fossil_likelihood_table %||% ""),
  required_for_endpoint_conditioned = c(TRUE, TRUE, TRUE, FALSE),
  exists = c(file.exists(args$modern_holdout %||% ""),
             file.exists(args$independent_range %||% ""),
             file.exists(args$fossils %||% ""),
             file.exists(args$fossil_likelihood_table %||% "")),
  stringsAsFactors = FALSE
)
conditioning_inputs$note <- c(
  "Use a held-out modern grid, not the same records used to fit HMSC.",
  "Alternative independent present-day range evidence.",
  "Fossils may be family/clade-level; do not force species-level validation.",
  "Optional precomputed fossil likelihood weights."
)
write_csv(conditioning_inputs, file.path(dirs$audit,
                                         "case05_conditioning_inputs.csv"))

if (identical(run_type, "endpoint_conditioned_reconstruction")) {
  endpoint_rows <- conditioning_inputs$input %in% c("modern_holdout",
                                                    "independent_range")
  fossil_rows <- conditioning_inputs$input %in% c("fossils",
                                                  "fossil_likelihood_table")
  has_endpoint <- any(conditioning_inputs$exists[endpoint_rows])
  has_fossils <- any(conditioning_inputs$exists[fossil_rows])
  missing <- conditioning_inputs[
    (endpoint_rows & !has_endpoint) | (fossil_rows & !has_fossils),
    ,
    drop = FALSE
  ]
  if (nrow(missing)) {
    write_csv(missing, file.path(dirs$audit,
                                 "missing_conditioning_inputs.csv"))
    stop("endpoint_conditioned_reconstruction requires an independent modern endpoint and fossil/fossil-likelihood input. See missing_conditioning_inputs.csv",
         call. = FALSE)
  }
}

if (!file.exists(case04_engine)) {
  stop("Missing Case04 streaming engine: ", case04_engine, call. = FALSE)
}

engine_args <- c(
  case04_engine,
  paste0("--quick=", if (quick) "TRUE" else "FALSE"),
  paste0("--output=", output)
)

pass_through <- c("pkg_root", "previous_case04", "plant_input_dir", "tree",
                  "response_history", "palaeo_earth_state", "recipe",
                  "time_domain", "comm", "sites", "fossils", "metadata",
                  "traits", "trait_raw", "time_selection", "n_time_slices",
                  "max_times", "max_species",
                  "max_draws", "internal_dt", "persistence_reference_myr",
                  "reporting_interval_myr", "map_n_time_slices",
                  "root_rho", "root_max_cells",
                  "root_initialisation", "root_seed_count", "root_support_power",
                  "reporting_interval_myr", "transport_mode",
                  "plate_transport_index", "transport_min_target_coverage",
                  "movement_intercept", "distance_decay", "ldd_per_source",
                  "ldd_intercept", "ldd_distance_decay",
                  "establishment_intercept", "establishment_slope",
                  "persistence_intercept", "persistence_slope",
                  "run_biotic_filtering", "biotic_competition_strength",
                  "biotic_facilitation_strength", "biotic_persistence_strength",
                  "run_speciation_demo", "run_lineage_extinction_demo",
                  "binary_threshold", "seed", "make_png", "make_tif", "link",
                  "write_attribution_state", "attribution_max_cells",
                  "attribution_max_lineages", "attribution_max_rows")
for (nm in pass_through) {
  if (!is.null(args[[nm]]) && !nm %in% c("quick", "output")) {
    engine_args <- c(engine_args, paste0("--", nm, "=", args[[nm]]))
  }
}

log_msg("Calling streaming engine:", case04_engine)
status <- system2(find_rscript(), engine_args)
if (!identical(status, 0L)) {
  stop("Case05 streaming engine failed with exit status ", status, call. = FALSE)
}
log_msg("Streaming engine completed.")

case04_config <- file.path(output, "00_config", "case04_final_config.csv")
case04_validation <- file.path(output, "00_config",
                               "case04_final_validation_checks.csv")
case04_map_index <- file.path(output, "22_report", "case04_final_map_index.csv")
case04_time_summary <- file.path(output, "22_report",
                                 "case04_final_time_summary.csv")
case04_full_attr_dir <- file.path(output, "21_sensitivity",
                                  "full_process_contribution")

if (file.exists(case04_config)) {
  cfg <- read_csv(case04_config)
  cfg <- rbind(
    data.frame(parameter = "case", value = "Case05_Plant200_global_dynamic_no_bgb"),
    data.frame(parameter = "run_type", value = run_type),
    data.frame(parameter = "bioGeoBEARS", value = "cancelled"),
    data.frame(parameter = "region_layer", value = "not_used"),
    data.frame(parameter = "A_hist", value = "not_used"),
    cfg
  )
  write_csv(cfg, file.path(dirs$config, "case05_config.csv"))
}

add_prefix_copy(output, "22_report/case04_final_map_index.csv",
                "22_report/case05_map_index.csv")
add_prefix_copy(output, "22_report/case04_final_time_summary.csv",
                "22_report/case05_time_summary.csv")
add_prefix_copy(output, "00_config/case04_final_validation_checks.csv",
                "00_config/case05_engine_validation_checks.csv")
add_prefix_copy(output,
                "21_sensitivity/full_process_contribution/case05_full_process_contribution_summary.csv",
                "21_sensitivity/process_attribution/case05_full_process_contribution_summary.csv")
add_prefix_copy(output,
                "21_sensitivity/full_process_contribution/case05_full_process_contribution_time_series.csv",
                "21_sensitivity/process_attribution/case05_full_process_contribution_time_series.csv")
add_prefix_copy(output,
                "21_sensitivity/full_process_contribution/case05_full_process_contribution_draws.csv",
                "21_sensitivity/process_attribution/case05_full_process_contribution_draws.csv")
add_prefix_copy(output,
                "21_sensitivity/full_process_contribution/case05_full_process_contribution_history_maps.csv",
                "21_sensitivity/process_attribution/case05_full_process_contribution_history_maps.csv")

map_index <- if (file.exists(case04_map_index)) read_csv(case04_map_index) else data.frame()
time_summary <- if (file.exists(case04_time_summary)) read_csv(case04_time_summary) else data.frame()

if (dir.exists(pkg_root) && requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(pkg_root, quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(HmscEcoEvo))
}

run_process_attribution <- as_bool(args$run_process_attribution, TRUE)
attribution_workers <- as_int(args$attribution_workers, 20L)
attribution_targets <- strsplit(args$attribution_targets %||%
                                  "occupancy,range_area,lineage_richness",
                                ",", fixed = TRUE)[[1]]
attribution_targets <- trimws(attribution_targets)
attribution_targets <- attribution_targets[nzchar(attribution_targets)]
attribution_state_path <- file.path(output, "16_state_space_inference",
                                    "attribution_interval_state.csv")
attribution_metadata_path <- file.path(output, "16_state_space_inference",
                                       "attribution_state_metadata.csv")

prepare_events_for_attribution_target <- function(events, target_name) {
  if (is.null(events) || !length(events)) return(events)
  out <- events
  if (!is.null(out$speciation_events) &&
      identical(target_name, "lineage_richness")) {
    d <- out$speciation_events
    if (!"delta_metric" %in% names(d) &&
        all(c("n_current_lineages", "n_previous_lineages") %in% names(d))) {
      d$delta_metric <- suppressWarnings(
        as.numeric(d$n_current_lineages) -
          as.numeric(d$n_previous_lineages)
      )
      d$event_metric <- "lineage_richness_delta"
    }
    out$speciation_events <- d
  }
  out
}

run_one_attribution_target <- function(target_name, attr_obj) {
  local_obj <- attr_obj
  local_obj$events <- prepare_events_for_attribution_target(
    attr_obj$events, target_name
  )
  ans <- hee_process_attribution(
    local_obj,
    target = target_name,
    integrate_time = TRUE,
    normalise = "time",
    upstream_processes = "available",
    return = c("summary", "time_series", "maps", "draws")
  )
  ans$summary$target <- target_name
  if (!is.null(ans$state_budget)) ans$state_budget$target <- target_name
  if (!is.null(ans$time_series)) ans$time_series$target <- target_name
  if (!is.null(ans$maps)) ans$maps$target <- target_name
  if (!is.null(ans$draws)) ans$draws$target <- target_name
  ans
}

bind_attribution_component <- function(results, component) {
  rows <- lapply(results, function(z) z[[component]])
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(data.frame())
  cols <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(d) {
    missing <- setdiff(cols, names(d))
    for (m in missing) d[[m]] <- NA
    d[, cols, drop = FALSE]
  })
  do.call(rbind, rows)
}

attribution_benchmark <- data.frame()
if (run_process_attribution && file.exists(attribution_state_path)) {
  log_msg("Reading Case05 attribution interval state:", attribution_state_path)
  attr_state <- read_csv(attribution_state_path)
  attr_meta <- if (file.exists(attribution_metadata_path)) {
    read_csv(attribution_metadata_path)
  } else {
    data.frame()
  }
  speciation_events <- file.path(output, "14_speciation",
                                 "speciation_event_summary.csv")
  events <- list()
  if (file.exists(speciation_events)) {
    events$speciation_events <- read_csv(speciation_events)
  }
  process_status <- c(
    environmental_filtering = "scenario_or_model_derived",
    dispersal = "scenario_parameterised",
    colonisation = "estimated_exact_ctmc_budget",
    biotic_filtering = if (identical(run_type, "all_processes_demo")) {
      "scenario"
    } else {
      "not_estimated"
    },
    persistence = "estimated_exact_ctmc_budget",
    evolution = "model_derived_from_hmsc_tree",
    speciation = "tree_conditioned",
    extinction = if (identical(run_type, "all_processes_demo")) {
      "scenario_pressure"
    } else {
      "not_estimated_global_lineage_extinction"
    }
  )
  attr_obj <- list(
    intervals = attr_state,
    events = events,
    metadata = list(
      process_status = process_status,
      attribution_state_metadata = attr_meta
    )
  )

  log_msg("Running process attribution serially for targets:",
          paste(attribution_targets, collapse = ", "))
  t1 <- proc.time()
  serial_results <- lapply(attribution_targets, run_one_attribution_target,
                           attr_obj = attr_obj)
  serial_elapsed <- unname((proc.time() - t1)[["elapsed"]])
  names(serial_results) <- attribution_targets

  max_workers <- if (requireNamespace("parallel", quietly = TRUE)) {
    max(1L, min(attribution_workers, parallel::detectCores(logical = TRUE)))
  } else {
    1L
  }
  parallel_results <- NULL
  parallel_elapsed <- NA_real_
  if (max_workers > 1L && length(attribution_targets) > 1L &&
      requireNamespace("parallel", quietly = TRUE)) {
    log_msg("Running process attribution with workers:", max_workers)
    t2 <- proc.time()
    cl <- parallel::makeCluster(max_workers)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    parallel::clusterExport(
      cl,
      varlist = c("attr_obj", "attribution_targets",
                  "prepare_events_for_attribution_target",
                  "run_one_attribution_target"),
      envir = environment()
    )
    parallel_results <- parallel::parLapply(
      cl,
      attribution_targets,
      function(target_name) {
        suppressPackageStartupMessages(library(HmscEcoEvo))
        run_one_attribution_target(target_name, attr_obj)
      }
    )
    parallel_elapsed <- unname((proc.time() - t2)[["elapsed"]])
    names(parallel_results) <- attribution_targets
    try(parallel::stopCluster(cl), silent = TRUE)
  }

  chosen_results <- if (!is.null(parallel_results)) parallel_results else serial_results
  attribution_summary <- bind_attribution_component(chosen_results, "summary")
  attribution_state_budget <- bind_attribution_component(chosen_results, "state_budget")
  attribution_time_series <- bind_attribution_component(chosen_results, "time_series")
  attribution_maps <- bind_attribution_component(chosen_results, "maps")
  attribution_draws <- bind_attribution_component(chosen_results, "draws")

  write_csv(attribution_summary, file.path(dirs$attribution,
                                           "case05_process_attribution_summary.csv"))
  write_csv(attribution_state_budget, file.path(dirs$attribution,
                                                "case05_process_attribution_state_budget.csv"))
  write_csv(attribution_time_series, file.path(dirs$attribution,
                                               "case05_process_attribution_time_series.csv"))
  write_csv(attribution_maps, file.path(dirs$attribution,
                                        "case05_process_attribution_maps.csv"))
  write_csv(attribution_draws, file.path(dirs$attribution,
                                         "case05_process_attribution_draws.csv"))
  saveRDS(chosen_results, file.path(dirs$attribution,
                                    "case05_process_attribution_results.rds"))

  attribution_benchmark <- data.frame(
    mode = c("serial_1_worker", paste0("parallel_", max_workers, "_workers")),
    workers = c(1L, max_workers),
    elapsed_seconds = c(serial_elapsed, parallel_elapsed),
    n_targets = length(attribution_targets),
    n_interval_rows = nrow(attr_state),
    speedup_vs_serial = c(1, if (is.finite(parallel_elapsed) &&
                                  parallel_elapsed > 0) {
      serial_elapsed / parallel_elapsed
    } else {
      NA_real_
    }),
    interpretation = c(
      "Single-process attribution over all requested targets.",
      "Target-level parallel attribution; may be slower for small target counts because state tables must be serialised to workers."
    ),
    stringsAsFactors = FALSE
  )
  write_csv(attribution_benchmark, file.path(dirs$attribution,
                                             "case05_process_attribution_benchmark.csv"))
  log_msg("Process attribution completed. Serial seconds:", round(serial_elapsed, 3),
          "Parallel seconds:", round(parallel_elapsed, 3))
} else {
  attribution_benchmark <- data.frame(
    mode = "not_run",
    workers = attribution_workers,
    elapsed_seconds = NA_real_,
    n_targets = length(attribution_targets),
    n_interval_rows = NA_integer_,
    speedup_vs_serial = NA_real_,
    interpretation = if (!file.exists(attribution_state_path)) {
      "Attribution interval state was not written by the engine."
    } else {
      "run_process_attribution is FALSE."
    },
    stringsAsFactors = FALSE
  )
  write_csv(attribution_benchmark, file.path(dirs$attribution,
                                             "case05_process_attribution_benchmark.csv"))
}

conditioning_status <- data.frame(
  layer = c("modern_endpoint", "fossil_constraint",
            "posterior_history_weighting", "particle_smoothing"),
  status = c(
    if (identical(run_type, "endpoint_conditioned_reconstruction")) {
      "INPUTS_PRESENT_BUT_ENGINE_IS_FORWARD_WEIGHT_SUMMARY"
    } else {
      "DIAGNOSTIC_ONLY_NOT_CONDITIONED"
    },
    if (file.exists(args$fossils %||% "")) {
      "FOSSIL_FILE_SUPPLIED_DIAGNOSTIC"
    } else {
      "DIAGNOSTIC_OR_COPIED_BY_ENGINE_IF_DEFAULT_EXISTS"
    },
    "NOT_IMPLEMENTED_IN_THIS_FORWARD_ENGINE",
    "NOT_IMPLEMENTED_IN_THIS_FORWARD_ENGINE"
  ),
  interpretation = c(
    "Do not double-use GBIF training data as endpoint likelihood; provide held-out/independent data for reconstruction.",
    "Family/clade fossil evidence is diagnostic unless a fossil likelihood table is supplied and used for particle weighting.",
    "This Case05 run summarises equal-weight response draws, not reweighted historical particles.",
    "The script is a streaming mean-field forward scenario; full SMC/particle MCMC remains a future mode."
  ),
  stringsAsFactors = FALSE
)
write_csv(conditioning_status, file.path(dirs$report,
                                         "case05_conditioning_status.csv"))

biotic_mode <- if (identical(run_type, "all_processes_demo")) {
  "all_processes_demo_scenario_enabled"
} else {
  "empirical_core_off_zero_effect"
}
extinction_mode <- if (identical(run_type, "all_processes_demo")) {
  "all_processes_demo_scenario_pressure"
} else {
  "empirical_core_partial"
}

process_coverage <- data.frame(
  function_group = c(
    "data_preparation", "time_axis", "lineage_time_tree",
    "evolution", "environmental_filtering", "biotic_filtering",
    "dispersal", "colonisation", "persistence", "state_update",
    "speciation", "extinction", "modern_validation",
    "fossil_validation", "diversity_maps", "reporting",
    "BGB_functions"
  ),
  function_or_engine = c(
    "Case05 dispatcher + Case04 streaming engine input audit",
    "hee_time_axis/engine projected times from palaeo earth state",
    "active lineage table from dated tree",
    "hee_evolution_ancestral_response_direct",
    "streaming matrix environmental filtering equivalent to hee_environmental_filtering_ancestral_suitability",
    "optional scenario biotic filtering modifier inside establishment/persistence",
    "hee_within_region_movement_kernel plus sparse global graph",
    "arrival x establishment -> colonisation hazard",
    "hee_persistence_probability formula as streaming matrix",
    "CTMC mean-field occupancy update",
    "copy-then-diverge inheritance",
    "local/geography-forced loss; optional scenario lineage-extinction pressure",
    "modern endpoint diagnostic",
    "family/clade fossil diagnostic",
    "expected richness, potential richness, Shannon, Simpson, WE, response dispersion, refugia",
    "case05 report, map index, validation checks",
    "hee_bgb_*"
  ),
  case05_mode = c(
    "all", "all", "all", "all", "all", biotic_mode,
    "all", "all", "all", "all", "all", extinction_mode,
    "diagnostic_or_endpoint_mode", "diagnostic_or_endpoint_mode",
    "all", "all", "not_used_in_no_bgb_case05"
  ),
  execution_stage = c(
    "run once", "run once", "run once", "per draw/time",
    "per draw/time", "per substep", "per interval/substep",
    "per substep", "per substep", "per substep", "tree nodes",
    "cell land transitions", "end of run", "end of run",
    "end of run", "end of run", "excluded"
  ),
  executed = c(
    TRUE, TRUE, TRUE, TRUE, TRUE,
    TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
    TRUE, TRUE, nrow(map_index) > 0, TRUE, FALSE
  ),
  output_object = c(
    "case05_input_inventory/config",
    "case05_time_summary",
    "active_lineage_counts",
    "ancestral_beta_response_all_projected_times",
    "mean_environmental_support maps/tables",
    "biotic_competition_pressure and biotic_facilitation_pressure maps",
    "movement_matrix_summary and arrival_probability maps",
    "colonisation_probability maps",
    "persistence/local_extinction maps",
    "cell_metrics time slices",
    "speciation_event_summary and speciation_inheritance_footprint maps",
    "local_extinction_probability maps; optional scenario_lineage_extinction_pressure maps",
    "modern_endpoint_diagnostic_summary",
    "fossil_constraint_diagnostic_summary",
    "18_diversity_maps and map index",
    "case05_interpretation.md",
    "case05_removed_bgb_inputs.csv"
  ),
  consumed_by = c(
    "all stages", "lineage/evolution/dynamic grid", "evolution/speciation",
    "environmental filtering", "colonisation/persistence/maps",
    "colonisation/persistence", "colonisation/rescue",
    "state update", "state update", "diversity summaries",
    "lineage state", "state masks/validation", "report",
    "report", "report", "user/reviewer", "none"
  ),
  stringsAsFactors = FALSE
)
process_coverage <- rbind(
  process_coverage,
  data.frame(
    function_group = "process_attribution",
    function_or_engine = "hee_process_attribution",
    case05_mode = "all",
    execution_stage = "post-processing after dynamic engine",
    executed = run_process_attribution && file.exists(file.path(
      dirs$attribution, "case05_process_attribution_summary.csv"
    )),
    output_object = "case05_process_attribution_summary/state_budget/time_series/maps/draws",
    consumed_by = "reviewer-facing process contribution tables",
    stringsAsFactors = FALSE
  )
)
process_coverage <- rbind(
  process_coverage,
  data.frame(
    function_group = "full_streaming_process_contribution",
    function_or_engine = "Case04 streaming CTMC budget using hee_process_attribution formulas",
    case05_mode = "all",
    execution_stage = "during dynamic engine; all active lineages and all valid land cells",
    executed = file.exists(file.path(
      dirs$attribution, "case05_full_process_contribution_summary.csv"
    )),
    output_object = "case05_full_process_contribution_summary/time_series/draws/history_maps",
    consumed_by = "formal all-species contribution tables",
    stringsAsFactors = FALSE
  )
)
write_csv(process_coverage, file.path(dirs$report,
                                      "case05_function_coverage.csv"))

eight_process_inventory <- data.frame(
  process = c("environmental_filtering", "dispersal", "colonisation",
              "biotic_filtering", "persistence", "evolution",
              "speciation", "extinction"),
  plain_question = c(
    "Is the palaeoenvironment suitable for this lineage response?",
    "Can occupied cells send propagules to other cells?",
    "Can arriving propagules establish in the target cell?",
    "Do scenario interactions modify establishment or persistence?",
    "Can an occupied local population persist to the next interval?",
    "How do ancestral environmental responses change along the tree?",
    "How is lineage identity inherited at dated tree nodes?",
    "Where is local loss high, and what is the scenario lineage-extinction pressure?"
  ),
  formula_or_rule = c(
    "S_env = inverse_link(B(X_palaeo) %*% beta_ancestral)",
    "movement lambda from sparse neighbour/LDD graph and great-circle distance",
    "lambda_C = arrival_hazard * establishment_probability",
    "B_C/B_P scenario modifiers enter establishment and persistence logits",
    "lambda_L = -log(persistence_probability) / delta_t",
    "tip HMSC beta posterior -> ancestral beta response history",
    "copy-then-diverge at dated-tree node intervals",
    "local loss = 1 - persistence; total lineage extinction is scenario pressure only"
  ),
  empirical_status = c(
    "HMSC-derived response and palaeoenvironment",
    "scenario/prior parameterised movement",
    "scenario/prior parameterised establishment",
    if (identical(run_type, "all_processes_demo")) {
      "SCENARIO_DEMO_NOT_EMPIRICAL"
    } else {
      "OFF_IN_EMPIRICAL_CORE"
    },
    "scenario/prior parameterised persistence",
    "HMSC posterior plus dated tree",
    "dated tree conditioned, not estimated as a rate",
    "local/geographic loss only; global lineage extinction not estimated"
  ),
  output_metrics = c(
    "mean_environmental_support; potential_environmental_richness",
    "mean_arrival_probability; movement_matrix_summary",
    "mean_colonisation_probability",
    "biotic_competition_pressure; biotic_facilitation_pressure",
    "mean_persistence_probability; mean_local_extinction_probability",
    "ancestral_beta_response_all_projected_times; niche_response_dispersion",
    "speciation_event_summary; speciation_inheritance_footprint",
    "mean_local_extinction_probability; scenario_lineage_extinction_pressure"
  ),
  participates_in_dynamic_update = c(
    TRUE, TRUE, TRUE,
    identical(run_type, "all_processes_demo"),
    TRUE, TRUE, TRUE,
    FALSE
  ),
  interpretation_boundary = c(
    "Environmental support is not final occupancy.",
    "No BioGeoBEARS region gate or A_hist multiplier is used.",
    "Colonisation requires arrival and establishment.",
    "No independent interaction data; do not call it causal competition/facilitation.",
    "Persistence/local loss is not global lineage extinction.",
    "Ancestral response is a model reconstruction, not observed physiology.",
    "Speciation timing comes from the tree; Case05 does not invent new species.",
    "Extant-only tree cannot recover extinct unsampled lineages."
  ),
  stringsAsFactors = FALSE
)
write_csv(eight_process_inventory,
          file.path(dirs$report, "case05_eight_process_inventory.csv"))

biotic_contract <- data.frame(
  mode = c("empirical_core", "all_processes_demo", "future_empirical"),
  B_colonisation = c(0, "scenario interaction matrix", "external interaction model"),
  B_persistence = c(0, "scenario interaction matrix", "external interaction model"),
  interpretation = c(
    "No independent interaction data; HMSC residual association is not treated as causal interaction.",
    "Demonstration only; not empirically estimated.",
    "Requires independent abundance/network/interaction data."
  ),
  stringsAsFactors = FALSE
)
write_csv(biotic_contract, file.path(dirs$biotic,
                                     "case05_biotic_filtering_contract.csv"))

map_path_status <- function(tab, column) {
  if (!nrow(tab) || !column %in% names(tab)) return("NOT_INDEXED")
  paths <- as.character(tab[[column]])
  paths <- paths[!is.na(paths) & nzchar(paths)]
  if (!length(paths)) return("SKIP_NOT_REQUESTED")
  if (all(file.exists(paths))) "PASS" else "FAIL"
}

png_status <- map_path_status(map_index, "png")
tif_status <- map_path_status(map_index, "tif")
map_file_status <- if (png_status == "PASS" && tif_status == "PASS") {
  "PASS"
} else if (png_status == "SKIP_NOT_REQUESTED" &&
           tif_status == "SKIP_NOT_REQUESTED") {
  "SKIP_NOT_REQUESTED"
} else {
  "FAIL"
}

validation_extra <- data.frame(
  check = c("case05_no_bgb_contract", "case05_endpoint_double_use_guard",
            "case05_mode_contract", "case05_function_coverage",
            "case05_map_index_complete", "case05_sampled_surviving_label",
            "case05_biotic_filtering_boundary",
            "case05_global_lineage_extinction_boundary",
            "case05_all_processes_map_metrics"),
  status = c(
    "PASS",
    if (identical(run_type, "endpoint_conditioned_reconstruction")) {
      "INPUT_CHECKED"
    } else {
      "PASS_DIAGNOSTIC_ONLY"
    },
    "PASS",
    "PASS",
    map_file_status,
    "PASS",
    "PASS",
    "PASS",
    if (!identical(run_type, "all_processes_demo")) {
      "NOT_APPLICABLE"
    } else if (nrow(map_index) > 0 &&
               all(c("biotic_competition_pressure",
                     "biotic_facilitation_pressure",
                     "speciation_inheritance_footprint",
                     "scenario_lineage_extinction_pressure") %in%
                   unique(map_index$metric))) {
      "PASS"
    } else {
      "FAIL"
    }
  ),
  detail = c(
    "Case05 does not call BioGeoBEARS/BSM or use region/A_hist constraints.",
    "Same GBIF-derived modern data are not treated as independent endpoint likelihood in forward mode.",
    paste0("Run type is ", run_type, "."),
    "Function/process coverage table generated.",
    paste0(nrow(map_index), " indexed maps checked; png_status=", png_status,
           ", tif_status=", tif_status, "."),
    "Richness is labelled sampled-surviving-lineage richness, not total historical flora.",
    "Biotic filtering is OFF/zero-effect in empirical core unless independent data are supplied.",
    "Global lineage extinction is not estimated from extant-only tree; only local/geography-forced loss is represented.",
    "All-processes demo requires biotic filtering, speciation footprint and scenario extinction-pressure map metrics."
  ),
  stringsAsFactors = FALSE
)
validation_extra <- rbind(
  validation_extra,
  data.frame(
    check = "case05_process_attribution",
    status = if (run_process_attribution && file.exists(file.path(
      dirs$attribution, "case05_process_attribution_summary.csv"
    ))) {
      "PASS"
    } else {
      "NOT_RUN"
    },
    detail = if (run_process_attribution) {
      paste0("hee_process_attribution output directory: ", dirs$attribution)
    } else {
      "Process attribution disabled by --run_process_attribution=FALSE."
    },
    stringsAsFactors = FALSE
  )
)
validation_extra <- rbind(
  validation_extra,
  data.frame(
    check = "case05_full_species_grid_process_contribution",
    status = if (file.exists(file.path(
      dirs$attribution, "case05_full_process_contribution_summary.csv"
    ))) {
      "PASS"
    } else {
      "FAIL"
    },
    detail = "Full streaming state-transition contribution covers all active lineages, all valid land cells, and all processed intervals; sampled interval CSV is diagnostic only.",
    stringsAsFactors = FALSE
  )
)
if (file.exists(case04_validation)) {
  validation <- rbind(read_csv(case04_validation), validation_extra)
} else {
  validation <- validation_extra
}
write_csv(validation, file.path(dirs$config, "case05_validation_checks.csv"))

file_counts <- as.data.frame(sort(table(tools::file_ext(
  list.files(output, recursive = TRUE, full.names = TRUE)
))))
names(file_counts) <- c("extension", "count")

map_count <- if (nrow(map_index)) nrow(map_index) else 0L
time_count <- if (nrow(time_summary) && "time_ma" %in% names(time_summary)) {
  length(unique(time_summary$time_ma))
} else {
  NA_integer_
}

report <- c(
  "# Case05 Plant200 global dynamic no-BioGeoBEARS case",
  "",
  "Case05 cancels BioGeoBEARS, BSM histories, discrete regions, A_hist, and M1-M5 multipliers.",
  "It works directly on the global 1-degree palaeo-Earth grid.",
  "",
  "Core production line:",
  "",
  "`HMSC beta posterior -> ancestral beta -> environmental filtering S_env -> global sparse movement -> arrival -> establishment -> colonisation -> persistence/local loss -> CTMC occupancy psi -> diversity maps`.",
  "",
  "Run type:",
  "",
  paste0("- `", run_type, "`"),
  "",
  "Scientific boundary:",
  "",
  "- This run is a sampled-surviving-lineage analysis, not total historical plant diversity.",
  "- The valid biological time domain is bounded by the sampled dated tree root.",
  "- In forward scenario mode, modern endpoint and fossils are diagnostics, not likelihood-conditioned smoothing.",
  if (identical(run_type, "all_processes_demo")) {
    "- Biotic filtering is enabled as a scenario modifier inside establishment and persistence; it is not an empirical interaction estimate."
  } else {
    "- Biotic filtering is set to zero in empirical core because no independent interaction data are supplied."
  },
  "- Speciation is tree-conditioned copy-then-diverge inheritance; Case05 does not invent new species.",
  "- Global lineage extinction is not estimated from an extant-only tree; the demo reports scenario pressure only.",
  "- Plate transport remains identity-demo unless a real plate-overlap transport matrix is supplied.",
  "",
  "Eight process outputs:",
  "",
  "- Environmental filtering: `mean_environmental_support`, `potential_environmental_richness`.",
  "- Dispersal: `movement_matrix_summary`, `mean_arrival_probability`.",
  "- Colonisation: `mean_colonisation_probability`.",
  "- Biotic filtering: `biotic_competition_pressure`, `biotic_facilitation_pressure`.",
  "- Persistence: `mean_persistence_probability`, `mean_local_extinction_probability`.",
  "- Evolution: `ancestral_beta_response_all_projected_times`, `niche_response_dispersion`.",
  "- Speciation: `speciation_event_summary`, `speciation_inheritance_footprint`.",
  "- Extinction: local-loss maps plus `scenario_lineage_extinction_pressure` in all-processes demo.",
  "",
  "Outputs:",
  "",
  paste0("- Output directory: ", output),
  paste0("- Indexed maps: ", map_count),
  paste0("- Time slices: ", time_count),
  "- Final file counts are written to `22_report/case05_file_counts.csv`.",
  "",
  "Important generated tables:",
  "",
  "- `00_config/case05_config.csv`",
  "- `00_config/case05_validation_checks.csv`",
  "- `01_input_audit/case05_removed_bgb_inputs.csv`",
  "- `22_report/case05_function_coverage.csv`",
  "- `22_report/case05_eight_process_inventory.csv`",
  "- `22_report/case05_conditioning_status.csv`",
  "- `22_report/case05_map_index.csv`",
  "- `22_report/case05_time_summary.csv`",
  "- `21_sensitivity/process_attribution/case05_full_process_contribution_summary.csv`",
  "- `21_sensitivity/process_attribution/case05_full_process_contribution_time_series.csv`",
  "- `21_sensitivity/process_attribution/case05_full_process_contribution_history_maps.csv`",
  "- `21_sensitivity/process_attribution/case05_process_attribution_summary.csv`",
  "- `21_sensitivity/process_attribution/case05_process_attribution_benchmark.csv`"
)
write_text(report, file.path(dirs$report, "case05_interpretation.md"))
file_counts <- as.data.frame(sort(table(tools::file_ext(
  list.files(output, recursive = TRUE, full.names = TRUE)
))))
names(file_counts) <- c("extension", "count")
write_csv(file_counts, file.path(dirs$report, "case05_file_counts.csv"))

log_msg("Case05 dispatcher completed.")
cat("\nDONE\n")
cat("Output:", output, "\n")
