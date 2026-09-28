#!/usr/bin/env Rscript

# Run the Case05 v2 HMSC-response ensemble as restartable, independent shards.
# Only complete shards are eligible for later GIS aggregation.

`%||%` <- function(x, y) if (is.null(x) || !length(x) || identical(x, "") || is.na(x)) y else x
parse_args <- function(args) {
  out <- list()
  for (arg in args) if (grepl("^--", arg)) {
    z <- sub("^--", "", arg); at <- regexpr("=", z, fixed = TRUE)
    if (at < 0L) out[[z]] <- TRUE else out[[substr(z, 1L, at - 1L)]] <- substr(z, at + 1L, nchar(z))
  }
  out
}
as_num <- function(x, default) { z <- suppressWarnings(as.numeric(x)[1L]); if (is.finite(z)) z else default }
as_int <- function(x, default) { as.integer(as_num(x, default)) }
safe_dir <- function(x) { if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE); normalizePath(x, winslash = "/", mustWork = FALSE) }
write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("data.table", quietly = TRUE)) data.table::fwrite(x, path) else utils::write.csv(x, path, row.names = FALSE)
  invisible(path)
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
file_arg <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
if (!length(script)) script <- "scripts/case05_v2_run_formal_parallel.R"
pkg_root <- normalizePath(cfg$pkg_root %||% file.path(dirname(normalizePath(script[[1L]], winslash = "/", mustWork = FALSE)), ".."), winslash = "/", mustWork = TRUE)
workspace <- normalizePath(file.path(pkg_root, ".."), winslash = "/", mustWork = TRUE)
draw_script <- normalizePath(file.path(pkg_root, "scripts", "case05_v2_run_carrier_draw.R"), winslash = "/", mustWork = TRUE)
beta_path <- normalizePath(cfg$beta_path %||% file.path(pkg_root, "derived_inputs", "case04_plant200_hmsc_beta_posterior_25draws_20260919", "case04_hmsc_beta_posterior_25draws.csv"), winslash = "/", mustWork = TRUE)
carrier_root <- normalizePath(cfg$carrier_root %||% file.path(pkg_root, "derived_inputs", "case05_v2_paleomap_h3_carrier_325_0Ma_1deg"), winslash = "/", mustWork = TRUE)
movement_root <- normalizePath(cfg$movement_root %||% file.path(carrier_root, "movement_cache_arias_topographic_rate0p4"), winslash = "/", mustWork = TRUE)
output <- safe_dir(cfg$output %||% file.path(workspace, "outputs", "HmscEcoEvo", paste0("case05_v2_plant200_arias_tardis_balanced_full25_all66_", format(Sys.time(), "%Y%m%d"))))
workers <- as_int(cfg$workers, 20L)
internal_dt_myr <- as_num(cfg$internal_dt_myr, 1)
root_max_carriers <- as_int(cfg$root_max_carriers, 240L)
root_occupancy <- as_num(cfg$root_occupancy, 0.8)
root_prior_mode <- tolower(cfg$root_prior_mode %||% "environment_weighted")
root_support_power <- as_num(cfg$root_support_power, 1)
establishment_intercept <- as_num(cfg$establishment_intercept, -0.25)
establishment_slope <- as_num(cfg$establishment_slope, 0.85)
half_life <- as_num(cfg$local_persistence_half_life_myr, 50)
persistence_slope <- as_num(cfg$persistence_slope, 0.75)
seed <- as_int(cfg$seed, 20260922L)
if (!file.exists(file.path(movement_root, "case05_v2_movement_cache_index.csv")) || workers < 1L || internal_dt_myr <= 0 || half_life <= 0) {
  stop("A completed movement cache and positive run settings are required.", call. = FALSE)
}

# A transparent baseline probability corresponding to one interpretable local
# occupancy half-life. It is a declared scenario parameter, not an HMSC fit.
baseline_persistence <- exp(-log(2) / half_life)
persistence_intercept <- stats::qlogis(baseline_persistence)
beta <- utils::read.csv(beta_path, check.names = FALSE, stringsAsFactors = FALSE)
draws <- unique(as.character(beta$response_draw))
if (!is.null(cfg$response_draw) && nzchar(cfg$response_draw)) {
  requested <- trimws(strsplit(cfg$response_draw, ",", fixed = TRUE)[[1L]])
  draws <- draws[draws %in% requested]
}
if (!length(draws)) stop("No response_draw rows were found in beta_path.", call. = FALSE)
workers <- min(workers, length(draws))
dirs <- list(config = safe_dir(file.path(output, "00_config")),
             logs = safe_dir(file.path(output, "01_run_logs")),
             shards = safe_dir(file.path(output, "02_posterior_draws")),
             summaries = safe_dir(file.path(output, "05_summaries")))
write_csv(data.frame(
  field = c("case_id", "process_mode", "n_response_draws", "n_workers", "earth_slices", "internal_dt_myr", "plate_carriage", "movement_kernel", "source_emigration_rate_per_myr", "root_prior_mode", "root_support_power", "local_persistence_half_life_myr", "baseline_persistence_probability_per_myr", "persistence_intercept", "terminal_anchor", "scientific_boundary"),
  value = c("case05_v2_plant200", "forward_mean_field_process_scenario", length(draws), workers, 66, internal_dt_myr, "stable_PALEOMAP_H3_carrier_identity", "Arias_spherical_normal_single_topographic_resistance_cost", 0.4, root_prior_mode, root_support_power, half_life, baseline_persistence, persistence_intercept, "none", "No independent endpoint likelihood or BioGeoBEARS regional history. These are process-constrained forward scenarios, not unique posterior palaeodistribution reconstructions."),
  stringsAsFactors = FALSE
), file.path(dirs$config, "case05_v2_formal_configuration.csv"))
write_csv(data.frame(
  process = c("Environmental filtering", "Dispersal", "Colonisation", "Biotic filtering", "Evolution", "Speciation", "Persistence", "Extinction"),
  implementation = c("HMSC posterior Beta -> ancestral direct-BM response -> palaeoenvironment; intercept excluded", "Arias spherical local kernel on sparse carrier graph; plate carriage is stable identity, not active movement", "arrival hazard x environmental establishment", "neutral zero in empirical core without independent interaction data", "direct conditional reconstruction of HMSC environmental responses along dated branches", "exact dated-tree nodes with copy-then-diverge inheritance", "environment-dependent local CTMC survival with declared half-life scenario", "local loss plus hard geographic loss; global lineage extinction is not identifiable from an extant-only tree"),
  stringsAsFactors = FALSE
), file.path(dirs$config, "case05_v2_process_responsibility.csv"))

shard_path <- function(draw) file.path(dirs$shards, paste0("draw_", draw))
complete <- function(draw) file.exists(file.path(shard_path(draw), "16_state_space_inference", "draw_completion.csv"))
pending <- draws[!vapply(draws, complete, logical(1))]
if (!length(pending)) {
  write_csv(data.frame(response_draw = draws, status = "already_complete"), file.path(dirs$summaries, "case05_v2_draw_status.csv"))
  message("All Case05 v2 shards are already complete: ", output)
  quit(save = "no", status = 0L)
}

rscript <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript)) rscript <- file.path(R.home("bin"), "Rscript")
run_one <- function(draw, settings) {
  out <- file.path(settings$output, "02_posterior_draws", paste0("draw_", draw))
  log <- file.path(settings$logs, paste0("draw_", draw, ".log"))
  args <- c(settings$draw_script, paste0("--response_draw=", draw), paste0("--output=", out), paste0("--movement_root=", settings$movement_root), paste0("--beta_path=", settings$beta_path), paste0("--internal_dt_myr=", settings$internal_dt_myr), paste0("--root_max_carriers=", settings$root_max_carriers), paste0("--root_occupancy=", settings$root_occupancy), paste0("--root_prior_mode=", settings$root_prior_mode), paste0("--root_support_power=", settings$root_support_power), paste0("--establishment_intercept=", settings$establishment_intercept), paste0("--establishment_slope=", settings$establishment_slope), paste0("--persistence_intercept=", settings$persistence_intercept), paste0("--persistence_slope=", settings$persistence_slope), "--save_all_carrier_states=true", paste0("--seed=", settings$seed + match(draw, settings$draws)))
  status <- suppressWarnings(system2(settings$rscript, args = args, stdout = log, stderr = log))
  data.frame(response_draw = draw, exit_status = status, completed = file.exists(file.path(out, "16_state_space_inference", "draw_completion.csv")), output = normalizePath(out, winslash = "/", mustWork = FALSE), log = normalizePath(log, winslash = "/", mustWork = FALSE), stringsAsFactors = FALSE)
}
settings <- list(output = output, logs = dirs$logs, draw_script = draw_script, movement_root = movement_root, beta_path = beta_path, internal_dt_myr = internal_dt_myr, root_max_carriers = root_max_carriers, root_occupancy = root_occupancy, root_prior_mode = root_prior_mode, root_support_power = root_support_power, establishment_intercept = establishment_intercept, establishment_slope = establishment_slope, persistence_intercept = persistence_intercept, persistence_slope = persistence_slope, seed = seed, draws = draws, rscript = rscript)
message("Running ", length(pending), " Case05 v2 response draws on ", workers, " workers.")
if (workers == 1L) {
  status <- lapply(pending, run_one, settings = settings)
} else {
  cl <- parallel::makeCluster(workers, outfile = "")
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, c("run_one", "settings"), envir = environment())
  status <- parallel::parLapply(cl, pending, function(draw) run_one(draw, settings))
}
status <- do.call(rbind, status)
write_csv(status, file.path(dirs$summaries, "case05_v2_draw_status.csv"))
if (!all(status$completed & status$exit_status == 0L)) stop("One or more Case05 v2 shards failed; inspect 01_run_logs.", call. = FALSE)
message("All Case05 v2 formal response draws completed: ", output)
