#!/usr/bin/env Rscript

# Screen explicit Case05 rate scenarios against the 0 Ma fixed-effect HMSC
# nowcast. This is a reduced, complete-grid rate-scale screen, not a formal
# PALEOMAP analysis or independent validation.

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

as_int <- function(x, default = NA_integer_) {
  y <- suppressWarnings(as.integer(x %||% default))
  if (!is.finite(y)) default else y
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
  } else {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
}
safe_id <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
scalar_arg <- function(name, value) paste0("--", name, "=", value)
bind_rows_fill <- function(x) {
  x <- x[!vapply(x, is.null, logical(1))]
  if (!length(x)) return(data.frame())
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(as.data.frame(data.table::rbindlist(x, fill = TRUE),
                         stringsAsFactors = FALSE))
  }
  all_names <- unique(unlist(lapply(x, names), use.names = FALSE))
  x <- lapply(x, function(z) {
    missing <- setdiff(all_names, names(z))
    for (nm in missing) z[[nm]] <- NA
    z[, all_names, drop = FALSE]
  })
  do.call(rbind, x)
}

cfg <- parse_args(args)
pkg_root <- normalizePath(cfg$pkg_root %||% "C:/Users/Google/Documents/HMSC-HIST/HmscEcoEvo",
                          winslash = "/", mustWork = TRUE)
rscript <- cfg$rscript %||% file.path(R.home("bin"), "Rscript.exe")
engine <- normalizePath(cfg$engine_script %||% file.path(
  pkg_root, "scripts", "case04_plant200_global_dynamic_no_bgb_final.R"
), winslash = "/", mustWork = TRUE)
merger <- normalizePath(cfg$merger_script %||% file.path(
  pkg_root, "scripts", "case04_merge_posterior_draw_shards.R"
), winslash = "/", mustWork = TRUE)
scorer <- normalizePath(cfg$endpoint_scorer %||% file.path(
  pkg_root, "scripts", "case05_score_modern_endpoint.R"
), winslash = "/", mustWork = TRUE)
response_history <- normalizePath(cfg$response_history %||% file.path(
  pkg_root, "derived_inputs", "case04_plant200_hmsc_beta_posterior_25draws_20260919",
  "case04_hmsc_beta_posterior_25draws.csv"
), winslash = "/", mustWork = TRUE)
slice_root <- file.path(pkg_root, "derived_inputs", "case04_plant200_palaeo_1deg_slices_20260919")
palaeo_slice_index <- normalizePath(cfg$palaeo_slice_index %||% file.path(
  slice_root, "case04_palaeo_earth_slice_index.csv"
), winslash = "/", mustWork = TRUE)
topographic_slice_index <- normalizePath(cfg$topographic_slice_index %||% file.path(
  slice_root, "case04_topography_slice_index.csv"
), winslash = "/", mustWork = TRUE)
output <- safe_dir(cfg$output %||% file.path(
  pkg_root, "outputs", "case05_endpoint_rate_screen_20260922"
))
logs <- safe_dir(file.path(output, "01_logs"))
screen_species <- as_int(cfg$screen_species, 20L)
screen_times <- as_int(cfg$screen_times, 5L)
response_draw <- cfg$response_draw %||% "chain01_iter0090"
seed <- as_int(cfg$seed, 20260922L)
if (screen_species < 2L || screen_times < 2L) {
  stop("screen_species and screen_times must both be >= 2.", call. = FALSE)
}

# These are broad rate-scale scenarios, not values inferred from the modern
# HMSC. The retained candidate must pass a full-times PALEOMAP pilot next.
default_candidates <- data.frame(
  candidate_id = c("legacy_default", "movement_reduced", "balanced_reduced",
                   "loss_increased", "conservative"),
  movement_intercept = c(-0.25, -1.50, -1.25, -0.75, -1.50),
  establishment_intercept = c(-0.25, -0.25, -0.75, -0.50, -1.00),
  persistence_intercept = c(1.25, 1.25, 0.75, 0.00, 0.25),
  root_rho = c(0.55, 0.55, 0.35, 0.35, 0.25),
  local_dispersal_scale_km = 200,
  establishment_slope = 0.85,
  persistence_slope = 0.75,
  stringsAsFactors = FALSE
)
candidate_csv <- cfg$candidate_csv %||% ""
candidates <- if (nzchar(candidate_csv)) {
  candidate_csv <- normalizePath(candidate_csv, winslash = "/", mustWork = TRUE)
  supplied <- read_csv(candidate_csv)
  required <- c("candidate_id", "movement_intercept", "establishment_intercept",
                "persistence_intercept", "root_rho")
  missing <- setdiff(required, names(supplied))
  if (length(missing)) {
    stop("candidate_csv is missing required columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  optional <- c("local_dispersal_scale_km", "establishment_slope",
                "persistence_slope")
  for (nm in optional) {
    if (!nm %in% names(supplied)) supplied[[nm]] <- default_candidates[[nm]][[1L]]
  }
  supplied <- supplied[, names(default_candidates), drop = FALSE]
  if (anyDuplicated(supplied$candidate_id) || any(!nzchar(as.character(supplied$candidate_id)))) {
    stop("candidate_csv candidate_id values must be non-empty and unique.", call. = FALSE)
  }
  numeric_cols <- setdiff(names(default_candidates), "candidate_id")
  for (nm in numeric_cols) supplied[[nm]] <- suppressWarnings(as.numeric(supplied[[nm]]))
  if (any(!is.finite(as.matrix(supplied[numeric_cols])))) {
    stop("candidate_csv contains non-finite rate or scale values.", call. = FALSE)
  }
  supplied
} else default_candidates
write_csv(candidates, file.path(output, "00_config", "endpoint_rate_screen_candidates.csv"))
writeLines(c(
  "Reduced complete-grid endpoint rate-scale screen.",
  "Uses identity transport and a reduced time/taxon subset only to reject implausible rate scales.",
  "It is not independent validation or a formal historical reconstruction.",
  if (nzchar(candidate_csv)) paste("Candidate table:", candidate_csv) else
    "Candidate table: built-in broad screen."
), file.path(output, "00_config", "scientific_boundary.txt"))

results <- vector("list", nrow(candidates))
for (i in seq_len(nrow(candidates))) {
  candidate <- candidates[i, , drop = FALSE]
  cid <- safe_id(candidate$candidate_id)
  shard <- file.path(output, "02_candidate_shards", cid, "posterior_shards",
                     paste0("draw_", safe_id(response_draw)))
  scenario <- file.path(output, "03_candidate_results", cid)
  log <- file.path(logs, paste0(cid, ".log"))
  args_engine <- c(
    engine, scalar_arg("pkg_root", pkg_root), "--quick=false",
    scalar_arg("response_history", response_history),
    scalar_arg("response_draw_ids", response_draw),
    scalar_arg("palaeo_slice_index", palaeo_slice_index),
    scalar_arg("topographic_slice_index", topographic_slice_index),
    scalar_arg("max_species", screen_species), "--time_selection=evenly_spaced",
    scalar_arg("n_time_slices", screen_times), "--transport_mode=identity_demo",
    "--dispersal_scheme=topographic_limited", "--ldd_per_source=0",
    "--run_biotic_filtering=false", "--run_speciation_demo=false",
    "--run_lineage_extinction_demo=false", "--make_png=false", "--make_tif=false",
    "--write_attribution_state=false", "--write_modern_tip_occupancy=true",
    "--metric_storage=rds_per_time", "--checkpoint_every_draw=true",
    "--resume_checkpoint=true", "--cache_movement_matrices=true",
    "--cache_movement_in_memory=false", "--internal_dt=5",
    "--max_hazard_per_step=0.10", "--persistence_reference_myr=5",
    "--reporting_interval_myr=1", "--root_initialisation=compact_shared_root_patch",
    "--root_seed_count=1", "--root_max_cells=600", "--root_support_power=2",
    scalar_arg("movement_intercept", candidate$movement_intercept),
    scalar_arg("local_dispersal_scale_km", candidate$local_dispersal_scale_km),
    scalar_arg("establishment_intercept", candidate$establishment_intercept),
    scalar_arg("establishment_slope", candidate$establishment_slope),
    scalar_arg("persistence_intercept", candidate$persistence_intercept),
    scalar_arg("persistence_slope", candidate$persistence_slope),
    scalar_arg("root_rho", candidate$root_rho), scalar_arg("seed", seed),
    scalar_arg("particle_seed", seed), scalar_arg("output", shard)
  )
  engine_status <- system2(rscript, args = shQuote(args_engine), stdout = log, stderr = log)
  merge_status <- score_status <- NA_integer_
  comparison_path <- file.path(scenario, "20_modern_endpoint_validation",
                               "dynamic_vs_hmsc_nowcast_endpoint_comparison.csv")
  gate_path <- file.path(scenario, "20_modern_endpoint_validation",
                         "modern_endpoint_calibration_gate.csv")
  if (identical(as.integer(engine_status), 0L)) {
    merge_status <- system2(rscript, args = shQuote(c(
      merger, scalar_arg("shard_root", dirname(shard)), scalar_arg("output", scenario)
    )), stdout = log, stderr = log)
  }
  if (identical(as.integer(merge_status), 0L)) {
    score_status <- system2(rscript, args = shQuote(c(
      scorer, scalar_arg("pkg_root", pkg_root), scalar_arg("scenario_output", scenario),
      "--data_role=training_diagnostic"
    )), stdout = log, stderr = log)
  }
  comparison <- if (file.exists(comparison_path)) read_csv(comparison_path) else data.frame()
  gate <- if (file.exists(gate_path)) read_csv(gate_path) else data.frame()
  result <- cbind(candidate, engine_status = as.integer(engine_status),
                  merge_status = as.integer(merge_status),
                  score_status = as.integer(score_status), scenario_output = scenario,
                  log = log, stringsAsFactors = FALSE)
  if (nrow(comparison)) result <- cbind(result, comparison[1L, , drop = FALSE])
  # Training cells were also used for the modular HMSC fit.  Their direct
  # agreement is retained as a diagnostic, but it cannot be treated as an
  # independent endpoint gate.  Required gates are dynamic-vs-nowcast checks;
  # a spatial holdout or external endpoint upgrades every row to PASS/FAIL.
  result$all_required_calibration_checks_pass <- nrow(gate) > 0L &&
    all(gate$status %in% c("PASS", "DIAGNOSTIC_ONLY"))
  results[[i]] <- result
  write_csv(bind_rows_fill(results[seq_len(i)]),
            file.path(output, "04_summary", "endpoint_rate_screen_summary_live.csv"))
}
summary <- bind_rows_fill(results)
write_csv(summary, file.path(output, "04_summary", "endpoint_rate_screen_summary.csv"))
ranking <- summary[order(abs(log(pmax(summary$dynamic_to_nowcast_prevalence_ratio, 1e-12))),
                         -summary$cell_richness_pearson, na.last = TRUE), , drop = FALSE]
ranking$screen_rank <- seq_len(nrow(ranking))
write_csv(ranking, file.path(output, "04_summary", "endpoint_rate_screen_ranking.csv"))
cat("Case05 endpoint rate-scale screen completed: ", output, "\n", sep = "")
