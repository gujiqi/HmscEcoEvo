#!/usr/bin/env Rscript

# Preflight calibration for the Case05 v3 forward process scenario.
#
# This script deliberately uses one retained HMSC response draw and the same
# Plant200 modern records that informed the current HMSC fit.  Its output is an
# internal *reasonableness diagnostic*, never an independent validation score,
# endpoint likelihood, or particle weight.  Its only job is to reject declared
# root, movement, colonisation, and local-persistence settings that produce an
# obviously collapsed or saturated present-day state before a full 25-draw run.

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
as_int <- function(x, default) { z <- as_num(x, default); if (is.finite(z)) as.integer(z) else default }
safe_dir <- function(x) { if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE); normalizePath(x, winslash = "/", mustWork = FALSE) }
write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("data.table", quietly = TRUE)) data.table::fwrite(x, path) else utils::write.csv(x, path, row.names = FALSE)
  invisible(path)
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
file_arg <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
if (!length(script)) script <- "scripts/case05_v3_calibrate_process_scenarios.R"
pkg_root <- normalizePath(cfg$pkg_root %||% file.path(dirname(normalizePath(script[[1L]], winslash = "/", mustWork = FALSE)), ".."), winslash = "/", mustWork = TRUE)
workspace <- normalizePath(file.path(pkg_root, ".."), winslash = "/", mustWork = TRUE)
rscript <- cfg$rscript %||% file.path(R.home("bin"), "Rscript.exe")
draw_script <- normalizePath(file.path(pkg_root, "scripts", "case05_v3_run_reference_calibrated_draw.R"), winslash = "/", mustWork = TRUE)
score_script <- normalizePath(file.path(pkg_root, "scripts", "case05_v3_score_training_diagnostic.R"), winslash = "/", mustWork = TRUE)
beta_path <- cfg$beta_path %||% file.path(pkg_root, "derived_inputs", "case04_plant200_hmsc_beta_posterior_25draws_20260919", "case04_hmsc_beta_posterior_25draws.csv")
draw_id <- cfg$response_draw %||% "chain01_iter0090"
output <- safe_dir(cfg$output %||% file.path(workspace, "outputs", "HmscEcoEvo", paste0("case05_v3_process_preflight_", format(Sys.time(), "%Y%m%d"))))
internal_dt_myr <- as_num(cfg$internal_dt_myr, 1)
seed <- as_int(cfg$seed, 20260923L)

if (!file.exists(rscript) || !file.exists(draw_script) || !file.exists(score_script) ||
    !file.exists(beta_path) || !is.finite(internal_dt_myr) || internal_dt_myr <= 0) {
  stop("Rscript, Case05 draw/score scripts, Beta posterior, and positive internal_dt_myr are required.", call. = FALSE)
}
beta <- utils::read.csv(beta_path, stringsAsFactors = FALSE)
if (!draw_id %in% as.character(beta$response_draw)) {
  stop("response_draw is absent from the HMSC Beta posterior: ", draw_id, call. = FALSE)
}

# These are deliberately broad prior/sensitivity settings, not estimates. All
# candidates use the same root-centre seed and the full 66-slice Earth history.
# Half-life is converted below to a one-Myr persistence probability so numerical
# substeps cannot alter the biological persistence scenario.
candidates <- data.frame(
  scenario_id = c(
    "compact_long_persistence",
    "broad_long_persistence",
    "broad_faster_spherical_dispersal",
    "broad_faster_higher_establishment"
  ),
  root_max_carriers = c(24L, 256L, 256L, 256L),
  root_occupancy = c(0.8, 0.8, 0.8, 0.8),
  local_persistence_half_life_myr = c(300, 300, 300, 300),
  persistence_slope = c(0.75, 0.75, 0.75, 0.75),
  movement_rate_multiplier = c(1, 1, 3, 3),
  establishment_intercept = c(-2, -2, -2, -0.75),
  establishment_slope = c(1, 1, 1, 1),
  root_prior_mode = "environment_weighted",
  root_support_power = 1,
  interpretation = c(
    "Narrow root prior; long declared local-persistence scenario.",
    "Broader root-prior sensitivity; same process rates.",
    "Broader root plus a threefold Arias-kernel emigration sensitivity.",
    "Broader root, faster movement, and a higher declared establishment baseline."
  ),
  stringsAsFactors = FALSE
)
write_csv(candidates, file.path(output, "00_config", "preflight_candidate_grid.csv"))
writeLines(c(
  "# Case05 v3 process preflight",
  "",
  "This is an internal diagnostic only. Modern Plant200 records informed both the HMSC posterior and the response-reference calibration.",
  "Plant200 sites are target-occurrence cells, so taxon zeros are non-detections rather than surveyed absences. The preflight writes presence-only spatial-support scores but cannot use prevalence or richness to calibrate occupancy mass.",
  "Candidates are not posterior draws, endpoint likelihoods, or empirical estimates of dispersal, colonisation, or persistence. A formal ensemble requires a surveyed or independently cross-fitted terminal-evidence manifest."
), file.path(output, "README.md"))

run_one <- function(row, i) {
  half_life <- as.numeric(row$local_persistence_half_life_myr)
  baseline_persistence <- exp(-log(2) / half_life)
  persistence_intercept <- stats::qlogis(baseline_persistence)
  scenario_out <- safe_dir(file.path(output, "01_candidate_runs", row$scenario_id))
  log_dir <- safe_dir(file.path(scenario_out, "01_run_logs"))
  draw_log <- file.path(log_dir, "draw.log")
  score_log <- file.path(log_dir, "training_score.log")
  args <- c(
    draw_script,
    paste0("--response_draw=", draw_id),
    paste0("--output=", scenario_out),
    paste0("--internal_dt_myr=", internal_dt_myr),
    paste0("--root_max_carriers=", row$root_max_carriers),
    paste0("--root_occupancy=", row$root_occupancy),
    paste0("--root_prior_mode=", row$root_prior_mode),
    paste0("--root_support_power=", row$root_support_power),
    paste0("--movement_rate_multiplier=", row$movement_rate_multiplier),
    paste0("--establishment_intercept=", row$establishment_intercept),
    paste0("--establishment_slope=", row$establishment_slope),
    paste0("--persistence_intercept=", persistence_intercept),
    paste0("--persistence_slope=", row$persistence_slope),
    "--condition_extant_tree_survival=false",
    "--save_all_carrier_states=true",
    paste0("--seed=", seed)
  )
  complete_before <- file.exists(file.path(scenario_out, "16_state_space_inference", "draw_completion.csv"))
  status_draw <- if (complete_before) {
    0L
  } else {
    suppressWarnings(system2(rscript, args = args, stdout = draw_log, stderr = draw_log))
  }
  complete <- file.exists(file.path(scenario_out, "16_state_space_inference", "draw_completion.csv"))
  status_score <- NA_integer_
  if (status_draw == 0L && complete) {
    status_score <- suppressWarnings(system2(
      rscript,
      args = c(score_script, paste0("--scenario_output=", scenario_out),
               "--max_output_distance_km=500"),
      stdout = score_log, stderr = score_log
    ))
  }
  diagnostic_dir <- file.path(scenario_out, "20_modern_training_diagnostic")
  observation_path <- file.path(diagnostic_dir, "observation_design_audit.csv")
  presence_path <- file.path(diagnostic_dir, "presence_only_training_diagnostic_summary.csv")
  binary_path <- file.path(diagnostic_dir, "training_diagnostic_summary.csv")
  observation <- if (file.exists(observation_path)) {
    utils::read.csv(observation_path, stringsAsFactors = FALSE)
  } else data.frame()
  presence_score <- if (file.exists(presence_path)) {
    utils::read.csv(presence_path, stringsAsFactors = FALSE)
  } else data.frame()
  interval_path <- file.path(scenario_out, "16_state_space_inference", "interval_audit.csv")
  audit <- if (file.exists(interval_path)) utils::read.csv(interval_path, stringsAsFactors = FALSE) else data.frame()
  score_value <- function(score, field) {
    if (!nrow(score) || !field %in% names(score)) return(NA_real_)
    suppressWarnings(as.numeric(score[[field]][[1L]]))
  }
  absence_eligible <- if (nrow(observation) &&
                           "bernoulli_absence_likelihood_eligible" %in% names(observation)) {
    isTRUE(as.logical(observation$bernoulli_absence_likelihood_eligible[[1L]]))
  } else FALSE
  # Older generated folders can retain a legacy binary score from before the
  # occurrence-only audit existed. Never read that stale table unless the
  # current observation-design audit explicitly permits confirmed absences.
  binary_score <- if (absence_eligible && file.exists(binary_path)) {
    utils::read.csv(binary_path, stringsAsFactors = FALSE)
  } else data.frame()
  observed_prevalence <- score_value(binary_score, "observed_prevalence")
  predicted_prevalence <- score_value(binary_score, "predicted_prevalence")
  observed_richness <- score_value(binary_score, "observed_mean_richness")
  predicted_richness <- score_value(binary_score, "predicted_mean_richness")
  presence_lift <- score_value(presence_score, "mean_log_support_lift_over_uniform")
  observation_design <- if (nrow(observation) && "observation_design" %in% names(observation)) {
    as.character(observation$observation_design[[1L]])
  } else NA_character_
  data.frame(
    scenario_id = row$scenario_id,
    draw_id = draw_id,
    draw_exit_status = status_draw,
    training_score_exit_status = status_score,
    draw_complete = complete,
    root_max_carriers = row$root_max_carriers,
    local_persistence_half_life_myr = half_life,
    movement_rate_multiplier = row$movement_rate_multiplier,
    establishment_intercept = row$establishment_intercept,
    baseline_persistence_probability_per_myr = baseline_persistence,
    observation_design = observation_design,
    bernoulli_absence_likelihood_eligible = absence_eligible,
    presence_only_mean_log_support_lift_over_uniform = presence_lift,
    observed_mean_prevalence = observed_prevalence,
    predicted_mean_prevalence = predicted_prevalence,
    prevalence_ratio = predicted_prevalence / observed_prevalence,
    observed_mean_site_richness = observed_richness,
    predicted_mean_site_richness = predicted_richness,
    richness_ratio = predicted_richness / observed_richness,
    min_expected_lineage_richness_over_time = if (nrow(audit)) min(audit$mean_expected_lineage_richness, na.rm = TRUE) else NA_real_,
    max_expected_lineage_richness_over_time = if (nrow(audit)) max(audit$mean_expected_lineage_richness, na.rm = TRUE) else NA_real_,
    min_survival_probability_mean_field = if (nrow(audit)) min(audit$min_lineage_survival_probability_mean_field, na.rm = TRUE) else NA_real_,
    output = scenario_out,
    stringsAsFactors = FALSE
  )
}

results <- lapply(seq_len(nrow(candidates)), function(i) run_one(candidates[i, , drop = FALSE], i))
results <- do.call(rbind, results)
# A deliberately modest numerical gate: identify values outside a broad
# one-order-of-magnitude envelope. It is not a goodness-of-fit claim.
results$numerical_stability_gate <- with(results,
  is.finite(min_expected_lineage_richness_over_time) &
    min_expected_lineage_richness_over_time > 1e-6
)
# A prevalence/richness gate is valid only for an actual surveyed
# presence--absence endpoint. Plant200 is a target-occurrence grid, so its
# presence-only score diagnoses spatial support but cannot calibrate absolute
# occupancy mass. Such rows deliberately cannot unlock a formal ensemble.
results$preflight_gate <- with(results,
  bernoulli_absence_likelihood_eligible &
  is.finite(prevalence_ratio) & is.finite(richness_ratio) &
    prevalence_ratio >= 0.1 & prevalence_ratio <= 10 &
    richness_ratio >= 0.1 & richness_ratio <= 10 &
    is.finite(min_expected_lineage_richness_over_time) &
    min_expected_lineage_richness_over_time > 1e-6
)
write_csv(results, file.path(output, "05_summaries", "case05_v3_preflight_results.csv"))
if (!any(results$preflight_gate)) {
  stop("No candidate has both a surveyed presence-absence endpoint and a passing numerical gate. Do not launch a formal Case05 ensemble. The presence-only diagnostic may still be inspected, but it cannot calibrate prevalence or richness.", call. = FALSE)
}
message("Case05 v3 process preflight complete: ", output)
