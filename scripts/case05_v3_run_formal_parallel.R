#!/usr/bin/env Rscript

# Run a preflight-approved Case05 v3 response ensemble as restartable shards.
#
# Important: this launcher intentionally refuses the former hard-coded
# "formal" settings.  A full ensemble is allowed only after the exact process
# settings have cleared a complete-grid preflight.  The default run mode is
# endpoint_conditioned, which additionally requires an external or spatially
# cross-fitted endpoint-evidence manifest.  Without that evidence, callers
# may explicitly request run_mode=forward_scenario; outputs are then labelled
# as unconditioned forward process scenarios, never palaeodistribution
# reconstructions.  Neither mode applies a terminal map anchor.

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
if (!length(script)) script <- "scripts/case05_v3_run_formal_parallel.R"
pkg_root <- normalizePath(cfg$pkg_root %||% file.path(dirname(normalizePath(script[[1L]], winslash = "/", mustWork = FALSE)), ".."), winslash = "/", mustWork = TRUE)
workspace <- normalizePath(file.path(pkg_root, ".."), winslash = "/", mustWork = TRUE)
draw_script <- normalizePath(file.path(pkg_root, "scripts", "case05_v3_run_reference_calibrated_draw.R"), winslash = "/", mustWork = TRUE)
beta_path <- normalizePath(cfg$beta_path %||% file.path(pkg_root, "derived_inputs", "case04_plant200_hmsc_beta_posterior_25draws_20260919", "case04_hmsc_beta_posterior_25draws.csv"), winslash = "/", mustWork = TRUE)
carrier_root <- normalizePath(cfg$carrier_root %||% file.path(pkg_root, "derived_inputs", "case05_v2_paleomap_h3_carrier_325_0Ma_1deg"), winslash = "/", mustWork = TRUE)
movement_root <- normalizePath(cfg$movement_root %||% file.path(pkg_root, "derived_inputs", "case05_v3_arias_target_landscape_cache"), winslash = "/", mustWork = TRUE)
reference_path <- normalizePath(cfg$reference_path %||% file.path(pkg_root, "derived_inputs", "case05_v3_reference_environmental_support", "ancestral_response_reference_thresholds.csv"), winslash = "/", mustWork = TRUE)
output <- safe_dir(cfg$output %||% file.path(workspace, "outputs", "HmscEcoEvo", paste0("case05_v3_plant200_arias_tardis_reference_calibrated_full25_all66_", format(Sys.time(), "%Y%m%d"))))
workers <- as_int(cfg$workers, 20L)
internal_dt_myr <- as_num(cfg$internal_dt_myr, 1)
run_mode <- match.arg(tolower(cfg$run_mode %||% "endpoint_conditioned"),
                      c("endpoint_conditioned", "forward_scenario"))
preflight_results <- cfg$preflight_results %||% ""
preflight_candidate_id <- cfg$preflight_candidate_id %||% ""
endpoint_evidence_manifest <- cfg$endpoint_evidence_manifest %||% ""
seed <- as_int(cfg$seed, 20260923L)

if (!nzchar(preflight_results) || !file.exists(preflight_results) ||
    !nzchar(preflight_candidate_id)) {
  stop(
    "A completed --preflight_results CSV and --preflight_candidate_id are required. ",
    "Run scripts/case05_v3_calibrate_process_scenarios.R first; do not use ",
    "hard-coded process defaults for a full ensemble.",
    call. = FALSE
  )
}
preflight <- utils::read.csv(preflight_results, stringsAsFactors = FALSE)
needed_preflight <- c(
  "scenario_id", "preflight_gate", "root_max_carriers", "local_persistence_half_life_myr",
  "movement_rate_multiplier", "establishment_intercept"
)
if (!all(needed_preflight %in% names(preflight))) {
  stop("preflight_results lacks required columns: ",
       paste(setdiff(needed_preflight, names(preflight)), collapse = ", "),
       call. = FALSE)
}
selected_preflight <- preflight[
  as.character(preflight$scenario_id) == preflight_candidate_id &
    !is.na(as.logical(preflight$preflight_gate)) &
    as.logical(preflight$preflight_gate), , drop = FALSE
]
if (nrow(selected_preflight) != 1L) {
  stop(
    "preflight_candidate_id must select exactly one preflight_gate=TRUE row. ",
    "A numerical preflight failure means Case05 cannot be expanded to 25 draws.",
    call. = FALSE
  )
}
if (run_mode == "endpoint_conditioned") {
  stop(
    "Case05 v3 implements an unconditioned forward mean-field scenario, ",
    "not an endpoint-conditioned solver. Do not pass an endpoint manifest to make ",
    "it formal. Run scripts/case05_v4_arias_tardis_formal_readiness.R first and ",
    "use the v4 particle/pruning workflow only after its required independent or ",
    "spatially cross-fitted endpoint evidence has passed audit. v3 may be invoked ",
    "only with --run_mode=forward_scenario and its outputs remain scenarios.",
    call. = FALSE
  )
}

# Exact settings come from the accepted complete-grid preflight row.  This
# prevents a later full run from silently changing any P2/Persistence scale.
root_max_carriers <- as_int(selected_preflight$root_max_carriers[[1L]], NA_integer_)
root_occupancy <- as_num(cfg$root_occupancy, 0.8)
root_prior_mode <- tolower(cfg$root_prior_mode %||% "environment_weighted")
root_support_power <- as_num(cfg$root_support_power, 1)
movement_rate_multiplier <- as_num(selected_preflight$movement_rate_multiplier[[1L]], NA_real_)
establishment_intercept <- as_num(selected_preflight$establishment_intercept[[1L]], NA_real_)
establishment_slope <- as_num(cfg$establishment_slope, 1)
half_life <- as_num(selected_preflight$local_persistence_half_life_myr[[1L]], NA_real_)
persistence_slope <- as_num(cfg$persistence_slope, 0.75)

movement_index <- file.path(movement_root, "case05_v3_movement_cache_index.csv")
if (!file.exists(movement_index) || !file.exists(reference_path) || workers < 1L ||
    internal_dt_myr <= 0 || half_life <= 0 || !is.finite(root_max_carriers) ||
    root_max_carriers < 1L || !is.finite(establishment_intercept) ||
    !is.finite(movement_rate_multiplier) || movement_rate_multiplier < 0) {
  stop("A complete v3 movement cache, ancestral reference table, and positive settings are required.", call. = FALSE)
}
cache <- utils::read.csv(movement_index, stringsAsFactors = FALSE)
if (nrow(cache) != 66L || any(cache$n_sources_truncated > 0L) || any(!file.exists(cache$file))) {
  stop("The Case05 v3 movement cache must contain all 66 validated, non-truncated slices.", call. = FALSE)
}
baseline_persistence <- exp(-log(2) / half_life)
persistence_intercept <- stats::qlogis(baseline_persistence)

beta <- utils::read.csv(beta_path, check.names = FALSE, stringsAsFactors = FALSE)
draws <- unique(as.character(beta$response_draw))
if (!is.null(cfg$response_draw) && nzchar(cfg$response_draw)) {
  requested <- trimws(strsplit(cfg$response_draw, ",", fixed = TRUE)[[1L]])
  draws <- draws[draws %in% requested]
}
if (!length(draws)) stop("No selected HMSC response draws were found.", call. = FALSE)
workers <- min(workers, length(draws))
dirs <- list(
  config = safe_dir(file.path(output, "00_config")),
  logs = safe_dir(file.path(output, "01_run_logs")),
  shards = safe_dir(file.path(output, "02_posterior_draws")),
  summaries = safe_dir(file.path(output, "05_summaries"))
)
write_csv(data.frame(
  field = c("case_id", "run_mode", "preflight_results", "preflight_candidate_id", "endpoint_evidence_manifest", "process_mode", "n_response_draws", "n_workers", "earth_slices",
            "internal_dt_myr", "plate_carriage", "movement_kernel", "candidate_tail_rule",
            "root_prior_mode", "root_max_carriers", "root_occupancy",
            "movement_rate_multiplier", "establishment_intercept", "establishment_slope",
            "local_persistence_half_life_myr", "baseline_persistence_probability_per_myr",
            "persistence_intercept", "extant_tree_survival_conditioning", "terminal_anchor", "endpoint_role", "scientific_boundary"),
  value = c(
    "case05_v3_plant200", run_mode, normalizePath(preflight_results, winslash = "/", mustWork = TRUE),
    preflight_candidate_id,
    if (nzchar(endpoint_evidence_manifest)) normalizePath(endpoint_evidence_manifest, winslash = "/", mustWork = TRUE) else "not_supplied",
    if (run_mode == "endpoint_conditioned") "endpoint_conditioned_process_ensemble" else "reference_calibrated_forward_process_scenario",
    length(draws), workers, nrow(cache), internal_dt_myr,
    "stable_PALEOMAP_H3_carrier_identity",
    "Arias_spherical_normal_great_circle_x_destination_landscape_weight",
    unique(cache$tail_radius_sigma)[[1L]],
    root_prior_mode, root_max_carriers, root_occupancy,
    movement_rate_multiplier, establishment_intercept, establishment_slope,
    half_life, baseline_persistence, persistence_intercept, "disabled_in_forward_scenario", "none",
    if (run_mode == "endpoint_conditioned") "external_or_cross_fitted_manifest_supplied" else "training_diagnostic_only; no independent terminal likelihood or particle reweighting",
    if (run_mode == "endpoint_conditioned") {
      "Endpoint evidence is supplied through the declared manifest. Validate its provenance and use before calling any output a palaeodistribution reconstruction."
    } else {
      "Modern Plant200 records informed both the HMSC posterior and the reference response threshold. Extant-tree survival is not reimposed at each numerical time step: zero-support trajectories must be rejected or weighted only in an explicit particle/endpoint workflow. This is an unconditioned forward mean-field process scenario, not an independently endpoint-conditioned posterior reconstruction of historical plant distributions."
    }
  ),
  stringsAsFactors = FALSE
), file.path(dirs$config, "case05_v3_formal_configuration.csv"))
write_csv(data.frame(
  process = c("Environmental filtering", "Dispersal", "Colonisation", "Biotic filtering",
              "Evolution", "Speciation", "Persistence", "Extinction"),
  implementation = c(
    "HMSC posterior Beta -> direct ancestral response -> modern-prevalence-calibrated relative environmental support; intercept excluded",
    "Arias spherical kernel: great-circle distance x destination cell area x one topographic-permeability landscape; plate carriage is stable identity, not active movement",
    "arrival hazard x conditional establishment; parameters are calibrated process-scenario settings",
    "zero in empirical core because no independent historical interaction network is available",
    "direct conditional reconstruction of HMSC environmental responses along dated branches",
    "exact dated-tree nodes with copy-then-diverge inheritance",
    "environment-dependent local CTMC survival with declared half-life",
    "local loss plus hard geographic loss; global lineage-extinction rates remain unidentifiable from an extant-only tree and are not repaired by per-step survival renormalisation"
  ),
  stringsAsFactors = FALSE
), file.path(dirs$config, "case05_v3_process_responsibility.csv"))

shard_path <- function(draw) file.path(dirs$shards, paste0("draw_", draw))
complete <- function(draw) file.exists(file.path(shard_path(draw), "16_state_space_inference", "draw_completion.csv"))
pending <- draws[!vapply(draws, complete, logical(1))]
if (!length(pending)) {
  write_csv(data.frame(response_draw = draws, status = "already_complete"),
            file.path(dirs$summaries, "case05_v3_draw_status.csv"))
  message("All Case05 v3 shards are already complete: ", output)
  quit(save = "no", status = 0L)
}

rscript <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript)) rscript <- file.path(R.home("bin"), "Rscript")
run_one <- function(draw, settings) {
  out <- file.path(settings$output, "02_posterior_draws", paste0("draw_", draw))
  log <- file.path(settings$logs, paste0("draw_", draw, ".log"))
  args <- c(
    settings$draw_script, paste0("--response_draw=", draw), paste0("--output=", out),
    paste0("--carrier_root=", settings$carrier_root), paste0("--movement_root=", settings$movement_root),
    paste0("--reference_path=", settings$reference_path), paste0("--beta_path=", settings$beta_path),
    paste0("--internal_dt_myr=", settings$internal_dt_myr),
    paste0("--root_max_carriers=", settings$root_max_carriers),
    paste0("--root_occupancy=", settings$root_occupancy),
    paste0("--root_prior_mode=", settings$root_prior_mode),
    paste0("--root_support_power=", settings$root_support_power),
    paste0("--movement_rate_multiplier=", settings$movement_rate_multiplier),
    paste0("--establishment_intercept=", settings$establishment_intercept),
    paste0("--establishment_slope=", settings$establishment_slope),
    paste0("--persistence_intercept=", settings$persistence_intercept),
    paste0("--persistence_slope=", settings$persistence_slope),
    "--condition_extant_tree_survival=false",
    "--save_all_carrier_states=true",
    paste0("--seed=", settings$seed + match(draw, settings$draws))
  )
  status <- suppressWarnings(system2(settings$rscript, args = args, stdout = log, stderr = log))
  data.frame(
    response_draw = draw, exit_status = status,
    completed = file.exists(file.path(out, "16_state_space_inference", "draw_completion.csv")),
    output = normalizePath(out, winslash = "/", mustWork = FALSE),
    log = normalizePath(log, winslash = "/", mustWork = FALSE),
    stringsAsFactors = FALSE
  )
}
settings <- list(
  output = output, logs = dirs$logs, draw_script = draw_script, carrier_root = carrier_root,
  movement_root = movement_root, reference_path = reference_path, beta_path = beta_path,
  internal_dt_myr = internal_dt_myr, root_max_carriers = root_max_carriers,
  root_occupancy = root_occupancy, root_prior_mode = root_prior_mode,
  root_support_power = root_support_power, movement_rate_multiplier = movement_rate_multiplier,
  establishment_intercept = establishment_intercept, establishment_slope = establishment_slope,
  persistence_intercept = persistence_intercept, persistence_slope = persistence_slope,
  seed = seed, draws = draws, rscript = rscript
)
message("Running ", length(pending), " Case05 v3 response draws on ", workers, " workers.")
if (workers == 1L) {
  status <- lapply(pending, run_one, settings = settings)
} else {
  cl <- parallel::makeCluster(workers, outfile = "")
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, c("run_one", "settings"), envir = environment())
  status <- parallel::parLapply(cl, pending, function(draw) run_one(draw, settings))
}
status <- do.call(rbind, status)
write_csv(status, file.path(dirs$summaries, "case05_v3_draw_status.csv"))
if (!all(status$completed & status$exit_status == 0L)) {
  stop("One or more Case05 v3 shards failed; inspect 01_run_logs.", call. = FALSE)
}
message("All Case05 v3 response draws completed: ", output)
