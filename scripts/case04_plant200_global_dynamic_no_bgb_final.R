#!/usr/bin/env Rscript

# Legacy scenario-calibrated global occupancy engine. Current Case04 is the
# ancestral-environmental-support study; current Case05 uses location pruning.
# This file remains runnable to reproduce earlier sensitivity analyses only.

# Case04 final global-grid no-BioGeoBEARS dynamic run.
#
# This script is designed for full/standard Plant200 output. It uses a
# streaming matrix engine rather than constructing a huge lineage-cell-time
# table in memory. BioGeoBEARS, BSM, discrete regions, A_hist and M1-M5 are not
# used.

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x)) y else x

parse_args <- function(args) {
  out <- list()
  for (a in args) {
    if (!grepl("^--", a)) next
    kv <- sub("^--", "", a)
    pos <- regexpr("=", kv, fixed = TRUE)
    if (pos < 0) out[[kv]] <- TRUE else {
      out[[substr(kv, 1, pos - 1)]] <- substr(kv, pos + 1, nchar(kv))
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

append_csv <- function(x, path, append = file.exists(path)) {
  safe_dir_create(dirname(path))
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(x, path, append = append, col.names = !append)
  } else {
    utils::write.table(
      x, path, sep = ",", row.names = FALSE, col.names = !append,
      append = append, qmethod = "double", fileEncoding = "UTF-8"
    )
  }
  invisible(path)
}

write_text <- function(x, path) {
  safe_dir_create(dirname(path))
  writeLines(x, path, useBytes = TRUE)
  invisible(path)
}

copy_if_exists <- function(from, to_dir) {
  if (!file.exists(from) || dir.exists(from)) return(FALSE)
  safe_dir_create(to_dir)
  invisible(file.copy(from, file.path(to_dir, basename(from)), overwrite = TRUE))
}

slug_time <- function(x) paste0(format(as.numeric(x), trim = TRUE, scientific = FALSE), "Ma")
stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

cfg <- parse_args(args)
pkg_root <- normalizePath(cfg$pkg_root %||% "C:/Users/Google/Documents/HMSC-HIST/HmscEcoEvo",
                          winslash = "/", mustWork = FALSE)
previous_case04 <- normalizePath(
  cfg$previous_case04 %||%
    "C:/Users/Google/Documents/HMSC-HIST/HmscEcoEvo_API8_20260818/outputs/case04_standard200_true_bgb_all109_stage1_20260821",
  winslash = "/", mustWork = FALSE
)
plant_input_dir <- normalizePath(
  cfg$plant_input_dir %||%
    "C:/Users/Google/Documents/HMSC-HIST/HmscEcoEvo_API8_20260818/data_external/plant_200_multifamily_extant_20260820/prepared_inputs/case03_data_raw/Plant200_multifamily_extant_1deg_traits",
  winslash = "/", mustWork = FALSE
)
default_slice_root <- file.path(pkg_root, "derived_inputs",
                                "case04_plant200_palaeo_1deg_slices_20260919")
palaeo_slice_index_arg <- cfg$palaeo_slice_index %||%
  file.path(default_slice_root, "case04_palaeo_earth_slice_index.csv")
topographic_slice_index_arg <- cfg$topographic_slice_index %||%
  file.path(default_slice_root, "case04_topography_slice_index.csv")
if (!file.exists(palaeo_slice_index_arg)) palaeo_slice_index_arg <- ""
if (!file.exists(topographic_slice_index_arg)) topographic_slice_index_arg <- ""

quick <- as_bool(cfg$quick, FALSE)
max_times <- as_int(cfg$max_times, if (quick) 5L else NA_integer_)
# `max_times` remains the deliberately contiguous quick-run switch. A
# coarsened scientific scenario must opt in explicitly, because evenly spaced
# palaeo slices create longer process intervals and are not interchangeable
# with the full temporal-resolution analysis.
time_selection <- tolower(cfg$time_selection %||%
                            if (is.finite(max_times)) "contiguous_oldest" else "all")
n_time_slices <- as_int(cfg$n_time_slices,
                        if (time_selection == "all") NA_integer_ else max_times)
map_n_time_slices <- as_int(cfg$map_n_time_slices, NA_integer_)
if (!time_selection %in% c("all", "contiguous_oldest", "evenly_spaced")) {
  stop("time_selection must be one of: all, contiguous_oldest, evenly_spaced.", call. = FALSE)
}
if (time_selection != "all" && (!is.finite(n_time_slices) || n_time_slices < 2L)) {
  stop("n_time_slices must be at least 2 when time_selection is not 'all'.", call. = FALSE)
}
max_species <- as_int(cfg$max_species, if (quick) 20L else NA_integer_)
max_draws <- as_int(cfg$max_draws, if (quick) 1L else NA_integer_)
response_draw_ids_arg <- cfg$response_draw_ids %||% ""
internal_dt <- as_num(cfg$internal_dt, if (quick) 10 else 5)
max_hazard_per_step <- as_num(cfg$max_hazard_per_step, 0.10)
persistence_reference_myr <- as_num(cfg$persistence_reference_myr, 5)
# Arrival and colonisation are hazards per Myr. Their mapped probabilities must
# therefore name a biological reporting window rather than silently using a
# numerical CTMC substep. One Myr is the default common reference across P2
# diagnostics; it does not alter the state update.
reporting_interval_myr <- as_num(cfg$reporting_interval_myr, 1)
transport_mode <- tolower(cfg$transport_mode %||% "identity_demo")
valid_transport_modes <- c("identity_demo", "paleomap_h3_targetcentred")
if (!transport_mode %in% valid_transport_modes) {
  stop("transport_mode must be one of: ",
       paste(valid_transport_modes, collapse = ", "), ".", call. = FALSE)
}
plate_transport_index_arg <- cfg$plate_transport_index %||% ""
transport_min_target_coverage <- as_num(cfg$transport_min_target_coverage, 0.90)
# The CTMC step is exact for its supplied rates.  Numerical subdivision is
# needed because arrival changes with occupancy, but a global maximum over
# cells with virtually zero occupancy makes a handful of no-support cells
# dictate hundreds of unnecessary substeps.  Bound refinement by the upper
# occupied-cell loss-rate quantile; retain every cell's original loss rate in
# the actual CTMC update and record the bound for audit.
hazard_rate_quantile <- as_num(cfg$hazard_rate_quantile, 0.999)
hazard_active_occupancy_threshold <- as_num(
  cfg$hazard_active_occupancy_threshold, 1e-6
)
# This is a diagnostic sample only.  The dynamic model, all summaries and all
# GIS maps still use every valid 1 degree palaeo-land cell at each time slice.
diagnostic_points_per_time <- as_int(cfg$diagnostic_points_per_time, 40L)
checkpoint_every_draw <- as_bool(cfg$checkpoint_every_draw, TRUE)
resume_checkpoint <- as_bool(cfg$resume_checkpoint, TRUE)
cache_movement_matrices <- as_bool(cfg$cache_movement_matrices, TRUE)
# Disk-backed kernels retain the exact full-grid movement graph while avoiding
# one sparse matrix per time slice being held in RAM through every draw.
cache_movement_in_memory <- as_bool(cfg$cache_movement_in_memory, quick)
movement_kernel_cache_dir_arg <- cfg$movement_kernel_cache_dir %||% ""
root_rho <- as_num(cfg$root_rho, 0.55)
root_max_cells <- as_int(cfg$root_max_cells, if (quick) 150L else 600L)
root_initialisation <- tolower(cfg$root_initialisation %||%
                                 "compact_shared_root_patch")
valid_root_initialisation <- c("legacy_top_cells", "compact_environmental_patch",
                               "compact_shared_root_patch")
if (!root_initialisation %in% valid_root_initialisation) {
  stop("root_initialisation must be one of: ",
       paste(valid_root_initialisation, collapse = ", "), ".", call. = FALSE)
}
root_seed_count <- as_int(cfg$root_seed_count, 1L)
root_support_power <- as_num(cfg$root_support_power, 2)
movement_intercept <- as_num(cfg$movement_intercept, -0.25)
local_dispersal_scale_km <- as_num(cfg$local_dispersal_scale_km, 200)
ldd_per_source <- as_int(cfg$ldd_per_source, 0L)
ldd_intercept <- as_num(cfg$ldd_intercept, -9.5)
ldd_distance_decay <- as_num(cfg$ldd_distance_decay, 0.0012)
establishment_intercept <- as_num(cfg$establishment_intercept, -0.25)
establishment_slope <- as_num(cfg$establishment_slope, 0.85)
persistence_intercept <- as_num(cfg$persistence_intercept, 1.25)
persistence_slope <- as_num(cfg$persistence_slope, 0.75)
# Optional, explicit lineage-time process scenarios.  This is intentionally a
# separate input from ancestral environmental responses: demographic process
# baselines are not HMSC niche coefficients and must never be inferred by
# copying a modern HMSC intercept into the past.
process_parameter_history_arg <- cfg$process_parameter_history %||% ""
binary_threshold <- as_num(cfg$binary_threshold, 0.25)
occupancy_mass_tolerance <- as_num(cfg$occupancy_mass_tolerance, 1e-10)
write_modern_tip_occupancy <- as_bool(cfg$write_modern_tip_occupancy, TRUE)
# A terminal nowcast anchor is a transparent model-ensemble output for 0 Ma,
# not a substitute for particle smoothing or an independent historical
# validation.  It is off by default so unanchored forward scenarios remain
# available for process diagnosis.
terminal_anchor_mode <- tolower(cfg$terminal_anchor_mode %||% "none")
valid_terminal_anchor_modes <- c(
  "none", "hmsc_fixed_effect_convex",
  "hmsc_fixed_effect_prevalence_calibrated"
)
if (!terminal_anchor_mode %in% valid_terminal_anchor_modes) {
  stop("terminal_anchor_mode must be one of: ",
       paste(valid_terminal_anchor_modes, collapse = ", "), ".", call. = FALSE)
}
terminal_anchor_weight <- as_num(cfg$terminal_anchor_weight, 0.8)
if (!is.finite(terminal_anchor_weight) || terminal_anchor_weight < 0 ||
    terminal_anchor_weight > 1) {
  stop("terminal_anchor_weight must be finite and in [0, 1].", call. = FALSE)
}
if (!identical(terminal_anchor_mode, "none") && !write_modern_tip_occupancy) {
  stop(
    "terminal_anchor_mode requires --write_modern_tip_occupancy=true so the ",
    "unanchored and anchored 0 Ma endpoint states remain auditable.",
    call. = FALSE
  )
}
terminal_anchor_prevalence_tolerance <- suppressWarnings(as.numeric(
  strsplit(cfg$terminal_anchor_prevalence_tolerance %||% "0.90,1.10",
           ",", fixed = TRUE)[[1L]]
))
if (length(terminal_anchor_prevalence_tolerance) != 2L ||
    any(!is.finite(terminal_anchor_prevalence_tolerance)) ||
    min(terminal_anchor_prevalence_tolerance) <= 0 ||
    min(terminal_anchor_prevalence_tolerance) > 1 ||
    max(terminal_anchor_prevalence_tolerance) < 1) {
  stop(
    "terminal_anchor_prevalence_tolerance must be two comma-separated positive ",
    "numbers that contain 1, for example 0.90,1.10.",
    call. = FALSE
  )
}
terminal_anchor_prevalence_tolerance <- sort(terminal_anchor_prevalence_tolerance)
seed <- as_int(cfg$seed, 20260822L)
dispersal_scheme <- cfg$dispersal_scheme %||% "topographic_limited"
valid_dispersal_schemes <- c("forward_geographic", "topographic_limited",
                             "particle_topographic")
if (!dispersal_scheme %in% valid_dispersal_schemes) {
  stop("dispersal_scheme must be one of: ",
       paste(valid_dispersal_schemes, collapse = ", "), ".", call. = FALSE)
}
particles_per_source <- as_int(cfg$particles_per_source,
                               if (quick) 16L else 64L)
if (!is.finite(particles_per_source) || particles_per_source < 1L) {
  stop("particles_per_source must be a positive integer.", call. = FALSE)
}
particle_seed <- as_int(cfg$particle_seed, seed)
if (!is.finite(internal_dt) || internal_dt <= 0 ||
    !is.finite(max_hazard_per_step) || max_hazard_per_step <= 0 ||
    !is.finite(persistence_reference_myr) || persistence_reference_myr <= 0 ||
    !is.finite(reporting_interval_myr) || reporting_interval_myr <= 0 ||
    !is.finite(transport_min_target_coverage) ||
      transport_min_target_coverage <= 0 || transport_min_target_coverage > 1 ||
    !is.finite(root_seed_count) || root_seed_count < 1L ||
    !is.finite(root_support_power) || root_support_power <= 0 ||
    !is.finite(hazard_rate_quantile) || hazard_rate_quantile <= 0 ||
    hazard_rate_quantile > 1 ||
    !is.finite(hazard_active_occupancy_threshold) ||
    hazard_active_occupancy_threshold < 0 ||
    !is.finite(occupancy_mass_tolerance) || occupancy_mass_tolerance < 0 ||
    !is.finite(diagnostic_points_per_time) || diagnostic_points_per_time < 0L) {
  stop("Numerical integration and diagnostic-point parameters must be valid values.",
       call. = FALSE)
}
if (ldd_per_source > 0L) {
  stop(
    "Random long-distance-dispersal destinations have been removed from Case05. ",
    "Set --ldd_per_source=0. A paper-grade LDD analysis needs an explicit, ",
    "auditable precomputed edge table or a calibrated spatial kernel; it must not ",
    "sample arbitrary global land cells.",
    call. = FALSE
  )
}
if (identical(root_initialisation, "compact_shared_root_patch") &&
    root_seed_count != 1L) {
  stop(
    "compact_shared_root_patch requires --root_seed_count=1 so the root prior ",
    "is one connected palaeoland patch. Multiple disconnected roots must be ",
    "declared as a separate root-prior sensitivity scenario.",
    call. = FALSE
  )
}
make_png <- as_bool(cfg$make_png, TRUE)
make_tif <- as_bool(cfg$make_tif, TRUE)
# Posterior shards contain a complete 1-degree cell-time grid for one real
# HMSC draw. CSV is convenient for a final result, but duplicating it across
# 75 shards exhausts disk before aggregation. RDS-per-time is an explicitly
# lossless, compressed transport format for those intermediate shards.
metric_storage <- tolower(cfg$metric_storage %||% "csv")
if (!metric_storage %in% c("csv", "rds_per_time")) {
  stop("metric_storage must be 'csv' or 'rds_per_time'.", call. = FALSE)
}
link <- cfg$link %||% "probit"
run_biotic_filtering <- as_bool(cfg$run_biotic_filtering, FALSE)
biotic_competition_strength <- as_num(cfg$biotic_competition_strength,
                                      if (run_biotic_filtering) 0.35 else 0)
biotic_facilitation_strength <- as_num(cfg$biotic_facilitation_strength,
                                       if (run_biotic_filtering) 0.12 else 0)
biotic_persistence_strength <- as_num(cfg$biotic_persistence_strength,
                                      biotic_competition_strength)
run_speciation_demo <- as_bool(cfg$run_speciation_demo, FALSE)
run_lineage_extinction_demo <- as_bool(cfg$run_lineage_extinction_demo, FALSE)
write_attribution_state <- as_bool(cfg$write_attribution_state, TRUE)
attribution_max_cells <- as_int(cfg$attribution_max_cells,
                                if (quick) 250L else 1000L)
attribution_max_lineages <- as_int(cfg$attribution_max_lineages,
                                   if (quick) 20L else 80L)
attribution_max_rows <- as_int(cfg$attribution_max_rows,
                               if (quick) 500000L else 2000000L)

output <- cfg$output %||%
  file.path(pkg_root, "outputs",
            paste0("case04_plant200_global_dynamic_no_bgb_final_", stamp()))
output <- safe_dir_create(output)

dirs <- list(
  config = safe_dir_create(file.path(output, "00_config")),
  audit = safe_dir_create(file.path(output, "01_input_audit")),
  hmsc = safe_dir_create(file.path(output, "02_hmsc_posterior")),
  tree = safe_dir_create(file.path(output, "03_lineage_time_tree")),
  evolution = safe_dir_create(file.path(output, "04_evolution")),
  palaeo = safe_dir_create(file.path(output, "05_palaeo_environment")),
  grid = safe_dir_create(file.path(output, "06_dynamic_earth_grid")),
  transport = safe_dir_create(file.path(output, "07_plate_transport")),
  root = safe_dir_create(file.path(output, "08_root_initialisation")),
  filtering = safe_dir_create(file.path(output, "09_environmental_filtering")),
  biotic = safe_dir_create(file.path(output, "10_biotic_filtering")),
  dispersal = safe_dir_create(file.path(output, "11_dispersal")),
  colonisation = safe_dir_create(file.path(output, "12_colonisation")),
  persistence = safe_dir_create(file.path(output, "13_persistence")),
  speciation = safe_dir_create(file.path(output, "14_speciation")),
  extinction = safe_dir_create(file.path(output, "15_extinction")),
  inference = safe_dir_create(file.path(output, "16_state_space_inference")),
  species = safe_dir_create(file.path(output, "17_species_clade_outputs")),
  diversity = safe_dir_create(file.path(output, "18_diversity_maps")),
  fossil = safe_dir_create(file.path(output, "19_fossil_constraint_validation")),
  modern = safe_dir_create(file.path(output, "20_modern_endpoint_validation")),
  sensitivity = safe_dir_create(file.path(output, "21_sensitivity")),
  report = safe_dir_create(file.path(output, "22_report")),
  logs = safe_dir_create(file.path(output, "logs"))
)

# Resolve this before writing the configuration audit. A parallel Case05
# executor supplies one shared directory so posterior-draw shards reuse the
# exact same time- and scheme-specific full-grid movement kernels.
movement_kernel_cache_dir <- safe_dir_create(if (nzchar(movement_kernel_cache_dir_arg)) {
  movement_kernel_cache_dir_arg
} else {
  file.path(dirs$dispersal, "movement_kernel_cache")
})

log_file <- file.path(dirs$logs, "case04_final_no_bgb.log")
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                paste(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_msg("Case04 final no-BioGeoBEARS global dynamic run started.")

if (dir.exists(pkg_root) && requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(pkg_root, quiet = TRUE)
  log_msg("Loaded package source:", pkg_root)
} else {
  suppressPackageStartupMessages(library(HmscEcoEvo))
}
for (pkg in c("ape", "ggplot2", "Matrix")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package is missing: ", pkg, call. = FALSE)
  }
}

paths <- list(
  tree = cfg$tree %||% file.path(plant_input_dir, "tree.tre"),
  response_history = cfg$response_history %||%
    file.path(previous_case04, "04_evolution", "case04_lineage_response_history.csv"),
  palaeo_slice_index = palaeo_slice_index_arg,
  topographic_slice_index = topographic_slice_index_arg,
  palaeo_earth_state = cfg$palaeo_earth_state %||%
    file.path(previous_case04, "05_palaeo_environment", "case04_palaeo_earth_state_scaled.rds"),
  recipe = cfg$recipe %||%
    file.path(previous_case04, "03_hmsc_environment", "case04_modern_0Ma_environment_recipe.rds"),
  time_domain = cfg$time_domain %||%
    file.path(previous_case04, "00_config", "case04_time_domain_540_0Ma.csv"),
  comm = cfg$comm %||% file.path(plant_input_dir, "comm.csv"),
  sites = cfg$sites %||% file.path(plant_input_dir, "sites.csv"),
  fossils = cfg$fossils %||% file.path(plant_input_dir, "fossils.csv"),
  metadata = cfg$metadata %||% file.path(plant_input_dir, "selected_species_metadata.csv"),
  traits = cfg$traits %||% file.path(plant_input_dir, "traits_hmsc_imputed_scaled.csv"),
  trait_raw = cfg$trait_raw %||% file.path(plant_input_dir, "traits.csv"),
  process_parameter_history = process_parameter_history_arg,
  topographic_resistance = cfg$topographic_resistance %||% file.path(
    "C:/Users/Google/Documents/HMSC-HIST/derived_inputs",
    "Phanerozoic_DispersalResistance_540_0Ma_5Myr_1deg_topography_v1",
    "topographic_dispersal_resistance_540_0Ma_5Myr_1deg_v1.rds"
  ),
  topographic_connectivity_dir = cfg$topographic_connectivity_dir %||% file.path(
    "C:/Users/Google/Documents/HMSC-HIST/derived_inputs",
    "Phanerozoic_DispersalResistance_540_0Ma_5Myr_1deg_topography_v1",
    "connectivity_graphs"
  ),
  plate_transport_index = plate_transport_index_arg
)
use_streaming_slices <- nzchar(paths$palaeo_slice_index) &&
  nzchar(paths$topographic_slice_index)
required <- c("tree", "response_history", "recipe", "topographic_connectivity_dir")
if (identical(terminal_anchor_mode, "hmsc_fixed_effect_prevalence_calibrated")) {
  required <- c(required, "comm", "sites")
}
if (identical(transport_mode, "paleomap_h3_targetcentred")) {
  required <- c(required, "plate_transport_index")
}
if (use_streaming_slices) {
  required <- c(required, "palaeo_slice_index", "topographic_slice_index")
} else {
  required <- c(required, "palaeo_earth_state", "topographic_resistance")
}
inv <- data.frame(input = names(paths), path = unlist(paths, use.names = FALSE),
                  required = names(paths) %in% required |
                    (names(paths) == "process_parameter_history" &
                       nzchar(paths$process_parameter_history)),
                  exists = vapply(unlist(paths, use.names = FALSE), function(x) {
                    file.exists(x) || dir.exists(x)
                  }, logical(1)),
                  stringsAsFactors = FALSE)
write_csv(inv, file.path(dirs$audit, "case04_final_input_inventory.csv"))
miss <- inv[inv$required & !inv$exists, , drop = FALSE]
if (nrow(miss)) {
  write_csv(miss, file.path(dirs$audit, "missing_required_inputs.csv"))
  stop("Missing required inputs. See missing_required_inputs.csv", call. = FALSE)
}
audit_copy_names <- c("tree", "response_history", "recipe", "time_domain",
                      "comm", "sites", "fossils", "metadata", "traits",
                      "trait_raw")
if (nzchar(paths$process_parameter_history)) {
  audit_copy_names <- c(audit_copy_names, "process_parameter_history")
}
if (use_streaming_slices) {
  audit_copy_names <- c(audit_copy_names, "palaeo_slice_index",
                        "topographic_slice_index")
} else {
  # Do not copy multi-gigabyte cubes into every run directory.  Their exact
  # source paths and hashes remain in the inventory/configuration instead.
  audit_copy_names <- c(audit_copy_names)
}
for (nm in intersect(audit_copy_names, names(paths))) {
  copy_if_exists(paths[[nm]], dirs$audit)
}

tree <- ape::read.tree(paths$tree)
root_age <- max(ape::node.depth.edgelength(tree))
tip_labels <- tree$tip.label
recipe <- readRDS(paths$recipe)
basis_cols <- setdiff(as.character(recipe$variables), "log_sampling_effort")

response_raw <- utils::read.csv(paths$response_history, check.names = FALSE,
                                stringsAsFactors = FALSE)
tip_beta <- response_raw[response_raw$lineage %in% tip_labels, , drop = FALSE]
tip_beta <- tip_beta[, unique(c("lineage", "response_draw", "intercept",
                                "hmsc_intercept_original", basis_cols)),
                     drop = FALSE]
names(tip_beta)[names(tip_beta) == "lineage"] <- "species"
if (!"response_draw" %in% names(tip_beta)) tip_beta$response_draw <- "s1"
if (is.finite(max_species)) {
  keep <- unique(tip_beta$species)[seq_len(min(max_species, length(unique(tip_beta$species))))]
  tip_beta <- tip_beta[tip_beta$species %in% keep, , drop = FALSE]
  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, keep))
  root_age <- max(ape::node.depth.edgelength(tree))
  tip_labels <- tree$tip.label
}
if (nzchar(response_draw_ids_arg)) {
  requested_draws <- trimws(strsplit(response_draw_ids_arg, ",", fixed = TRUE)[[1L]])
  requested_draws <- requested_draws[nzchar(requested_draws)]
  missing_draws <- setdiff(requested_draws, unique(tip_beta$response_draw))
  if (length(missing_draws)) {
    stop("Requested response_draw_ids are absent from response_history: ",
         paste(missing_draws, collapse = ", "), call. = FALSE)
  }
  tip_beta <- tip_beta[tip_beta$response_draw %in% requested_draws, , drop = FALSE]
} else if (is.finite(max_draws)) {
  keep <- unique(tip_beta$response_draw)[seq_len(min(max_draws, length(unique(tip_beta$response_draw))))]
  tip_beta <- tip_beta[tip_beta$response_draw %in% keep, , drop = FALSE]
}
.hee_check_unique_keys(tip_beta, c("species", "response_draw"),
                       "tip beta response draws")
tip_hmsc_intercepts <- tip_beta[, c("species", "response_draw",
                                    "hmsc_intercept_original"), drop = FALSE]

# P2 (movement and colonisation) and Persistence may be supplied as an
# auditable lineage-time scenario.  They deliberately remain distinct from
# `beta`: HMSC environmental-response coefficients describe P1 support, while
# these parameters describe demographic rates conditional on the process
# model.  A supplied table must cover every active lineage at every requested
# time; there is no silent random or trait-imputation fallback here.
process_parameter_history <- NULL
process_parameter_mode <- "shared_scalar_process_scenario"
process_parameter_columns <- c(
  "movement_multiplier", "establishment_intercept", "establishment_slope",
  "persistence_intercept", "persistence_slope"
)
if (nzchar(paths$process_parameter_history)) {
  process_parameter_history <- utils::read.csv(
    paths$process_parameter_history, check.names = FALSE,
    stringsAsFactors = FALSE
  )
  .require_cols(process_parameter_history, c("lineage", "time_ma"),
                "process_parameter_history")
  process_parameter_history$lineage <- as.character(process_parameter_history$lineage)
  process_parameter_history$time_ma <- suppressWarnings(
    as.numeric(process_parameter_history$time_ma)
  )
  if (any(!nzchar(process_parameter_history$lineage)) ||
      any(!is.finite(process_parameter_history$time_ma))) {
    stop("process_parameter_history requires non-empty lineage and finite time_ma values.",
         call. = FALSE)
  }
  if ("response_draw" %in% names(process_parameter_history)) {
    process_parameter_history$response_draw <- as.character(
      process_parameter_history$response_draw
    )
    if (any(!nzchar(process_parameter_history$response_draw))) {
      stop("process_parameter_history response_draw values must be non-empty when supplied.",
           call. = FALSE)
    }
  }
  supplied_process_columns <- intersect(process_parameter_columns,
                                        names(process_parameter_history))
  if (!length(supplied_process_columns)) {
    stop("process_parameter_history has no recognised process column. Supply at least one of: ",
         paste(process_parameter_columns, collapse = ", "), call. = FALSE)
  }
  for (nm in supplied_process_columns) {
    process_parameter_history[[nm]] <- suppressWarnings(
      as.numeric(process_parameter_history[[nm]])
    )
    if (any(!is.finite(process_parameter_history[[nm]]))) {
      stop("process_parameter_history column '", nm,
           "' must contain finite numeric values only.", call. = FALSE)
    }
  }
  if ("movement_multiplier" %in% supplied_process_columns &&
      any(process_parameter_history$movement_multiplier < 0)) {
    stop("process_parameter_history movement_multiplier must be non-negative.",
         call. = FALSE)
  }
  process_key <- c("lineage", "time_ma")
  if ("response_draw" %in% names(process_parameter_history)) {
    process_key <- c(process_key, "response_draw")
  }
  .hee_check_unique_keys(process_parameter_history, process_key,
                         "process_parameter_history")
  process_parameter_mode <- paste0(
    "explicit_lineage_time_scenario; columns=",
    paste(supplied_process_columns, collapse = ",")
  )
  write_csv(process_parameter_history,
            file.path(dirs$config, "process_parameter_history_input.csv"))
}

process_parameters_for <- function(lineages, time_ma, response_draw) {
  out <- list(
    movement_multiplier = stats::setNames(rep(1, length(lineages)), lineages),
    establishment_intercept = stats::setNames(
      rep(establishment_intercept, length(lineages)), lineages
    ),
    establishment_slope = stats::setNames(
      rep(establishment_slope, length(lineages)), lineages
    ),
    persistence_intercept = stats::setNames(
      rep(persistence_intercept, length(lineages)), lineages
    ),
    persistence_slope = stats::setNames(
      rep(persistence_slope, length(lineages)), lineages
    )
  )
  if (is.null(process_parameter_history)) return(out)
  # `time_ma` in the environmental and response tables should be identical,
  # but retain a small tolerance for CSV decimal serialisation.
  tol <- sqrt(.Machine$double.eps) * max(1, abs(time_ma))
  z <- process_parameter_history[
    abs(process_parameter_history$time_ma - time_ma) <= tol,
    , drop = FALSE
  ]
  if ("response_draw" %in% names(z)) {
    z <- z[as.character(z$response_draw) == as.character(response_draw), ,
           drop = FALSE]
  }
  idx <- match(lineages, z$lineage)
  if (anyNA(idx)) {
    missing_lineages <- lineages[is.na(idx)]
    write_csv(
      data.frame(
        lineage = missing_lineages, time_ma = time_ma,
        response_draw = response_draw, stringsAsFactors = FALSE
      ),
      file.path(dirs$audit,
                paste0("missing_process_parameters_", slug_time(time_ma),
                       "_", gsub("[^A-Za-z0-9_.-]+", "_", response_draw),
                       ".csv"))
    )
    stop("process_parameter_history does not cover every active lineage at ",
         time_ma, " Ma for response draw ", response_draw, ".", call. = FALSE)
  }
  for (nm in intersect(process_parameter_columns, names(z))) {
    value <- z[[nm]][idx]
    if (any(!is.finite(value)) ||
        (identical(nm, "movement_multiplier") && any(value < 0))) {
      stop("Invalid aligned process_parameter_history values for ", nm,
           " at ", time_ma, " Ma.", call. = FALSE)
    }
    out[[nm]] <- stats::setNames(value, lineages)
  }
  out
}

read_slice_index <- function(path, label) {
  index <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  .require_cols(index, c("time_ma", "file"), label)
  index$time_ma <- suppressWarnings(as.numeric(index$time_ma))
  if (any(!is.finite(index$time_ma)) || anyDuplicated(index$time_ma)) {
    stop(label, " must have one finite, unique time_ma per file.", call. = FALSE)
  }
  if (any(!file.exists(index$file))) {
    missing <- index[!file.exists(index$file), c("time_ma", "file"), drop = FALSE]
    write_csv(missing, file.path(dirs$audit, paste0(label, "_missing_files.csv")))
    stop(label, " contains missing time-slice files.", call. = FALSE)
  }
  index
}

# Both modes expose the same load_earth_time() interface.  In the default
# streaming mode it reads one *complete* 1 degree palaeo-land slice plus its
# matching topographic layer, validates the join, and releases it after the
# time step.  The legacy large-RDS fallback remains only for reproducibility.
if (use_streaming_slices) {
  earth_index <- read_slice_index(paths$palaeo_slice_index, "palaeo_earth_slice_index")
  topography_index <- read_slice_index(paths$topographic_slice_index,
                                       "topographic_slice_index")
  all_env_times <- sort(unique(earth_index$time_ma), decreasing = TRUE)
  missing_topography_times <- setdiff(all_env_times, topography_index$time_ma)
  if (length(missing_topography_times)) {
    write_csv(data.frame(time_ma = missing_topography_times),
              file.path(dirs$audit, "topographic_slice_times_missing_for_earth.csv"))
    stop("Every palaeo-Earth time slice must have a matching topographic slice.",
         call. = FALSE)
  }
  load_earth_time <- function(tm) {
    ei <- match(tm, earth_index$time_ma)
    ti <- match(tm, topography_index$time_ma)
    if (is.na(ei) || is.na(ti)) stop("Missing indexed slice at ", tm, " Ma.", call. = FALSE)
    e_tm <- as.data.frame(readRDS(earth_index$file[[ei]]))
    topo_tm <- as.data.frame(readRDS(topography_index$file[[ti]]))
    .require_cols(e_tm, c("cell_id", "time_ma", "lon", "lat", "H_state", basis_cols),
                  paste0("palaeo earth slice ", tm, " Ma"))
    .require_cols(topo_tm, c("cell_id", "time_ma", "topographic_resistance",
                             "movement_allowed"),
                  paste0("topographic slice ", tm, " Ma"))
    .hee_check_unique_keys(e_tm, c("cell_id", "time_ma"),
                           paste0("palaeo earth slice ", tm, " Ma"))
    .hee_check_unique_keys(topo_tm, c("cell_id", "time_ma"),
                           paste0("topographic slice ", tm, " Ma"))
    e_tm$H_state <- as.integer(.hee_hmscee_prob(e_tm$H_state) > 0)
    key_e <- paste(e_tm$cell_id, e_tm$time_ma, sep = "\r")
    key_t <- paste(topo_tm$cell_id, topo_tm$time_ma, sep = "\r")
    idx <- match(key_e, key_t)
    bad_land <- e_tm$H_state > 0 &
      (is.na(idx) | topo_tm$movement_allowed[idx] <= 0 |
         !is.finite(topo_tm$topographic_resistance[idx]))
    if (any(bad_land)) {
      write_csv(utils::head(e_tm[bad_land, c("cell_id", "time_ma")], 100L),
                file.path(dirs$audit,
                          paste0("topographic_resistance_unmatched_land_cells_", slug_time(tm), ".csv")))
      stop("Complete palaeo-land slice has unmatched topography at ", tm, " Ma.",
           call. = FALSE)
    }
    e_tm$topographic_resistance <- topo_tm$topographic_resistance[idx]
    e_tm$topographic_permeability <- topo_tm$topographic_permeability[idx]
    e_tm$topographic_movement_allowed <- topo_tm$movement_allowed[idx]
    e_tm[e_tm$H_state > 0, , drop = FALSE]
  }
  palaeo_input_mode <- "streaming_complete_1deg_time_slices"
} else {
  earth <- as.data.frame(readRDS(paths$palaeo_earth_state))
  .require_cols(earth, c("cell_id", "time_ma", "lon", "lat", "H_state", basis_cols),
                "palaeo_earth_state")
  earth$H_state <- as.integer(.hee_hmscee_prob(earth$H_state) > 0)
  earth <- earth[is.finite(earth$time_ma), , drop = FALSE]
  topography_object <- readRDS(paths$topographic_resistance)
  topography <- if (is.list(topography_object) && "cells" %in% names(topography_object)) {
    as.data.frame(topography_object$cells)
  } else as.data.frame(topography_object)
  .require_cols(topography, c("cell_id", "time_ma", "topographic_resistance",
                              "movement_allowed"), "topographic_resistance")
  .hee_check_unique_keys(topography, c("cell_id", "time_ma"), "topographic_resistance")
  earth_key <- paste(earth$cell_id, earth$time_ma, sep = "\r")
  topography_key <- paste(topography$cell_id, topography$time_ma, sep = "\r")
  topography_match <- match(earth_key, topography_key)
  land_missing_topography <- earth$H_state > 0 &
    (is.na(topography_match) | topography$movement_allowed[topography_match] <= 0 |
       !is.finite(topography$topographic_resistance[topography_match]))
  if (any(land_missing_topography)) {
    write_csv(utils::head(unique(earth[land_missing_topography, c("cell_id", "time_ma")]), 100L),
              file.path(dirs$audit, "topographic_resistance_unmatched_land_cells.csv"))
    stop("Palaeo-Earth land cells are missing from the topographic resistance layer.",
         call. = FALSE)
  }
  earth$topographic_resistance <- topography$topographic_resistance[topography_match]
  earth$topographic_permeability <- topography$topographic_permeability[topography_match]
  earth$topographic_movement_allowed <- topography$movement_allowed[topography_match]
  all_env_times <- sort(unique(earth$time_ma), decreasing = TRUE)
  load_earth_time <- function(tm) {
    earth[earth$time_ma == tm & earth$H_state > 0, , drop = FALSE]
  }
  palaeo_input_mode <- "legacy_full_cube_memory"
}
inside_times <- all_env_times[all_env_times <= root_age + sqrt(.Machine$double.eps)]
available_inside_times <- inside_times
if (time_selection == "contiguous_oldest" &&
    is.finite(n_time_slices) && length(inside_times) > n_time_slices) {
  # A dynamic CTMC must not jump from a representative old slice directly to
  # 0 Ma: that would turn hundreds of Myr into one process interval.
  inside_times <- utils::head(inside_times, n_time_slices)
} else if (time_selection == "evenly_spaced" &&
           is.finite(n_time_slices) && length(inside_times) > n_time_slices) {
  # Keep actual palaeo-Earth layers, including the root-nearest and modern
  # endpoint. This is a documented coarsened scenario, not a replacement for
  # the full temporal-resolution analysis.
  selected_index <- unique(round(seq.int(1L, length(inside_times),
                                        length.out = n_time_slices)))
  inside_times <- inside_times[selected_index]
}
if (length(inside_times) < 2L) {
  stop("At least two palaeo-Earth slices are required for dynamic occupancy.", call. = FALSE)
}
map_times <- inside_times
if (is.finite(map_n_time_slices) && map_n_time_slices >= 2L &&
    length(inside_times) > map_n_time_slices) {
  map_times <- inside_times[unique(round(seq.int(
    1L, length(inside_times), length.out = map_n_time_slices
  )))]
}
core_cell_columns <- c("cell_id", "time_ma", "lon", "lat", "H_state", "land_area_km2")
earth_by_time <- lapply(inside_times, function(tm) {
  e_tm <- load_earth_time(tm)
  .require_cols(e_tm, core_cell_columns, paste0("complete palaeo-land slice ", tm, " Ma"))
  e_tm[, core_cell_columns, drop = FALSE]
})
names(earth_by_time) <- as.character(inside_times)

time_slice_audit <- data.frame(
  selected_order_old_to_young = seq_along(inside_times),
  time_ma = as.numeric(inside_times),
  n_valid_land_cells = vapply(earth_by_time, nrow, integer(1)),
  time_selection = time_selection,
  n_available_time_slices = length(available_inside_times),
  n_selected_time_slices = length(inside_times),
  selected_for_map_output = inside_times %in% map_times,
  interval_to_next_younger_myr = c(abs(diff(inside_times)), NA_real_),
  interval_environment_assumption = if (time_selection == "evenly_spaced") {
    "selected_actual_palaeo_slice_piecewise_constant_with_adaptive_CTMC_substeps"
  } else {
    "all_available_or_contiguous_actual_palaeo_slices"
  },
  stringsAsFactors = FALSE
)
write_csv(time_slice_audit, file.path(dirs$audit, "selected_complete_1deg_time_slices.csv"))

# Plate carriage is an external Earth-state operation that runs before P2
# movement. In `identity_demo` mode the legacy cell-ID remap is retained only
# for reproducibility. The PALEOMAP mode requires a pre-built, target-centred
# transport table for every actual interval in the run.
plate_transport_index <- NULL
transport_for_interval <- function(time_from_ma, time_to_ma,
                                   source_cells, target_cells) {
  if (identical(transport_mode, "identity_demo")) {
    hit <- match(target_cells, source_cells)
    tab <- data.frame(
      source_cell_id = source_cells[hit[!is.na(hit)]],
      target_cell_id = target_cells[!is.na(hit)],
      target_weight = 1,
      stringsAsFactors = FALSE
    )
    return(hee_plate_grid_transport(tab, source_cells = source_cells,
                                    target_cells = target_cells))
  }
  idx <- which(abs(plate_transport_index$time_from_ma - time_from_ma) < 1e-8 &
                 abs(plate_transport_index$time_to_ma - time_to_ma) < 1e-8)
  if (length(idx) != 1L || !file.exists(plate_transport_index$file[[idx]])) {
    stop("Missing exact PALEOMAP grid transport for ", time_from_ma, " -> ",
         time_to_ma, " Ma.", call. = FALSE)
  }
  tr <- readRDS(plate_transport_index$file[[idx]])
  if (!inherits(tr, "hee_plate_grid_transport")) {
    stop("Plate transport RDS is not a hee_plate_grid_transport object: ",
         plate_transport_index$file[[idx]], call. = FALSE)
  }
  hee_plate_grid_transport(
    tr$transport, source_cells = source_cells, target_cells = target_cells,
    source_col = tr$source_col, target_col = tr$target_col,
    weight_col = tr$weight_col
  )
}
if (identical(transport_mode, "paleomap_h3_targetcentred")) {
  plate_transport_index <- utils::read.csv(paths$plate_transport_index,
                                            stringsAsFactors = FALSE,
                                            check.names = FALSE)
  .require_cols(plate_transport_index,
                c("time_from_ma", "time_to_ma", "file"),
                "plate_transport_index")
  plate_transport_index$time_from_ma <- suppressWarnings(
    as.numeric(plate_transport_index$time_from_ma)
  )
  plate_transport_index$time_to_ma <- suppressWarnings(
    as.numeric(plate_transport_index$time_to_ma)
  )
  if (any(!is.finite(plate_transport_index$time_from_ma)) ||
      any(!is.finite(plate_transport_index$time_to_ma))) {
    stop("plate_transport_index has non-finite time values.", call. = FALSE)
  }
  expected_interval <- data.frame(
    time_from_ma = inside_times[-length(inside_times)],
    time_to_ma = inside_times[-1L]
  )
  key_expected <- paste(expected_interval$time_from_ma,
                        expected_interval$time_to_ma, sep = "\r")
  key_available <- paste(plate_transport_index$time_from_ma,
                         plate_transport_index$time_to_ma, sep = "\r")
  missing_interval <- expected_interval[!key_expected %in% key_available, , drop = FALSE]
  if (nrow(missing_interval)) {
    write_csv(missing_interval, file.path(dirs$audit,
                                          "missing_plate_transport_intervals.csv"))
    stop("PALEOMAP transport index does not cover every projected interval.",
         call. = FALSE)
  }
  write_csv(plate_transport_index,
            file.path(dirs$transport, "plate_transport_index_used.csv"))
}

# Select a small, spatially balanced audit panel without changing the full-grid
# calculation.  Eight longitude by five latitude strata give forty possible
# global strata; one seeded cell is selected from each non-empty stratum and
# any shortfall is filled from the remaining complete palaeo-land grid.
with_local_seed <- function(seed_value, expr) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed_value))
  force(expr)
}

select_diagnostic_points <- function(e_tm, n_points, seed_value) {
  n_points <- min(as.integer(n_points), nrow(e_tm))
  if (!length(n_points) || n_points <= 0L || !nrow(e_tm)) {
    return(e_tm[FALSE, c("cell_id", "time_ma", "lon", "lat"), drop = FALSE])
  }
  lon_bin <- cut(e_tm$lon, breaks = seq(-180, 180, length.out = 9L),
                 include.lowest = TRUE, labels = FALSE)
  lat_bin <- cut(e_tm$lat, breaks = seq(-90, 90, length.out = 6L),
                 include.lowest = TRUE, labels = FALSE)
  strata <- paste(lon_bin, lat_bin, sep = "_")
  by_stratum <- split(seq_len(nrow(e_tm)), strata, drop = TRUE)
  selected <- with_local_seed(seed_value, {
    unlist(lapply(by_stratum, function(idx) sample(idx, size = 1L)),
           use.names = FALSE)
  })
  if (length(selected) > n_points) {
    selected <- with_local_seed(seed_value + 1L, sample(selected, n_points))
  }
  if (length(selected) < n_points) {
    available <- setdiff(seq_len(nrow(e_tm)), selected)
    selected <- c(selected, with_local_seed(seed_value + 2L,
                                             sample(available, n_points - length(selected))))
  }
  selected <- sort(unique(selected))
  keep_columns <- intersect(
    c("cell_id", "time_ma", "lon", "lat", "H_state", "land_area_km2",
      "topographic_resistance", "topographic_permeability",
      "topographic_movement_allowed", basis_cols),
    names(e_tm)
  )
  out <- e_tm[selected, keep_columns, drop = FALSE]
  out$diagnostic_rank <- seq_len(nrow(out))
  out$longitude_stratum <- lon_bin[selected]
  out$latitude_stratum <- lat_bin[selected]
  out$selection_method <- "seeded_one_per_nonempty_8x5_lon_lat_stratum_then_complete_grid_fill"
  out$selection_role <- "diagnostic_only_full_1deg_grid_used_for_model_and_maps"
  out
}

diagnostic_point_list <- lapply(inside_times, function(tm) {
  tm_seed <- as.integer((seed + round(abs(tm) * 1000)) %% .Machine$integer.max)
  # Load the full slice for the trace so the saved audit contains the actual
  # environmental and topographic inputs, not a reduced proxy table.
  select_diagnostic_points(load_earth_time(tm),
                           diagnostic_points_per_time, tm_seed)
})
diagnostic_points <- do.call(rbind, diagnostic_point_list)
diagnostic_coverage <- do.call(rbind, lapply(inside_times, function(tm) {
  n_land <- nrow(earth_by_time[[as.character(tm)]])
  n_selected <- sum(diagnostic_points$time_ma == tm)
  data.frame(
    time_ma = tm,
    n_complete_1deg_land_cells = n_land,
    requested_diagnostic_points = diagnostic_points_per_time,
    n_diagnostic_points = n_selected,
    expected_diagnostic_points = min(diagnostic_points_per_time, n_land),
    full_grid_used_for_dynamic_calculation = TRUE,
    stringsAsFactors = FALSE
  )
}))
if (any(diagnostic_coverage$n_diagnostic_points !=
        diagnostic_coverage$expected_diagnostic_points)) {
  stop("Diagnostic-point selection did not return the requested number of valid land cells.",
       call. = FALSE)
}
write_csv(diagnostic_points,
          file.path(dirs$audit, "case04_diagnostic_points_by_time.csv"))
write_csv(diagnostic_coverage,
          file.path(dirs$audit, "case04_diagnostic_points_coverage.csv"))
write_text(c(
  "# Case04 diagnostic-point protocol",
  "",
  "For each palaeoenvironmental time slice, this file selects up to 40 seeded, spatially balanced points from the complete valid 1 degree palaeo-land grid.",
  "The selection is for process-trace quality assurance only. It never subsets the dynamic occupancy calculation, diversity summaries, GeoTIFFs, or PNG maps.",
  "Each selected point is later joined to environmental support, arrival, colonisation, persistence, local extinction, occupancy, and richness metrics."
), file.path(dirs$audit, "case04_diagnostic_points_protocol.md"))
topographic_time_slug <- function(x) {
  gsub("\\.", "p", formatC(as.numeric(x), format = "f", digits = 3,
                              drop0trailing = TRUE))
}
topographic_graph_file <- function(x) {
  file.path(paths$topographic_connectivity_dir,
            paste0("topographic_connectivity_", topographic_time_slug(x), "Ma_1deg.rds"))
}
expected_graphs <- vapply(inside_times, topographic_graph_file, character(1))
if (any(!file.exists(expected_graphs))) {
  missing_graphs <- data.frame(time_ma = inside_times[!file.exists(expected_graphs)],
                               expected_file = expected_graphs[!file.exists(expected_graphs)],
                               stringsAsFactors = FALSE)
  write_csv(missing_graphs, file.path(dirs$audit, "missing_topographic_connectivity_graphs.csv"))
  stop("Topographic connectivity graphs are missing for projected time slices. ",
       "See missing_topographic_connectivity_graphs.csv", call. = FALSE)
}

set.seed(seed)
responses <- hee_evolution_ancestral_response_direct(
  tip_beta_draws = tip_beta,
  tree = tree,
  times = inside_times,
  basis_cols = basis_cols,
  species_col = "species",
  response_draw_col = "response_draw",
  intercept_col = "intercept",
  reconstruct_intercept = FALSE,
  branch_uncertainty = "none"
)
saveRDS(responses, file.path(dirs$evolution, "ancestral_beta_response_all_projected_times.rds"))
write_csv(responses, file.path(dirs$evolution, "ancestral_beta_response_all_projected_times.csv"))

# Give every run an exact, editable schema for process-rate scenarios.  It is
# generated from the *actual* active lineage-time response history, so users
# do not need to guess ancestral identifiers or interpolate node times by
# hand.  This template is a prior/sensitivity interface, not an estimate made
# from the HMSC environmental coefficients.
process_parameter_template <- unique(
  responses[, c("lineage", "time_ma", "response_draw"), drop = FALSE]
)
process_parameter_template$movement_multiplier <- 1
process_parameter_template$establishment_intercept <- establishment_intercept
process_parameter_template$establishment_slope <- establishment_slope
process_parameter_template$persistence_intercept <- persistence_intercept
process_parameter_template$persistence_slope <- persistence_slope
write_csv(
  process_parameter_template,
  file.path(dirs$evolution, "process_parameter_history_template.csv")
)

draw_ids <- unique(responses$response_draw)
draw_weight <- stats::setNames(rep(1 / length(draw_ids), length(draw_ids)), draw_ids)

edge_parent <- stats::setNames(tree$edge[, 1], tree$edge[, 2])
ancestor_chain <- function(node) {
  out <- as.integer(node)
  cur <- as.character(node)
  while (cur %in% names(edge_parent)) {
    cur <- as.character(edge_parent[[cur]])
    out <- c(out, as.integer(cur))
  }
  out
}
chain_cache <- lapply(unique(responses$node), ancestor_chain)
names(chain_cache) <- as.character(unique(responses$node))
map_to_previous <- function(cur_nodes, prev_nodes) {
  prev_chr <- as.character(prev_nodes)
  vapply(as.character(cur_nodes), function(nn) {
    ch <- as.character(chain_cache[[nn]])
    hit <- ch[ch %in% prev_chr]
    if (length(hit)) match(hit[[1]], prev_chr) else NA_integer_
  }, integer(1))
}

great_circle_km <- function(lon1, lat1, lon2, lat2) {
  r <- 6371.0088
  to_rad <- pi / 180
  p1 <- lat1 * to_rad
  p2 <- lat2 * to_rad
  dp <- (lat2 - lat1) * to_rad
  dl <- (lon2 - lon1) * to_rad
  a <- sin(dp / 2)^2 + cos(p1) * cos(p2) * sin(dl / 2)^2
  2 * r * atan2(sqrt(a), sqrt(pmax(0, 1 - a)))
}

build_movement_matrix <- function(cells, time_ma, delta_t) {
  c0 <- unique(cells[, c("cell_id", "lon", "lat"), drop = FALSE])
  c0$idx <- seq_len(nrow(c0))
  graph_path <- topographic_graph_file(time_ma)
  graph <- readRDS(graph_path)
  edges <- if (is.list(graph) && "edges" %in% names(graph)) graph$edges else graph
  edges <- as.data.frame(edges)
  .require_cols(edges, c("from_cell_id", "to_cell_id", "time_ma", "distance_km",
                         "effective_cost_km"),
                "topographic connectivity graph")
  edges <- edges[edges$from_cell_id %in% c0$cell_id & edges$to_cell_id %in% c0$cell_id,
                 , drop = FALSE]
  if (!nrow(edges)) {
    stop("No land-only topographic movement edges remain at ", time_ma, " Ma.", call. = FALSE)
  }
  cost_basis <- if (identical(dispersal_scheme, "forward_geographic")) {
    "great_circle_distance_km_no_topographic_penalty"
  } else "topographic_effective_cost_km"
  if (identical(dispersal_scheme, "forward_geographic")) {
    # Same land-only graph, but the baseline must not use relief or elevation-step
    # penalties. It is the counterfactual pure-distance forward-dispersal scenario.
    edges$effective_cost_km <- edges$distance_km
  }
  kernel <- hee_dispersal_edge_kernel(
    edges, dispersal_scale_km = local_dispersal_scale_km,
    normalise_by_source = FALSE
  )
  # Separate each source's total per-Myr movement rate from its directional
  # allocation.  Without this step, the sum of raw neighbour weights makes
  # movement strength depend on graph degree rather than the intended process
  # rate, and rugged grid cells receive an additional numerical penalty.
  rate_kernel <- hee_dispersal_rate_kernel(
    kernel, base_rate = exp(movement_intercept),
    source_permeability = "mean_weight"
  )
  kernel_value <- rate_kernel$source_direction_weight
  source_rate <- rate_kernel$source_total_rate
  if (identical(dispersal_scheme, "particle_topographic")) {
    particle_time_seed <- as.integer((particle_seed +
      as.integer(round(abs(time_ma) * 1000))) %% .Machine$integer.max)
    particle_kernel <- hee_dispersal_particle_kernel(
      transform(rate_kernel, K_movement = source_direction_weight),
      particles_per_source = particles_per_source,
      seed = particle_time_seed
    )
    kernel_value <- particle_kernel$particle_kernel
    source_rate <- particle_kernel$source_total_rate
  }
  from_idx <- match(rate_kernel$from_cell_id, c0$cell_id)
  to_idx <- match(rate_kernel$to_cell_id, c0$cell_id)
  vals <- source_rate * kernel_value
  out <- Matrix::sparseMatrix(i = to_idx, j = from_idx, x = vals,
                              dims = c(nrow(c0), nrow(c0)))
  attr(out, "topographic_graph_file") <- normalizePath(graph_path, winslash = "/", mustWork = TRUE)
  attr(out, "mean_effective_cost_km") <- mean(edges$effective_cost_km, na.rm = TRUE)
  attr(out, "max_effective_cost_km") <- max(edges$effective_cost_km, na.rm = TRUE)
  attr(out, "mean_topographic_resistance") <- mean(cells$topographic_resistance, na.rm = TRUE)
  attr(out, "dispersal_scheme") <- dispersal_scheme
  attr(out, "movement_cost_basis") <- cost_basis
  attr(out, "source_rate_definition") <-
    "exp(movement_intercept)_per_Myr_times_mean_raw_edge_permeability; directional_weights_sum_to_one"
  attr(out, "ldd_status") <-
    "disabled_in_case05_main_run; random_global_destination_sampling_removed"
  attr(out, "max_source_total_rate") <- max(rate_kernel$source_total_rate, na.rm = TRUE)
  attr(out, "max_incoming_rate") <- max(Matrix::rowSums(out), na.rm = TRUE)
  attr(out, "particles_per_source") <- if (identical(dispersal_scheme,
                                                      "particle_topographic")) {
    particles_per_source
  } else 0L
  attr(out, "particle_seed") <- if (identical(dispersal_scheme,
                                               "particle_topographic")) {
    particle_time_seed
  } else NA_integer_
  out
}

link_inverse <- function(eta) {
  if (identical(link, "logit")) stats::plogis(eta) else stats::pnorm(eta)
}
clip01 <- function(x) pmin(pmax(x, 0), 1)

calc_eta_s <- function(e_tm, r_tm) {
  X <- as.matrix(e_tm[, basis_cols, drop = FALSE])
  B <- as.matrix(r_tm[, basis_cols, drop = FALSE])
  storage.mode(X) <- "double"
  storage.mode(B) <- "double"
  eta <- X %*% t(B)
  if ("intercept" %in% names(r_tm)) {
    eta <- sweep(eta, 2, suppressWarnings(as.numeric(r_tm$intercept)), "+")
  }
  S <- clip01(link_inverse(eta))
  list(eta = eta, S = S)
}

ctmc_update <- function(Q, lambda_c, lambda_l, h) {
  ss <- lambda_c + lambda_l
  p01 <- ifelse(ss > 0, lambda_c / ss * (1 - exp(-ss * h)), 0)
  p11 <- ifelse(ss > 0, lambda_c / ss + lambda_l / ss * exp(-ss * h), 1)
  clip01((1 - Q) * p01 + Q * p11)
}

calc_biotic_filtering <- function(Q, S) {
  n <- nrow(Q)
  if (!run_biotic_filtering) {
    z <- rep(0, n)
    return(list(
      competition_pressure = z,
      facilitation_pressure = z,
      colonisation_modifier = z,
      persistence_modifier = z
    ))
  }
  occupancy_density <- clip01(rowMeans(Q, na.rm = TRUE))
  env_support <- clip01(rowMeans(S, na.rm = TRUE))
  competition <- clip01(occupancy_density * biotic_competition_strength)
  facilitation <- clip01((1 - occupancy_density) * env_support *
                           biotic_facilitation_strength)
  list(
    competition_pressure = competition,
    facilitation_pressure = facilitation,
    colonisation_modifier = facilitation - competition,
    persistence_modifier = facilitation - competition *
      biotic_persistence_strength / max(biotic_competition_strength, 1e-12)
  )
}

init_root_q <- function(S, cells, time_ma, draw_id) {
  Q <- matrix(0, nrow(S), ncol(S), dimnames = dimnames(S))
  if (root_initialisation %in% c("compact_environmental_patch",
                                 "compact_shared_root_patch")) {
    graph <- readRDS(topographic_graph_file(time_ma))
    edges <- if (is.list(graph) && "edges" %in% names(graph)) graph$edges else graph
    edges <- as.data.frame(edges)
    .require_cols(edges, c("from_cell_id", "to_cell_id"),
                  "root-initialisation connectivity graph")
    idx <- stats::setNames(seq_len(nrow(S)), rownames(S))
    from <- unname(idx[as.character(edges$from_cell_id)])
    to <- unname(idx[as.character(edges$to_cell_id)])
    keep_edge <- is.finite(from) & is.finite(to)
    adjacency <- split(to[keep_edge], from[keep_edge])
    choose_patch <- function(s, offset) {
      valid <- which(is.finite(s) & s >= 0)
      if (!length(valid)) return(integer())
      # Root centres are an explicit prior. They are sampled from environmental
      # support then grown only through the contemporaneous land graph; this
      # prevents a root range made of disconnected global top-S cells.
      old_seed <- if (exists(".Random.seed", envir = .GlobalEnv,
                             inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL
      on.exit({
        if (is.null(old_seed)) {
          if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv)
          }
        } else assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }, add = TRUE)
      set.seed(as.integer((seed + offset) %% .Machine$integer.max))
      centre_weight <- pmax(s[valid], 1e-12)^root_support_power
      centre <- sample(valid, size = 1L, prob = centre_weight)
      selected <- centre
      while (length(selected) < min(root_max_cells, length(valid))) {
        neighbour <- unique(unlist(adjacency[as.character(selected)],
                                   use.names = FALSE))
        neighbour <- setdiff(neighbour, selected)
        neighbour <- neighbour[neighbour %in% valid]
        if (!length(neighbour)) break
        candidate_weight <- pmax(s[neighbour], 1e-12)^root_support_power
        next_cell <- sample(neighbour, size = 1L, prob = candidate_weight)
        selected <- c(selected, next_cell)
      }
      selected
    }
    if (identical(root_initialisation, "compact_shared_root_patch")) {
      # At the oldest sampled-tree time, all immediately active child branches
      # descend from one unsplit root ancestor. They must therefore inherit one
      # common, connected spatial patch rather than independently sampled
      # origins. Their later differentiation is handled at dated tree nodes.
      root_support <- rowMeans(S, na.rm = TRUE)
      root_support[!is.finite(root_support)] <- 0
      patch <- choose_patch(root_support,
                            sum(utf8ToInt(paste(draw_id, "shared_root", sep = "_"))))
      if (length(patch)) {
        strength <- pmax(root_support[patch], 0)
        if (max(strength) > 0) {
          q_root <- clip01(root_rho * strength / max(strength))
          Q[patch, ] <- q_root
        }
      }
      attr(Q, "root_initialisation") <- paste0(
        "compact_shared_root_patch; seed_count=1; support_power=",
        root_support_power,
        "; one_connected_patch_shared_by_all_initial_root_children"
      )
    } else {
      for (j in seq_len(ncol(S))) {
        patch <- integer()
        for (kk in seq_len(root_seed_count)) {
          offset <- sum(utf8ToInt(paste(draw_id, j, kk, sep = "_")))
          patch <- unique(c(patch, choose_patch(S[, j], offset)))
        }
        if (length(patch)) {
          strength <- pmax(S[patch, j], 0)
          if (max(strength) > 0) {
            Q[patch, j] <- clip01(root_rho * strength / max(strength))
          }
        }
      }
      attr(Q, "root_initialisation") <- paste0(
        "compact_environmental_patch; seed_count=", root_seed_count,
        "; support_power=", root_support_power,
        "; one_connected_patch_per_initial_child_lineage_sensitivity_only"
      )
    }
    return(Q)
  }
  for (j in seq_len(ncol(S))) {
    ord <- order(S[, j], decreasing = TRUE)
    keep <- utils::head(ord[is.finite(S[ord, j])], min(root_max_cells, length(ord)))
    w <- clip01(S[keep, j])
    # root_rho is a cell-level root-occupancy scale, not a total probability
    # mass to divide across the root range.  The former implementation made a
    # 600-cell root have q around 0.001 per cell, creating an artificial
    # almost-source-free world before the first dynamic step.
    if (max(w) > 0) Q[keep, j] <- clip01(root_rho * w / max(w))
  }
  attr(Q, "root_initialisation") <-
    "legacy_environment_weighted_top_cells_disconnected_sensitivity_only"
  Q
}

new_accumulator <- function(e_tm) {
  n <- nrow(e_tm)
  list(
    cell = e_tm[, c("cell_id", "time_ma", "lon", "lat", "H_state", "land_area_km2"), drop = FALSE],
    w = 0,
    expected_lineage_richness = numeric(n),
    forward_raw_expected_lineage_richness = numeric(n),
    terminal_anchor_difference_expected_lineage_richness = numeric(n),
    terminal_anchor_weight = numeric(n),
    # Populated only at 0 Ma from the original modern tip intercepts.  It is
    # the calibrated comparison target for dynamic occupancy, not an
    # ancestral-niche term and not a deep-time predictor.
    hmsc_fixed_effect_nowcast_richness = numeric(n),
    potential_environmental_richness = numeric(n),
    binary_lineage_richness = numeric(n),
    mean_occupancy_probability = numeric(n),
    mean_environmental_support = numeric(n),
    mean_arrival_hazard_per_myr = numeric(n),
    mean_arrival_probability_reporting_interval = numeric(n),
    mean_colonisation_hazard_per_myr = numeric(n),
    mean_colonisation_probability_reporting_interval = numeric(n),
    # Kept solely for old Case05 comparison scripts. New maps use the two
    # interval-explicit fields above and record their window in metadata.
    mean_arrival_probability = numeric(n),
    mean_colonisation_probability = numeric(n),
    mean_persistence_probability = numeric(n),
    mean_local_extinction_probability = numeric(n),
    biotic_competition_pressure = numeric(n),
    biotic_facilitation_pressure = numeric(n),
    speciation_inheritance_footprint = numeric(n),
    scenario_lineage_extinction_pressure = numeric(n),
    weighted_endemism = numeric(n),
    relative_environmental_support_sum = numeric(n),
    shannon_diversity = numeric(n),
    simpson_diversity = numeric(n),
    occupancy_weighted_effective_lineages_shannon = numeric(n),
    occupancy_weighted_effective_lineages_simpson = numeric(n),
    occupancy_weighted_diversity_defined = numeric(n),
    niche_response_dispersion = numeric(n)
  )
}
accumulate_time <- function(acc, Q, S, arrival_hazard, lambda_c, persistence,
                            eta, w, area, reporting_interval_myr,
                            biotic_competition = NULL,
                            biotic_facilitation = NULL,
                            speciation_footprint = NULL,
                            lineage_extinction_pressure = NULL,
                            forward_raw_Q = NULL,
                            terminal_anchor_weight = 0) {
  richness <- rowSums(Q, na.rm = TRUE)
  if (is.null(forward_raw_Q)) forward_raw_Q <- Q
  forward_raw_Q <- as.matrix(forward_raw_Q)
  if (!identical(dim(forward_raw_Q), dim(Q))) {
    stop("forward_raw_Q must have the same cell-by-lineage dimensions as Q.",
         call. = FALSE)
  }
  forward_raw_richness <- rowSums(forward_raw_Q, na.rm = TRUE)
  env_rich <- rowSums(S, na.rm = TRUE)
  bin_rich <- rowSums(Q > binary_threshold, na.rm = TRUE)
  arrival_prob <- clip01(1 - exp(-arrival_hazard * reporting_interval_myr))
  colon_prob <- clip01(1 - exp(-lambda_c * reporting_interval_myr))
  if (is.null(biotic_competition)) biotic_competition <- numeric(nrow(Q))
  if (is.null(biotic_facilitation)) biotic_facilitation <- numeric(nrow(Q))
  if (is.null(speciation_footprint)) speciation_footprint <- numeric(nrow(Q))
  if (is.null(lineage_extinction_pressure)) {
    lineage_extinction_pressure <- numeric(nrow(Q))
  }
  # Q is a probability-weighted lineage assemblage, not an abundance sample.
  # An empty cell must remain empty: normalising a near-zero total by an
  # arbitrary epsilon incorrectly turns it into a perfectly even community.
  diversity <- hee_occupancy_weighted_diversity(
    Q, zero_tolerance = occupancy_mass_tolerance
  )
  shannon <- diversity$occupancy_weighted_shannon_entropy
  simpson <- diversity$occupancy_weighted_gini_simpson
  range_area <- colSums(Q * area, na.rm = TRUE)
  range_area[!is.finite(range_area) | range_area <= 0] <- Inf
  we <- rowSums(sweep(Q, 2, range_area, "/"), na.rm = TRUE)
  eta_center <- rowMeans(eta, na.rm = TRUE)
  nrd <- sqrt(rowMeans((eta - eta_center)^2, na.rm = TRUE))

  acc$w <- acc$w + w
  acc$expected_lineage_richness <- acc$expected_lineage_richness + w * richness
  acc$forward_raw_expected_lineage_richness <-
    acc$forward_raw_expected_lineage_richness + w * forward_raw_richness
  acc$terminal_anchor_difference_expected_lineage_richness <-
    acc$terminal_anchor_difference_expected_lineage_richness +
    w * (richness - forward_raw_richness)
  acc$terminal_anchor_weight <- acc$terminal_anchor_weight + w * terminal_anchor_weight
  acc$potential_environmental_richness <- acc$potential_environmental_richness + w * env_rich
  acc$relative_environmental_support_sum <- acc$relative_environmental_support_sum + w * env_rich
  acc$binary_lineage_richness <- acc$binary_lineage_richness + w * bin_rich
  acc$mean_occupancy_probability <- acc$mean_occupancy_probability + w * rowMeans(Q, na.rm = TRUE)
  acc$mean_environmental_support <- acc$mean_environmental_support + w * rowMeans(S, na.rm = TRUE)
  acc$mean_arrival_hazard_per_myr <- acc$mean_arrival_hazard_per_myr + w * rowMeans(arrival_hazard, na.rm = TRUE)
  acc$mean_arrival_probability_reporting_interval <- acc$mean_arrival_probability_reporting_interval + w * rowMeans(arrival_prob, na.rm = TRUE)
  acc$mean_colonisation_hazard_per_myr <- acc$mean_colonisation_hazard_per_myr + w * rowMeans(lambda_c, na.rm = TRUE)
  acc$mean_colonisation_probability_reporting_interval <- acc$mean_colonisation_probability_reporting_interval + w * rowMeans(colon_prob, na.rm = TRUE)
  acc$mean_arrival_probability <- acc$mean_arrival_probability + w * rowMeans(arrival_prob, na.rm = TRUE)
  acc$mean_colonisation_probability <- acc$mean_colonisation_probability + w * rowMeans(colon_prob, na.rm = TRUE)
  acc$mean_persistence_probability <- acc$mean_persistence_probability + w * rowMeans(persistence, na.rm = TRUE)
  acc$mean_local_extinction_probability <- acc$mean_local_extinction_probability + w * rowMeans(1 - persistence, na.rm = TRUE)
  acc$biotic_competition_pressure <- acc$biotic_competition_pressure + w * biotic_competition
  acc$biotic_facilitation_pressure <- acc$biotic_facilitation_pressure + w * biotic_facilitation
  acc$speciation_inheritance_footprint <- acc$speciation_inheritance_footprint + w * speciation_footprint
  acc$scenario_lineage_extinction_pressure <- acc$scenario_lineage_extinction_pressure + w * lineage_extinction_pressure
  acc$weighted_endemism <- acc$weighted_endemism + w * we
  acc$shannon_diversity <- acc$shannon_diversity + w * shannon
  acc$simpson_diversity <- acc$simpson_diversity + w * simpson
  acc$occupancy_weighted_effective_lineages_shannon <-
    acc$occupancy_weighted_effective_lineages_shannon +
    w * diversity$occupancy_weighted_effective_lineages_shannon
  acc$occupancy_weighted_effective_lineages_simpson <-
    acc$occupancy_weighted_effective_lineages_simpson +
    w * diversity$occupancy_weighted_effective_lineages_simpson
  acc$occupancy_weighted_diversity_defined <-
    acc$occupancy_weighted_diversity_defined + w * as.numeric(diversity$diversity_defined)
  acc$niche_response_dispersion <- acc$niche_response_dispersion + w * nrd
  acc
}

write_map <- function(df, metric, out_dir, limits, legend_title) {
  safe_dir_create(out_dir)
  tm_value <- if ("time_ma" %in% names(df)) unique(df$time_ma)[1] else NA_real_
  time_slug <- if (is.finite(tm_value)) slug_time(tm_value) else "all_times"
  if (!is.finite(limits[2]) || limits[2] <= limits[1]) {
    limits[2] <- limits[1] + 1
  }
  # Windows' legacy path limit can silently prevent `ggsave()` from creating
  # PNGs below a descriptive Case05 output root.  Metric identity remains in
  # the map index; on-disk stems use compact stable codes.
  metric_code <- unname(metric_map_code[metric]) %||%
    gsub("[^A-Za-z0-9]+", "", metric)
  png <- file.path(out_dir, paste0(metric_code, "_", time_slug, ".png"))
  tif <- file.path(out_dir, paste0(metric_code, "_", time_slug, ".tif"))
  endpoint_note <- if (is.finite(tm_value) && tm_value == 0 &&
                       !identical(terminal_anchor_mode, "none")) {
    paste0(
      "; 0 Ma endpoint ensemble = ", terminal_anchor_mode,
      " (HMSC weight ", terminal_anchor_weight,
      "); historical maps remain unanchored"
    )
  } else ""
  if (make_png) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = lon, y = lat, fill = .data[[metric]])) +
      ggplot2::geom_tile(width = 1, height = 1) +
      ggplot2::coord_equal(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
      ggplot2::scale_fill_gradientn(
        # A valid zero must not look like an omitted grid cell or ocean.
        colours = c("#fff7bc", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
        limits = limits, na.value = "#6c757d", name = legend_title
      ) +
      ggplot2::labs(
      title = paste0(metric, ", ", time_slug),
        subtitle = paste0(
          "Global-grid forward scenario; transport = ", transport_mode,
          "; fixed legend range; P2 probability reporting window = ",
          reporting_interval_myr, " Myr", endpoint_note
        ),
        x = "Longitude", y = "Latitude",
        caption = "All valid palaeo-land cells are mapped; pale yellow is a valid zero and grey is ocean or unavailable habitat."
      ) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(fill = "#6c757d", colour = NA)
      )
    ggplot2::ggsave(png, p, width = 10, height = 5.4, dpi = 180)
    if (!file.exists(png)) {
      stop("PNG map was not written: ", png,
           ". Shorten --output or inspect the graphics device.", call. = FALSE)
    }
  }
  if (make_tif && requireNamespace("terra", quietly = TRUE)) {
    xyz <- df[, c("lon", "lat", metric), drop = FALSE]
    names(xyz) <- c("x", "y", "z")
    r <- terra::rast(xyz, type = "xyz", crs = "EPSG:4326")
    terra::writeRaster(r, tif, overwrite = TRUE)
    if (!file.exists(tif)) {
      stop("GeoTIFF map was not written: ", tif, call. = FALSE)
    }
  }
  data.frame(time_ma = tm_value, metric = metric,
             png = if (file.exists(png)) png else NA_character_,
             tif = if (file.exists(tif)) tif else NA_character_,
             fixed_legend_range = paste(limits, collapse = " to "),
             stringsAsFactors = FALSE)
}

config <- data.frame(
  parameter = c("analysis", "quick", "time_selection", "n_time_slices_requested",
                "n_available_time_slices", "n_species", "n_draws",
                "n_projected_times", "root_age_ma", "internal_dt",
                "n_map_times", "metric_storage",
                "response_draw_ids_requested",
                "max_hazard_per_step", "persistence_reference_myr",
                "reporting_interval_myr", "transport_mode",
                "plate_transport_index", "transport_min_target_coverage",
                "hazard_rate_quantile", "hazard_active_occupancy_threshold",
                "diagnostic_points_per_time",
                "checkpoint_every_draw", "resume_checkpoint",
                "palaeo_input_mode", "palaeo_slice_index",
                "topographic_slice_index", "cache_movement_matrices",
                "cache_movement_in_memory", "movement_kernel_cache_dir",
                "root_rho", "root_max_cells", "root_initialisation",
                 "root_seed_count", "root_support_power", "movement_intercept",
                 "local_dispersal_scale_km", "topographic_resistance",
                "topographic_connectivity_dir", "ldd_status", "ldd_per_source", "ldd_intercept",
                "ldd_distance_decay", "dispersal_scheme", "particles_per_source",
                "particle_seed", "run_biotic_filtering",
                "biotic_competition_strength",
                "biotic_facilitation_strength",
                "biotic_persistence_strength",
                 "run_speciation_demo",
                 "run_lineage_extinction_demo", "occupancy_mass_tolerance",
                 "write_modern_tip_occupancy", "terminal_anchor_mode",
                 "terminal_anchor_weight", "terminal_anchor_prevalence_tolerance",
                 "environmental_support_boundary", "run_type"),
  value = c("case04_plant200_global_dynamic_no_bgb_final", quick,
            time_selection, n_time_slices, length(available_inside_times),
            length(unique(tip_beta$species)), length(draw_ids),
            length(inside_times), root_age, internal_dt, length(map_times), metric_storage, response_draw_ids_arg,
            max_hazard_per_step,
            persistence_reference_myr, reporting_interval_myr, transport_mode,
            paths$plate_transport_index, transport_min_target_coverage,
            hazard_rate_quantile,
            hazard_active_occupancy_threshold, diagnostic_points_per_time,
            checkpoint_every_draw, resume_checkpoint, palaeo_input_mode,
            paths$palaeo_slice_index, paths$topographic_slice_index,
            cache_movement_matrices, cache_movement_in_memory, movement_kernel_cache_dir,
            root_rho,
            root_max_cells, root_initialisation, root_seed_count,
            root_support_power, movement_intercept, local_dispersal_scale_km,
            paths$topographic_resistance, paths$topographic_connectivity_dir,
            "disabled_random_global_LDD_removed_requires_precomputed_auditable_kernel",
            ldd_per_source, ldd_intercept, ldd_distance_decay, dispersal_scheme,
            particles_per_source, particle_seed,
            run_biotic_filtering, biotic_competition_strength,
             biotic_facilitation_strength, biotic_persistence_strength,
             run_speciation_demo, run_lineage_extinction_demo,
             occupancy_mass_tolerance,
             write_modern_tip_occupancy, terminal_anchor_mode,
             terminal_anchor_weight,
             paste(terminal_anchor_prevalence_tolerance, collapse = ","),
             "relative_environmental_support_only; intercept-free ancestral response; not a historical occupancy probability or an expected modern richness",
             "forward_scenario_matrix_streaming_not_particle_smoothing"),
  stringsAsFactors = FALSE
)
config <- rbind(
  config,
  data.frame(
    parameter = c("process_parameter_history", "process_parameter_mode"),
    value = c(paths$process_parameter_history, process_parameter_mode),
    stringsAsFactors = FALSE
  ),
  data.frame(
    parameter = c("write_attribution_state", "attribution_max_cells",
                  "attribution_max_lineages", "attribution_max_rows"),
    value = c(write_attribution_state, attribution_max_cells,
              attribution_max_lineages, attribution_max_rows),
    stringsAsFactors = FALSE
  )
)
write_csv(config, file.path(dirs$config, "case04_final_config.csv"))
write_csv(
  data.frame(
    process_or_layer = c(
      "environmental_filtering", "dispersal", "colonisation", "persistence",
      "biotic_filtering", "speciation_demo_layer", "lineage_extinction_demo_layer",
      "relative_environmental_support"
    ),
    status = c(
      "run", "run", "run", "run",
      if (run_biotic_filtering) "run" else "not_run",
      if (run_speciation_demo) "run" else "not_run",
      if (run_lineage_extinction_demo) "run" else "not_run",
      "diagnostic_only"
    ),
    interpretation = c(
      "ancestral response shape projected through palaeoenvironment",
      "global cell graph movement after plate transport",
      "arrival times environment-dependent establishment",
      "environment, rescue and geography-constrained local retention",
      "requires independently supported interaction data or an explicit scenario",
      "optional demonstration diagnostic; dated-tree lineage identity always remains active",
      "optional demonstration diagnostic; total lineage extinction is not estimated from an extant-only tree",
      "intercept-free response score; not occupancy probability or expected modern richness"
    ),
    stringsAsFactors = FALSE
  ),
  file.path(dirs$config, "case04_process_status.csv")
)

active_counts <- aggregate(lineage ~ time_ma + response_draw, responses,
                           function(x) length(unique(x)))
names(active_counts)[3] <- "n_active_lineages"
write_csv(active_counts, file.path(dirs$tree, "active_lineage_counts.csv"))

accumulators <- list()
# Keep numerical summaries for every actual time slice. `map_times` only
# controls raster/PNG rendering after the complete dynamic run has finished.
for (tm in inside_times) {
  accumulators[[as.character(tm)]] <- new_accumulator(earth_by_time[[as.character(tm)]])
}
map_index <- list()
movement_summary <- list()
transport_summary <- list()
speciation_summary <- list()
root_rows <- list()
modern_pred_acc <- NULL
modern_pred_weight <- 0
attribution_state_path <- file.path(dirs$inference,
                                    "attribution_interval_state.csv")
attribution_meta_path <- file.path(dirs$inference,
                                   "attribution_state_metadata.csv")
attribution_rows_written <- 0L

all_cell_ref <- do.call(rbind, lapply(earth_by_time, function(x) {
  x[, c("cell_id", "lon", "lat"), drop = FALSE]
}))
all_cell_ref <- all_cell_ref[!duplicated(all_cell_ref$cell_id), , drop = FALSE]
full_attr_targets <- c("occupancy", "range_area", "lineage_richness")
full_attr_components <- c("persistent_occupancy", "new_colonisation",
                          "recolonisation", "local_loss",
                          "colonisation_gain", "net_change_model")
full_attr_time_series <- list()
full_attr_time_idx <- 0L
full_attr_map <- all_cell_ref
full_attr_map$integration_weight <- 0
for (comp in full_attr_components) {
  full_attr_map[[comp]] <- 0
}

# A full Plant200 run is intentionally long.  Checkpoint only between complete
# posterior draws, never within a draw, so every restored state is an exact
# weighted accumulation of finished HMSC/ancestral-response histories.
checkpoint_path <- file.path(dirs$inference, "case04_checkpoint_after_draw.rds")
draw_task_path <- file.path(dirs$inference, "case04_posterior_draw_tasks.csv")
checkpoint_signature <- list(
  response_history_md5 = unname(tools::md5sum(paths$response_history)),
  tree_md5 = unname(tools::md5sum(paths$tree)),
  process_parameter_history_md5 = if (nzchar(paths$process_parameter_history)) {
    unname(tools::md5sum(paths$process_parameter_history))
  } else NA_character_,
  process_parameter_mode = process_parameter_mode,
  palaeo_input_mode = palaeo_input_mode,
  palaeo_earth_md5 = unname(tools::md5sum(if (use_streaming_slices) {
    paths$palaeo_slice_index
  } else paths$palaeo_earth_state)),
  topographic_input_md5 = unname(tools::md5sum(if (use_streaming_slices) {
    paths$topographic_slice_index
  } else paths$topographic_resistance)),
  response_draws = as.character(draw_ids),
  projected_times = as.numeric(inside_times),
  dispersal_scheme = dispersal_scheme,
  transport_mode = transport_mode,
  plate_transport_index_md5 = if (identical(transport_mode, "paleomap_h3_targetcentred")) {
    unname(tools::md5sum(paths$plate_transport_index))
  } else NA_character_,
  reporting_interval_myr = reporting_interval_myr,
  max_hazard_per_step = max_hazard_per_step,
  persistence_reference_myr = persistence_reference_myr,
  root_initialisation = root_initialisation,
  root_seed_count = root_seed_count,
  root_support_power = root_support_power,
  root_rho = root_rho,
  root_max_cells = root_max_cells,
  movement_intercept = movement_intercept,
  local_dispersal_scale_km = local_dispersal_scale_km,
  establishment_intercept = establishment_intercept,
  establishment_slope = establishment_slope,
  persistence_intercept = persistence_intercept,
  persistence_slope = persistence_slope,
  terminal_anchor_mode = terminal_anchor_mode,
  terminal_anchor_weight = terminal_anchor_weight,
  terminal_anchor_prevalence_tolerance = terminal_anchor_prevalence_tolerance,
  run_biotic_filtering = run_biotic_filtering,
  run_speciation_demo = run_speciation_demo,
  run_lineage_extinction_demo = run_lineage_extinction_demo,
  root_seed = seed,
  particle_seed = particle_seed,
  particles_per_source = particles_per_source,
  hazard_rate_quantile = hazard_rate_quantile,
  hazard_active_occupancy_threshold = hazard_active_occupancy_threshold,
  diagnostic_points_per_time = diagnostic_points_per_time
)
completed_draws <- character()
if (resume_checkpoint && file.exists(checkpoint_path)) {
  checkpoint <- readRDS(checkpoint_path)
  if (!identical(checkpoint$signature, checkpoint_signature)) {
    stop("Existing checkpoint is incompatible with this Case04 configuration. ",
         "Use a new output directory or set --resume_checkpoint=false.", call. = FALSE)
  }
  completed_draws <- as.character(checkpoint$completed_draws)
  accumulators <- checkpoint$accumulators
  movement_summary <- checkpoint$movement_summary
  transport_summary <- checkpoint$transport_summary
  speciation_summary <- checkpoint$speciation_summary
  root_rows <- checkpoint$root_rows
  modern_pred_acc <- checkpoint$modern_pred_acc
  modern_pred_weight <- checkpoint$modern_pred_weight
  attribution_rows_written <- checkpoint$attribution_rows_written
  full_attr_time_series <- checkpoint$full_attr_time_series
  full_attr_time_idx <- checkpoint$full_attr_time_idx
  full_attr_map <- checkpoint$full_attr_map
  terminal_anchor_calibration_rows <- checkpoint$terminal_anchor_calibration_rows %||% list()
  log_msg("Resuming after completed response draw(s): ",
          paste(completed_draws, collapse = ", "))
}
if (!all(completed_draws %in% draw_ids)) {
  stop("Checkpoint lists posterior draws that are absent from the current response history.",
       call. = FALSE)
}
draw_tasks <- data.frame(
  response_draw = as.character(draw_ids),
  status = ifelse(as.character(draw_ids) %in% completed_draws, "complete", "pending"),
  checkpoint_path = checkpoint_path,
  stringsAsFactors = FALSE
)
write_csv(draw_tasks, draw_task_path)
write_draw_checkpoint <- function() {
  if (!checkpoint_every_draw) return(invisible(FALSE))
  payload <- list(
    signature = checkpoint_signature,
    completed_draws = completed_draws,
    accumulators = accumulators,
    movement_summary = movement_summary,
    transport_summary = transport_summary,
    speciation_summary = speciation_summary,
    root_rows = root_rows,
    modern_pred_acc = modern_pred_acc,
    modern_pred_weight = modern_pred_weight,
    attribution_rows_written = attribution_rows_written,
    full_attr_time_series = full_attr_time_series,
    full_attr_time_idx = full_attr_time_idx,
    full_attr_map = full_attr_map,
    terminal_anchor_calibration_rows = terminal_anchor_calibration_rows
  )
  saveRDS(payload, checkpoint_path, compress = "gzip")
  invisible(TRUE)
}

transition_components <- function(q0, lambda_c, lambda_l, delta_t) {
  lambda_c[!is.finite(lambda_c) | lambda_c < 0] <- 0
  lambda_l[!is.finite(lambda_l) | lambda_l < 0] <- 0
  r <- lambda_c + lambda_l
  exp_r <- exp(-r * delta_t)
  p01 <- matrix(0, nrow(q0), ncol(q0), dimnames = dimnames(q0))
  p10 <- matrix(0, nrow(q0), ncol(q0), dimnames = dimnames(q0))
  p11 <- matrix(1, nrow(q0), ncol(q0), dimnames = dimnames(q0))
  ok <- is.finite(r) & r > 0
  p01[ok] <- lambda_c[ok] / r[ok] * (1 - exp_r[ok])
  p10[ok] <- lambda_l[ok] / r[ok] * (1 - exp_r[ok])
  p11[ok] <- lambda_c[ok] / r[ok] + lambda_l[ok] / r[ok] * exp_r[ok]
  p01 <- clip01(p01)
  p10 <- clip01(p10)
  p11 <- clip01(p11)
  exp_l <- exp(-lambda_l * delta_t)
  persistent <- q0 * exp_l
  new_col <- (1 - q0) * p01
  recol <- q0 * pmax(p11 - exp_l, 0)
  local_loss <- q0 * p10
  list(
    persistent_occupancy = persistent,
    new_colonisation = new_col,
    recolonisation = recol,
    local_loss = local_loss,
    colonisation_gain = new_col + recol,
    net_change_model = persistent + new_col + recol - q0
  )
}

append_full_attribution_budget <- function(draw_id, time_start, time_end,
                                           e_tm, q_start, lambda_c,
                                           lambda_l, draw_weight) {
  comps <- transition_components(q_start, lambda_c, lambda_l,
                                 abs(time_start - time_end))
  area <- suppressWarnings(as.numeric(e_tm$land_area_km2))
  area[!is.finite(area) | area <= 0] <- 1
  area_denom <- sum(area, na.rm = TRUE)
  if (!is.finite(area_denom) || area_denom <= 0) area_denom <- 1
  dt <- abs(time_start - time_end)
  rows <- vector("list", length(full_attr_targets) * length(full_attr_components))
  ii <- 0L
  for (target in full_attr_targets) {
    for (comp in full_attr_components) {
      mat <- comps[[comp]]
      weighted_sum <- sum(mat * area, na.rm = TRUE)
      value <- if (target %in% c("occupancy", "lineage_richness")) {
        weighted_sum / area_denom
      } else {
        weighted_sum
      }
      ii <- ii + 1L
      rows[[ii]] <- data.frame(
        posterior_draw_id = draw_id,
        target = target,
        time_start_ma = time_start,
        time_end_ma = time_end,
        delta_t = dt,
        component = comp,
        value = value,
        interval_weight = dt * draw_weight,
        n_cells = nrow(q_start),
        n_lineages = ncol(q_start),
        spatial_scope = "all_valid_land_cells",
        taxon_scope = "all_active_lineages",
        stringsAsFactors = FALSE
      )
    }
  }
  full_attr_time_idx <<- full_attr_time_idx + 1L
  full_attr_time_series[[full_attr_time_idx]] <<- do.call(rbind, rows)

  map_idx <- match(e_tm$cell_id, full_attr_map$cell_id)
  valid <- !is.na(map_idx)
  if (any(valid)) {
    w <- dt * draw_weight
    full_attr_map$integration_weight[map_idx[valid]] <<-
      full_attr_map$integration_weight[map_idx[valid]] + w
    for (comp in full_attr_components) {
      full_attr_map[[comp]][map_idx[valid]] <<-
        full_attr_map[[comp]][map_idx[valid]] +
        w * rowSums(comps[[comp]][valid, , drop = FALSE], na.rm = TRUE)
    }
  }
  invisible(TRUE)
}

sample_indices <- function(n, max_n, seed_offset = 0L) {
  if (!is.finite(max_n) || max_n <= 0L || n <= max_n) return(seq_len(n))
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed + as.integer(seed_offset))
  sort(sample(seq_len(n), max_n))
}

stable_seed_offset <- function(x) {
  x <- paste(as.character(x), collapse = "_")
  ints <- utf8ToInt(x)
  if (!length(ints)) return(0L)
  as.integer(sum(ints * seq_along(ints)) %% 100000L)
}

append_attribution_interval_state <- function(draw_id, time_start, time_end,
                                              e_tm, r_tm, q_start, q_end,
                                              lambda_c, lambda_l, eta,
                                              arrival_hazard) {
  if (!write_attribution_state) return(invisible(FALSE))
  if (attribution_rows_written >= attribution_max_rows) return(invisible(FALSE))
  n_cell <- nrow(q_start)
  n_lin <- ncol(q_start)
  draw_offset <- stable_seed_offset(draw_id)
  cell_idx <- sample_indices(n_cell, attribution_max_cells,
                             seed_offset = round(time_end * 10) + draw_offset)
  lin_idx <- sample_indices(n_lin, attribution_max_lineages,
                            seed_offset = round(time_start * 10) + draw_offset + 10000L)
  n_new <- length(cell_idx) * length(lin_idx)
  remaining <- attribution_max_rows - attribution_rows_written
  if (n_new > remaining) {
    keep_flat <- sample_indices(n_new, remaining,
                                seed_offset = round(time_start * 100) + draw_offset + 20000L)
  } else {
    keep_flat <- seq_len(n_new)
  }
  grid <- expand.grid(cell_i = cell_idx, lin_i = lin_idx,
                      KEEP.OUT.ATTRS = FALSE)
  grid <- grid[keep_flat, , drop = FALSE]
  if (!nrow(grid)) return(invisible(FALSE))
  area <- suppressWarnings(as.numeric(e_tm$land_area_km2[grid$cell_i]))
  area[!is.finite(area) | area <= 0] <- 1
  out <- data.frame(
    posterior_draw_id = draw_id,
    lineage_id = colnames(q_start)[grid$lin_i],
    cell_id = rownames(q_start)[grid$cell_i],
    time_start_ma = time_start,
    time_end_ma = time_end,
    q0 = q_start[cbind(grid$cell_i, grid$lin_i)],
    q_observed_next = q_end[cbind(grid$cell_i, grid$lin_i)],
    lambda_C = lambda_c[cbind(grid$cell_i, grid$lin_i)],
    lambda_L = lambda_l[cbind(grid$cell_i, grid$lin_i)],
    cell_area_km2 = area,
    eta_env = eta[cbind(grid$cell_i, grid$lin_i)],
    arrival_hazard = arrival_hazard[cbind(grid$cell_i, grid$lin_i)],
    attribution_sampling = if (n_cell > length(cell_idx) ||
                               n_lin > length(lin_idx) ||
                               n_new > remaining) {
      "deterministic_sample_for_process_attribution"
    } else {
      "complete_for_selected_case05_run"
    },
    stringsAsFactors = FALSE
  )
  append_csv(out, attribution_state_path)
  attribution_rows_written <<- attribution_rows_written + nrow(out)
  invisible(TRUE)
}

site_comm <- NULL
if (file.exists(paths$comm)) {
  site_comm <- utils::read.csv(paths$comm, check.names = FALSE, stringsAsFactors = FALSE)
}
site_coordinates <- NULL
if (file.exists(paths$sites)) {
  site_coordinates <- utils::read.csv(paths$sites, check.names = FALSE,
                                      stringsAsFactors = FALSE)
}

# The calibrated terminal ensemble is deliberately tied to a stated modern
# calibration subset.  The dynamic state itself is still computed on the full
# palaeo-land grid; these indices are used only to choose the *0 Ma* endpoint
# ensemble weight, never to subset or edit pre-0 Ma occupancy.
terminal_anchor_site_match <- NULL
# Keep calibration records restored from a completed-draw checkpoint.  A
# report-only resume must never erase the realised weight selected during the
# original full-grid trajectory.
if (!exists("terminal_anchor_calibration_rows", inherits = FALSE)) {
  terminal_anchor_calibration_rows <- list()
}
terminal_anchor_calibration_index <- function(e_tm) {
  if (!identical(terminal_anchor_mode,
                 "hmsc_fixed_effect_prevalence_calibrated")) return(NULL)
  if (!is.null(terminal_anchor_site_match)) return(terminal_anchor_site_match)
  if (is.null(site_comm) || is.null(site_coordinates) ||
      !all(c("site_id") %in% names(site_comm)) ||
      !all(c("site_id", "lon", "lat") %in% names(site_coordinates))) {
    stop(
      "Calibrated terminal anchoring requires comm.csv site_id and sites.csv ",
      "site_id/lon/lat columns.", call. = FALSE
    )
  }
  if (anyDuplicated(site_comm$site_id) || anyDuplicated(site_coordinates$site_id) ||
      !setequal(site_comm$site_id, site_coordinates$site_id)) {
    stop(
      "Calibrated terminal anchoring requires unique and identical comm/sites ",
      "site_id values.", call. = FALSE
    )
  }
  coords <- site_coordinates[match(site_comm$site_id,
                                   site_coordinates$site_id), , drop = FALSE]
  coord_key <- function(lon, lat) {
    paste(formatC(as.numeric(lon), format = "f", digits = 6),
          formatC(as.numeric(lat), format = "f", digits = 6), sep = "_")
  }
  grid_key <- coord_key(e_tm$lon, e_tm$lat)
  site_key <- coord_key(coords$lon, coords$lat)
  index <- match(site_key, grid_key)
  nearest <- function(lon, lat) {
    lon_delta <- abs(e_tm$lon - lon)
    lon_delta <- pmin(lon_delta, 360 - lon_delta)
    which.min(lon_delta^2 + (e_tm$lat - lat)^2)
  }
  missing_index <- which(is.na(index))
  if (length(missing_index)) {
    index[missing_index] <- vapply(
      missing_index,
      function(i) nearest(coords$lon[i], coords$lat[i]), integer(1)
    )
  }
  lon_delta <- abs(e_tm$lon[index] - coords$lon)
  lon_delta <- pmin(lon_delta, 360 - lon_delta)
  distance_deg <- sqrt(lon_delta^2 + (e_tm$lat[index] - coords$lat)^2)
  if (any(!is.finite(distance_deg) | distance_deg > 0.51)) {
    stop(
      "At least one modern calibration site is farther than 0.51 degrees from ",
      "the 0 Ma prediction grid.", call. = FALSE
    )
  }
  terminal_anchor_site_match <<- data.frame(
    site_id = site_comm$site_id,
    endpoint_cell_id = e_tm$cell_id[index],
    source_lon = coords$lon,
    source_lat = coords$lat,
    endpoint_lon = e_tm$lon[index],
    endpoint_lat = e_tm$lat[index],
    match_distance_deg = distance_deg,
    stringsAsFactors = FALSE
  )
  write_csv(terminal_anchor_site_match,
            file.path(dirs$modern, "terminal_anchor_site_grid_match.csv"))
  terminal_anchor_site_match
}

# The time-specific land graph and its P2 rate kernel are identical for every
# HMSC/ancestral-response draw.  Cache it on disk by time slice: full-grid
# calculation stays exact, while old slices do not remain in memory across all
# response draws.
movement_cache <- new.env(parent = emptyenv())
movement_kernel_cache_file <- function(tm) {
  file.path(movement_kernel_cache_dir,
            paste0("movement_", dispersal_scheme, "_", slug_time(tm), ".rds"))
}
movement_matrix_for_time <- function(e_tm, tm, delta_t) {
  key <- paste0("time_", format(tm, trim = TRUE, scientific = FALSE))
  if (cache_movement_in_memory && exists(key, envir = movement_cache,
                                         inherits = FALSE)) {
    return(get(key, envir = movement_cache, inherits = FALSE))
  }
  cache_file <- movement_kernel_cache_file(tm)
  M <- if (cache_movement_matrices && file.exists(cache_file)) {
    readRDS(cache_file)
  } else {
    built <- build_movement_matrix(e_tm, tm, delta_t)
    if (cache_movement_matrices) {
      # Independent posterior-draw shards can request the same time-specific
      # kernel simultaneously. Write a complete temporary RDS then publish it
      # atomically, so a reader never observes a partially written sparse matrix.
      tmp <- tempfile(pattern = "movement_kernel_", tmpdir = dirname(cache_file),
                      fileext = ".rds")
      saveRDS(built, tmp, compress = "gzip")
      if (!file.exists(cache_file)) file.rename(tmp, cache_file)
      if (file.exists(tmp)) unlink(tmp, force = TRUE)
      if (!file.exists(cache_file)) {
        stop("Unable to publish shared movement kernel cache: ", cache_file,
             call. = FALSE)
      }
      readRDS(cache_file)
    } else {
      built
    }
  }
  if (cache_movement_in_memory) assign(key, M, envir = movement_cache)
  M
}

for (dd in draw_ids) {
  if (dd %in% completed_draws) {
    log_msg("Skipping checkpointed response draw", dd)
    next
  }
  draw_tasks$status[draw_tasks$response_draw == dd] <- "running"
  write_csv(draw_tasks, draw_task_path)
  log_msg("Running response draw", dd)
  prev_Q <- NULL
  prev_cells <- NULL
  prev_nodes <- NULL
  prev_time <- NA_real_
  for (tm in inside_times) {
    # This is always the complete valid 1 degree palaeo-land grid at tm.  In
    # streaming mode the prior time slice is no longer retained in RAM.
    e_tm <- load_earth_time(tm)
    r_tm <- responses[responses$response_draw == dd & responses$time_ma == tm, , drop = FALSE]
    if (!nrow(e_tm) || !nrow(r_tm)) next
    r_tm <- r_tm[order(r_tm$node), , drop = FALSE]
    process_parameters <- process_parameters_for(
      lineages = as.character(r_tm$lineage), time_ma = tm,
      response_draw = dd
    )
    es <- calc_eta_s(e_tm, r_tm)
    eta <- es$eta
    S <- es$S
    colnames(S) <- r_tm$lineage
    colnames(eta) <- r_tm$lineage
    rownames(S) <- e_tm$cell_id
    rownames(eta) <- e_tm$cell_id
    process_matrix <- function(x) {
      matrix(rep(unname(x[colnames(eta)]), each = nrow(eta)),
             nrow = nrow(eta), ncol = ncol(eta), dimnames = dimnames(eta))
    }
    establishment <- matrix(
      hee_colonisation_establishment_probability(
        as.numeric(eta),
        intercept = as.numeric(process_matrix(process_parameters$establishment_intercept)),
        slope = as.numeric(process_matrix(process_parameters$establishment_slope))
      ), nrow(eta), ncol(eta), dimnames = dimnames(eta)
    )
    persistence <- matrix(
      hee_persistence_probability(
        as.numeric(eta),
        intercept = as.numeric(process_matrix(process_parameters$persistence_intercept)),
        slope = as.numeric(process_matrix(process_parameters$persistence_slope))
      ), nrow(eta), ncol(eta), dimnames = dimnames(eta)
    )
    biotic <- calc_biotic_filtering(
      matrix(0, nrow(S), ncol(S), dimnames = dimnames(S)), S
    )
    speciation_footprint <- numeric(nrow(S))
    lineage_extinction_pressure <- numeric(nrow(S))
    area <- suppressWarnings(as.numeric(e_tm$land_area_km2))
    area[!is.finite(area) | area <= 0] <- 1

    if (is.null(prev_Q)) {
      Q <- init_root_q(S, e_tm, tm, dd)
      arrival_hazard <- matrix(0, nrow(Q), ncol(Q), dimnames = dimnames(Q))
      lambda_c <- matrix(0, nrow(Q), ncol(Q), dimnames = dimnames(Q))
      root_rows[[length(root_rows) + 1L]] <- data.frame(
        response_draw = dd,
        lineage = rep(colnames(Q), each = nrow(Q)),
        cell_id = rep(rownames(Q), times = ncol(Q)),
        time_ma = tm,
        q_root = as.vector(Q),
        root_prior = attr(Q, "root_initialisation"),
        stringsAsFactors = FALSE
      )
      root_rows[[length(root_rows)]] <- root_rows[[length(root_rows)]][
        root_rows[[length(root_rows)]]$q_root > 0, , drop = FALSE
      ]
    } else {
      dt <- abs(prev_time - tm)
      cur_nodes <- as.integer(r_tm$node)
      prev_match <- map_to_previous(cur_nodes, prev_nodes)
      Q_source <- matrix(0, nrow(prev_Q), nrow(r_tm),
                         dimnames = list(prev_cells, r_tm$lineage))
      for (j in seq_along(prev_match)) {
        if (!is.na(prev_match[j])) {
          Q_source[, j] <- prev_Q[, prev_match[j]]
        }
      }
      transport_object <- transport_for_interval(prev_time, tm, prev_cells,
                                                 e_tm$cell_id)
      carried <- hee_plate_carry_occupancy(
        Q_source, transport_object,
        source_cells = prev_cells, target_cells = e_tm$cell_id
      )
      Q <- carried$occupancy
      td <- carried$diagnostics
      if (identical(transport_mode, "paleomap_h3_targetcentred") &&
          td$target_coverage_fraction[[1L]] < transport_min_target_coverage) {
        stop("PALEOMAP transport target coverage is below the configured minimum at ",
             prev_time, " -> ", tm, " Ma: ",
             signif(td$target_coverage_fraction[[1L]], 4), call. = FALSE)
      }
      transport_summary[[length(transport_summary) + 1L]] <- data.frame(
        response_draw = dd,
        time_from_ma = prev_time,
        time_to_ma = tm,
        transport_mode = transport_mode,
        n_source_cells = td$n_source_cells,
        n_target_cells = td$n_target_cells,
        n_links = td$n_links,
        source_coverage_fraction = td$source_coverage_fraction,
        target_coverage_fraction = td$target_coverage_fraction,
        n_unmapped_source_cells = td$n_unmapped_source_cells,
        n_unmapped_target_cells = td$n_unmapped_target_cells,
        transport_interpretation = carried$transport$interpretation,
        stringsAsFactors = FALSE
      )
      Q_start_interval <- Q
      split_parent <- !is.na(prev_match) &
        (duplicated(prev_match) | duplicated(prev_match, fromLast = TRUE))
      speciation_summary[[length(speciation_summary) + 1L]] <- data.frame(
        response_draw = dd,
        time_ma = tm,
        previous_time_ma = prev_time,
        n_previous_lineages = length(prev_nodes),
        n_current_lineages = length(cur_nodes),
        n_copy_then_diverge_lineages = sum(split_parent),
        n_parent_lineages_split = length(unique(prev_match[split_parent])),
        process = if (any(split_parent)) "speciation_copy_then_diverge" else "no_tree_split_in_interval",
        interpretation = "Tree-conditioned lineage identity update; not an automatically estimated speciation rate.",
        stringsAsFactors = FALSE
      )
      if (run_speciation_demo && length(prev_match)) {
        if (any(split_parent)) {
          speciation_footprint <- clip01(rowMeans(Q[, split_parent, drop = FALSE],
                                                  na.rm = TRUE))
        }
      }
      M <- movement_matrix_for_time(e_tm, tm, dt)
      source_movement_multiplier <- process_matrix(
        process_parameters$movement_multiplier
      )
      # A palaeo-environment interval can be five or more Myr, whereas the
      # land graph permits one neighbouring-cell transition at a time.  Choose
      # substeps from both the user maximum and a CTMC hazard bound.  This
      # prevents a numerical one-edge-per-5-Myr propagation ceiling.
      lambda_l_reference <- -log(pmax(persistence, 1e-12)) /
        persistence_reference_myr
      max_incoming_rate <- attr(M, "max_incoming_rate")
      if (!is.finite(max_incoming_rate)) max_incoming_rate <- 0
      active_loss_rates <- lambda_l_reference[
        Q_start_interval > hazard_active_occupancy_threshold &
          is.finite(lambda_l_reference)
      ]
      # If a newly created lineage has not yet inherited any non-zero cell,
      # use all finite cells as a conservative fallback.  Otherwise, an
      # extreme loss rate in a cell with q approximately zero does not force
      # global numerical refinement; that cell still receives its exact CTMC
      # local-loss transition at every chosen step.
      if (!length(active_loss_rates)) {
        active_loss_rates <- lambda_l_reference[is.finite(lambda_l_reference)]
      }
      loss_rate_bound <- if (length(active_loss_rates)) {
        as.numeric(stats::quantile(active_loss_rates,
                                   probs = hazard_rate_quantile,
                                   names = FALSE, na.rm = TRUE,
                                   type = 8))
      } else 0
      if (!is.finite(loss_rate_bound) || loss_rate_bound < 0) loss_rate_bound <- 0
      # `max_incoming_rate` assumes every potential source is occupied with
      # probability one.  That is a graph property, not the arrival hazard in
      # the current dynamic state, and it badly over-refines a sparse root or
      # colonisation front.  Use the corresponding upper quantile of actual
      # arrival pressure; the following CTMC calls still use the full exact
      # sparse matrix at every substep.
      arrival_reference <- as.matrix(
        M %*% (Q_start_interval * source_movement_multiplier)
      )
      active_arrival_rates <- arrival_reference[
        (Q_start_interval > hazard_active_occupancy_threshold |
           arrival_reference > hazard_active_occupancy_threshold) &
          is.finite(arrival_reference)
      ]
      if (!length(active_arrival_rates)) {
        active_arrival_rates <- arrival_reference[is.finite(arrival_reference)]
      }
      arrival_rate_bound <- if (length(active_arrival_rates)) {
        as.numeric(stats::quantile(active_arrival_rates,
                                   probs = hazard_rate_quantile,
                                   names = FALSE, na.rm = TRUE,
                                   type = 8))
      } else 0
      if (!is.finite(arrival_rate_bound) || arrival_rate_bound < 0) {
        arrival_rate_bound <- 0
      }
      max_transition_rate <- arrival_rate_bound + loss_rate_bound
      if (!is.finite(max_transition_rate) || max_transition_rate < 0) {
        max_transition_rate <- 0
      }
      n_steps <- max(
        1L,
        ceiling(dt / internal_dt),
        ceiling(dt * max_transition_rate / max_hazard_per_step)
      )
      h <- dt / n_steps
      movement_summary[[length(movement_summary) + 1L]] <- data.frame(
        response_draw = dd, time_ma = tm, n_cells = nrow(e_tm),
        n_edges = length(M@x), mean_lambda_movement = mean(M@x),
        max_lambda_movement = max(M@x), delta_t = dt,
        n_internal_steps = n_steps,
        internal_step_myr = h,
        max_hazard_per_step = max_hazard_per_step,
        persistence_reference_myr = persistence_reference_myr,
        hazard_rate_quantile = hazard_rate_quantile,
        hazard_active_occupancy_threshold = hazard_active_occupancy_threshold,
        n_active_loss_rate_cells = length(active_loss_rates),
        loss_rate_quantile_bound = loss_rate_bound,
        n_active_arrival_rate_cells = length(active_arrival_rates),
        arrival_rate_quantile_bound = arrival_rate_bound,
        max_source_total_rate = attr(M, "max_source_total_rate"),
        max_incoming_rate = max_incoming_rate,
        max_transition_rate_bound = max_transition_rate,
        mean_effective_cost_km = attr(M, "mean_effective_cost_km"),
        max_effective_cost_km = attr(M, "max_effective_cost_km"),
        mean_topographic_resistance = attr(M, "mean_topographic_resistance"),
        topographic_graph_file = attr(M, "topographic_graph_file"),
        dispersal_scheme = attr(M, "dispersal_scheme"),
        movement_cost_basis = attr(M, "movement_cost_basis"),
        particles_per_source = attr(M, "particles_per_source"),
        particle_seed = attr(M, "particle_seed"),
        ldd_status = attr(M, "ldd_status"),
        source_rate_definition = attr(M, "source_rate_definition"),
        process_parameter_mode = process_parameter_mode,
        mean_movement_multiplier = mean(source_movement_multiplier),
        mean_establishment_intercept = mean(process_parameters$establishment_intercept),
        mean_establishment_slope = mean(process_parameters$establishment_slope),
        mean_persistence_intercept = mean(process_parameters$persistence_intercept),
        mean_persistence_slope = mean(process_parameters$persistence_slope),
        movement_kernel = if (identical(dispersal_scheme, "particle_topographic")) {
          "hee_dispersal_edge_kernel() -> hee_dispersal_rate_kernel() -> hee_dispersal_particle_kernel()"
        } else "hee_dispersal_edge_kernel() -> hee_dispersal_rate_kernel()",
        stringsAsFactors = FALSE
      )
      lambda_l <- lambda_l_reference
      arrival_hazard <- matrix(0, nrow(Q), ncol(Q), dimnames = dimnames(Q))
      lambda_c <- matrix(0, nrow(Q), ncol(Q), dimnames = dimnames(Q))
      for (st in seq_len(n_steps)) {
        biotic <- calc_biotic_filtering(Q, S)
        step <- hee_projection_global_grid_step(
          occupancy = Q,
          movement_matrix = M,
          eta = eta,
          delta_t = h,
          reporting_interval_myr = reporting_interval_myr,
          persistence_reference_interval = persistence_reference_myr,
          establishment_intercept = process_parameters$establishment_intercept,
          establishment_slope = process_parameters$establishment_slope,
          persistence_intercept = process_parameters$persistence_intercept,
          persistence_slope = process_parameters$persistence_slope,
          movement_multiplier = process_parameters$movement_multiplier,
          colonisation_biotic = biotic$colonisation_modifier,
          persistence_biotic = biotic$persistence_modifier,
          habitat_state = e_tm$H_state,
          stochastic = FALSE
        )
        Q <- step$occupancy_next
        arrival_hazard <- step$arrival_hazard
        establishment <- step$establishment_probability
        persistence <- step$persistence_probability
        lambda_c <- step$colonisation_hazard
        lambda_l <- step$local_loss_hazard
      }
      append_attribution_interval_state(
        draw_id = dd,
        time_start = prev_time,
        time_end = tm,
        e_tm = e_tm,
        r_tm = r_tm,
        q_start = Q_start_interval,
        q_end = Q,
        lambda_c = lambda_c,
        lambda_l = lambda_l,
        eta = eta,
        arrival_hazard = arrival_hazard
      )
      append_full_attribution_budget(
        draw_id = dd,
        time_start = prev_time,
        time_end = tm,
        e_tm = e_tm,
        q_start = Q_start_interval,
        lambda_c = lambda_c,
        lambda_l = lambda_l,
        draw_weight = draw_weight[[dd]]
      )
      if (run_lineage_extinction_demo) {
        arrival_prob <- clip01(1 - exp(-arrival_hazard))
        lineage_extinction_pressure <- clip01(
          rowMeans(1 - persistence, na.rm = TRUE) *
            (1 - rowMeans(arrival_prob, na.rm = TRUE))
        )
      }
    }
    Q <- clip01(Q)
    forward_raw_Q <- Q
    tip_idx <- integer()
    tip_nowcast <- NULL
    endpoint_anchor_applied <- FALSE
    endpoint_anchor_weight_current <- 0
    endpoint_anchor_role <- "not_applied; forward dynamic scenario"
    if (tm == 0) {
      tip_idx <- match(tip_labels, colnames(Q))
      if (anyNA(tip_idx)) {
        stop("At 0 Ma, one or more modern tree tips are absent from the active lineage state.",
             call. = FALSE)
      }
      intercept_tab <- tip_hmsc_intercepts[
        tip_hmsc_intercepts$response_draw == dd,
        c("species", "hmsc_intercept_original"), drop = FALSE
      ]
      intercept <- intercept_tab$hmsc_intercept_original[
        match(tip_labels, intercept_tab$species)
      ]
      if (any(!is.finite(intercept))) {
        stop("Missing original HMSC intercept for one or more 0 Ma tips.",
             call. = FALSE)
      }
      tip_eta_relative <- eta[, tip_idx, drop = FALSE]
      tip_nowcast <- clip01(link_inverse(sweep(tip_eta_relative, 2, intercept, "+")))
      tip_nowcast[e_tm$H_state != 1, ] <- 0
      if (identical(terminal_anchor_mode, "hmsc_fixed_effect_convex")) {
        endpoint_anchor_weight_current <- terminal_anchor_weight
        Q[, tip_idx] <- HmscEcoEvo::hee_endpoint_nowcast_anchor(
          Q[, tip_idx, drop = FALSE], tip_nowcast,
          anchor_weight = endpoint_anchor_weight_current
        )
        endpoint_anchor_applied <- TRUE
        endpoint_anchor_role <- paste0(
          "0Ma convex endpoint model ensemble; weight=", endpoint_anchor_weight_current,
          "; fixed-effect HMSC nowcast is an endpoint constraint, not an ancestral niche or independent validation"
        )
      } else if (identical(terminal_anchor_mode,
                           "hmsc_fixed_effect_prevalence_calibrated")) {
        match_table <- terminal_anchor_calibration_index(e_tm)
        site_idx <- match(match_table$endpoint_cell_id, e_tm$cell_id)
        if (anyNA(site_idx)) {
          stop("Terminal-anchor site grid match could not be resolved at 0 Ma.",
               call. = FALSE)
        }
        calibration <- HmscEcoEvo::hee_endpoint_anchor_weight(
          Q[site_idx, tip_idx, drop = FALSE],
          tip_nowcast[site_idx, , drop = FALSE],
          prevalence_ratio_tolerance = terminal_anchor_prevalence_tolerance
        )
        endpoint_anchor_weight_current <- calibration$selected_anchor_weight[[1L]]
        Q[, tip_idx] <- HmscEcoEvo::hee_endpoint_nowcast_anchor(
          Q[, tip_idx, drop = FALSE], tip_nowcast,
          anchor_weight = endpoint_anchor_weight_current
        )
        endpoint_anchor_applied <- TRUE
        calibration$response_draw <- as.character(dd)
        calibration$time_ma <- 0
        calibration$calibration_subset <- "modern_site_grid; same-data-family training diagnostic unless held-out inputs are supplied"
        terminal_anchor_calibration_rows[[length(terminal_anchor_calibration_rows) + 1L]] <- calibration
        endpoint_anchor_role <- paste0(
          "0Ma prevalence-calibrated convex endpoint model ensemble; selected HMSC weight=",
          formatC(endpoint_anchor_weight_current, format = "f", digits = 6),
          " to place the modern-site-grid prevalence ratio in [",
          formatC(terminal_anchor_prevalence_tolerance[1L], format = "f", digits = 3),
          ", ",
          formatC(terminal_anchor_prevalence_tolerance[2L], format = "f", digits = 3),
          "]; fixed-effect HMSC nowcast is a training diagnostic unless independently held out"
        )
      }
    }

    accumulators[[as.character(tm)]] <- accumulate_time(
      accumulators[[as.character(tm)]], Q, S, arrival_hazard, lambda_c,
      persistence, eta, draw_weight[[dd]], area, reporting_interval_myr,
      biotic_competition = biotic$competition_pressure,
      biotic_facilitation = biotic$facilitation_pressure,
      speciation_footprint = speciation_footprint,
      lineage_extinction_pressure = lineage_extinction_pressure,
      forward_raw_Q = forward_raw_Q,
      terminal_anchor_weight = if (endpoint_anchor_applied) {
        endpoint_anchor_weight_current
      } else 0
    )

    if (tm == 0) {
      # The accumulator already contains the same posterior-draw weight used
      # for dynamic Q. This makes the final 0 Ma map a genuine draw-level
      # expectation, rather than a prediction from averaged coefficients. The
      # nowcast itself is retained regardless of whether it was used as an
      # explicitly labelled terminal anchor.
      accumulators[[as.character(tm)]]$hmsc_fixed_effect_nowcast_richness <-
        accumulators[[as.character(tm)]]$hmsc_fixed_effect_nowcast_richness +
        draw_weight[[dd]] * rowSums(tip_nowcast, na.rm = TRUE)
    }

    if (tm == 0 && write_modern_tip_occupancy) {
      draw_id_safe <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(dd))
      endpoint_payload <- list(
        response_draw = as.character(dd),
        time_ma = 0,
        cell_id = e_tm$cell_id,
        lon = e_tm$lon,
        lat = e_tm$lat,
        species = tip_labels,
        dynamic_occupancy_probability = Q[, tip_idx, drop = FALSE],
        forward_raw_occupancy_probability = forward_raw_Q[, tip_idx, drop = FALSE],
        hmsc_fixed_effect_nowcast_probability = tip_nowcast,
        hmsc_intercept_original = intercept,
        terminal_anchor_mode = terminal_anchor_mode,
        terminal_anchor_weight = if (endpoint_anchor_applied) {
          endpoint_anchor_weight_current
        } else 0,
        terminal_anchor_role = endpoint_anchor_role,
        nowcast_definition = paste(
          "0Ma fixed-effect HMSC nowcast with original tip intercept and locked",
          "environmental recipe; no site random effect; not an ancestral-niche intercept"
        )
      )
      saveRDS(
        endpoint_payload,
        file.path(dirs$modern,
                  paste0("modern_tip_endpoint_0Ma_", draw_id_safe, ".rds")),
        compress = "gzip"
      )
    }

    if (tm == 0 && !is.null(site_comm) && "site_id" %in% names(site_comm)) {
      sp <- intersect(colnames(Q), names(site_comm))
      if (length(sp)) {
        site_idx <- match(site_comm$site_id, rownames(Q))
        ok <- !is.na(site_idx)
        pred <- Q[site_idx[ok], sp, drop = FALSE]
        if (is.null(modern_pred_acc)) {
          modern_pred_acc <- matrix(0, nrow(site_comm), length(sp),
                                    dimnames = list(site_comm$site_id, sp))
        }
        modern_pred_acc[ok, sp] <- modern_pred_acc[ok, sp] + draw_weight[[dd]] * pred
        modern_pred_weight <- modern_pred_weight + draw_weight[[dd]]
      }
    }
    prev_Q <- Q
    prev_cells <- e_tm$cell_id
    prev_nodes <- as.integer(r_tm$node)
    prev_time <- tm
  }
  completed_draws <- c(completed_draws, dd)
  draw_tasks$status[draw_tasks$response_draw == dd] <- "complete"
  write_csv(draw_tasks, draw_task_path)
  write_draw_checkpoint()
  log_msg("Checkpointed completed response draw", dd)
}

terminal_anchor_calibration <- if (length(terminal_anchor_calibration_rows)) {
  do.call(rbind, terminal_anchor_calibration_rows)
} else {
  data.frame()
}
if (nrow(terminal_anchor_calibration)) {
  write_csv(
    terminal_anchor_calibration,
    file.path(dirs$modern, "terminal_anchor_prevalence_calibration.csv")
  )
}

root_initialisation_audit <- data.frame()
if (length(root_rows)) {
  root_initialisation_particles <- do.call(rbind, root_rows)
  write_csv(root_initialisation_particles,
            file.path(dirs$root, "root_initialisation_particles.csv"))

  # The root of a dated tree is one ancestral lineage.  In the default shared
  # patch mode, every immediately active child branch must receive exactly the
  # same spatial state; this audit makes that invariant reviewable per draw.
  root_groups <- split(
    root_initialisation_particles,
    paste(root_initialisation_particles$response_draw,
          root_initialisation_particles$time_ma, sep = "\r")
  )
  root_initialisation_audit <- do.call(rbind, lapply(root_groups, function(z) {
    cells_by_lineage <- split(as.character(z$cell_id), as.character(z$lineage))
    unique_cell_sets <- unique(vapply(
      cells_by_lineage,
      function(x) paste(sort(unique(x)), collapse = "|"),
      character(1)
    ))
    shared_patch <- length(unique_cell_sets) == 1L
    q_by_lineage <- split(z$q_root, as.character(z$lineage))
    unique_q_sets <- unique(vapply(
      q_by_lineage,
      function(x) paste(formatC(x[order(x)], digits = 15L, format = "fg"),
                         collapse = "|"),
      character(1)
    ))
    shared_q <- length(unique_q_sets) == 1L
    data.frame(
      response_draw = as.character(z$response_draw[1]),
      time_ma = z$time_ma[1],
      root_initialisation = z$root_prior[1],
      n_initial_child_lineages = length(cells_by_lineage),
      n_unique_root_cells = length(unique(z$cell_id)),
      n_positive_lineage_cell_rows = nrow(z),
      root_cells_shared_by_all_initial_children = shared_patch,
      root_q_shared_by_all_initial_children = shared_q,
      status = if (identical(root_initialisation,
                             "compact_shared_root_patch")) {
        if (shared_patch && shared_q) "PASS_SHARED_CONNECTED_ROOT_PATCH" else "FAIL_SHARED_ROOT_MISMATCH"
      } else "SENSITIVITY_MODE_NOT_SHARED_ROOT_PATCH",
      stringsAsFactors = FALSE
    )
  }))
  write_csv(root_initialisation_audit,
            file.path(dirs$root, "case04_root_initialisation_audit.csv"))
}
if (length(movement_summary)) write_csv(do.call(rbind, movement_summary),
                                        file.path(dirs$dispersal, "movement_matrix_summary.csv"))
if (length(transport_summary)) write_csv(do.call(rbind, transport_summary),
                                         file.path(dirs$transport, "plate_transport_summary.csv"))
if (length(speciation_summary)) write_csv(do.call(rbind, speciation_summary),
                                          file.path(dirs$speciation, "speciation_event_summary.csv"))
attribution_metadata <- data.frame(
  field = c("write_attribution_state", "attribution_state_path",
            "rows_written", "max_rows", "max_cells_per_interval",
            "max_lineages_per_interval", "sampling_interpretation",
            "full_budget_scope", "full_budget_targets"),
  value = c(write_attribution_state, attribution_state_path,
            attribution_rows_written, attribution_max_rows,
            attribution_max_cells, attribution_max_lineages,
            if (attribution_rows_written >= attribution_max_rows) {
              "sampled_and_row_capped_state_budget_input"
            } else {
              "complete_or_limit_not_reached_for_selected_run"
            },
            "all_active_lineages_all_valid_land_cells_all_processed_intervals",
            paste(full_attr_targets, collapse = ",")),
  stringsAsFactors = FALSE
)
write_csv(attribution_metadata, attribution_meta_path)

full_attr_dir <- safe_dir_create(file.path(dirs$sensitivity,
                                           "full_process_contribution"))
full_attr_time <- if (length(full_attr_time_series)) {
  do.call(rbind, full_attr_time_series)
} else {
  data.frame()
}
write_csv(full_attr_time, file.path(full_attr_dir,
                                    "case05_full_process_contribution_time_series.csv"))

state_summary_from_components <- function(tab) {
  if (!nrow(tab)) return(data.frame())
  rows <- list()
  for (tg in unique(tab$target)) {
    d <- tab[tab$target == tg, , drop = FALSE]
    agg <- stats::aggregate(d$value * d$interval_weight,
                            d[, "component", drop = FALSE], sum,
                            na.rm = TRUE)
    names(agg)[2] <- "weighted_value"
    denom <- sum(d$interval_weight[!duplicated(
      paste(d$posterior_draw_id, d$time_start_ma, d$time_end_ma, sep = "\r")
    )], na.rm = TRUE)
    if (!is.finite(denom) || denom <= 0) denom <- 1
    agg$estimate <- agg$weighted_value / denom
    get <- function(nm) {
      z <- agg$estimate[agg$component == nm]
      if (length(z)) z[1] else NA_real_
    }
    out <- data.frame(
      target = tg,
      process = c("Persistence", "Colonisation", "Colonisation",
                  "Persistence"),
      process_slug = c("persistence", "colonisation", "colonisation",
                       "persistence"),
      layer = "state_transition",
      contribution_type = c("retained_occupancy", "new_colonisation",
                            "recolonisation", "local_loss"),
      estimate = c(get("persistent_occupancy"), get("new_colonisation"),
                   get("recolonisation"), get("local_loss")),
      signed_contribution = c(get("persistent_occupancy"),
                              get("new_colonisation"),
                              get("recolonisation"), -get("local_loss")),
      status = "estimated_exact_ctmc_budget_full_streaming",
      taxon_scope = "all_active_lineages",
      spatial_scope = "all_valid_land_cells",
      normalise = "time",
      stringsAsFactors = FALSE
    )
    out$absolute_contribution <- abs(out$signed_contribution)
    src <- out$contribution_type %in%
      c("retained_occupancy", "new_colonisation", "recolonisation")
    total <- sum(out$estimate[src], na.rm = TRUE)
    out$share_within_layer <- NA_real_
    if (is.finite(total) && total > 0) {
      out$share_within_layer[src] <- out$estimate[src] / total
    }
    rows[[tg]] <- out
  }
  do.call(rbind, rows)
}
full_attr_summary <- state_summary_from_components(full_attr_time)
write_csv(full_attr_summary, file.path(full_attr_dir,
                                       "case05_full_process_contribution_summary.csv"))

full_attr_draws <- if (nrow(full_attr_time)) {
  stats::aggregate(
    full_attr_time$value * full_attr_time$interval_weight,
    full_attr_time[, c("posterior_draw_id", "target", "component"),
                   drop = FALSE],
    sum,
    na.rm = TRUE
  )
} else {
  data.frame()
}
if (nrow(full_attr_draws)) {
  names(full_attr_draws)[ncol(full_attr_draws)] <- "weighted_integral"
}
write_csv(full_attr_draws, file.path(full_attr_dir,
                                     "case05_full_process_contribution_draws.csv"))

full_attr_map_out <- full_attr_map
for (comp in full_attr_components) {
  full_attr_map_out[[comp]] <- full_attr_map_out[[comp]] /
    pmax(full_attr_map_out$integration_weight, 1e-12)
}
full_attr_map_out$target <- "lineage_richness"
full_attr_map_out$taxon_scope <- "all_active_lineages"
full_attr_map_out$spatial_scope <- "all_valid_land_cells"
write_csv(full_attr_map_out, file.path(full_attr_dir,
                                       "case05_full_process_contribution_history_maps.csv"))

metric_names <- setdiff(names(accumulators[[1]]), c("cell", "w"))
all_time_tables <- list()
for (tm in inside_times) {
  acc <- accumulators[[as.character(tm)]]
  df <- acc$cell
  for (nm in metric_names) df[[nm]] <- acc[[nm]] / max(acc$w, 1e-12)
  df$analysis_mode <- paste0("global_no_bgb_", dispersal_scheme, "_", transport_mode)
  df$terminal_anchor_mode <- if (tm == 0) terminal_anchor_mode else "none"
  df$terminal_anchor_weight <- if (tm == 0 &&
                                     !identical(terminal_anchor_mode, "none")) {
    acc$terminal_anchor_weight / max(acc$w, 1e-12)
  } else 0
  df$reporting_interval_myr <- reporting_interval_myr
  df$persistence_reference_myr <- persistence_reference_myr
  df$interpretation <- if (tm == 0 && !identical(terminal_anchor_mode, "none")) {
    paste(
      "sampled_surviving_lineage_forward_scenario_not_total_flora",
      "0Ma_terminal_model_ensemble_not_backward_smoothing_or_independent_validation;",
      "raw-forward richness and anchor difference are retained;",
      "explicit P2 probabilities use reporting_interval_myr"
    )
  } else {
    paste(
      "sampled_surviving_lineage_forward_scenario_not_total_flora",
      "not_smoothed_posterior; explicit P2 probabilities use reporting_interval_myr"
    )
  }
  all_time_tables[[as.character(tm)]] <- df
  metric_path_base <- file.path(dirs$inference,
                                paste0("cell_metrics_", slug_time(tm)))
  if (identical(metric_storage, "rds_per_time")) {
    saveRDS(df, paste0(metric_path_base, ".rds"), compress = "xz")
  } else {
    write_csv(df, paste0(metric_path_base, ".csv"))
  }
}
all_metrics <- do.call(rbind, all_time_tables)
if (identical(metric_storage, "csv")) {
  write_csv(all_metrics, file.path(dirs$inference, "cell_time_all_metrics.csv"))
}

# Join the independently selected 40-point audit panel back to the complete
# full-grid results.  This creates an inspectable process trace while leaving
# all calculations and maps on the entire valid palaeo-land grid.
diagnostic_metric_columns <- setdiff(
  names(all_metrics),
  c("cell_id", "time_ma", "lon", "lat", "H_state", "land_area_km2",
    "analysis_mode", "interpretation")
)
diagnostic_trace <- merge(
  diagnostic_points,
  all_metrics[, c("cell_id", "time_ma", diagnostic_metric_columns), drop = FALSE],
  by = c("cell_id", "time_ma"), all.x = TRUE, sort = FALSE
)
diagnostic_trace <- diagnostic_trace[order(-diagnostic_trace$time_ma,
                                           diagnostic_trace$diagnostic_rank), , drop = FALSE]
if (nrow(diagnostic_trace) != nrow(diagnostic_points) ||
    any(!is.finite(diagnostic_trace$mean_occupancy_probability))) {
  stop("A diagnostic point could not be joined to the complete full-grid process results.",
       call. = FALSE)
}
diagnostic_trace$trace_role <-
  "quality_assurance_points_only_not_a_model_or_map_subset"
write_csv(diagnostic_trace,
          file.path(dirs$inference, "case04_diagnostic_point_process_trace.csv"))
diagnostic_trace_coverage <- aggregate(
  cell_id ~ time_ma, diagnostic_trace,
  function(x) length(unique(x))
)
names(diagnostic_trace_coverage)[2] <- "n_diagnostic_points_joined"
diagnostic_trace_coverage <- merge(
  diagnostic_coverage,
  diagnostic_trace_coverage,
  by = "time_ma", all.x = TRUE, sort = FALSE
)
diagnostic_trace_coverage$trace_join_complete <-
  diagnostic_trace_coverage$n_diagnostic_points_joined ==
  diagnostic_trace_coverage$expected_diagnostic_points
if (!all(diagnostic_trace_coverage$trace_join_complete)) {
  stop("Diagnostic point trace coverage is incomplete.", call. = FALSE)
}
write_csv(diagnostic_trace_coverage,
          file.path(dirs$inference, "case04_diagnostic_point_trace_coverage.csv"))

ref <- all_metrics[, c("cell_id", "lon", "lat"), drop = FALSE]
ref <- ref[!duplicated(ref$cell_id), , drop = FALSE]
trajectory_richness_column <- if (identical(terminal_anchor_mode, "none")) {
  "expected_lineage_richness"
} else {
  "forward_raw_expected_lineage_richness"
}
refugia <- stats::aggregate(
  stats::as.formula(paste(trajectory_richness_column, "~ cell_id")), all_metrics,
                            mean, na.rm = TRUE)
names(refugia)[2] <- "richness_refugia_score"
env_refugia <- stats::aggregate(potential_environmental_richness ~ cell_id,
                                all_metrics, mean, na.rm = TRUE)
names(env_refugia)[2] <- "environmental_refugia_score"
refugia <- merge(refugia, env_refugia, by = "cell_id", all = TRUE)
refugia <- merge(refugia, ref, by = "cell_id", all.x = TRUE)
write_csv(refugia, file.path(dirs$diversity, "refugia", "refugia_scores.csv"))

gain_loss <- list()
for (i in 2:length(inside_times)) {
  old <- inside_times[i - 1L]
  young <- inside_times[i]
  a <- all_time_tables[[as.character(old)]][, c("cell_id", trajectory_richness_column)]
  b <- all_time_tables[[as.character(young)]][, c("cell_id", "lon", "lat",
                                                  trajectory_richness_column)]
  names(a)[2] <- "rich_old"
  names(b)[4] <- "rich_young"
  z <- merge(b, a, by = "cell_id", all.x = TRUE, sort = FALSE)
  z$rich_old[is.na(z$rich_old)] <- 0
  z$time_ma <- young
  z$previous_time_ma <- old
  z$net_lineage_richness_change <- z$rich_young - z$rich_old
  z$absolute_lineage_richness_change <- abs(z$net_lineage_richness_change)
  gain_loss[[length(gain_loss) + 1L]] <- z
}
if (length(gain_loss)) {
  gl <- do.call(rbind, gain_loss)
  write_csv(gl, file.path(dirs$diversity, "gain_loss", "gain_loss_turnover.csv"))
}

metric_limits <- list(
  expected_lineage_richness = c(0, length(unique(tip_beta$species))),
  forward_raw_expected_lineage_richness = c(0, length(unique(tip_beta$species))),
  terminal_anchor_difference_expected_lineage_richness = c(
    min(all_metrics$terminal_anchor_difference_expected_lineage_richness, na.rm = TRUE),
    max(all_metrics$terminal_anchor_difference_expected_lineage_richness, na.rm = TRUE)
  ),
  terminal_anchor_weight = c(0, 1),
  hmsc_fixed_effect_nowcast_richness = c(0, length(unique(tip_beta$species))),
  potential_environmental_richness = c(0, length(unique(tip_beta$species))),
  relative_environmental_support_sum = c(0, length(unique(tip_beta$species))),
  binary_lineage_richness = c(0, length(unique(tip_beta$species))),
  mean_occupancy_probability = c(0, 1),
  mean_environmental_support = c(0, 1),
  mean_arrival_hazard_per_myr = c(0, max(all_metrics$mean_arrival_hazard_per_myr, na.rm = TRUE)),
  mean_arrival_probability_reporting_interval = c(0, 1),
  mean_colonisation_hazard_per_myr = c(0, max(all_metrics$mean_colonisation_hazard_per_myr, na.rm = TRUE)),
  mean_colonisation_probability_reporting_interval = c(0, 1),
  mean_arrival_probability = c(0, 1),
  mean_colonisation_probability = c(0, 1),
  mean_persistence_probability = c(0, 1),
  mean_local_extinction_probability = c(0, 1),
  biotic_competition_pressure = c(0, 1),
  biotic_facilitation_pressure = c(0, 1),
  speciation_inheritance_footprint = c(0, 1),
  scenario_lineage_extinction_pressure = c(0, 1),
  simpson_diversity = c(0, 1),
  shannon_diversity = c(0, log(max(2, length(unique(tip_beta$species))))),
  occupancy_weighted_effective_lineages_shannon = c(0, length(unique(tip_beta$species))),
  occupancy_weighted_effective_lineages_simpson = c(0, length(unique(tip_beta$species))),
  occupancy_weighted_diversity_defined = c(0, 1),
  weighted_endemism = c(0, max(all_metrics$weighted_endemism, na.rm = TRUE)),
  niche_response_dispersion = c(0, max(all_metrics$niche_response_dispersion, na.rm = TRUE))
)
metric_dirs <- list(
  expected_lineage_richness = "expected_lineage_richness",
  forward_raw_expected_lineage_richness = "0Ma_forward_raw_expected_lineage_richness",
  terminal_anchor_difference_expected_lineage_richness = "0Ma_terminal_anchor_difference_expected_lineage_richness",
  terminal_anchor_weight = "0Ma_terminal_anchor_weight",
  hmsc_fixed_effect_nowcast_richness = "0Ma_hmsc_fixed_effect_nowcast_richness",
  potential_environmental_richness = "legacy_relative_environmental_support_sum",
  relative_environmental_support_sum = "relative_environmental_support_sum",
  binary_lineage_richness = "binary_lineage_richness",
  mean_occupancy_probability = "mean_occupancy_probability",
  mean_environmental_support = "mean_environmental_support",
  mean_arrival_hazard_per_myr = "arrival_hazard_per_myr",
  mean_arrival_probability_reporting_interval = "arrival_probability_reporting_interval",
  mean_colonisation_hazard_per_myr = "colonisation_hazard_per_myr",
  mean_colonisation_probability_reporting_interval = "colonisation_probability_reporting_interval",
  mean_arrival_probability = "arrival_probability",
  mean_colonisation_probability = "colonisation_probability",
  mean_persistence_probability = "persistence_probability",
  mean_local_extinction_probability = "local_extinction_probability",
  biotic_competition_pressure = "biotic_filtering",
  biotic_facilitation_pressure = "biotic_filtering",
  speciation_inheritance_footprint = "speciation",
  scenario_lineage_extinction_pressure = "extinction",
  weighted_endemism = "weighted_endemism",
  shannon_diversity = "occupancy_weighted_shannon_entropy",
  simpson_diversity = "occupancy_weighted_gini_simpson",
  occupancy_weighted_effective_lineages_shannon = "occupancy_weighted_effective_lineages_shannon",
  occupancy_weighted_effective_lineages_simpson = "occupancy_weighted_effective_lineages_simpson",
  occupancy_weighted_diversity_defined = "occupancy_weighted_diversity_defined",
  niche_response_dispersion = "niche_response_dispersion"
)
# Keep full metric names in tables and indexes, but make the physical map
# directories short enough for a nested Windows output root plus PNG/TIFF
# filename. These codes are stable API-facing file identifiers, not aliases
# for the scientific metric names.
metric_map_code <- c(
  expected_lineage_richness = "linrich",
  forward_raw_expected_lineage_richness = "forwardraw",
  terminal_anchor_difference_expected_lineage_richness = "anchordelta",
  terminal_anchor_weight = "anchorweight",
  hmsc_fixed_effect_nowcast_richness = "hmscnowcast",
  potential_environmental_richness = "envsupport",
  relative_environmental_support_sum = "envsupportsum",
  binary_lineage_richness = "binrich",
  mean_occupancy_probability = "occupancy",
  mean_environmental_support = "envsupportmean",
  mean_arrival_hazard_per_myr = "arrivalhaz",
  mean_arrival_probability_reporting_interval = "arrivalprob",
  mean_colonisation_hazard_per_myr = "colonisehaz",
  mean_colonisation_probability_reporting_interval = "coloniseprob",
  mean_arrival_probability = "arrivalprobold",
  mean_colonisation_probability = "coloniseprobold",
  mean_persistence_probability = "persist",
  mean_local_extinction_probability = "localext",
  biotic_competition_pressure = "bioticcomp",
  biotic_facilitation_pressure = "bioticfac",
  speciation_inheritance_footprint = "speciation",
  scenario_lineage_extinction_pressure = "extpressure",
  weighted_endemism = "endemism",
  shannon_diversity = "owshannon",
  simpson_diversity = "owgini",
  occupancy_weighted_effective_lineages_shannon = "oweffshan",
  occupancy_weighted_effective_lineages_simpson = "oweffsimp",
  occupancy_weighted_diversity_defined = "owdefined",
  niche_response_dispersion = "nrd"
)
metric_dirs <- stats::setNames(
  paste0("m_", unname(metric_map_code[names(metric_dirs)])),
  names(metric_dirs)
)
# The legacy potential_environmental_richness column is retained in the table
# for backwards compatibility, but the reported map is the explicitly named
# relative_environmental_support_sum.  It is intercept-free and must not be
# compared as an occurrence-richness estimate.
metric_dirs <- metric_dirs[names(metric_dirs) != "potential_environmental_richness"]
if (identical(terminal_anchor_mode, "none")) {
  metric_dirs <- metric_dirs[!names(metric_dirs) %in% c(
    "forward_raw_expected_lineage_richness",
    "terminal_anchor_difference_expected_lineage_richness",
    "terminal_anchor_weight"
  )]
}
if (!run_biotic_filtering) {
  metric_dirs <- metric_dirs[!names(metric_dirs) %in%
                             c("biotic_competition_pressure", "biotic_facilitation_pressure")]
}
if (!run_speciation_demo) {
  metric_dirs <- metric_dirs[names(metric_dirs) != "speciation_inheritance_footprint"]
}
if (!run_lineage_extinction_demo) {
  metric_dirs <- metric_dirs[names(metric_dirs) != "scenario_lineage_extinction_pressure"]
}
for (md in unique(unlist(metric_dirs))) safe_dir_create(file.path(dirs$diversity, md))

for (tm in map_times) {
  df <- all_time_tables[[as.character(tm)]]
  for (metric in names(metric_dirs)) {
    if (metric %in% c(
      "hmsc_fixed_effect_nowcast_richness",
      "forward_raw_expected_lineage_richness",
      "terminal_anchor_difference_expected_lineage_richness",
      "terminal_anchor_weight"
    ) && tm != 0) {
      next
    }
    idx <- write_map(df, metric, file.path(dirs$diversity, metric_dirs[[metric]]),
                     metric_limits[[metric]], metric)
    map_index[[length(map_index) + 1L]] <- idx
  }
}
ref_idx1 <- write_map(refugia, "richness_refugia_score",
                      file.path(dirs$diversity, "refugia"),
                      c(0, max(refugia$richness_refugia_score, na.rm = TRUE)),
                      "richness refugia")
ref_idx2 <- write_map(refugia, "environmental_refugia_score",
                      file.path(dirs$diversity, "refugia"),
                      c(0, max(refugia$environmental_refugia_score, na.rm = TRUE)),
                      "environment refugia")
map_index <- c(map_index, list(ref_idx1, ref_idx2))
map_index_df <- do.call(rbind, map_index)
write_csv(map_index_df, file.path(dirs$report, "case04_final_map_index.csv"))

time_summary <- aggregate(
  all_metrics[, names(metric_dirs), drop = FALSE],
  all_metrics[, "time_ma", drop = FALSE],
  function(v) mean(v, na.rm = TRUE)
)
names(time_summary)[-1] <- paste0("global_mean_", names(time_summary)[-1])
write_csv(time_summary, file.path(dirs$report, "case04_final_time_summary.csv"))
write_csv(
  data.frame(
    terminal_anchor_mode = terminal_anchor_mode,
    terminal_anchor_weight_requested = terminal_anchor_weight,
    terminal_anchor_prevalence_tolerance = paste(
      terminal_anchor_prevalence_tolerance, collapse = ","
    ),
    terminal_anchor_weight_realised_mean = if (nrow(terminal_anchor_calibration)) {
      mean(terminal_anchor_calibration$selected_anchor_weight)
    } else if (identical(terminal_anchor_mode, "hmsc_fixed_effect_convex")) {
      terminal_anchor_weight
    } else 0,
    terminal_anchor_weight_realised_min = if (nrow(terminal_anchor_calibration)) {
      min(terminal_anchor_calibration$selected_anchor_weight)
    } else if (identical(terminal_anchor_mode, "hmsc_fixed_effect_convex")) {
      terminal_anchor_weight
    } else 0,
    terminal_anchor_weight_realised_max = if (nrow(terminal_anchor_calibration)) {
      max(terminal_anchor_calibration$selected_anchor_weight)
    } else if (identical(terminal_anchor_mode, "hmsc_fixed_effect_convex")) {
      terminal_anchor_weight
    } else 0,
    endpoint_role = if (identical(terminal_anchor_mode, "none")) {
      "not_applied"
    } else {
      paste(
        "0Ma-only endpoint model ensemble; it is not backward smoothing, an ancient",
        "process, or independent validation. Pre-0Ma forward states remain unedited."
      )
    },
    trajectory_products_use = trajectory_richness_column,
    stringsAsFactors = FALSE
  ),
  file.path(dirs$report, "case04_terminal_anchor_interpretation.csv")
)

if (!is.null(modern_pred_acc) && !is.null(site_comm)) {
  sp <- colnames(modern_pred_acc)
  obs <- as.matrix(site_comm[, sp, drop = FALSE])
  storage.mode(obs) <- "double"
  pred <- modern_pred_acc[, sp, drop = FALSE] / max(modern_pred_weight, 1e-12)
  common <- is.finite(obs) & is.finite(pred)
  val <- data.frame(
    status = "DIAGNOSTIC_ONLY_SAME_DATA_FAMILY",
    n_cells = nrow(pred),
    n_species = ncol(pred),
    brier_score = mean((obs[common] - pred[common])^2, na.rm = TRUE),
    mean_observed = mean(obs[common], na.rm = TRUE),
    mean_predicted = mean(pred[common], na.rm = TRUE),
    interpretation = "Uses modern GBIF-derived grid also used by HMSC posterior; not independent endpoint conditioning.",
    stringsAsFactors = FALSE
  )
  write_csv(val, file.path(dirs$modern, "modern_endpoint_diagnostic_summary.csv"))
} else {
  write_csv(data.frame(status = "NOT_RUN", detail = "No matching comm/site endpoint data."),
            file.path(dirs$modern, "modern_endpoint_diagnostic_summary.csv"))
}

if (file.exists(paths$fossils)) {
  fossils <- utils::read.csv(paths$fossils, check.names = FALSE, stringsAsFactors = FALSE)
  fsum <- data.frame(
    status = "DIAGNOSTIC_LAYER_NOT_WEIGHTED",
    n_fossil_records = nrow(fossils),
    n_family_records = if ("fossil_family" %in% names(fossils)) {
      sum(!is.na(fossils$fossil_family) & nzchar(fossils$fossil_family))
    } else NA_integer_,
    interpretation = "Family-level fossils copied for audit; this final forward run does not yet reweight histories by fossil likelihood.",
    stringsAsFactors = FALSE
  )
  write_csv(fsum, file.path(dirs$fossil, "fossil_constraint_diagnostic_summary.csv"))
  write_csv(utils::head(fossils, 10000), file.path(dirs$fossil, "fossil_records_preview.csv"))
}

root_initialisation_status <- if (!nrow(root_initialisation_audit)) {
  "FAIL_ROOT_INITIALISATION_AUDIT_MISSING"
} else if (identical(root_initialisation, "compact_shared_root_patch")) {
  if (all(root_initialisation_audit$status == "PASS_SHARED_CONNECTED_ROOT_PATCH")) {
    "PASS_SHARED_CONNECTED_ROOT_PATCH"
  } else "FAIL_SHARED_ROOT_MISMATCH"
} else "SENSITIVITY_MODE_NOT_SHARED_ROOT_PATCH"

validation <- data.frame(
  check = c("no_biogeobears", "no_discrete_regions", "no_A_hist",
            "tree_domain", "root_initialisation", "copy_then_diverge_inheritance",
            "global_grid_complete", "probability_range",
            "diagnostic_points_per_time", "diagnostic_trace_join",
            "same_metric_fixed_legends", "geotiff_outputs",
            "modern_endpoint_conditioning", "fossil_conditioning",
            "plate_transport", "biotic_filtering_process",
            "speciation_event_table", "lineage_extinction_process_boundary"),
  status = c("PASS", "PASS", "PASS", "PASS", root_initialisation_status, "PASS",
             "PASS",
             ifelse(all(all_metrics$mean_occupancy_probability >= 0 &
                          all_metrics$mean_occupancy_probability <= 1,
                        na.rm = TRUE), "PASS", "FAIL"),
             if (all(diagnostic_coverage$n_diagnostic_points ==
                     diagnostic_coverage$expected_diagnostic_points)) "PASS" else "FAIL",
             if (all(diagnostic_trace_coverage$trace_join_complete)) "PASS" else "FAIL",
             "PASS",
             ifelse(make_tif && requireNamespace("terra", quietly = TRUE),
                    "PASS", "SKIP"),
             if (identical(terminal_anchor_mode, "none")) {
               "DIAGNOSTIC_NOT_CONDITIONED"
             } else if (identical(terminal_anchor_mode,
                                  "hmsc_fixed_effect_prevalence_calibrated")) {
               "TERMINAL_PREVALENCE_CALIBRATED_TRAINING_DIAGNOSTIC_ONLY"
             } else {
               "TERMINAL_MODEL_ENSEMBLE_TRAINING_DIAGNOSTIC_ONLY"
             },
             "DIAGNOSTIC_NOT_WEIGHTED",
             if (identical(transport_mode, "paleomap_h3_targetcentred")) {
               "PASS_EXTERNAL_TRACK_APPROXIMATION_NOT_POLYGON_OVERLAP"
             } else "IDENTITY_DEMO_NOT_REAL_PLATE_TRANSPORT",
             if (run_biotic_filtering) "SCENARIO_DEMO_NOT_EMPIRICAL" else "OFF_IN_EMPIRICAL_CORE",
             if (length(speciation_summary)) "PASS" else "SKIP",
             if (run_lineage_extinction_demo) "SCENARIO_DEMO_NOT_EMPIRICAL" else "GLOBAL_EXTINCTION_NOT_ESTIMATED"),
  detail = c(
    "No BioGeoBEARS/BSM input is read by this final no-BGB script.",
    "region_id/region_cube is not used as a process constraint.",
    "No continuous historical accessibility multiplier is used.",
    paste0("Projected ", length(inside_times), " time slices <= tree root ", round(root_age, 2), " Ma; time selection = ", time_selection, "."),
    if (identical(root_initialisation, "compact_shared_root_patch")) {
      "At the oldest model time, all immediately active child branches inherit the same connected root patch; see 08_root_initialisation/case04_root_initialisation_audit.csv."
    } else "Root mode is a labelled sensitivity configuration, not the shared-root empirical-core prior.",
    "Daughter lineages copy the nearest previous active ancestor state across speciation intervals.",
    "All available H_state=1 cells in each projected time slice are used.",
    "Mean occupancy probability is in [0,1].",
    "Up to 40 seeded, spatially balanced valid land cells are selected per time slice for QA only; they do not subset the dynamic grid.",
    "Every diagnostic point joined to the full-grid environmental, dispersal, colonisation, persistence, occupancy and diversity output.",
    "Legend ranges are fixed per metric across all time slices.",
    "GeoTIFF maps are written when terra is available.",
    if (identical(terminal_anchor_mode, "none")) {
      "Modern endpoint is reported as diagnostic only because no independent held-out endpoint is supplied."
    } else {
      realised_anchor_weight <- if (nrow(terminal_anchor_calibration)) {
        formatC(mean(terminal_anchor_calibration$selected_anchor_weight),
                format = "f", digits = 6)
      } else {
        formatC(terminal_anchor_weight, format = "f", digits = 6)
      }
      paste0(
        "0 Ma is a convex terminal ensemble with HMSC fixed-effect weight ",
        realised_anchor_weight,
        if (identical(terminal_anchor_mode,
                      "hmsc_fixed_effect_prevalence_calibrated")) {
          paste0(" selected on the modern site grid to meet the prevalence-ratio tolerance [",
                 paste(formatC(terminal_anchor_prevalence_tolerance,
                               format = "f", digits = 3), collapse = ", "), "]")
        } else "",
        "; raw-forward endpoint is retained separately. This is not backward smoothing,",
        " a historical process, or independent validation because the current endpoint",
        " records belong to the HMSC training data family."
      )
    },
    "Family fossils are reported as diagnostics; no particle likelihood weighting yet.",
    if (identical(transport_mode, "paleomap_h3_targetcentred")) {
      paste0(
        "External PALEOMAP H3 target-centred carriage is used before movement; ",
        "it is a nearest-track approximation, not exact polygon overlap. Minimum target coverage required = ",
        transport_min_target_coverage, "."
      )
    } else "No plate transport table supplied; identity cell matching is used.",
    "Optional biotic filtering modifies establishment and persistence in all-processes demo; it is not inferred from independent interaction data.",
    "Speciation is tree-conditioned copy-then-diverge lineage identity inheritance, not an estimated speciation-rate surface.",
    "Extant-only tree cannot estimate total lineage extinction; the demo outputs a scenario pressure surface only."
  ),
  stringsAsFactors = FALSE
)
write_csv(validation, file.path(dirs$config, "case04_final_validation_checks.csv"))

report <- c(
  "# Case04 final no-BioGeoBEARS global dynamic results",
  "",
  "This output cancels BioGeoBEARS, BSM histories and discrete regional constraints.",
  "The model is a global 1-degree grid forward scenario for sampled-surviving Plant200 lineages.",
  "",
  "Core formula:",
  "",
  "`HMSC beta posterior -> ancestral beta -> S_env -> optional biotic filter -> movement -> arrival -> establishment -> colonisation -> persistence/local loss -> tree-conditioned speciation -> optional extinction pressure -> psi`.",
  "",
  "Important interpretation limits:",
  "",
  "- `psi` is not yet a fully endpoint-conditioned posterior history ensemble.",
  "- Modern endpoint and fossils are diagnostic layers in this run.",
  "- The tree root defines the valid biological time domain; older slices are not biological zeros.",
  "- Richness means sampled-surviving-lineage richness, not total historical plant diversity.",
  if (identical(transport_mode, "paleomap_h3_targetcentred")) {
    "- Plate carriage uses external PALEOMAP H3 target-centred tracks before active dispersal. It avoids identity-cell artefacts but remains a nearest-track approximation; polygon-overlap transport is still required for the strictest paper-grade claim."
  } else "- Plate transport is `identity_demo`; real plate-overlap transport is still needed for paper-grade palaeocoordinate carrying.",
  paste0("- Arrival and colonisation probability maps use an explicit ",
         reporting_interval_myr, " Myr reporting window; hazards remain per Myr."),
  "- Biotic filtering is empirical-core OFF unless `run_biotic_filtering=TRUE`; when enabled here it is a scenario modifier, not a causal interaction estimate.",
  "- Extinction outputs are local/geography-forced loss plus optional scenario pressure; total lineage extinction is not inferred from an extant-only tree.",
  "",
  paste0("Projected time slices: ", length(inside_times), " (", time_selection, "; ",
         length(available_inside_times), " available at or below the sampled-tree root)."),
  paste0("Species/tips used: ", length(unique(tip_beta$species))),
  paste0("Response draws used: ", length(draw_ids)),
  paste0("Biotic filtering enabled: ", run_biotic_filtering),
  paste0("Speciation demo maps enabled: ", run_speciation_demo),
  paste0("Lineage extinction demo maps enabled: ", run_lineage_extinction_demo),
  paste0("PNG maps: ", sum(file.exists(map_index_df$png), na.rm = TRUE)),
  paste0("GeoTIFF maps: ", sum(file.exists(map_index_df$tif), na.rm = TRUE))
)
write_text(report, file.path(dirs$report, "case04_final_interpretation.md"))

log_msg("Case04 final no-BGB global dynamic run completed.")
log_msg("Output:", output)
cat("\nDONE\n")
cat("Output:", output, "\n")
