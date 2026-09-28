#!/usr/bin/env Rscript

# Score one completed Case05 v3 forward-process draw at 0 Ma. Plant200 rows
# are target-occurrence cells, so their zeros are non-detections for the other
# taxa, not confirmed absences. This script therefore writes a presence-only
# spatial-support diagnostic by default. A Bernoulli prevalence/Brier score is
# produced only when an explicit surveyed-presence-absence design is supplied.

`%||%` <- function(x, y) if (is.null(x) || !length(x) || identical(x, "") || is.na(x)) y else x
parse_args <- function(args) {
  out <- list()
  for (arg in args) if (grepl("^--", arg)) {
    z <- sub("^--", "", arg); at <- regexpr("=", z, fixed = TRUE)
    if (at < 0L) out[[z]] <- TRUE else out[[substr(z, 1L, at - 1L)]] <- substr(z, at + 1L, nchar(z))
  }
  out
}
safe_dir <- function(x) { if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE); normalizePath(x, winslash = "/", mustWork = FALSE) }
write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("data.table", quietly = TRUE)) data.table::fwrite(x, path) else utils::write.csv(x, path, row.names = FALSE)
  invisible(path)
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
file_arg <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
if (!length(script)) script <- "scripts/case05_v3_score_training_diagnostic.R"
pkg_root <- normalizePath(cfg$pkg_root %||% file.path(dirname(normalizePath(script[[1L]], winslash = "/", mustWork = FALSE)), ".."), winslash = "/", mustWork = TRUE)
workspace <- normalizePath(file.path(pkg_root, ".."), winslash = "/", mustWork = TRUE)
scenario_output <- normalizePath(cfg$scenario_output %||% stop("--scenario_output is required.", call. = FALSE), winslash = "/", mustWork = TRUE)
carrier_root <- normalizePath(cfg$carrier_root %||% file.path(pkg_root, "derived_inputs", "case05_v2_paleomap_h3_carrier_325_0Ma_1deg"), winslash = "/", mustWork = TRUE)
earth_root <- normalizePath(cfg$earth_root %||% file.path(pkg_root, "derived_inputs", "case04_plant200_palaeo_1deg_slices_20260919"), winslash = "/", mustWork = TRUE)
plant_input_dir <- normalizePath(cfg$plant_input_dir %||% file.path(workspace, "HmscEcoEvo_API8_20260818", "data_external", "plant_200_multifamily_extant_20260820", "prepared_inputs", "case03_data_raw", "Plant200_multifamily_extant_1deg_traits"), winslash = "/", mustWork = TRUE)
max_output_distance_km <- suppressWarnings(as.numeric(cfg$max_output_distance_km %||% 500)[1L])
interpolation_k <- suppressWarnings(as.integer(cfg$interpolation_k %||% 4L)[1L])
if (!is.finite(max_output_distance_km) || max_output_distance_km <= 0 || !is.finite(interpolation_k) || interpolation_k < 1L) {
  stop("max_output_distance_km and interpolation_k must be positive.", call. = FALSE)
}
if (!requireNamespace("pkgload", quietly = TRUE)) stop("pkgload is required.", call. = FALSE)
pkgload::load_all(pkg_root, quiet = TRUE)

carrier_index <- utils::read.csv(file.path(carrier_root, "case05_v2_plate_carrier_index.csv"), stringsAsFactors = FALSE)
carrier_index$time_ma <- as.numeric(carrier_index$time_ma)
carrier_path <- carrier_index$file[which.min(abs(carrier_index$time_ma - 0))]
earth_index <- utils::read.csv(file.path(earth_root, "case04_palaeo_earth_slice_index.csv"), stringsAsFactors = FALSE)
earth_index$time_ma <- as.numeric(earth_index$time_ma)
earth_path <- earth_index$file[which.min(abs(earth_index$time_ma - 0))]
state_path <- file.path(scenario_output, "16_state_space_inference", "carrier_states", "occupancy_0Ma.rds")
if (!file.exists(state_path)) stop("Missing completed 0 Ma state: ", state_path, call. = FALSE)

carrier_blob <- readRDS(carrier_path)
earth <- as.data.frame(readRDS(earth_path), stringsAsFactors = FALSE)
state <- readRDS(state_path)
occupancy <- as.matrix(state$occupancy)
if (any(!is.finite(occupancy)) || any(occupancy < 0 | occupancy > 1)) stop("0 Ma occupancy must be finite probabilities in [0, 1].", call. = FALSE)
operator <- HmscEcoEvo::hee_plate_carrier_output_weights(
  carrier_grid = carrier_blob$carriers,
  target_grid = earth,
  k = interpolation_k,
  max_output_distance_km = max_output_distance_km
)
grid_probability <- HmscEcoEvo::hee_plate_carrier_apply_output_weights(occupancy, operator)

comm <- utils::read.csv(file.path(plant_input_dir, "comm.csv"), check.names = FALSE, stringsAsFactors = FALSE)
sites <- utils::read.csv(file.path(plant_input_dir, "sites.csv"), stringsAsFactors = FALSE)
if (!all(c("site_id", "lon", "lat") %in% names(sites)) || !"site_id" %in% names(comm)) {
  stop("comm.csv and sites.csv require site_id; sites.csv also requires lon and lat.", call. = FALSE)
}
sites <- sites[match(comm$site_id, sites$site_id), , drop = FALSE]
if (anyNA(sites$site_id)) stop("comm site_id values must match sites.csv exactly.", call. = FALSE)
species <- setdiff(names(comm), "site_id")
if (!setequal(species, colnames(grid_probability))) stop("0 Ma state species do not match comm.csv.", call. = FALSE)
grid_probability <- grid_probability[, species, drop = FALSE]
coord_key <- function(lon, lat) paste(formatC(as.numeric(lon), format = "f", digits = 6), formatC(as.numeric(lat), format = "f", digits = 6), sep = "_")
site_idx <- match(coord_key(sites$lon, sites$lat), coord_key(earth$lon, earth$lat))
if (anyNA(site_idx)) stop("Some modern sites are absent from the 0 Ma Earth grid.", call. = FALSE)
observed <- as.matrix(comm[, species, drop = FALSE]); storage.mode(observed) <- "double"
rownames(observed) <- comm$site_id
covered <- rowSums(is.finite(grid_probability[site_idx, , drop = FALSE])) == length(species)
if (!any(covered)) stop("No complete 0 Ma carrier-output predictions cover training sites.", call. = FALSE)
observation_audit <- HmscEcoEvo::hee_audit_observation_design(
  observed,
  design_label = "Plant200 GBIF target-occurrence cells; a row exists only when at least one selected taxon was recorded"
)
# Place known occurrence cells onto the full 0 Ma Earth grid. Unknown cells
# remain zero *records*, never biological absences, and the score below is
# conditional on the observed occurrence locations.
full_occurrence <- matrix(0, nrow = nrow(earth), ncol = length(species),
                          dimnames = list(as.character(earth$cell_id), species))
full_occurrence[site_idx[covered], ] <- observed[covered, , drop = FALSE]
presence_score <- HmscEcoEvo::hee_endpoint_presence_only_score(
  grid_probability[, species, drop = FALSE], full_occurrence,
  data_role = "training_diagnostic"
)
binary_score <- NULL
if (isTRUE(observation_audit$bernoulli_absence_likelihood_eligible)) {
  probability <- grid_probability[site_idx[covered], , drop = FALSE]
  rownames(probability) <- comm$site_id[covered]
  observed_binary <- observed[covered, , drop = FALSE]
  binary_score <- HmscEcoEvo::hee_modern_endpoint_score(
    probability, observed_binary, data_role = "training_diagnostic"
  )
}
coverage <- data.frame(
  n_training_sites = nrow(comm), n_complete_covered_sites = sum(covered),
  covered_site_fraction = mean(covered), interpolation_k = interpolation_k,
  max_output_distance_km = max_output_distance_km,
  carrier_output_coverage_fraction = operator$diagnostics$output_coverage_fraction,
  stringsAsFactors = FALSE
)
out <- safe_dir(file.path(scenario_output, "20_modern_training_diagnostic"))
write_csv(coverage, file.path(out, "carrier_to_training_site_coverage.csv"))
write_csv(observation_audit, file.path(out, "observation_design_audit.csv"))
write_csv(presence_score$summary,
          file.path(out, "presence_only_training_diagnostic_summary.csv"))
write_csv(presence_score$taxa,
          file.path(out, "presence_only_training_diagnostic_by_species.csv"))
if (is.null(binary_score)) {
  write_csv(data.frame(
    status = "NOT_RUN_UNCONFIRMED_ABSENCES",
    reason = "Plant200 rows are GBIF target-occurrence cells; zeros are non-detections, not surveyed absences.",
    stringsAsFactors = FALSE
  ), file.path(out, "bernoulli_endpoint_score_status.csv"))
} else {
  write_csv(binary_score$summary, file.path(out, "training_diagnostic_summary.csv"))
  write_csv(binary_score$species, file.path(out, "training_diagnostic_by_species.csv"))
  write_csv(binary_score$cells, file.path(out, "training_diagnostic_by_site.csv"))
}
writeLines(c(
  "# Case05 v3 modern training diagnostic",
  "",
  "Plant200 rows are GBIF target-occurrence cells. A zero for one selected taxon in a row where another taxon was observed is a non-detection, not a confirmed absence.",
  "presence_only_training_diagnostic_summary.csv therefore scores whether known occurrences fall in predicted spatial support relative to a uniform available-cell reference. It neither estimates prevalence nor corrects GBIF effort bias.",
  "The same records informed the HMSC posterior and its response reference, so this remains a training diagnostic, not an independent endpoint likelihood or particle weight.",
  "Carrier-to-grid interpolation is an output operation. Uncovered locations are excluded and counted in carrier_to_training_site_coverage.csv; they are never recoded as biological absences."
), file.path(out, "README.md"))
message("Case05 v3 training diagnostic complete: ", out)
