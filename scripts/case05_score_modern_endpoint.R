#!/usr/bin/env Rscript

# Score the merged 0 Ma dynamic endpoint against the fixed-effect HMSC nowcast
# and an explicitly declared modern observation role. This script does not
# retrofit a historical trajectory: it provides the gate that a forward
# scenario must pass before its deep-time maps are called calibrated.

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

safe_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

read_csv <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    as.data.frame(data.table::fread(path, data.table = FALSE))
  } else utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

write_csv <- function(x, path) {
  safe_dir(dirname(path))
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(x, path)
  } else utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
  invisible(path)
}

cfg <- parse_args(args)
pkg_root <- normalizePath(
  cfg$pkg_root %||% "C:/Users/Google/Documents/HMSC-HIST/HmscEcoEvo",
  winslash = "/", mustWork = TRUE
)
scenario_output <- normalizePath(cfg$scenario_output %||% "",
                                 winslash = "/", mustWork = FALSE)
if (!dir.exists(scenario_output)) {
  stop("Use --scenario_output=<merged Case05 scenario directory>.", call. = FALSE)
}
default_input <- file.path(
  "C:/Users/Google/Documents/HMSC-HIST/HmscEcoEvo_API8_20260818",
  "data_external/plant_200_multifamily_extant_20260820/prepared_inputs",
  "case03_data_raw/Plant200_multifamily_extant_1deg_traits"
)
comm_path <- normalizePath(cfg$comm %||% file.path(default_input, "comm.csv"),
                           winslash = "/", mustWork = FALSE)
sites_path <- normalizePath(cfg$sites %||% file.path(default_input, "sites.csv"),
                            winslash = "/", mustWork = FALSE)
data_role <- cfg$data_role %||% "training_diagnostic"
if (!data_role %in% c("training_diagnostic", "spatial_holdout", "independent_validation")) {
  stop("data_role must be training_diagnostic, spatial_holdout or independent_validation.",
       call. = FALSE)
}
min_richness_correlation <- suppressWarnings(as.numeric(
  cfg$min_richness_correlation %||% 0.50
))
min_species_range_correlation <- suppressWarnings(as.numeric(
  cfg$min_species_range_correlation %||% 0.30
))
max_prevalence_ratio <- suppressWarnings(as.numeric(
  cfg$max_prevalence_ratio %||% 2
))
max_site_match_distance_deg <- suppressWarnings(as.numeric(
  cfg$max_site_match_distance_deg %||% 0.51
))
if (!is.finite(min_richness_correlation) ||
    !is.finite(min_species_range_correlation) ||
    min_species_range_correlation < -1 || min_species_range_correlation > 1 ||
    !is.finite(max_prevalence_ratio) ||
    max_prevalence_ratio < 1 || !is.finite(max_site_match_distance_deg) ||
    max_site_match_distance_deg <= 0) {
  stop("Endpoint gate thresholds are invalid.", call. = FALSE)
}

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(pkg_root, quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(HmscEcoEvo))
}
endpoint_path <- file.path(
  scenario_output, "20_modern_endpoint_validation",
  "modern_tip_endpoint_0Ma_posterior_mean.rds"
)
if (!file.exists(endpoint_path)) {
  stop("Missing merged endpoint payload: ", endpoint_path, call. = FALSE)
}
missing <- c(comm = comm_path, sites = sites_path)
missing <- missing[!file.exists(missing)]
if (length(missing)) {
  audit_dir <- safe_dir(file.path(scenario_output, "20_modern_endpoint_validation"))
  write_csv(data.frame(input = names(missing), path = unname(missing),
                       stringsAsFactors = FALSE),
            file.path(audit_dir, "missing_modern_endpoint_inputs.csv"))
  stop("Missing modern endpoint inputs. See missing_modern_endpoint_inputs.csv",
       call. = FALSE)
}

endpoint <- readRDS(endpoint_path)
comm <- read_csv(comm_path)
sites <- read_csv(sites_path)
if (!"site_id" %in% names(comm) || !"site_id" %in% names(sites) ||
    !all(c("lon", "lat") %in% names(sites))) {
  stop("comm.csv and sites.csv must contain site_id; sites.csv must also contain lon and lat.",
       call. = FALSE)
}
if (anyDuplicated(comm$site_id) || anyDuplicated(sites$site_id)) {
  stop("site_id must be unique in comm.csv and sites.csv.", call. = FALSE)
}
if (!setequal(comm$site_id, sites$site_id)) {
  stop("comm.csv and sites.csv have non-identical site_id sets.", call. = FALSE)
}
species <- as.character(endpoint$species)
missing_species <- setdiff(species, names(comm))
if (length(missing_species)) {
  stop("Modern community is missing endpoint species: ",
       paste(head(missing_species, 10L), collapse = ", "), call. = FALSE)
}
# `site_id` identifies the modern sampling table, whereas `cell_id` identifies
# the palaeoenvironment grid.  Match their spatial locations explicitly; do
# not rely on their unrelated naming conventions.
sites <- sites[match(comm$site_id, sites$site_id), , drop = FALSE]
coord_key <- function(lon, lat) {
  paste(formatC(as.numeric(lon), format = "f", digits = 6),
        formatC(as.numeric(lat), format = "f", digits = 6), sep = "_")
}
endpoint_key <- coord_key(endpoint$lon, endpoint$lat)
site_key <- coord_key(sites$lon, sites$lat)
site_index <- match(site_key, endpoint_key)

# A raster may have undergone harmless floating-point serialization.  Recover
# from that only by a nearest-grid lookup, then record the actual distance.
nearest_index <- function(lon, lat, grid_lon, grid_lat) {
  lon_delta <- abs(grid_lon - lon)
  lon_delta <- pmin(lon_delta, 360 - lon_delta)
  which.min(lon_delta^2 + (grid_lat - lat)^2)
}
fallback <- which(is.na(site_index))
if (length(fallback)) {
  site_index[fallback] <- vapply(fallback, function(i) {
    nearest_index(sites$lon[i], sites$lat[i], endpoint$lon, endpoint$lat)
  }, integer(1))
}
lon_delta <- abs(endpoint$lon[site_index] - sites$lon)
lon_delta <- pmin(lon_delta, 360 - lon_delta)
match_distance_deg <- sqrt(lon_delta^2 + (endpoint$lat[site_index] - sites$lat)^2)
site_match <- data.frame(
  site_id = comm$site_id,
  source_lon = sites$lon,
  source_lat = sites$lat,
  endpoint_cell_id = endpoint$cell_id[site_index],
  endpoint_lon = endpoint$lon[site_index],
  endpoint_lat = endpoint$lat[site_index],
  match_method = ifelse(site_key %in% endpoint_key, "exact_coordinate", "nearest_grid"),
  match_distance_deg = match_distance_deg,
  stringsAsFactors = FALSE
)
audit_dir <- safe_dir(file.path(scenario_output, "20_modern_endpoint_validation"))
write_csv(site_match, file.path(audit_dir, "modern_site_to_endpoint_grid_match.csv"))
if (any(!is.finite(match_distance_deg) | match_distance_deg > max_site_match_distance_deg)) {
  write_csv(site_match[!is.finite(match_distance_deg) |
                       match_distance_deg > max_site_match_distance_deg, , drop = FALSE],
            file.path(audit_dir, "modern_sites_exceed_endpoint_match_threshold.csv"))
  stop("Some modern sites exceed the endpoint-grid matching threshold. See modern_sites_exceed_endpoint_match_threshold.csv",
       call. = FALSE)
}
observed <- as.matrix(comm[, species, drop = FALSE])
storage.mode(observed) <- "double"
rownames(observed) <- comm$site_id
dynamic <- endpoint$dynamic_occupancy_probability[site_index, species, drop = FALSE]
nowcast <- endpoint$hmsc_fixed_effect_nowcast_probability[site_index, species, drop = FALSE]
has_forward_raw <- "forward_raw_occupancy_probability" %in% names(endpoint)
forward_raw <- if (has_forward_raw) {
  endpoint$forward_raw_occupancy_probability[site_index, species, drop = FALSE]
} else NULL
rownames(dynamic) <- rownames(nowcast) <- comm$site_id
colnames(dynamic) <- colnames(nowcast) <- species
if (has_forward_raw) {
  rownames(forward_raw) <- comm$site_id
  colnames(forward_raw) <- species
}

dynamic_score <- hee_modern_endpoint_score(
  dynamic, observed, data_role = data_role
)
nowcast_score <- hee_modern_endpoint_score(
  nowcast, observed, data_role = data_role
)
forward_raw_score <- if (has_forward_raw) {
  hee_modern_endpoint_score(forward_raw, observed, data_role = data_role)
} else NULL
matrix_comparison <- function(dynamic, nowcast) {
  dynamic <- as.matrix(dynamic)
  nowcast <- as.matrix(nowcast)
  dyn_richness <- rowSums(dynamic, na.rm = TRUE)
  now_richness <- rowSums(nowcast, na.rm = TRUE)
  dyn_range <- colSums(dynamic, na.rm = TRUE)
  now_range <- colSums(nowcast, na.rm = TRUE)
  safe_cor <- function(x, y, method = "pearson") {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3L || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) {
      return(NA_real_)
    }
    stats::cor(x[ok], y[ok], method = method)
  }
  data.frame(
    dynamic_prevalence = mean(dynamic),
    hmsc_nowcast_prevalence = mean(nowcast),
    dynamic_to_nowcast_prevalence_ratio = if (mean(nowcast) > 0) {
      mean(dynamic) / mean(nowcast)
    } else NA_real_,
    cell_richness_pearson = safe_cor(dyn_richness, now_richness),
    cell_richness_spearman = safe_cor(dyn_richness, now_richness, "spearman"),
    species_range_pearson = safe_cor(dyn_range, now_range),
    species_range_spearman = safe_cor(dyn_range, now_range, "spearman"),
    n_cells = nrow(dynamic),
    n_species = ncol(dynamic),
    stringsAsFactors = FALSE
  )
}
dynamic_vs_nowcast <- matrix_comparison(dynamic, nowcast)
forward_raw_vs_nowcast <- if (has_forward_raw) {
  matrix_comparison(forward_raw, nowcast)
} else NULL
score <- rbind(
  transform(dynamic_score$summary, prediction = "dynamic_endpoint"),
  if (has_forward_raw) {
    transform(forward_raw_score$summary, prediction = "forward_raw_endpoint")
  },
  transform(nowcast_score$summary, prediction = "hmsc_fixed_effect_nowcast")
)
score <- score[, c("prediction", setdiff(names(score), "prediction")), drop = FALSE]
dyn_prev <- score$predicted_prevalence[score$prediction == "dynamic_endpoint"]
now_prev <- score$predicted_prevalence[score$prediction == "hmsc_fixed_effect_nowcast"]
ratio <- if (is.finite(dyn_prev) && is.finite(now_prev) && now_prev > 0) dyn_prev / now_prev else NA_real_
rich_cor <- score$richness_correlation[score$prediction == "dynamic_endpoint"]
nowcast_observed_rich_cor <- score$richness_correlation[
  score$prediction == "hmsc_fixed_effect_nowcast"
]
gate <- data.frame(
  check = c("endpoint_matrix_complete", "dynamic_vs_nowcast_prevalence_ratio",
            "dynamic_vs_nowcast_cell_richness_correlation",
            "dynamic_vs_nowcast_species_range_correlation",
            "dynamic_endpoint_observed_richness_correlation",
            "hmsc_nowcast_observed_richness_correlation"),
  value = c(TRUE, ratio, dynamic_vs_nowcast$cell_richness_pearson,
            dynamic_vs_nowcast$species_range_pearson, rich_cor,
            nowcast_observed_rich_cor),
  threshold = c(NA_character_,
                paste0("[", 1 / max_prevalence_ratio, ", ", max_prevalence_ratio, "]"),
                paste0(">=", min_richness_correlation),
                paste0(">=", min_species_range_correlation),
                paste0(">=", min_richness_correlation),
                paste0(">=", min_richness_correlation)),
  status = c(
    "PASS",
    if (is.finite(ratio) && ratio >= 1 / max_prevalence_ratio && ratio <= max_prevalence_ratio) "PASS" else "FAIL",
    if (is.finite(dynamic_vs_nowcast$cell_richness_pearson) &&
          dynamic_vs_nowcast$cell_richness_pearson >= min_richness_correlation) "PASS" else "FAIL",
    if (is.finite(dynamic_vs_nowcast$species_range_pearson) &&
          dynamic_vs_nowcast$species_range_pearson >= min_species_range_correlation) "PASS" else "FAIL",
    if (data_role == "training_diagnostic") "DIAGNOSTIC_ONLY" else if (
      is.finite(rich_cor) && rich_cor >= min_richness_correlation
    ) "PASS" else "FAIL",
    if (data_role == "training_diagnostic") "DIAGNOSTIC_ONLY" else if (
      is.finite(nowcast_observed_rich_cor) &&
        nowcast_observed_rich_cor >= min_richness_correlation
    ) "PASS" else "FAIL"
  ),
  stringsAsFactors = FALSE
)
gate$scientific_role <- if (data_role == "training_diagnostic") {
  "diagnostic_only; do_not_call_independent_validation"
} else {
  "endpoint_calibration_or_validation"
}

out_dir <- safe_dir(file.path(scenario_output, "20_modern_endpoint_validation"))
write_csv(score, file.path(out_dir, "modern_endpoint_score_summary.csv"))
comparison <- rbind(
  transform(dynamic_vs_nowcast, prediction = "dynamic_endpoint"),
  if (has_forward_raw) {
    transform(forward_raw_vs_nowcast, prediction = "forward_raw_endpoint")
  }
)
comparison <- comparison[, c("prediction", setdiff(names(comparison), "prediction")),
                         drop = FALSE]
write_csv(comparison,
          file.path(out_dir, "dynamic_vs_hmsc_nowcast_endpoint_comparison.csv"))
write_csv(gate, file.path(out_dir, "modern_endpoint_calibration_gate.csv"))
write_csv(rbind(
  transform(dynamic_score$species, prediction = "dynamic_endpoint"),
  if (has_forward_raw) {
    transform(forward_raw_score$species, prediction = "forward_raw_endpoint")
  },
  transform(nowcast_score$species, prediction = "hmsc_fixed_effect_nowcast")
), file.path(out_dir, "modern_endpoint_species_scores.csv"))
write_csv(rbind(
  transform(dynamic_score$cells, prediction = "dynamic_endpoint"),
  if (has_forward_raw) {
    transform(forward_raw_score$cells, prediction = "forward_raw_endpoint")
  },
  transform(nowcast_score$cells, prediction = "hmsc_fixed_effect_nowcast")
), file.path(out_dir, "modern_endpoint_cell_scores.csv"))
write_csv(data.frame(
  endpoint_payload = normalizePath(endpoint_path, winslash = "/", mustWork = TRUE),
  comm = normalizePath(comm_path, winslash = "/", mustWork = TRUE),
  sites = normalizePath(sites_path, winslash = "/", mustWork = TRUE),
  data_role = data_role,
  n_sites = nrow(observed),
  n_species = ncol(observed),
  n_response_draws = endpoint$n_response_draws,
  terminal_anchor_mode = endpoint$terminal_anchor_mode %||% "none",
  terminal_anchor_weight = endpoint$terminal_anchor_weight %||% 0,
  forward_raw_endpoint_available = has_forward_raw,
  stringsAsFactors = FALSE
), file.path(out_dir, "modern_endpoint_score_inputs.csv"))

cat("Case05 modern endpoint score complete: ", out_dir, "\n", sep = "")
