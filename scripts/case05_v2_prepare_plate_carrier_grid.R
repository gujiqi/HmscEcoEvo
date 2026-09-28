#!/usr/bin/env Rscript

# Prepare the plate-carrier state space used by Case05 v2.
#
# Stable PALEOMAP H3 identities carry occupancy through plate motion. Complete
# 1 degree palaeo-land grids provide environment/habitat at each time and are
# used again only for map rendering. This intentionally replaces the legacy
# target-centred 1 degree to 1 degree carriage approximation.

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || identical(x, "") || is.na(x)) y else x
}

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!grepl("^--", arg)) next
    x <- sub("^--", "", arg)
    at <- regexpr("=", x, fixed = TRUE)
    if (at < 0L) out[[x]] <- TRUE else {
      out[[substr(x, 1L, at - 1L)]] <- substr(x, at + 1L, nchar(x))
    }
  }
  out
}

as_num <- function(x, default) {
  z <- suppressWarnings(as.numeric(x)[1L])
  if (is.finite(z)) z else default
}

safe_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

write_csv <- function(x, path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(x, path)
  } else {
    utils::write.csv(x, path, row.names = FALSE)
  }
  invisible(path)
}

ma_slug <- function(x) sub("\\.$", "", formatC(x, format = "f", digits = 6,
                                                   drop0trailing = TRUE))

args <- parse_args(commandArgs(trailingOnly = TRUE))
script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
if (!length(script_file)) script_file <- "scripts/case05_v2_prepare_plate_carrier_grid.R"
script_dir <- dirname(normalizePath(script_file[[1L]], winslash = "/", mustWork = FALSE))
pkg_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
workspace <- normalizePath(file.path(pkg_root, ".."), winslash = "/", mustWork = TRUE)

slice_index <- args$palaeo_slice_index %||% file.path(
  pkg_root, "derived_inputs", "case04_plant200_palaeo_1deg_slices_20260919",
  "case04_palaeo_earth_slice_index.csv"
)
tracks_path <- args$plate_tracks %||% file.path(
  workspace, "data", "raw", "plate_reconstruction", "palaeoverse_gridded_recon_zenodo",
  "PALEOMAP.RDS"
)
output <- safe_dir(args$output %||% file.path(
  pkg_root, "derived_inputs", "case05_v2_paleomap_h3_carrier_325_0Ma_1deg"
))
max_match_distance_km <- as_num(args$max_match_distance_km, 250)
max_output_distance_km <- as_num(args$max_output_distance_km, 300)
min_output_coverage <- as_num(args$min_output_coverage, 0.995)
start_ma <- as_num(args$start_ma, 325)

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("pkgload is required to load the current package source.", call. = FALSE)
}
pkgload::load_all(pkg_root, quiet = TRUE)
if (!file.exists(slice_index) || !file.exists(tracks_path)) {
  stop("Missing palaeo slice index or PALEOMAP carrier-track table.", call. = FALSE)
}
index <- utils::read.csv(slice_index, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("time_ma", "file") %in% names(index))) {
  stop("palaeo_slice_index must contain time_ma and file.", call. = FALSE)
}
index$time_ma <- suppressWarnings(as.numeric(index$time_ma))
index <- index[is.finite(index$time_ma) & index$time_ma <= start_ma + 1e-8, , drop = FALSE]
index <- index[order(index$time_ma, decreasing = TRUE), , drop = FALSE]
if (nrow(index) < 2L || anyDuplicated(index$time_ma) || any(!file.exists(index$file))) {
  stop("palaeo_slice_index has duplicate ages, fewer than two usable slices, or a missing RDS.",
       call. = FALSE)
}
tracks <- readRDS(tracks_path)
if (!is.data.frame(tracks) || !"h3" %in% names(tracks)) {
  stop("PALEOMAP.RDS must be a data frame with stable h3 IDs.", call. = FALSE)
}

track_at_time <- function(time_ma) {
  suffix <- if (abs(time_ma) < 1e-8) "" else paste0("_", as.integer(round(time_ma)))
  cols <- paste0(c("lng", "lat"), suffix)
  if (!all(cols %in% names(tracks))) {
    stop("PALEOMAP carrier table has no coordinates at ", time_ma, " Ma.", call. = FALSE)
  }
  data.frame(track_id = as.character(tracks$h3),
             lon = suppressWarnings(as.numeric(tracks[[cols[[1L]]]])),
             lat = suppressWarnings(as.numeric(tracks[[cols[[2L]]]])),
             stringsAsFactors = FALSE)
}

carrier_dir <- safe_dir(file.path(output, "carrier_slices"))
rows <- vector("list", nrow(index))
for (ii in seq_len(nrow(index))) {
  time_ma <- index$time_ma[[ii]]
  earth <- as.data.frame(readRDS(index$file[[ii]]))
  required <- c("cell_id", "lon", "lat", "H_state", "land_area_km2")
  if (!all(required %in% names(earth))) {
    stop("Earth slice at ", time_ma, " Ma lacks: ",
         paste(setdiff(required, names(earth)), collapse = ", "), call. = FALSE)
  }
  matched <- HmscEcoEvo::hee_plate_carrier_match_grid(
    carriers = track_at_time(time_ma), grid = earth, time_ma = time_ma,
    max_match_distance_km = max_match_distance_km
  )
  active <- matched$track_id[matched$active_land]
  probe <- matrix(1, nrow = length(active), ncol = 1L,
                  dimnames = list(active, "coverage_probe"))
  rendered <- HmscEcoEvo::hee_plate_carrier_to_grid(
    carrier_values = probe, carrier_grid = matched, target_grid = earth,
    max_output_distance_km = max_output_distance_km
  )
  diag <- attr(matched, "diagnostics")
  diag$output_coverage_fraction <- rendered$diagnostics$output_coverage_fraction
  diag$output_nearest_carrier_km_q95 <- rendered$diagnostics$nearest_carrier_km_q95
  diag$file <- normalizePath(file.path(carrier_dir, paste0("carrier_", ma_slug(time_ma), "Ma.rds")),
                             winslash = "/", mustWork = FALSE)
  if (diag$output_coverage_fraction < min_output_coverage) {
    stop("Carrier-to-grid coverage below threshold at ", time_ma, " Ma: ",
         format(diag$output_coverage_fraction, digits = 5), " < ", min_output_coverage,
         ". Do not run Case05 v2 with an incomplete carrier map.", call. = FALSE)
  }
  saveRDS(list(carriers = matched, output_coverage = rendered$coverage,
               diagnostics = diag), diag$file, compress = "gzip")
  rows[[ii]] <- diag
  message("Prepared carrier grid at ", time_ma, " Ma; active carriers ",
          diag$n_active_land_carriers, "; 1 degree coverage ",
          formatC(diag$output_coverage_fraction, format = "f", digits = 4))
}
summary <- do.call(rbind, rows)
write_csv(summary, file.path(output, "case05_v2_plate_carrier_index.csv"))
write_csv(data.frame(
  parameter = c("plate_tracks", "palaeo_slice_index", "state_space",
                "environment_role", "output_role", "max_match_distance_km",
                "max_output_distance_km", "min_output_coverage",
                "scientific_boundary"),
  value = c(normalizePath(tracks_path, winslash = "/"),
            normalizePath(slice_index, winslash = "/"),
            "stable PALEOMAP H3 carrier IDs with time-specific palaeo coordinates",
            "complete 1 degree palaeo-land grid sampled onto active carriers",
            "carrier state interpolated to complete 1 degree palaeo-land grid only for display and GIS output",
            max_match_distance_km, max_output_distance_km, min_output_coverage,
            "Plate carriage is a coordinate/land-unit transfer. It is not active dispersal, a colonisation multiplier, or an independent probability term."),
  stringsAsFactors = FALSE
), file.path(output, "case05_v2_plate_carrier_metadata.csv"))
writeLines(c(
  "# Case05 v2 plate-carrier preparation",
  "",
  "The state space consists of stable PALEOMAP H3 carrier IDs, each with reconstructed coordinates at every time slice.",
  "At each time, the complete 1 degree palaeo-land grid supplies environmental covariates and habitat availability to the carriers.",
  "Carrier states are interpolated back to the complete 1 degree grid only after process updates, for GIS products. This output interpolation is not a second movement model.",
  "A time slice fails preparation when full-land output coverage is below the declared threshold; no target-centred fallback is used."
), file.path(output, "README.md"))
message("Case05 v2 carrier preparation complete: ", output)
