#!/usr/bin/env Rscript

# Aggregate complete Case05 v2 posterior shards and write full 1 degree GIS
# products.  Process states remain on stable PALEOMAP carrier identities; this
# script performs output-only interpolation to the complete palaeo-land grid.
# It never reuses that interpolation as biological plate transport or movement.

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
as_flag <- function(x, default = FALSE) tolower(as.character(x %||% default)) %in% c("true", "t", "1", "yes")
safe_dir <- function(x) { if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE); normalizePath(x, winslash = "/", mustWork = FALSE) }
write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("data.table", quietly = TRUE)) data.table::fwrite(x, path) else utils::write.csv(x, path, row.names = FALSE)
  invisible(path)
}
time_slug <- function(x) formatC(as.numeric(x), format = "f", digits = 6, drop0trailing = TRUE)

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
file_arg <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
if (!length(script)) script <- "scripts/case05_v2_aggregate_render_maps.R"
pkg_root <- normalizePath(cfg$pkg_root %||% file.path(dirname(normalizePath(script[[1L]], winslash = "/", mustWork = FALSE)), ".."), winslash = "/", mustWork = TRUE)
workspace <- normalizePath(file.path(pkg_root, ".."), winslash = "/", mustWork = TRUE)
output <- normalizePath(cfg$output %||% stop("--output is required.", call. = FALSE), winslash = "/", mustWork = TRUE)
carrier_root <- normalizePath(cfg$carrier_root %||% file.path(pkg_root, "derived_inputs", "case05_v2_paleomap_h3_carrier_325_0Ma_1deg"), winslash = "/", mustWork = TRUE)
earth_index_path <- normalizePath(cfg$earth_index %||% file.path(pkg_root, "derived_inputs", "case04_plant200_palaeo_1deg_slices_20260919", "case04_palaeo_earth_slice_index.csv"), winslash = "/", mustWork = TRUE)
max_output_distance_km <- as_num(cfg$max_output_distance_km, 1000)
interpolation_k <- as_int(cfg$interpolation_k, 4L)
allow_partial <- as_flag(cfg$allow_partial, FALSE)
render <- as_flag(cfg$render, TRUE)
write_geotiff <- as_flag(cfg$write_geotiff, TRUE)
if (!requireNamespace("pkgload", quietly = TRUE) || !requireNamespace("terra", quietly = TRUE)) {
  stop("pkgload and terra are required for Case05 v2 aggregation.", call. = FALSE)
}
pkgload::load_all(pkg_root, quiet = TRUE)
if (!is.finite(max_output_distance_km) || max_output_distance_km <= 0 || interpolation_k < 1L) {
  stop("max_output_distance_km and interpolation_k must be positive.", call. = FALSE)
}

dirs <- list(
  shards = file.path(output, "02_posterior_draws"),
  aggregate = safe_dir(file.path(output, "02_posterior_aggregate")),
  maps = safe_dir(file.path(output, "03_shared_scale_maps")),
  summaries = safe_dir(file.path(output, "05_summaries"))
)
shards <- list.dirs(dirs$shards, recursive = FALSE, full.names = TRUE)
complete <- shards[file.exists(file.path(shards, "16_state_space_inference", "draw_completion.csv"))]
if (!length(complete)) stop("No completed posterior draw shards are available.", call. = FALSE)
expected <- if (file.exists(file.path(output, "00_config", "case05_v2_formal_configuration.csv"))) {
  x <- utils::read.csv(file.path(output, "00_config", "case05_v2_formal_configuration.csv"), stringsAsFactors = FALSE)
  as_int(x$value[x$field == "n_response_draws"], length(complete))
} else length(complete)
if (!allow_partial && length(complete) != expected) {
  stop("Only ", length(complete), " of ", expected, " posterior draw shards are complete. Re-run after all shards finish, or set --allow_partial=true for a diagnostic aggregation.", call. = FALSE)
}
draw_ids <- sub("^draw_", "", basename(complete))

carrier_index <- utils::read.csv(file.path(carrier_root, "case05_v2_plate_carrier_index.csv"), stringsAsFactors = FALSE)
earth_index <- utils::read.csv(earth_index_path, stringsAsFactors = FALSE)
carrier_index$time_ma <- suppressWarnings(as.numeric(carrier_index$time_ma))
earth_index$time_ma <- suppressWarnings(as.numeric(earth_index$time_ma))
carrier_index <- carrier_index[order(carrier_index$time_ma, decreasing = TRUE), , drop = FALSE]
earth_index <- earth_index[match(carrier_index$time_ma, earth_index$time_ma), , drop = FALSE]
if (any(!is.finite(earth_index$time_ma)) || any(!file.exists(carrier_index$file)) || any(!file.exists(earth_index$file))) {
  stop("Carrier and complete Earth-grid indexes must contain the same valid times and files.", call. = FALSE)
}
default_map_times <- c(325, 300, 275, 250, 225, 200, 175, 150, 125, 100,
                       80, 65, 50, 40, 30, 20, 15, 10, 5, 0)
map_times <- if (!is.null(cfg$map_times) && nzchar(cfg$map_times)) {
  suppressWarnings(as.numeric(trimws(strsplit(cfg$map_times, ",", fixed = TRUE)[[1L]])))
} else default_map_times
map_times <- sort(unique(map_times[is.finite(map_times)]), decreasing = TRUE)
map_index_rows <- match(map_times, carrier_index$time_ma)
if (anyNA(map_index_rows)) {
  stop("Requested map_times are not present in the palaeoenvironment index: ",
       paste(map_times[is.na(map_index_rows)], collapse = ", "), call. = FALSE)
}
map_index_rows <- sort(unique(as.integer(map_index_rows)))

state_metrics <- c(
  "expected_sampled_surviving_lineage_richness",
  "potential_environmental_support_richness",
  "occupancy_weighted_shannon",
  "occupancy_weighted_simpson"
)
process_metrics <- c(
  "mean_arrival_probability_1myr",
  "mean_colonisation_probability_1myr",
  "mean_persistence_probability"
)
metric_meta <- data.frame(
  metric = c(state_metrics, process_metrics, "carrier_output_coverage", "nearest_carrier_distance_km", "diversity_defined"),
  label = c(
    "Expected sampled-surviving-lineage richness",
    "Potential ancestral environmental-support richness",
    "Occupancy-weighted Shannon lineage diversity",
    "Occupancy-weighted Simpson lineage diversity",
    "Mean arrival probability (per 1 Myr reporting interval)",
    "Mean colonisation probability (per 1 Myr reporting interval)",
    "Mean persistence probability (per 1 Myr reporting interval)",
    "Carrier output coverage", "Nearest active carrier distance", "Diversity-defined mask"
  ),
  unit = c("expected active lineages", "expected active lineages", "entropy", "Gini-Simpson", "probability", "probability", "probability", "0/1", "km", "0/1"),
  role = c(rep("state", length(state_metrics)), rep("interval_process", length(process_metrics)), "coverage", "coverage", "state"),
  boundary = c(
    "Process-constrained scenario summary; not all historical plant richness.",
    "Environmental support only; not historical occupancy probability.",
    "Defined only where expected lineage richness reaches the configured threshold.",
    "Defined only where expected lineage richness reaches the configured threshold.",
    "Interval process diagnostic on the older time slice, not a final occupancy probability.",
    "Interval process diagnostic on the older time slice, not a final occupancy probability.",
    "Interval process diagnostic on the older time slice, not a lineage survival estimate.",
    "Output interpolation support; zero means unresolved/no nearby carrier, not biological absence.",
    "GIS interpolation diagnostic, not a movement distance.",
    "Shows where diversity summaries are scientifically defined."
  ),
  stringsAsFactors = FALSE
)

read_draw_table <- function(shard, dir_name, file_name) {
  path <- file.path(shard, "16_state_space_inference", dir_name, file_name)
  if (!file.exists(path)) stop("Completed shard is missing ", path, call. = FALSE)
  readRDS(path)
}
summary_stats <- function(x) {
  if (!is.matrix(x)) x <- as.matrix(x)
  data.frame(mean = rowMeans(x, na.rm = TRUE), sd = apply(x, 1L, stats::sd, na.rm = TRUE),
             q025 = apply(x, 1L, stats::quantile, probs = .025, na.rm = TRUE, names = FALSE),
             q975 = apply(x, 1L, stats::quantile, probs = .975, na.rm = TRUE, names = FALSE),
             stringsAsFactors = FALSE)
}
matrix_from_draws <- function(tables, metric, ids) {
  out <- vapply(tables, function(z) {
    v <- as.numeric(z[[metric]][match(ids, z$track_id)])
    if (any(!is.finite(v))) v[!is.finite(v)] <- NA_real_
    v
  }, numeric(length(ids)))
  if (is.null(dim(out))) out <- matrix(out, ncol = 1L)
  out
}
aggregate_slice <- function(ii) {
  time_ma <- carrier_index$time_ma[[ii]]
  slug <- time_slug(time_ma)
  cache <- file.path(dirs$aggregate, paste0("state_aggregate_", slug, "Ma.rds"))
  if (file.exists(cache)) return(readRDS(cache))
  tables <- lapply(complete, read_draw_table, dir_name = "carrier_metrics", file_name = paste0("carrier_metrics_", slug, "Ma.rds"))
  ids <- as.character(tables[[1L]]$track_id)
  if (anyDuplicated(ids) || any(!vapply(tables, function(z) setequal(ids, z$track_id), logical(1)))) {
    stop("Carrier identifiers differ among posterior shards at ", time_ma, " Ma.", call. = FALSE)
  }
  base <- tables[[1L]][match(ids, tables[[1L]]$track_id), c("track_id", "time_ma", "paleo_lon", "paleo_lat", "active_land", "n_active_lineages"), drop = FALSE]
  for (metric in state_metrics) {
    ss <- summary_stats(matrix_from_draws(tables, metric, ids))
    names(ss) <- paste(metric, names(ss), sep = "__")
    base <- cbind(base, ss)
  }
  defined <- matrix_from_draws(tables, "diversity_defined", ids)
  base$diversity_defined__mean <- rowMeans(defined, na.rm = TRUE)
  saveRDS(base, cache, compress = "gzip")
  invisible(base)
}
aggregate_interval <- function(ii) {
  if (ii >= nrow(carrier_index)) return(invisible(NULL))
  older <- carrier_index$time_ma[[ii]]; younger <- carrier_index$time_ma[[ii + 1L]]
  cache <- file.path(dirs$aggregate, paste0("interval_aggregate_", time_slug(older), "_to_", time_slug(younger), "Ma.rds"))
  if (file.exists(cache)) return(readRDS(cache))
  stem <- paste0("carrier_interval_process_", time_slug(older), "_to_", time_slug(younger), "Ma.rds")
  tables <- lapply(complete, read_draw_table, dir_name = "carrier_interval_process", file_name = stem)
  ids <- as.character(tables[[1L]]$track_id)
  if (anyDuplicated(ids) || any(!vapply(tables, function(z) setequal(ids, z$track_id), logical(1)))) {
    stop("Carrier identifiers differ among interval-process shards at ", older, " -> ", younger, " Ma.", call. = FALSE)
  }
  base <- tables[[1L]][match(ids, tables[[1L]]$track_id), c("track_id", "time_from_ma", "time_to_ma", "paleo_lon", "paleo_lat", "active_land", "n_internal_process_steps"), drop = FALSE]
  for (metric in process_metrics) {
    ss <- summary_stats(matrix_from_draws(tables, metric, ids))
    names(ss) <- paste(metric, names(ss), sep = "__")
    base <- cbind(base, ss)
  }
  saveRDS(base, cache, compress = "gzip")
  invisible(base)
}

message("Aggregating ", length(complete), " complete Case05 v2 response draws across ", nrow(carrier_index), " time slices.")
invisible(lapply(seq_len(nrow(carrier_index)), aggregate_slice))
invisible(lapply(seq_len(nrow(carrier_index) - 1L), aggregate_interval))

# HMSC environmental support is evaluated directly on the complete 1 degree
# palaeoenvironment grid. It is intentionally not interpolated from the
# plate-carrier state space: only dynamic occupancy needs the carrier lattice.
direct_support_dir <- safe_dir(file.path(dirs$aggregate, "direct_environmental_support"))
direct_support_for <- function(ii) {
  time_ma <- carrier_index$time_ma[[ii]]; slug <- time_slug(time_ma)
  cache <- file.path(direct_support_dir, paste0("direct_support_", slug, "Ma.rds"))
  if (file.exists(cache)) return(readRDS(cache))
  earth <- as.data.frame(readRDS(earth_index$file[[ii]]))
  environment_exclude <- c("time_ma", "cell_id", "lon", "lat", "H_state", "land_area_km2")
  draw_support <- vector("list", length(complete))
  for (dd in seq_along(complete)) {
    shard <- complete[[dd]]
    state <- readRDS(file.path(shard, "16_state_space_inference", "carrier_states",
                               paste0("occupancy_", slug, "Ma.rds")))
    lineages <- as.character(state$lineages)
    response <- utils::read.csv(file.path(shard, "04_evolution",
                                           "ancestral_beta_response_master_time.csv"),
                                check.names = FALSE, stringsAsFactors = FALSE)
    axis <- intersect(setdiff(names(earth), environment_exclude), names(response))
    if (!length(axis)) stop("No shared environmental response axes at ", time_ma, " Ma.", call. = FALSE)
    take <- response[abs(suppressWarnings(as.numeric(response$time_ma)) - time_ma) <= 1e-8 &
                       response$lineage %in% lineages, , drop = FALSE]
    missing <- setdiff(lineages, take$lineage)
    if (length(missing)) {
      # Exact dated-tree nodes are queried on the younger side in the runner.
      # Reuse that explicit convention for full-grid environmental projection.
      epsilon_rows <- response[response$lineage %in% missing &
                                 response$time_ma < time_ma &
                                 response$time_ma >= time_ma - 1e-3, , drop = FALSE]
      if (nrow(epsilon_rows)) {
        epsilon_rows <- do.call(rbind, lapply(split(epsilon_rows, epsilon_rows$lineage), function(z) {
          z[which.max(z$time_ma), , drop = FALSE]
        }))
        take <- rbind(take, epsilon_rows)
      }
    }
    if (!setequal(lineages, take$lineage) || anyDuplicated(take$lineage)) {
      stop("Full-grid ancestral response is missing active lineages at ", time_ma,
           " Ma in draw ", basename(shard), call. = FALSE)
    }
    take <- take[match(lineages, take$lineage), , drop = FALSE]
    X <- as.matrix(earth[, axis, drop = FALSE]); B <- as.matrix(take[, axis, drop = FALSE])
    storage.mode(X) <- "double"; storage.mode(B) <- "double"
    eta <- X %*% t(B)
    complete_environment <- rowSums(is.finite(eta)) == ncol(eta)
    eta[!is.finite(eta)] <- NA_real_
    draw_support[[dd]] <- rowSums(stats::pnorm(eta), na.rm = TRUE)
    draw_support[[dd]][!complete_environment] <- NA_real_
  }
  mat <- do.call(cbind, draw_support)
  ss <- summary_stats(mat); names(ss) <- paste0("potential_environmental_support_richness__", names(ss))
  out <- cbind(earth[, c("cell_id", "time_ma", "lon", "lat", "H_state"), drop = FALSE], ss)
  attr(out, "axes") <- axis
  attr(out, "n_response_draws") <- length(complete)
  saveRDS(out, cache, compress = "gzip")
  out
}
direct_support <- lapply(map_index_rows, direct_support_for)
names(direct_support) <- time_slug(carrier_index$time_ma[map_index_rows])

# Fixed legends are calculated from all time-slice posterior means.  The 99.5th
# percentile controls display saturation only; raw GeoTIFF values are not
# clipped. This gives time-comparable colour scales without one extreme cell
# making an entire historical atlas visually blank.
scale_rows <- list(); si <- 0L
for (metric in c(state_metrics, process_metrics)) {
  pattern <- if (metric %in% state_metrics) "^state_aggregate_" else "^interval_aggregate_"
  files <- list.files(dirs$aggregate, pattern = pattern, full.names = TRUE)
  for (stat in c("mean", "sd", "q025", "q975")) {
    vals <- if (metric == "potential_environmental_support_richness") {
      unlist(lapply(direct_support, function(z) z[[paste0(metric, "__", stat)]]), use.names = FALSE)
    } else {
      unlist(lapply(files, function(f) readRDS(f)[[paste0(metric, "__", stat)]]), use.names = FALSE)
    }
    vals <- vals[is.finite(vals)]
    if (!length(vals)) next
    si <- si + 1L
    upper <- as.numeric(stats::quantile(vals, .995, names = FALSE))
    if (!is.finite(upper) || upper <= 0) upper <- max(vals, 1e-8)
    scale_rows[[si]] <- data.frame(metric = metric, statistic = stat, display_min = 0, display_max = upper,
                                   display_rule = "global_q995_of_this_summary_across_all_times",
                                   raw_values_clipped = FALSE, stringsAsFactors = FALSE)
  }
}
# Deterministic diagnostic layers also receive fixed legends across every age.
scale_rows[[length(scale_rows) + 1L]] <- data.frame(metric = "carrier_output_coverage", statistic = "deterministic", display_min = 0, display_max = 1,
                                                     display_rule = "fixed_probability_scale", raw_values_clipped = FALSE)
scale_rows[[length(scale_rows) + 1L]] <- data.frame(metric = "nearest_carrier_distance_km", statistic = "deterministic", display_min = 0, display_max = max_output_distance_km,
                                                     display_rule = "fixed_declared_interpolation_radius", raw_values_clipped = FALSE)
scale_rows[[length(scale_rows) + 1L]] <- data.frame(metric = "diversity_defined", statistic = "posterior_mean", display_min = 0, display_max = 1,
                                                     display_rule = "fixed_probability_scale", raw_values_clipped = FALSE)
display_scales <- do.call(rbind, scale_rows)
write_csv(display_scales, file.path(dirs$summaries, "case05_v2_display_scale_catalog.csv"))

make_raster <- function(grid, values) {
  r <- terra::rast(nrows = 180L, ncols = 360L, xmin = -180, xmax = 180, ymin = -90, ymax = 90, crs = "EPSG:4326")
  target <- terra::cellFromXY(r, cbind(grid$lon, grid$lat))
  vals <- rep(NA_real_, terra::ncell(r)); vals[target] <- as.numeric(values)
  terra::values(r) <- vals
  r
}
render_map <- function(r, png, title, unit, lim, divergent = FALSE) {
  dir.create(dirname(png), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(png, width = 1800, height = 1000, res = 160)
  on.exit(grDevices::dev.off(), add = TRUE)
  pal <- if (divergent) grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(128L) else grDevices::hcl.colors(128L, "YlOrRd", rev = FALSE)
  terra::plot(r, col = pal, range = lim, colNA = "grey84", main = title, axes = TRUE,
              plg = list(title = unit, cex = .78), xlab = "Palaeolongitude (degrees)", ylab = "Palaeolatitude (degrees)")
  graphics::mtext("Grey = ocean, unavailable habitat, or output beyond the declared carrier-interpolation radius.", side = 1L, line = 3, cex = .66)
}
output_map_rows <- list(); oi <- 0L
interpolate_metric <- function(carrier_tab, metric_column, output_operator) {
  value <- as.matrix(carrier_tab[[metric_column]])
  rownames(value) <- carrier_tab$track_id; colnames(value) <- metric_column
  # Conditional diversity is undefined below the configured expected-richness
  # threshold. Zeros are used only for interpolation mechanics and the full-grid
  # result is masked immediately after using the interpolated richness summary.
  value[!is.finite(value)] <- 0
  out <- as.numeric(HmscEcoEvo::hee_plate_carrier_apply_output_weights(
    carrier_values = value, output_weights = output_operator
  )[, 1L])
  list(values = out, coverage = output_operator$coverage,
       diagnostics = output_operator$diagnostics)
}
write_product <- function(values, grid, metric, stat, time_ma = NA_real_, time_from_ma = NA_real_, time_to_ma = NA_real_, coverage, diagnostics, mask = NULL) {
  if (!is.null(mask)) values[!mask] <- NA_real_
  meta <- metric_meta[match(metric, metric_meta$metric), , drop = FALSE]
  subdir <- safe_dir(file.path(dirs$maps, metric, stat))
  tag <- if (is.finite(time_ma)) paste0(time_slug(time_ma), "Ma") else paste0(time_slug(time_from_ma), "_to_", time_slug(time_to_ma), "Ma")
  tif <- file.path(subdir, paste0(metric, "_", stat, "_", tag, ".tif"))
  png <- file.path(subdir, paste0(metric, "_", stat, "_", tag, ".png"))
  r <- make_raster(grid, values)
  if (write_geotiff) terra::writeRaster(r, tif, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
  scale <- display_scales[display_scales$metric == metric & display_scales$statistic == stat, , drop = FALSE]
  finite_values <- values[is.finite(values)]
  lim <- if (nrow(scale)) {
    c(scale$display_min[[1L]], scale$display_max[[1L]])
  } else if (length(finite_values)) {
    range(finite_values)
  } else {
    c(0, 1)
  }
  if (length(lim) != 2L || any(!is.finite(lim)) || lim[[1L]] == lim[[2L]]) lim <- c(0, 1)
  title <- paste0("Case05 v2: ", meta$label, " (", stat, "); ", tag, "; forward process scenario")
  if (render) render_map(r, png, title, meta$unit, lim)
  oi <<- oi + 1L
  output_map_rows[[oi]] <<- data.frame(
    metric = metric, statistic = stat, time_ma = time_ma, time_from_ma = time_from_ma, time_to_ma = time_to_ma,
    GeoTIFF = if (write_geotiff) normalizePath(tif, winslash = "/", mustWork = FALSE) else NA_character_,
    PNG = if (render) normalizePath(png, winslash = "/", mustWork = FALSE) else NA_character_,
    unit = meta$unit, role = meta$role, boundary = meta$boundary,
    CRS = "EPSG:4326", resolution_degrees = 1, n_posterior_draws = length(complete),
    output_coverage_fraction = diagnostics$output_coverage_fraction,
    interpolation_k = interpolation_k, max_output_distance_km = max_output_distance_km,
    stringsAsFactors = FALSE
  )
  invisible(r)
}

for (ii in map_index_rows) {
  time_ma <- carrier_index$time_ma[[ii]]; slug <- time_slug(time_ma)
  carrier_blob <- readRDS(carrier_index$file[[ii]])
  earth <- as.data.frame(readRDS(earth_index$file[[ii]]))
  output_operator <- HmscEcoEvo::hee_plate_carrier_output_weights(
    carrier_grid = carrier_blob$carriers, target_grid = earth,
    k = interpolation_k, max_output_distance_km = max_output_distance_km
  )
  tab <- readRDS(file.path(dirs$aggregate, paste0("state_aggregate_", slug, "Ma.rds")))
  direct_env <- direct_support[[time_slug(time_ma)]]
  richness <- interpolate_metric(tab, "expected_sampled_surviving_lineage_richness__mean", output_operator)
  # All state summaries are interpolated from the carrier-state result to the
  # complete 1 degree palaeo-land grid. Diversity is deliberately masked where
  # expected richness is below 0.05: there is no meaningful Shannon/Simpson
  # index for a nearly empty reconstructed assemblage.
  for (metric in state_metrics) for (stat in c("mean", "sd", "q025", "q975")) {
    direct_environment <- metric == "potential_environmental_support_richness"
    got <- if (direct_environment) {
      list(values = as.numeric(direct_env[[paste0(metric, "__", stat)]]),
           coverage = data.frame(cell_id = earth$cell_id, output_covered = earth$H_state > 0,
                                 nearest_carrier_distance_km = NA_real_),
           diagnostics = data.frame(output_coverage_fraction = 1))
    } else if (metric == "expected_sampled_surviving_lineage_richness" && stat == "mean") richness else
      interpolate_metric(tab, paste0(metric, "__", stat), output_operator)
    mask <- NULL
    if (metric %in% c("occupancy_weighted_shannon", "occupancy_weighted_simpson")) mask <- is.finite(richness$values) & richness$values >= .05
    write_product(got$values, earth, metric, stat, time_ma = time_ma, coverage = got$coverage, diagnostics = got$diagnostics, mask = mask)
  }
  # Coverage diagnostics must use the exact same output radius/operator as the
  # metric maps. The carrier-preparation cache may contain a deliberately more
  # conservative screening radius, which is retained in its own metadata but
  # must not be confused with this GIS product's actual support domain.
  cov <- rep(NA_real_, nrow(earth)); cov[match(output_operator$coverage$cell_id, earth$cell_id)] <- as.numeric(output_operator$coverage$output_covered)
  distance <- rep(NA_real_, nrow(earth)); distance[match(output_operator$coverage$cell_id, earth$cell_id)] <- as.numeric(output_operator$coverage$nearest_carrier_distance_km)
  write_product(cov, earth, "carrier_output_coverage", "deterministic", time_ma = time_ma, coverage = richness$coverage, diagnostics = richness$diagnostics)
  write_product(distance, earth, "nearest_carrier_distance_km", "deterministic", time_ma = time_ma, coverage = richness$coverage, diagnostics = richness$diagnostics)
  defined <- interpolate_metric(tab, "diversity_defined__mean", output_operator)
  write_product(defined$values, earth, "diversity_defined", "posterior_mean", time_ma = time_ma, coverage = defined$coverage, diagnostics = defined$diagnostics)
  if (ii < nrow(carrier_index)) {
    younger <- carrier_index$time_ma[[ii + 1L]]
    proc <- readRDS(file.path(dirs$aggregate, paste0("interval_aggregate_", slug, "_to_", time_slug(younger), "Ma.rds")))
    for (metric in process_metrics) for (stat in c("mean", "sd", "q025", "q975")) {
      got <- interpolate_metric(proc, paste0(metric, "__", stat), output_operator)
      write_product(got$values, earth, metric, stat, time_from_ma = time_ma, time_to_ma = younger, coverage = got$coverage, diagnostics = got$diagnostics)
    }
  }
  message("Rendered Case05 v2 full-grid products for ", time_ma, " Ma.")
}

map_index <- do.call(rbind, output_map_rows)
write_csv(map_index, file.path(dirs$summaries, "case05_v2_map_index.csv"))
write_csv(data.frame(
  check = c("complete_posterior_draws", "stable_carrier_state", "full_1_degree_gis_grid", "terminal_anchor", "BioGeoBEARS_or_region_multiplier", "diversity_mask", "coverage_written"),
  passed = c(length(complete) == expected, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  value = c(length(complete), "PALEOMAP_H3_track_id", "EPSG:4326; 360x180", "none", "none", "expected_richness >= 0.05", "carrier coverage and nearest distance GeoTIFFs"),
  stringsAsFactors = FALSE
), file.path(dirs$summaries, "case05_v2_map_quality_gates.csv"))
writeLines(c(
  "# Case05 v2 full-grid map products",
  "",
  "The biological process was updated on stable PALEOMAP H3 carrier identities. Each GeoTIFF/PNG is an output-only inverse-distance interpolation to the complete 1 degree palaeo-land grid.",
  "No target-centred raster copying, terminal convex anchor, BioGeoBEARS multiplier, or regional-history multiplier is used.",
  "Grey cells are ocean, unavailable habitat, or cells farther than the declared carrier-output radius. Consult carrier_output_coverage and nearest_carrier_distance_km before interpreting a map.",
  "Expected sampled-surviving-lineage richness is a process-constrained forward-scenario result for the sampled extant tree, not a reconstruction of all historical plants.",
  "Arrival, colonisation, and persistence are interval diagnostics attached to their older time slice. They must not be read as final occupancy maps.",
  "Shannon and Simpson maps are masked when expected lineage richness is below 0.05."
), file.path(dirs$maps, "README.md"))
message("Case05 v2 aggregation and full-grid GIS rendering complete: ", output)
