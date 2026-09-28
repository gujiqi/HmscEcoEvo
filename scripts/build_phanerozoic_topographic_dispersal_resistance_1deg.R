#!/usr/bin/env Rscript

# Build a reproducible 1-degree Phanerozoic topographic dispersal input.
#
# This script reads the raw, unscaled continuous palaeo-elevation and land-mask
# slices, not the standardised HMSC predictor table. It writes a land-only
# resistance GeoTIFF, a sparse local connectivity graph, and a compact Case04/
# Case05-ready cell-time layer for every available time slice. Climate is
# deliberately absent: it belongs to environmental filtering rather than the
# structural movement resistance surface.

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || identical(x, "")) y else x
}
parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) next
    pair <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    out[[pair[[1L]]]] <- paste(pair[-1L], collapse = "=")
  }
  out
}

cli <- parse_args(args)
default_env <- file.path(
  "C:/Users/Google/Documents/HMSC-HIST",
  "Phanerozoic_Environment_540_0Ma_5Myr_1deg_landharmonized_v2",
  "08_rds_export",
  "Phanerozoic_Environment_540_0Ma_5Myr_Level2_landharmonized_v2_compact_float32.rds"
)
default_output <- file.path(
  "C:/Users/Google/Documents/HMSC-HIST/derived_inputs",
  "Phanerozoic_DispersalResistance_540_0Ma_5Myr_1deg_topography_v1"
)
env_rds <- normalizePath(cli$env_rds %||% default_env, winslash = "/", mustWork = TRUE)
output <- normalizePath(cli$output %||% default_output, winslash = "/", mustWork = FALSE)
relief_ref_m <- as.numeric(cli$relief_ref_m %||% 500)
step_ref_m <- as.numeric(cli$step_ref_m %||% 500)
omega_topo <- as.numeric(cli$omega_topo %||% 1)
neighbours <- as.integer(cli$neighbours %||% 8L)
weight_relief <- as.numeric(cli$weight_relief %||% 0.5)
weight_step <- as.numeric(cli$weight_step %||% 0.5)
overwrite <- identical(tolower(cli$overwrite %||% "false"), "true")

if (!requireNamespace("terra", quietly = TRUE)) {
  stop("Package 'terra' is required to write GeoTIFF resistance layers.", call. = FALSE)
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required to write the representative map panel.", call. = FALSE)
}

script_arg <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))][1L])
script_path <- normalizePath(ifelse(is.na(script_arg),
                                    "scripts/build_phanerozoic_topographic_dispersal_resistance_1deg.R",
                                    script_arg),
                             winslash = "/", mustWork = FALSE)
pkg_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
source_files <- list.files(file.path(pkg_root, "R"), pattern = "\\.R$", full.names = TRUE)
for (file in source_files) source(file, local = .GlobalEnv)

if (dir.exists(output) && length(list.files(output, all.files = TRUE, no.. = TRUE)) && !overwrite) {
  stop("Output exists and is non-empty: ", output,
       ". Use --overwrite=true only for a deliberate replacement.", call. = FALSE)
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
dirs <- list(
  geotiff = file.path(output, "geotiff"),
  cell_tables = file.path(output, "cell_tables"),
  graphs = file.path(output, "connectivity_graphs"),
  metadata = file.path(output, "metadata"),
  figures = file.path(output, "figures")
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

cube <- readRDS(env_rds)
needed <- c("elevation_m", "land_mask_dem", "land_area_km2")
available <- cube$list_variables()
absent <- setdiff(needed, available)
if (length(absent)) stop("Raw environment RDS lacks: ", paste(absent, collapse = ", "), call. = FALSE)
times <- sort(unique(as.numeric(cube$age_ma)), decreasing = TRUE)
if (!length(times) || any(!is.finite(times))) stop("No finite time slices found in env_rds.", call. = FALSE)
if (!is.null(cli$times) && nzchar(cli$times)) {
  requested_times <- suppressWarnings(as.numeric(strsplit(cli$times, ",", fixed = TRUE)[[1L]]))
  if (!length(requested_times) || any(!is.finite(requested_times))) {
    stop("--times must be a comma-separated list of finite ages in Ma.", call. = FALSE)
  }
  unmatched <- setdiff(requested_times, times)
  if (length(unmatched)) {
    stop("Requested ages are unavailable in env_rds: ",
         paste(unmatched, collapse = ", "), call. = FALSE)
  }
  times <- sort(unique(requested_times), decreasing = TRUE)
}

time_slug <- function(x) {
  gsub("\\.", "p", formatC(x, format = "f", digits = 3, drop0trailing = TRUE))
}
cell_id_from_coordinates <- function(lon, lat) {
  # Match the canonical `hee_make_paleo_grid()` convention used by Case04:
  # explicit signs prevent positive and negative coordinates from collapsing
  # into visually similar but non-identical cross-time keys.
  sprintf("lon%+.3f_lat%+.3f", as.numeric(lon), as.numeric(lat))
}
make_raster <- function(cells, value_col) {
  x <- cells[, c("lon", "lat", value_col), drop = FALSE]
  names(x)[[3L]] <- "value"
  r <- terra::rast(x, type = "xyz", crs = "EPSG:4326")
  names(r) <- value_col
  r
}
write_metadata <- function(path, values) {
  utils::write.csv(values, path, row.names = FALSE, na = "")
}

cell_index <- vector("list", length(times))
map_index <- vector("list", length(times))
graph_index <- vector("list", length(times))
diagnostic_index <- vector("list", length(times))
plot_cells <- list()
resistance_range <- c(Inf, -Inf)

message("Building 1-degree topographic resistance for ", length(times), " time slices.")
for (i in seq_along(times)) {
  tm <- times[[i]]
  message(sprintf("[%03d/%03d] %s Ma", i, length(times), tm))
  slice <- as.data.frame(cube$get_slice(tm, variables = needed))
  names(slice)[names(slice) == "age_ma"] <- "time_ma"
  required <- c("time_ma", "lon", "lat", "elevation_m", "land_mask_dem", "land_area_km2")
  if (!all(required %in% names(slice))) {
    stop("Slice ", tm, " Ma does not expose expected columns.", call. = FALSE)
  }
  slice$cell_id <- cell_id_from_coordinates(slice$lon, slice$lat)
  surface <- hee_dispersal_resistance_surface(
    landscape = slice,
    elevation_col = "elevation_m", land_col = "land_mask_dem",
    cell_col = "cell_id", time_col = "time_ma", lon_col = "lon", lat_col = "lat",
    neighbours = neighbours, resolution_deg = 1,
    relief_ref_m = relief_ref_m, absolute_elevation = FALSE,
    missing_elevation = "barrier"
  )
  graph <- hee_dispersal_connectivity_graph(
    surface, cell_col = "cell_id", time_col = "time_ma", lon_col = "lon", lat_col = "lat",
    elevation_col = "elevation_m", neighbours = neighbours, resolution_deg = 1,
    step_ref_m = step_ref_m, weight_relief = weight_relief,
    weight_step = weight_step, omega_topo = omega_topo
  )
  cells <- surface$cells
  cells$analysis_mode <- "derived_topographic_structural_resistance"
  cells$source_env_rds <- env_rds
  cells$resolution_deg <- 1
  cells$crs <- "EPSG:4326"
  cells$units_topographic_resistance <- "dimensionless_log_relief_index"
  cells$units_local_relief_m <- "m"

  slug <- time_slug(tm)
  cell_file <- file.path(dirs$cell_tables, paste0("topographic_resistance_", slug, "Ma_1deg.rds"))
  graph_file <- file.path(dirs$graphs, paste0("topographic_connectivity_", slug, "Ma_1deg.rds"))
  tif_file <- file.path(dirs$geotiff, paste0("topographic_resistance_", slug, "Ma_1deg.tif"))
  # gzip keeps one-slice graph reads practical for the dynamic occupancy
  # workflow while avoiding the long CPU cost of xz on 109 graph files.
  saveRDS(cells, cell_file, compress = "gzip")
  saveRDS(graph, graph_file, compress = "gzip")
  resistance_raster <- make_raster(cells, "topographic_resistance")
  terra::writeRaster(
    resistance_raster, tif_file, overwrite = TRUE,
    wopt = list(gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
  )

  keep <- cells[, c("cell_id", "time_ma", "lon", "lat", "elevation_m",
                    "land_mask_dem", "land_area_km2", "local_relief_m",
                    "topographic_resistance", "topographic_permeability",
                    "movement_allowed", "terrain_missing"), drop = FALSE]
  cell_index[[i]] <- keep
  valid_resistance <- keep$topographic_resistance[is.finite(keep$topographic_resistance)]
  if (length(valid_resistance)) {
    resistance_range[[1L]] <- min(resistance_range[[1L]], min(valid_resistance))
    resistance_range[[2L]] <- max(resistance_range[[2L]], max(valid_resistance))
  }
  map_index[[i]] <- data.frame(
    filename = normalizePath(tif_file, winslash = "/", mustWork = TRUE),
    metric = "topographic_resistance", time_ma = tm,
    unit = "dimensionless_log_relief_index", crs = "EPSG:4326",
    resolution_deg = 1, source = env_rds,
    analysis_mode = "derived_topographic_structural_resistance",
    stringsAsFactors = FALSE
  )
  graph_index[[i]] <- data.frame(
    filename = normalizePath(graph_file, winslash = "/", mustWork = TRUE),
    metric = "local_topographic_connectivity_graph", time_ma = tm,
    edge_count = nrow(graph$edges), source = env_rds,
    analysis_mode = "derived_topographic_structural_resistance",
    stringsAsFactors = FALSE
  )
  diagnostic_index[[i]] <- data.frame(
    time_ma = tm,
    n_total_cells = nrow(cells),
    n_land_cells = sum(cells$land_mask_dem > 0, na.rm = TRUE),
    n_movement_cells = sum(cells$movement_allowed > 0, na.rm = TRUE),
    n_missing_land_elevation = sum(cells$terrain_missing, na.rm = TRUE),
    n_local_edges = nrow(graph$edges),
    min_resistance = if (length(valid_resistance)) min(valid_resistance) else NA_real_,
    max_resistance = if (length(valid_resistance)) max(valid_resistance) else NA_real_,
    stringsAsFactors = FALSE
  )
  if (tm %in% c(0, 65, 200, 325, 540)) {
    plot_cells[[as.character(tm)]] <- keep[, c("lon", "lat", "time_ma", "topographic_resistance"), drop = FALSE]
  }
}

all_cells <- do.call(rbind, cell_index)
combined_file <- file.path(output, "topographic_dispersal_resistance_540_0Ma_5Myr_1deg_v1.rds")
saveRDS(
  list(
    cells = all_cells,
    metadata = list(
      source_env_rds = env_rds,
      time_direction = "time_ma: larger values are older",
      grid_resolution_deg = 1,
      crs = "EPSG:4326",
      local_relief_formula = "median(abs(z_neighbour - z_cell)) over valid land 8-neighbours",
      resistance_formula = "log1p(local_relief_m / relief_ref_m)",
      edge_barrier_formula = "weight_relief * mean(cell resistance) + weight_step * log1p(elevation_step_m / step_ref_m)",
      edge_cost_formula = "distance_km * exp(omega_topo * barrier_index)",
      relief_ref_m = relief_ref_m, step_ref_m = step_ref_m,
      weight_relief = weight_relief, weight_step = weight_step,
      omega_topo = omega_topo, neighbours = neighbours,
      ocean_rule = "Ocean and land cells with missing elevation have no local movement edge.",
      scientific_boundary = "Generic structural P2 resistance scenario, not a calibrated taxon-specific dispersal, occurrence, or gene-flow estimate. Climate suitability, plate transport, and rare long-distance dispersal are separate components."
    )
  ), combined_file, compress = "gzip"
)

map_index <- do.call(rbind, map_index)
graph_index <- do.call(rbind, graph_index)
diagnostics <- do.call(rbind, diagnostic_index)
write_metadata(file.path(output, "topographic_dispersal_resistance_map_index.csv"), map_index)
write_metadata(file.path(output, "topographic_connectivity_graph_index.csv"), graph_index)
write_metadata(file.path(output, "topographic_dispersal_resistance_diagnostics.csv"), diagnostics)

settings <- data.frame(
  field = c("env_rds", "combined_cell_time_rds", "time_slices", "time_ma_direction",
            "crs", "resolution_deg", "neighbours", "relief_ref_m", "step_ref_m",
            "weight_relief", "weight_step", "omega_topo", "ocean_rule",
            "absolute_elevation_used", "scientific_boundary"),
  value = c(env_rds, normalizePath(combined_file, winslash = "/", mustWork = TRUE),
            length(times), "larger values are older", "EPSG:4326", 1, neighbours,
            relief_ref_m, step_ref_m, weight_relief, weight_step, omega_topo,
            "no local movement edge", FALSE,
            "Generic topographic structural-resistance input to the dispersal process; not calibrated realised dispersal or gene flow."),
  stringsAsFactors = FALSE
)
write_metadata(file.path(dirs$metadata, "topographic_dispersal_resistance_metadata.csv"), settings)
write_metadata(file.path(dirs$metadata, "surface_function_metadata.csv"), surface$metadata)
write_metadata(file.path(dirs$metadata, "connectivity_function_metadata.csv"), graph$metadata)

if (length(plot_cells) && all(is.finite(resistance_range))) {
  panel <- do.call(rbind, plot_cells)
  panel <- panel[is.finite(panel$topographic_resistance), , drop = FALSE]
  panel$time_label <- paste0(panel$time_ma, " Ma")
  p <- ggplot2::ggplot(panel, ggplot2::aes(x = lon, y = lat, fill = topographic_resistance)) +
    ggplot2::geom_tile(width = 1, height = 1) +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::facet_wrap(~time_label, ncol = 2) +
    ggplot2::scale_fill_viridis_c(
      option = "C", limits = resistance_range, oob = scales::squish,
      name = "Topographic\nresistance\n(log relief)"
    ) +
    ggplot2::labs(
      title = "Continuous 1-degree topographic dispersal resistance",
      subtitle = "Land-only structural P2 input; identical colour range across time slices",
      x = "Longitude (degrees)", y = "Latitude (degrees)"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
  ggplot2::ggsave(file.path(dirs$figures, "topographic_resistance_representative_times.png"),
                  p, width = 12, height = 8, dpi = 220)
}

message("Completed ", length(times), " time slices.")
message("Combined Case04/Case05 layer: ", normalizePath(combined_file, winslash = "/", mustWork = TRUE))
message("GeoTIFF index: ", normalizePath(file.path(output, "topographic_dispersal_resistance_map_index.csv"), winslash = "/", mustWork = TRUE))
