#!/usr/bin/env Rscript

# Build sparse, time-specific Arias spherical movement kernels on stable
# PALEOMAP H3 carriers. Topographic resistance is read from the complete 1
# degree grid and attached to carriers before the kernel is built.

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
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(x, path)
  } else {
    utils::write.csv(x, path, row.names = FALSE)
  }
  invisible(path)
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
file_arg <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
if (!length(script)) script <- "scripts/case05_v2_prepare_movement_cache.R"
pkg_root <- normalizePath(file.path(dirname(normalizePath(script[[1L]], winslash = "/", mustWork = FALSE)), ".."), winslash = "/", mustWork = TRUE)
carrier_root <- cfg$carrier_root %||% file.path(pkg_root, "derived_inputs", "case05_v2_paleomap_h3_carrier_325_0Ma_1deg")
topo_index_path <- cfg$topographic_slice_index %||% file.path(pkg_root, "derived_inputs", "case04_plant200_palaeo_1deg_slices_20260919", "case04_topography_slice_index.csv")
output <- safe_dir(cfg$output %||% file.path(carrier_root, "movement_cache_arias_topographic"))
n_neighbours <- as_int(cfg$n_neighbours, 8L)
max_edge_distance_km <- as_num(cfg$max_edge_distance_km, 550)
omega_topography <- as_num(cfg$omega_topography, 1)
diffusion_variance <- as_num(cfg$diffusion_variance_rad2_per_myr, 0.0002)
source_emigration <- as_num(cfg$source_emigration_rate_per_myr, 0.05)

if (!requireNamespace("pkgload", quietly = TRUE) || !requireNamespace("RANN", quietly = TRUE)) {
  stop("pkgload and RANN are required.", call. = FALSE)
}
pkgload::load_all(pkg_root, quiet = TRUE)
carrier_index_path <- file.path(carrier_root, "case05_v2_plate_carrier_index.csv")
if (!file.exists(carrier_index_path) || !file.exists(topo_index_path)) {
  stop("Missing carrier or topographic index.", call. = FALSE)
}
carrier_index <- utils::read.csv(carrier_index_path, stringsAsFactors = FALSE)
topo_index <- utils::read.csv(topo_index_path, stringsAsFactors = FALSE)
topo_index$time_ma <- suppressWarnings(as.numeric(topo_index$time_ma))
if (any(!is.finite(carrier_index$time_ma)) || any(!file.exists(carrier_index$file)) ||
    any(!carrier_index$time_ma %in% topo_index$time_ma)) {
  stop("Carrier and topographic slices must have aligned valid time values.", call. = FALSE)
}

xyz <- function(lon, lat) {
  rad <- pi / 180; lon <- as.numeric(lon) * rad; lat <- as.numeric(lat) * rad
  cbind(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat))
}
great_circle_km <- function(lon1, lat1, lon2, lat2) {
  r <- pi / 180; p1 <- lat1 * r; p2 <- lat2 * r
  dp <- (lat2 - lat1) * r; dl <- (lon2 - lon1) * r
  a <- sin(dp / 2)^2 + cos(p1) * cos(p2) * sin(dl / 2)^2
  2 * 6371.0088 * atan2(sqrt(a), sqrt(pmax(0, 1 - a)))
}

cache_dir <- safe_dir(file.path(output, "matrices"))
rows <- vector("list", nrow(carrier_index))
for (ii in seq_len(nrow(carrier_index))) {
  time_ma <- carrier_index$time_ma[[ii]]
  carrier_blob <- readRDS(carrier_index$file[[ii]])
  carrier <- as.data.frame(carrier_blob$carriers)
  topo_file <- topo_index$file[[match(time_ma, topo_index$time_ma)]]
  topo <- as.data.frame(readRDS(topo_file))
  required_topo <- c("cell_id", "topographic_resistance", "topographic_permeability", "movement_allowed")
  if (!all(required_topo %in% names(topo))) stop("Topographic slice lacks required columns at ", time_ma, " Ma.")
  hit <- match(carrier$matched_grid_cell_id, topo$cell_id)
  carrier$topographic_resistance <- suppressWarnings(as.numeric(topo$topographic_resistance[hit]))
  carrier$topographic_permeability <- suppressWarnings(as.numeric(topo$topographic_permeability[hit]))
  carrier$movement_allowed <- as.integer(carrier$active_land & is.finite(carrier$topographic_resistance) &
                                            is.finite(carrier$topographic_permeability) &
                                            topo$movement_allowed[hit] > 0)
  active <- which(carrier$movement_allowed > 0L)
  if (length(active) <= n_neighbours) stop("Too few active carriers at ", time_ma, " Ma.")
  active_xyz <- xyz(carrier$paleo_lon[active], carrier$paleo_lat[active])
  nn <- RANN::nn2(active_xyz, query = active_xyz,
                  k = min(n_neighbours + 1L, length(active)))
  from <- rep(active, each = ncol(nn$nn.idx) - 1L)
  to <- active[as.vector(t(nn$nn.idx[, -1L, drop = FALSE]))]
  keep <- from != to
  dist <- great_circle_km(carrier$paleo_lon[from], carrier$paleo_lat[from],
                          carrier$paleo_lon[to], carrier$paleo_lat[to])
  keep <- keep & is.finite(dist) & dist > 0 & dist <= max_edge_distance_km
  from <- from[keep]; to <- to[keep]; dist <- dist[keep]
  mean_resistance <- (carrier$topographic_resistance[from] + carrier$topographic_resistance[to]) / 2
  edge <- data.frame(
    from_cell_id = carrier$track_id[from], to_cell_id = carrier$track_id[to],
    distance_km = dist,
    effective_cost_km = dist * exp(omega_topography * mean_resistance),
    stringsAsFactors = FALSE
  )
  cells <- data.frame(cell_id = carrier$track_id, lon = carrier$paleo_lon,
                      lat = carrier$paleo_lat, time_ma = time_ma,
                      movement_allowed = carrier$movement_allowed,
                      movement_landscape_weight = pmax(pmin(carrier$topographic_permeability, 1), 0),
                      cell_measure_weight = 1, stringsAsFactors = FALSE)
  kernel <- HmscEcoEvo::hee_dispersal_spherical_kernel(
    cells = cells, edges = edge, delta_t_myr = 1,
    diffusion_variance_rad2_per_myr = diffusion_variance,
    source_emigration_rate_per_myr = source_emigration,
    lineage = "carrier_shared", time_ma = time_ma,
    edge_cost_col = "effective_cost_km", spatial_domain = "global_explicit"
  )
  movement <- HmscEcoEvo::hee_dispersal_transition_matrix(kernel,
    cell_order = carrier$track_id, value_col = "movement_rate_per_myr", sparse = TRUE)
  cache_file <- file.path(cache_dir, paste0("movement_", formatC(time_ma, format = "f", digits = 6, drop0trailing = TRUE), "Ma.rds"))
  saveRDS(list(movement = movement, carrier = carrier,
               kernel_diagnostics = kernel$diagnostics,
               kernel_metadata = kernel$metadata), cache_file, compress = "gzip")
  rows[[ii]] <- data.frame(
    time_ma = time_ma, file = normalizePath(cache_file, winslash = "/", mustWork = FALSE),
    n_total_carriers = nrow(carrier), n_active_carriers = length(active),
    n_edges = nrow(kernel$edges), n_active_sources = kernel$diagnostics$n_active_sources,
    output_coverage_fraction = carrier_blob$diagnostics$output_coverage_fraction,
    source_emigration_rate_per_myr = source_emigration,
    diffusion_variance_rad2_per_myr = diffusion_variance,
    omega_topography = omega_topography,
    stringsAsFactors = FALSE
  )
  message("Cached movement at ", time_ma, " Ma: ", nrow(kernel$edges), " edges.")
}
out <- do.call(rbind, rows)
write_csv(out, file.path(output, "case05_v2_movement_cache_index.csv"))
write_csv(data.frame(
  field = c("kernel", "topographic_role", "plate_carriage_role", "tardis_role", "scientific_boundary"),
  value = c(
    "Arias spherical normal kernel on time-specific stable PALEOMAP H3 carrier coordinates",
    "topographic resistance enters once as a resistance-expanded edge distance",
    "stable carrier identity transports state with plate motion before active movement; no target-centred grid carriage",
    "time-ordered least-cost paths are post-hoc route diagnostics, not a second multiplier",
    "movement, establishment and persistence parameters are declared process scenarios unless independently calibrated"
  ), stringsAsFactors = FALSE
), file.path(output, "case05_v2_movement_metadata.csv"))
message("Case05 v2 movement cache complete: ", output)
