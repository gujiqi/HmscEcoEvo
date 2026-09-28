#!/usr/bin/env Rscript

# Build the Case05 v3 active-dispersal cache.
#
# The active P2 transition is deliberately narrow in scope:
# K(i -> j) is an Arias-style spherical diffusion allocation weighted once
# by the destination landscape measure (cell area x permeability).
#
# It does NOT put topographic least-cost distance into the same kernel. A
# TARDIS-style least-cost graph is a separate route diagnostic, constructed
# after this cache, so the same topography is not counted twice as both a
# target penalty and a distance penalty.

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
if (!length(script)) script <- "scripts/case05_v3_prepare_arias_target_landscape_cache.R"
pkg_root <- normalizePath(file.path(dirname(normalizePath(script[[1L]], winslash = "/", mustWork = FALSE)), ".."), winslash = "/", mustWork = TRUE)
carrier_root <- cfg$carrier_root %||% file.path(pkg_root, "derived_inputs", "case05_v2_paleomap_h3_carrier_325_0Ma_1deg")
topo_index_path <- cfg$topographic_slice_index %||% file.path(pkg_root, "derived_inputs", "case04_plant200_palaeo_1deg_slices_20260919", "case04_topography_slice_index.csv")
output <- safe_dir(cfg$output %||% file.path(carrier_root, "movement_cache_arias_target_landscape"))
candidate_neighbours <- as_int(cfg$candidate_neighbours, 128L)
max_edge_distance_km <- as_num(cfg$max_edge_distance_km, 450)
diffusion_variance <- as_num(cfg$diffusion_variance_rad2_per_myr, 0.0002)
source_emigration <- as_num(cfg$source_emigration_rate_per_myr, 0.05)
strict_radius_search <- as_flag(cfg$strict_radius_search, TRUE)
requested_times <- if (!is.null(cfg$time_ma) && nzchar(cfg$time_ma)) suppressWarnings(as.numeric(trimws(strsplit(cfg$time_ma, ",", fixed = TRUE)[[1L]]))) else numeric()

if (!requireNamespace("pkgload", quietly = TRUE) || !requireNamespace("RANN", quietly = TRUE)) {
  stop("pkgload and RANN are required.", call. = FALSE)
}
if (candidate_neighbours < 2L || !is.finite(max_edge_distance_km) || max_edge_distance_km <= 0 ||
    !is.finite(diffusion_variance) || diffusion_variance <= 0 ||
    !is.finite(source_emigration) || source_emigration < 0) {
  stop("Invalid candidate graph or diffusion settings.", call. = FALSE)
}
pkgload::load_all(pkg_root, quiet = TRUE)
carrier_index_path <- file.path(carrier_root, "case05_v2_plate_carrier_index.csv")
if (!file.exists(carrier_index_path) || !file.exists(topo_index_path)) {
  stop("Missing carrier or topographic slice index.", call. = FALSE)
}
carrier_index <- utils::read.csv(carrier_index_path, stringsAsFactors = FALSE)
topo_index <- utils::read.csv(topo_index_path, stringsAsFactors = FALSE)
carrier_index$time_ma <- suppressWarnings(as.numeric(carrier_index$time_ma))
topo_index$time_ma <- suppressWarnings(as.numeric(topo_index$time_ma))
if (any(!is.finite(carrier_index$time_ma)) || any(!file.exists(carrier_index$file)) ||
    any(!carrier_index$time_ma %in% topo_index$time_ma)) {
  stop("Carrier and topographic indexes must contain aligned valid slices.", call. = FALSE)
}
if (length(requested_times)) {
  carrier_index <- carrier_index[carrier_index$time_ma %in% requested_times, , drop = FALSE]
  if (!nrow(carrier_index) || !setequal(carrier_index$time_ma, requested_times)) {
    stop("Every requested --time_ma value must exist in the carrier index.", call. = FALSE)
  }
}
carrier_index <- carrier_index[order(carrier_index$time_ma, decreasing = TRUE), , drop = FALSE]

xyz <- function(lon, lat) {
  rad <- pi / 180
  lon <- as.numeric(lon) * rad; lat <- as.numeric(lat) * rad
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
  required_topo <- c("cell_id", "topographic_permeability", "movement_allowed")
  if (!all(required_topo %in% names(topo))) {
    stop("Topographic slice lacks required columns at ", time_ma, " Ma.", call. = FALSE)
  }
  hit <- match(carrier$matched_grid_cell_id, topo$cell_id)
  carrier$topographic_permeability <- suppressWarnings(as.numeric(topo$topographic_permeability[hit]))
  carrier$movement_allowed <- as.integer(
    carrier$active_land & is.finite(carrier$topographic_permeability) &
      topo$movement_allowed[hit] > 0 & is.finite(carrier$land_area_km2) & carrier$land_area_km2 > 0
  )
  active <- which(carrier$movement_allowed > 0L)
  if (length(active) <= candidate_neighbours) {
    stop("Too few active carriers at ", time_ma, " Ma for candidate_neighbours.", call. = FALSE)
  }
  active_xyz <- xyz(carrier$paleo_lon[active], carrier$paleo_lat[active])
  # This is a sparse numerical integration of the continuous spherical
  # kernel. The requested k is deliberately checked: if its outer neighbour
  # is still inside the declared tail radius, the candidate graph is not
  # sufficiently complete and a formal cache must fail rather than silently
  # omit high-weight targets.
  nn <- RANN::nn2(active_xyz, query = active_xyz,
                  k = min(candidate_neighbours + 1L, length(active)))
  from <- rep(active, each = ncol(nn$nn.idx) - 1L)
  to <- active[as.vector(t(nn$nn.idx[, -1L, drop = FALSE]))]
  distance <- great_circle_km(carrier$paleo_lon[from], carrier$paleo_lat[from],
                              carrier$paleo_lon[to], carrier$paleo_lat[to])
  keep <- from != to & is.finite(distance) & distance > 0 & distance <= max_edge_distance_km
  edge <- data.frame(
    from_cell_id = carrier$track_id[from[keep]],
    to_cell_id = carrier$track_id[to[keep]],
    distance_km = distance[keep],
    stringsAsFactors = FALSE
  )
  kth_distance <- great_circle_km(
    carrier$paleo_lon[active], carrier$paleo_lat[active],
    carrier$paleo_lon[active[nn$nn.idx[, ncol(nn$nn.idx)]]],
    carrier$paleo_lat[active[nn$nn.idx[, ncol(nn$nn.idx)]]]
  )
  n_sources_truncated <- sum(is.finite(kth_distance) & kth_distance < max_edge_distance_km - 1e-8)
  if (strict_radius_search && n_sources_truncated > 0L) {
    stop("Candidate-neighbour search is truncated inside the ", max_edge_distance_km,
         " km tail radius at ", time_ma, " Ma for ", n_sources_truncated,
         " sources. Increase --candidate_neighbours; do not treat this as a formal kernel.", call. = FALSE)
  }
  cells <- data.frame(
    cell_id = carrier$track_id, lon = carrier$paleo_lon, lat = carrier$paleo_lat,
    time_ma = time_ma, movement_allowed = carrier$movement_allowed,
    movement_landscape_weight = pmax(pmin(carrier$topographic_permeability, 1), 0),
    cell_area_km2 = carrier$land_area_km2, stringsAsFactors = FALSE
  )
  landscape <- HmscEcoEvo::hee_dispersal_landscape_weights(
    cells, hard_mask_col = "movement_allowed",
    movement_weight_col = "movement_landscape_weight",
    cell_area_col = "cell_area_km2"
  )
  kernel <- HmscEcoEvo::hee_dispersal_spherical_kernel(
    cells = landscape, edges = edge, delta_t_myr = 1,
    diffusion_variance_rad2_per_myr = diffusion_variance,
    source_emigration_rate_per_myr = source_emigration,
    lineage = "carrier_shared", time_ma = time_ma,
    # Critical: no effective-cost column. Topography enters exactly once as
    # a target landscape weight; least-cost paths remain a separate TARDIS
    # diagnostic layer.
    edge_cost_col = NULL, spatial_domain = "global_explicit"
  )
  movement <- HmscEcoEvo::hee_dispersal_transition_matrix(
    kernel, cell_order = carrier$track_id, value_col = "movement_rate_per_myr", sparse = TRUE
  )
  cache_file <- file.path(cache_dir, paste0("movement_", time_slug(time_ma), "Ma.rds"))
  saveRDS(list(
    movement = movement, carrier = carrier, landscape = landscape,
    kernel_diagnostics = kernel$diagnostics, kernel_metadata = kernel$metadata,
    candidate_diagnostics = data.frame(
      candidate_neighbours = candidate_neighbours, max_edge_distance_km = max_edge_distance_km,
      n_active_sources = length(active), n_sources_truncated = n_sources_truncated,
      kth_distance_km_q05 = as.numeric(stats::quantile(kth_distance, .05)),
      kth_distance_km_q50 = as.numeric(stats::quantile(kth_distance, .50)),
      kth_distance_km_q95 = as.numeric(stats::quantile(kth_distance, .95)),
      stringsAsFactors = FALSE
    )
  ), cache_file, compress = "gzip")
  sigma_km <- sqrt(diffusion_variance) * 6371.0088
  rows[[ii]] <- data.frame(
    time_ma = time_ma, file = normalizePath(cache_file, winslash = "/", mustWork = FALSE),
    n_total_carriers = nrow(carrier), n_active_carriers = length(active),
    n_edges = nrow(kernel$edges), n_active_sources = kernel$diagnostics$n_active_sources,
    n_sources_truncated = n_sources_truncated,
    output_coverage_fraction = carrier_blob$diagnostics$output_coverage_fraction,
    source_emigration_rate_per_myr = source_emigration,
    diffusion_variance_rad2_per_myr = diffusion_variance,
    sigma_km_per_sqrt_myr = sigma_km,
    tail_radius_sigma = max_edge_distance_km / sigma_km,
    stringsAsFactors = FALSE
  )
  message("Cached v3 Arias target-landscape kernel at ", time_ma, " Ma: ", nrow(kernel$edges), " edges.")
}
out <- do.call(rbind, rows)
write_csv(out, file.path(output, "case05_v3_movement_cache_index.csv"))
write_csv(data.frame(
  field = c("active_kernel", "distance_term", "landscape_term", "topography_count",
            "candidate_tail", "plate_carriage_role", "tardis_role", "scientific_boundary"),
  value = c(
    "Arias-style spherical-normal target-by-source transition on time-specific PALEOMAP H3 carriers",
    "great-circle distance only; no least-cost or resistance-expanded distance in the active kernel",
    "destination area x topographic-permeability, normalized by source",
    "one: topography is a declared destination-permeability surface only",
    paste0(max_edge_distance_km, " km (", formatC(max_edge_distance_km / (sqrt(diffusion_variance) * 6371.0088), digits = 3, format = "f"), " Gaussian SD); strict candidate-radius check=", strict_radius_search),
    "stable carrier identity provides plate carriage before active movement; it is not organismal dispersal",
    "TARDIS-style time-ordered least-cost paths are separate route diagnostics and are not multiplied into movement, colonisation, or persistence",
    "source-emigration, establishment and persistence parameters remain declared process scenarios unless calibrated against independent endpoint or fossil evidence"
  ), stringsAsFactors = FALSE
), file.path(output, "case05_v3_movement_metadata.csv"))
message("Case05 v3 Arias target-landscape cache complete: ", output)

