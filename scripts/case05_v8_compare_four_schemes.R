#!/usr/bin/env Rscript

# Compare four matched terrestrial location-density scenarios. These fields are
# location distributions, not cell occupancy or historical richness.
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1L]] else
  "C:/Users/Google/Documents/HMSC-HIST/outputs/HmscEcoEvo"
ids <- c(spherical_land = "case05_v8_full200_spherical_land_with_ldd_20260926",
  spherical_no_ldd = "case05_v8_full200_spherical_land_20260926",
  terrain_resistance = "case05_v8_full200_terrain_resistance_20260926",
  particle_topography = "case05_v8_full200_particle_topography_20260926")
if (!requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("terra", quietly = TRUE)) stop("ggplot2 and terra are required.")
files <- file.path(root, ids, "05_shared_scale_maps", "case05_v8_map_values.rds")
if (!all(file.exists(files))) stop("All four complete Case05 v8 map-value files are required.")
quality <- lapply(file.path(root, ids, "06_quality", "case05_v8_quality_gates.csv"),
  function(path) utils::read.csv(path, stringsAsFactors = FALSE))
if (any(vapply(quality, function(x) !all(x$pass), logical(1)))) {
  stop("At least one scheme failed quality gates.")
}
fields <- lapply(files, readRDS)
names(fields) <- names(ids)
times <- names(fields[[1L]]$location_density)
if (!all(vapply(fields, function(x) identical(names(x$location_density), times), logical(1)))) {
  stop("Scenario display times differ.")
}
out <- file.path(root, "case05_v8_four_scheme_comparison_20260926")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
pair <- list(c("terrain_resistance", "spherical_land"),
  c("particle_topography", "terrain_resistance"),
  c("particle_topography", "spherical_land"),
  c("spherical_land", "spherical_no_ldd"))
land <- lapply(seq_along(times), function(k) fields[[1L]]$landscape_weight[[k]] > 0)
if (any(vapply(seq_along(times), function(k) {
  !identical(land[[k]], fields[[2L]]$landscape_weight[[k]] > 0) ||
    !identical(land[[k]], fields[[3L]]$landscape_weight[[k]] > 0) ||
    !identical(land[[k]], fields[[4L]]$landscape_weight[[k]] > 0)
}, logical(1)))) stop("Scenario land masks differ.")
summary <- do.call(rbind, lapply(seq_along(times), function(k) {
  do.call(rbind, lapply(pair, function(p) {
    a <- fields[[p[[1L]]]]$location_density[[k]]
    b <- fields[[p[[2L]]]]$location_density[[k]]
    if (length(a) != 4050L || length(b) != 4050L ||
        any(a[!land[[k]]] != 0) || any(b[!land[[k]]] != 0)) {
      stop("Invalid terrestrial density at ", times[[k]], " Ma.")
    }
    data.frame(time_ma = as.numeric(times[[k]]), numerator = p[[1L]],
      denominator = p[[2L]], active_lineage_mass = sum(a),
      total_variation_per_lineage = 0.5 * sum(abs(a - b)) / sum(a),
      maximum_absolute_cell_difference = max(abs(a - b)),
      signed_land_difference = sum(a - b),
      land_cells = sum(land[[k]]), marine_cells = sum(!land[[k]]))
  }))
}))
utils::write.csv(summary, file.path(out, "four_scheme_difference_by_time.csv"), row.names = FALSE)

selected <- intersect(c(0, 20, 65, 100, 200, 300), as.numeric(times))
map_data <- list()
for (k in which(as.numeric(times) %in% selected)) {
  for (p in pair) {
    x <- fields[[p[[1L]]]]$location_density[[k]] -
      fields[[p[[2L]]]]$location_density[[k]]
    x[!land[[k]]] <- NA_real_
    map_data[[paste(p, collapse = "_minus_")]][[times[[k]]]] <- x
  }
}
common_limit <- max(abs(unlist(map_data)), na.rm = TRUE)
grid <- expand.grid(lon = seq(-178, 178, by = 4), lat = seq(-88, 88, by = 4))
index <- list()
for (name in names(map_data)) {
  folder <- file.path(out, name); dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  for (key in names(map_data[[name]])) {
    tab <- grid; tab$difference <- map_data[[name]][[key]]
    png <- file.path(folder, paste0("location_density_difference_", key, "Ma.png"))
    tif <- file.path(folder, paste0("location_density_difference_", key, "Ma.tif"))
    figure <- ggplot2::ggplot(tab, ggplot2::aes(lon, lat)) +
      ggplot2::geom_tile(ggplot2::aes(fill = difference), colour = NA) +
      ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7",
        high = "#b2182b", midpoint = 0, limits = c(-common_limit, common_limit),
        na.value = "white", name = "Density difference") +
      ggplot2::coord_equal(expand = FALSE) +
      ggplot2::labs(title = paste(name, key, "Ma"),
        subtitle = "Same terrestrial mask and symmetric colour scale across all comparisons",
        caption = "Lineage location-density difference, not occupancy or richness.",
        x = "Palaeogeographic longitude", y = "Palaeogeographic latitude") +
      ggplot2::theme_minimal(base_size = 12)
    ggplot2::ggsave(png, figure, width = 12, height = 6.4, dpi = 150)
    raster <- terra::rast(ncols = 90, nrows = 45, xmin = -180, xmax = 180,
      ymin = -90, ymax = 90, crs = "EPSG:4326")
    pixels <- terra::cellFromXY(raster, as.matrix(tab[, c("lon", "lat")]))
    value <- rep(NA_real_, terra::ncell(raster)); value[pixels] <- tab$difference
    terra::values(raster) <- value
    terra::writeRaster(raster, tif, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
    index[[length(index) + 1L]] <- data.frame(comparison = name,
      time_ma = as.numeric(key), png = png, geotiff = tif,
      scale_min = -common_limit, scale_max = common_limit)
  }
}
index <- do.call(rbind, index)
utils::write.csv(index, file.path(out, "four_scheme_map_index.csv"), row.names = FALSE)
stopifnot(all(file.exists(index$png)), all(file.exists(index$geotiff)))
message("Four-scheme comparison complete: ", out)
