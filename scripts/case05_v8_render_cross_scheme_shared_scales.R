#!/usr/bin/env Rscript

# Re-render every completed Case05 v8 metric on one metric-specific scale
# shared across all four schemes and all 66 display times.
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1L]] else
  "C:/Users/Google/Documents/HMSC-HIST/outputs/HmscEcoEvo"
ids <- c(spherical_no_ldd = "case05_v8_full200_spherical_land_20260926",
  spherical_land = "case05_v8_full200_spherical_land_with_ldd_20260926",
  terrain_resistance = "case05_v8_full200_terrain_resistance_20260926",
  particle_topography = "case05_v8_full200_particle_topography_20260926")
paths <- file.path(root, ids, "05_shared_scale_maps", "case05_v8_map_values.rds")
if (!all(file.exists(paths))) stop("All four complete Case05 v8 outputs are required.")
for (package in c("ggplot2", "scales", "terra")) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Missing package: ", package)
}
fields <- lapply(paths, readRDS); names(fields) <- names(ids)
times <- names(fields[[1L]]$location_density)
if (!all(vapply(fields, function(x) identical(names(x$location_density), times), logical(1)))) {
  stop("Time grids differ across schemes.")
}
metrics <- c("01_density", "02_log_density", "03_ancestral_eta", "04_land_fraction",
  "05_support_intensity", "06_support_shannon", "07_support_effective")
source_name <- c("location_density", "location_density_log10",
  "environmental_linear_predictor", "landscape_weight",
  "terminal_calibrated_lineage_support_intensity", "support_mixture_shannon",
  "support_mixture_effective_lineages")
names(source_name) <- metrics
value_at <- function(object, metric, k) {
  x <- as.numeric(object[[source_name[[metric]]]][[k]])
  if (length(x) != 4050L) stop("Metric has wrong grid length: ", metric)
  x
}
limits <- lapply(metrics, function(metric) {
  value <- unlist(lapply(fields, function(object) {
    unlist(lapply(seq_along(times), function(k) value_at(object, metric, k)))
  }), use.names = FALSE)
  value <- value[is.finite(value)]
  if (!length(value)) stop("No finite values for ", metric)
  if (metric == "03_ancestral_eta") {
    extent <- max(abs(value)); c(-extent, extent)
  } else if (metric == "02_log_density") {
    range(value)
  } else if (metric == "04_land_fraction") {
    c(0, 1)
  } else if (metric == "07_support_effective") {
    c(1, max(value))
  } else c(0, max(value))
})
names(limits) <- metrics
output <- file.path(root, "case05_v8_four_scheme_comparison_20260926",
  "07_cross_scheme_shared_scale_maps")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
grid <- expand.grid(lon = seq(-178, 178, by = 4), lat = seq(-88, 88, by = 4))
index <- vector("list", length(ids) * length(times) * length(metrics))
at <- 0L
for (scheme in names(fields)) {
  object <- fields[[scheme]]
  for (k in seq_along(times)) {
    land <- value_at(object, "04_land_fraction", k) > 0
    for (metric in metrics) {
      x <- value_at(object, metric, k)
      x[!land] <- NA_real_
      tab <- grid; tab$value <- x
      folder <- file.path(output, scheme, metric)
      dir.create(folder, recursive = TRUE, showWarnings = FALSE)
      filename <- paste0(metric, "_", times[[k]], "Ma")
      png <- file.path(folder, paste0(filename, ".png"))
      tif <- file.path(folder, paste0(filename, ".tif"))
      range <- limits[[metric]]
      scale <- if (metric == "03_ancestral_eta") {
        ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7",
          high = "#b2182b", midpoint = 0, limits = range,
          oob = scales::squish, na.value = "white", name = metric)
      } else ggplot2::scale_fill_viridis_c(option = "C", limits = range,
        oob = scales::squish, na.value = "white", name = metric)
      figure <- ggplot2::ggplot(tab, ggplot2::aes(lon, lat)) +
        ggplot2::geom_tile(ggplot2::aes(fill = value), colour = NA) + scale +
        ggplot2::coord_equal(expand = FALSE) +
        ggplot2::labs(title = paste(metric, "|", times[[k]], "Ma |", scheme),
          subtitle = "One legend range per metric across all four schemes and 66 time slices",
          x = "Palaeogeographic longitude", y = "Palaeogeographic latitude",
          caption = "Terrestrial location support or diagnostic; not occupancy or historical richness.") +
        ggplot2::theme_minimal(base_size = 11)
      ggplot2::ggsave(png, figure, width = 12, height = 6.4, dpi = 150)
      raster <- terra::rast(ncols = 90, nrows = 45, xmin = -180, xmax = 180,
        ymin = -90, ymax = 90, crs = "EPSG:4326")
      cells <- terra::cellFromXY(raster, as.matrix(tab[, c("lon", "lat")]))
      values <- rep(NA_real_, terra::ncell(raster)); values[cells] <- x
      terra::values(raster) <- values; names(raster) <- metric
      terra::writeRaster(raster, tif, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
      at <- at + 1L
      index[[at]] <- data.frame(scheme = scheme, metric = metric,
        time_ma = as.numeric(times[[k]]), png = png, geotiff = tif,
        scale_min = range[[1L]], scale_max = range[[2L]])
    }
  }
}
index <- do.call(rbind, index)
utils::write.csv(index, file.path(output, "cross_scheme_shared_scale_map_index.csv"),
  row.names = FALSE)
stopifnot(nrow(index) == 1848L, all(file.exists(index$png)),
  all(file.exists(index$geotiff)))
for (metric in metrics) {
  rows <- index[index$metric == metric, , drop = FALSE]
  stopifnot(length(unique(rows$scale_min)) == 1L,
    length(unique(rows$scale_max)) == 1L)
}
message("Shared-scale Case05 v8 atlas complete: ", output)
