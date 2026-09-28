#!/usr/bin/env Rscript

# Four time-resolved quantitative map plates from saved Case05/Case07 outputs.
# Figures defend distinct claims: location, predictive movement, support-mixture
# diversity, and ecological refugium candidates. None is occupancy richness.
stopifnot(all(vapply(c("ggplot2", "ragg", "svglite", "terra", "scales"),
                     requireNamespace, logical(1), quietly = TRUE)))

root <- Sys.getenv("HMSCEE_ROOT", "")
if (!nzchar(root)) root <- normalizePath(file.path(getwd(), ".."), winslash = "/")
archive <- file.path(root, "outputs", "HmscEcoEvo")
out <- file.path(archive, "manuscript_results_six_process_20260927")
fig_dir <- file.path(out, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
case05 <- file.path(archive, "case05_v8_four_scheme_comparison_20260926")
case07 <- file.path(archive, "case07_patch_evidence_ladder_20260926")
field_path <- file.path(archive, "case05_v8_full200_terrain_resistance_20260926",
                        "05_shared_scale_maps", "case05_v8_map_values.rds")
index_path <- file.path(case05, "07_cross_scheme_shared_scale_maps",
                        "cross_scheme_shared_scale_map_index.csv")
movement_index_path <- file.path(case07, "07_dispersal_source_corridor_maps",
                                 "plant200_source_corridor_map_index.csv")
refuge_index_path <- file.path(case07, "06_real_time_maps",
                               "case07_real_time_map_index.csv")
diagnostic_path <- file.path(case07, "06_real_time_maps",
                             "case07_real_time_diagnostics.csv")
required <- c(field_path, index_path, movement_index_path,
              refuge_index_path, diagnostic_path)
if (!all(file.exists(required))) stop("Missing saved Case05/Case07 input: ",
                                     paste(required[!file.exists(required)], collapse = "; "))

field <- readRDS(field_path)
shared_index <- utils::read.csv(index_path, stringsAsFactors = FALSE)
movement_index <- utils::read.csv(movement_index_path, stringsAsFactors = FALSE)
refuge_index <- utils::read.csv(refuge_index_path, stringsAsFactors = FALSE)
diagnostic <- utils::read.csv(diagnostic_path, stringsAsFactors = FALSE)
times <- c(300, 200, 100, 65, 20, 5)
refuge_times <- c(245, 170, 125, 55, 25, 20)
focal <- "Camellia_sinensis"
grid <- expand.grid(lon = seq(-178, 178, by = 4), lat = seq(-88, 88, by = 4))
stopifnot(nrow(grid) == 4050L, all(as.character(times) %in% names(field$location_density)))

metric_limits <- function(metric) {
  z <- shared_index[shared_index$metric == metric, c("scale_min", "scale_max")]
  if (!nrow(z) || length(unique(z$scale_min)) != 1L ||
      length(unique(z$scale_max)) != 1L) {
    stop("Expected one all-scheme, all-time scale for ", metric)
  }
  c(z$scale_min[[1L]], z$scale_max[[1L]])
}

land_for <- function(time_ma) {
  at <- match(as.character(time_ma), names(field$location_density))
  z <- as.numeric(field$landscape_weight[[at]]) > 0
  if (length(z) != nrow(grid)) stop("Land mask grid mismatch.")
  z
}

base_theme <- ggplot2::theme_minimal(base_size = 12, base_family = "Arial") +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    panel.background = ggplot2::element_rect(fill = "white", colour = NA),
    strip.text = ggplot2::element_text(face = "bold", colour = "#243B50", size = 11),
    plot.title = ggplot2::element_text(face = "bold", colour = "#243B50", size = 17),
    plot.subtitle = ggplot2::element_text(colour = "#627681", size = 11),
    plot.caption = ggplot2::element_text(colour = "#53656F", size = 9),
    axis.title = ggplot2::element_text(colour = "#243B50", size = 10),
    axis.text = ggplot2::element_text(size = 8),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(size = 10),
    legend.text = ggplot2::element_text(size = 9),
    plot.margin = ggplot2::margin(8, 12, 8, 8)
  )

save_plate <- function(plot, id) {
  ggplot2::ggsave(file.path(fig_dir, paste0(id, ".png")), plot,
                  width = 12.8, height = 8.6, dpi = 400,
                  device = ragg::agg_png, bg = "white")
  ggplot2::ggsave(file.path(fig_dir, paste0(id, ".svg")), plot,
                  width = 12.8, height = 8.6,
                  device = svglite::svglite, bg = "white")
}

map_frame <- function(time_ma, value, mark_zero_as_unavailable = FALSE) {
  if (length(value) != nrow(grid)) stop("Metric grid mismatch at ", time_ma)
  land <- land_for(time_ma)
  tab <- grid
  tab$time_ma <- time_ma
  panel_label <- function(t) paste0(t, " Ma")
  tab$panel <- factor(panel_label(time_ma),
                      levels = vapply(times, panel_label, character(1)))
  tab$land <- land
  tab$value <- as.numeric(value)
  tab$value[!land] <- NA_real_
  if (mark_zero_as_unavailable) tab$value[!is.na(tab$value) & tab$value <= 0] <- NA_real_
  tab
}

numeric_plate <- function(tab, title, subtitle, legend_title, limits,
                          caption, breaks = NULL, palette = "C") {
  ggplot2::ggplot(tab, ggplot2::aes(lon, lat)) +
    ggplot2::geom_tile(data = tab[tab$land, , drop = FALSE],
                       fill = "#E8ECEA", width = 4, height = 4) +
    ggplot2::geom_tile(ggplot2::aes(fill = value), width = 4, height = 4) +
    ggplot2::scale_fill_viridis_c(option = palette, limits = limits,
                                  breaks = breaks, oob = scales::squish,
                                  na.value = NA, name = legend_title,
                                  guide = ggplot2::guide_colourbar(barwidth = 12,
                                                                   barheight = 0.55)) +
    ggplot2::facet_wrap(~panel, ncol = 3, drop = FALSE) +
    ggplot2::coord_fixed(expand = FALSE) +
    ggplot2::scale_x_continuous(breaks = c(-120, 0, 120), limits = c(-180, 180)) +
    ggplot2::scale_y_continuous(breaks = c(-60, 0, 60), limits = c(-90, 90)) +
    ggplot2::labs(title = title, subtitle = subtitle,
                  x = "Palaeolongitude", y = "Palaeolatitude",
                  caption = caption) + base_theme
}

# Tree-conditioned geographic locations under one declared movement kernel.
location <- do.call(rbind, lapply(times, function(time_ma) {
  z <- as.numeric(field$location_density[[as.character(time_ma)]])
  tab <- map_frame(time_ma, log10(pmax(z, 1e-12)))
  tab$value[tab$land & z <= 0] <- NA_real_
  tab
}))
fig_location <- numeric_plate(
  location, "Tree-conditioned lineage locations through time",
  "Plant200 | terrain-resistance movement | common display scale across all six slices",
  expression(log[10] * " location density"), metric_limits("02_log_density"),
  "Log location density is not occupancy. The last comparable modelled slice is 5 Ma; the record-conditioned 0 Ma endpoint is audited separately.",
  breaks = c(-12, -9, -6, -3, 0))
save_plate(fig_location, "Fig11_lineage_locations_six_times")

# A kernel prediction, not inferred source-to-target posterior movement.
select_index <- function(index, time_ma, metric, focal_tip = NULL) {
  rows <- index[index$time_ma == time_ma & index$metric == metric, , drop = FALSE]
  if (!is.null(focal_tip)) rows <- rows[rows$focal_tip == focal_tip, , drop = FALSE]
  if (nrow(rows) != 1L || !file.exists(rows$geotiff[[1L]])) {
    stop("Missing unique indexed map for ", metric, " at ", time_ma, " Ma")
  }
  rows
}
raster_values <- function(path) {
  tab <- terra::as.data.frame(terra::rast(path), xy = TRUE, na.rm = FALSE)
  names(tab)[1:3] <- c("lon", "lat", "value")
  key <- paste(tab$lon, tab$lat)
  z <- tab$value[match(paste(grid$lon, grid$lat), key)]
  if (length(z) != nrow(grid)) stop("Raster coordinate grid mismatch: ", path)
  z
}
movement <- do.call(rbind, lapply(times, function(time_ma) {
  row <- select_index(movement_index, time_ma, "potential_land_corridors")
  raw <- raster_values(row$geotiff[[1L]])
  tab <- map_frame(time_ma, log10(pmax(raw, 1e-10)))
  tab$value[tab$land & (!is.finite(raw) | raw <= 0)] <- NA_real_
  tab
}))
movement_max <- unique(movement_index$scale_max[movement_index$metric == "potential_land_corridors"])
if (length(movement_max) != 1L) stop("Corridor scales differ across time.")
fig_movement <- numeric_plate(
  movement, "Potential local movement changes with the palaeogeographic stage",
  "Plant200 | forward-predictive terrain kernel | per 1 Myr, not posterior-used routes",
  expression(log[10] * " predictive mass / Myr"),
  c(-10, log10(movement_max)),
  "Potential local-link throughput is not a realized route. Panels end at 5 Ma; the record-conditioned 0 Ma endpoint is audited separately.",
  breaks = c(-10, -8, -6, -4, -2))
save_plate(fig_movement, "Fig12_predictive_movement_six_times")

# exp(Shannon) of a terminal-calibrated lineage-support mixture.
diversity <- do.call(rbind, lapply(times, function(time_ma) {
  key <- as.character(time_ma)
  z <- field$support_mixture_effective_lineages[[key]]
  z[field$terminal_calibrated_lineage_support_intensity[[key]] <=
      .Machine$double.eps] <- NA_real_
  tab <- map_frame(time_ma, z, mark_zero_as_unavailable = TRUE)
  tab
}))
fig_diversity <- numeric_plate(
  diversity, "Effective diversity of the lineage-support mixture",
  "Plant200 | terminal-calibrated support fields | one scale across the entire 66-time atlas",
  "Effective support\nlineages", metric_limits("07_support_effective"),
  "Support-mixture entropy is not species richness. Panels end at 5 Ma; grey land has undefined or sub-precision support, not confirmed absence.",
  breaks = c(1, 10, 20, 30, 40))
save_plate(fig_diversity, "Fig13_support_mixture_diversity_six_times")

# Candidate refugia: 1 = opportunity-only; 2 = same candidate with modelled
# geographic-location mass above the declared threshold. Neither is survival.
refugia <- do.call(rbind, lapply(refuge_times, function(time_ma) {
  row <- select_index(refuge_index, time_ma, "candidate_refugia", focal)
  if (row$status[[1L]] != "evaluated") stop("Selected refuge boundary slice is not evaluable.")
  raw <- raster_values(row$geotiff[[1L]])
  tab <- grid
  tab$time_ma <- time_ma
  tab$panel <- factor(paste0(time_ma, " Ma | ", row$active_ancestor[[1L]]),
    levels = paste0(refuge_times, " Ma | ",
      vapply(refuge_times, function(t) select_index(refuge_index, t,
                      "candidate_refugia", focal)$active_ancestor[[1L]], character(1))))
  tab$land <- land_for(time_ma)
  tab$class <- factor(as.integer(raw), levels = 1:2,
    labels = c("Ecological candidate", "Candidate + modelled location mass"))
  tab$class[!tab$land] <- NA
  tab
}))
fig_refugia <- ggplot2::ggplot(refugia, ggplot2::aes(lon, lat)) +
  ggplot2::geom_tile(data = refugia[refugia$land, , drop = FALSE],
                     fill = "#E8ECEA", width = 4, height = 4) +
  ggplot2::geom_tile(ggplot2::aes(fill = class), width = 4, height = 4) +
  ggplot2::scale_fill_manual(values = c("Ecological candidate" = "#19866D",
    "Candidate + modelled location mass" = "#D58A28"),
    drop = FALSE, na.translate = FALSE, na.value = NA, name = "Evidence class") +
  ggplot2::facet_wrap(~panel, ncol = 3, drop = FALSE) +
  ggplot2::coord_fixed(expand = FALSE) +
  ggplot2::scale_x_continuous(breaks = c(-120, 0, 120), limits = c(-180, 180)) +
  ggplot2::scale_y_continuous(breaks = c(-60, 0, 60), limits = c(-90, 90)) +
  ggplot2::labs(title = "Candidate ecological refugia along one ancestral path",
    subtitle = "Path ending in Camellia sinensis | strip labels identify the active ancestral lineage",
    x = "Palaeolongitude", y = "Palaeolatitude",
    caption = "Patch continuity through modelled contractions, not confirmed population persistence. Ocean is white; other land is grey.") +
  base_theme
save_plate(fig_refugia, "Fig14_refugia_candidates_six_times")

# The abrupt 5-to-0 Ma support drop is caused by exact terminal-record
# conditioning, so make the endpoint-data footprint explicit.
endpoint_times <- c(5, 0)
endpoint <- do.call(rbind, lapply(endpoint_times, function(time_ma) {
  tab <- grid
  tab$land <- land_for(time_ma)
  tab$time_ma <- time_ma
  tab$panel <- factor(if (time_ma == 5) "5 Ma | modelled location field"
                      else "0 Ma | recorded terminal bins",
                      levels = c("5 Ma | modelled location field",
                                 "0 Ma | recorded terminal bins"))
  positive <- field$location_density[[as.character(time_ma)]] > 0
  tab$signal <- ifelse(positive & tab$land,
    if (time_ma == 5) "Positive modelled location" else "Recorded endpoint support", NA_character_)
  tab
}))
endpoint_counts <- do.call(rbind, lapply(endpoint_times, function(time_ma) {
  land <- land_for(time_ma)
  density <- field$location_density[[as.character(time_ma)]]
  data.frame(time_ma = time_ma, land_cells = sum(land),
             positive_location_cells = sum(land & density > 0),
             zero_location_cells = sum(land & density == 0))
}))
utils::write.csv(endpoint_counts, file.path(out, "terminal_endpoint_coverage_audit.csv"),
                 row.names = FALSE)
fig_endpoint <- ggplot2::ggplot(endpoint, ggplot2::aes(lon, lat)) +
  ggplot2::geom_tile(data = endpoint[endpoint$land, , drop = FALSE],
                     fill = "#E8ECEA", width = 4, height = 4) +
  ggplot2::geom_tile(ggplot2::aes(fill = signal), width = 4, height = 4) +
  ggplot2::scale_fill_manual(values = c("Positive modelled location" = "#13846F",
                                   "Recorded endpoint support" = "#D28C24"),
                             na.value = NA, na.translate = FALSE,
                             name = "Evidence class") +
  ggplot2::facet_wrap(~panel, ncol = 2, drop = FALSE) +
  ggplot2::coord_fixed(expand = FALSE) +
  ggplot2::scale_x_continuous(breaks = c(-120, 0, 120), limits = c(-180, 180)) +
  ggplot2::scale_y_continuous(breaks = c(-60, 0, 60), limits = c(-90, 90)) +
  ggplot2::labs(title = "The record-only terminal constraint creates a spatial discontinuity",
    subtitle = sprintf("5 Ma: %s/%s land cells with location mass; 0 Ma: %s/%s recorded endpoint cells",
      endpoint_counts$positive_location_cells[[1L]], endpoint_counts$land_cells[[1L]],
      endpoint_counts$positive_location_cells[[2L]], endpoint_counts$land_cells[[2L]]),
    x = "Palaeolongitude", y = "Palaeolatitude",
    caption = "Grey land at 0 Ma is unrepresented by terminal records, not demonstrated absence or an inferred extinction.") +
  base_theme
save_plate(fig_endpoint, "Fig15_terminal_endpoint_coverage_audit")

summaries <- do.call(rbind, lapply(times, function(time_ma) {
  key <- as.character(time_ma)
  land <- land_for(time_ma)
  den <- as.numeric(field$location_density[[key]])
  eff <- as.numeric(field$support_mixture_effective_lineages[[key]])
  eff[field$terminal_calibrated_lineage_support_intensity[[key]] <=
      .Machine$double.eps] <- NA_real_
  data.frame(time_ma = time_ma,
    endpoint_type = if (time_ma == 0) "exact_record_conditioning" else "modelled_location",
    active_lineages = ncol(field$lineage_location_density[[key]]),
    land_cells = sum(land),
    max_location_density = max(den[land], na.rm = TRUE),
    max_effective_support_lineages = max(eff[land], na.rm = TRUE),
    n_land_with_defined_support_diversity = sum(is.finite(eff[land])))
}))
utils::write.csv(summaries, file.path(out, "time_map_source_summary.csv"), row.names = FALSE)

new_manifest <- data.frame(
  figure_id = c("Fig11_lineage_locations_six_times",
                "Fig12_predictive_movement_six_times",
                "Fig13_support_mixture_diversity_six_times",
                "Fig14_refugia_candidates_six_times",
                "Fig15_terminal_endpoint_coverage_audit"),
  source = c(field_path, movement_index_path, field_path, refuge_index_path,
             field_path),
  evidence_class = c("archived_1.0.3_full200_location_field",
                     "real_inputs_forward_predictive_not_posterior",
                     "archived_1.0.3_support_mixture_with_1.0.4_precision_mask",
                     "corrected_plant200_environmental_candidate",
                     "archived_1.0.3_exact_terminal_record_conditioning_audit"),
  interpretation = c("Tree-conditioned lineage-location density, log display, not occupancy",
                     "Kernel-predicted local-link throughput, not historical posterior routes",
                     "Effective number of lineage-support fields; sub-precision cells masked; not historical richness",
                     "Patch-history candidate refugia, not verified population persistence",
                     "5-to-0 Ma coverage collapse comes from record-only terminal conditioning, not biological extinction"),
  stringsAsFactors = FALSE)
manifest_path <- file.path(out, "figure_source_manifest.csv")
old_manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
old_manifest <- old_manifest[!old_manifest$figure_id %in% new_manifest$figure_id, , drop = FALSE]
utils::write.csv(rbind(old_manifest, new_manifest), manifest_path, row.names = FALSE)

stopifnot(all(file.exists(file.path(fig_dir, paste0(new_manifest$figure_id, ".png")))),
          all(file.exists(file.path(fig_dir, paste0(new_manifest$figure_id, ".svg")))),
          all(summaries$land_cells > 0),
          all(summaries$active_lineages > 0),
          all(summaries$max_effective_support_lineages <=
              summaries$active_lineages + 1e-9))
message("Saved five time-resolved and endpoint-audit figure plates: ", fig_dir)
