#!/usr/bin/env Rscript

# Rebuild only the manuscript panels changed in the six-figure revision.
root <- Sys.getenv("HMSCEE_ROOT", "C:/Users/Google/Documents/HMSC-HIST")
source_dir <- file.path(root, "HmscEcoEvo", "docs", "manuscript",
                        "MEE_full_manuscript_20260924", "source_data")
out <- file.path(root, "outputs", "HmscEcoEvo",
                 "manuscript_results_six_process_20260927")
fig_dir <- file.path(out, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(requireNamespace("ggplot2", quietly = TRUE),
          requireNamespace("patchwork", quietly = TRUE),
          requireNamespace("ragg", quietly = TRUE),
          requireNamespace("svglite", quietly = TRUE))
library(ggplot2)
library(patchwork)

save_panel <- function(plot, stem, width, height) {
  ggsave(file.path(fig_dir, paste0(stem, ".png")), plot,
         width = width, height = height, units = "in", dpi = 400,
         device = ragg::agg_png, bg = "white")
  ggsave(file.path(fig_dir, paste0(stem, ".svg")), plot,
         width = width, height = height, units = "in",
         device = svglite::svglite, bg = "white")
}

workflow_archive <- file.path(fig_dir, "Fig01_conceptual_workflow_user.png")
workflow_attachment <- "C:/Users/Google/AppData/Local/Temp/codex-clipboard-953b7146-e3f9-4b88-821a-ecf1d71f2740.png"
if (!file.exists(workflow_archive)) {
  stopifnot(file.exists(workflow_attachment))
  stopifnot(file.copy(workflow_attachment, workflow_archive))
}

e <- readRDS(file.path(source_dir, "figure_inputs.rds"))
design <- e$design
temporal <- e$temporal
grid <- data.frame(cell = sprintf("g4_%04d", 1:4050),
                   lon = rep(seq(-178, 178, 4), 45),
                   lat = rep(seq(-88, 88, 4), each = 90))
carrier_root <- file.path(root, "HmscEcoEvo", "derived_inputs",
                          "case05_v2_paleomap_h3_carrier_325_0Ma_1deg",
                          "movement_cache_arias_topographic_rate0p4")
carrier_index <- read.csv(file.path(carrier_root,
                                    "case05_v2_movement_cache_index.csv"))
carrier_file <- carrier_index$file[match(0, carrier_index$time_ma)]
stopifnot(length(carrier_file) == 1, !is.na(carrier_file), file.exists(carrier_file))
carrier <- readRDS(carrier_file)$carrier
carrier <- carrier[is.finite(carrier$paleo_lon) & is.finite(carrier$paleo_lat) &
                   !is.na(carrier$active_land) & carrier$active_land > 0 &
                   is.finite(carrier$land_area_km2) & carrier$land_area_km2 > 0, ]
bin <- pmin(89L, pmax(0L, floor((carrier$paleo_lon + 180) / 4))) +
  90L * pmin(44L, pmax(0L, floor((carrier$paleo_lat + 90) / 4))) + 1L
grid$land <- seq_len(nrow(grid)) %in% unique(bin)

split <- read.csv(file.path(e$v5, "01_modern_observation_design",
                                 "modern_4deg_spatial_holdout_split.csv"),
                  stringsAsFactors = FALSE)
stopifnot(nrow(split) == 554, all(split$split %in% c("training", "holdout")))
prev <- data.frame(species = colnames(design$Y),
                   occupied_cells = colSums(design$Y > 0))
stopifnot(nrow(prev) == 200, nrow(temporal) > 0)

base_theme <- theme_bw(base_size = 12, base_family = "Arial") +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        plot.tag = element_text(face = "bold", size = 15))
p_map <- ggplot(grid, aes(lon, lat)) +
  geom_tile(data = grid[grid$land, , drop = FALSE], fill = "#DBE0E0") +
  geom_contour(aes(z = as.numeric(land)), breaks = .5,
               colour = "#69777A", linewidth = .18) +
  geom_point(data = split, aes(colour = split), size = .95) +
  scale_colour_manual(values = c(training = "#17685B", holdout = "#B5472E")) +
  coord_fixed(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  labs(title = "Modern four-degree model cells", x = "Longitude",
       y = "Latitude", colour = NULL) +
  base_theme + theme(legend.position = "bottom")
p_tree <- ggplot(temporal, aes(time_ma, active_lineages)) +
  geom_step(colour = "#573F76", linewidth = .9) +
  scale_x_reverse() +
  labs(title = "Lineages in the sampled tree",
       x = "Time before present (Ma)", y = "Active lineages") + base_theme
p_records <- ggplot(prev, aes(occupied_cells)) +
  geom_histogram(binwidth = 3, fill = "#557890", colour = "white") +
  labs(title = "Uneven modern recording",
       x = "Recorded four-degree cells per species", y = "Species") + base_theme
figure2 <- (p_map / (p_tree + p_records)) +
  plot_layout(heights = c(1.35, 1)) + plot_annotation(tag_levels = "a")
save_panel(figure2, "Fig02_input_coverage_three_panels", 10.5, 7.5)

noanalog <- read.csv(file.path(root, "HmscEcoEvo", "outputs",
                                "six_process_case06_full_20260927", "data",
                                "noanalog_diagnostic_by_time.csv"))
stopifnot(all(c("time_ma", "noanalog_area_fraction") %in% names(noanalog)))
p_extrap <- ggplot(noanalog, aes(time_ma, noanalog_area_fraction)) +
  geom_area(fill = "#D58323", alpha = .22) +
  geom_line(linewidth = 1.2, colour = "#D58323") +
  scale_x_reverse(limits = c(325, 0), breaks = c(325, 300, 200, 100, 0)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(title = "Land outside the modern training envelope",
       subtitle = "Environmental extrapolation is a diagnostic, not absence",
       x = "Time before present (Ma)", y = "Area-weighted fraction") +
  theme_minimal(base_size = 14, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())
save_panel(p_extrap, "Supp_environmental_extrapolation", 8.5, 4.5)

cat("Revised figures saved in", fig_dir, "\n")
