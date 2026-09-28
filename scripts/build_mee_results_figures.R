#!/usr/bin/env Rscript

# Manuscript figures from archived Plant200 and current six-process results.
root <- Sys.getenv("HMSCEE_ROOT", "C:/Users/Google/Documents/HMSC-HIST")
archived <- file.path(root, "outputs", "HmscEcoEvo")
current <- file.path(root, "HmscEcoEvo", "outputs")
out <- file.path(archived, "manuscript_results_six_process_20260927")
fig_dir <- file.path(out, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(requireNamespace("ggplot2", quietly = TRUE),
          requireNamespace("patchwork", quietly = TRUE),
          requireNamespace("ragg", quietly = TRUE),
          requireNamespace("svglite", quietly = TRUE))
library(ggplot2)
library(patchwork)

pal <- c(blue = "#2766A3", teal = "#16836F", orange = "#D58323",
         red = "#B9473A", ink = "#243B50", grey = "#728493")
theme_set(theme_minimal(base_size = 13, base_family = "Arial") +
            theme(panel.grid.minor = element_blank(),
                  panel.grid.major.x = element_line(colour = "#E5E9EC"),
                  panel.grid.major.y = element_line(colour = "#E5E9EC"),
                  plot.title = element_text(face = "bold", colour = pal[["ink"]], size = 16),
                  plot.subtitle = element_text(colour = pal[["grey"]], size = 11),
                  axis.title = element_text(colour = pal[["ink"]]),
                  legend.position = "bottom"))

read_checked <- function(path) {
  if (!file.exists(path)) stop("Missing source data: ", path)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
save_figure <- function(plot, name, width = 11, height = 7) {
  png <- file.path(fig_dir, paste0(name, ".png"))
  svg <- file.path(fig_dir, paste0(name, ".svg"))
  ggplot2::ggsave(png, plot, width = width, height = height, dpi = 400,
                  device = ragg::agg_png, bg = "white")
  ggplot2::ggsave(svg, plot, width = width, height = height,
                  device = svglite::svglite, bg = "white")
  c(png = png, svg = svg)
}
manifest <- list()
register <- function(id, source, evidence, interpretation) {
  manifest[[length(manifest) + 1L]] <<- data.frame(
    figure_id = id, source = source, evidence_class = evidence,
    interpretation = interpretation, stringsAsFactors = FALSE)
}

case06 <- file.path(current, "six_process_case06_full_20260927")
case07 <- file.path(current, "six_process_case07_full_20260927")
case05 <- file.path(archived, "case05_v8_four_scheme_comparison_20260926")
case05_land <- file.path(archived, "case05_v8_full200_spherical_land_with_ldd_20260926")

# Figure 1: tree-conditioned domain and climatic extrapolation are distinct.
movement_path <- file.path(case07, "07_dispersal_source_corridor_maps",
                          "plant200_movement_diagnostics.csv")
noanalog_path <- file.path(case06, "data", "noanalog_diagnostic_by_time.csv")
movement <- read_checked(movement_path)
noanalog <- read_checked(noanalog_path)
p_lineage <- ggplot(movement, aes(time_ma, active_lineages)) +
  geom_line(linewidth = 1.1, colour = pal[["blue"]]) +
  geom_point(data = subset(movement, time_ma %in% c(0, 100, 200, 300, 325)),
             size = 2.4, colour = pal[["blue"]]) +
  scale_x_reverse(limits = c(325, 0), breaks = c(325, 300, 200, 100, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
  labs(title = "A  Active sampled lineages", y = "Number of branches", x = NULL)
p_noanalog <- ggplot(noanalog, aes(time_ma, noanalog_area_fraction)) +
  geom_area(fill = pal[["orange"]], alpha = 0.22) +
  geom_line(linewidth = 1.1, colour = pal[["orange"]]) +
  scale_x_reverse(limits = c(325, 0), breaks = c(325, 300, 200, 100, 0)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(title = "B  Land outside the modern training envelope",
       y = "Area-weighted fraction", x = "Time before present (Ma)")
save_figure((p_lineage / p_noanalog) +
              plot_annotation(title = "The biological and environmental domains change differently",
                subtitle = "Plant200, 66 native time slices; environmental extrapolation is a diagnostic, not absence"),
            "Fig01_lineage_domain_and_extrapolation", height = 8.2)
register("Fig01_lineage_domain_and_extrapolation", paste(movement_path, noanalog_path, sep = " | "),
         "real_input_descriptive", "Tree-defined active branches and external no-analog diagnostic")

# Figure 2: alternative evolutionary models and boundary-limited fits.
fit_path <- file.path(case06, "data", "model_comparison_by_axis_draw.csv")
fit <- read_checked(fit_path)
fit$delta_aicc <- as.numeric(fit$delta_aicc)
fit$axis <- gsub("_", " ", fit$axis)
fit$boundary <- fit$boundary_status != "none"
summ <- aggregate(delta_aicc ~ axis + model, fit, stats::median)
bound <- aggregate(boundary ~ axis + model, fit, mean)
summ <- merge(summ, bound, by = c("axis", "model"))
summ$model <- factor(summ$model, levels = c("BM", "OU", "EB"))
p_fit <- ggplot(summ, aes(model, reorder(axis, delta_aicc, FUN = max))) +
  geom_tile(aes(fill = pmin(delta_aicc, 250)), colour = "white", linewidth = 0.9) +
  geom_text(aes(label = sprintf("%.0f%s", delta_aicc, ifelse(boundary == 1, "*", ""))),
            size = 3.3, colour = pal[["ink"]]) +
  scale_fill_gradient(low = "#E1F0EC", high = "#D9813B", name = "Median delta AICc\n(display clipped at 250)") +
  labs(title = "Ancestral-response model fits",
       subtitle = "Five HMSC draws for each of eight axes; * all five estimates at a parameter bound",
       x = "Evolution model", y = "Environmental-response axis") +
  theme(panel.grid = element_blank(), legend.position = "right")
save_figure(p_fit, "Fig02_evolution_model_comparison", width = 10, height = 6.5)
register("Fig02_evolution_model_comparison", fit_path, "fitted_model_comparison",
         "AICc rankings are conditional on the selected HMSC draws; boundary fits limit inference")

# Figure 3: response-model choice changes ancient support, not the modern endpoint.
diff_path <- file.path(case06, "data", "model_difference_summary.csv")
diff <- read_checked(diff_path)
diff$comparison <- factor(diff$comparison, levels = c("OU - BM", "EB - BM"))
diff$metric <- factor(diff$metric,
                      levels = c("focal_path_support", "mean_active_lineage_support"),
                      labels = c("Focal ancestral path", "Mean active-lineage support"))
p_diff <- ggplot(diff, aes(time_ma, area_weighted_mean_abs_difference,
                            colour = comparison, group = comparison)) +
  geom_line(linewidth = 1.1) + geom_point(size = 2.1) +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  scale_x_reverse(breaks = c(300, 200, 100, 0)) +
  scale_colour_manual(values = c("OU - BM" = pal[["teal"]], "EB - BM" = pal[["orange"]])) +
  labs(title = "Model choice changes historical environmental support",
       subtitle = "Area-weighted absolute differences over land; the shared modern endpoint is constrained",
       x = "Time before present (Ma)", y = "Mean absolute support difference", colour = NULL)
save_figure(p_diff, "Fig03_evolution_support_sensitivity", width = 11, height = 5.4)
register("Fig03_evolution_support_sensitivity", diff_path, "derived_environmental_support",
         "Model sensitivity, not evidence of realized historical occupancy")

# Figure 4: four archived 200-tip movement settings on a matched 4-degree grid.
movement_diff_path <- file.path(case05, "four_scheme_difference_by_time.csv")
md <- read_checked(movement_diff_path)
md$comparison <- paste(md$numerator, "versus", md$denominator)
lab <- c("spherical_land versus spherical_no_ldd" = "Rare land jumps added",
         "terrain_resistance versus spherical_land" = "Terrain resistance",
         "particle_topography versus terrain_resistance" = "Finite-propagule particle",
         "particle_topography versus spherical_land" = "Particle versus spherical")
md$comparison <- factor(md$comparison, levels = names(lab), labels = unname(lab))
p_movement <- ggplot(md, aes(time_ma, total_variation_per_lineage,
                              colour = comparison, group = comparison)) +
  geom_line(linewidth = 1.05) +
  scale_x_reverse(limits = c(325, 0), breaks = c(325, 300, 200, 100, 0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values = c(pal[["red"]], pal[["blue"]],
                                 pal[["teal"]], pal[["orange"]])) +
  labs(title = "Movement assumptions separate ancient location fields",
       subtitle = "Total variation per active lineage; all four schemes share the conditioned modern endpoint",
       x = "Time before present (Ma)", y = "Total variation distance", colour = NULL)
save_figure(p_movement, "Fig04_movement_scheme_sensitivity", width = 11, height = 6)
register("Fig04_movement_scheme_sensitivity", movement_diff_path,
         "archived_1.0.3_full200_scenario_comparison",
         "Differences among conditional location distributions, not occupancy or empirical movement rates")

# Figure 5: candidate ecological refugia depend on declared contraction threshold.
patch_path <- file.path(case07, "01_real_plant200", "focal_path_summary.csv")
patch <- read_checked(patch_path)
focal <- gsub("_", " ", patch$focal_tip)
rs <- do.call(rbind, lapply(seq_len(nrow(patch)), function(i) data.frame(
  focal_tip = focal[[i]], threshold = c("10%", "25%"),
  contraction_intervals = c(patch$contraction_intervals_10pct[[i]],
                             patch$contraction_intervals_25pct[[i]]),
  candidate_patches = c(patch$candidate_refugia_10pct[[i]],
                        patch$candidate_refugia_25pct[[i]]))))
long <- rbind(data.frame(rs[c("focal_tip", "threshold")],
                         metric = "Contraction intervals", count = rs$contraction_intervals),
              data.frame(rs[c("focal_tip", "threshold")],
                         metric = "Candidate patches", count = rs$candidate_patches))
long$threshold <- factor(long$threshold, levels = c("10%", "25%"))
p_refugia <- ggplot(long, aes(focal_tip, count, fill = threshold)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("10%" = pal[["blue"]], "25%" = pal[["orange"]])) +
  labs(title = "Candidate refugia depend on the contraction rule",
       subtitle = "Three focal ancestral paths; counts describe supported patch histories, not confirmed populations",
       x = NULL, y = "Count", fill = "Area contraction") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_figure(p_refugia, "Fig05_candidate_refugia_thresholds", width = 11, height = 6)
register("Fig05_candidate_refugia_thresholds", patch_path,
         "real_input_descriptive_candidate", "Threshold-sensitive ecological candidates, not demographic persistence")

# Curated existing map and hmscHist plates. Copying does not alter pixels.
curated <- list(
  Fig06_modern_holdout = file.path(archived, "case05_data_driven_framework",
                                     "modern_spatial_holdout_endpoint.png"),
  Fig07_ancestral_response_paths = file.path(case06, "figures",
                                              "02_focal_ancestral_beta_trajectories.png"),
  Fig08_ancient_environmental_support = file.path(case06, "figures",
                                                   "06_mean_active_lineage_support_four_times.png"),
  Fig09_terrestrial_location_100Ma = file.path(case05_land,
    "05_shared_scale_maps", "01_location_density", "01_location_density_100Ma.png"),
  Fig10_noanalog_100Ma = file.path(case05, "10_patch_history_results_corrected",
                                  "maps", "noanalog_training_envelope_100Ma.png"),
  Fig11_candidate_refugia_history = file.path(case07, "04_figures",
                                             "02_real_candidate_refugia.png"),
  Fig12_potential_land_corridors_100Ma = file.path(case07,
    "07_dispersal_source_corridor_maps", "potential_land_corridors",
    "potential_land_corridors_100Ma.png"))
for (id in names(curated)) {
  from <- curated[[id]]
  if (!file.exists(from)) stop("Missing image source: ", from)
  if (!file.copy(from, file.path(fig_dir, paste0(id, ".png")), overwrite = TRUE))
    stop("Could not copy image source: ", from)
  kind <- if (grepl("Fig09|Fig10", id)) "archived_1.0.3_full200_diagnostic" else
    if (grepl("Fig12", id)) "forward_predictive_not_posterior" else
      "real_input_diagnostic"
  register(id, from, kind, "Curated unchanged PNG; see Results legend for interpretation")
}
hist_dir <- file.path(root, "HmscEcoEvo", "inst", "extdata", "ecoevo_cases",
                      "outputs", "case07_validation_dashboard")
for (name in c("niche_summary", "phylo_signal", "evo_transition",
               "trait_mediation", "gamma_evolution", "population_evolution",
               "validation_dashboard")) {
  id <- paste0("Supp_hmscHist_", name)
  from <- file.path(hist_dir, paste0(name, ".png"))
  if (!file.exists(from)) stop("Missing hmscHist demonstration plate: ", from)
  if (!file.copy(from, file.path(fig_dir, paste0(id, ".png")), overwrite = TRUE))
    stop("Could not copy hmscHist plate: ", from)
  register(id, from, "bundled_synthetic_demonstration",
           "Interface demonstration only; not a Plant200 fitted-history result")
}

manifest <- do.call(rbind, manifest)
utils::write.csv(manifest, file.path(out, "figure_source_manifest.csv"), row.names = FALSE)
cat("Figures:", nrow(manifest), "\nOutput:", out, "\n")
