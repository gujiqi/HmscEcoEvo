#!/usr/bin/env Rscript

args <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)][1])
script_dir <- if (!is.na(file_arg) && nzchar(file_arg)) dirname(normalizePath(file_arg)) else getwd()
pkg_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = TRUE)

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(pkg_root, quiet = TRUE)
} else {
  library(HmscEcoEvo)
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required to save deep-time workflow plots.", call. = FALSE)
}
if (!requireNamespace("ape", quietly = TRUE)) {
  stop("Package 'ape' is required to draw dated phylogenies.", call. = FALSE)
}

pal <- c(
  navy = "#163A5F",
  blue = "#2F6DB5",
  teal = "#2D7F7B",
  green = "#2D7F5E",
  gold = "#D8A03D",
  orange = "#C76D2A",
  red = "#B23A3A",
  purple = "#6B4FA3",
  grey = "#6E7781"
)

theme_case <- function() {
  ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.18, colour = "grey88"),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey70"),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = 10),
      legend.key.height = grid::unit(3.5, "mm"),
      legend.key.width = grid::unit(4.5, "mm")
    )
}

out_root <- file.path(script_dir, "outputs")
if (dir.exists(out_root)) unlink(out_root, recursive = TRUE, force = TRUE)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

write_csv <- function(x, path) {
  utils::write.csv(as.data.frame(x), path, row.names = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

save_rds <- function(x, path) {
  saveRDS(x, path, compress = "xz")
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

save_plot_bundle <- function(plot, base_file, width = 8.2, height = 5.2, dpi = 320) {
  png_file <- paste0(base_file, ".png")
  pdf_file <- paste0(base_file, ".pdf")
  ggplot2::ggsave(png_file, plot = plot, width = width, height = height,
                  dpi = dpi, bg = "white")
  ggplot2::ggsave(pdf_file, plot = plot, width = width, height = height,
                  bg = "white")
  c(normalizePath(png_file, winslash = "/", mustWork = FALSE),
    normalizePath(pdf_file, winslash = "/", mustWork = FALSE))
}

save_tree_bundle <- function(tree, species_origin, base_file, title) {
  png_file <- paste0(base_file, ".png")
  pdf_file <- paste0(base_file, ".pdf")
  draw <- function() {
    old <- par(mar = c(4, 1, 3, 1))
    on.exit(par(old), add = TRUE)
    ape::plot.phylo(tree, direction = "rightwards", cex = 0.65,
                    label.offset = 4, no.margin = FALSE)
    title(main = title, cex.main = 0.95)
    axis(1, at = seq(0, 540, by = 90), labels = 540 - seq(0, 540, by = 90),
         cex.axis = 0.7)
    mtext("Time before present (Ma)", side = 1, line = 2.4, cex = 0.75)
  }
  grDevices::png(png_file, width = 8.2, height = 5.2, units = "in", res = 320, bg = "white")
  draw()
  grDevices::dev.off()
  grDevices::pdf(pdf_file, width = 8.2, height = 5.2)
  draw()
  grDevices::dev.off()
  c(normalizePath(png_file, winslash = "/", mustWork = FALSE),
    normalizePath(pdf_file, winslash = "/", mustWork = FALSE))
}

long_matrix <- function(M, row_name = "row", col_name = "column", value_name = "value") {
  d <- as.data.frame(as.table(M), stringsAsFactors = FALSE)
  names(d) <- c(row_name, col_name, value_name)
  d
}

summarise_richness_time <- function(richness) {
  rows <- lapply(split(richness$expected_richness, richness$time_ma), function(z) {
    data.frame(mean = mean(z, na.rm = TRUE),
               median = stats::median(z, na.rm = TRUE),
               q95 = as.numeric(stats::quantile(z, 0.95, na.rm = TRUE)),
               max = max(z, na.rm = TRUE))
  })
  out <- do.call(rbind, rows)
  out$time_ma <- as.numeric(rownames(out))
  out[, c("time_ma", "mean", "median", "q95", "max")]
}

summarise_predictors <- function(paleo_grid, variables) {
  rows <- lapply(variables, function(v) {
    pieces <- lapply(split(paleo_grid[[v]], paleo_grid$time_ma), function(x) {
      data.frame(mean = mean(x, na.rm = TRUE),
                 lwr = as.numeric(stats::quantile(x, 0.05, na.rm = TRUE)),
                 upr = as.numeric(stats::quantile(x, 0.95, na.rm = TRUE)))
    })
    z <- do.call(rbind, pieces)
    z$time_ma <- as.numeric(rownames(z))
    data.frame(variable = v, time_ma = z$time_ma,
               mean = z$mean, lwr = z$lwr, upr = z$upr)
  })
  do.call(rbind, rows)
}

summarise_turnover <- function(turnover) {
  stats::aggregate(turnover ~ time_to_ma, turnover, mean, na.rm = TRUE)
}

workflow_steps <- function(case, x) {
  rows <- data.frame(
    case = case,
    step = c(
      "load authoritative 4-degree Phanerozoic RDS",
      "index and check full 540-0 Ma timecube",
      "construct dated 540 Ma phylogeny and species origins",
      "lock modern predictor recipe",
      "check leakage-safe exogenous predictors",
      "project Beta through all 109 time slices",
      "apply phylogenetic existence mask",
      "apply extrapolation/accessibility/dynamic weights",
      "calculate richness",
      "calculate turnover",
      "calculate refugia",
      "validate modern prediction",
      "validate fossil or track evidence",
      "summarize BioGeoBEARS-style accessibility/events",
      "write source-backed figures and tables"
    ),
    status = "ready",
    stringsAsFactors = FALSE
  )
  rows$n_rows <- c(
    length(x$times),
    if (!is.null(x$timecube_index)) nrow(x$timecube_index) else length(x$times),
    length(x$species_origin),
    length(x$recipe$variables),
    length(x$variables),
    nrow(x$suitability),
    nrow(x$projection),
    nrow(x$projection),
    nrow(x$richness),
    nrow(x$turnover),
    nrow(x$refugia),
    nrow(x$modern_validation$evaluation),
    if (!is.null(x$fossil_validation)) nrow(x$fossil_validation) else nrow(x$projection),
    nrow(x$bgb$compare),
    NA_integer_
  )
  rows
}

plot_workflow_steps <- function(steps, title) {
  steps$step <- factor(steps$step, levels = rev(steps$step))
  steps$plot_rows <- ifelse(is.na(steps$n_rows), 1, pmax(steps$n_rows, 1))
  ggplot2::ggplot(steps, ggplot2::aes(x = status, y = step)) +
    ggplot2::geom_point(ggplot2::aes(size = plot_rows), color = pal["green"], alpha = 0.85) +
    ggplot2::scale_size_continuous(range = c(2, 8), labels = scales::comma) +
    ggplot2::labs(x = NULL, y = NULL, size = "Rows/items", title = title) +
    theme_case()
}

plot_beta <- function(beta, title) {
  B <- beta[, setdiff(colnames(beta), "(Intercept)"), drop = FALSE]
  d <- long_matrix(B, "species", "predictor", "beta")
  ggplot2::ggplot(d, ggplot2::aes(x = predictor, y = species, fill = beta)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.18) +
    ggplot2::scale_fill_gradient2(low = pal["blue"], mid = "white", high = pal["red"],
                                  midpoint = 0) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Beta", title = title) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)) +
    theme_case()
}

plot_validation <- function(validation, title) {
  ev <- validation$evaluation
  d <- data.frame(metric = names(ev), value = as.numeric(ev[1, ]))
  fill_values <- stats::setNames(unname(c(pal["blue"], pal["green"], pal["orange"], pal["purple"])),
                                 d$metric)
  ggplot2::ggplot(d, ggplot2::aes(x = metric, y = value, fill = metric)) +
    ggplot2::geom_col(width = 0.65, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", value)), vjust = -0.25, size = 2.7) +
    ggplot2::scale_fill_manual(values = fill_values) +
    ggplot2::labs(x = NULL, y = "Value", title = title) +
    ggplot2::coord_cartesian(ylim = c(0, max(1, max(d$value, na.rm = TRUE) * 1.15))) +
    theme_case()
}

plot_bgb <- function(bgb, title) {
  events <- bgb$events
  event_time <- stats::aggregate(event_prob ~ time_ma + event_type, events, mean, na.rm = TRUE)
  ggplot2::ggplot() +
    ggplot2::geom_col(data = bgb$compare,
                      ggplot2::aes(x = reorder(model, AICc_weight), y = AICc_weight),
                      fill = pal["navy"], alpha = 0.85, width = 0.65) +
    ggplot2::geom_line(data = event_time,
                       ggplot2::aes(x = scales::rescale(time_ma, to = c(1, nrow(bgb$compare))),
                                    y = event_prob, color = event_type, group = event_type),
                       linewidth = 0.8, inherit.aes = FALSE) +
    ggplot2::geom_point(data = event_time,
                        ggplot2::aes(x = scales::rescale(time_ma, to = c(1, nrow(bgb$compare))),
                                     y = event_prob, color = event_type),
                        size = 1.3, inherit.aes = FALSE) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Model / rescaled event time", y = "AICc weight or event probability",
                  color = "Event", title = title) +
    theme_case()
}

run_case01 <- function() {
  case <- "case01_global540_20"
  x <- readRDS(file.path(script_dir, case, "case01_global540_20.rds"))
  out_dir <- file.path(out_root, "c01")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  files <- c(
    save_rds(x, file.path(out_dir, "case_input_and_all_results.rds")),
    save_rds(x$suitability, file.path(out_dir, "suitability_full_540_0Ma.rds")),
    save_rds(x$projection, file.path(out_dir, "projection_full_540_0Ma.rds")),
    save_rds(x$richness, file.path(out_dir, "richness_full_540_0Ma.rds")),
    save_rds(x$turnover, file.path(out_dir, "turnover_full_540_0Ma.rds")),
    save_rds(x$refugia, file.path(out_dir, "refugia_full_540_0Ma.rds"))
  )

  rich_time <- summarise_richness_time(x$richness)
  turnover_time <- summarise_turnover(x$turnover)
  predictor_time <- summarise_predictors(x$paleo_grid, x$variables)
  risk_time <- stats::aggregate(cbind(extrapolation_risk, extrapolation_weight) ~ time_ma,
                                x$extrapolation, mean, na.rm = TRUE)
  land_area <- stats::aggregate(land_area_km2 ~ time_ma, x$paleo_grid, sum, na.rm = TRUE)
  fossil_time <- stats::aggregate(cbind(predicted_probability, hit = as.numeric(hit)) ~ matched_time_ma,
                                  x$fossil_validation, mean, na.rm = TRUE)
  fossil_time_plot <- fossil_time[is.finite(fossil_time$predicted_probability), , drop = FALSE]
  steps <- workflow_steps(case, x)
  source_summary <- data.frame(
    source = x$source_level,
    n_time_slices = length(x$times),
    oldest_ma = max(x$times),
    youngest_ma = min(x$times),
    n_variables = length(x$variables),
    n_species = length(x$species),
    n_projection_rows = nrow(x$projection),
    n_richness_rows = nrow(x$richness),
    n_turnover_rows = nrow(x$turnover)
  )

  files <- c(files,
             write_csv(source_summary, file.path(out_dir, "source_and_size_summary.csv")),
             write_csv(data.frame(species = x$species, clade = x$clade[x$species],
                                  origin_ma = x$species_origin[x$species]),
                       file.path(out_dir, "species_origin.csv")),
             write_csv(x$beta, file.path(out_dir, "beta.csv")),
             write_csv(x$traits, file.path(out_dir, "traits.csv")),
             write_csv(x$timecube_index, file.path(out_dir, "timecube_index_full.csv")),
             write_csv(predictor_time, file.path(out_dir, "summary_predictors_by_time.csv")),
             write_csv(risk_time, file.path(out_dir, "summary_extrapolation_by_time.csv")),
             write_csv(land_area, file.path(out_dir, "summary_land_area_by_time.csv")),
             write_csv(rich_time, file.path(out_dir, "summary_richness_by_time.csv")),
             write_csv(turnover_time, file.path(out_dir, "summary_turnover_by_time.csv")),
             write_csv(x$refugia, file.path(out_dir, "refugia_summary.csv")),
             write_csv(x$fossils, file.path(out_dir, "simulated_fossils_full_axis.csv")),
             write_csv(x$fossil_validation, file.path(out_dir, "fossil_validation.csv")),
             write_csv(fossil_time, file.path(out_dir, "summary_fossil_validation_by_time.csv")),
             write_csv(x$modern_validation$evaluation, file.path(out_dir, "modern_validation_metrics.csv")),
             write_csv(x$bgb$compare, file.path(out_dir, "bgb_model_compare.csv")),
             write_csv(x$bgb$events, file.path(out_dir, "bgb_events.csv")),
             write_csv(steps, file.path(out_dir, "workflow_steps.csv")))

  p_source <- ggplot2::ggplot(land_area, ggplot2::aes(x = time_ma, y = land_area_km2 / 1e6)) +
    ggplot2::geom_area(fill = pal["teal"], alpha = 0.35) +
    ggplot2::geom_line(color = pal["navy"], linewidth = 0.7) +
    ggplot2::scale_x_reverse() +
    ggplot2::labs(x = "Time (Ma)", y = "Land area (million km2)",
                  title = "Authoritative 4-degree Phanerozoic source: full 540-0 Ma coverage") +
    theme_case()

  p_pred <- ggplot2::ggplot(predictor_time, ggplot2::aes(x = time_ma, y = mean)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr), fill = pal["blue"], alpha = 0.18) +
    ggplot2::geom_line(color = pal["navy"], linewidth = 0.55) +
    ggplot2::scale_x_reverse() +
    ggplot2::facet_wrap(~ variable, scales = "free_y", ncol = 2) +
    ggplot2::labs(x = "Time (Ma)", y = "Land-cell mean with 5-95% range",
                  title = "Paleoenvironmental predictors used by the HMSC projection") +
    theme_case()

  p_risk <- ggplot2::ggplot(risk_time, ggplot2::aes(x = time_ma)) +
    ggplot2::geom_line(ggplot2::aes(y = extrapolation_risk), color = pal["red"], linewidth = 0.75) +
    ggplot2::geom_line(ggplot2::aes(y = extrapolation_weight), color = pal["green"], linewidth = 0.75) +
    ggplot2::scale_x_reverse() +
    ggplot2::labs(x = "Time (Ma)", y = "Mean risk / weight",
                  title = "Extrapolation risk propagated through all time slices") +
    theme_case()

  p_rich <- ggplot2::ggplot(rich_time, ggplot2::aes(x = time_ma)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = median, ymax = q95), fill = pal["blue"], alpha = 0.18) +
    ggplot2::geom_line(ggplot2::aes(y = mean), color = pal["navy"], linewidth = 0.85) +
    ggplot2::geom_line(ggplot2::aes(y = max), color = pal["orange"], linewidth = 0.55, linetype = 2) +
    ggplot2::scale_x_reverse() +
    ggplot2::labs(x = "Time (Ma)", y = "Expected richness",
                  title = "Global expected richness across the full Phanerozoic") +
    theme_case()

  map_times <- c(540, 480, 420, 360, 300, 240, 180, 120, 60, 0)
  map_data <- x$richness[x$richness$time_ma %in% map_times, ]
  p_maps <- ggplot2::ggplot(map_data, ggplot2::aes(x = lon, y = lat, color = expected_richness)) +
    ggplot2::geom_point(size = 0.28) +
    ggplot2::facet_wrap(~ time_ma, labeller = ggplot2::label_both, ncol = 5) +
    ggplot2::scale_color_gradientn(colors = c(pal["navy"], "#F8F7F0", pal["red"])) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "Longitude", y = "Latitude", color = "Expected\nrichness",
                  title = "Global richness maps at representative Phanerozoic slices") +
    theme_case()

  turnover_time_plot <- turnover_time[is.finite(turnover_time$turnover), , drop = FALSE]
  p_turn <- ggplot2::ggplot(turnover_time_plot, ggplot2::aes(x = time_to_ma, y = turnover)) +
    ggplot2::geom_line(color = pal["purple"], linewidth = 0.75) +
    ggplot2::geom_point(color = pal["purple"], size = 1.1) +
    ggplot2::scale_x_reverse() +
    ggplot2::labs(x = "Interval end (Ma)", y = "Mean turnover",
                  title = "Community turnover between every adjacent 5 Myr slice") +
    theme_case()

  p_ref <- ggplot2::ggplot(x$refugia, ggplot2::aes(x = lon, y = lat, color = refugia_score)) +
    ggplot2::geom_point(size = 0.45) +
    ggplot2::scale_color_gradient2(low = pal["navy"], mid = "#F8F7F0", high = pal["red"],
                                   midpoint = 0) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "Longitude", y = "Latitude", color = "Refugia\nscore",
                  title = "Long-term refugia potential from 540-0 Ma richness stability") +
    theme_case()

  p_fossil <- ggplot2::ggplot(fossil_time_plot, ggplot2::aes(x = matched_time_ma)) +
    ggplot2::geom_line(ggplot2::aes(y = predicted_probability), color = pal["navy"], linewidth = 0.7) +
    ggplot2::geom_point(ggplot2::aes(y = predicted_probability, color = hit), size = 1.6) +
    ggplot2::scale_x_reverse() +
    ggplot2::scale_color_gradient(low = pal["red"], high = pal["green"]) +
    ggplot2::labs(x = "Matched time (Ma)", y = "Mean fossil-match probability / hit rate",
                  color = "Hit rate", title = "Simulated fossil validation across the full time axis") +
    theme_case()

  files <- c(files,
             save_tree_bundle(x$tree, x$species_origin, file.path(out_dir, "plot_dated_phylogeny_540Ma"),
                              "Case 01 dated phylogeny with unequal species origins"),
             save_plot_bundle(plot_beta(x$beta, "Case 01 HMSC Beta used for deep-time projection"),
                              file.path(out_dir, "plot_beta_heatmap")),
             save_plot_bundle(p_source, file.path(out_dir, "plot_source_land_area_full_axis")),
             save_plot_bundle(p_pred, file.path(out_dir, "plot_predictor_trajectories_full_axis"), width = 9, height = 7),
             save_plot_bundle(p_risk, file.path(out_dir, "plot_extrapolation_risk_full_axis")),
             save_plot_bundle(p_rich, file.path(out_dir, "plot_richness_full_axis")),
             save_plot_bundle(p_maps, file.path(out_dir, "plot_representative_richness_maps"), width = 10, height = 5.8),
             save_plot_bundle(p_turn, file.path(out_dir, "plot_turnover_full_axis")),
             save_plot_bundle(p_ref, file.path(out_dir, "plot_refugia_map")),
             save_plot_bundle(p_fossil, file.path(out_dir, "plot_fossil_validation_full_axis")),
             save_plot_bundle(plot_validation(x$modern_validation, "Modern validation of simulated HMSC projection"),
                              file.path(out_dir, "plot_modern_validation")),
             save_plot_bundle(plot_bgb(x$bgb, "BioGeoBEARS-style model weights and event probabilities"),
                              file.path(out_dir, "plot_bgb_accessibility_events")),
             save_plot_bundle(plot_workflow_steps(steps, "Case 01 full workflow completion"),
                              file.path(out_dir, "plot_workflow_steps")))

  data.frame(case = case, output_dir = normalizePath(out_dir, winslash = "/"),
             n_files = length(files), n_time_slices = length(x$times),
             n_projection_rows = nrow(x$projection), stringsAsFactors = FALSE)
}

run_case02 <- function() {
  case <- "case02_tracks540_20"
  x <- readRDS(file.path(script_dir, case, "case02_tracks540_20.rds"))
  out_dir <- file.path(out_root, "c02")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  files <- c(
    save_rds(x, file.path(out_dir, "case_input_and_all_results.rds")),
    save_rds(x$projection, file.path(out_dir, "projection_tracks_full_540_0Ma.rds")),
    save_rds(x$track_environment, file.path(out_dir, "track_environment_full_540_0Ma.rds")),
    save_rds(x$richness, file.path(out_dir, "richness_tracks_full_540_0Ma.rds")),
    save_rds(x$turnover, file.path(out_dir, "turnover_tracks_full_540_0Ma.rds")),
    save_rds(x$refugia, file.path(out_dir, "refugia_tracks_full_540_0Ma.rds"))
  )

  track_env_long <- reshape(
    x$track_environment[, c("track_id", "region", "time_ma", x$variables), drop = FALSE],
    varying = x$variables, v.names = "value", timevar = "variable",
    times = x$variables, direction = "long"
  )
  rownames(track_env_long) <- NULL
  rich_time <- stats::aggregate(expected_richness ~ time_ma, x$richness, mean, na.rm = TRUE)
  turnover_time <- summarise_turnover(x$turnover)
  turnover_time_plot <- turnover_time[is.finite(turnover_time$turnover), , drop = FALSE]
  component_time <- stats::aggregate(cbind(accessibility, dynamic_weight, extrapolation_weight,
                                           phylo_existence) ~ time_ma,
                                     x$projection, mean, na.rm = TRUE)
  acc_region <- stats::aggregate(accessibility ~ region + time_ma, x$projection, mean, na.rm = TRUE)
  refugia_region <- merge(x$refugia, unique(x$tracks[, c("track_id", "region")]),
                          by = "track_id", all.x = TRUE, sort = FALSE)
  steps <- workflow_steps(case, x)
  source_summary <- data.frame(
    source = x$source_level,
    n_time_slices = length(x$times),
    oldest_ma = max(x$times),
    youngest_ma = min(x$times),
    n_tracks = length(unique(x$tracks$track_id)),
    n_species = length(x$species),
    n_projection_rows = nrow(x$projection),
    n_richness_rows = nrow(x$richness),
    n_turnover_rows = nrow(x$turnover)
  )

  files <- c(files,
             write_csv(source_summary, file.path(out_dir, "source_and_size_summary.csv")),
             write_csv(data.frame(species = x$species, clade = x$clade[x$species],
                                  origin_ma = x$species_origin[x$species]),
                       file.path(out_dir, "species_origin.csv")),
             write_csv(x$beta, file.path(out_dir, "beta.csv")),
             write_csv(x$traits, file.path(out_dir, "traits.csv")),
             write_csv(x$tracks, file.path(out_dir, "plate_corrected_tracks_full_axis.csv")),
             write_csv(x$track_environment, file.path(out_dir, "track_environment_full_axis.csv")),
             write_csv(track_env_long, file.path(out_dir, "track_environment_long.csv")),
             write_csv(x$projection, file.path(out_dir, "track_projection_full_axis.csv")),
             write_csv(x$richness, file.path(out_dir, "track_richness_full_axis.csv")),
             write_csv(x$turnover, file.path(out_dir, "track_turnover_full_axis.csv")),
             write_csv(refugia_region, file.path(out_dir, "track_refugia_summary.csv")),
             write_csv(component_time, file.path(out_dir, "summary_component_weights_by_time.csv")),
             write_csv(acc_region, file.path(out_dir, "summary_accessibility_by_region_time.csv")),
             write_csv(x$modern_validation$evaluation, file.path(out_dir, "modern_validation_metrics.csv")),
             write_csv(x$bgb$compare, file.path(out_dir, "bgb_model_compare.csv")),
             write_csv(x$bgb$events, file.path(out_dir, "bgb_events.csv")),
             write_csv(steps, file.path(out_dir, "workflow_steps.csv")))

  p_tracks <- ggplot2::ggplot(x$tracks, ggplot2::aes(x = paleo_lon, y = paleo_lat,
                                                     color = region, group = region)) +
    ggplot2::geom_path(linewidth = 0.75, alpha = 0.9) +
    ggplot2::geom_point(ggplot2::aes(size = 540 - time_ma), alpha = 0.45) +
    ggplot2::coord_equal(xlim = c(-180, 180), ylim = c(-90, 90)) +
    ggplot2::labs(x = "Paleo-longitude", y = "Paleo-latitude",
                  size = "Elapsed Myr", color = "Region",
                  title = "Plate-corrected tracks spanning 540-0 Ma") +
    theme_case()

  p_env <- ggplot2::ggplot(track_env_long, ggplot2::aes(x = time_ma, y = value,
                                                       color = region, group = region)) +
    ggplot2::geom_line(linewidth = 0.55) +
    ggplot2::scale_x_reverse() +
    ggplot2::facet_wrap(~ variable, scales = "free_y", ncol = 2) +
    ggplot2::labs(x = "Time (Ma)", y = "Nearest-land paleoenvironment",
                  color = "Region", title = "Track environmental histories from the 4-degree RDS") +
    theme_case()

  p_rich <- ggplot2::ggplot(x$richness, ggplot2::aes(x = time_ma, y = expected_richness,
                                                     color = region, group = track_id)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::scale_x_reverse() +
    ggplot2::labs(x = "Time (Ma)", y = "Expected richness",
                  color = "Region", title = "Track richness through every 5 Myr slice") +
    theme_case()

  turnover_plot <- x$turnover[is.finite(x$turnover$turnover), , drop = FALSE]
  p_turn <- ggplot2::ggplot(turnover_plot, ggplot2::aes(x = time_to_ma, y = turnover,
                                                       color = cell_id, group = cell_id)) +
    ggplot2::geom_line(linewidth = 0.65) +
    ggplot2::scale_x_reverse() +
    ggplot2::labs(x = "Interval end (Ma)", y = "Turnover",
                  color = "Track", title = "Track turnover between adjacent 5 Myr slices") +
    theme_case()

  comp_long <- reshape(component_time, varying = c("accessibility", "dynamic_weight",
                                                   "extrapolation_weight", "phylo_existence"),
                       v.names = "value", timevar = "component",
                       times = c("accessibility", "dynamic", "extrapolation", "phylogeny"),
                       direction = "long")
  rownames(comp_long) <- NULL
  p_comp <- ggplot2::ggplot(comp_long, ggplot2::aes(x = time_ma, y = value,
                                                   color = component)) +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::scale_x_reverse() +
    ggplot2::labs(x = "Time (Ma)", y = "Mean component weight",
                  color = "Component", title = "Accessibility, dynamic filtering, extrapolation and phylogeny") +
    theme_case()

  p_acc <- ggplot2::ggplot(acc_region, ggplot2::aes(x = time_ma, y = region, fill = accessibility)) +
    ggplot2::geom_tile() +
    ggplot2::scale_x_reverse() +
    ggplot2::scale_fill_gradientn(colors = c(pal["navy"], "#F8F7F0", pal["red"])) +
    ggplot2::labs(x = "Time (Ma)", y = NULL, fill = "Mean\naccessibility",
                  title = "Model-averaged accessibility by region and time") +
    theme_case()

  p_ref <- ggplot2::ggplot(refugia_region, ggplot2::aes(x = reorder(track_id, refugia_score),
                                                        y = refugia_score, fill = region)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Track", y = "Refugia score", fill = "Region",
                  title = "Track refugia from long-term richness and stability") +
    theme_case()

  files <- c(files,
             save_tree_bundle(x$tree, x$species_origin, file.path(out_dir, "plot_dated_phylogeny_540Ma"),
                              "Case 02 dated phylogeny with unequal species origins"),
             save_plot_bundle(plot_beta(x$beta, "Case 02 HMSC Beta used for track projection"),
                              file.path(out_dir, "plot_beta_heatmap")),
             save_plot_bundle(p_tracks, file.path(out_dir, "plot_plate_corrected_tracks")),
             save_plot_bundle(p_env, file.path(out_dir, "plot_track_environment_full_axis"), width = 9, height = 7),
             save_plot_bundle(p_rich, file.path(out_dir, "plot_track_richness_full_axis")),
             save_plot_bundle(p_turn, file.path(out_dir, "plot_track_turnover_full_axis")),
             save_plot_bundle(p_comp, file.path(out_dir, "plot_projection_component_weights")),
             save_plot_bundle(p_acc, file.path(out_dir, "plot_accessibility_heatmap")),
             save_plot_bundle(p_ref, file.path(out_dir, "plot_track_refugia")),
             save_plot_bundle(plot_validation(x$modern_validation, "Modern validation of track-case projection"),
                              file.path(out_dir, "plot_modern_validation")),
             save_plot_bundle(plot_bgb(x$bgb, "BioGeoBEARS-style model weights and event probabilities"),
                              file.path(out_dir, "plot_bgb_accessibility_events")),
             save_plot_bundle(plot_workflow_steps(steps, "Case 02 full workflow completion"),
                              file.path(out_dir, "plot_workflow_steps")))

  data.frame(case = case, output_dir = normalizePath(out_dir, winslash = "/"),
             n_files = length(files), n_time_slices = length(x$times),
             n_projection_rows = nrow(x$projection), stringsAsFactors = FALSE)
}

notes_dir <- file.path(pkg_root, "inst", "extdata", "deep_time_workflow_notes")
coverage <- data.frame(word_functions = NA_integer_, exported_hee = NA_integer_, missing = NA_integer_)
if (dir.exists(notes_dir)) {
  txt_files <- list.files(notes_dir, pattern = "\\.txt$", full.names = TRUE)
  txt <- paste(unlist(lapply(txt_files, readLines, warn = FALSE)), collapse = "\n")
  hits <- gregexpr("\\bhee_[A-Za-z][A-Za-z0-9_]*\\b", txt, perl = TRUE)
  funcs <- unique(regmatches(txt, hits)[[1]])
  funcs <- funcs[vapply(funcs, function(z) length(gregexpr("hee_", z, fixed = TRUE)[[1]]) == 1, logical(1))]
  exports <- getNamespaceExports("HmscEcoEvo")
  coverage <- data.frame(word_functions = length(funcs),
                         exported_hee = sum(grepl("^hee_", exports)),
                         missing = length(setdiff(funcs, exports)))
}

summary <- rbind(run_case01(), run_case02())
write_csv(summary, file.path(out_root, "deep_time_case_run_summary.csv"))
write_csv(coverage, file.path(out_root, "word_function_coverage_runtime_check.csv"))
cat("Ran ", nrow(summary), " full 540-0 Ma deep-time cases and saved outputs under:\n",
    normalizePath(out_root, winslash = "/"), "\n", sep = "")
