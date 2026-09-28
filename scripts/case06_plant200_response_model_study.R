#!/usr/bin/env Rscript

# Case06: matched BM/OU/EB ancestral-response comparison on real Plant200 data.
# All three models use the same HMSC tip Beta draws, dated tree and palaeo-Earth.
# Conditional branch means isolate model assumptions from bridge Monte Carlo.

root <- "C:/Users/Google/Documents/HMSC-HIST"
base <- file.path(root, "outputs", "HmscEcoEvo")
output_arg <- commandArgs(trailingOnly = TRUE)
output_arg <- output_arg[startsWith(output_arg, "--output=")]
if (length(output_arg) > 1L) stop("Supply --output only once.", call. = FALSE)
out <- if (length(output_arg)) {
  sub("^--output=", "", output_arg[[1L]])
} else {
  file.path(base, "case06_plant200_ancestral_response_BM_OU_EB")
}
fig_dir <- file.path(out, "figures")
data_dir <- file.path(out, "data")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
for (pkg in c("HmscEcoEvo", "ape", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Install required package: ", pkg, call. = FALSE)
  }
}
if (utils::packageVersion("HmscEcoEvo") < package_version("1.0.4") ||
    sum(HmscEcoEvo::hee_core_process_catalog()$layer == "core_process") != 6L) {
  stop("Case06 requires HmscEcoEvo >= 1.0.4 with the six-process API.",
       call. = FALSE)
}

source_case <- file.path(base, "case05_v8_full200_terrain_resistance_20260926")
input_tree <- file.path(root, "HmscEcoEvo_API8_20260818", "data_external",
  "plant_200_multifamily_extant_20260820", "prepared_inputs",
  "case03_data_raw", "Plant200_multifamily_extant_1deg_traits", "tree.tre")
beta_file <- file.path(base,
  "case05_plant200_4deg_effortaware_formal25_all66_20260923",
  "02_hmsc_4deg", "hmsc_environmental_beta_posterior_draws.rds")
fit_file <- file.path(base,
  "case05_plant200_4deg_effortaware_formal25_all66_20260923",
  "02_hmsc_4deg", "case05_v5_hmsc_4deg_target_group_model.rds")
cache_root <- file.path(root, "HmscEcoEvo", "derived_inputs",
  "case05_v2_paleomap_h3_carrier_325_0Ma_1deg",
  "movement_cache_arias_topographic_rate0p4")
cache_index <- utils::read.csv(file.path(cache_root,
  "case05_v2_movement_cache_index.csv"))
times <- sort(unique(cache_index$time_ma), decreasing = TRUE)
if (length(times) != 66L || !isTRUE(all.equal(times, seq(325, 0, by = -5)))) {
  stop("Case06 requires the 66 native 325-0 Ma palaeo-Earth slices.")
}

selection <- utils::read.csv(file.path(source_case, "01_input_audit",
  "species_selection.csv"))
tips <- selection$species[selection$selected]
tree <- ape::reorder.phylo(ape::keep.tip(ape::read.tree(input_tree), tips),
                           "cladewise")
n_tip <- length(tree$tip.label)
if (n_tip != 200L) stop("Case06 requires the saved Plant200 tip selection.")
depth <- ape::node.depth.edgelength(tree)
root_age <- max(depth[seq_len(n_tip)])
tip_lag <- root_age - depth[seq_len(n_tip)]
if (max(abs(tip_lag)) > 1e-4) stop("Tree is not near-ultrametric.")
terminal_edge <- match(seq_len(n_tip), tree$edge[, 2L])
tree$edge.length[terminal_edge] <- tree$edge.length[terminal_edge] + tip_lag
root_age <- max(ape::node.depth.edgelength(tree)[seq_len(n_tip)])
if (root_age > 325 + 1e-4) {
  tree$edge.length <- tree$edge.length * (325 / root_age)
}

axes <- c("MAT_pohl_C", "MAP_pohl_mm_yr", "elevation_m",
  "T_seasonality_pohl_sd", "P_seasonality_pohl_sd",
  "moisture_availability_index_z", "wetland_potential_index_z",
  "relief_3x3_m")
draw_ids <- unique(utils::read.csv(file.path(source_case,
  "04_phylogeographic_pruning", "pruning_log_likelihoods.csv"))$response_draw)
tip_beta <- readRDS(beta_file)
tip_beta <- tip_beta[tip_beta$lineage %in% tree$tip.label &
  tip_beta$response_draw %in% draw_ids,
  c("lineage", "response_draw", axes), drop = FALSE]
names(tip_beta)[names(tip_beta) == "lineage"] <- "species"
if (nrow(tip_beta) != n_tip * length(draw_ids) || length(draw_ids) != 5L) {
  stop("The saved five-draw HMSC tip Beta table is incomplete.")
}
modern_fit <- readRDS(fit_file)
training_x <- as.matrix(modern_fit$XData[, axes, drop = FALSE])
rm(modern_fit)
training_bounds <- apply(training_x, 2L, stats::quantile,
  probs = c(.01, .99), na.rm = TRUE)

models <- c("BM", "OU", "EB")
palette <- c(BM = "#2764A4", OU = "#16856A", EB = "#D68416")
responses <- setNames(vector("list", length(models)), models)
fit_tables <- list()
for (model in models) {
  file <- file.path(data_dir, paste0("ancestral_beta_conditional_", model, ".rds"))
  if (file.exists(file) && Sys.getenv("HMSCEE_CASE06_FORCE") != "1") {
    responses[[model]] <- readRDS(file)
  } else {
    message("Case06 ancestral response: ", model)
    responses[[model]] <- HmscEcoEvo::hee_evolution_ancestral_response_direct(
      tip_beta_draws = tip_beta, tree = tree, times = times,
      basis_cols = axes, evolution_model = model,
      branch_uncertainty = "none", reconstruct_intercept = FALSE)
    saveRDS(responses[[model]], file, compress = "gzip")
  }
  fit <- attr(responses[[model]], "model_comparison")
  if (is.null(fit) || nrow(fit) != length(draw_ids) * length(axes)) {
    stop("Missing model-fit diagnostics for ", model)
  }
  fit_tables[[model]] <- fit
}
fit_table <- do.call(rbind, fit_tables)
fit_table$group <- paste(fit_table$response_draw, fit_table$axis)
fit_table$delta_aicc <- fit_table$aicc - ave(fit_table$aicc, fit_table$group, FUN = min)
fit_table$selected_by_aicc <- fit_table$delta_aicc < 1e-8
utils::write.csv(fit_table, file.path(data_dir, "model_comparison_by_axis_draw.csv"),
                 row.names = FALSE)

parent <- setNames(tree$edge[, 1L], as.character(tree$edge[, 2L]))
focal <- c("Ilex_aquifolium", "Camellia_sinensis", "Jasminum_officinale")
if (!all(focal %in% tree$tip.label)) stop("Focal taxa are absent from the tree.")
tip_path <- function(tip) {
  node <- match(tip, tree$tip.label)
  path <- integer()
  while (!is.na(node)) {
    path <- c(path, node)
    next_node <- parent[as.character(node)]
    node <- if (is.na(next_node)) NA_integer_ else as.integer(next_node)
  }
  path
}
paths <- setNames(lapply(focal, tip_path), focal)
trajectory <- list()
for (model in models) {
  x <- responses[[model]]
  for (tip in focal) {
    z <- x[x$node %in% paths[[tip]],
           c("time_ma", "response_draw", "node", axes[1:2]), drop = FALSE]
    if (nrow(z) != length(times) * length(draw_ids)) {
      stop("Focal ancestral path is incomplete: ", tip, ", ", model)
    }
    for (axis in axes[1:2]) {
      trajectory[[length(trajectory) + 1L]] <- data.frame(
        model = model, tip = tip, axis = axis,
        time_ma = z$time_ma, response_draw = z$response_draw,
        ancestor_node = z$node, beta = z[[axis]])
    }
  }
}
trajectory <- do.call(rbind, trajectory)
utils::write.csv(trajectory, file.path(data_dir,
  "three_focal_tip_paths_beta_325_0Ma.csv"), row.names = FALSE)

grid <- data.frame(cell_id = sprintf("g4_%04d", seq_len(4050L)),
  lon = rep(seq(-178, 178, by = 4), times = 45),
  lat = rep(seq(-88, 88, by = 4), each = 90))
bin4 <- function(lon, lat) {
  ix <- pmin(89L, pmax(0L, floor((lon + 180) / 4)))
  iy <- pmin(44L, pmax(0L, floor((lat + 90) / 4)))
  as.integer(iy * 90L + ix + 1L)
}
slice_from_carrier <- function(carrier) {
  land <- !is.na(carrier$active_land) & carrier$active_land &
    is.finite(carrier$land_area_km2) & carrier$land_area_km2 > 0 &
    is.finite(carrier$paleo_lon) & is.finite(carrier$paleo_lat)
  x <- carrier[land, , drop = FALSE]
  x$bin <- bin4(x$paleo_lon, x$paleo_lat)
  area <- numeric(nrow(grid))
  total <- rowsum(x$land_area_km2, x$bin, reorder = FALSE)
  area[as.integer(rownames(total))] <- total[, 1L]
  env <- matrix(NA_real_, nrow(grid), length(axes),
                dimnames = list(NULL, axes))
  complete <- stats::complete.cases(x[, axes, drop = FALSE])
  x <- x[complete, , drop = FALSE]
  if (nrow(x)) {
    den <- rowsum(x$land_area_km2, x$bin, reorder = FALSE)
    ids <- as.integer(rownames(den))
    for (axis in axes) {
      num <- rowsum(x$land_area_km2 * x[[axis]], x$bin,
                     reorder = FALSE)
      env[ids, axis] <- num[, 1L] / den[, 1L]
    }
  }
  list(area = area, env = env, valid = area > 0 &
         stats::complete.cases(env))
}

selected_times <- c(300, 200, 100, 0)
map_rows <- list()
global_rows <- list()
noanalog_rows <- list()
for (tm in times) {
  file <- cache_index$file[match(tm, cache_index$time_ma)]
  slice <- slice_from_carrier(readRDS(file)$carrier)
  valid <- slice$valid
  if (!any(valid)) stop("No complete land environment at ", tm, " Ma")
  X <- slice$env[valid, axes, drop = FALSE]
  a <- slice$area[valid]
  noanalog <- rowSums(sweep(X, 2L, training_bounds[1L, ], "<") |
    sweep(X, 2L, training_bounds[2L, ], ">")) > 0
  noanalog_rows[[length(noanalog_rows) + 1L]] <- data.frame(
    time_ma = tm, n_supported_land_cells = length(a),
    noanalog_area_fraction = sum(a[noanalog]) / sum(a))
  if (tm %% 25 == 0) message("Case06 palaeoenvironment: ", tm, " Ma")
  for (model in models) {
    response <- responses[[model]]
    draw_mean <- matrix(NA_real_, length(a), length(draw_ids))
    focal_mean <- matrix(NA_real_, length(a), length(draw_ids))
    n_active <- integer(length(draw_ids))
    for (di in seq_along(draw_ids)) {
      rr <- response[response$time_ma == tm &
        response$response_draw == draw_ids[di], , drop = FALSE]
      if (!nrow(rr)) stop("Missing active responses at ", tm, " Ma")
      beta <- as.matrix(rr[, axes, drop = FALSE])
      draw_mean[, di] <- rowMeans(stats::pnorm(X %*% t(beta)))
      n_active[di] <- nrow(rr)
      focal_row <- rr[rr$node %in% paths[[focal[1L]]], , drop = FALSE]
      if (nrow(focal_row) != 1L) {
        stop("Focal lineage is not unique at ", tm, " Ma")
      }
      focal_mean[, di] <- as.numeric(stats::pnorm(
        X %*% as.numeric(focal_row[1L, axes])))
    }
    if (length(unique(n_active)) != 1L) {
      stop("Active-lineage count changes between HMSC draws.")
    }
    global_rows[[length(global_rows) + 1L]] <- data.frame(
      time_ma = tm, model = model, response_draw = draw_ids,
      n_active_sampled_lineages = n_active,
      area_weighted_mean_support = colSums(draw_mean * a) / sum(a),
      focal_tip = focal[1L],
      focal_area_weighted_support = colSums(focal_mean * a) / sum(a))
    if (tm %in% selected_times) {
      for (metric in c("mean_active_lineage_support", "focal_path_support")) {
        state <- if (metric == "mean_active_lineage_support") draw_mean else
          focal_mean
        value <- rep(NA_real_, nrow(grid))
        value[valid] <- rowMeans(state)
        map_rows[[length(map_rows) + 1L]] <- data.frame(
          grid, time_ma = tm, model = model, metric = metric,
          support = value, land_exists = slice$area > 0,
          land_area_km2 = slice$area,
          noanalog_flag = rep(NA, nrow(grid)))
        map_rows[[length(map_rows)]]$noanalog_flag[valid] <- noanalog
      }
    }
  }
}
global_summary <- do.call(rbind, global_rows)
noanalog_summary <- do.call(rbind, noanalog_rows)
maps <- do.call(rbind, map_rows)
utils::write.csv(global_summary, file.path(data_dir,
  "global_environmental_support_by_time_model_draw.csv"), row.names = FALSE)
utils::write.csv(noanalog_summary, file.path(data_dir,
  "noanalog_diagnostic_by_time.csv"), row.names = FALSE)
saveRDS(maps, file.path(data_dir, "selected_time_support_maps_4deg.rds"),
        compress = "gzip")

bm <- maps[maps$model == "BM", c("metric", "time_ma", "cell_id", "support")]
other <- maps[maps$model != "BM", , drop = FALSE]
other$baseline <- bm$support[match(paste(other$metric, other$time_ma,
                                         other$cell_id),
                                   paste(bm$metric, bm$time_ma, bm$cell_id))]
other$difference <- other$support - other$baseline
groups <- split(other, interaction(other$metric, other$time_ma,
                                   other$model, drop = TRUE))
difference_summary <- do.call(rbind, lapply(groups, function(z) {
  valid <- is.finite(z$difference) & is.finite(z$land_area_km2) &
    z$land_area_km2 > 0
  d <- abs(z$difference[valid])
  a <- z$land_area_km2[valid]
  data.frame(metric = as.character(z$metric[1L]), time_ma = z$time_ma[1L],
    comparison = paste0(as.character(z$model[1L]), " - BM"),
    n_land_cells = sum(valid),
    area_weighted_mean_abs_difference = sum(d * a) / sum(a),
    max_abs_difference = max(d),
    land_area_fraction_difference_gt_0p05 = sum(a[d > .05]) / sum(a))
}))
rownames(difference_summary) <- NULL
utils::write.csv(difference_summary, file.path(data_dir,
  "model_difference_summary.csv"), row.names = FALSE)
diff_limit <- max(abs(other$difference), na.rm = TRUE)
if (!is.finite(diff_limit) || diff_limit <= 0) diff_limit <- 1e-6
saveRDS(other, file.path(data_dir, "selected_time_model_difference_maps_4deg.rds"),
        compress = "gzip")

theme_case <- ggplot2::theme_minimal(base_size = 16) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", size = 15),
    plot.title = ggplot2::element_text(face = "bold", size = 20),
    plot.subtitle = ggplot2::element_text(size = 12),
    axis.text = ggplot2::element_text(size = 11),
    legend.title = ggplot2::element_text(size = 12))
save_plot <- function(plot, name, width, height) {
  ggplot2::ggsave(file.path(fig_dir, name), plot = plot,
    width = width, height = height, units = "in", dpi = 220, bg = "white")
}
axis_labels <- c(MAT_pohl_C = "Mean annual temperature",
  MAP_pohl_mm_yr = "Annual precipitation")

fit_plot <- ggplot2::ggplot(fit_table,
  ggplot2::aes(x = delta_aicc, y = axis, colour = model,
               shape = boundary_status)) +
  ggplot2::geom_point(position = ggplot2::position_jitter(height = .12),
                      size = 3, alpha = .8) +
  ggplot2::scale_colour_manual(values = palette) +
  ggplot2::scale_x_continuous(trans = "log1p") +
  ggplot2::labs(title = "Model comparison across eight HMSC response axes",
    subtitle = "Five posterior draws per axis; points on parameter bounds are flagged",
    x = expression(Delta * "AICc (log1p display scale)"),
    y = NULL, colour = "Evolution model", shape = "Fit status") +
  theme_case
save_plot(fit_plot, "01_model_fit_AICc_and_boundaries.png", 14, 8)

trajectory$axis_label <- unname(axis_labels[trajectory$axis])
trajectory$tip_label <- gsub("_", " ", trajectory$tip)
trajectory$model <- factor(trajectory$model, levels = models)
trajectory_mean <- stats::aggregate(beta ~ model + tip_label + axis_label +
  time_ma, trajectory, mean)
trajectory_plot <- ggplot2::ggplot() +
  ggplot2::geom_line(data = trajectory,
    ggplot2::aes(x = time_ma, y = beta, colour = model,
                 group = interaction(model, response_draw)), alpha = .16,
    linewidth = .35) +
  ggplot2::geom_line(data = trajectory_mean,
    ggplot2::aes(x = time_ma, y = beta, colour = model,
                 linetype = model), linewidth = 1.15) +
  ggplot2::facet_grid(axis_label ~ tip_label, scales = "free_y") +
  ggplot2::scale_colour_manual(values = palette) +
  ggplot2::scale_linetype_manual(values = c(BM = "solid", OU = "solid",
                                           EB = "longdash")) +
  ggplot2::scale_x_reverse(breaks = c(300, 200, 100, 0)) +
  ggplot2::labs(title = "Three models, the same modern HMSC responses",
    subtitle = "Focal tip-to-root paths; faint lines are five posterior draws",
    x = "Time before present (Ma)", y = "Ancestral response coefficient",
    colour = "Model") + theme_case
save_plot(trajectory_plot, "02_focal_ancestral_beta_trajectories.png", 17, 10)

root_rows <- list()
for (model in models) {
  x <- responses[[model]]
  z <- x[x$lineage_type == "root", , drop = FALSE]
  for (axis in axes) {
    root_rows[[length(root_rows) + 1L]] <- data.frame(
      model = model, axis = axis, response_draw = z$response_draw,
      beta = z[[axis]])
  }
}
root_beta <- do.call(rbind, root_rows)
utils::write.csv(root_beta, file.path(data_dir, "root_beta_by_model_draw.csv"),
                 row.names = FALSE)
root_plot <- ggplot2::ggplot(root_beta,
  ggplot2::aes(x = beta, y = axis, colour = model)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = .4, colour = "grey65") +
  ggplot2::geom_point(position = ggplot2::position_jitterdodge(
    jitter.width = .08, jitter.height = 0, dodge.width = .65),
    size = 2.4, alpha = .8) +
  ggplot2::scale_colour_manual(values = palette) +
  ggplot2::labs(title = "Root ancestral response at 325 Ma",
    subtitle = "Each point is one HMSC posterior draw; axes are fitted-model inputs",
    x = "Response coefficient", y = NULL, colour = "Model") + theme_case
save_plot(root_plot, "03_root_beta_model_comparison.png", 14, 8)

global_mean <- stats::aggregate(cbind(area_weighted_mean_support,
  focal_area_weighted_support) ~ time_ma + model, global_summary, mean)
trend_plot <- ggplot2::ggplot(global_mean,
  ggplot2::aes(x = time_ma, y = area_weighted_mean_support,
               colour = model, linetype = model)) +
  ggplot2::geom_line(linewidth = 1.25) +
  ggplot2::scale_colour_manual(values = palette) +
  ggplot2::scale_linetype_manual(values = c(BM = "solid", OU = "solid",
                                           EB = "longdash")) +
  ggplot2::scale_x_reverse(breaks = c(325, 300, 250, 200, 150, 100, 50, 0)) +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  ggplot2::labs(title = "Global land-weighted environmental support",
    subtitle = "Mean across active sampled-surviving lineages; not occupancy or richness",
    x = "Time before present (Ma)", y = "Mean support index",
    colour = "Model") + theme_case
save_plot(trend_plot, "04_global_environmental_support_through_time.png", 14, 6.5)

map_time_label <- function(x) factor(paste0(x, " Ma"),
                                    levels = paste0(selected_times, " Ma"))
maps$time_label <- map_time_label(maps$time_ma)
maps$model <- factor(maps$model, levels = models)
other$time_label <- map_time_label(other$time_ma)
other$comparison <- factor(paste0(other$model, " - BM"),
                           levels = c("OU - BM", "EB - BM"))
for (metric in c("mean_active_lineage_support", "focal_path_support")) {
  z <- maps[maps$metric == metric, , drop = FALSE]
  title <- if (metric == "focal_path_support") {
    paste("Environmental support along the", gsub("_", " ", focal[1L]),
          "ancestral path")
  } else "Mean environmental support across active lineages"
  p <- ggplot2::ggplot(z, ggplot2::aes(lon, lat, fill = support)) +
    ggplot2::geom_raster() + ggplot2::facet_grid(time_label ~ model) +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), na.value = "white") +
    ggplot2::labs(title = title,
      subtitle = "Same 4-degree palaeo-Earth and fixed [0,1] scale; ocean is blank",
      x = "Palaeolongitude", y = "Palaeolatitude", fill = "Support") + theme_case
  save_plot(p, paste0(if (metric == "focal_path_support") "05" else "06",
                      "_", metric, "_four_times.png"), 16, 13)
  d <- other[other$metric == metric, , drop = FALSE]
  p_diff <- ggplot2::ggplot(d, ggplot2::aes(lon, lat, fill = difference)) +
    ggplot2::geom_raster() + ggplot2::facet_grid(time_label ~ comparison) +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::scale_fill_gradient2(low = "#2D5DA8", mid = "#E5E9E7",
      high = "#C84536", midpoint = 0,
      limits = c(-diff_limit, diff_limit), na.value = "white") +
    ggplot2::labs(title = paste("Model difference:", title),
      subtitle = "A single symmetric legend across all times and both metrics",
      x = "Palaeolongitude", y = "Palaeolatitude",
      fill = "Support difference") + theme_case
  save_plot(p_diff,
    paste0(if (metric == "focal_path_support") "07" else "08",
           "_", metric, "_difference_four_times.png"), 14, 13)
}

valid_land <- maps$land_exists & !is.na(maps$support)
wrong_ocean <- any(!maps$land_exists & is.finite(maps$support))
present <- maps[maps$time_ma == 0, , drop = FALSE]
present_bm <- present$support[present$model == "BM"]
present_ou <- present$support[present$model == "OU"]
present_eb <- present$support[present$model == "EB"]
present_equal <- isTRUE(all.equal(present_bm, present_ou,
  tolerance = 1e-10)) && isTRUE(all.equal(present_bm, present_eb,
  tolerance = 1e-10))
qc <- data.frame(
  check = c("three_real_fitted_models", "five_shared_hmsc_draws",
    "sixty_six_native_earth_slices", "all_conditional_beta_finite",
    "all_support_in_unit_interval", "ocean_has_no_support",
    "modern_tip_endpoint_equal_across_models", "aicc_boundaries_reported"),
  pass = c(length(models) == 3L, length(draw_ids) == 5L,
    length(times) == 66L,
    all(vapply(responses, function(x) all(is.finite(
      as.matrix(x[, axes, drop = FALSE]))), logical(1))),
    all(maps$support[valid_land] >= 0 & maps$support[valid_land] <= 1),
    !wrong_ocean, present_equal,
    all(c("boundary_status", "delta_aicc") %in% names(fit_table))),
  stringsAsFactors = FALSE)
utils::write.csv(qc, file.path(out, "case06_quality_gates.csv"),
                 row.names = FALSE)
if (!all(qc$pass)) stop("Case06 quality gate failure: ",
                       paste(qc$check[!qc$pass], collapse = ", "))
message("Case06 complete: ", out)
