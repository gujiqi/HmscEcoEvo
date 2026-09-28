#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
} else {
  normalizePath("HmscEcoEvo/inst/extdata/ecoevo_cases/run_all_case_plots.R", mustWork = TRUE)
}
case_root <- dirname(script_path)
pkg_dir <- normalizePath(file.path(case_root, "..", "..", ".."), mustWork = TRUE)

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(pkg_dir, quiet = TRUE)
} else {
  library(HmscEcoEvo)
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required to save diagnostic plots.", call. = FALSE)
}

output_root <- file.path(case_root, "outputs")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

case_dirs <- list.dirs(case_root, recursive = FALSE, full.names = TRUE)
case_dirs <- case_dirs[grepl("^case[0-9]{2}_", basename(case_dirs))]

assert_nonempty <- function(x, label) {
  ok <- FALSE
  if (is.data.frame(x) || is.matrix(x)) ok <- nrow(x) > 0 && ncol(x) > 0
  if (is.array(x)) ok <- length(x) > 0
  if (is.list(x) && !is.data.frame(x)) ok <- length(x) > 0
  if (!ok) stop(label, " is empty.", call. = FALSE)
  invisible(TRUE)
}

check_metrics <- function(metrics, case_name) {
  checks <- list(
    niche = c("beta_table", "beta_summary", "response_curves", "magnitude", "pca_scores"),
    phylo = c("global_signal", "correlogram", "scale_dependent_signal", "local_signal",
              "node_signal", "clade_conservatism", "residual_phylogenetic_map"),
    evo = c("branch_shifts", "shift_counts", "shift_support", "multivariate_shift",
            "rate_metrics", "variance_metrics", "dtt", "mdi", "early_burst",
            "convergence", "peak_reuse", "phylo_beta_distance",
            "evolutionary_integration", "trait_niche_covariance"),
    trait = c("gamma", "trait_explained_r2", "history_explained_r2",
              "combined_explained_r2", "residual_rho", "signal_decomposition",
              "missing_trait_risk", "trait_phylogenetic_redundancy",
              "trait_omission_sensitivity", "residual_beta_tree_heatmap"),
    gamma = c("gamma", "clade_specific_gamma", "regime_specific_gamma",
              "gamma_shift_probability", "gamma_sign_flip_index",
              "gamma_variance_among_clades", "trait_function_turnover_index",
              "pairwise_turnover_matrix", "trait_axis_specialization",
              "trait_function_network"),
    population = c("population_beta", "within_species_niche_divergence",
                   "local_adaptation_contribution", "plasticity_contribution",
                   "LA_PL_ratio", "variance_partition", "genetic_niche_signal",
                   "reaction_norms", "GxE", "eco_evolutionary_feedback_path",
                   "beta_t_shift_rate", "trait_evolution_contribution",
                   "community_driven_selection_index", "feedback_strength",
                   "eco_evo_turnover_ratio", "lagged_niche_shift",
                   "interaction_mediated_evolution"),
    validation = c("posterior_uncertainty_propagation", "tree_uncertainty_sensitivity",
                   "trait_omission_sensitivity", "spatial_confounding_check",
                   "environment_omission_check", "prior_sensitivity",
                   "simulation_recovery", "block_cross_validation",
                   "posterior_predictive_check", "mcmc_diagnostics", "trace")
  )
  for (module in names(checks)) {
    for (nm in checks[[module]]) {
      assert_nonempty(metrics[[module]][[nm]], paste(case_name, module, nm))
    }
  }
  invisible(TRUE)
}

print_plot_object <- function(plot_obj) {
  if (inherits(plot_obj, "hmsc_ecoevo_nature_page")) {
    print.hmsc_ecoevo_nature_page <- get("print.hmsc_ecoevo_nature_page",
                                         envir = asNamespace("HmscEcoEvo"))
    print.hmsc_ecoevo_nature_page(plot_obj)
  } else {
    print(plot_obj)
  }
}

save_plot_object <- function(plot_obj, base_file, width = 12, height = 8.4, dpi = 220) {
  pdf_file <- paste0(base_file, ".pdf")
  grDevices::pdf(pdf_file, width = width, height = height, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  print_plot_object(plot_obj)
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)

  files <- pdf_file
  if (inherits(plot_obj, "hmsc_ecoevo_nature_page")) {
    png_file <- paste0(base_file, ".png")
    grDevices::png(png_file, width = width, height = height, units = "in",
                   res = dpi, bg = "white")
    on.exit(grDevices::dev.off(), add = TRUE)
    print_plot_object(plot_obj)
    grDevices::dev.off()
    on.exit(NULL, add = FALSE)
    files <- c(files, png_file)
  } else if (inherits(plot_obj, "hmsc_ecoevo_plot_list")) {
    for (i in seq_along(plot_obj$plots)) {
      panel_file <- sprintf("%s_panel%02d.png", base_file, i)
      ggplot2::ggsave(panel_file, plot = plot_obj$plots[[i]], width = width,
                      height = height, dpi = dpi, bg = "white")
      files <- c(files, panel_file)
    }
  } else {
    png_file <- paste0(base_file, ".png")
    ggplot2::ggsave(png_file, plot = plot_obj, width = width,
                    height = height, dpi = dpi, bg = "white")
    files <- c(files, png_file)
  }
  files
}

run_case <- function(case_dir) {
  case_name <- basename(case_dir)
  evo <- readRDS(file.path(case_dir, "evo.rds"))
  axes <- evo$metadata$analysis_axes
  linear_axes <- if (is.null(evo$metadata$linear_axes)) {
    axes[!grepl("_sq$|_squared$|\\^2$", axes)]
  } else {
    evo$metadata$linear_axes
  }
  quadratic <- evo$metadata$quadratic_map
  if (is.null(quadratic)) {
    q_axes <- paste0(linear_axes, "_sq")
    quadratic <- data.frame(linear = linear_axes[q_axes %in% axes],
                            quadratic = q_axes[q_axes %in% axes],
                            stringsAsFactors = FALSE)
  }
  case_out <- file.path(output_root, case_name)
  if (dir.exists(case_out)) unlink(case_out, recursive = TRUE)
  dir.create(case_out, recursive = TRUE, showWarnings = FALSE)

  metrics <- list(
    niche = calc_niche_metrics(evo, axes = axes, quadratic = quadratic),
    phylo = calc_phylo_signal_metrics(evo, axes = linear_axes),
    evo = calc_evo_transition_metrics(evo, axes = linear_axes, p_shift_threshold = 0.5),
    trait = calc_trait_mediation_metrics(evo, axes = linear_axes),
    gamma = calc_gamma_evolution_metrics(evo, axes = linear_axes),
    population = calc_population_evolution_metrics(evo),
    validation = calc_validation_metrics(evo)
  )
  check_metrics(metrics, case_name)
  saveRDS(metrics, file.path(case_out, "metrics.rds"))

  plots <- list(
    niche_summary = plot_niche_summary(metrics$niche),
    phylo_signal = plot_phylo_signal(metrics$phylo),
    evo_transition = plot_evo_transition(metrics$evo),
    trait_mediation = plot_trait_mediation(metrics$trait),
    gamma_evolution = plot_gamma_evolution(metrics$gamma),
    population_evolution = plot_population_evolution(metrics$population),
    validation_dashboard = plot_validation_dashboard(metrics$validation)
  )

  plot_files <- unlist(Map(function(nm, obj) {
    save_plot_object(obj, file.path(case_out, nm))
  }, names(plots), plots), use.names = FALSE)

  data.frame(
    case = case_name,
    metrics_file = file.path(case_out, "metrics.rds"),
    plot_file_count = length(plot_files),
    output_dir = case_out,
    stringsAsFactors = FALSE
  )
}

summary <- do.call(rbind, lapply(case_dirs, run_case))
utils::write.csv(summary, file.path(output_root, "case_run_summary.csv"), row.names = FALSE)
cat("Ran ", nrow(summary), " cases and saved plots under:\n", output_root, "\n", sep = "")
