#' Calculate validation and robustness diagnostics
#'
#' Computes and standardizes validation outputs: posterior uncertainty
#' propagation, tree uncertainty sensitivity, trait omission sensitivity,
#' spatial confounding checks, environment omission checks, prior sensitivity,
#' simulation recovery, block cross-validation, posterior predictive checks, and
#' MCMC diagnostics. Expensive checks should be supplied as precomputed tables.
#'
#' @param evo A `hmsc_ecoevo` object.
#' @param precomputed Optional validation tables.
#' @param ... Additional arguments passed to [as_hmsc_ecoevo()].
#' @return A `hmsc_validation_metrics` object.
#' @export
calc_validation_metrics <- function(evo,
                                    precomputed = NULL,
                                    ...) {
  evo <- as_hmsc_ecoevo(evo, ...)
  pre <- utils::modifyList(evo$precomputed, .as_precomputed_list(precomputed))
  diagnostics <- list(messages = character(), warnings = character())

  posterior_uncertainty <- .posterior_uncertainty_propagation(evo)
  tree_sensitivity <- .pre_or_na(pre$tree_uncertainty_sensitivity,
    data.frame(metric = "tree_uncertainty_sensitivity", value = NA_real_,
               source = "requires_precomputed_tree_set"))
  trait_omission <- if (!is.null(pre$trait_omission_sensitivity)) {
    as.data.frame(pre$trait_omission_sensitivity)
  } else calc_trait_mediation_metrics(evo)$trait_omission_sensitivity
  spatial <- .spatial_confounding_check(evo$XData, evo$coords)
  environment_omission <- .pre_or_na(pre$environment_omission_check,
    data.frame(metric = "environment_omission_check", value = NA_real_,
               source = "requires_precomputed_model_comparison"))
  prior_sensitivity <- .pre_or_na(pre$prior_sensitivity,
    data.frame(metric = "prior_sensitivity", value = NA_real_,
               source = "requires_precomputed_prior_sensitivity"))
  simulation_recovery <- .pre_or_na(pre$simulation_recovery,
    data.frame(true_value = NA_real_, estimated_value = NA_real_,
               metric = "simulation_recovery",
               source = "requires_precomputed_simulation"))
  block_cv <- .pre_or_na(pre$block_cross_validation,
    data.frame(block = NA_character_, metric = NA_character_, value = NA_real_,
               source = "requires_precomputed_block_cv"))
  ppc <- .pre_or_na(pre$posterior_predictive_check,
    data.frame(metric = "posterior_predictive_check", value = NA_real_,
               source = "requires_precomputed_ppc"))
  mcmc <- .mcmc_diagnostics(evo, pre)
  trace <- .trace_table(evo, pre)

  for (msg in c(
    "Tree uncertainty sensitivity requires precomputed results over a tree set.",
    "Environment omission, prior sensitivity, simulation recovery, block CV, and PPC are model-comparison or simulation workflows and should be supplied through precomputed tables when not already available."
  )) diagnostics$messages <- c(diagnostics$messages, msg)

  out <- list(
    posterior_uncertainty_propagation = posterior_uncertainty,
    tree_uncertainty_sensitivity = tree_sensitivity,
    trait_omission_sensitivity = trait_omission,
    spatial_confounding_check = spatial,
    environment_omission_check = environment_omission,
    prior_sensitivity = prior_sensitivity,
    simulation_recovery = simulation_recovery,
    block_cross_validation = block_cv,
    posterior_predictive_check = ppc,
    mcmc_diagnostics = mcmc,
    trace = trace,
    diagnostics = diagnostics,
    source = evo
  )
  class(out) <- "hmsc_validation_metrics"
  out
}

#' Plot validation and robustness diagnostics
#'
#' Returns the seventh HmscEcoEvo diagnostic figure: simulation recovery,
#' posterior uncertainty propagation, tree uncertainty sensitivity, model
#' sensitivity checks, block cross-validation, and posterior predictive/MCMC
#' diagnostics.
#'
#' @param x A `hmsc_ecoevo` or `hmsc_validation_metrics` object.
#' @param style Figure style. `"nature"` returns a single full-page
#'   Nature-style diagnostic plate; `"classic"` returns the original compact
#'   ggplot/patchwork summary.
#' @param ... Arguments passed to [calc_validation_metrics()] when needed.
#' @return A `ggplot`, `patchwork`, or printable plot-list object.
#' @export
plot_validation_dashboard <- function(x, style = c("nature", "classic"), ...) {
  style <- match.arg(style)
  if (identical(style, "nature")) return(.plot_validation_dashboard_nature(x, ...))
  .gg_ok()
  m <- if (inherits(x, "hmsc_validation_metrics")) x else calc_validation_metrics(x, ...)
  S <- m$simulation_recovery
  U <- m$posterior_uncertainty_propagation
  T <- m$tree_uncertainty_sensitivity
  P <- m$prior_sensitivity
  E <- m$environment_omission_check
  CV <- m$block_cross_validation
  PPC <- m$posterior_predictive_check
  M <- m$mcmc_diagnostics
  TO <- m$trait_omission_sensitivity
  SP <- m$spatial_confounding_check
  TR <- m$trace

  p1 <- if (all(c("true_value", "estimated_value") %in% names(S))) {
    ggplot2::ggplot(S, ggplot2::aes(true_value, estimated_value)) +
      ggplot2::geom_point(color = "#2b6cb0", na.rm = TRUE) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
      ggplot2::labs(title = "A. Simulation recovery",
                    x = "True value", y = "Estimated value") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p2 <- if (nrow(U) > 0) {
    ggplot2::ggplot(U, ggplot2::aes(metric, value, fill = axis)) +
      ggplot2::geom_col(position = "dodge", width = 0.7, na.rm = TRUE) +
      ggplot2::labs(title = "B. Posterior uncertainty propagation",
                    x = "Metric", y = "Uncertainty") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p3 <- if (nrow(T) > 0) {
    ggplot2::ggplot(T, ggplot2::aes(metric, value, fill = metric)) +
      ggplot2::geom_col(width = 0.7, na.rm = TRUE) +
      ggplot2::labs(title = "C. Tree uncertainty sensitivity",
                    x = "Metric", y = "Value") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     legend.position = "none")
  } else NULL

  sens <- rbind(.standard_metric_table(P, "prior_sensitivity"),
                .standard_metric_table(E, "environment_omission"))
  p4 <- if (nrow(sens) > 0) {
    ggplot2::ggplot(sens, ggplot2::aes(metric, value, fill = source_group)) +
      ggplot2::geom_col(position = "dodge", width = 0.7, na.rm = TRUE) +
      ggplot2::labs(title = "D. Model sensitivity checks",
                    x = "Metric", y = "Value") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p5 <- if (nrow(CV) > 0 && all(c("block", "value") %in% names(CV))) {
    ggplot2::ggplot(CV, ggplot2::aes(block, value, fill = metric)) +
      ggplot2::geom_col(position = "dodge", width = 0.7, na.rm = TRUE) +
      ggplot2::labs(title = "E. Block cross-validation",
                    x = "Block", y = "Value") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  ppc_mcmc <- rbind(.standard_metric_table(PPC, "posterior_predictive_check"),
                    .standard_metric_table(M, "mcmc"))
  p6 <- if (nrow(ppc_mcmc) > 0) {
    ggplot2::ggplot(ppc_mcmc, ggplot2::aes(metric, value, fill = source_group)) +
      ggplot2::geom_col(position = "dodge", width = 0.7, na.rm = TRUE) +
      ggplot2::facet_wrap(~source_group, scales = "free_y") +
      ggplot2::labs(title = "F. PPC and MCMC diagnostics",
                    x = "Metric", y = "Value") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  } else NULL

  p7 <- if (nrow(TO) > 0 && all(c("omitted_trait", "delta_r2") %in% names(TO))) {
    ggplot2::ggplot(TO, ggplot2::aes(omitted_trait, delta_r2, fill = axis)) +
      ggplot2::geom_col(position = "dodge", width = 0.7, na.rm = TRUE) +
      ggplot2::labs(title = "G. Trait omission sensitivity",
                    x = "Omitted trait", y = "Delta R2") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  } else NULL

  p8 <- if (nrow(SP) > 0 && all(c("variable", "r2_with_coords") %in% names(SP))) {
    ggplot2::ggplot(SP, ggplot2::aes(variable, r2_with_coords, fill = variable)) +
      ggplot2::geom_col(width = 0.7, na.rm = TRUE) +
      ggplot2::labs(title = "H. Spatial confounding check",
                    x = "Predictor", y = "R2 with coordinates") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     legend.position = "none")
  } else NULL

  p9 <- if (nrow(TR) > 0 && all(c("draw", "parameter", "value") %in% names(TR))) {
    ggplot2::ggplot(TR, ggplot2::aes(draw, value, color = parameter)) +
      ggplot2::geom_line(linewidth = 0.5, alpha = 0.8, na.rm = TRUE) +
      ggplot2::labs(title = "I. Trace diagnostics",
                    x = "Draw", y = "Value") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  .combine_plots(list(p1, p2, p3, p4, p5, p6, p7, p8, p9),
                 title = "Validation and Robustness", ncol = 2)
}

.posterior_uncertainty_propagation <- function(evo) {
  draws <- evo$beta_draws
  if (is.null(draws)) {
    return(data.frame(metric = "beta_posterior_sd", axis = evo$axes,
                      value = NA_real_, source = "requires_beta_draws",
                      stringsAsFactors = FALSE))
  }
  out <- list()
  k <- 1L
  for (axis in evo$axes) {
    z <- draws[, , axis, drop = FALSE]
    out[[k]] <- data.frame(metric = "mean_beta_posterior_sd", axis = axis,
                           value = mean(apply(z[, , 1, drop = FALSE], 2, stats::sd, na.rm = TRUE), na.rm = TRUE),
                           source = "beta_draws", stringsAsFactors = FALSE)
    k <- k + 1L
  }
  mag <- apply(draws, 1, function(m) {
    mm <- matrix(m, nrow = length(evo$species), ncol = length(evo$axes))
    sqrt(rowSums(mm^2, na.rm = TRUE))
  })
  out[[k]] <- data.frame(metric = "mean_magnitude_posterior_sd", axis = "all",
                         value = mean(apply(mag, 1, stats::sd, na.rm = TRUE), na.rm = TRUE),
                         source = "beta_draws", stringsAsFactors = FALSE)
  do.call(rbind, out)
}

.pre_or_na <- function(x, fallback) {
  if (is.null(x)) fallback else as.data.frame(x)
}

.spatial_confounding_check <- function(XData, coords) {
  if (is.null(XData) || is.null(coords) || ncol(coords) < 2) {
    return(data.frame(variable = NA_character_, r2_with_coords = NA_real_,
                      source = "requires_XData_and_coords", stringsAsFactors = FALSE))
  }
  X <- XData[vapply(XData, is.numeric, logical(1))]
  C <- coords[vapply(coords, is.numeric, logical(1))]
  if (ncol(X) == 0 || ncol(C) < 2) return(data.frame())
  out <- lapply(names(X), function(v) {
    data.frame(variable = v, r2_with_coords = .safe_lm_r2(X[[v]], C),
               source = "lm_variable_on_coordinates", stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

.mcmc_diagnostics <- function(evo, pre) {
  if (!is.null(pre$mcmc_diagnostics)) return(as.data.frame(pre$mcmc_diagnostics))
  draws <- evo$beta_draws
  if (is.null(draws)) {
    return(data.frame(metric = c("ESS", "Rhat"), value = NA_real_,
                      parameter = "beta", source = "requires_draws_or_precomputed",
                      stringsAsFactors = FALSE))
  }
  mat <- matrix(draws, nrow = dim(draws)[1])
  colnames(mat) <- paste0("beta_", seq_len(ncol(mat)))
  ess <- if (.soft_pkg("coda")) {
    coda::effectiveSize(coda::mcmc(mat))
  } else rep(NA_real_, ncol(mat))
  # `coda::effectiveSize()` may legitimately return numeric(0) for a
  # degenerate or too-short chain.  Preserve that diagnostic boundary rather
  # than recycling a non-empty parameter vector into a zero-row ESS result.
  if (!length(ess)) {
    return(data.frame(
      metric = "ESS", parameter = "beta", value = NA_real_,
      source = "coda_effectiveSize_unavailable_for_degenerate_chain",
      stringsAsFactors = FALSE
    ))
  }
  parameter <- names(ess)
  if (is.null(parameter) || length(parameter) != length(ess)) {
    parameter <- colnames(mat)
  }
  data.frame(metric = rep("ESS", length(ess)), parameter = parameter,
             value = as.numeric(ess), source = rep("coda_effectiveSize", length(ess)),
             stringsAsFactors = FALSE)
}

.trace_table <- function(evo, pre) {
  if (!is.null(pre$trace)) return(as.data.frame(pre$trace))
  draws <- evo$beta_draws
  if (is.null(draws)) return(data.frame())
  n <- min(3, dim(draws)[2] * dim(draws)[3])
  mat <- matrix(draws, nrow = dim(draws)[1])
  out <- lapply(seq_len(n), function(i) {
    data.frame(draw = seq_len(nrow(mat)), parameter = paste0("beta_", i),
               value = mat[, i], stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

.standard_metric_table <- function(x, group) {
  if (is.null(x) || nrow(x) == 0) return(data.frame())
  x <- as.data.frame(x)
  if (!"metric" %in% names(x)) x$metric <- group
  if (!"value" %in% names(x)) {
    val <- names(x)[vapply(x, is.numeric, logical(1))][1]
    if (is.na(val)) x$value <- NA_real_ else x$value <- x[[val]]
  }
  x$source_group <- group
  x[, c("metric", "value", "source_group"), drop = FALSE]
}
