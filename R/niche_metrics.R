#' Calculate niche metrics from HMSC Beta
#'
#' Computes the foundational Beta-derived niche metrics: coefficient table,
#' posterior support and uncertainty, response curves, curvature-based optimum
#' and Gaussian-equivalent breadth, vector magnitude and direction, pairwise
#' cosine similarity and angle, and Beta-PCA. Optimum and breadth are reported
#' only for valid concave quadratic responses. For `eta = b0 + b1*x + b2*x^2`,
#' the optimum is `-b1/(2*b2)` and the Gaussian-equivalent breadth is
#' `sqrt(-1/(2*b2))`.
#'
#' @param evo A `hmsc_ecoevo` object or arguments coercible with
#'   [as_hmsc_ecoevo()].
#' @param quadratic Optional data.frame with columns `linear` and `quadratic`
#'   naming first- and second-order terms. If omitted, common quadratic names
#'   such as `temp_sq`, `temp2`, and `I(temp^2)` are detected.
#' @param x_grid Numeric grid used for response curves.
#' @param axes Optional subset of niche axes.
#' @param weights Optional named axis weights for Beta-vector magnitude.
#' @param distr Optional model distribution/link label. Defaults to
#'   `evo$distr`.
#' @param include_intercept Logical; include intercept-like columns in vector
#'   magnitude/PCA summaries. Defaults to `FALSE`.
#' @param intercept_axis Optional name of an intercept column. If omitted,
#'   common names such as `(Intercept)` and `intercept` are detected.
#' @param ... Additional arguments passed to [as_hmsc_ecoevo()] if needed.
#' @return A `hmsc_niche_metrics` object.
#' @export
#' @examples
#' beta <- matrix(c(0.3, -0.2, 0.1, 0.4), nrow = 2,
#'   dimnames = list(c("sp1", "sp2"), c("temp", "precip")))
#' evo <- hmsc_ecoevo(beta = beta)
#' calc_niche_metrics(evo)
calc_niche_metrics <- function(evo,
                               quadratic = NULL,
                               x_grid = seq(-2, 2, length.out = 101),
                               axes = NULL,
                               weights = NULL,
                               distr = NULL,
                               include_intercept = FALSE,
                               intercept_axis = NULL,
                               ...) {
  evo <- as_hmsc_ecoevo(evo, ...)
  B_all <- evo$beta
  all_draws <- evo$beta_draws
  if (is.null(intercept_axis)) {
    hit <- colnames(B_all)[.is_intercept_axis(colnames(B_all))]
    intercept_axis <- if (length(hit) > 0) hit[1] else NULL
  }
  if (is.null(axes)) {
    axes <- if (include_intercept) colnames(B_all) else .non_intercept_axes(colnames(B_all))
  }
  if (length(axes) == 0) stop("No niche axes remain after intercept filtering.", call. = FALSE)
  B <- B_all
  if (!is.null(axes)) {
    missing <- setdiff(axes, colnames(B_all))
    if (length(missing) > 0) stop("Unknown niche axis/axes: ", paste(missing, collapse = ", "), call. = FALSE)
    B <- B_all[, axes, drop = FALSE]
  }
  draws <- all_draws
  if (!is.null(draws)) draws <- draws[, rownames(B), colnames(B), drop = FALSE]
  distr <- distr %||% evo$distr
  diagnostics <- list(messages = character(), warnings = character())

  beta_table <- .matrix_to_long(B, value_name = "beta",
                                row_name = "species", col_name = "axis")
  beta_summary <- if (!is.null(draws)) {
    .posterior_summary_array(draws, value_name = "beta")
  } else {
    .support_from_matrix(B, value_name = "beta",
                         row_name = "species", col_name = "axis")
  }

  if (is.null(weights)) {
    w <- rep(1, ncol(B))
    names(w) <- colnames(B)
  } else {
    w <- as.numeric(weights[colnames(B)])
    names(w) <- colnames(B)
    w[is.na(w)] <- 1
  }
  magnitude <- sqrt(rowSums(sweep(B^2, 2, w, `*`), na.rm = TRUE))
  direction <- B
  nz <- magnitude > 0 & !is.na(magnitude)
  direction[nz, ] <- direction[nz, , drop = FALSE] / magnitude[nz]
  direction[!nz, ] <- NA_real_
  cosine <- tcrossprod(direction)
  cosine <- pmax(pmin(cosine, 1), -1)
  angle <- acos(cosine)

  magnitude_table <- data.frame(
    species = rownames(B),
    beta_magnitude = magnitude,
    mean_axis_uncertainty = NA_real_,
    stringsAsFactors = FALSE
  )
  if (!is.null(draws)) {
    mag_draws <- apply(draws, 1, function(m) {
      mm <- matrix(m, nrow = nrow(B), ncol = ncol(B),
                   dimnames = list(rownames(B), colnames(B)))
      sqrt(rowSums(sweep(mm^2, 2, w, `*`), na.rm = TRUE))
    })
    if (is.vector(mag_draws)) mag_draws <- matrix(mag_draws, nrow = nrow(B))
    magnitude_table$mean <- rowMeans(mag_draws, na.rm = TRUE)
    magnitude_table$median <- apply(mag_draws, 1, stats::median, na.rm = TRUE)
    magnitude_table$lwr <- apply(mag_draws, 1, stats::quantile, 0.025, na.rm = TRUE)
    magnitude_table$upr <- apply(mag_draws, 1, stats::quantile, 0.975, na.rm = TRUE)
    magnitude_table$sd <- apply(mag_draws, 1, stats::sd, na.rm = TRUE)
    magnitude_table$mean_axis_uncertainty <- stats::aggregate(
      beta_summary$sd, list(species = beta_summary$species), mean, na.rm = TRUE
    )$x
  }

  pairs <- .axis_pairs(colnames(B), quadratic = quadratic)
  optimum <- data.frame()
  response_curves <- data.frame()
  used_linear <- character()
  if (nrow(pairs) > 0) {
    for (i in seq_len(nrow(pairs))) {
      lin <- pairs$linear[i]
      quad <- pairs$quadratic[i]
      if (!all(c(lin, quad) %in% colnames(B))) next
      used_linear <- c(used_linear, lin, quad)
      b1 <- B[, lin]
      b2 <- B[, quad]
      valid <- b2 < 0
      opt <- ifelse(valid, -b1 / (2 * b2), NA_real_)
      breadth <- ifelse(valid, sqrt(-1 / (2 * b2)), NA_real_)
      optimum <- rbind(optimum, data.frame(
        species = rownames(B),
        response_axis = pairs$axis[i],
        linear_axis = lin,
        quadratic_axis = quad,
        beta_linear = b1,
        beta_quadratic = b2,
        optimum = opt,
        optimum_type = ifelse(valid, "maximum", "not_defined"),
        breadth_index = breadth,
        stringsAsFactors = FALSE
      ))
      response_curves <- rbind(response_curves, .build_response_curves(
        B_all = B_all, draws = all_draws, species = rownames(B), x_grid = x_grid,
        response_axis = pairs$axis[i], linear_axis = lin, quadratic_axis = quad,
        intercept_axis = intercept_axis, distr = distr, curve_type = "quadratic"
      ))
    }
  }
  linear_axes <- setdiff(colnames(B), used_linear)
  if (length(linear_axes) == 0 && nrow(pairs) == 0) {
    diagnostics$messages <- c(diagnostics$messages,
      "No quadratic axis pairs were detected; optimum and breadth are unavailable.")
  }
  if (length(linear_axes) > 0) {
    if (nrow(pairs) == 0) diagnostics$messages <- c(diagnostics$messages,
      "No quadratic axis pairs were detected; optimum and breadth are unavailable.")
    for (axis in linear_axes) {
      response_curves <- rbind(response_curves, .build_response_curves(
        B_all = B_all, draws = all_draws, species = rownames(B), x_grid = x_grid,
        response_axis = axis, linear_axis = axis, quadratic_axis = NULL,
        intercept_axis = intercept_axis, distr = distr, curve_type = "linear"
      ))
    }
  }

  pca <- NULL
  pca_scores <- data.frame()
  pca_loadings <- data.frame()
  pca_variance <- data.frame()
  if (nrow(B) >= 2 && ncol(B) >= 2) {
    pca <- tryCatch(stats::prcomp(B, center = TRUE, scale. = TRUE),
                    error = function(e) NULL)
    if (!is.null(pca)) {
      pca_scores <- data.frame(species = rownames(B), pca$x, check.names = FALSE)
      pca_loadings <- data.frame(axis = rownames(pca$rotation), pca$rotation,
                                 check.names = FALSE)
      pca_variance <- data.frame(
        component = paste0("PC", seq_along(pca$sdev)),
        variance = pca$sdev^2,
        proportion = pca$sdev^2 / sum(pca$sdev^2),
        stringsAsFactors = FALSE
      )
    }
  } else {
    diagnostics$messages <- c(diagnostics$messages,
      "Beta-PCA needs at least two species and two niche axes.")
  }

  out <- list(
    beta = B,
    beta_table = beta_table,
    beta_summary = beta_summary,
    response_curves = response_curves,
    optimum_breadth = optimum,
    magnitude = magnitude_table,
    direction = direction,
    cosine_similarity = cosine,
    angle = angle,
    pca = pca,
    pca_scores = pca_scores,
    pca_loadings = pca_loadings,
    pca_variance = pca_variance,
    quadratic_pairs = pairs,
    diagnostics = diagnostics,
    source = evo
  )
  class(out) <- "hmsc_niche_metrics"
  out
}

#' Plot niche summary diagnostics
#'
#' Builds the first HmscEcoEvo diagnostic figure: Beta heatmap, posterior
#' support/uncertainty, response curves, optimum-breadth scatter, Beta-PCA, and
#' vector magnitude/uncertainty. The function returns a `patchwork` object when
#' `patchwork` is installed, otherwise a printable plot-list object.
#'
#' @param x A `hmsc_ecoevo` or `hmsc_niche_metrics` object.
#' @param style Figure style. `"nature"` returns a single full-page
#'   Nature-style diagnostic plate; `"classic"` returns the original compact
#'   ggplot/patchwork summary.
#' @param ... Arguments passed to [calc_niche_metrics()] when `x` is not already
#'   a metrics object.
#' @return A `ggplot`, `patchwork`, or printable plot-list object.
#' @export
plot_niche_summary <- function(x, style = c("nature", "classic"), ...) {
  style <- match.arg(style)
  if (identical(style, "nature")) return(.plot_niche_summary_nature(x, ...))
  .gg_ok()
  metrics <- if (inherits(x, "hmsc_niche_metrics")) x else calc_niche_metrics(x, ...)
  B <- metrics$beta_table
  S <- metrics$beta_summary
  P <- metrics$pca_scores
  L <- metrics$pca_loadings
  M <- metrics$magnitude
  O <- metrics$optimum_breadth
  R <- metrics$response_curves

  p1 <- ggplot2::ggplot(B, ggplot2::aes(axis, species, fill = beta)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(low = "#2b6cb0", mid = "white",
                                  high = "#c53030", midpoint = 0,
                                  na.value = "grey85") +
    ggplot2::labs(title = "A. Species-by-environment Beta",
                  x = "Niche axis", y = "Species", fill = "Beta") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  p2 <- ggplot2::ggplot(S, ggplot2::aes(beta, species, color = axis)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lwr, xmax = upr),
                            height = 0.15, na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(alpha = support), size = 1.8) +
    ggplot2::scale_alpha_continuous(limits = c(0, 1), na.value = 0.45) +
    ggplot2::labs(title = "B. Beta support and uncertainty",
                  x = "Posterior mean", y = "Species", alpha = "Support") +
    ggplot2::theme_minimal(base_size = 10)

  p3 <- if (nrow(R) > 0) {
    p <- ggplot2::ggplot(R, ggplot2::aes(x, response, group = species, color = species))
    if (all(c("response_lwr", "response_upr") %in% names(R))) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = response_lwr, ymax = response_upr, fill = species),
        alpha = 0.10, color = NA, show.legend = FALSE, na.rm = TRUE
      )
    }
    p +
      ggplot2::geom_line(alpha = 0.8, linewidth = 0.6) +
      ggplot2::facet_wrap(~response_axis, scales = "free_y") +
      ggplot2::labs(title = "C. Niche response curves",
                    x = "Standardized environmental gradient", y = "Expected response") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(legend.position = "none")
  } else NULL

  p4 <- if (nrow(O) > 0) {
    ggplot2::ggplot(O, ggplot2::aes(optimum, breadth_index, color = response_axis)) +
      ggplot2::geom_point(ggplot2::aes(shape = optimum_type), size = 2, na.rm = TRUE) +
      ggplot2::labs(title = "D. Niche optimum vs breadth",
                    x = "Optimum", y = "Curvature-based breadth") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p5 <- if (nrow(P) > 0 && all(c("PC1", "PC2") %in% names(P))) {
    p <- ggplot2::ggplot(P, ggplot2::aes(PC1, PC2, label = species)) +
      ggplot2::geom_point(size = 2, color = "#2f855a") +
      ggplot2::geom_text(vjust = -0.6, size = 2.6, check_overlap = TRUE) +
      ggplot2::labs(title = "E. Multivariate Beta-PCA",
                    x = "PC1", y = "PC2") +
      ggplot2::theme_minimal(base_size = 10)
    if (nrow(L) > 0 && all(c("PC1", "PC2") %in% names(L))) {
      scale <- max(abs(P[, c("PC1", "PC2")]), na.rm = TRUE)
      p <- p + ggplot2::geom_segment(data = L,
          ggplot2::aes(x = 0, y = 0, xend = PC1 * scale, yend = PC2 * scale),
          inherit.aes = FALSE, arrow = ggplot2::arrow(length = grid::unit(0.12, "in")),
          color = "#4a5568") +
        ggplot2::geom_text(data = L,
          ggplot2::aes(x = PC1 * scale, y = PC2 * scale, label = axis),
          inherit.aes = FALSE, color = "#4a5568", size = 2.5)
    }
    p
  } else NULL

  p6 <- ggplot2::ggplot(M, ggplot2::aes(stats::reorder(species, beta_magnitude), beta_magnitude)) +
    ggplot2::geom_col(fill = "#805ad5", width = 0.7) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lwr, ymax = upr), width = 0.15,
                           na.rm = TRUE, color = "grey30") +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "F. Beta-vector magnitude",
                  x = "Species", y = "Magnitude") +
    ggplot2::theme_minimal(base_size = 10)

  .combine_plots(list(p1, p2, p3, p4, p5, p6),
                 title = "Niche Summary", ncol = 2)
}

.build_response_curves <- function(B_all, draws = NULL, species, x_grid,
                                   response_axis, linear_axis,
                                   quadratic_axis = NULL,
                                   intercept_axis = NULL,
                                   distr = NULL,
                                   curve_type = "linear") {
  B_all <- as.matrix(B_all)
  species <- intersect(species, rownames(B_all))
  has_axis <- function(axis) !is.null(axis) && axis %in% colnames(B_all)
  draw_names <- if (!is.null(draws)) dimnames(draws) else NULL
  can_use_draws <- !is.null(draws) &&
    !is.null(draw_names[[2]]) && !is.null(draw_names[[3]]) &&
    linear_axis %in% draw_names[[3]]
  out <- vector("list", length(species))
  for (i in seq_along(species)) {
    sp <- species[i]
    b0 <- if (has_axis(intercept_axis)) B_all[sp, intercept_axis] else 0
    b1 <- B_all[sp, linear_axis]
    b2 <- if (has_axis(quadratic_axis)) B_all[sp, quadratic_axis] else 0
    eta <- b0 + b1 * x_grid + b2 * x_grid^2
    response <- .link_response(eta, distr = distr)
    tab <- data.frame(
      species = sp,
      response_axis = response_axis,
      linear_axis = linear_axis,
      quadratic_axis = quadratic_axis %||% NA_character_,
      x = x_grid,
      eta = eta,
      response = response,
      eta_lwr = NA_real_,
      eta_upr = NA_real_,
      response_lwr = NA_real_,
      response_upr = NA_real_,
      curve_type = curve_type,
      stringsAsFactors = FALSE
    )
    if (can_use_draws && sp %in% draw_names[[2]]) {
      b0d <- if (!is.null(intercept_axis) && intercept_axis %in% draw_names[[3]]) {
        draws[, sp, intercept_axis]
      } else {
        rep(0, dim(draws)[1])
      }
      b1d <- draws[, sp, linear_axis]
      b2d <- if (!is.null(quadratic_axis) && quadratic_axis %in% draw_names[[3]]) {
        draws[, sp, quadratic_axis]
      } else {
        rep(0, dim(draws)[1])
      }
      eta_draws <- b0d + outer(b1d, x_grid, `*`) + outer(b2d, x_grid^2, `*`)
      q_eta <- apply(eta_draws, 2, stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE)
      response_draws <- .link_response(eta_draws, distr = distr)
      q_response <- apply(response_draws, 2, stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE)
      tab$eta_lwr <- q_eta[1, ]
      tab$eta_upr <- q_eta[2, ]
      tab$response_lwr <- q_response[1, ]
      tab$response_upr <- q_response[2, ]
    }
    out[[i]] <- tab
  }
  if (length(out) == 0) data.frame() else do.call(rbind, out)
}
