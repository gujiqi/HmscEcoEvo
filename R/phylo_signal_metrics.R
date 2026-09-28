#' Calculate phylogenetic signal diagnostics for Beta-derived niche traits
#'
#' Computes global and localized phylogenetic signal summaries for HMSC Beta.
#' HMSC rho can be extracted from the `hmsc_ecoevo` object. Pagel's lambda,
#' Blomberg's K, and Abouheif's Cmean are read from `precomputed` when supplied
#' or delegated to optional comparative packages when available. Descriptive
#' Moran's I, correlogram, scale-dependent signal, LIPA/local Moran's I,
#' node-level signal, clade conservatism index, and residual Beta maps are
#' computed directly when `ape` and a phylogeny are available. When optional
#' `phylosignal` and `phylobase` are installed, Abouheif's Cmean and the
#' correlogram are delegated to those comparative-method implementations.
#'
#' @param evo A `hmsc_ecoevo` object.
#' @param axes Optional subset of Beta axes.
#' @param clade Optional clade vector/table; defaults to `evo$clade`.
#' @param bins Number of phylogenetic-distance bins for correlograms.
#' @param correlogram_permutations Number of label permutations used to obtain
#'   descriptive two-sided p-values for the fallback Moran's I correlogram.
#'   Set to 0 to skip permutation p-values. Precomputed or phylosignal
#'   correlograms are used as supplied.
#' @param precomputed Optional list overriding `evo$precomputed`.
#' @param ... Additional arguments passed to [as_hmsc_ecoevo()].
#' @return A `hmsc_phylo_signal_metrics` object.
#' @export
calc_phylo_signal_metrics <- function(evo,
                                      axes = NULL,
                                      clade = NULL,
                                      bins = 5,
                                      correlogram_permutations = 99,
                                      precomputed = NULL,
                                      ...) {
  evo <- as_hmsc_ecoevo(evo, ...)
  B <- evo$beta
  if (!is.null(axes)) B <- B[, axes, drop = FALSE]
  pre <- utils::modifyList(evo$precomputed, .as_precomputed_list(precomputed))
  diagnostics <- list(messages = character(), warnings = character())
  phylo <- evo$phylo
  clade <- clade %||% evo$clade

  if (is.null(phylo)) {
    diagnostics$warnings <- c(diagnostics$warnings,
      "No phylogeny supplied; phylogenetic signal metrics are returned as precomputed values or NA.")
  }

  rho_table <- .rho_to_table(evo$rho, evo$rho_draws, axes = colnames(B))
  if (is.null(rho_table) && !is.null(pre$rho)) rho_table <- as.data.frame(pre$rho)

  global <- data.frame()
  for (axis in colnames(B)) {
    moran_ci <- if (is.null(phylo)) NULL else
      .phylo_moran_draw_summary(evo$beta_draws, axis, phylo, rownames(B))
    z <- B[, axis]
    moran <- if (is.null(phylo)) NA_real_ else .phylo_moran_i(z, phylo)
    global <- rbind(global, data.frame(
      axis = axis,
      metric = "Moran_I",
      estimate = if (is.null(moran_ci)) moran else moran_ci["mean"],
      lwr = if (is.null(moran_ci)) NA_real_ else moran_ci["lwr"],
      upr = if (is.null(moran_ci)) NA_real_ else moran_ci["upr"],
      source = if (is.null(moran_ci)) "descriptive_inverse_distance" else "beta_draws_inverse_distance",
      stringsAsFactors = FALSE
    ))
  }
  if (!is.null(rho_table) && nrow(rho_table) > 0) {
    rho_out <- rho_table
    rho_out$metric <- "Hmsc_rho"
    rho_out$source <- "Hmsc_or_user"
    rho_out$estimate <- rho_out$mean
    global <- rbind(global, rho_out[, c("axis", "metric", "estimate", "lwr", "upr", "source")])
  }

  pkg_tables <- .optional_phylo_package_metrics(B, phylo, pre, diagnostics)
  global <- rbind(global, pkg_tables$global)
  diagnostics <- pkg_tables$diagnostics

  correlogram <- if (is.null(phylo)) data.frame() else {
    if (!is.null(pre$phylogenetic_correlogram)) {
      as.data.frame(pre$phylogenetic_correlogram)
    } else {
      pc <- .phylosignal_correlogram(B, phylo)
      if (nrow(pc) > 0) pc else .correlogram_table(B, phylo, bins = bins,
                                                   permutations = correlogram_permutations)
    }
  }
  scale_signal <- .scale_signal_table(correlogram)

  lipa <- data.frame()
  if (!is.null(phylo)) {
    for (axis in colnames(B)) {
      li <- .local_moran(B[, axis], phylo)
      if (!is.null(names(li))) li <- li[rownames(B)]
      lipa <- rbind(lipa, data.frame(
        species = rownames(B),
        axis = axis,
        local_moran_i = li,
        local_signal_strength = abs(li),
        stringsAsFactors = FALSE
      ))
    }
  }

  node_signal <- if (is.null(phylo)) data.frame() else .node_signal_table(B, phylo)
  cci <- .clade_cci(B, clade)
  residual_map <- .residual_beta_map(B, traits = evo$traits,
                                     history = evo$history_species)

  out <- list(
    global_signal = global,
    rho = rho_table,
    correlogram = correlogram,
    scale_dependent_signal = scale_signal,
    local_signal = lipa,
    node_signal = node_signal,
    clade_conservatism = cci,
    residual_phylogenetic_map = residual_map,
    diagnostics = diagnostics,
    source = evo
  )
  class(out) <- "hmsc_phylo_signal_metrics"
  out
}

#' Plot phylogenetic signal diagnostics
#'
#' Returns the second HmscEcoEvo diagnostic figure: global phylogenetic signal,
#' phylogenetic correlogram, scale-dependent signal and CCI, local LIPA,
#' node-level signal, and residual Beta map.
#'
#' @param x A `hmsc_ecoevo` or `hmsc_phylo_signal_metrics` object.
#' @param style Figure style. `"nature"` returns a single full-page
#'   Nature-style diagnostic plate; `"classic"` returns the original compact
#'   ggplot/patchwork summary.
#' @param ... Arguments passed to [calc_phylo_signal_metrics()] when needed.
#' @return A `ggplot`, `patchwork`, or printable plot-list object.
#' @export
plot_phylo_signal <- function(x, style = c("nature", "classic"), ...) {
  style <- match.arg(style)
  if (identical(style, "nature")) return(.plot_phylo_signal_nature(x, ...))
  .gg_ok()
  m <- if (inherits(x, "hmsc_phylo_signal_metrics")) x else calc_phylo_signal_metrics(x, ...)
  G <- m$global_signal
  C <- m$correlogram
  S <- m$scale_dependent_signal
  L <- m$local_signal
  N <- m$node_signal
  R <- m$residual_phylogenetic_map
  CCI <- m$clade_conservatism

  p1 <- ggplot2::ggplot(G, ggplot2::aes(estimate, metric, color = axis)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey70") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lwr, xmax = upr),
                            height = 0.15, na.rm = TRUE) +
    ggplot2::geom_point(size = 2, na.rm = TRUE) +
    ggplot2::labs(title = "A. Global phylogenetic signal",
                  x = "Estimate", y = "Metric") +
    ggplot2::theme_minimal(base_size = 10)

  p2 <- if (nrow(C) > 0) {
    p <- ggplot2::ggplot(C, ggplot2::aes(phylo_distance, similarity, color = axis))
    if (all(c("lwr", "upr") %in% names(C))) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lwr, ymax = upr, fill = axis),
        alpha = 0.10, color = NA, show.legend = FALSE, na.rm = TRUE
      )
    }
    p +
      ggplot2::geom_line(linewidth = 0.7, na.rm = TRUE) +
      ggplot2::geom_point(size = 1.8, na.rm = TRUE) +
      ggplot2::labs(title = "B. Phylogenetic correlogram",
                    x = "Phylogenetic distance", y = "Moran's I") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p3a <- if (nrow(S) > 0) {
    ggplot2::ggplot(S, ggplot2::aes(scale, signal, fill = axis)) +
      ggplot2::geom_col(position = "dodge", width = 0.7, na.rm = TRUE) +
      ggplot2::labs(title = "C. Scale-dependent signal",
                    x = "Scale", y = "Signal") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL
  p3b <- if (nrow(CCI) > 0) {
    ggplot2::ggplot(CCI, ggplot2::aes(stats::reorder(clade, CCI), CCI)) +
      ggplot2::geom_col(fill = "#2b6cb0", width = 0.7, na.rm = TRUE) +
      ggplot2::coord_flip() +
      ggplot2::labs(title = "C2. Clade conservatism index",
                    x = "Clade", y = "CCI") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p4 <- if (nrow(L) > 0) {
    ggplot2::ggplot(L, ggplot2::aes(axis, species, fill = local_moran_i)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2) +
      ggplot2::scale_fill_gradient2(low = "#2b6cb0", mid = "white",
                                    high = "#c53030", midpoint = 0,
                                    na.value = "grey85") +
      ggplot2::labs(title = "D. Tip-level local Moran's I",
                    x = "Axis", y = "Species", fill = "Local I") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p5 <- if (nrow(N) > 0) {
    ggplot2::ggplot(N, ggplot2::aes(node, node_signal, color = axis, size = n_species)) +
      ggplot2::geom_point(alpha = 0.8, na.rm = TRUE) +
      ggplot2::labs(title = "E. Internal node-level signal",
                    x = "Node", y = "Signal", size = "Tips") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p6 <- if (nrow(R) > 0) {
    ggplot2::ggplot(R, ggplot2::aes(axis, species, fill = residual_beta)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2) +
      ggplot2::scale_fill_gradient2(low = "#2b6cb0", mid = "white",
                                    high = "#c53030", midpoint = 0,
                                    na.value = "grey85") +
      ggplot2::labs(title = "F. Residual Beta map",
                    x = "Axis", y = "Species", fill = "Residual") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  .combine_plots(list(p1, p2, p3a, p3b, p4, p5, p6),
                 title = "Phylogenetic Signal Localization", ncol = 2)
}

.rho_to_table <- function(rho, rho_draws = NULL, axes) {
  if (!is.null(rho_draws)) {
    if (is.matrix(rho_draws) || is.data.frame(rho_draws)) {
      d <- as.data.frame(rho_draws)
      if (is.null(names(d))) names(d) <- axes[seq_len(ncol(d))]
      out <- lapply(names(d), function(ax) {
        ci <- .posterior_ci_vector(d[[ax]])
        data.frame(axis = ax, mean = ci["mean"], median = ci["median"],
                   lwr = ci["lwr"], upr = ci["upr"], sd = ci["sd"],
                   support = ci["support"], stringsAsFactors = FALSE)
      })
      return(do.call(rbind, out))
    }
    if (is.vector(rho_draws)) {
      ci <- .posterior_ci_vector(rho_draws)
      return(data.frame(axis = "rho", mean = ci["mean"], median = ci["median"],
                        lwr = ci["lwr"], upr = ci["upr"], sd = ci["sd"],
                        support = ci["support"], stringsAsFactors = FALSE))
    }
  }
  if (is.null(rho)) return(NULL)
  if (is.data.frame(rho)) {
    if (!"axis" %in% names(rho)) rho$axis <- axes[seq_len(nrow(rho))]
    if (!"mean" %in% names(rho) && "estimate" %in% names(rho)) rho$mean <- rho$estimate
    if (!"lwr" %in% names(rho)) rho$lwr <- NA_real_
    if (!"upr" %in% names(rho)) rho$upr <- NA_real_
    return(rho)
  }
  if (is.matrix(rho)) {
    vals <- as.numeric(rho)
    return(data.frame(axis = rep(axes, length.out = length(vals)),
                      mean = vals, lwr = NA_real_, upr = NA_real_,
                      stringsAsFactors = FALSE))
  }
  vals <- as.numeric(rho)
  data.frame(axis = rep(axes, length.out = length(vals)),
             mean = vals, lwr = NA_real_, upr = NA_real_,
             stringsAsFactors = FALSE)
}

.optional_phylo_package_metrics <- function(B, phylo, pre, diagnostics) {
  global <- data.frame()
  add_pre <- function(tab, metric, source) {
    if (is.null(tab)) return(data.frame())
    tab <- as.data.frame(tab)
    if (!"axis" %in% names(tab)) tab$axis <- colnames(B)[seq_len(nrow(tab))]
    if (!"estimate" %in% names(tab)) {
      val <- setdiff(names(tab), c("axis", "lwr", "upr", "source"))[1]
      tab$estimate <- tab[[val]]
    }
    if (!"lwr" %in% names(tab)) tab$lwr <- NA_real_
    if (!"upr" %in% names(tab)) tab$upr <- NA_real_
    tab$metric <- metric
    tab$source <- source
    tab[, c("axis", "metric", "estimate", "lwr", "upr", "source")]
  }
  global <- rbind(global, add_pre(pre$pagel_lambda, "Pagel_lambda", "precomputed"))
  global <- rbind(global, add_pre(pre$blomberg_k, "Blomberg_K", "precomputed"))
  global <- rbind(global, add_pre(pre$abouheif_cmean, "Abouheif_Cmean", "precomputed"))

  if (!is.null(phylo) && is.null(pre$pagel_lambda) && .soft_pkg("phytools")) {
    for (axis in colnames(B)) {
      val <- tryCatch(phytools::phylosig(phylo, B[phylo$tip.label, axis], method = "lambda"),
                      error = function(e) NA_real_)
      est <- .phylosig_estimate(val, "lambda")
      global <- rbind(global, data.frame(axis = axis, metric = "Pagel_lambda",
        estimate = est, lwr = NA_real_, upr = NA_real_,
        source = "phytools::phylosig", stringsAsFactors = FALSE))
    }
  } else if (is.null(pre$pagel_lambda)) {
    diagnostics$messages <- c(diagnostics$messages,
      "Pagel's lambda requires precomputed$pagel_lambda or the optional phytools package.")
    global <- rbind(global, data.frame(axis = colnames(B), metric = "Pagel_lambda",
      estimate = NA_real_, lwr = NA_real_, upr = NA_real_,
      source = "requires_precomputed_or_phytools", stringsAsFactors = FALSE))
  }

  if (!is.null(phylo) && is.null(pre$blomberg_k) && .soft_pkg("phytools")) {
    for (axis in colnames(B)) {
      val <- tryCatch(phytools::phylosig(phylo, B[phylo$tip.label, axis], method = "K"),
                      error = function(e) NA_real_)
      est <- .phylosig_estimate(val, "K")
      global <- rbind(global, data.frame(axis = axis, metric = "Blomberg_K",
        estimate = est, lwr = NA_real_, upr = NA_real_,
        source = "phytools::phylosig", stringsAsFactors = FALSE))
    }
  } else if (is.null(pre$blomberg_k)) {
    diagnostics$messages <- c(diagnostics$messages,
      "Blomberg's K requires precomputed$blomberg_k or the optional phytools package.")
    global <- rbind(global, data.frame(axis = colnames(B), metric = "Blomberg_K",
      estimate = NA_real_, lwr = NA_real_, upr = NA_real_,
      source = "requires_precomputed_or_phytools", stringsAsFactors = FALSE))
  }

  if (is.null(pre$abouheif_cmean)) {
    cmean <- .phylosignal_global_metric(B, phylo, method = "Cmean")
    if (nrow(cmean) > 0) {
      global <- rbind(global, cmean)
    } else {
      diagnostics$messages <- c(diagnostics$messages,
        "Abouheif's Cmean requires precomputed$abouheif_cmean or optional packages phylosignal and phylobase.")
      global <- rbind(global, data.frame(axis = colnames(B), metric = "Abouheif_Cmean",
        estimate = NA_real_, lwr = NA_real_, upr = NA_real_,
        source = "requires_precomputed_or_phylosignal", stringsAsFactors = FALSE))
    }
  }
  list(global = global, diagnostics = diagnostics)
}

.phylosig_estimate <- function(x, field) {
  if (is.null(x)) return(NA_real_)
  if (is.atomic(x) && length(x) > 0L) return(as.numeric(x[1]))
  if (is.list(x)) {
    candidates <- unique(c(field, toupper(field), tolower(field), "estimate", "value"))
    for (nm in candidates) {
      if (!is.null(x[[nm]]) && is.atomic(x[[nm]]) && length(x[[nm]]) > 0L) {
        return(as.numeric(x[[nm]][1]))
      }
    }
  }
  NA_real_
}

.make_phylo4d <- function(B, phylo) {
  if (is.null(phylo) || !.soft_pkg("phylobase")) return(NULL)
  if (!all(rownames(B) %in% phylo$tip.label)) return(NULL)
  p4 <- tryCatch(phylobase::phylo4(phylo), error = function(e) NULL)
  if (is.null(p4)) return(NULL)
  data <- as.data.frame(B[phylo$tip.label, , drop = FALSE])
  tryCatch(phylobase::phylo4d(p4, tip.data = data), error = function(e) NULL)
}

.phylosignal_global_metric <- function(B, phylo, method = "Cmean") {
  if (is.null(phylo) || !.soft_pkg("phylosignal") || !.soft_pkg("phylobase")) {
    return(data.frame())
  }
  p4d <- .make_phylo4d(B, phylo)
  if (is.null(p4d)) return(data.frame())
  sig <- tryCatch(phylosignal::phyloSignal(p4d, methods = method, reps = 0),
                  error = function(e) NULL)
  if (is.null(sig) || is.null(sig$stat) || !method %in% names(sig$stat)) {
    return(data.frame())
  }
  data.frame(
    axis = rownames(sig$stat),
    metric = if (identical(method, "Cmean")) "Abouheif_Cmean" else method,
    estimate = as.numeric(sig$stat[[method]]),
    lwr = NA_real_,
    upr = NA_real_,
    source = paste0("phylosignal::phyloSignal:", method),
    stringsAsFactors = FALSE
  )
}

.phylosignal_correlogram <- function(B, phylo, n.points = 25, ci.bs = 100) {
  if (is.null(phylo) || !.soft_pkg("phylosignal") || !.soft_pkg("phylobase")) {
    return(data.frame())
  }
  p4d <- .make_phylo4d(B, phylo)
  if (is.null(p4d)) return(data.frame())
  out <- list()
  k <- 1L
  for (axis in colnames(B)) {
    pc <- tryCatch(phylosignal::phyloCorrelogram(
      p4d, trait = axis, dist.phylo = "patristic",
      n.points = n.points, ci.bs = ci.bs
    ), error = function(e) NULL)
    if (is.null(pc) || is.null(pc$res) || ncol(pc$res) < 4) next
    out[[k]] <- data.frame(
      distance_bin = seq_len(nrow(pc$res)),
      phylo_distance = pc$res[, 1],
      lwr = pc$res[, 2],
      upr = pc$res[, 3],
      moran_i = pc$res[, 4],
      beta_distance = NA_real_,
      axis = axis,
      similarity = pc$res[, 4],
      source = "phylosignal::phyloCorrelogram",
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
  if (length(out) == 0) data.frame() else do.call(rbind, out)
}

.scale_signal_table <- function(correlogram) {
  if (is.null(correlogram) || nrow(correlogram) == 0) return(data.frame())
  out <- lapply(split(correlogram, correlogram$axis), function(d) {
    ord <- order(d$phylo_distance)
    d <- d[ord, , drop = FALSE]
    n <- nrow(d)
    scale <- rep(c("shallow", "intermediate", "deep"), length.out = n)
    data.frame(axis = d$axis, scale = scale, signal = d$similarity,
               phylo_distance = d$phylo_distance, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

.phylo_moran_draw_summary <- function(beta_draws, axis, phylo, species) {
  if (is.null(beta_draws) || is.null(phylo) || !axis %in% dimnames(beta_draws)[[3]]) {
    return(NULL)
  }
  species <- species[species %in% dimnames(beta_draws)[[2]]]
  if (length(species) < 3) return(NULL)
  vals <- vapply(seq_len(dim(beta_draws)[1]), function(i) {
    z <- beta_draws[i, species, axis]
    names(z) <- species
    .phylo_moran_i(z, phylo)
  }, numeric(1))
  if (all(is.na(vals))) return(NULL)
  .posterior_ci_vector(vals)
}

.residual_beta_map <- function(B, traits = NULL, history = NULL) {
  predictors <- NULL
  if (!is.null(traits)) predictors <- traits
  if (!is.null(history)) predictors <- if (is.null(predictors)) history else cbind(predictors, history)
  out <- list()
  k <- 1L
  for (axis in colnames(B)) {
    y <- B[, axis]
    residual <- y - mean(y, na.rm = TRUE)
    if (!is.null(predictors) && ncol(predictors) > 0) {
      d <- .align_species_frame(predictors, rownames(B), "predictors")
      keep <- vapply(d, function(z) length(unique(z[!is.na(z)])) > 1, logical(1))
      d <- d[, keep, drop = FALSE]
      if (ncol(d) > 0) {
        d <- d[vapply(d, is.numeric, logical(1))]
      }
      if (ncol(d) > 0) {
        for (nm in names(d)) {
          if (anyNA(d[[nm]])) d[[nm]][is.na(d[[nm]])] <- mean(d[[nm]], na.rm = TRUE)
        }
        n_obs <- sum(!is.na(y))
        if (ncol(d) >= max(2, n_obs - 3)) {
          pc <- tryCatch(stats::prcomp(d, center = TRUE, scale. = TRUE),
                         error = function(e) NULL)
          if (!is.null(pc)) {
            n_pc <- min(3, ncol(pc$x), max(1, floor((n_obs - 2) / 2)))
            d <- as.data.frame(pc$x[, seq_len(n_pc), drop = FALSE])
            names(d) <- paste0("predictor_PC", seq_len(ncol(d)))
          }
        }
        fit <- tryCatch(stats::lm(y ~ ., data = d), error = function(e) NULL)
        if (!is.null(fit)) residual <- stats::residuals(fit)
      }
    }
    out[[k]] <- data.frame(species = names(residual), axis = axis,
                           residual_beta = as.numeric(residual),
                           stringsAsFactors = FALSE)
    k <- k + 1L
  }
  do.call(rbind, out)
}
