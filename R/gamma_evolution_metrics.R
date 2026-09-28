#' Calculate Gamma / trait-function evolution diagnostics
#'
#' Computes diagnostics for evolution in the trait-to-niche relationship
#' represented by Gamma. Global Gamma is taken from HMSC when available.
#' Clade- and regime-specific Gamma are estimated descriptively as within-group
#' regressions of Beta on traits unless precomputed tables are supplied. Gamma
#' shift probabilities require external model output and are not inferred by
#' HmscEcoEvo.
#'
#' @param evo A `hmsc_ecoevo` object.
#' @param axes Optional subset of Beta axes.
#' @param clade Optional clade vector/table; defaults to `evo$clade`.
#' @param regime Optional regime vector/table; defaults to `evo$regime`.
#' @param network_threshold Minimum absolute Gamma value for network edges.
#' @param precomputed Optional list containing `clade_gamma`, `regime_gamma`,
#'   and/or `gamma_shift_probability`.
#' @param ... Additional arguments passed to [as_hmsc_ecoevo()].
#' @return A `hmsc_gamma_evolution_metrics` object.
#' @export
calc_gamma_evolution_metrics <- function(evo,
                                         axes = NULL,
                                         clade = NULL,
                                         regime = NULL,
                                         network_threshold = 0,
                                         precomputed = NULL,
                                         ...) {
  evo <- as_hmsc_ecoevo(evo, ...)
  B <- evo$beta
  if (!is.null(axes)) B <- B[, axes, drop = FALSE]
  pre <- utils::modifyList(evo$precomputed, .as_precomputed_list(precomputed))
  clade <- clade %||% evo$clade
  regime <- regime %||% evo$regime
  diagnostics <- list(messages = character(), warnings = character())

  gamma <- .gamma_summary(evo$gamma, evo$gamma_draws)
  if (nrow(gamma) > 0 && "axis" %in% names(gamma)) {
    gamma <- gamma[gamma$axis %in% colnames(B), , drop = FALSE]
  }
  traits <- evo$traits

  clade_gamma <- if (!is.null(pre$clade_gamma)) {
    as.data.frame(pre$clade_gamma)
  } else .group_gamma(B, traits, clade, "clade")
  if (nrow(clade_gamma) == 0) diagnostics$messages <- c(diagnostics$messages,
    "Clade-specific Gamma requires traits and a clade vector, or precomputed$clade_gamma.")

  regime_gamma <- if (!is.null(pre$regime_gamma)) {
    as.data.frame(pre$regime_gamma)
  } else .group_gamma(B, traits, regime, "regime")
  if (nrow(regime_gamma) == 0) diagnostics$messages <- c(diagnostics$messages,
    "Regime-specific Gamma requires traits and a regime vector, or precomputed$regime_gamma.")

  shift_prob <- if (!is.null(pre$gamma_shift_probability)) {
    as.data.frame(pre$gamma_shift_probability)
  } else {
    diagnostics$messages <- c(diagnostics$messages,
      "Gamma-shift probability requires precomputed$gamma_shift_probability from an external model.")
    .gamma_shift_placeholder(gamma)
  }
  branch_shift_prob <- if (!is.null(pre$gamma_branch_shift_probability)) {
    as.data.frame(pre$gamma_branch_shift_probability)
  } else if ("branch" %in% names(shift_prob)) {
    stats::aggregate(gamma_shift_probability ~ branch, shift_prob, mean, na.rm = TRUE)
  } else {
    data.frame(branch = integer(), gamma_shift_probability = numeric(),
               source = character())
  }

  sign_flip <- .gamma_sign_flip(gamma, clade_gamma)
  gamma_var <- .gamma_variance_among_groups(clade_gamma)
  turnover <- .gamma_turnover(clade_gamma)
  specialization <- .trait_axis_specialization(gamma)
  network <- .trait_function_network(gamma, threshold = network_threshold)

  out <- list(
    gamma = gamma,
    clade_specific_gamma = clade_gamma,
    regime_specific_gamma = regime_gamma,
    gamma_shift_probability = shift_prob,
    gamma_branch_shift_probability = branch_shift_prob,
    gamma_sign_flip_index = sign_flip,
    gamma_variance_among_clades = gamma_var,
    trait_function_turnover_index = turnover$index,
    pairwise_turnover_matrix = turnover$matrix,
    trait_axis_specialization = specialization,
    trait_function_network = network,
    diagnostics = diagnostics,
    source = evo
  )
  class(out) <- "hmsc_gamma_evolution_metrics"
  out
}

#' Plot Gamma evolution diagnostics
#'
#' Returns the fifth HmscEcoEvo diagnostic figure: clade-specific Gamma,
#' regime-specific Gamma, Gamma-shift probability, sign-flip and variance
#' diagnostics, trait-function network, and pairwise turnover matrix.
#'
#' @param x A `hmsc_ecoevo` or `hmsc_gamma_evolution_metrics` object.
#' @param style Figure style. `"nature"` returns a single full-page
#'   Nature-style diagnostic plate; `"classic"` returns the original compact
#'   ggplot/patchwork summary.
#' @param ... Arguments passed to [calc_gamma_evolution_metrics()] when needed.
#' @return A `ggplot`, `patchwork`, or printable plot-list object.
#' @export
plot_gamma_evolution <- function(x, style = c("nature", "classic"), ...) {
  style <- match.arg(style)
  if (identical(style, "nature")) return(.plot_gamma_evolution_nature(x, ...))
  .gg_ok()
  m <- if (inherits(x, "hmsc_gamma_evolution_metrics")) x else calc_gamma_evolution_metrics(x, ...)
  C <- m$clade_specific_gamma
  R <- m$regime_specific_gamma
  S <- m$gamma_shift_probability
  F <- m$gamma_sign_flip_index
  V <- m$gamma_variance_among_clades
  N <- m$trait_function_network
  T <- m$trait_function_turnover_index
  A <- m$trait_axis_specialization

  p1 <- if (nrow(C) > 0) {
    ggplot2::ggplot(C, ggplot2::aes(axis, trait, fill = gamma)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2) +
      ggplot2::facet_wrap(~clade) +
      ggplot2::scale_fill_gradient2(low = "#2b6cb0", mid = "white",
                                    high = "#c53030", midpoint = 0) +
      ggplot2::labs(title = "A. Clade-specific Gamma",
                    x = "Axis", y = "Trait", fill = "Gamma") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p2 <- if (nrow(R) > 0) {
    ggplot2::ggplot(R, ggplot2::aes(axis, trait, fill = gamma)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2) +
      ggplot2::facet_wrap(~regime) +
      ggplot2::scale_fill_gradient2(low = "#2b6cb0", mid = "white",
                                    high = "#c53030", midpoint = 0) +
      ggplot2::labs(title = "B. Regime-specific Gamma",
                    x = "Axis", y = "Trait", fill = "Gamma") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p3 <- if (nrow(S) > 0) {
    ggplot2::ggplot(S, ggplot2::aes(axis, trait, fill = gamma_shift_probability)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2) +
      ggplot2::scale_fill_gradient(low = "white", high = "#c53030", na.value = "grey85") +
      ggplot2::labs(title = "C. Gamma-shift probability",
                    x = "Axis", y = "Trait", fill = "P shift") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p4 <- if (nrow(F) > 0) {
    ggplot2::ggplot(F, ggplot2::aes(axis, gamma_sign_flip_index, fill = trait)) +
      ggplot2::geom_col(position = "dodge", width = 0.7, na.rm = TRUE) +
      ggplot2::labs(title = "D. Gamma sign-flip index",
                    x = "Axis", y = "Sign-flip index") +
      ggplot2::theme_minimal(base_size = 10)
  } else if (nrow(V) > 0) {
    ggplot2::ggplot(V, ggplot2::aes(axis, gamma_variance_among_clades, fill = trait)) +
      ggplot2::geom_col(position = "dodge", width = 0.7, na.rm = TRUE) +
      ggplot2::labs(title = "D. Gamma variance among clades",
                    x = "Axis", y = "Variance") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p5 <- if (nrow(N) > 0) {
    ggplot2::ggplot(N, ggplot2::aes(trait, axis, fill = weight)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2) +
      ggplot2::scale_fill_gradient2(low = "#2b6cb0", mid = "white",
                                    high = "#c53030", midpoint = 0) +
      ggplot2::labs(title = "E. Trait-function network",
                    x = "Trait", y = "Axis", fill = "Weight") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  } else NULL

  p6 <- if (nrow(T) > 0) {
    ggplot2::ggplot(T, ggplot2::aes(group1, group2, fill = turnover)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2) +
      ggplot2::scale_fill_gradient(low = "white", high = "#2b6cb0", na.value = "grey85") +
      ggplot2::labs(title = "F. Trait-function turnover",
                    x = "Group", y = "Group", fill = "Turnover") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p7 <- if (nrow(V) > 0) {
    ggplot2::ggplot(V, ggplot2::aes(axis, gamma_variance_among_clades, fill = trait)) +
      ggplot2::geom_col(position = "dodge", width = 0.7, na.rm = TRUE) +
      ggplot2::labs(title = "G. Gamma variance among clades",
                    x = "Axis", y = "Variance") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p8 <- if (nrow(A) > 0) {
    ggplot2::ggplot(A, ggplot2::aes(axis, trait, fill = trait_axis_specialization)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2) +
      ggplot2::scale_fill_gradient(low = "white", high = "#2f855a", na.value = "grey85") +
      ggplot2::labs(title = "H. Trait-axis specialization",
                    x = "Axis", y = "Trait", fill = "Share") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  .combine_plots(list(p1, p2, p3, p4, p5, p6, p7, p8),
                 title = "Trait-Function Evolution", ncol = 2)
}

.group_gamma <- function(B, traits, group, group_name) {
  if (is.null(traits) || is.null(group)) return(data.frame())
  group <- if (group_name == "clade") .align_clade(group, rownames(B)) else .align_regime(group, rownames(B))
  out <- list()
  k <- 1L
  for (g in unique(group)) {
    sp <- names(group)[group == g]
    if (length(sp) <= 2) next
    tr <- traits[sp, , drop = FALSE]
    for (axis in colnames(B)) {
      co <- .safe_lm_coefs(B[sp, axis], tr)
      if (is.null(co)) next
      co$axis <- axis
      co$gamma <- co$estimate
      co[[group_name]] <- g
      names(co)[names(co) == "term"] <- "trait"
      out[[k]] <- co[, c(group_name, "trait", "axis", "gamma"), drop = FALSE]
      k <- k + 1L
    }
  }
  if (length(out) == 0) data.frame() else do.call(rbind, out)
}

.gamma_shift_placeholder <- function(gamma) {
  if (is.null(gamma) || nrow(gamma) == 0) return(data.frame())
  unique(gamma[, c("trait", "axis")]) |>
    transform(gamma_shift_probability = NA_real_,
              source = "requires_precomputed_gamma_shift_probability")
}

.gamma_sign_flip <- function(global_gamma, clade_gamma) {
  if (nrow(global_gamma) == 0 || nrow(clade_gamma) == 0) return(data.frame())
  g <- global_gamma[, c("trait", "axis", "gamma")]
  names(g)[3] <- "global_gamma"
  d <- merge(clade_gamma, g, by = c("trait", "axis"), all.x = TRUE)
  d$flip <- sign(d$gamma) != sign(d$global_gamma)
  out <- stats::aggregate(d$flip, list(trait = d$trait, axis = d$axis), mean, na.rm = TRUE)
  names(out)[3] <- "gamma_sign_flip_index"
  out
}

.gamma_variance_among_groups <- function(clade_gamma) {
  if (nrow(clade_gamma) == 0) return(data.frame())
  out <- stats::aggregate(clade_gamma$gamma,
    list(trait = clade_gamma$trait, axis = clade_gamma$axis),
    stats::var, na.rm = TRUE)
  names(out)[3] <- "gamma_variance_among_clades"
  out
}

.gamma_turnover <- function(clade_gamma) {
  if (nrow(clade_gamma) == 0 || !"clade" %in% names(clade_gamma)) {
    return(list(index = data.frame(), matrix = matrix(numeric(0), 0, 0)))
  }
  clades <- unique(clade_gamma$clade)
  keys <- unique(paste(clade_gamma$trait, clade_gamma$axis, sep = "::"))
  M <- matrix(0, nrow = length(clades), ncol = length(keys),
              dimnames = list(clades, keys))
  for (i in seq_len(nrow(clade_gamma))) {
    M[clade_gamma$clade[i], paste(clade_gamma$trait[i], clade_gamma$axis[i], sep = "::")] <- clade_gamma$gamma[i]
  }
  D <- as.matrix(stats::dist(M))
  idx <- data.frame(group1 = rownames(D)[row(D)], group2 = colnames(D)[col(D)],
                    turnover = as.vector(D), stringsAsFactors = FALSE)
  idx <- idx[idx$group1 != idx$group2, , drop = FALSE]
  list(index = idx, matrix = D)
}

.trait_axis_specialization <- function(gamma) {
  if (nrow(gamma) == 0) return(data.frame())
  out <- lapply(split(gamma, gamma$trait), function(d) {
    total <- sum(abs(d$gamma), na.rm = TRUE)
    spec <- if (is.na(total) || total == 0) rep(NA_real_, nrow(d)) else abs(d$gamma) / total
    data.frame(trait = d$trait,
               axis = d$axis,
               trait_axis_specialization = spec,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

.trait_function_network <- function(gamma, threshold = 0) {
  if (nrow(gamma) == 0) return(data.frame())
  out <- gamma[, c("trait", "axis", "gamma", "support"), drop = FALSE]
  out$weight <- out$gamma
  out$edge <- abs(out$weight) >= threshold
  out
}
