#' Calculate trait, history, and hidden-trait mediation diagnostics
#'
#' Quantifies how measured traits and historical predictors explain Beta-derived
#' niche traits, and how much residual phylogenetic structure remains. TMNS,
#' HMNS, THMNS, and HPNS are descriptive decomposition indices; they should not
#' be interpreted as causal mediation without a causal design.
#'
#' @param evo A `hmsc_ecoevo` object.
#' @param axes Optional subset of Beta axes.
#' @param precomputed Optional list containing `trait_omission_sensitivity`,
#'   `history_explained_r2`, or other externally computed sensitivity tables.
#' @param ... Additional arguments passed to [as_hmsc_ecoevo()].
#' @return A `hmsc_trait_mediation_metrics` object.
#' @export
calc_trait_mediation_metrics <- function(evo,
                                         axes = NULL,
                                         precomputed = NULL,
                                         ...) {
  evo <- as_hmsc_ecoevo(evo, ...)
  B <- evo$beta
  if (!is.null(axes)) B <- B[, axes, drop = FALSE]
  pre <- utils::modifyList(evo$precomputed, .as_precomputed_list(precomputed))
  diagnostics <- list(messages = character(), warnings = character())

  gamma_summary <- .gamma_summary(evo$gamma, evo$gamma_draws)
  if (nrow(gamma_summary) > 0 && "axis" %in% names(gamma_summary)) {
    gamma_summary <- gamma_summary[gamma_summary$axis %in% colnames(B), , drop = FALSE]
  }
  traits <- evo$traits
  history <- evo$history_species

  trait_r2 <- .r2_table(B, traits, "trait")
  if (nrow(trait_r2) == 0) diagnostics$messages <- c(diagnostics$messages,
    "Trait-explained R2 requires species-level traits.")

  history_r2 <- if (!is.null(pre$history_explained_r2)) {
    as.data.frame(pre$history_explained_r2)
  } else .r2_table(B, history, "history")
  if (nrow(history_r2) == 0) diagnostics$messages <- c(diagnostics$messages,
    "History-explained R2 requires species-level historical predictors, such as TrData_history.")

  combined <- if (!is.null(traits) || !is.null(history)) {
    pred <- if (is.null(traits)) history else if (is.null(history)) traits else cbind(traits, history)
    .r2_table(B, pred, "trait_history")
  } else data.frame()

  residual_map <- .residual_beta_map(B, traits = traits, history = history)
  residual_rho <- .residual_phylo_signal(residual_map, evo$phylo)
  decomposition <- .mediation_decomposition(B, trait_r2, history_r2, combined, residual_rho)
  missing_trait_risk <- .missing_trait_risk(decomposition)
  trait_redundancy <- .trait_phylo_redundancy(traits, evo$phylo)
  trait_omission <- if (!is.null(pre$trait_omission_sensitivity)) {
    as.data.frame(pre$trait_omission_sensitivity)
  } else .trait_omission_table(B, traits)

  out <- list(
    gamma = gamma_summary,
    trait_explained_r2 = trait_r2,
    history_explained_r2 = history_r2,
    combined_explained_r2 = combined,
    residual_rho = residual_rho,
    signal_decomposition = decomposition,
    TMNS = decomposition[decomposition$component == "TMNS", , drop = FALSE],
    HMNS = decomposition[decomposition$component == "HMNS", , drop = FALSE],
    THMNS = decomposition[decomposition$component == "THMNS", , drop = FALSE],
    HPNS = decomposition[decomposition$component == "HPNS", , drop = FALSE],
    missing_trait_risk = missing_trait_risk,
    trait_phylogenetic_redundancy = trait_redundancy,
    trait_omission_sensitivity = trait_omission,
    residual_beta_tree_heatmap = residual_map,
    diagnostics = diagnostics,
    source = evo
  )
  class(out) <- "hmsc_trait_mediation_metrics"
  out
}

#' Plot trait/history mediation diagnostics
#'
#' Returns the fourth HmscEcoEvo diagnostic figure: Gamma matrix,
#' trait/history/hidden-signal decomposition, missing-trait risk, trait
#' phylogenetic redundancy, trait-omission sensitivity, and residual Beta map.
#'
#' @param x A `hmsc_ecoevo` or `hmsc_trait_mediation_metrics` object.
#' @param style Figure style. `"nature"` returns a single full-page
#'   Nature-style diagnostic plate; `"classic"` returns the original compact
#'   ggplot/patchwork summary.
#' @param ... Arguments passed to [calc_trait_mediation_metrics()] when needed.
#' @return A `ggplot`, `patchwork`, or printable plot-list object.
#' @export
plot_trait_mediation <- function(x, style = c("nature", "classic"), ...) {
  style <- match.arg(style)
  if (identical(style, "nature")) return(.plot_trait_mediation_nature(x, ...))
  .gg_ok()
  m <- if (inherits(x, "hmsc_trait_mediation_metrics")) x else calc_trait_mediation_metrics(x, ...)
  G <- m$gamma
  D <- m$signal_decomposition
  M <- m$missing_trait_risk
  R <- m$trait_phylogenetic_redundancy
  O <- m$trait_omission_sensitivity
  H <- m$residual_beta_tree_heatmap
  TR2 <- m$trait_explained_r2
  HR2 <- m$history_explained_r2
  CR2 <- m$combined_explained_r2

  p1 <- if (nrow(G) > 0) {
    ggplot2::ggplot(G, ggplot2::aes(axis, trait, fill = gamma)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2) +
      ggplot2::scale_fill_gradient2(low = "#2b6cb0", mid = "white",
                                    high = "#c53030", midpoint = 0) +
      ggplot2::labs(title = "A. Gamma trait-by-Beta effect matrix",
                    x = "Beta axis", y = "Trait", fill = "Gamma") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p2 <- if (nrow(D) > 0) {
    ggplot2::ggplot(D, ggplot2::aes(axis, value, fill = component)) +
      ggplot2::geom_col(position = "stack", width = 0.7, na.rm = TRUE) +
      ggplot2::labs(title = "B. Phylogenetic signal decomposition",
                    x = "Axis", y = "Descriptive share") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p3 <- if (nrow(M) > 0) {
    ggplot2::ggplot(M, ggplot2::aes(axis, missing_trait_risk_index, fill = axis)) +
      ggplot2::geom_col(width = 0.7, na.rm = TRUE) +
      ggplot2::labs(title = "C. Missing trait risk index",
                    x = "Axis", y = "Risk") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(legend.position = "none")
  } else NULL

  p4 <- if (nrow(R) > 0) {
    ggplot2::ggplot(R, ggplot2::aes(trait, redundancy, fill = moran_i)) +
      ggplot2::geom_col(width = 0.7, na.rm = TRUE) +
      ggplot2::labs(title = "D. Trait phylogenetic redundancy",
                    x = "Trait", y = "Redundancy") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  } else NULL

  p5 <- if (nrow(O) > 0) {
    ggplot2::ggplot(O, ggplot2::aes(omitted_trait, delta_r2, fill = axis)) +
      ggplot2::geom_col(position = "dodge", width = 0.7, na.rm = TRUE) +
      ggplot2::labs(title = "E. Trait omission sensitivity",
                    x = "Omitted trait", y = "Delta R2") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  } else NULL

  p6 <- if (nrow(H) > 0) {
    ggplot2::ggplot(H, ggplot2::aes(axis, species, fill = residual_beta)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2) +
      ggplot2::scale_fill_gradient2(low = "#2b6cb0", mid = "white",
                                    high = "#c53030", midpoint = 0) +
      ggplot2::labs(title = "F. Residual Beta tree heatmap",
                    x = "Axis", y = "Species", fill = "Residual") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  r2_parts <- list(
    if (nrow(TR2) > 0) data.frame(axis = TR2$axis, predictor_set = "trait", r2 = TR2$r2,
                                  stringsAsFactors = FALSE) else NULL,
    if (nrow(HR2) > 0) data.frame(axis = HR2$axis, predictor_set = "history", r2 = HR2$r2,
                                  stringsAsFactors = FALSE) else NULL,
    if (nrow(CR2) > 0) data.frame(axis = CR2$axis, predictor_set = "trait_history", r2 = CR2$r2,
                                  stringsAsFactors = FALSE) else NULL
  )
  r2_parts <- r2_parts[!vapply(r2_parts, is.null, logical(1))]
  r2_long <- if (length(r2_parts) > 0) do.call(rbind, r2_parts) else data.frame()
  p7 <- if (nrow(r2_long) > 0) {
    ggplot2::ggplot(r2_long, ggplot2::aes(axis, r2, fill = predictor_set)) +
      ggplot2::geom_col(position = "dodge", width = 0.7, na.rm = TRUE) +
      ggplot2::labs(title = "G. Trait/history explained R2",
                    x = "Axis", y = "R2") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  .combine_plots(list(p1, p2, p3, p4, p5, p6, p7),
                 title = "Trait Mediation and Hidden Traits", ncol = 2)
}

.gamma_summary <- function(gamma, gamma_draws = NULL) {
  if (!is.null(gamma_draws)) {
    out <- .posterior_summary_array(gamma_draws, value_name = "gamma")
    names(out)[names(out) == "species"] <- "trait"
    return(out)
  }
  if (is.null(gamma)) return(data.frame())
  .support_from_matrix(gamma, value_name = "gamma",
                       row_name = "trait", col_name = "axis")
}

.r2_table <- function(B, predictors, label) {
  if (is.null(predictors) || ncol(predictors) == 0) return(data.frame())
  out <- lapply(colnames(B), function(axis) {
    data.frame(axis = axis, predictor_set = label,
               r2 = .safe_lm_r2(B[, axis], predictors),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

.residual_phylo_signal <- function(residual_map, phylo) {
  if (is.null(phylo) || nrow(residual_map) == 0) return(data.frame())
  out <- lapply(split(residual_map, residual_map$axis), function(d) {
    z <- setNames(d$residual_beta, d$species)
    data.frame(axis = d$axis[1],
               residual_rho = .phylo_moran_i(z, phylo),
               metric = "residual_Moran_I",
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

.mediation_decomposition <- function(B, trait_r2, history_r2, combined, residual_rho) {
  axes <- colnames(B)
  get <- function(tab, col, axis, default = 0) {
    if (is.null(tab) || nrow(tab) == 0 || !axis %in% tab$axis) return(default)
    val <- tab[match(axis, tab$axis), col]
    ifelse(is.na(val), default, val)
  }
  out <- list()
  k <- 1L
  for (axis in axes) {
    tr <- get(trait_r2, "r2", axis)
    hr <- get(history_r2, "r2", axis)
    cr <- get(combined, "r2", axis, max(tr, hr, na.rm = TRUE))
    rr <- abs(get(residual_rho, "residual_rho", axis))
    vals <- c(TMNS = tr, HMNS = hr, THMNS = cr, HPNS = rr * (1 - cr))
    vals[!is.finite(vals)] <- NA_real_
    out[[k]] <- data.frame(axis = axis, component = names(vals), value = as.numeric(vals),
                           stringsAsFactors = FALSE)
    k <- k + 1L
  }
  do.call(rbind, out)
}

.missing_trait_risk <- function(decomposition) {
  if (nrow(decomposition) == 0) return(data.frame())
  hp <- decomposition[decomposition$component == "HPNS", c("axis", "value")]
  th <- decomposition[decomposition$component == "THMNS", c("axis", "value")]
  names(hp)[2] <- "HPNS"
  names(th)[2] <- "THMNS"
  out <- merge(hp, th, by = "axis", all = TRUE)
  out$missing_trait_risk_index <- out$HPNS * (1 - out$THMNS)
  out
}

.trait_phylo_redundancy <- function(traits, phylo) {
  if (is.null(traits) || ncol(traits) == 0) return(data.frame())
  num <- traits[vapply(traits, is.numeric, logical(1))]
  if (ncol(num) == 0) return(data.frame())
  cmat <- stats::cor(num, use = "pairwise.complete.obs")
  out <- lapply(names(num), function(tr) {
    z <- num[[tr]]
    names(z) <- rownames(num)
    moran <- if (is.null(phylo)) NA_real_ else .phylo_moran_i(z, phylo)
    redundancy <- mean(abs(cmat[tr, setdiff(colnames(cmat), tr)]), na.rm = TRUE)
    data.frame(trait = tr, moran_i = moran,
               redundancy = ifelse(is.nan(redundancy), NA_real_, redundancy),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

.trait_omission_table <- function(B, traits) {
  if (is.null(traits) || ncol(traits) == 0) return(data.frame())
  base <- .r2_table(B, traits, "trait")
  out <- list()
  k <- 1L
  for (trait in names(traits)) {
    reduced <- traits[, setdiff(names(traits), trait), drop = FALSE]
    red <- .r2_table(B, reduced, "trait_minus_one")
    for (axis in colnames(B)) {
      b <- base$r2[match(axis, base$axis)]
      r <- if (nrow(red) == 0 || !"r2" %in% names(red)) NA_real_ else red$r2[match(axis, red$axis)]
      out[[k]] <- data.frame(axis = axis, omitted_trait = trait,
                             full_r2 = b, reduced_r2 = r,
                             delta_r2 = b - r,
                             stringsAsFactors = FALSE)
      k <- k + 1L
    }
  }
  do.call(rbind, out)
}
