#' Calculate evolutionary transition and diversification diagnostics
#'
#' Computes descriptive Beta-transition metrics and standardizes optional
#' model-based results supplied by the user. HmscEcoEvo does not fit bayou,
#' SURFACE, OUwie, mvMORPH, BM/OU/EB, or adaptive-peak models internally. Pass
#' those results through `precomputed` to populate branch shift probabilities,
#' rate-shift probabilities, DTT/MDI, early-burst parameters, convergence
#' regimes, peak reuse, and distance-corrected Beta similarity.
#'
#' @param evo A `hmsc_ecoevo` object.
#' @param axes Optional subset of Beta axes.
#' @param shift_threshold Descriptive outlier threshold used to count shifts
#'   when model-based `P_shift` is unavailable.
#' @param p_shift_threshold Posterior probability threshold used to count
#'   branch niche shifts when `P_shift` is available.
#' @param integration_threshold Absolute Beta-axis correlation threshold for the
#'   evolutionary integration network.
#' @param precomputed Optional list of tables: `branch_shifts`, `rate_shifts`,
#'   `dtt`, `early_burst`, `convergence`, `peaks`, `integration_network`.
#' @param ... Additional arguments passed to [as_hmsc_ecoevo()].
#' @return A `hmsc_evo_transition_metrics` object.
#' @export
calc_evo_transition_metrics <- function(evo,
                                        axes = NULL,
                                        shift_threshold = 2,
                                        p_shift_threshold = 0.5,
                                        integration_threshold = 0.5,
                                        precomputed = NULL,
                                        ...) {
  evo <- as_hmsc_ecoevo(evo, ...)
  B <- evo$beta
  if (!is.null(axes)) B <- B[, axes, drop = FALSE]
  phylo <- evo$phylo
  pre <- utils::modifyList(evo$precomputed, .as_precomputed_list(precomputed))
  diagnostics <- list(messages = character(), warnings = character())

  branch <- .branch_transition_table(B, phylo)
  if (!is.null(pre$branch_shifts)) {
    branch <- .merge_precomputed_branch(branch, pre$branch_shifts)
  } else {
    diagnostics$messages <- c(diagnostics$messages,
      "Model-based branch shift probability P_shift requires precomputed$branch_shifts.")
    if (nrow(branch) > 0) branch$P_shift <- NA_real_
  }
  if (nrow(branch) > 0) {
    if (!"branch_contrast_outlier_score" %in% names(branch)) {
      branch$branch_contrast_outlier_score <- NA_real_
    }
    if (!"delta_theta" %in% names(branch)) branch$delta_theta <- NA_real_
    if (!"axis_specific_shift" %in% names(branch)) branch$axis_specific_shift <- branch$delta_theta
    if (!"P_shift" %in% names(branch)) branch$P_shift <- NA_real_
    branch$descriptive_shift <- !is.na(branch$branch_contrast_outlier_score) &
      branch$branch_contrast_outlier_score >= shift_threshold
    branch$model_shift <- !is.na(branch$P_shift) & branch$P_shift >= p_shift_threshold
    has_model_shift <- any(!is.na(branch$P_shift))
    branch$shift_flag <- if (has_model_shift) branch$model_shift else branch$descriptive_shift
  }

  shift_counts <- if (nrow(branch) > 0) {
    stats::aggregate(branch$shift_flag, list(axis = branch$axis), sum, na.rm = TRUE)
  } else data.frame(axis = colnames(B), x = NA_real_)
  names(shift_counts)[2] <- "number_of_niche_shifts"
  shift_support <- .shift_support_table(branch)
  multivariate_shift <- .multivariate_shift_table(branch)
  shift_count_posterior <- .shift_count_posterior_table(branch, pre$shift_count_posterior)

  rates <- .rate_table(B, phylo, evo$clade)
  if (!is.null(pre$rate_shifts)) {
    rates <- .merge_precomputed_rate(rates, pre$rate_shifts)
  } else {
    diagnostics$messages <- c(diagnostics$messages,
      "P_rate_shift requires precomputed$rate_shifts from an external rate-shift model.")
    if (nrow(rates) > 0) rates$P_rate_shift <- NA_real_
  }
  branch_rate <- .branch_rate_shift_table(phylo, pre$branch_rate_shifts)

  variance <- .variance_transition_table(B, evo$clade)
  disparity <- .disparity_table(B)
  dtt <- .optional_table(pre$dtt, name = "precomputed$dtt")
  if (is.null(dtt)) {
    dtt <- data.frame(time = NA_real_, disparity = NA_real_, axis = NA_character_,
                      source = "requires_precomputed_dtt", stringsAsFactors = FALSE)
    diagnostics$messages <- c(diagnostics$messages,
      "DTT and MDI require precomputed$dtt from an external disparity-through-time analysis.")
  }
  mdi <- if (!is.null(pre$mdi)) as.data.frame(pre$mdi) else {
    data.frame(axis = colnames(B), MDI = NA_real_,
               source = "requires_precomputed_dtt", stringsAsFactors = FALSE)
  }
  eb <- if (!is.null(pre$early_burst)) as.data.frame(pre$early_burst) else {
    diagnostics$messages <- c(diagnostics$messages,
      "Early-burst parameters require precomputed$early_burst from an external BM/OU/EB model.")
    data.frame(axis = colnames(B), early_burst_parameter = NA_real_,
               source = "requires_precomputed_early_burst", stringsAsFactors = FALSE)
  }

  pairwise <- .pairwise_beta_phylo(B, phylo)
  pairwise$beta_cosine_similarity <- .pairwise_from_matrix(
    calc_niche_metrics(evo)$cosine_similarity, pairwise$species1, pairwise$species2
  )
  distance_corrected_similarity <- .distance_corrected_similarity_table(pairwise)
  convergence <- .convergence_table(pairwise, pre$convergence)
  peaks <- .peak_reuse_table(pre$peaks, rownames(B))
  adaptive_regimes <- if (!is.null(pre$peaks)) as.data.frame(pre$peaks) else data.frame()
  integration <- .integration_network(B, threshold = integration_threshold)
  trait_niche_covariance <- .trait_niche_covariance_table(B, evo$traits)

  out <- list(
    branch_shifts = branch,
    shift_counts = shift_counts,
    shift_support = shift_support,
    shift_count_posterior = shift_count_posterior,
    multivariate_shift = multivariate_shift,
    rate_metrics = rates,
    branch_rate_shifts = branch_rate,
    variance_metrics = variance,
    disparity = disparity,
    dtt = dtt,
    mdi = mdi,
    early_burst = eb,
    convergence = convergence,
    peak_reuse = peaks,
    adaptive_regimes = adaptive_regimes,
    phylo_beta_distance = if (nrow(convergence) > 0) convergence else pairwise,
    distance_corrected_similarity = distance_corrected_similarity,
    beta_cosine_similarity = calc_niche_metrics(evo)$cosine_similarity,
    evolutionary_integration = integration$network,
    modularity = integration$modules,
    trait_niche_covariance = trait_niche_covariance,
    diagnostics = diagnostics,
    source = evo
  )
  class(out) <- "hmsc_evo_transition_metrics"
  out
}

#' Plot evolutionary transition diagnostics
#'
#' Returns the third HmscEcoEvo diagnostic figure: branch-level shift summaries,
#' shift magnitudes, rates, variance/disparity/DTT, phylogenetic distance versus
#' Beta distance, and multivariate integration/modularity.
#'
#' @param x A `hmsc_ecoevo` or `hmsc_evo_transition_metrics` object.
#' @param style Figure style. `"nature"` returns a single full-page
#'   Nature-style diagnostic plate; `"classic"` returns the original compact
#'   ggplot/patchwork summary.
#' @param ... Arguments passed to [calc_evo_transition_metrics()] when needed.
#' @return A `ggplot`, `patchwork`, or printable plot-list object.
#' @export
plot_evo_transition <- function(x, style = c("nature", "classic"), ...) {
  style <- match.arg(style)
  if (identical(style, "nature")) return(.plot_evo_transition_nature(x, ...))
  .gg_ok()
  m <- if (inherits(x, "hmsc_evo_transition_metrics")) x else calc_evo_transition_metrics(x, ...)
  B <- m$branch_shifts
  R <- m$rate_metrics
  V <- m$variance_metrics
  D <- m$dtt
  P <- m$phylo_beta_distance
  I <- m$evolutionary_integration
  SC <- m$shift_counts
  SS <- m$shift_support
  MS <- m$multivariate_shift
  EB <- m$early_burst
  MDI <- m$mdi
  COV <- m$trait_niche_covariance

  p1 <- if (nrow(B) > 0) {
    ggplot2::ggplot(B, ggplot2::aes(branch, branch_contrast_outlier_score, fill = axis)) +
      ggplot2::geom_col(position = "dodge", width = 0.7, na.rm = TRUE) +
      ggplot2::labs(title = "A. Branch-level niche shift score",
                    x = "Branch", y = "Outlier score") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank())
  } else NULL

  p2 <- if (nrow(B) > 0) {
    ggplot2::ggplot(B, ggplot2::aes(delta_theta, P_shift, color = axis)) +
      ggplot2::geom_point(na.rm = TRUE) +
      ggplot2::labs(title = "B. Shift magnitude and P_shift",
                    x = "Delta theta / branch delta", y = "P_shift") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p3 <- if (nrow(R) > 0) {
    ggplot2::ggplot(R, ggplot2::aes(axis, axis_specific_rate, fill = clade)) +
      ggplot2::geom_col(position = "dodge", width = 0.7, na.rm = TRUE) +
      ggplot2::labs(title = "C. Axis-specific evolutionary rate",
                    x = "Axis", y = "Rate") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p4 <- if (nrow(V) > 0) {
    ggplot2::ggplot(V, ggplot2::aes(clade, variance_ratio, fill = axis)) +
      ggplot2::geom_col(position = "dodge", width = 0.7, na.rm = TRUE) +
      ggplot2::geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
      ggplot2::labs(title = "D. Niche expansion/contraction",
                    x = "Clade", y = "Variance ratio") +
      ggplot2::theme_minimal(base_size = 10)
  } else if (nrow(D) > 0) {
    ggplot2::ggplot(D, ggplot2::aes(time, disparity, color = axis)) +
      ggplot2::geom_line(na.rm = TRUE) +
      ggplot2::labs(title = "D. Disparity through time",
                    x = "Time", y = "Disparity") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p5 <- if (nrow(P) > 0) {
    ggplot2::ggplot(P, ggplot2::aes(phylo_distance, beta_distance)) +
      ggplot2::geom_point(color = "#2b6cb0", alpha = 0.75, na.rm = TRUE) +
      ggplot2::geom_smooth(method = "lm", se = FALSE, color = "#c53030", na.rm = TRUE) +
      ggplot2::labs(title = "E. Phylogenetic distance vs Beta distance",
                    x = "Phylogenetic distance", y = "Beta distance") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p6 <- if (nrow(I) > 0) {
    ggplot2::ggplot(I, ggplot2::aes(axis1, axis2, fill = correlation)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2) +
      ggplot2::scale_fill_gradient2(low = "#2b6cb0", mid = "white",
                                    high = "#c53030", midpoint = 0) +
      ggplot2::labs(title = "F. Evolutionary integration network",
                    x = "Axis", y = "Axis", fill = "r") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p7 <- if (nrow(D) > 0 && all(c("time", "disparity") %in% names(D))) {
    ggplot2::ggplot(D, ggplot2::aes(time, disparity, color = axis)) +
      ggplot2::geom_line(linewidth = 0.7, na.rm = TRUE) +
      ggplot2::geom_point(size = 1.6, na.rm = TRUE) +
      ggplot2::labs(title = "G. Disparity through time",
                    x = "Relative time", y = "Disparity") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p8 <- if (nrow(SC) > 0) {
    d <- merge(SC, SS, by = "axis", all.x = TRUE)
    ggplot2::ggplot(d, ggplot2::aes(axis, number_of_niche_shifts, fill = max_P_shift)) +
      ggplot2::geom_col(width = 0.7, na.rm = TRUE) +
      ggplot2::scale_fill_gradient(low = "grey85", high = "#c53030", na.value = "grey70") +
      ggplot2::labs(title = "H. Number of shifts and max P_shift",
                    x = "Axis", y = "Number of shifts", fill = "Max P") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p9 <- if (nrow(MS) > 0) {
    ggplot2::ggplot(MS, ggplot2::aes(branch, multivariate_shift_magnitude, fill = max_P_shift)) +
      ggplot2::geom_col(width = 0.7, na.rm = TRUE) +
      ggplot2::scale_fill_gradient(low = "grey85", high = "#805ad5", na.value = "grey70") +
      ggplot2::labs(title = "I. Multivariate branch shift magnitude",
                    x = "Branch", y = "Multivariate shift", fill = "Max P") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank())
  } else NULL

  eb_mdi_parts <- list(
    if (nrow(EB) > 0 && "early_burst_parameter" %in% names(EB)) {
      data.frame(axis = EB$axis, metric = "early_burst_parameter",
                 value = EB$early_burst_parameter, stringsAsFactors = FALSE)
    } else NULL,
    if (nrow(MDI) > 0 && "MDI" %in% names(MDI)) {
      data.frame(axis = MDI$axis, metric = "MDI",
                 value = MDI$MDI, stringsAsFactors = FALSE)
    } else NULL
  )
  eb_mdi_parts <- eb_mdi_parts[!vapply(eb_mdi_parts, is.null, logical(1))]
  eb_mdi <- if (length(eb_mdi_parts) > 0) do.call(rbind, eb_mdi_parts) else data.frame()
  p10 <- if (nrow(eb_mdi) > 0) {
    ggplot2::ggplot(eb_mdi, ggplot2::aes(axis, value, fill = metric)) +
      ggplot2::geom_col(position = "dodge", width = 0.7, na.rm = TRUE) +
      ggplot2::labs(title = "J. Early burst and MDI",
                    x = "Axis", y = "Value") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p11 <- if (nrow(COV) > 0) {
    ggplot2::ggplot(COV, ggplot2::aes(axis, trait, fill = correlation)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2) +
      ggplot2::scale_fill_gradient2(low = "#2b6cb0", mid = "white",
                                    high = "#c53030", midpoint = 0,
                                    na.value = "grey85") +
      ggplot2::labs(title = "K. Trait-niche covariance",
                    x = "Axis", y = "Trait", fill = "r") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  .combine_plots(list(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11),
                 title = "Evolutionary Transition and Diversification", ncol = 2)
}

.branch_transition_table <- function(B, phylo) {
  if (is.null(phylo)) return(data.frame())
  .require_pkg("ape", "branch transition diagnostics")
  ntip <- length(phylo$tip.label)
  branch_ids <- seq_len(nrow(phylo$edge))
  out <- list()
  k <- 1L
  for (axis in colnames(B)) {
    z <- B[phylo$tip.label, axis]
    ace <- tryCatch(ape::ace(z, phylo, type = "continuous")$ace,
                    error = function(e) rep(NA_real_, phylo$Nnode))
    states <- c(setNames(z, seq_len(ntip)), setNames(ace, (ntip + 1):(ntip + phylo$Nnode)))
    delta <- states[as.character(phylo$edge[, 2])] - states[as.character(phylo$edge[, 1])]
    score <- abs(delta) / stats::sd(delta, na.rm = TRUE)
    score[!is.finite(score)] <- NA_real_
    out[[k]] <- data.frame(
      branch = branch_ids,
      parent = phylo$edge[, 1],
      child = phylo$edge[, 2],
      axis = axis,
      delta_theta = as.numeric(delta),
      axis_specific_shift = as.numeric(delta),
      branch_contrast_outlier_score = score,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
  do.call(rbind, out)
}

.merge_precomputed_branch <- function(branch, pre) {
  pre <- as.data.frame(pre)
  if ("branch_shift_probability" %in% names(pre) && !"P_shift" %in% names(pre)) {
    names(pre)[names(pre) == "branch_shift_probability"] <- "P_shift"
  }
  if (nrow(branch) == 0) return(pre)
  key <- intersect(c("branch", "axis"), names(pre))
  if (length(key) == 0) return(pre)
  out <- merge(branch, pre, by = key, all.x = TRUE, suffixes = c("", ".pre"))
  for (nm in c("P_shift", "delta_theta", "axis_specific_shift",
               "branch_contrast_outlier_score")) {
    pre_nm <- paste0(nm, ".pre")
    if (pre_nm %in% names(out)) {
      out[[nm]] <- ifelse(!is.na(out[[pre_nm]]), out[[pre_nm]], out[[nm]])
      out[[pre_nm]] <- NULL
    }
  }
  out
}

.shift_support_table <- function(branch) {
  if (nrow(branch) == 0 || !"axis" %in% names(branch)) return(data.frame())
  if (!"P_shift" %in% names(branch)) branch$P_shift <- NA_real_
  if (!"shift_flag" %in% names(branch)) branch$shift_flag <- NA
  out <- lapply(split(branch, branch$axis), function(z) {
    p <- z$P_shift[!is.na(z$P_shift)]
    data.frame(
      axis = z$axis[1],
      mean_P_shift = if (length(p) > 0) mean(p) else NA_real_,
      max_P_shift = if (length(p) > 0) max(p) else NA_real_,
      sum_P_shift = if (length(p) > 0) sum(p) else NA_real_,
      shift_count = sum(z$shift_flag, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

.shift_count_posterior_table <- function(branch, pre = NULL, n_draws = 1000) {
  if (!is.null(pre)) {
    pre <- as.data.frame(pre)
    if (!"number_of_niche_shifts" %in% names(pre)) {
      val <- names(pre)[vapply(pre, is.numeric, logical(1))][1]
      if (!is.na(val)) pre$number_of_niche_shifts <- pre[[val]]
    }
    if (!"draw" %in% names(pre)) pre$draw <- seq_len(nrow(pre))
    return(pre)
  }
  if (nrow(branch) == 0 || !"P_shift" %in% names(branch) ||
      all(is.na(branch$P_shift))) {
    return(data.frame(draw = integer(), number_of_niche_shifts = numeric(),
                      source = character()))
  }
  branch_prob <- stats::aggregate(P_shift ~ branch, branch, max, na.rm = TRUE)
  p <- pmin(pmax(branch_prob$P_shift, 0), 1)
  p <- p[is.finite(p)]
  if (length(p) == 0) {
    return(data.frame(draw = integer(), number_of_niche_shifts = numeric(),
                      source = character()))
  }
  pmf <- .poisson_binomial_pmf(p)
  u <- (seq_len(n_draws) - 0.5) / n_draws
  counts <- findInterval(u, c(0, cumsum(pmf)), rightmost.closed = TRUE)
  counts <- pmin(pmax(counts - 1L, 0L), length(pmf) - 1L)
  data.frame(
    draw = seq_len(n_draws),
    number_of_niche_shifts = counts,
    source = "deterministic_poisson_binomial_quantiles_from_branch_P_shift",
    stringsAsFactors = FALSE
  )
}

.poisson_binomial_pmf <- function(p) {
  p <- pmin(pmax(as.numeric(p), 0), 1)
  p <- p[is.finite(p)]
  if (length(p) == 0) return(1)
  prob <- c(1)
  for (pi in p) {
    prob <- c(prob * (1 - pi), 0) + c(0, prob * pi)
  }
  prob / sum(prob)
}

.branch_rate_shift_table <- function(phylo, pre = NULL) {
  if (!is.null(pre)) return(as.data.frame(pre))
  if (is.null(phylo)) {
    return(data.frame(branch = integer(), P_rate_shift = numeric(),
                      rate_ratio = numeric(), source = character()))
  }
  meta <- data.frame(branch = seq_len(nrow(phylo$edge)),
                     parent = phylo$edge[, 1],
                     child = phylo$edge[, 2],
                     branch_length = phylo$edge.length,
                     P_rate_shift = NA_real_,
                     rate_ratio = NA_real_,
                     source = "requires_precomputed_branch_rate_shifts",
                     stringsAsFactors = FALSE)
  meta
}

.multivariate_shift_table <- function(branch) {
  if (nrow(branch) == 0 || !"branch" %in% names(branch) ||
      !"delta_theta" %in% names(branch)) {
    return(data.frame())
  }
  if (!"P_shift" %in% names(branch)) branch$P_shift <- NA_real_
  if (!"shift_flag" %in% names(branch)) branch$shift_flag <- NA
  out <- lapply(split(branch, branch$branch), function(z) {
    p <- z$P_shift[!is.na(z$P_shift)]
    data.frame(
      branch = z$branch[1],
      parent = if ("parent" %in% names(z)) z$parent[1] else NA,
      child = if ("child" %in% names(z)) z$child[1] else NA,
      multivariate_shift_magnitude = sqrt(sum(z$delta_theta^2, na.rm = TRUE)),
      max_P_shift = if (length(p) > 0) max(p) else NA_real_,
      n_axis_shifts = sum(z$shift_flag, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

.trait_niche_covariance_table <- function(B, traits) {
  if (is.null(traits)) {
    return(data.frame(trait = character(), axis = character(),
                      covariance = numeric(), correlation = numeric(),
                      n_species = integer(), source = character()))
  }
  traits <- .align_species_frame(traits, rownames(B), "traits")
  keep <- vapply(traits, is.numeric, logical(1))
  traits <- traits[, keep, drop = FALSE]
  if (ncol(traits) == 0) return(data.frame())
  out <- list()
  k <- 1L
  for (tr in colnames(traits)) {
    for (axis in colnames(B)) {
      ok <- !is.na(traits[[tr]]) & !is.na(B[, axis])
      out[[k]] <- data.frame(
        trait = tr,
        axis = axis,
        covariance = if (sum(ok) > 1) stats::cov(traits[[tr]][ok], B[ok, axis]) else NA_real_,
        correlation = if (sum(ok) > 2) stats::cor(traits[[tr]][ok], B[ok, axis]) else NA_real_,
        n_species = sum(ok),
        source = "descriptive_species_trait_beta_covariance",
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  do.call(rbind, out)
}

.rate_table <- function(B, phylo, clade = NULL) {
  out <- list()
  k <- 1L
  for (axis in colnames(B)) {
    rate <- if (!is.null(phylo)) {
      z <- B[phylo$tip.label, axis]
      pic <- tryCatch(ape::pic(z, phylo), error = function(e) NA_real_)
      mean(pic^2, na.rm = TRUE)
    } else stats::var(B[, axis], na.rm = TRUE)
    out[[k]] <- data.frame(axis = axis, clade = "global",
      axis_specific_rate = rate, rate_ratio = 1, stringsAsFactors = FALSE)
    k <- k + 1L
  }
  if (!is.null(clade)) {
    clade <- .align_clade(clade, rownames(B))
    for (cl in unique(clade)) {
      sp <- names(clade)[clade == cl]
      if (length(sp) < 2) next
      for (axis in colnames(B)) {
        cr <- stats::var(B[sp, axis], na.rm = TRUE)
        gr <- stats::var(B[, axis], na.rm = TRUE)
        out[[k]] <- data.frame(axis = axis, clade = cl,
          axis_specific_rate = cr,
          rate_ratio = ifelse(is.na(gr) || gr == 0, NA_real_, cr / gr),
          stringsAsFactors = FALSE)
        k <- k + 1L
      }
    }
  }
  do.call(rbind, out)
}

.merge_precomputed_rate <- function(rates, pre) {
  pre <- as.data.frame(pre)
  if (nrow(rates) == 0) return(pre)
  key <- intersect(c("axis", "clade"), names(pre))
  if (length(key) == 0) return(pre)
  out <- merge(rates, pre, by = key, all.x = TRUE, suffixes = c("", ".pre"))
  for (nm in c("axis_specific_rate", "rate_ratio", "P_rate_shift")) {
    pre_nm <- paste0(nm, ".pre")
    if (pre_nm %in% names(out)) {
      out[[nm]] <- ifelse(!is.na(out[[pre_nm]]), out[[pre_nm]], out[[nm]])
      out[[pre_nm]] <- NULL
    }
  }
  out
}

.variance_transition_table <- function(B, clade) {
  if (is.null(clade)) return(data.frame())
  clade <- .align_clade(clade, rownames(B))
  out <- list()
  k <- 1L
  for (cl in unique(clade)) {
    sp <- names(clade)[clade == cl]
    if (length(sp) < 2) next
    for (axis in colnames(B)) {
      cv <- stats::var(B[sp, axis], na.rm = TRUE)
      gv <- stats::var(B[, axis], na.rm = TRUE)
      ratio <- ifelse(is.na(gv) || gv == 0, NA_real_, cv / gv)
      out[[k]] <- data.frame(
        clade = cl,
        axis = axis,
        clade_niche_variance = cv,
        global_niche_variance = gv,
        variance_ratio = ratio,
        niche_expansion_index = pmax(ratio - 1, 0),
        niche_contraction_index = pmax(1 - ratio, 0),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  if (length(out) == 0) data.frame() else do.call(rbind, out)
}

.disparity_table <- function(B) {
  data.frame(
    metric = c("total_multivariate_disparity", paste0("axis_disparity:", colnames(B))),
    value = c(sum(diag(stats::var(B, na.rm = TRUE)), na.rm = TRUE),
              apply(B, 2, stats::var, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}

.pairwise_from_matrix <- function(M, sp1, sp2) {
  vapply(seq_along(sp1), function(i) M[sp1[i], sp2[i]], numeric(1))
}

.convergence_table <- function(pairwise, pre = NULL) {
  if (!is.null(pre)) return(as.data.frame(pre))
  if (nrow(pairwise) == 0) return(data.frame())
  pd <- pairwise$phylo_distance
  bd <- pairwise$beta_distance
  high_phylo <- pd >= stats::quantile(pd, 0.75, na.rm = TRUE)
  low_beta <- bd <= stats::quantile(bd, 0.25, na.rm = TRUE)
  pairwise$convergence_score <- as.numeric(high_phylo & low_beta)
  pairwise$source <- "descriptive_far_phylo_low_beta_distance"
  pairwise
}

.distance_corrected_similarity_table <- function(pairwise) {
  empty <- data.frame(
    species1 = character(),
    species2 = character(),
    phylo_distance = numeric(),
    beta_cosine_similarity = numeric(),
    expected_beta_cosine_similarity = numeric(),
    distance_corrected_similarity = numeric(),
    source = character(),
    stringsAsFactors = FALSE
  )
  if (nrow(pairwise) == 0) return(empty)
  out <- pairwise[, intersect(c("species1", "species2", "phylo_distance",
                                "beta_distance", "beta_cosine_similarity"),
                              names(pairwise)), drop = FALSE]
  if (!"beta_cosine_similarity" %in% names(out)) out$beta_cosine_similarity <- NA_real_
  if (!"phylo_distance" %in% names(out)) out$phylo_distance <- NA_real_
  out$expected_beta_cosine_similarity <- NA_real_
  out$distance_corrected_similarity <- NA_real_
  ok <- is.finite(out$phylo_distance) & is.finite(out$beta_cosine_similarity)
  enough <- sum(ok) >= 3 && length(unique(out$phylo_distance[ok])) >= 2
  if (!enough) {
    out$source <- "insufficient_phylogenetic_distance_gradient"
    return(out)
  }
  fit <- tryCatch(stats::lm(beta_cosine_similarity ~ phylo_distance,
                            data = out[ok, , drop = FALSE]),
                  error = function(e) NULL)
  if (is.null(fit)) {
    out$source <- "distance_correction_failed"
    return(out)
  }
  out$expected_beta_cosine_similarity[ok] <-
    stats::predict(fit, newdata = out[ok, , drop = FALSE])
  out$distance_corrected_similarity <-
    out$beta_cosine_similarity - out$expected_beta_cosine_similarity
  out$source <- "descriptive_lm_residual_cosine_on_phylogenetic_distance_not_causal"
  out
}

.peak_reuse_table <- function(peaks, species) {
  if (is.null(peaks)) {
    return(data.frame(peak_reuse_index = NA_real_,
                      source = "requires_precomputed_peaks",
                      stringsAsFactors = FALSE))
  }
  peaks <- as.data.frame(peaks)
  .require_cols(peaks, c("species", "peak"), "precomputed$peaks")
  peaks <- peaks[peaks$species %in% species, , drop = FALSE]
  tab <- table(peaks$peak)
  denom <- choose(sum(tab), 2)
  reuse <- ifelse(denom == 0, NA_real_, sum(choose(tab, 2)) / denom)
  data.frame(peak_reuse_index = as.numeric(reuse),
             n_peaks = length(tab),
             source = "precomputed_peaks", stringsAsFactors = FALSE)
}

.integration_network <- function(B, threshold = 0.5) {
  C <- stats::cor(B, use = "pairwise.complete.obs")
  diag(C) <- NA_real_
  net <- data.frame(
    axis1 = rownames(C)[row(C)],
    axis2 = colnames(C)[col(C)],
    correlation = as.vector(C),
    stringsAsFactors = FALSE
  )
  net <- net[net$axis1 != net$axis2, , drop = FALSE]
  net$edge <- abs(net$correlation) >= threshold
  modules <- data.frame(axis = colnames(B), module = NA_integer_, stringsAsFactors = FALSE)
  if (ncol(B) >= 2) {
    d <- stats::as.dist(1 - abs(C))
    hc <- tryCatch(stats::hclust(d), error = function(e) NULL)
    if (!is.null(hc)) modules$module <- stats::cutree(hc, k = min(2, ncol(B)))
  }
  list(network = net, modules = modules)
}
