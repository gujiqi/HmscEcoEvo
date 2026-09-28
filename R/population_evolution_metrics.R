#' Calculate population-level and dynamic eco-evolutionary diagnostics
#'
#' Computes optional population-level and time-series diagnostics: population
#' Beta, within-species niche divergence, local adaptation and plasticity
#' contributions, LA:PL ratio, genetic niche signal, reaction norms, G x E,
#' eco-evolutionary feedback paths, population trait shifts, trait-evolution
#' contribution, and eco-evo turnover ratio. Feedback-path summaries are
#' descriptive unless the user supplies a causal or mechanistic model output.
#'
#' @param evo A `hmsc_ecoevo` object.
#' @param precomputed Optional list of population or dynamic diagnostic tables.
#' @param ... Additional arguments passed to [as_hmsc_ecoevo()].
#' @return A `hmsc_population_evolution_metrics` object.
#' @export
calc_population_evolution_metrics <- function(evo,
                                              precomputed = NULL,
                                              ...) {
  evo <- as_hmsc_ecoevo(evo, ...)
  pre <- utils::modifyList(evo$precomputed, .as_precomputed_list(precomputed))
  diagnostics <- list(messages = character(), warnings = character())
  pop <- evo$population
  ts <- evo$timeseries

  if (is.null(pop)) diagnostics$messages <- c(diagnostics$messages,
    "Population-level Beta diagnostics require population data with species, population, axis, and beta.")

  population_beta <- pop %||% data.frame(
    species = character(), population = character(), axis = character(), beta = numeric()
  )
  divergence <- .within_species_divergence(population_beta)
  partition <- .local_plasticity_partition(population_beta, pre)
  genetic <- .genetic_niche_signal(population_beta, pre)
  reaction <- .reaction_norms(population_beta, pre)
  gxe <- .gxe_table(population_beta, pre)
  trait_shift <- .population_trait_shift(population_beta, pre)
  dynamic <- .dynamic_feedback(ts, pre)

  out <- list(
    population_beta = population_beta,
    within_species_niche_divergence = divergence,
    local_adaptation_contribution = partition[, intersect(c("species", "axis", "local_adaptation_contribution"), names(partition)), drop = FALSE],
    plasticity_contribution = partition[, intersect(c("species", "axis", "plasticity_contribution"), names(partition)), drop = FALSE],
    LA_PL_ratio = partition[, intersect(c("species", "axis", "LA_PL_ratio"), names(partition)), drop = FALSE],
    variance_partition = partition,
    genetic_niche_signal = genetic,
    reaction_norms = reaction,
    GxE = gxe,
    population_trait_shift = trait_shift,
    eco_evolutionary_feedback_path = dynamic$feedback_path,
    beta_t_shift_rate = dynamic$beta_t_shift_rate,
    trait_evolution_contribution = dynamic$trait_evolution_contribution,
    community_driven_selection_index = dynamic$community_driven_selection_index,
    feedback_strength = dynamic$feedback_strength,
    eco_evo_turnover_ratio = dynamic$eco_evo_turnover_ratio,
    lagged_niche_shift = dynamic$lagged_niche_shift,
    interaction_mediated_evolution = dynamic$interaction_mediated_evolution,
    diagnostics = diagnostics,
    source = evo
  )
  class(out) <- "hmsc_population_evolution_metrics"
  out
}

#' Plot population and dynamic eco-evolution diagnostics
#'
#' Returns the sixth HmscEcoEvo diagnostic figure: species-level versus
#' population-level Beta, within-species divergence, local-adaptation/plasticity
#' partitioning, LA:PL ratio, genetic niche signal/reaction norms, and feedback
#' diagnostics.
#'
#' @param x A `hmsc_ecoevo` or `hmsc_population_evolution_metrics` object.
#' @param style Figure style. `"nature"` returns a single full-page
#'   Nature-style diagnostic plate; `"classic"` returns the original compact
#'   ggplot/patchwork summary.
#' @param ... Arguments passed to [calc_population_evolution_metrics()] when
#'   needed.
#' @return A `ggplot`, `patchwork`, or printable plot-list object.
#' @export
plot_population_evolution <- function(x, style = c("nature", "classic"), ...) {
  style <- match.arg(style)
  if (identical(style, "nature")) return(.plot_population_evolution_nature(x, ...))
  .gg_ok()
  m <- if (inherits(x, "hmsc_population_evolution_metrics")) x else calc_population_evolution_metrics(x, ...)
  P <- m$population_beta
  D <- m$within_species_niche_divergence
  V <- m$variance_partition
  G <- m$genetic_niche_signal
  R <- m$reaction_norms
  F <- m$feedback_strength
  X <- m$GxE
  T <- m$trait_evolution_contribution
  ETR <- m$eco_evo_turnover_ratio
  PTS <- m$population_trait_shift

  p1 <- if (nrow(P) > 0) {
    ggplot2::ggplot(P, ggplot2::aes(axis, beta, color = population)) +
      ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.08), na.rm = TRUE) +
      ggplot2::facet_wrap(~species) +
      ggplot2::labs(title = "A. Population-level Beta",
                    x = "Axis", y = "Population Beta") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p2 <- if (nrow(D) > 0) {
    ggplot2::ggplot(D, ggplot2::aes(stats::reorder(species, within_species_niche_divergence),
                                    within_species_niche_divergence, fill = axis)) +
      ggplot2::geom_col(position = "dodge", width = 0.7, na.rm = TRUE) +
      ggplot2::coord_flip() +
      ggplot2::labs(title = "B. Within-species niche divergence",
                    x = "Species", y = "Mean population Beta distance") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p3 <- if (nrow(V) > 0) {
    vv <- V[, intersect(c("species", "axis", "local_adaptation_contribution", "plasticity_contribution"), names(V)), drop = FALSE]
    if (ncol(vv) >= 4) {
      vv <- stats::reshape(vv, varying = c("local_adaptation_contribution", "plasticity_contribution"),
                           v.names = "value", timevar = "component",
                           times = c("local_adaptation", "plasticity"),
                           direction = "long")
      ggplot2::ggplot(vv, ggplot2::aes(axis, value, fill = component)) +
        ggplot2::geom_col(position = "stack", width = 0.7, na.rm = TRUE) +
        ggplot2::facet_wrap(~species) +
        ggplot2::labs(title = "C. Niche-difference variance partition",
                      x = "Axis", y = "Contribution") +
        ggplot2::theme_minimal(base_size = 10)
    } else NULL
  } else NULL

  p4 <- if (nrow(V) > 0 && "LA_PL_ratio" %in% names(V)) {
    ggplot2::ggplot(V, ggplot2::aes(axis, LA_PL_ratio, fill = species)) +
      ggplot2::geom_col(position = "dodge", width = 0.7, na.rm = TRUE) +
      ggplot2::geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
      ggplot2::labs(title = "D. Local adaptation to plasticity ratio",
                    x = "Axis", y = "LA:PL") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p5 <- if (nrow(R) > 0) {
    ggplot2::ggplot(R, ggplot2::aes(environment, response, color = genotype, group = genotype)) +
      ggplot2::geom_point(na.rm = TRUE) +
      ggplot2::geom_smooth(method = "lm", se = FALSE, na.rm = TRUE) +
      ggplot2::facet_wrap(~species) +
      ggplot2::labs(title = "E. Reaction norms / GxE",
                    x = "Environment", y = "Response") +
      ggplot2::theme_minimal(base_size = 10)
  } else if (nrow(G) > 0) {
    ggplot2::ggplot(G, ggplot2::aes(axis, genetic_niche_signal, fill = species)) +
      ggplot2::geom_col(position = "dodge", width = 0.7, na.rm = TRUE) +
      ggplot2::labs(title = "E. Genetic niche signal",
                    x = "Axis", y = "Signal") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p6 <- if (nrow(F) > 0) {
    ggplot2::ggplot(F, ggplot2::aes(path, feedback_strength, fill = path)) +
      ggplot2::geom_col(width = 0.7, na.rm = TRUE) +
      ggplot2::labs(title = "F. Eco-evolutionary feedback diagnostics",
                    x = "Path", y = "Feedback strength") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     legend.position = "none")
  } else NULL

  p7 <- if (nrow(G) > 0) {
    ggplot2::ggplot(G, ggplot2::aes(axis, genetic_niche_signal, fill = species)) +
      ggplot2::geom_col(position = "dodge", width = 0.7, na.rm = TRUE) +
      ggplot2::labs(title = "G. Genetic niche signal",
                    x = "Axis", y = "Signal") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  p8 <- if (nrow(X) > 0 && "GxE" %in% names(X)) {
    ggplot2::ggplot(X, ggplot2::aes(stats::reorder(species, GxE), GxE, fill = GxE)) +
      ggplot2::geom_col(width = 0.7, na.rm = TRUE) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_gradient2(low = "#2b6cb0", mid = "white",
                                    high = "#c53030", midpoint = 0,
                                    na.value = "grey85") +
      ggplot2::labs(title = "H. GxE interaction coefficient",
                    x = "Species", y = "GxE") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(legend.position = "none")
  } else NULL

  p9 <- if (nrow(T) > 0 && all(c("trait", "contribution") %in% names(T))) {
    ggplot2::ggplot(T, ggplot2::aes(trait, contribution, fill = trait)) +
      ggplot2::geom_col(width = 0.7, na.rm = TRUE) +
      ggplot2::labs(title = "I. Trait evolution contribution",
                    x = "Trait", y = "Contribution") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(legend.position = "none")
  } else NULL

  p10 <- if (nrow(ETR) > 0) {
    d <- ETR
    if (!"metric" %in% names(d)) d$metric <- "eco_evo_turnover_ratio"
    if (!"value" %in% names(d)) d$value <- d[[names(d)[vapply(d, is.numeric, logical(1))][1]]]
    ggplot2::ggplot(d, ggplot2::aes(metric, value, fill = metric)) +
      ggplot2::geom_col(width = 0.7, na.rm = TRUE) +
      ggplot2::labs(title = "J. Eco-evolutionary turnover ratio",
                    x = "Metric", y = "Value") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(legend.position = "none")
  } else NULL

  p11 <- if (nrow(PTS) > 0 && "population_trait_shift" %in% names(PTS)) {
    ggplot2::ggplot(PTS, ggplot2::aes(trait, population_trait_shift, fill = species)) +
      ggplot2::geom_col(position = "dodge", width = 0.7, na.rm = TRUE) +
      ggplot2::labs(title = "K. Population trait shift",
                    x = "Trait", y = "Mean trait shift") +
      ggplot2::theme_minimal(base_size = 10)
  } else NULL

  .combine_plots(list(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11),
                 title = "Population-Level and Dynamic Eco-Evolution", ncol = 2)
}

.within_species_divergence <- function(pop) {
  if (is.null(pop) || nrow(pop) == 0) return(data.frame())
  out <- list()
  k <- 1L
  for (sp in unique(pop$species)) {
    d <- pop[pop$species == sp, , drop = FALSE]
    for (axis in unique(d$axis)) {
      z <- d$beta[d$axis == axis]
      div <- if (length(z) < 2) NA_real_ else mean(abs(stats::dist(matrix(z, ncol = 1))), na.rm = TRUE)
      out[[k]] <- data.frame(species = sp, axis = axis,
                             within_species_niche_divergence = div,
                             stringsAsFactors = FALSE)
      k <- k + 1L
    }
  }
  do.call(rbind, out)
}

.local_plasticity_partition <- function(pop, pre) {
  if (!is.null(pre$local_plasticity_partition)) return(as.data.frame(pre$local_plasticity_partition))
  if (is.null(pop) || nrow(pop) == 0) {
    return(data.frame(species = character(), axis = character(),
      local_adaptation_contribution = numeric(), plasticity_contribution = numeric(),
      LA_PL_ratio = numeric()))
  }
  if (!all(c("local_adaptation", "plasticity") %in% names(pop))) {
    return(data.frame(species = unique(pop$species), axis = NA_character_,
      local_adaptation_contribution = NA_real_, plasticity_contribution = NA_real_,
      LA_PL_ratio = NA_real_, source = "requires_columns_or_precomputed",
      stringsAsFactors = FALSE))
  }
  out <- stats::aggregate(cbind(local_adaptation, plasticity) ~ species + axis,
                          data = pop, mean, na.rm = TRUE)
  out$local_adaptation_contribution <- out$local_adaptation
  out$plasticity_contribution <- out$plasticity
  out$LA_PL_ratio <- out$local_adaptation / out$plasticity
  out$LA_PL_ratio[!is.finite(out$LA_PL_ratio)] <- NA_real_
  out
}

.genetic_niche_signal <- function(pop, pre) {
  if (!is.null(pre$genetic_niche_signal)) return(as.data.frame(pre$genetic_niche_signal))
  if (is.null(pop) || nrow(pop) == 0 || !all(c("genetic_value", "beta") %in% names(pop))) return(data.frame())
  out <- list()
  k <- 1L
  for (sp in unique(pop$species)) {
    d1 <- pop[pop$species == sp, , drop = FALSE]
    for (axis in unique(d1$axis)) {
      d <- d1[d1$axis == axis, , drop = FALSE]
      out[[k]] <- data.frame(species = sp, axis = axis,
        genetic_niche_signal = stats::cor(d$genetic_value, d$beta, use = "pairwise.complete.obs"),
        stringsAsFactors = FALSE)
      k <- k + 1L
    }
  }
  do.call(rbind, out)
}

.reaction_norms <- function(pop, pre) {
  if (!is.null(pre$reaction_norms)) return(as.data.frame(pre$reaction_norms))
  if (is.null(pop) || nrow(pop) == 0 ||
      !all(c("environment", "response", "genotype") %in% names(pop))) return(data.frame())
  pop[, c("species", "population", "genotype", "environment", "response"), drop = FALSE]
}

.gxe_table <- function(pop, pre) {
  if (!is.null(pre$GxE)) return(as.data.frame(pre$GxE))
  if (is.null(pop) || nrow(pop) == 0 ||
      !all(c("environment", "response", "genotype") %in% names(pop))) return(data.frame())
  out <- list()
  k <- 1L
  for (sp in unique(pop$species)) {
    d <- pop[pop$species == sp, , drop = FALSE]
    fit <- tryCatch(stats::lm(response ~ environment * genotype, data = d), error = function(e) NULL)
    if (is.null(fit)) next
    co <- suppressWarnings(stats::coef(summary(fit)))
    idx <- grepl("environment:genotype|genotype.*:environment", rownames(co))
    out[[k]] <- data.frame(species = sp,
      GxE = ifelse(any(idx), mean(co[idx, "Estimate"], na.rm = TRUE), NA_real_),
      stringsAsFactors = FALSE)
    k <- k + 1L
  }
  if (length(out) == 0) data.frame() else do.call(rbind, out)
}

.population_trait_shift <- function(pop, pre) {
  if (!is.null(pre$population_trait_shift)) return(as.data.frame(pre$population_trait_shift))
  empty <- data.frame(species = character(), trait = character(),
                      population_trait_shift = numeric(),
                      pairwise_trait_divergence = numeric(),
                      n_population = integer(), source = character(),
                      stringsAsFactors = FALSE)
  if (is.null(pop) || nrow(pop) == 0) return(empty)
  if ("trait_shift" %in% names(pop)) {
    d <- pop
    if (!"trait" %in% names(d)) d$trait <- "trait_shift"
    out <- stats::aggregate(trait_shift ~ species + trait, data = d, mean, na.rm = TRUE)
    names(out)[names(out) == "trait_shift"] <- "population_trait_shift"
    n <- stats::aggregate(population ~ species + trait, data = d, function(z) length(unique(z)))
    names(n)[names(n) == "population"] <- "n_population"
    out <- merge(out, n, by = c("species", "trait"), all.x = TRUE)
    out$pairwise_trait_divergence <- NA_real_
    out$source <- "population_trait_shift_column"
    return(out[, c("species", "trait", "population_trait_shift",
                   "pairwise_trait_divergence", "n_population", "source")])
  }
  if (!"trait_value" %in% names(pop)) return(empty)
  d <- pop
  if (!"trait" %in% names(d)) d$trait <- "trait_value"
  d <- unique(d[, intersect(c("species", "population", "trait", "trait_value",
                              "reference_trait_value", "species_trait_mean"),
                            names(d)), drop = FALSE])
  out <- lapply(split(d, list(d$species, d$trait), drop = TRUE), function(z) {
    if (nrow(z) == 0) return(NULL)
    ref <- NULL
    if ("reference_trait_value" %in% names(z) && any(!is.na(z$reference_trait_value))) {
      ref <- z$reference_trait_value
    } else if ("species_trait_mean" %in% names(z) && any(!is.na(z$species_trait_mean))) {
      ref <- z$species_trait_mean
    }
    shift <- if (!is.null(ref)) {
      mean(abs(z$trait_value - ref), na.rm = TRUE)
    } else {
      mean(abs(z$trait_value - mean(z$trait_value, na.rm = TRUE)), na.rm = TRUE)
    }
    div <- if (sum(!is.na(z$trait_value)) >= 2) {
      mean(abs(stats::dist(matrix(z$trait_value, ncol = 1))), na.rm = TRUE)
    } else NA_real_
    data.frame(
      species = z$species[1],
      trait = z$trait[1],
      population_trait_shift = shift,
      pairwise_trait_divergence = div,
      n_population = length(unique(z$population)),
      source = if (!is.null(ref)) "population_trait_value_from_reference" else
        "population_trait_value_deviation_from_species_mean",
      stringsAsFactors = FALSE
    )
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0) empty else do.call(rbind, out)
}

.dynamic_feedback <- function(ts, pre) {
  get_pre <- function(name, cols = character()) {
    if (!is.null(pre[[name]])) return(as.data.frame(pre[[name]]))
    data.frame()
  }
  feedback_path <- get_pre("eco_evolutionary_feedback_path")
  beta_rate <- get_pre("beta_t_shift_rate")
  trait_contrib <- get_pre("trait_evolution_contribution")
  selection <- get_pre("community_driven_selection_index")
  feedback <- get_pre("feedback_strength")
  turnover <- get_pre("eco_evo_turnover_ratio")
  lagged <- get_pre("lagged_niche_shift")
  interaction <- get_pre("interaction_mediated_evolution")

  if (!is.null(ts) && nrow(ts) > 1) {
    if (all(c("species", "axis", "time", "beta") %in% names(ts)) && nrow(beta_rate) == 0) {
      beta_rate <- do.call(rbind, lapply(split(ts, list(ts$species, ts$axis), drop = TRUE), function(d) {
        d <- d[order(d$time), , drop = FALSE]
        data.frame(species = d$species[1], axis = d$axis[1],
                   beta_t_shift_rate = mean(abs(diff(d$beta)) / diff(d$time), na.rm = TRUE),
                   stringsAsFactors = FALSE)
      }))
    }
    if (all(c("ecological_turnover", "evolutionary_turnover") %in% names(ts)) && nrow(turnover) == 0) {
      turnover <- data.frame(eco_evo_turnover_ratio =
        mean(ts$evolutionary_turnover, na.rm = TRUE) / mean(ts$ecological_turnover, na.rm = TRUE),
        source = "descriptive_timeseries", stringsAsFactors = FALSE)
    }
    if (all(c("community_metric", "beta") %in% names(ts)) && nrow(feedback) == 0) {
      feedback <- data.frame(path = "community_to_beta",
        feedback_strength = stats::cor(ts$community_metric, ts$beta, use = "pairwise.complete.obs"),
        source = "descriptive_correlation_not_causal", stringsAsFactors = FALSE)
    }
  }

  list(
    feedback_path = feedback_path,
    beta_t_shift_rate = beta_rate,
    trait_evolution_contribution = trait_contrib,
    community_driven_selection_index = selection,
    feedback_strength = feedback,
    eco_evo_turnover_ratio = turnover,
    lagged_niche_shift = lagged,
    interaction_mediated_evolution = interaction
  )
}
