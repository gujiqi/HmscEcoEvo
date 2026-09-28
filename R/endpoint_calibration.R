#' Mean-field occupancy-composition summaries with a valid zero-occupancy state
#'
#' Computes entropy and concentration summaries from a non-negative
#' cell-by-lineage occupancy matrix. These are **mean-field occupancy
#' composition summaries**, not posterior expectations of observed abundance
#' diversity indices. They summarize the relative distribution of expected
#' lineage occupancy mass within one cell. In general,
#' `diversity(E[Z]) != E[diversity(Z)]`; posterior Shannon, Simpson, PD, and
#' functional diversity therefore require binary or abundance community draws
#' before each draw-level index is calculated and summarized. Rows whose total
#' occupancy mass is at or below `zero_tolerance` are explicitly assigned zero
#' diversity; this avoids the invalid result in which an empty cell is
#' normalized into a perfectly even community.
#'
#' @param occupancy A numeric matrix with cells in rows and lineages or species
#'   in columns. Values must be finite and non-negative.
#' @param zero_tolerance Non-negative total occupancy mass below which a row is
#'   treated as unoccupied.
#'
#' @return A data frame with total occupancy mass, clearly labelled mean-field
#'   composition summaries, backwards-compatible `occupancy_weighted_*`
#'   aliases, and a `diversity_defined` flag. Row names are retained as
#'   `cell_id` when present.
#' @export
#'
#' @examples
#' x <- rbind(c(0.5, 0.5), c(0, 0))
#' hee_occupancy_weighted_diversity(x)
hee_occupancy_weighted_diversity <- function(occupancy,
                                             zero_tolerance = 1e-10) {
  if (is.data.frame(occupancy)) occupancy <- as.matrix(occupancy)
  if (!is.matrix(occupancy) || !is.numeric(occupancy)) {
    stop("occupancy must be a numeric matrix or data frame.", call. = FALSE)
  }
  if (length(zero_tolerance) != 1L || !is.finite(zero_tolerance) ||
      zero_tolerance < 0) {
    stop("zero_tolerance must be one finite non-negative number.", call. = FALSE)
  }
  if (any(!is.finite(occupancy))) {
    stop("occupancy must not contain NA, NaN or infinite values.", call. = FALSE)
  }
  if (any(occupancy < 0)) {
    stop("occupancy must be non-negative.", call. = FALSE)
  }

  occupancy_mass <- rowSums(occupancy)
  defined <- occupancy_mass > zero_tolerance
  relative <- matrix(0, nrow = nrow(occupancy), ncol = ncol(occupancy))
  if (any(defined)) {
    relative[defined, ] <- occupancy[defined, , drop = FALSE] /
      occupancy_mass[defined]
  }
  entropy <- -rowSums(ifelse(relative > 0, relative * log(relative), 0))
  concentration <- rowSums(relative * relative)
  gini_simpson <- 1 - concentration
  effective_shannon <- exp(entropy)
  effective_simpson <- ifelse(concentration > 0, 1 / concentration, 0)

  entropy[!defined] <- 0
  gini_simpson[!defined] <- 0
  effective_shannon[!defined] <- 0
  effective_simpson[!defined] <- 0
  out <- data.frame(
    occupancy_mass = occupancy_mass,
    mean_field_compositional_shannon_entropy = entropy,
    mean_field_compositional_gini_simpson = gini_simpson,
    mean_field_effective_lineages_shannon = effective_shannon,
    mean_field_effective_lineages_simpson = effective_simpson,
    # Retained so legacy scripts remain runnable. New reports must use the
    # mean_field_* names and must not call them posterior diversity.
    occupancy_weighted_shannon_entropy = entropy,
    occupancy_weighted_gini_simpson = gini_simpson,
    occupancy_weighted_effective_lineages_shannon = effective_shannon,
    occupancy_weighted_effective_lineages_simpson = effective_simpson,
    diversity_defined = defined,
    scientific_boundary = ifelse(
      defined,
      "mean_field_occupancy_composition_not_posterior_draw_diversity",
      "empty_or_negligible_mean_field_occupancy"
    ),
    stringsAsFactors = FALSE
  )
  if (!is.null(rownames(occupancy))) out$cell_id <- rownames(occupancy)
  out
}

#' Score a modern endpoint prediction against supplied validation observations
#'
#' This helper evaluates a probability matrix against a modern endpoint data
#' set. The caller must state whether observations are independent, spatially
#' held out, or only a training diagnostic. It does not turn unobserved records
#' into detection-corrected absences and does not make a historical claim.
#'
#' @param probability Numeric cell-by-species probability matrix in `[0, 1]`.
#' @param observed Numeric or logical cell-by-species matrix of 0/1 endpoint
#'   observations with identical dimensions and, where supplied, identical
#'   dimnames.
#' @param data_role One of `"independent_validation"`, `"spatial_holdout"`, or
#'   `"training_diagnostic"`.
#' @param epsilon Small clipping constant used only for log scores and logits.
#'
#' @return A list containing one-row `summary`, per-species `species`, and
#'   per-cell `cells` score tables.
#' @export
#'
#' @examples
#' p <- rbind(c(0.8, 0.2), c(0.1, 0.6))
#' y <- rbind(c(1, 0), c(0, 1))
#' hee_modern_endpoint_score(p, y, data_role = "spatial_holdout")$summary
hee_modern_endpoint_score <- function(
    probability,
    observed,
    data_role = c("independent_validation", "spatial_holdout", "training_diagnostic"),
    epsilon = 1e-8) {
  data_role <- match.arg(data_role)
  if (is.data.frame(probability)) probability <- as.matrix(probability)
  if (is.data.frame(observed)) observed <- as.matrix(observed)
  if (!is.matrix(probability) || !is.numeric(probability) ||
      !is.matrix(observed) || !is.numeric(observed)) {
    stop("probability and observed must be numeric matrices or data frames.",
         call. = FALSE)
  }
  if (!identical(dim(probability), dim(observed))) {
    stop("probability and observed must have identical dimensions.", call. = FALSE)
  }
  if (!is.null(colnames(probability)) && !is.null(colnames(observed)) &&
      !identical(colnames(probability), colnames(observed))) {
    stop("probability and observed species columns must be identically ordered.",
         call. = FALSE)
  }
  if (!is.null(rownames(probability)) && !is.null(rownames(observed)) &&
      !identical(rownames(probability), rownames(observed))) {
    stop("probability and observed cell rows must be identically ordered.",
         call. = FALSE)
  }
  if (any(!is.finite(probability)) || any(probability < 0 | probability > 1)) {
    stop("probability must contain finite values in [0, 1].", call. = FALSE)
  }
  if (any(!is.finite(observed)) || any(!(observed %in% c(0, 1)))) {
    stop("observed must contain finite binary 0/1 values.", call. = FALSE)
  }
  if (length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0 ||
      epsilon >= 0.5) {
    stop("epsilon must be one finite number in (0, 0.5).", call. = FALSE)
  }

  clipped <- pmin(pmax(probability, epsilon), 1 - epsilon)
  log_score <- observed * log(clipped) + (1 - observed) * log(1 - clipped)
  brier <- (observed - probability)^2
  species <- data.frame(
    species = if (is.null(colnames(probability))) {
      paste0("species_", seq_len(ncol(probability)))
    } else colnames(probability),
    n_cells = nrow(probability),
    observed_prevalence = colMeans(observed),
    predicted_prevalence = colMeans(probability),
    brier_score = colMeans(brier),
    mean_log_score = colMeans(log_score),
    stringsAsFactors = FALSE
  )
  cells <- data.frame(
    cell_id = if (is.null(rownames(probability))) {
      paste0("cell_", seq_len(nrow(probability)))
    } else rownames(probability),
    observed_richness = rowSums(observed),
    predicted_expected_richness = rowSums(probability),
    cell_brier_score = rowMeans(brier),
    cell_mean_log_score = rowMeans(log_score),
    stringsAsFactors = FALSE
  )
  pred_vector <- as.vector(probability)
  obs_vector <- as.vector(observed)
  fit <- tryCatch(
    stats::glm(obs_vector ~ stats::qlogis(pmin(pmax(pred_vector, epsilon), 1 - epsilon)),
               family = stats::binomial()),
    error = function(e) NULL
  )
  co <- if (is.null(fit)) c(NA_real_, NA_real_) else stats::coef(fit)
  richness_cor <- if (stats::sd(cells$observed_richness) == 0 ||
                      stats::sd(cells$predicted_expected_richness) == 0) {
    NA_real_
  } else {
    stats::cor(cells$observed_richness, cells$predicted_expected_richness)
  }
  summary <- data.frame(
    data_role = data_role,
    n_cells = nrow(probability),
    n_species = ncol(probability),
    observed_prevalence = mean(obs_vector),
    predicted_prevalence = mean(pred_vector),
    brier_score = mean(brier),
    mean_log_score = mean(log_score),
    calibration_intercept = unname(co[[1L]]),
    calibration_slope = unname(co[[2L]]),
    observed_mean_richness = mean(cells$observed_richness),
    predicted_mean_richness = mean(cells$predicted_expected_richness),
    richness_correlation = richness_cor,
    stringsAsFactors = FALSE
  )
  list(summary = summary, species = species, cells = cells)
}

#' Reweight complete historical particles by a modern endpoint observation
#'
#' Calculates an explicitly declared endpoint likelihood for a list of complete
#' cell-by-species dynamic-occupancy particles. This is useful when uncertain
#' root locations, plate-transport alternatives, or process scenarios have
#' been forward-simulated and an **independent** modern endpoint is available.
#' It does not edit a particle's historical states, turn the endpoint into an
#' ancestral niche, or make non-detections detection-corrected absences.
#'
#' With `data_role = "training_diagnostic"`, the returned weights are a
#' diagnostic, not independent posterior weights: the same modern observations
#' may already have informed the HMSC component. A value of `likelihood_power`
#' below one is a transparent power-likelihood sensitivity device, not a remedy
#' for duplicated data use.
#'
#' @param particles A named list of numeric cell-by-species probability
#'   matrices, all with identical dimensions and aligned dimnames.
#' @param observed Numeric or logical binary cell-by-species endpoint matrix.
#' @param data_role One of `"independent_validation"`, `"spatial_holdout"`, or
#'   `"training_diagnostic"`.
#' @param likelihood_power Non-negative exponent applied to each particle's
#'   Bernoulli log score before normalization. Use `1` for a properly specified
#'   independent likelihood; choose and report another value only as a
#'   sensitivity analysis.
#' @param epsilon Small clipping constant used only for log scores.
#'
#' @return A data frame with particle log score, normalized and squared weights,
#'   ensemble effective sample size, and the scientific role of the endpoint.
#' @export
#'
#' @examples
#' p <- list(a = matrix(c(0.8, 0.2), ncol = 1),
#'           b = matrix(c(0.3, 0.7), ncol = 1))
#' y <- matrix(c(1, 0), ncol = 1)
#' hee_endpoint_particle_weights(p, y, data_role = "spatial_holdout")
hee_endpoint_particle_weights <- function(
    particles,
    observed,
    data_role = c("independent_validation", "spatial_holdout", "training_diagnostic"),
    likelihood_power = 1,
    epsilon = 1e-8) {
  data_role <- match.arg(data_role)
  if (!is.list(particles) || !length(particles)) {
    stop("particles must be a non-empty named list of probability matrices.",
         call. = FALSE)
  }
  if (is.null(names(particles)) || any(!nzchar(names(particles))) ||
      anyDuplicated(names(particles))) {
    stop("particles must have unique, non-empty names.", call. = FALSE)
  }
  if (is.data.frame(observed)) observed <- as.matrix(observed)
  if (!is.matrix(observed) || !is.numeric(observed) ||
      any(!is.finite(observed)) || any(!(observed %in% c(0, 1)))) {
    stop("observed must be a finite binary numeric matrix.", call. = FALSE)
  }
  if (length(likelihood_power) != 1L || !is.finite(likelihood_power) ||
      likelihood_power < 0) {
    stop("likelihood_power must be one finite non-negative number.",
         call. = FALSE)
  }
  if (length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0 ||
      epsilon >= 0.5) {
    stop("epsilon must be one finite number in (0, 0.5).", call. = FALSE)
  }

  score_one <- function(probability, id) {
    if (is.data.frame(probability)) probability <- as.matrix(probability)
    if (!is.matrix(probability) || !is.numeric(probability) ||
        !identical(dim(probability), dim(observed)) ||
        any(!is.finite(probability)) || any(probability < 0 | probability > 1)) {
      stop("Particle '", id,
           "' must be a finite [0, 1] matrix with observed dimensions.",
           call. = FALSE)
    }
    if (!is.null(rownames(probability)) && !is.null(rownames(observed)) &&
        !identical(rownames(probability), rownames(observed))) {
      stop("Particle '", id, "' cell rows are not aligned to observed.",
           call. = FALSE)
    }
    if (!is.null(colnames(probability)) && !is.null(colnames(observed)) &&
        !identical(colnames(probability), colnames(observed))) {
      stop("Particle '", id, "' species columns are not aligned to observed.",
           call. = FALSE)
    }
    clipped <- pmin(pmax(probability, epsilon), 1 - epsilon)
    sum(observed * log(clipped) + (1 - observed) * log1p(-clipped))
  }
  log_score <- vapply(names(particles), function(id) score_one(particles[[id]], id),
                      numeric(1))
  log_weight <- likelihood_power * log_score
  centered <- log_weight - max(log_weight)
  weight <- exp(centered)
  weight <- weight / sum(weight)
  out <- data.frame(
    particle_id = names(particles),
    log_endpoint_score = unname(log_score),
    likelihood_power = likelihood_power,
    log_weight_unnormalized = unname(log_weight),
    weight = unname(weight),
    squared_weight = unname(weight^2),
    data_role = data_role,
    scientific_role = if (data_role == "training_diagnostic") {
      "diagnostic_reweighting_not_independent_posterior"
    } else "endpoint_conditioned_particle_weight",
    stringsAsFactors = FALSE
  )
  out$ensemble_effective_sample_size <- 1 / sum(out$weight^2)
  out
}
