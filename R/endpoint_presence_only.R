#' Audit whether a community matrix supplies confirmed absences
#'
#' A matrix assembled from occurrence records commonly has one row for every
#' cell in which at least one focal taxon was observed.  Zeros for the other
#' taxa in those rows are then non-detections, not confirmed absences.  This
#' helper makes that distinction explicit before an endpoint score is used to
#' calibrate a dynamic-distribution model.
#'
#' @param observed Numeric or logical cell-by-taxon matrix. Values must be
#'   finite non-negative counts or binary records.
#' @param survey_complete Optional logical vector of length `nrow(observed)`.
#'   `TRUE` means the row came from a survey in which all focal taxa could have
#'   been detected. When omitted, no confirmed-absence claim is made.
#' @param design_label Optional user-supplied description of the observation
#'   design.
#'
#' @return A one-row data frame describing zero semantics and whether a
#'   Bernoulli presence--absence endpoint likelihood is eligible.
#' @export
#'
#' @examples
#' y <- rbind(c(1, 0), c(0, 1))
#' hee_audit_observation_design(y)
hee_audit_observation_design <- function(observed,
                                         survey_complete = NULL,
                                         design_label = NULL) {
  if (is.data.frame(observed)) observed <- as.matrix(observed)
  if (!is.matrix(observed) || !is.numeric(observed) || !nrow(observed) ||
      !ncol(observed) || any(!is.finite(observed)) || any(observed < 0)) {
    stop("observed must be a non-empty finite, non-negative numeric matrix.",
         call. = FALSE)
  }
  if (is.null(survey_complete)) {
    survey_complete <- rep(FALSE, nrow(observed))
    survey_complete_supplied <- FALSE
  } else {
    if (!is.logical(survey_complete) || length(survey_complete) != nrow(observed) ||
        anyNA(survey_complete)) {
      stop("survey_complete must be a non-missing logical vector with one value per row.",
           call. = FALSE)
    }
    survey_complete_supplied <- TRUE
  }
  row_records <- rowSums(observed > 0)
  all_rows_have_focal_record <- all(row_records > 0)
  all_rows_survey_complete <- all(survey_complete)
  absence_eligible <- all_rows_survey_complete
  inferred_design <- if (absence_eligible) {
    "surveyed_presence_absence"
  } else if (all_rows_have_focal_record) {
    "target_occurrence_grid_non_detections_not_absences"
  } else {
    "observation_design_unknown_non_detections_not_absences"
  }
  data.frame(
    n_cells = nrow(observed),
    n_taxa = ncol(observed),
    n_zero_record_cells = sum(row_records == 0),
    all_rows_have_focal_record = all_rows_have_focal_record,
    survey_complete_supplied = survey_complete_supplied,
    all_rows_survey_complete = all_rows_survey_complete,
    bernoulli_absence_likelihood_eligible = absence_eligible,
    observation_design = inferred_design,
    design_label = if (is.null(design_label) || !length(design_label)) {
      NA_character_
    } else as.character(design_label)[1L],
    scientific_boundary = if (absence_eligible) {
      "Confirmed surveyed absences may be used in a Bernoulli endpoint likelihood, subject to detection modelling."
    } else {
      "Zeros must not be interpreted as confirmed absences; use a presence-only or explicitly detection-aware endpoint analysis."
    },
    stringsAsFactors = FALSE
  )
}

#' Score predicted spatial support using presence-only endpoint records
#'
#' This computes the conditional point-pattern score of known occurrences
#' under a non-negative predicted spatial intensity.  For taxon \eqn{j}, with
#' observed cell counts \eqn{n_{ij}} and optional effort/exposure \eqn{e_i},
#' the normalized sampling intensity is
#' \deqn{\pi_{ij}=\frac{e_i p_{ij}}{\sum_u e_u p_{uj}}}
#' and the score is \eqn{\sum_i n_{ij}\log(\pi_{ij})}.  It uses no record as
#' a confirmed absence.  With `effort = NULL`, the reference is a uniform
#' available-cell sampling surface; this is a diagnostic, not a correction for
#' GBIF sampling bias.
#'
#' @param probability Numeric cell-by-taxon matrix of finite values in
#'   `[0, 1]`. Rows must be all candidate modern cells, not only occupied
#'   occurrence cells.
#' @param occurrence Non-negative numeric cell-by-taxon matrix of occurrence
#'   counts or binary records, aligned to `probability`.
#' @param effort Optional finite non-negative cell vector giving a declared
#'   observation-effort or exposure surface. It is normalized internally and
#'   must be positive in cells containing records.
#' @param data_role One of `"independent_validation"`, `"spatial_holdout"`, or
#'   `"training_diagnostic"`. The role is recorded; it does not turn reused
#'   HMSC training records into independent evidence.
#' @param occurrence_representation Either `"binary_grid"` (the default) or
#'   `"record_count"`. `"binary_grid"` counts at most one known occurrence of
#'   each taxon per grid cell and is the defensible default for aggregated GBIF
#'   data. `"record_count"` gives repeated records within a cell repeated
#'   likelihood weight; it should be used only with an observation model that
#'   justifies that sampling-unit interpretation.
#' @param epsilon Positive lower bound used only to avoid `log(0)` in a
#'   declared numerical likelihood.
#'
#' @return A list containing a one-row `summary` and a per-taxon `taxa` table.
#' @export
#'
#' @examples
#' p <- rbind(c(0.8, 0.2), c(0.2, 0.8), c(0.1, 0.1))
#' y <- rbind(c(1, 0), c(0, 1), c(0, 0))
#' hee_endpoint_presence_only_score(
#'   p, y, data_role = "training_diagnostic"
#' )$summary
hee_endpoint_presence_only_score <- function(
    probability,
    occurrence,
    effort = NULL,
    data_role = c("independent_validation", "spatial_holdout", "training_diagnostic"),
    occurrence_representation = c("binary_grid", "record_count"),
    epsilon = 1e-12) {
  data_role <- match.arg(data_role)
  occurrence_representation <- match.arg(occurrence_representation)
  if (is.data.frame(probability)) probability <- as.matrix(probability)
  if (is.data.frame(occurrence)) occurrence <- as.matrix(occurrence)
  if (!is.matrix(probability) || !is.numeric(probability) ||
      !is.matrix(occurrence) || !is.numeric(occurrence) ||
      !identical(dim(probability), dim(occurrence))) {
    stop("probability and occurrence must be numeric matrices with identical dimensions.",
         call. = FALSE)
  }
  if (any(!is.finite(probability)) || any(probability < 0 | probability > 1) ||
      any(!is.finite(occurrence)) || any(occurrence < 0)) {
    stop("probability must be in [0, 1] and occurrence must be finite and non-negative.",
         call. = FALSE)
  }
  if (!is.null(rownames(probability)) && !is.null(rownames(occurrence)) &&
      !identical(rownames(probability), rownames(occurrence))) {
    stop("probability and occurrence cell rows are not aligned.", call. = FALSE)
  }
  if (!is.null(colnames(probability)) && !is.null(colnames(occurrence)) &&
      !identical(colnames(probability), colnames(occurrence))) {
    stop("probability and occurrence taxon columns are not aligned.", call. = FALSE)
  }
  if (length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0 ||
      epsilon >= 1) {
    stop("epsilon must be one finite number in (0, 1).", call. = FALSE)
  }
  if (is.null(effort)) {
    if (data_role != "training_diagnostic") {
      stop(
        "Independent or spatial-holdout presence-only endpoint scoring requires ",
        "a declared effort/exposure surface. A uniform available-cell surface is ",
        "only permitted for training_diagnostic.",
        call. = FALSE
      )
    }
    effort <- rep(1, nrow(probability))
    effort_source <- "uniform_available_cell_reference"
  } else {
    if (!is.numeric(effort) || length(effort) != nrow(probability) ||
        any(!is.finite(effort)) || any(effort < 0)) {
      stop("effort must be a finite non-negative numeric vector with one value per cell.",
           call. = FALSE)
    }
    effort_source <- "user_supplied_effort_surface"
  }
  occurrence_raw <- occurrence
  occurrence <- if (occurrence_representation == "binary_grid") {
    (occurrence > 0) * 1
  } else {
    occurrence
  }
  positive_cells <- rowSums(occurrence) > 0
  if (any(effort <= 0 & positive_cells)) {
    stop("effort must be positive in every cell that contains an occurrence record.",
         call. = FALSE)
  }
  taxon <- if (is.null(colnames(probability))) {
    paste0("taxon_", seq_len(ncol(probability)))
  } else colnames(probability)
  n_records_raw <- colSums(occurrence_raw)
  n_records <- colSums(occurrence)
  n_positive_cells <- colSums(occurrence > 0)
  support_mass <- colSums(probability * effort)
  availability <- sum(effort > 0)
  if (!is.finite(availability) || availability < 1L) {
    stop("effort must define at least one available cell.", call. = FALSE)
  }
  conditional_log_score <- mean_log_support <- uniform_log_score <-
    relative_log_support <- rep(NA_real_, ncol(probability))
  for (j in seq_len(ncol(probability))) {
    if (n_records[[j]] <= 0 || !is.finite(support_mass[[j]]) ||
        support_mass[[j]] <= 0) next
    pi_j <- pmax(effort * probability[, j] / support_mass[[j]], epsilon)
    conditional_log_score[[j]] <- sum(occurrence[, j] * log(pi_j))
    mean_log_support[[j]] <- conditional_log_score[[j]] / n_records[[j]]
    uniform_log_score[[j]] <- -log(availability)
    relative_log_support[[j]] <- mean_log_support[[j]] - uniform_log_score[[j]]
  }
  taxa <- data.frame(
    taxon = taxon,
    occurrence_representation = occurrence_representation,
    n_occurrence_records_raw = as.numeric(n_records_raw),
    n_occurrence_units = as.numeric(n_records),
    n_positive_grid_cells = as.numeric(n_positive_cells),
    predicted_effort_weighted_mass = as.numeric(support_mass),
    conditional_log_score = conditional_log_score,
    mean_log_support_per_record = mean_log_support,
    uniform_mean_log_score = uniform_log_score,
    mean_log_support_lift_over_uniform = relative_log_support,
    stringsAsFactors = FALSE
  )
  usable <- is.finite(taxa$mean_log_support_per_record) &
    taxa$n_occurrence_units > 0
  summary <- data.frame(
    data_role = data_role,
    endpoint_type = "presence_only_conditional_spatial_support",
    effort_source = effort_source,
    occurrence_representation = occurrence_representation,
    n_cells = nrow(probability),
    n_taxa = ncol(probability),
    n_taxa_with_records_and_support = sum(usable),
    n_occurrence_records_raw = sum(n_records_raw),
    n_occurrence_units = sum(n_records),
    n_positive_grid_cells = sum(n_positive_cells),
    mean_log_support_per_record = if (any(usable)) {
      stats::weighted.mean(taxa$mean_log_support_per_record[usable],
                           taxa$n_occurrence_units[usable])
    } else NA_real_,
    mean_log_support_lift_over_uniform = if (any(usable)) {
      stats::weighted.mean(taxa$mean_log_support_lift_over_uniform[usable],
                           taxa$n_occurrence_units[usable])
    } else NA_real_,
    scientific_boundary = if (data_role == "training_diagnostic") {
      "Training-record presence-only diagnostic; it is not an independent endpoint likelihood or posterior particle weight. Binary-grid scoring avoids duplicate within-cell records but does not remove spatial sampling bias."
    } else {
      "Presence-only conditional spatial-support score with a declared effort surface; absence, prevalence, detection, and residual spatial dependence are not estimated."
    },
    stringsAsFactors = FALSE
  )
  list(summary = summary, taxa = taxa)
}

#' Reweight complete dynamic-history particles with presence-only endpoint data
#'
#' This is the particle analogue of [hee_endpoint_presence_only_score()].  It
#' scores every complete modern occupancy map against occurrence locations as a
#' conditional spatial point pattern and converts the resulting scores to
#' normalized particle weights.  Crucially, no unrecorded cell is treated as a
#' confirmed absence.  It is therefore suitable for an independently sampled
#' occurrence endpoint or for a spatially cross-fitted endpoint, but it does
#' **not** make GBIF records that already informed HMSC into independent
#' evidence.
#'
#' Let `d` index a complete particle and let ψᵈᵢⱼ be its predicted
#' terminal occupancy.  With occurrence counts nᵢⱼ and declared effort
#' eᵢ, the score used here is
#' \deqn{\log L_d=\sum_{j,i}n_{ij}\log\left(
#' \frac{e_i\psi^{(d)}_{ij}}{\sum_u e_u\psi^{(d)}_{uj}}\right).}
#' The normalized weight is proportional to `exp(likelihood_power * log L_d)`.
#' This conditional score identifies spatial allocation among available cells;
#' it does not estimate prevalence, detection, or absolute abundance.
#'
#' @param particles A named list of finite cell-by-taxon matrices with
#'   identical dimensions and aligned dimnames. By default each must be a
#'   binary terminal occupancy state from one dynamic-history particle.
#' @param occurrence Non-negative cell-by-taxon occurrence counts or binary
#'   records, aligned to each particle.
#' @param effort Optional finite non-negative cell-level observation effort or
#'   exposure surface. It must be positive wherever a record occurs.
#' @param data_role One of `"independent_validation"`, `"spatial_holdout"`, or
#'   `"training_diagnostic"`. Only the first two may be described as an
#'   endpoint-conditioned ensemble.
#' @param occurrence_representation Passed to
#'   [hee_endpoint_presence_only_score()]. The default aggregates repeated
#'   records to a unique taxon-by-grid occurrence before weighting particles.
#' @param particle_representation Either `"binary_state"` (the default) or
#'   `"probability_map"`. Only binary states can be described as dynamic
#'   posterior particles. A continuous `"probability_map"` may be scored as a
#'   sensitivity/scenario map but is not a particle posterior.
#' @param likelihood_power Non-negative exponent on the conditional endpoint
#'   log score. A value below one is a disclosed sensitivity device, not a cure
#'   for reused observations.
#' @param epsilon Positive numerical lower bound passed to
#'   [hee_endpoint_presence_only_score()].
#'
#' @return A data frame with per-particle scores, normalized weights, ensemble
#'   effective sample size, and an explicit scientific-role label.
#' @export
#'
#' @examples
#' p <- list(
#'   focused = rbind(c(1, 0), c(0, 1), c(0, 0)),
#'   diffuse = rbind(c(1, 1), c(1, 1), c(1, 1))
#' )
#' y <- rbind(c(1, 0), c(0, 1), c(0, 0))
#' hee_endpoint_presence_only_particle_weights(
#'   p, y, data_role = "training_diagnostic"
#' )
hee_endpoint_presence_only_particle_weights <- function(
    particles,
    occurrence,
    effort = NULL,
    data_role = c("independent_validation", "spatial_holdout", "training_diagnostic"),
    occurrence_representation = c("binary_grid", "record_count"),
    particle_representation = c("binary_state", "probability_map"),
    likelihood_power = 1,
    epsilon = 1e-12) {
  data_role <- match.arg(data_role)
  occurrence_representation <- match.arg(occurrence_representation)
  particle_representation <- match.arg(particle_representation)
  if (!is.list(particles) || !length(particles) || is.null(names(particles)) ||
      any(!nzchar(names(particles))) || anyDuplicated(names(particles))) {
    stop("particles must be a non-empty named list with unique particle IDs.",
         call. = FALSE)
  }
  if (is.data.frame(occurrence)) occurrence <- as.matrix(occurrence)
  if (!is.matrix(occurrence) || !is.numeric(occurrence) ||
      any(!is.finite(occurrence)) || any(occurrence < 0)) {
    stop("occurrence must be a finite non-negative numeric matrix.", call. = FALSE)
  }
  likelihood_power <- suppressWarnings(as.numeric(likelihood_power)[1L])
  if (!is.finite(likelihood_power) || likelihood_power < 0) {
    stop("likelihood_power must be one finite non-negative number.", call. = FALSE)
  }

  score_one <- function(probability, particle_id) {
    if (is.data.frame(probability)) probability <- as.matrix(probability)
    if (!is.matrix(probability) || !is.numeric(probability) ||
        !identical(dim(probability), dim(occurrence))) {
      stop("Particle '", particle_id,
           "' must be a numeric matrix aligned to occurrence.", call. = FALSE)
    }
    if (any(!is.finite(probability)) || any(probability < 0 | probability > 1)) {
      stop("Particle '", particle_id,
           "' must contain finite values in [0, 1].", call. = FALSE)
    }
    if (particle_representation == "binary_state" &&
        any(!(probability %in% c(0, 1)))) {
      stop("Particle '", particle_id,
           "' is not a binary terminal occupancy state. Use particle_representation='probability_map' only for a non-posterior scenario-map sensitivity.",
           call. = FALSE)
    }
    if (!is.null(rownames(probability)) && !is.null(rownames(occurrence)) &&
        !identical(rownames(probability), rownames(occurrence))) {
      stop("Particle '", particle_id, "' cell rows are not aligned to occurrence.",
           call. = FALSE)
    }
    if (!is.null(colnames(probability)) && !is.null(colnames(occurrence)) &&
        !identical(colnames(probability), colnames(occurrence))) {
      stop("Particle '", particle_id, "' taxon columns are not aligned to occurrence.",
           call. = FALSE)
    }
    hee_endpoint_presence_only_score(
      probability = probability, occurrence = occurrence, effort = effort,
      data_role = data_role,
      occurrence_representation = occurrence_representation,
      epsilon = epsilon
    )
  }

  scores <- lapply(names(particles), function(id) score_one(particles[[id]], id))
  summaries <- do.call(rbind, lapply(scores, `[[`, "summary"))
  log_score <- summaries$mean_log_support_per_record * summaries$n_occurrence_units
  if (any(!is.finite(log_score))) {
    stop("At least one particle has no finite conditional endpoint score.",
         call. = FALSE)
  }
  log_weight <- likelihood_power * log_score
  centered <- log_weight - max(log_weight)
  weight <- exp(centered)
  weight <- weight / sum(weight)
  data.frame(
    particle_id = names(particles),
    conditional_log_endpoint_score = unname(log_score),
    mean_log_support_per_record = summaries$mean_log_support_per_record,
    mean_log_support_lift_over_uniform = summaries$mean_log_support_lift_over_uniform,
    occurrence_representation = occurrence_representation,
    particle_representation = particle_representation,
    likelihood_power = likelihood_power,
    weight = unname(weight),
    squared_weight = unname(weight^2),
    ensemble_effective_sample_size = 1 / sum(weight^2),
    data_role = data_role,
    scientific_role = if (data_role == "training_diagnostic") {
      "diagnostic_only_not_an_independent_posterior_particle_weight"
    } else if (particle_representation == "binary_state") {
      "presence_only_endpoint_conditioned_binary_particle_weight"
    } else {
      "scenario_probability_map_reweighting_not_particle_posterior"
    },
    scientific_boundary = if (data_role == "training_diagnostic") {
      "The occurrence data are reused diagnostic information; these weights must not be called a posterior reconstruction."
    } else {
      if (particle_representation == "binary_state") {
        "Weights condition binary terminal particle states by spatial allocation among available cells only; prevalence, absences, imperfect detection, and residual spatial dependence are not estimated."
      } else {
        "Continuous maps are reweighted only as a scenario sensitivity. They are not binary dynamic-history particles and must not be called a posterior palaeodistribution ensemble."
      }
    },
    stringsAsFactors = FALSE
  )
}
