#' Update global-grid dynamic occupancy for one interval
#'
#' This is the matrix implementation of a *calibrated or explicitly
#' scenario-labelled* no-BioGeoBEARS grid-level update. It keeps process roles
#' explicit: `movement_matrix` supplies P2 dispersal/arrival, `eta` supplies
#' P1 environmental filtering, optional biotic terms supply P3 scenario or
#' externally measured effects, and the CTMC transition combines colonisation
#' and local persistence. It does not apply geographic accessibility
#' multipliers, regional history gates, or extrapolation weights.
#'
#' It is not the default empirical phylogeographic reconstruction for a
#' one-time modern occurrence matrix. Such data do not identify historical
#' establishment and persistence logits. When those parameters have no
#' independent dynamic calibration, use [hee_dispersal_location_transition()],
#' [hee_dispersal_time_ordered_pruning()], and separately report HMSC-derived
#' ancestral environmental support rather than interpreting a CTMC scenario as
#' historical occupancy.
#'
#' @param occupancy Numeric cell-by-lineage occupancy probabilities, or binary
#'   occupancy states when `stochastic = TRUE`.
#' @param movement_matrix Square target-by-source movement-rate matrix. Sparse
#'   matrices are supported.
#' @param eta Cell-by-lineage environmental linear predictor.
#' @param delta_t Positive interval duration in Myr.
#' @param reporting_interval_myr Positive duration in Myr used only to
#'   translate instantaneous arrival and colonisation hazards into explicitly
#'   labelled reporting probabilities. It defaults to `delta_t`. It has no
#'   effect on the CTMC state transition. Supplying a common value (for
#'   example, one Myr) is useful when maps from several processes are compared.
#' @param persistence_reference_interval Positive duration in Myr represented
#'   by `persistence_intercept` and `persistence_slope`. It defaults to
#'   `delta_t` for backward compatibility. When an environmental interval is
#'   integrated with several numerical substeps, keep this value fixed at the
#'   biological reference interval and vary only `delta_t`; otherwise reducing
#'   a numerical step would inadvertently change the local-loss hazard.
#' @param establishment_intercept,establishment_slope Parameters of the
#'   conditional establishment logit. Each may be a scalar, a named cell
#'   vector, a named lineage vector, or a cell-by-lineage matrix. Named vectors
#'   are aligned to the occupancy dimnames, preventing accidental recycling of
#'   lineage-specific process baselines across cells.
#' @param persistence_intercept,persistence_slope Parameters of the
#'   persistence logit in the same scalar, vector, or matrix forms.
#' @param movement_multiplier Optional non-negative scalar, named cell or
#'   lineage vector, or cell-by-lineage matrix applied to occupied **source**
#'   cells before the movement matrix. It supports an explicitly supplied
#'   lineage-specific dispersal scenario without changing the target-side
#'   movement graph or confusing the multiplier with environmental suitability.
#' @param colonisation_biotic,persistence_biotic Optional scalar or
#'   cell-by-lineage additive logit effects. They must originate from an
#'   external interaction model or be labelled as scenarios.
#' @param habitat_state Optional cell-level or cell-by-lineage hard habitat
#'   state. Zeros force occupancy to zero after the transition.
#' @param stochastic Draw binary next states when `TRUE`; otherwise return the
#'   deterministic CTMC expectation.
#' @param seed Optional random seed used only for `stochastic = TRUE`.
#' @return A list with `occupancy_next`, arrival and transition quantities.
#' @export
#'
#' @examples
#' q <- matrix(c(1, 0), ncol = 1)
#' m <- matrix(c(0, 0.2, 0, 0), nrow = 2)
#' hee_projection_global_grid_step(q, m, matrix(0, 2, 1), delta_t = 1)
hee_projection_global_grid_step <- function(occupancy,
                                            movement_matrix,
                                            eta,
                                            delta_t,
                                            reporting_interval_myr = delta_t,
                                            persistence_reference_interval = delta_t,
                                            establishment_intercept = -0.25,
                                            establishment_slope = 0.85,
                                            persistence_intercept = 1.25,
                                            persistence_slope = 0.75,
                                            movement_multiplier = 1,
                                            colonisation_biotic = 0,
                                            persistence_biotic = 0,
                                            habitat_state = NULL,
                                            stochastic = FALSE,
                                            seed = NULL) {
  q <- as.matrix(occupancy)
  eta <- as.matrix(eta)
  storage.mode(q) <- "double"
  storage.mode(eta) <- "double"
  if (!identical(dim(q), dim(eta))) {
    stop("occupancy and eta must have identical cell-by-lineage dimensions.",
         call. = FALSE)
  }
  if (any(!is.finite(q)) || any(q < 0 | q > 1)) {
    stop("occupancy must be finite and within [0, 1].", call. = FALSE)
  }
  delta_t <- suppressWarnings(as.numeric(delta_t)[1L])
  if (!is.finite(delta_t) || delta_t <= 0) {
    stop("delta_t must be a positive finite duration in Myr.", call. = FALSE)
  }
  reporting_interval_myr <- suppressWarnings(
    as.numeric(reporting_interval_myr)[1L]
  )
  if (!is.finite(reporting_interval_myr) || reporting_interval_myr <= 0) {
    stop("reporting_interval_myr must be a positive finite duration in Myr.",
         call. = FALSE)
  }
  # A cached Matrix::dgCMatrix can be restored before Matrix's S4 methods are
  # registered in a fresh Rscript session. Register them before dim() and %*%
  # are used, so sparse production movement caches behave like in-memory ones.
  if (methods::is(movement_matrix, "Matrix")) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Matrix is required to use a sparse movement_matrix.", call. = FALSE)
    }
  }
  persistence_reference_interval <- suppressWarnings(
    as.numeric(persistence_reference_interval)[1L]
  )
  if (!is.finite(persistence_reference_interval) ||
      persistence_reference_interval <= 0) {
    stop("persistence_reference_interval must be a positive finite duration in Myr.",
         call. = FALSE)
  }
  md <- dim(movement_matrix)
  if (length(md) != 2L || md[[1L]] != nrow(q) || md[[2L]] != nrow(q)) {
    stop("movement_matrix must be square with one row and column per cell.",
         call. = FALSE)
  }
  as_effect_matrix <- function(x, name) {
    if (length(x) == 1L) return(matrix(as.numeric(x), nrow(q), ncol(q)))
    # On a square grid, a lineage vector can have the same length as the
    # number of cells. Names disambiguate the intended process dimension.
    if (!is.null(names(x)) && !is.null(colnames(q))) {
      column_hit <- match(colnames(q), names(x))
      if (!anyNA(column_hit)) {
        value <- as.numeric(x)[column_hit]
        z <- matrix(rep(value, each = nrow(q)), nrow(q), ncol(q),
                    dimnames = dimnames(q))
        z[!is.finite(z)] <- 0
        return(z)
      }
    }
    if (!is.null(names(x)) && !is.null(rownames(q))) {
      row_hit <- match(rownames(q), names(x))
      if (!anyNA(row_hit)) {
        value <- as.numeric(x)[row_hit]
        z <- matrix(rep(value, ncol(q)), nrow(q), ncol(q),
                    dimnames = dimnames(q))
        z[!is.finite(z)] <- 0
        return(z)
      }
    }
    if (length(x) == nrow(q)) {
      value <- as.numeric(x)
      z <- matrix(rep(value, ncol(q)), nrow(q), ncol(q),
                  dimnames = dimnames(q))
      z[!is.finite(z)] <- 0
      return(z)
    }
    if (length(x) == ncol(q)) {
      value <- as.numeric(x)
      z <- matrix(rep(value, each = nrow(q)), nrow(q), ncol(q),
                  dimnames = dimnames(q))
      z[!is.finite(z)] <- 0
      return(z)
    }
    z <- as.matrix(x)
    if (!identical(dim(z), dim(q))) {
      stop(name, " must be scalar, cell-level, lineage-level, or have occupancy dimensions.",
           call. = FALSE)
    }
    storage.mode(z) <- "double"
    z[!is.finite(z)] <- 0
    z
  }
  b_col <- as_effect_matrix(colonisation_biotic, "colonisation_biotic")
  b_per <- as_effect_matrix(persistence_biotic, "persistence_biotic")
  est_intercept <- as_effect_matrix(establishment_intercept,
                                    "establishment_intercept")
  est_slope <- as_effect_matrix(establishment_slope,
                                "establishment_slope")
  per_intercept <- as_effect_matrix(persistence_intercept,
                                    "persistence_intercept")
  per_slope <- as_effect_matrix(persistence_slope, "persistence_slope")
  source_movement_multiplier <- as_effect_matrix(
    movement_multiplier, "movement_multiplier"
  )
  if (any(!is.finite(source_movement_multiplier)) ||
      any(source_movement_multiplier < 0)) {
    stop("movement_multiplier must be finite and non-negative.", call. = FALSE)
  }
  arrival <- as.matrix(movement_matrix %*% (q * source_movement_multiplier))
  arrival[!is.finite(arrival) | arrival < 0] <- 0

  est <- hee_colonisation_establishment_probability(
    as.numeric(eta + b_col), intercept = as.numeric(est_intercept),
    slope = as.numeric(est_slope)
  )
  establishment <- matrix(est, nrow(q), ncol(q), dimnames = dimnames(q))
  per <- hee_persistence_probability(
    as.numeric(eta + b_per), intercept = as.numeric(per_intercept),
    slope = as.numeric(per_slope)
  )
  persistence <- matrix(per, nrow(q), ncol(q), dimnames = dimnames(q))
  lambda_colonisation <- arrival * establishment
  lambda_loss <- -log(pmax(persistence, .Machine$double.eps)) /
    persistence_reference_interval
  total_rate <- lambda_colonisation + lambda_loss
  p01 <- ifelse(total_rate > 0,
                lambda_colonisation / total_rate *
                  (1 - exp(-total_rate * delta_t)), 0)
  p11 <- ifelse(total_rate > 0,
                lambda_colonisation / total_rate +
                  lambda_loss / total_rate * exp(-total_rate * delta_t), 1)
  p_next <- pmin(pmax((1 - q) * p01 + q * p11, 0), 1)

  if (isTRUE(stochastic)) {
    if (!is.null(seed)) {
      seed <- suppressWarnings(as.integer(seed)[1L])
      if (!is.finite(seed)) stop("seed must be finite when supplied.", call. = FALSE)
      old_seed <- if (exists(".Random.seed", envir = .GlobalEnv,
                             inherits = FALSE)) {
        get(".Random.seed", envir = .GlobalEnv)
      } else NULL
      on.exit({
        if (is.null(old_seed)) {
          if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv)
          }
        } else {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        }
      }, add = TRUE)
      set.seed(seed)
    }
    q_next <- matrix(stats::rbinom(length(p_next), size = 1L,
                                   prob = as.numeric(p_next)),
                     nrow(q), ncol(q), dimnames = dimnames(q))
  } else {
    q_next <- p_next
  }
  if (!is.null(habitat_state)) {
    h <- as_effect_matrix(habitat_state, "habitat_state")
    q_next[!is.finite(h) | h <= 0] <- 0
  }
  list(
    occupancy_next = q_next,
    arrival_hazard = arrival,
    arrival_probability_delta_t = pmin(
      pmax(1 - exp(-arrival * delta_t), 0), 1
    ),
    arrival_probability_reporting_interval = pmin(
      pmax(1 - exp(-arrival * reporting_interval_myr), 0), 1
    ),
    establishment_probability = establishment,
    colonisation_hazard = lambda_colonisation,
    colonisation_probability_delta_t = pmin(
      pmax(1 - exp(-lambda_colonisation * delta_t), 0), 1
    ),
    colonisation_probability_reporting_interval = pmin(
      pmax(1 - exp(-lambda_colonisation * reporting_interval_myr), 0), 1
    ),
    # Compatibility aliases.  New analysis code should use one of the
    # interval-explicit fields above; the old names were ambiguous whenever
    # the numerical CTMC step was not exactly one Myr.
    arrival_probability = pmin(pmax(1 - exp(-arrival * delta_t), 0), 1),
    colonisation_probability = pmin(
      pmax(1 - exp(-lambda_colonisation * delta_t), 0), 1
    ),
    persistence_probability = persistence,
    source_movement_multiplier = source_movement_multiplier,
    local_extinction_probability = 1 - persistence,
    local_loss_hazard = lambda_loss,
    delta_t_myr = delta_t,
    reporting_interval_myr = reporting_interval_myr,
    persistence_reference_interval = persistence_reference_interval,
    transition_p01 = p01,
    transition_p11 = p11,
    stochastic = isTRUE(stochastic)
  )
}
