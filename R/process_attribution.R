#' Attribute HmscEE dynamic results to the eight core processes
#'
#' `hee_process_attribution()` implements a history-integrated hierarchical
#' attribution framework for HmscEE dynamic occupancy outputs. It deliberately
#' does not create a single flat "six process percentage" table, because the
#' eight named HmscEE processes operate at different causal levels.
#'
#' The returned object separates three non-overlapping layers:
#'
#' * state transitions: colonisation and persistence are decomposed exactly
#'   from the continuous-time occupancy transition model;
#' * upstream mechanisms: environmental filtering, dispersal, biotic filtering
#'   and evolution are attributed with counterfactual Shapley values when the
#'   object provides a simulator;
#' * lineage events: speciation and lineage extinction are summarised from
#'   event jumps or supplied counterfactual event effects.
#'
#' Processes that are not identifiable from the supplied object are reported as
#' `NA` with an explicit status. They are not silently converted to zero.
#'
#' @param object A HmscEE dynamic result, a list containing a state table, or a
#'   state table. The state table must contain `time_ma`, `lineage_id`,
#'   `cell_id`, occupancy probability `q` (or an accepted alias), and the
#'   colonisation and local-loss hazards (`lambda_C`, `lambda_L`, or accepted
#'   aliases) for exact state-transition budgets.
#' @param target Target metric. One of `"occupancy"`, `"range_area"`,
#'   `"lineage_richness"`, `"phylogenetic_diversity"`,
#'   `"functional_diversity"`, or `"turnover"`, or a custom function.
#' @param target_args List of target-specific settings. Important entries are
#'   `aggregation`, `cell_weights`, `group_metric`, `spatial_summary`,
#'   `branch_lengths`, `traits`, `approximation`, `capacity`, and
#'   `target_function`.
#' @param taxa Optional lineage/taxon filter. A character vector is matched to
#'   `lineage_id`. A list may contain `ids`.
#' @param cells Optional cell filter. A vector is matched to `cell_id`. A list
#'   may contain `ids`, `weights`, and `availability`.
#' @param time Optional time filter. `NULL` uses all times, a single number
#'   requires an exact time, two numbers select the closed interval between
#'   them, and a list may contain `from`, `to`, `at`, and `snap`.
#' @param integrate_time Logical; if `TRUE`, return history-integrated
#'   summaries using actual interval durations.
#' @param time_weight `"duration"`, `"equal"`, a numeric vector, or a function
#'   of `(time_start, time_end, metadata)`.
#' @param normalise `"none"`, `"time"`, or `"capacity"`.
#' @param posterior_draws Optional posterior draw selector. `NULL` uses all
#'   available draws; a vector selects draw IDs; a single integer smaller than
#'   the number of draws samples that many draws with `seed`.
#' @param upstream_processes `"available"` or a character vector drawn from
#'   `"environmental_filtering"`, `"dispersal"`, `"biotic_filtering"`, and
#'   `"evolution"` (short aliases `"environment"`, `"biotic"` are accepted).
#' @param reference Counterfactual reference definition passed to the simulator.
#' @param event_effect `"instantaneous"` or `"counterfactual"` for lineage
#'   events. The current implementation uses supplied event deltas for both and
#'   labels unavailable counterfactual event effects as not estimated.
#' @param predictive_support Logical; if `TRUE`, include predictive-support
#'   records when supplied by the object.
#' @param exact_shapley Logical; use exact coalition enumeration when possible.
#' @param n_permutations Number of Monte Carlo permutations when
#'   `exact_shapley = FALSE`.
#' @param common_random_numbers Logical; when simulating coalitions, reuse the
#'   same base seed for paired process switches.
#' @param seed Integer seed used for draw sampling and counterfactual runs.
#' @param parallel Logical placeholder for future parallel execution. The base
#'   implementation is serial to keep package checks deterministic.
#' @param workers Optional worker count recorded in metadata.
#' @param output_dir Optional directory recorded in metadata; heavy simulators
#'   may use it for their own checkpointing.
#' @param return Character vector controlling returned components. Supported
#'   entries are `"summary"`, `"time_series"`, `"maps"`, and `"draws"`.
#' @param ... Reserved for future extensions.
#'
#' @return An object of class `hee_process_attribution`, a list with `summary`,
#'   `state_budget`, `time_series`, `maps`, `draws`, `events`,
#'   `predictive_support`, and `metadata`.
#' @export
#'
#' @examples
#' state <- data.frame(
#'   posterior_draw_id = 1,
#'   lineage_id = "sp1",
#'   cell_id = "c1",
#'   time_ma = c(1, 0),
#'   q = c(0.25, 0.40),
#'   lambda_C = c(0.3, 0.3),
#'   lambda_L = c(0.1, 0.1),
#'   cell_area_km2 = 10
#' )
#' hee_process_attribution(state, target = "occupancy")
hee_process_attribution <- function(
    object,
    target = "occupancy",
    target_args = list(),
    taxa = NULL,
    cells = NULL,
    time = NULL,
    integrate_time = TRUE,
    time_weight = "duration",
    normalise = "none",
    posterior_draws = NULL,
    upstream_processes = "available",
    reference = NULL,
    event_effect = "instantaneous",
    predictive_support = FALSE,
    exact_shapley = TRUE,
    n_permutations = 128L,
    common_random_numbers = TRUE,
    seed = 1L,
    parallel = TRUE,
    workers = NULL,
    output_dir = NULL,
    return = c("summary", "time_series"),
    ...) {
  if (is.function(target)) {
    target_name <- "custom"
    target_args$target_function <- target
  } else {
    target <- match.arg(
      as.character(target)[1],
      c("occupancy", "range_area", "lineage_richness",
        "phylogenetic_diversity", "functional_diversity", "turnover")
    )
    target_name <- target
  }
  normalise <- match.arg(as.character(normalise)[1],
                         c("none", "time", "capacity"))
  event_effect <- match.arg(as.character(event_effect)[1],
                            c("instantaneous", "counterfactual"))
  return <- unique(as.character(return))
  bad_return <- setdiff(return, c("summary", "time_series", "maps", "draws"))
  if (length(bad_return) > 0L) {
    stop("Unsupported return component(s): ",
         paste(bad_return, collapse = ", "), call. = FALSE)
  }
  if (is.null(target_args)) target_args <- list()
  if (!is.list(target_args)) {
    stop("target_args must be a list.", call. = FALSE)
  }
  seed <- as.integer(seed)[1]
  if (!is.finite(seed)) seed <- 1L
  n_permutations <- as.integer(n_permutations)[1]
  if (!is.finite(n_permutations) || n_permutations < 1L) {
    stop("n_permutations must be a positive integer.", call. = FALSE)
  }

  raw_state <- .hee_pa_state_table(object)
  interval_input <- .hee_pa_is_interval_table(raw_state)
  if (interval_input) {
    state <- NULL
    interval <- .hee_pa_standardise_interval(raw_state)
    interval <- .hee_pa_filter_interval(interval, taxa = taxa, cells = cells,
                                        time = time)
    interval <- .hee_pa_select_draws_interval(interval, posterior_draws,
                                              seed = seed)
    interval <- .hee_pa_prepare_target_interval(interval, target_name,
                                                target_args)
    if (nrow(interval) == 0L) {
      stop("No interval rows remain after taxa/cells/time/draw filtering.",
           call. = FALSE)
    }
  } else {
    state <- .hee_pa_standardise_state(raw_state)
    state <- .hee_pa_filter_state(state, taxa = taxa, cells = cells, time = time)
    state <- .hee_pa_select_draws(state, posterior_draws, seed = seed)
    state <- .hee_pa_prepare_target_state(state, target_name, target_args)
    interval <- NULL
    if (nrow(state) == 0L) {
      stop("No state rows remain after taxa/cells/time/draw filtering.",
           call. = FALSE)
    }
  }

  metadata <- list(
    target = target_name,
    target_args = target_args,
    integrate_time = isTRUE(integrate_time),
    time_weight = if (is.character(time_weight)) time_weight[1] else
      class(time_weight)[1],
    normalise = normalise,
    event_effect = event_effect,
    exact_shapley = isTRUE(exact_shapley),
    n_permutations = n_permutations,
    common_random_numbers = isTRUE(common_random_numbers),
    seed = seed,
    parallel_requested = isTRUE(parallel),
    workers = workers,
    output_dir = output_dir,
    n_state_rows = if (interval_input) nrow(interval) else nrow(state),
    n_draws = if (interval_input) {
      length(unique(interval$posterior_draw_id))
    } else {
      length(unique(state$posterior_draw_id))
    },
    n_lineages = if (interval_input) {
      length(unique(interval$lineage_id))
    } else {
      length(unique(state$lineage_id))
    },
    n_cells = if (interval_input) {
      length(unique(interval$cell_id))
    } else {
      length(unique(state$cell_id))
    },
    time_range_ma = if (interval_input) {
      range(c(interval$time_start_ma, interval$time_end_ma), na.rm = TRUE)
    } else {
      range(state$time_ma, na.rm = TRUE)
    },
    interval_input = interval_input
  )

  state_budget <- if (interval_input) {
    .hee_pa_state_budget_from_interval(
      interval,
      target = target_name,
      target_args = target_args,
      integrate_time = integrate_time,
      time_weight = time_weight,
      normalise = normalise,
      metadata = metadata
    )
  } else {
    .hee_pa_state_budget(
      state,
      target = target_name,
      target_args = target_args,
      integrate_time = integrate_time,
      time_weight = time_weight,
      normalise = normalise,
      metadata = metadata
    )
  }

  upstream <- .hee_pa_upstream_summary(
    object = object,
    state = if (interval_input) interval else state,
    target = target_name,
    target_args = target_args,
    upstream_processes = upstream_processes,
    reference = reference,
    exact_shapley = exact_shapley,
    n_permutations = n_permutations,
    common_random_numbers = common_random_numbers,
    seed = seed,
    metadata = metadata
  )

  events <- .hee_pa_event_summary(
    object = object,
    target = target_name,
    target_args = target_args,
    time = time,
    event_effect = event_effect
  )

  support <- .hee_pa_predictive_support(object, predictive_support)

  summary <- .hee_pa_bind_summary(
    state_budget$summary,
    upstream$summary,
    events$summary,
    support$summary
  )

  out <- list(
    target = target_name,
    summary = summary,
    state_budget = state_budget$state_budget,
    time_series = .hee_pa_requested("time_series", return,
                                    state_budget$time_series),
    maps = .hee_pa_requested("maps", return, state_budget$maps),
    draws = .hee_pa_requested("draws", return,
                              .hee_pa_bind_rows(list(state_budget$draws,
                                                     upstream$draws))),
    events = events$events,
    predictive_support = support$predictive_support,
    metadata = metadata
  )
  class(out) <- c("hee_process_attribution", "list")
  .hee_pa_validate_result(out)
  out
}

#' @export
print.hee_process_attribution <- function(x, ...) {
  cat("HmscEE process attribution\n")
  cat("Target:", x$target, "\n")
  if (!is.null(x$metadata$n_draws)) {
    cat("Draws:", x$metadata$n_draws,
        " Lineages:", x$metadata$n_lineages,
        " Cells:", x$metadata$n_cells, "\n")
  }
  if (!is.null(x$summary) && nrow(x$summary) > 0L) {
    cols <- intersect(c("process", "layer", "contribution_type",
                        "estimate", "share_within_layer", "status"),
                      names(x$summary))
    print(utils::head(x$summary[, cols, drop = FALSE], 12L),
          row.names = FALSE)
    if (nrow(x$summary) > 12L) {
      cat("... ", nrow(x$summary) - 12L, " more rows\n", sep = "")
    }
  }
  invisible(x)
}

.hee_pa_requested <- function(name, requested, value) {
  if (name %in% requested) value else NULL
}

.hee_pa_state_table <- function(object) {
  if (is.data.frame(object)) return(as.data.frame(object))
  if (!is.list(object)) {
    stop("object must be a data.frame or a list-like HmscEE dynamic result.",
         call. = FALSE)
  }
  candidates <- c("intervals", "interval_state", "attribution_state",
                  "state", "state_table", "trajectory", "occupancy",
                  "dynamic_state", "posterior_state")
  for (nm in candidates) {
    if (!is.null(object[[nm]]) && is.data.frame(object[[nm]])) {
      return(as.data.frame(object[[nm]]))
    }
    if (!is.null(object[[nm]]) && is.list(object[[nm]])) {
      nested <- object[[nm]]
      for (nm2 in c("q", "state", "table", "trajectory")) {
        if (!is.null(nested[[nm2]]) && is.data.frame(nested[[nm2]])) {
          return(as.data.frame(nested[[nm2]]))
        }
      }
    }
  }
  stop("object does not contain a state table. Expected object$state, ",
       "object$state_table, object$trajectory, or a data.frame input.",
       call. = FALSE)
}

.hee_pa_standardise_state <- function(x) {
  x <- as.data.frame(x)
  x <- .rename_first_existing(x, "posterior_draw_id",
                              c("draw_id", "state_draw_id", "particle_id",
                                "posterior_draw", "response_draw"))
  x <- .rename_first_existing(x, "lineage_id",
                              c("lineage", "species", "taxon", "taxon_id"))
  x <- .rename_first_existing(x, "cell_id",
                              c("cell", "paleo_cell_id", "grid_id"))
  x <- .rename_first_existing(x, "q",
                              c("occupancy_probability", "probability",
                                "psi", "p", "mean_occupancy_probability"))
  x <- .rename_first_existing(x, "lambda_C",
                              c("colonisation_hazard", "colonization_hazard",
                                "lambda_c", "colonisation_rate",
                                "colonization_rate", "lambda_colonisation"))
  x <- .rename_first_existing(x, "lambda_L",
                              c("local_loss_hazard", "local_extinction_hazard",
                                "lambda_l", "loss_hazard",
                                "lambda_local_loss", "local_loss_rate"))
  x <- .rename_first_existing(x, "cell_area_km2",
                              c("cell_area", "area_km2", "land_area_km2",
                                "area"))
  x <- .rename_first_existing(x, "z",
                              c("state", "occupied", "occupancy_state"))
  .require_cols(x, c("time_ma", "lineage_id", "cell_id", "q"),
                "state")
  if (!"posterior_draw_id" %in% names(x)) x$posterior_draw_id <- 1L
  if (!"lambda_C" %in% names(x)) x$lambda_C <- NA_real_
  if (!"lambda_L" %in% names(x)) x$lambda_L <- NA_real_
  if (!"cell_area_km2" %in% names(x)) x$cell_area_km2 <- 1
  x$posterior_draw_id <- as.character(x$posterior_draw_id)
  x$lineage_id <- as.character(x$lineage_id)
  x$cell_id <- as.character(x$cell_id)
  x$time_ma <- suppressWarnings(as.numeric(x$time_ma))
  x$q <- suppressWarnings(as.numeric(x$q))
  x$lambda_C <- suppressWarnings(as.numeric(x$lambda_C))
  x$lambda_L <- suppressWarnings(as.numeric(x$lambda_L))
  x$cell_area_km2 <- suppressWarnings(as.numeric(x$cell_area_km2))
  x$cell_area_km2[!is.finite(x$cell_area_km2) |
                    x$cell_area_km2 < 0] <- 1
  if ("land_exists" %in% names(x)) {
    land <- suppressWarnings(as.numeric(x$land_exists))
    x <- x[is.na(land) | land > 0, , drop = FALSE]
  }
  bad <- !is.finite(x$time_ma) | !is.finite(x$q)
  if (any(bad)) x <- x[!bad, , drop = FALSE]
  if (nrow(x) == 0L) {
    stop("state has no finite time_ma and q rows.", call. = FALSE)
  }
  x$q <- pmin(pmax(x$q, 0), 1)
  x
}

.hee_pa_is_interval_table <- function(x) {
  nms <- names(as.data.frame(x))
  any(c("time_start_ma", "time_from_ma", "older_time_ma") %in% nms) &&
    any(c("time_end_ma", "time_to_ma", "younger_time_ma") %in% nms) &&
    any(c("q0", "q_start", "occupancy_start", "q_before") %in% nms)
}

.hee_pa_standardise_interval <- function(x) {
  x <- as.data.frame(x)
  x <- .rename_first_existing(x, "posterior_draw_id",
                              c("draw_id", "response_draw", "particle_id",
                                "state_draw_id", "posterior_draw"))
  x <- .rename_first_existing(x, "lineage_id",
                              c("lineage", "species", "taxon", "taxon_id"))
  x <- .rename_first_existing(x, "cell_id",
                              c("cell", "paleo_cell_id", "grid_id"))
  x <- .rename_first_existing(x, "time_start_ma",
                              c("time_from_ma", "older_time_ma",
                                "previous_time_ma", "start_ma"))
  x <- .rename_first_existing(x, "time_end_ma",
                              c("time_to_ma", "younger_time_ma",
                                "current_time_ma", "end_ma"))
  x <- .rename_first_existing(x, "q0",
                              c("q_start", "occupancy_start", "q_before",
                                "previous_q"))
  x <- .rename_first_existing(x, "q_observed_next",
                              c("q1", "q_end", "occupancy_end", "q_after",
                                "current_q", "q_next"))
  x <- .rename_first_existing(x, "lambda_C",
                              c("colonisation_hazard", "colonization_hazard",
                                "lambda_c", "colonisation_rate",
                                "colonization_rate", "lambda_colonisation"))
  x <- .rename_first_existing(x, "lambda_L",
                              c("local_loss_hazard", "local_extinction_hazard",
                                "lambda_l", "loss_hazard",
                                "lambda_local_loss", "local_loss_rate"))
  x <- .rename_first_existing(x, "cell_area_km2",
                              c("cell_area", "area_km2", "land_area_km2",
                                "area"))
  .require_cols(x, c("time_start_ma", "time_end_ma", "lineage_id",
                     "cell_id", "q0", "lambda_C", "lambda_L"),
                "interval_state")
  if (!"posterior_draw_id" %in% names(x)) x$posterior_draw_id <- 1L
  if (!"cell_area_km2" %in% names(x)) x$cell_area_km2 <- 1
  if (!"q_observed_next" %in% names(x)) x$q_observed_next <- NA_real_
  x$posterior_draw_id <- as.character(x$posterior_draw_id)
  x$lineage_id <- as.character(x$lineage_id)
  x$cell_id <- as.character(x$cell_id)
  for (cc in c("time_start_ma", "time_end_ma", "q0", "q_observed_next",
               "lambda_C", "lambda_L", "cell_area_km2")) {
    x[[cc]] <- suppressWarnings(as.numeric(x[[cc]]))
  }
  x$delta_t <- abs(x$time_start_ma - x$time_end_ma)
  bad <- !is.finite(x$time_start_ma) | !is.finite(x$time_end_ma) |
    !is.finite(x$delta_t) | x$delta_t <= 0 | !is.finite(x$q0)
  if (any(bad)) x <- x[!bad, , drop = FALSE]
  if (nrow(x) == 0L) {
    stop("interval_state has no finite positive-duration rows.",
         call. = FALSE)
  }
  swap <- x$time_start_ma < x$time_end_ma
  if (any(swap)) {
    old <- x$time_start_ma[swap]
    x$time_start_ma[swap] <- x$time_end_ma[swap]
    x$time_end_ma[swap] <- old
  }
  x$time_mid_ma <- (x$time_start_ma + x$time_end_ma) / 2
  x$q0 <- pmin(pmax(x$q0, 0), 1)
  x$q_observed_next <- pmin(pmax(x$q_observed_next, 0), 1)
  x$lambda_C[!is.finite(x$lambda_C) | x$lambda_C < 0] <- NA_real_
  x$lambda_L[!is.finite(x$lambda_L) | x$lambda_L < 0] <- NA_real_
  x$cell_area_km2[!is.finite(x$cell_area_km2) |
                    x$cell_area_km2 < 0] <- 1
  rownames(x) <- NULL
  x
}

.hee_pa_filter_interval <- function(interval, taxa = NULL, cells = NULL,
                                    time = NULL) {
  if (!is.null(taxa)) {
    ids <- if (is.list(taxa)) taxa$ids else taxa
    if (is.null(ids)) stop("taxa list must contain ids.", call. = FALSE)
    interval <- interval[interval$lineage_id %in% as.character(ids),
                         , drop = FALSE]
  }
  if (!is.null(cells)) {
    ids <- if (is.list(cells)) cells$ids else cells
    if (is.null(ids)) stop("cells list must contain ids.", call. = FALSE)
    interval <- interval[interval$cell_id %in% as.character(ids),
                         , drop = FALSE]
  }
  if (!is.null(time)) {
    snap <- "exact"
    if (is.list(time)) {
      snap <- if (!is.null(time$snap)) as.character(time$snap)[1] else "exact"
      if (!is.null(time$at)) time <- time$at else time <- c(time$from, time$to)
    }
    all_times <- sort(unique(c(interval$time_start_ma, interval$time_end_ma)),
                      decreasing = TRUE)
    if (length(time) == 1L) {
      t0 <- suppressWarnings(as.numeric(time))[1]
      if (!any(abs(all_times - t0) < 1e-8)) {
        if (identical(snap, "nearest")) {
          t0 <- all_times[which.min(abs(all_times - t0))]
        } else {
          stop("Requested time is not an exact interval boundary. ",
               "Use time = list(at = ..., snap = 'nearest') for nearest ",
               "matching.", call. = FALSE)
        }
      }
      interval <- interval[abs(interval$time_start_ma - t0) < 1e-8 |
                             abs(interval$time_end_ma - t0) < 1e-8,
                           , drop = FALSE]
    } else if (length(time) >= 2L) {
      tr <- range(suppressWarnings(as.numeric(time[1:2])), na.rm = TRUE)
      interval <- interval[
        interval$time_start_ma >= tr[1] & interval$time_end_ma <= tr[2],
        , drop = FALSE
      ]
    }
  }
  rownames(interval) <- NULL
  interval
}

.hee_pa_select_draws_interval <- function(interval, posterior_draws = NULL,
                                          seed = 1L) {
  ids <- unique(interval$posterior_draw_id)
  if (is.null(posterior_draws)) return(interval)
  if (is.list(posterior_draws)) {
    if (!is.null(posterior_draws$ids)) {
      keep <- as.character(posterior_draws$ids)
    } else if (!is.null(posterior_draws$n)) {
      set.seed(seed)
      n <- min(as.integer(posterior_draws$n)[1], length(ids))
      keep <- sample(ids, n)
    } else {
      return(interval)
    }
  } else if (length(posterior_draws) == 1L && is.numeric(posterior_draws) &&
             as.integer(posterior_draws)[1] < length(ids) &&
             as.integer(posterior_draws)[1] > 0L) {
    set.seed(seed)
    keep <- sample(ids, as.integer(posterior_draws)[1])
  } else {
    keep <- as.character(posterior_draws)
  }
  interval[interval$posterior_draw_id %in% keep, , drop = FALSE]
}

.hee_pa_prepare_target_interval <- function(interval, target, target_args) {
  if (target == "phylogenetic_diversity" &&
      !"branch_length" %in% names(interval)) {
    bl <- target_args$branch_lengths
    if (is.null(bl)) {
      stop("target = 'phylogenetic_diversity' requires a branch_length ",
           "column or target_args$branch_lengths.", call. = FALSE)
    }
    interval <- .hee_pa_add_branch_lengths(interval, bl)
  }
  if (target %in% c("functional_diversity", "turnover")) {
    stop("Interval-state attribution currently supports exact state budgets ",
         "for occupancy, range_area, lineage_richness, and branch-weighted ",
         "phylogenetic_diversity. Use full joint state draws for ", target,
         ".", call. = FALSE)
  }
  interval
}

.hee_pa_filter_state <- function(state, taxa = NULL, cells = NULL, time = NULL) {
  if (!is.null(taxa)) {
    ids <- if (is.list(taxa)) taxa$ids else taxa
    if (is.null(ids)) stop("taxa list must contain ids.", call. = FALSE)
    state <- state[state$lineage_id %in% as.character(ids), , drop = FALSE]
  }
  if (!is.null(cells)) {
    ids <- if (is.list(cells)) cells$ids else cells
    if (is.null(ids)) stop("cells list must contain ids.", call. = FALSE)
    state <- state[state$cell_id %in% as.character(ids), , drop = FALSE]
  }
  if (!is.null(time)) {
    snap <- "exact"
    if (is.list(time)) {
      snap <- if (!is.null(time$snap)) as.character(time$snap)[1] else "exact"
      if (!is.null(time$at)) {
        time <- time$at
      } else {
        time <- c(time$from, time$to)
      }
    }
    tt <- sort(unique(state$time_ma), decreasing = TRUE)
    if (length(time) == 1L) {
      t0 <- suppressWarnings(as.numeric(time))[1]
      if (!any(abs(tt - t0) < 1e-8)) {
        if (identical(snap, "nearest")) {
          t0 <- tt[which.min(abs(tt - t0))]
        } else {
          stop("Requested time is not an exact time_ma in object. ",
               "Use time = list(at = ..., snap = 'nearest') for nearest ",
               "matching.", call. = FALSE)
        }
      }
      state <- state[abs(state$time_ma - t0) < 1e-8, , drop = FALSE]
    } else if (length(time) >= 2L) {
      tr <- range(suppressWarnings(as.numeric(time[1:2])), na.rm = TRUE)
      state <- state[state$time_ma >= tr[1] & state$time_ma <= tr[2],
                     , drop = FALSE]
    }
  }
  rownames(state) <- NULL
  state
}

.hee_pa_select_draws <- function(state, posterior_draws = NULL, seed = 1L) {
  ids <- unique(state$posterior_draw_id)
  if (is.null(posterior_draws)) return(state)
  if (is.list(posterior_draws)) {
    if (!is.null(posterior_draws$ids)) {
      keep <- as.character(posterior_draws$ids)
    } else if (!is.null(posterior_draws$n)) {
      set.seed(seed)
      n <- min(as.integer(posterior_draws$n)[1], length(ids))
      keep <- sample(ids, n)
    } else {
      return(state)
    }
  } else {
    pd <- posterior_draws
    if (length(pd) == 1L && is.numeric(pd) &&
        as.integer(pd)[1] < length(ids) && as.integer(pd)[1] > 0L) {
      set.seed(seed)
      keep <- sample(ids, as.integer(pd)[1])
    } else {
      keep <- as.character(pd)
    }
  }
  state[state$posterior_draw_id %in% keep, , drop = FALSE]
}

.hee_pa_prepare_target_state <- function(state, target, target_args) {
  if (target %in% c("phylogenetic_diversity", "functional_diversity",
                    "turnover")) {
    if (!"z" %in% names(state)) {
      approximation <- target_args$approximation %||% "none"
      if (!identical(approximation, "independent_bernoulli")) {
        stop(target, " requires joint binary state draws in column z. ",
             "Use target_args = list(approximation = ",
             "'independent_bernoulli') only for an explicit approximation.",
             call. = FALSE)
      }
      state$z <- as.numeric(state$q >= (target_args$threshold %||% 0.5))
      warning("Using independent Bernoulli/threshold approximation for ",
              target, "; this is not a strict joint posterior-state metric.",
              call. = FALSE)
    }
    state$z <- suppressWarnings(as.numeric(state$z))
    state$z[!is.finite(state$z)] <- 0
    state$z <- as.numeric(state$z > 0)
  }
  if (target == "phylogenetic_diversity") {
    if (!"branch_length" %in% names(state)) {
      bl <- target_args$branch_lengths
      if (is.null(bl)) {
        stop("target = 'phylogenetic_diversity' requires a branch_length ",
             "column or target_args$branch_lengths.", call. = FALSE)
      }
      state <- .hee_pa_add_branch_lengths(state, bl)
    }
  }
  state
}

.hee_pa_add_branch_lengths <- function(state, branch_lengths) {
  if (is.data.frame(branch_lengths)) {
    bl <- as.data.frame(branch_lengths)
    bl <- .rename_first_existing(bl, "lineage_id",
                                 c("lineage", "species", "taxon"))
    bl <- .rename_first_existing(bl, "branch_length",
                                 c("length", "branch_length_myr",
                                   "pd_weight"))
    .require_cols(bl, c("lineage_id", "branch_length"), "branch_lengths")
    bl$lineage_id <- as.character(bl$lineage_id)
    if ("time_ma" %in% names(bl)) {
      bl$time_ma <- suppressWarnings(as.numeric(bl$time_ma))
      state <- merge(state, bl[, c("lineage_id", "time_ma", "branch_length"),
                               drop = FALSE],
                     by = c("lineage_id", "time_ma"), all.x = TRUE,
                     sort = FALSE)
    } else {
      state <- merge(state, bl[, c("lineage_id", "branch_length"),
                               drop = FALSE],
                     by = "lineage_id", all.x = TRUE, sort = FALSE)
    }
  } else {
    vals <- suppressWarnings(as.numeric(branch_lengths))
    names(vals) <- names(branch_lengths)
    state$branch_length <- vals[match(state$lineage_id, names(vals))]
  }
  state$branch_length <- suppressWarnings(as.numeric(state$branch_length))
  state$branch_length[!is.finite(state$branch_length) |
                        state$branch_length < 0] <- 0
  state
}

.hee_pa_state_budget_from_interval <- function(interval, target, target_args,
                                               integrate_time, time_weight,
                                               normalise, metadata) {
  if (is.null(interval) || nrow(interval) == 0L) {
    empty <- .hee_pa_empty_summary("state_transition",
                                   "insufficient_time_points")
    return(list(summary = empty, state_budget = interval,
                time_series = interval, maps = NULL, draws = NULL))
  }
  interval <- .hee_pa_add_transition_components(interval)
  interval <- .hee_pa_add_target_weights(interval, target, target_args)
  interval$interval_weight <- .hee_pa_interval_weights(
    interval, time_weight, metadata
  )
  interval$capacity <- .hee_pa_interval_capacity(
    interval, target, target_args, normalise
  )

  metric_components <- .hee_pa_component_metrics(
    interval, target, integrate_time, normalise
  )
  draw_components <- .hee_pa_draw_component_metrics(
    interval, target, integrate_time, normalise
  )
  time_series <- .hee_pa_time_series_components(interval, target)
  maps <- .hee_pa_map_components(interval)

  summary <- .hee_pa_state_summary(metric_components)
  list(summary = summary, state_budget = metric_components,
       time_series = time_series, maps = maps, draws = draw_components)
}

.hee_pa_state_budget <- function(state, target, target_args,
                                 integrate_time, time_weight, normalise,
                                 metadata) {
  interval <- .hee_pa_intervals(state)
  .hee_pa_state_budget_from_interval(
    interval = interval,
    target = target,
    target_args = target_args,
    integrate_time = integrate_time,
    time_weight = time_weight,
    normalise = normalise,
    metadata = metadata
  )
}

.hee_pa_intervals <- function(state) {
  key <- paste(state$posterior_draw_id, state$lineage_id, state$cell_id,
               sep = "\r")
  rows <- split(seq_len(nrow(state)), key)
  out <- lapply(rows, function(idx) {
    d <- state[idx, , drop = FALSE]
    d <- d[order(d$time_ma, decreasing = TRUE), , drop = FALSE]
    if (nrow(d) < 2L) return(NULL)
    older <- d[-nrow(d), , drop = FALSE]
    younger <- d[-1L, , drop = FALSE]
    delta_t <- older$time_ma - younger$time_ma
    keep <- is.finite(delta_t) & delta_t > 0
    if (!any(keep)) return(NULL)
    out <- data.frame(
      posterior_draw_id = older$posterior_draw_id[keep],
      lineage_id = older$lineage_id[keep],
      cell_id = older$cell_id[keep],
      time_start_ma = older$time_ma[keep],
      time_end_ma = younger$time_ma[keep],
      time_mid_ma = (older$time_ma[keep] + younger$time_ma[keep]) / 2,
      delta_t = delta_t[keep],
      q0 = older$q[keep],
      q_observed_next = younger$q[keep],
      lambda_C = older$lambda_C[keep],
      lambda_L = older$lambda_L[keep],
      cell_area_km2 = older$cell_area_km2[keep],
      stringsAsFactors = FALSE
    )
    if ("branch_length" %in% names(older)) {
      out$branch_length <- suppressWarnings(as.numeric(older$branch_length[keep]))
      out$branch_length[!is.finite(out$branch_length) |
                          out$branch_length < 0] <- 0
    }
    out
  })
  .hee_pa_bind_rows(out)
}

.hee_pa_add_transition_components <- function(interval) {
  interval$lambda_C[!is.finite(interval$lambda_C) |
                      interval$lambda_C < 0] <- NA_real_
  interval$lambda_L[!is.finite(interval$lambda_L) |
                      interval$lambda_L < 0] <- NA_real_
  missing_rate <- is.na(interval$lambda_C) | is.na(interval$lambda_L)
  if (any(missing_rate)) {
    stop("Exact state-transition attribution requires finite non-negative ",
         "lambda_C and lambda_L for every interval.", call. = FALSE)
  }
  r <- interval$lambda_C + interval$lambda_L
  dt <- interval$delta_t
  exp_r <- exp(-r * dt)
  interval$P01 <- 0
  interval$P10 <- 0
  interval$P11 <- 1
  ok <- is.finite(r) & r > 0
  interval$P01[ok] <- interval$lambda_C[ok] / r[ok] * (1 - exp_r[ok])
  interval$P10[ok] <- interval$lambda_L[ok] / r[ok] * (1 - exp_r[ok])
  interval$P11[ok] <- interval$lambda_C[ok] / r[ok] +
    interval$lambda_L[ok] / r[ok] * exp_r[ok]
  interval$P11[!ok] <- 1
  interval$P01 <- pmin(pmax(interval$P01, 0), 1)
  interval$P10 <- pmin(pmax(interval$P10, 0), 1)
  interval$P11 <- pmin(pmax(interval$P11, 0), 1)
  interval$persistent_occupancy <- interval$q0 * exp(-interval$lambda_L * dt)
  interval$new_colonisation <- (1 - interval$q0) * interval$P01
  interval$recolonisation <- interval$q0 *
    pmax(interval$P11 - exp(-interval$lambda_L * dt), 0)
  interval$local_loss <- interval$q0 * interval$P10
  interval$colonisation_gain <- interval$new_colonisation +
    interval$recolonisation
  interval$q_next_model <- interval$persistent_occupancy +
    interval$new_colonisation + interval$recolonisation
  interval$net_change_model <- interval$q_next_model - interval$q0
  interval
}

.hee_pa_add_target_weights <- function(interval, target, target_args) {
  interval$target_weight <- 1
  interval$target_denominator_weight <- 1
  if (target %in% c("occupancy", "lineage_richness")) {
    cell_weights <- target_args$cell_weights %||% "area"
    if (identical(cell_weights, "area")) {
      interval$target_weight <- interval$cell_area_km2
      interval$target_denominator_weight <- interval$cell_area_km2
    }
  }
  if (target == "range_area") {
    interval$target_weight <- interval$cell_area_km2
    interval$target_denominator_weight <- 1
  }
  if (target == "phylogenetic_diversity") {
    if (!"branch_length" %in% names(interval)) {
      stop("target = 'phylogenetic_diversity' requires branch_length for ",
           "state-transition attribution.", call. = FALSE)
    }
    interval$target_weight <- interval$cell_area_km2 * interval$branch_length
    interval$target_denominator_weight <- interval$cell_area_km2
  }
  interval
}

.hee_pa_interval_weights <- function(interval, time_weight, metadata) {
  if (is.character(time_weight)) {
    time_weight <- match.arg(time_weight[1], c("duration", "equal"))
    if (identical(time_weight, "duration")) return(interval$delta_t)
    return(rep(1, nrow(interval)))
  }
  if (is.numeric(time_weight)) {
    tw <- as.numeric(time_weight)
    groups <- unique(paste(interval$time_start_ma, interval$time_end_ma,
                           sep = "\r"))
    if (length(tw) != length(groups)) {
      stop("Numeric time_weight must have one value per unique interval.",
           call. = FALSE)
    }
    key <- paste(interval$time_start_ma, interval$time_end_ma, sep = "\r")
    out <- tw[match(key, groups)]
    out[!is.finite(out) | out < 0] <- 0
    return(out)
  }
  if (is.function(time_weight)) {
    out <- mapply(time_weight, interval$time_start_ma, interval$time_end_ma,
                  MoreArgs = list(metadata = metadata))
    out <- suppressWarnings(as.numeric(out))
    out[!is.finite(out) | out < 0] <- 0
    return(out)
  }
  stop("time_weight must be 'duration', 'equal', numeric, or a function.",
       call. = FALSE)
}

.hee_pa_interval_capacity <- function(interval, target, target_args, normalise) {
  if (!identical(normalise, "capacity")) return(rep(1, nrow(interval)))
  cap <- target_args$capacity
  if (!is.null(cap)) {
    if (is.function(cap)) {
      out <- mapply(cap, interval$time_start_ma, interval$time_end_ma,
                    interval$cell_id, interval$lineage_id)
      out <- suppressWarnings(as.numeric(out))
      out[!is.finite(out) | out <= 0] <- NA_real_
      return(out)
    }
    if (is.numeric(cap) && length(cap) == 1L) {
      return(rep(as.numeric(cap), nrow(interval)))
    }
    if (is.character(cap) && length(cap) == 1L && cap %in% names(interval)) {
      out <- suppressWarnings(as.numeric(interval[[cap]]))
      out[!is.finite(out) | out <= 0] <- NA_real_
      return(out)
    }
  }
  if (identical(target, "occupancy")) return(rep(1, nrow(interval)))
  if (identical(target, "lineage_richness")) {
    n <- length(unique(interval$lineage_id))
    return(rep(max(n, 1), nrow(interval)))
  }
  stop("normalise = 'capacity' for target '", target,
       "' requires target_args$capacity.", call. = FALSE)
}

.hee_pa_component_metrics <- function(interval, target, integrate_time,
                                      normalise) {
  components <- c("persistent_occupancy", "new_colonisation",
                  "recolonisation", "local_loss", "colonisation_gain",
                  "net_change_model")
  rows <- lapply(components, function(comp) {
    estimate <- .hee_pa_integrate_component(interval, comp, target,
                                            integrate_time, normalise)
    data.frame(component = comp, estimate = estimate,
               stringsAsFactors = FALSE)
  })
  .hee_pa_bind_rows(rows)
}

.hee_pa_draw_component_metrics <- function(interval, target, integrate_time,
                                           normalise) {
  rows <- lapply(split(interval, interval$posterior_draw_id, drop = TRUE),
                 function(d) {
    metrics <- .hee_pa_component_metrics(d, target, integrate_time, normalise)
    metrics$posterior_draw_id <- d$posterior_draw_id[1]
    metrics
  })
  .hee_pa_bind_rows(rows)
}

.hee_pa_integrate_component <- function(interval, comp, target,
                                        integrate_time, normalise) {
  if (!comp %in% names(interval)) return(NA_real_)
  keys <- unique(paste(interval$posterior_draw_id, interval$time_start_ma,
                       interval$time_end_ma, sep = "\r"))
  vals <- lapply(keys, function(key) {
    d <- interval[paste(interval$posterior_draw_id, interval$time_start_ma,
                        interval$time_end_ma, sep = "\r") == key,
                  , drop = FALSE]
    num <- sum(d[[comp]] * d$target_weight, na.rm = TRUE)
    den <- 1
    if (target %in% c("occupancy", "lineage_richness")) {
      den <- sum(unique(d[, c("cell_id", "target_denominator_weight"),
                          drop = FALSE])$target_denominator_weight,
                 na.rm = TRUE)
      if (!is.finite(den) || den <= 0) den <- 1
    }
    value <- num / den
    w <- mean(d$interval_weight, na.rm = TRUE)
    cap <- mean(d$capacity, na.rm = TRUE)
    data.frame(value = value, interval_weight = w, capacity = cap,
               stringsAsFactors = FALSE)
  })
  tab <- .hee_pa_bind_rows(vals)
  if (is.null(tab) || nrow(tab) == 0L) return(NA_real_)
  if (!isTRUE(integrate_time)) {
    return(mean(tab$value, na.rm = TRUE))
  }
  out <- sum(tab$value * tab$interval_weight, na.rm = TRUE)
  if (identical(normalise, "time")) {
    denom <- sum(tab$interval_weight, na.rm = TRUE)
    if (denom > 0) out <- out / denom
  } else if (identical(normalise, "capacity")) {
    denom <- sum(tab$capacity * tab$interval_weight, na.rm = TRUE)
    if (denom > 0) out <- out / denom
  }
  out
}

.hee_pa_time_series_components <- function(interval, target) {
  components <- c("persistent_occupancy", "new_colonisation",
                  "recolonisation", "local_loss", "colonisation_gain",
                  "net_change_model")
  keys <- unique(paste(interval$time_start_ma, interval$time_end_ma,
                       sep = "\r"))
  rows <- list()
  for (key in keys) {
    d <- interval[paste(interval$time_start_ma, interval$time_end_ma,
                        sep = "\r") == key, , drop = FALSE]
    for (comp in components) {
      num <- sum(d[[comp]] * d$target_weight, na.rm = TRUE)
      den <- 1
      if (target %in% c("occupancy", "lineage_richness")) {
        den <- sum(unique(d[, c("cell_id", "target_denominator_weight"),
                            drop = FALSE])$target_denominator_weight,
                   na.rm = TRUE)
        if (!is.finite(den) || den <= 0) den <- 1
      }
      rows[[length(rows) + 1L]] <- data.frame(
        time_start_ma = d$time_start_ma[1],
        time_end_ma = d$time_end_ma[1],
        delta_t = d$delta_t[1],
        component = comp,
        value = num / den,
        stringsAsFactors = FALSE
      )
    }
  }
  .hee_pa_bind_rows(rows)
}

.hee_pa_map_components <- function(interval) {
  components <- c("persistent_occupancy", "new_colonisation",
                  "recolonisation", "local_loss", "colonisation_gain",
                  "net_change_model")
  rows <- list()
  for (comp in components) {
    val <- stats::aggregate(
      interval[[comp]] * interval$target_weight * interval$interval_weight,
      interval[, c("cell_id", "time_start_ma", "time_end_ma"), drop = FALSE],
      sum, na.rm = TRUE
    )
    names(val)[ncol(val)] <- "value"
    val$component <- comp
    rows[[comp]] <- val
  }
  .hee_pa_bind_rows(rows)
}

.hee_pa_state_summary <- function(metric_components) {
  get <- function(nm) {
    x <- metric_components$estimate[metric_components$component == nm]
    if (length(x) == 0L) NA_real_ else x[1]
  }
  rows <- data.frame(
    process = c("Retained occupancy", "New occupancy", "Recolonisation",
                "Local loss"),
    process_slug = rep("state_transition", 4),
    layer = rep("state_transition", 4),
    contribution_type = c("retained_occupancy", "new_colonisation",
                          "recolonisation", "local_loss"),
    estimate = c(get("persistent_occupancy"), get("new_colonisation"),
                 get("recolonisation"), get("local_loss")),
    signed_contribution = c(get("persistent_occupancy"),
                            get("new_colonisation"),
                            get("recolonisation"),
                            -get("local_loss")),
    status = "estimated_exact_ctmc_budget",
    interpretation = c(
      "Occupancy retained from already occupied cells.",
      "Previously empty cells occupied by new colonisation.",
      "Cells initially occupied, lost, and recolonised within the interval.",
      "Local loss from previously occupied cells; this is not lineage extinction."
    ),
    stringsAsFactors = FALSE
  )
  rows$absolute_contribution <- abs(rows$signed_contribution)
  source_rows <- rows$contribution_type %in%
    c("retained_occupancy", "new_colonisation", "recolonisation")
  denom <- sum(rows$estimate[source_rows], na.rm = TRUE)
  rows$share_within_layer <- NA_real_
  if (is.finite(denom) && denom > 0) {
    rows$share_within_layer[source_rows] <- rows$estimate[source_rows] / denom
  }
  rows
}

.hee_pa_upstream_summary <- function(object, state, target, target_args,
                                     upstream_processes, reference,
                                     exact_shapley, n_permutations,
                                     common_random_numbers, seed, metadata) {
  requested <- .hee_pa_upstream_processes(upstream_processes)
  status <- .hee_pa_process_status(object, requested)
  simulator <- .hee_pa_simulator(object)
  active <- requested[status[requested] != "not_estimated" &
                        status[requested] != "fixed_off"]
  if (is.null(simulator) || length(active) == 0L) {
    rows <- lapply(requested, function(p) {
      st <- status[[p]]
      if (is.null(st) || is.na(st)) st <- "not_estimated"
      if (p == "biotic_filtering" && st == "not_estimated") {
        st <- "not_estimated_no_external_biotic_evidence"
      } else if (is.null(simulator)) {
        st <- paste0(st, "_no_counterfactual_simulator")
      }
      .hee_pa_upstream_row(p, NA_real_, NA_real_, st,
                           "Counterfactual Shapley attribution is not available.")
    })
    return(list(summary = .hee_pa_bind_rows(rows), draws = NULL))
  }
  if (isTRUE(exact_shapley)) {
    shp <- .hee_pa_exact_shapley(
      simulator = simulator, processes = active, target = target,
      target_args = target_args, reference = reference, seed = seed,
      metadata = metadata
    )
  } else {
    shp <- .hee_pa_mc_shapley(
      simulator = simulator, processes = active, target = target,
      target_args = target_args, reference = reference, seed = seed,
      n_permutations = n_permutations,
      common_random_numbers = common_random_numbers,
      metadata = metadata
    )
  }
  rows <- lapply(requested, function(p) {
    if (p %in% shp$process) {
      z <- shp[shp$process == p, , drop = FALSE]
      .hee_pa_upstream_row(p, z$estimate[1], z$mcse[1],
                           "estimated_counterfactual_shapley",
                           "Model-implied upstream contribution; not a direct causal proof.")
    } else {
      st <- status[[p]]
      if (is.null(st) || is.na(st)) st <- "not_estimated"
      .hee_pa_upstream_row(p, NA_real_, NA_real_, st,
                           "Process was not included in Shapley coalitions.")
    }
  })
  list(summary = .hee_pa_bind_rows(rows), draws = shp)
}

.hee_pa_upstream_processes <- function(upstream_processes) {
  allowed <- c("environmental_filtering", "dispersal",
               "biotic_filtering", "evolution")
  if (identical(upstream_processes, "available")) return(allowed)
  x <- as.character(upstream_processes)
  aliases <- c(environment = "environmental_filtering",
               environmental = "environmental_filtering",
               biotic = "biotic_filtering")
  x <- ifelse(x %in% names(aliases), aliases[x], x)
  bad <- setdiff(x, allowed)
  if (length(bad) > 0L) {
    stop("upstream_processes may only include environmental_filtering, ",
         "dispersal, biotic_filtering, and evolution. Invalid: ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
  unique(x)
}

.hee_pa_process_status <- function(object, processes) {
  out <- stats::setNames(rep("data_estimated", length(processes)), processes)
  if (!is.list(object)) return(out)
  ps <- object$metadata$process_status %||% object$process_status
  if (is.null(ps)) {
    out["biotic_filtering"] <- "not_estimated"
    return(out)
  }
  if (is.data.frame(ps)) {
    pcol <- intersect(c("process", "process_slug", "name"), names(ps))[1]
    scol <- intersect(c("status", "process_status"), names(ps))[1]
    if (!is.na(pcol) && !is.na(scol)) {
      idx <- match(processes, as.character(ps[[pcol]]))
      hit <- !is.na(idx)
      out[hit] <- as.character(ps[[scol]][idx[hit]])
    }
  } else if (is.list(ps) || is.atomic(ps)) {
    vals <- unlist(ps)
    hit <- intersect(processes, names(vals))
    out[hit] <- as.character(vals[hit])
  }
  out[is.na(out) | out == ""] <- "not_estimated"
  out
}

.hee_pa_simulator <- function(object) {
  if (!is.list(object)) return(NULL)
  sim <- object$simulator %||% object$simulate
  if (is.function(sim)) sim else NULL
}

.hee_pa_exact_shapley <- function(simulator, processes, target, target_args,
                                  reference, seed, metadata) {
  n <- length(processes)
  coalitions <- .hee_pa_coalitions(processes)
  values <- lapply(seq_len(nrow(coalitions)), function(i) {
    switches <- as.logical(coalitions[i, processes, drop = TRUE])
    names(switches) <- processes
    value <- .hee_pa_simulate_value(simulator, switches, reference, seed,
                                    target, target_args, metadata)
    data.frame(coalition_id = i, coalitions[i, processes, drop = FALSE],
               value = value, stringsAsFactors = FALSE)
  })
  values <- .hee_pa_bind_rows(values)
  rows <- lapply(processes, function(p) {
    phi <- 0
    others <- setdiff(processes, p)
    for (i in seq_len(nrow(values))) {
      if (isTRUE(values[[p]][i])) next
      s_size <- sum(as.logical(values[i, others, drop = TRUE]))
      with_p <- values
      for (q in processes) {
        with_p <- with_p[with_p[[q]] == values[[q]][i], , drop = FALSE]
      }
      key <- values[i, processes, drop = FALSE]
      key[[p]] <- TRUE
      j <- rep(TRUE, nrow(values))
      for (q in processes) j <- j & values[[q]] == key[[q]]
      weight <- factorial(s_size) * factorial(n - s_size - 1) / factorial(n)
      phi <- phi + weight * (values$value[j][1] - values$value[i])
    }
    data.frame(process = p, estimate = phi, mcse = NA_real_,
               stringsAsFactors = FALSE)
  })
  shp <- .hee_pa_bind_rows(rows)
  attr(shp, "coalition_values") <- values
  shp
}

.hee_pa_mc_shapley <- function(simulator, processes, target, target_args,
                               reference, seed, n_permutations,
                               common_random_numbers, metadata) {
  set.seed(seed)
  contrib <- matrix(0, nrow = n_permutations, ncol = length(processes))
  colnames(contrib) <- processes
  for (i in seq_len(n_permutations)) {
    ord <- sample(processes)
    switches <- stats::setNames(rep(FALSE, length(processes)), processes)
    seed_i <- if (isTRUE(common_random_numbers)) seed else seed + i
    current <- .hee_pa_simulate_value(simulator, switches, reference, seed_i,
                                      target, target_args, metadata)
    for (p in ord) {
      switches[p] <- TRUE
      next_value <- .hee_pa_simulate_value(simulator, switches, reference,
                                           seed_i, target, target_args,
                                           metadata)
      contrib[i, p] <- next_value - current
      current <- next_value
    }
  }
  data.frame(process = colnames(contrib),
             estimate = colMeans(contrib),
             mcse = apply(contrib, 2, stats::sd) / sqrt(n_permutations),
             stringsAsFactors = FALSE)
}

.hee_pa_coalitions <- function(processes) {
  grid <- expand.grid(stats::setNames(rep(list(c(FALSE, TRUE)),
                                         length(processes)), processes),
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  grid
}

.hee_pa_simulate_value <- function(simulator, switches, reference, seed,
                                   target, target_args, metadata) {
  args <- list(process_switches = switches, reference = reference, seed = seed)
  ans <- tryCatch(do.call(simulator, args), error = function(e) {
    args2 <- list(draw_id = NULL, process_switches = switches,
                  reference = reference, seed = seed)
    do.call(simulator, args2)
  })
  if (is.numeric(ans) && length(ans) == 1L) return(as.numeric(ans))
  if (is.list(ans) && !is.null(ans$target_value)) {
    return(as.numeric(ans$target_value)[1])
  }
  st <- .hee_pa_standardise_state(.hee_pa_state_table(ans))
  val <- .hee_pa_target_values(st, target, target_args)
  if (nrow(val) == 0L) return(NA_real_)
  mean(val$value, na.rm = TRUE)
}

.hee_pa_target_values <- function(state, target, target_args) {
  if (!is.null(target_args$target_function) &&
      is.function(target_args$target_function)) {
    value <- target_args$target_function(state, target_args)
    return(data.frame(value = as.numeric(value)[1], stringsAsFactors = FALSE))
  }
  if (target == "occupancy") {
    state$w <- if (identical(target_args$cell_weights %||% "area", "area")) {
      state$cell_area_km2
    } else 1
    rows <- stats::aggregate(state$q * state$w,
                             state[, c("posterior_draw_id", "time_ma"),
                                   drop = FALSE],
                             sum, na.rm = TRUE)
    den <- stats::aggregate(state$w,
                            state[, c("posterior_draw_id", "time_ma"),
                                  drop = FALSE],
                            sum, na.rm = TRUE)
    rows$value <- rows$x / pmax(den$x, .Machine$double.eps)
    rows$x <- NULL
    return(rows)
  }
  if (target == "range_area") {
    rows <- stats::aggregate(state$q * state$cell_area_km2,
                             state[, c("posterior_draw_id", "time_ma"),
                                   drop = FALSE],
                             sum, na.rm = TRUE)
    names(rows)[ncol(rows)] <- "value"
    return(rows)
  }
  if (target == "lineage_richness") {
    richness <- stats::aggregate(state$q,
                                 state[, c("posterior_draw_id", "time_ma",
                                           "cell_id"),
                                       drop = FALSE],
                                 sum, na.rm = TRUE)
    names(richness)[ncol(richness)] <- "richness"
    area <- unique(state[, c("posterior_draw_id", "time_ma", "cell_id",
                             "cell_area_km2"), drop = FALSE])
    richness <- merge(richness, area,
                      by = c("posterior_draw_id", "time_ma", "cell_id"),
                      all.x = TRUE, sort = FALSE)
    rows <- stats::aggregate(richness$richness * richness$cell_area_km2,
                             richness[, c("posterior_draw_id", "time_ma"),
                                      drop = FALSE],
                             sum, na.rm = TRUE)
    den <- stats::aggregate(richness$cell_area_km2,
                            richness[, c("posterior_draw_id", "time_ma"),
                                     drop = FALSE],
                            sum, na.rm = TRUE)
    rows$value <- rows$x / pmax(den$x, .Machine$double.eps)
    rows$x <- NULL
    return(rows)
  }
  if (target == "phylogenetic_diversity") {
    rows <- stats::aggregate(state$z * state$branch_length,
                             state[, c("posterior_draw_id", "time_ma",
                                       "cell_id"),
                                   drop = FALSE],
                             sum, na.rm = TRUE)
    names(rows)[ncol(rows)] <- "pd"
    area <- unique(state[, c("posterior_draw_id", "time_ma", "cell_id",
                             "cell_area_km2"), drop = FALSE])
    rows <- merge(rows, area,
                  by = c("posterior_draw_id", "time_ma", "cell_id"),
                  all.x = TRUE, sort = FALSE)
    out <- stats::aggregate(rows$pd * rows$cell_area_km2,
                            rows[, c("posterior_draw_id", "time_ma"),
                                 drop = FALSE],
                            sum, na.rm = TRUE)
    den <- stats::aggregate(rows$cell_area_km2,
                            rows[, c("posterior_draw_id", "time_ma"),
                                 drop = FALSE],
                            sum, na.rm = TRUE)
    out$value <- out$x / pmax(den$x, .Machine$double.eps)
    out$x <- NULL
    return(out)
  }
  if (target == "functional_diversity") {
    return(.hee_pa_fd_values(state, target_args))
  }
  if (target == "turnover") {
    return(.hee_pa_turnover_values(state, target_args))
  }
  data.frame(value = NA_real_)
}

.hee_pa_fd_values <- function(state, target_args) {
  traits <- target_args$traits
  if (is.null(traits)) {
    stop("target = 'functional_diversity' requires target_args$traits ",
         "unless a custom target_function is supplied.", call. = FALSE)
  }
  tr <- as.data.frame(traits)
  tr <- .rename_first_existing(tr, "lineage_id",
                               c("lineage", "species", "taxon"))
  .require_cols(tr, "lineage_id", "traits")
  tr$lineage_id <- as.character(tr$lineage_id)
  trait_cols <- setdiff(names(tr), c("lineage_id", "time_ma"))
  trait_cols <- trait_cols[vapply(tr[trait_cols], is.numeric, logical(1))]
  if (length(trait_cols) == 0L) {
    stop("target_args$traits must contain at least one numeric trait column.",
         call. = FALSE)
  }
  for (cc in trait_cols) {
    z <- suppressWarnings(as.numeric(tr[[cc]]))
    if (isTRUE(target_args$standardise %||% TRUE)) {
      s <- stats::sd(z, na.rm = TRUE)
      if (is.finite(s) && s > 0) z <- (z - mean(z, na.rm = TRUE)) / s
    }
    tr[[cc]] <- z
  }
  state <- merge(state, tr, by = "lineage_id", all.x = TRUE, sort = FALSE)
  groups <- split(state, paste(state$posterior_draw_id, state$time_ma,
                               state$cell_id, sep = "\r"), drop = TRUE)
  vals <- lapply(groups, function(d) {
    occ <- d[d$z > 0, , drop = FALSE]
    val <- 0
    if (nrow(occ) > 1L) {
      mat <- as.matrix(occ[, trait_cols, drop = FALSE])
      keep <- stats::complete.cases(mat)
      mat <- mat[keep, , drop = FALSE]
      if (nrow(mat) > 1L) {
        dst <- as.matrix(stats::dist(mat))
        n <- nrow(dst)
        w <- rep(1 / n, n)
        val <- sum(outer(w, w) * dst)
      }
    }
    data.frame(posterior_draw_id = d$posterior_draw_id[1],
               time_ma = d$time_ma[1], cell_id = d$cell_id[1],
               cell_area_km2 = d$cell_area_km2[1], fd = val,
               stringsAsFactors = FALSE)
  })
  cell_vals <- .hee_pa_bind_rows(vals)
  out <- stats::aggregate(cell_vals$fd * cell_vals$cell_area_km2,
                          cell_vals[, c("posterior_draw_id", "time_ma"),
                                    drop = FALSE],
                          sum, na.rm = TRUE)
  den <- stats::aggregate(cell_vals$cell_area_km2,
                          cell_vals[, c("posterior_draw_id", "time_ma"),
                                    drop = FALSE],
                          sum, na.rm = TRUE)
  out$value <- out$x / pmax(den$x, .Machine$double.eps)
  out$x <- NULL
  out
}

.hee_pa_turnover_values <- function(state, target_args) {
  groups <- split(state, paste(state$posterior_draw_id, state$cell_id,
                               sep = "\r"), drop = TRUE)
  rows <- lapply(groups, function(d) {
    times <- sort(unique(d$time_ma), decreasing = TRUE)
    if (length(times) < 2L) return(NULL)
    z <- lapply(times, function(t) {
      as.character(d$lineage_id[d$time_ma == t & d$z > 0])
    })
    out <- lapply(seq_len(length(times) - 1L), function(i) {
      a <- z[[i]]
      b <- z[[i + 1L]]
      u <- union(a, b)
      inter <- intersect(a, b)
      total <- if (length(u) == 0L) 0 else 1 - length(inter) / length(u)
      gain <- length(setdiff(b, a))
      loss <- length(setdiff(a, b))
      val <- switch(target_args$component %||% "total",
                    gain = gain, loss = loss, total = total, total)
      data.frame(posterior_draw_id = d$posterior_draw_id[1],
                 time_ma = times[i], time_end_ma = times[i + 1L],
                 cell_id = d$cell_id[1],
                 cell_area_km2 = d$cell_area_km2[1],
                 value_cell = val, stringsAsFactors = FALSE)
    })
    .hee_pa_bind_rows(out)
  })
  tab <- .hee_pa_bind_rows(rows)
  if (is.null(tab)) return(data.frame(value = NA_real_))
  out <- stats::aggregate(tab$value_cell * tab$cell_area_km2,
                          tab[, c("posterior_draw_id", "time_ma"),
                              drop = FALSE],
                          sum, na.rm = TRUE)
  den <- stats::aggregate(tab$cell_area_km2,
                          tab[, c("posterior_draw_id", "time_ma"),
                              drop = FALSE],
                          sum, na.rm = TRUE)
  out$value <- out$x / pmax(den$x, .Machine$double.eps)
  out$x <- NULL
  out
}

.hee_pa_upstream_row <- function(process_slug, estimate, mcse, status,
                                 interpretation) {
  nm <- c(environmental_filtering = "Environmental filtering",
          dispersal = "Dispersal",
          biotic_filtering = "Biotic filtering",
          evolution = "Evolution")
  data.frame(
    process = unname(nm[process_slug]),
    process_slug = process_slug,
    layer = "upstream_mechanism",
    contribution_type = "shapley_counterfactual",
    estimate = estimate,
    signed_contribution = estimate,
    absolute_contribution = abs(estimate),
    mcse = mcse,
    share_within_layer = NA_real_,
    status = status,
    interpretation = interpretation,
    stringsAsFactors = FALSE
  )
}

.hee_pa_event_summary <- function(object, target, target_args, time,
                                  event_effect) {
  events <- .hee_pa_events(object)
  if (is.null(events) || nrow(events) == 0L) {
    summary <- .hee_pa_bind_rows(list(
      .hee_pa_event_row("Speciation", "speciation", NA_real_,
                        "not_estimated_no_event_log"),
      .hee_pa_event_row("Extinction", "extinction", NA_real_,
                        "not_estimated_no_event_log")
    ))
    return(list(summary = summary, events = events))
  }
  events <- .hee_pa_filter_events(events, time)
  spec <- .hee_pa_event_estimate(events, "speciation")
  ext <- .hee_pa_event_estimate(events, "lineage_extinction")
  earth <- .hee_pa_event_estimate(events, "earth_forcing")
  summary <- .hee_pa_bind_rows(list(
    .hee_pa_event_row("Speciation", "speciation", spec$estimate,
                      spec$status),
    .hee_pa_event_row("Extinction", "extinction", ext$estimate,
                      ext$status),
    .hee_pa_external_row("Earth forcing", "external_earth_forcing",
                         earth$estimate, earth$status)
  ))
  list(summary = summary, events = events)
}

.hee_pa_events <- function(object) {
  if (!is.list(object)) return(NULL)
  ev <- object$events %||% object$event_log
  if (is.null(ev)) return(NULL)
  if (is.data.frame(ev)) {
    out <- as.data.frame(ev)
  } else if (is.list(ev)) {
    rows <- list()
    if (!is.null(ev$speciation_events)) {
      d <- as.data.frame(ev$speciation_events)
      d$event_type <- d$event_type %||% "speciation"
      rows$speciation_events <- d
    }
    if (!is.null(ev$lineage_extinction_events)) {
      d <- as.data.frame(ev$lineage_extinction_events)
      d$event_type <- d$event_type %||% "lineage_extinction"
      rows$lineage_extinction_events <- d
    }
    if (!is.null(ev$earth_events)) {
      d <- as.data.frame(ev$earth_events)
      d$event_type <- d$event_type %||% "earth_forcing"
      rows$earth_events <- d
    }
    out <- .hee_pa_bind_rows(rows)
  } else {
    return(NULL)
  }
  if (is.null(out) || nrow(out) == 0L) return(NULL)
  out <- .rename_first_existing(out, "event_type", c("type", "event"))
  out <- .rename_first_existing(out, "time_ma", c("age_ma", "age"))
  out <- .rename_first_existing(out, "delta_metric",
                                c("target_delta", "delta", "contribution"))
  out$event_type <- tolower(as.character(out$event_type))
  out$event_type[grepl("spec", out$event_type)] <- "speciation"
  out$event_type[grepl("extinct|lineage.*loss", out$event_type)] <-
    "lineage_extinction"
  out$event_type[grepl("earth|land|habitat|submer|emerg", out$event_type)] <-
    "earth_forcing"
  out
}

.hee_pa_filter_events <- function(events, time) {
  if (is.null(events) || !("time_ma" %in% names(events)) || is.null(time)) {
    return(events)
  }
  t <- if (is.list(time)) c(time$from, time$to, time$at) else time
  t <- suppressWarnings(as.numeric(t))
  t <- t[is.finite(t)]
  if (length(t) == 0L) return(events)
  if (length(t) == 1L) return(events[abs(events$time_ma - t) < 1e-8,
                                    , drop = FALSE])
  tr <- range(t)
  events[events$time_ma >= tr[1] & events$time_ma <= tr[2], , drop = FALSE]
}

.hee_pa_event_estimate <- function(events, event_type) {
  d <- events[events$event_type == event_type, , drop = FALSE]
  if (nrow(d) == 0L) {
    return(list(estimate = NA_real_, status = "not_estimated_no_events"))
  }
  if ("delta_metric" %in% names(d)) {
    val <- suppressWarnings(as.numeric(d$delta_metric))
    return(list(estimate = sum(val, na.rm = TRUE),
                status = "estimated_event_delta"))
  }
  before_col <- intersect(c("metric_before", "target_before", "before"),
                          names(d))[1]
  after_col <- intersect(c("metric_after", "target_after", "after"),
                         names(d))[1]
  if (!is.na(before_col) && !is.na(after_col)) {
    val <- suppressWarnings(as.numeric(d[[after_col]]) -
                              as.numeric(d[[before_col]]))
    return(list(estimate = sum(val, na.rm = TRUE),
                status = "estimated_event_jump"))
  }
  list(estimate = NA_real_, status = "not_estimated_missing_event_effect")
}

.hee_pa_event_row <- function(process, process_slug, estimate, status) {
  data.frame(
    process = process,
    process_slug = process_slug,
    layer = "lineage_event",
    contribution_type = "event_jump",
    estimate = estimate,
    signed_contribution = estimate,
    absolute_contribution = abs(estimate),
    share_within_layer = NA_real_,
    status = status,
    interpretation = if (process_slug == "speciation") {
      "Tree/node event contribution; target dependent."
    } else {
      "Lineage extinction event contribution; local loss is a separate state-transition term."
    },
    stringsAsFactors = FALSE
  )
}

.hee_pa_external_row <- function(process, process_slug, estimate, status) {
  data.frame(
    process = process,
    process_slug = process_slug,
    layer = "external_earth",
    contribution_type = "external_event_delta",
    estimate = estimate,
    signed_contribution = estimate,
    absolute_contribution = abs(estimate),
    share_within_layer = NA_real_,
    status = status,
    interpretation = "External Earth forcing is a driver layer, not a biological process.",
    stringsAsFactors = FALSE
  )
}

.hee_pa_predictive_support <- function(object, predictive_support) {
  if (!isTRUE(predictive_support)) {
    return(list(summary = NULL, predictive_support = NULL))
  }
  ps <- if (is.list(object)) object$predictive_support %||%
    object$support else NULL
  if (is.null(ps)) {
    row <- data.frame(
      process = "Predictive support",
      process_slug = "predictive_support",
      layer = "predictive_support",
      contribution_type = "heldout_or_fossil_support",
      estimate = NA_real_,
      signed_contribution = NA_real_,
      absolute_contribution = NA_real_,
      share_within_layer = NA_real_,
      status = "not_estimated_no_predictive_support",
      interpretation = "Predictive support is reported separately from process contribution.",
      stringsAsFactors = FALSE
    )
    return(list(summary = row, predictive_support = NULL))
  }
  ps <- as.data.frame(ps)
  row <- data.frame(
    process = "Predictive support",
    process_slug = "predictive_support",
    layer = "predictive_support",
    contribution_type = "heldout_or_fossil_support",
    estimate = if ("estimate" %in% names(ps)) {
      mean(suppressWarnings(as.numeric(ps$estimate)), na.rm = TRUE)
    } else NA_real_,
    signed_contribution = NA_real_,
    absolute_contribution = NA_real_,
    share_within_layer = NA_real_,
    status = "reported_separately",
    interpretation = "Predictive support audits contribution; it is not a causal effect.",
    stringsAsFactors = FALSE
  )
  list(summary = row, predictive_support = ps)
}

.hee_pa_bind_summary <- function(...) {
  out <- .hee_pa_bind_rows(list(...))
  if (is.null(out)) return(data.frame())
  if ("layer" %in% names(out) && "absolute_contribution" %in% names(out)) {
    for (ly in unique(out$layer)) {
      idx <- out$layer == ly & is.na(out$share_within_layer)
      denom <- sum(out$absolute_contribution[out$layer == ly], na.rm = TRUE)
      if (is.finite(denom) && denom > 0) {
        out$share_within_layer[idx] <- out$absolute_contribution[idx] / denom
      }
    }
  }
  rownames(out) <- NULL
  out
}

.hee_pa_empty_summary <- function(layer, status) {
  data.frame(process = NA_character_, process_slug = NA_character_,
             layer = layer, contribution_type = NA_character_,
             estimate = NA_real_, signed_contribution = NA_real_,
             absolute_contribution = NA_real_,
             share_within_layer = NA_real_,
             status = status, interpretation = NA_character_,
             stringsAsFactors = FALSE)
}

.hee_pa_bind_rows <- function(x) {
  if (length(x) == 0L) return(NULL)
  if (!is.list(x) || is.data.frame(x)) x <- list(x)
  x <- x[!vapply(x, is.null, logical(1))]
  if (length(x) == 0L) return(NULL)
  cols <- unique(unlist(lapply(x, names), use.names = FALSE))
  out <- lapply(x, function(d) {
    d <- as.data.frame(d)
    missing <- setdiff(cols, names(d))
    for (m in missing) d[[m]] <- NA
    d[, cols, drop = FALSE]
  })
  do.call(rbind, out)
}

.hee_pa_validate_result <- function(out) {
  s <- out$summary
  if (!is.null(s) && nrow(s) > 0L &&
      any(s$status == "not_estimated" & !is.na(s$estimate))) {
    stop("Internal attribution validation failed: not_estimated rows must ",
         "not contain numeric estimates.", call. = FALSE)
  }
  invisible(TRUE)
}
