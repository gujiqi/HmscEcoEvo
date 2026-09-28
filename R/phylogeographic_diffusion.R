# Landscape-explicit phylogeographic diffusion -------------------------------

.hee_phygeo_normalise <- function(x, label) {
  x <- suppressWarnings(as.numeric(x))
  if (any(!is.finite(x) | x < 0)) {
    stop(label, " must be finite and non-negative.", call. = FALSE)
  }
  total <- sum(x)
  if (!is.finite(total) || total <= 0) {
    stop(label, " must have positive total mass.", call. = FALSE)
  }
  x / total
}

.hee_phygeo_child_key <- function(tree, node) {
  if (node <= length(tree$tip.label)) tree$tip.label[[node]] else {
    paste0("node_", node)
  }
}

#' Enforce terrestrial states after diffusion and plate carriage
#'
#' Restricts an already-composed location transition to land at its older and
#' younger endpoints. Probability assigned to marine targets is sent to an
#' explicit absorbing sink, not silently renormalised onto land. Cross-sea
#' land-to-land jumps remain possible when they are present in the input
#' dispersal transition. The sink is geography loss, not biological extinction.
#'
#' @param transition Named target-by-source column-stochastic matrix.
#' @param source_land,target_land Logical vectors named by location state (or in
#'   matrix order) for the older and younger endpoints.
#' @param sink_state Name of the added absorbing state.
#' @return A column-stochastic sparse transition with one additional sink row
#'   and column, plus mass-loss diagnostics.
#' @export
hee_dispersal_terrestrial_sink_transition <- function(
    transition, source_land, target_land,
    sink_state = "__marine_or_lost_land__") {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix is required for terrestrial sink transitions.", call. = FALSE)
  }
  P <- if (inherits(transition, "hee_dispersal_location_transition")) {
    transition$transition
  } else transition
  if (is.null(dim(P)) || nrow(P) != ncol(P) ||
      is.null(rownames(P)) || is.null(colnames(P)) ||
      !identical(rownames(P), colnames(P))) {
    stop("transition must be a square matrix with identical ordered state names.", call. = FALSE)
  }
  states <- colnames(P)
  if (sink_state %in% states) stop("sink_state duplicates a location state.", call. = FALSE)
  align_land <- function(x, label) {
    if (!is.null(names(x))) {
      if (!all(states %in% names(x))) stop(label, " lacks states.", call. = FALSE)
      x <- x[states]
    }
    if (length(x) != length(states) || anyNA(x)) {
      stop(label, " must contain one non-missing value per state.", call. = FALSE)
    }
    as.logical(x)
  }
  old <- align_land(source_land, "source_land")
  young <- align_land(target_land, "target_land")
  if (!methods::is(P, "sparseMatrix")) P <- Matrix::Matrix(P, sparse = TRUE)
  entries <- Matrix::summary(P)
  if (any(!is.finite(entries$x) | entries$x < 0) ||
      max(abs(Matrix::colSums(P) - 1)) > 1e-8) {
    stop("transition must be finite, non-negative, and source-normalised.", call. = FALSE)
  }
  keep <- old[entries$j] & young[entries$i] & entries$x > 0
  n <- length(states)
  retained <- numeric(n)
  if (any(keep)) {
    sums <- rowsum(entries$x[keep], entries$j[keep], reorder = FALSE)
    retained[as.integer(rownames(sums))] <- sums[, 1L]
  }
  lost <- pmax(0, 1 - retained)
  destination <- c(entries$i[keep], rep(n + 1L, n), n + 1L)
  source <- c(entries$j[keep], seq_len(n), n + 1L)
  weight <- c(entries$x[keep], lost, 1)
  names_all <- c(states, sink_state)
  out <- Matrix::sparseMatrix(i = destination, j = source, x = weight,
    dims = c(n + 1L, n + 1L), dimnames = list(names_all, names_all))
  error <- max(abs(Matrix::colSums(out) - 1))
  if (error > 1e-8) stop("Terrestrial sink transition lost column mass.", call. = FALSE)
  list(transition = out, sink_state = sink_state,
    diagnostics = data.frame(source_land_cells = sum(old),
      target_land_cells = sum(young),
      mean_land_source_sink_probability = if (any(old)) mean(lost[old]) else NA_real_,
      max_source_normalisation_error = error))
}

#' Add an explicit rare land-to-land jump tail to a location transition
#'
#' A low-rate jump connects occupied terrestrial source cells to a declared
#' sparse set of terrestrial landing cells. It does not create marine ancestry
#' and does not infer the jump rate from present-only community observations.
#' Use this as a documented dispersal sensitivity, not a calibrated fact.
#'
#' @param transition Named target-by-source column-stochastic transition.
#' @param land_cells Logical vector over the states at the movement stage.
#' @param hub_cells Land-state IDs available as rare-jump destinations.
#' @param jump_rate_per_myr Non-negative per-Myr jump rate.
#' @param delta_t_myr Positive stage duration.
#' @return Sparse column-stochastic transition on the original state set.
#' @export
hee_dispersal_rare_land_jump_transition <- function(
    transition, land_cells, hub_cells, jump_rate_per_myr, delta_t_myr) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix is required for rare land jumps.", call. = FALSE)
  }
  P <- if (inherits(transition, "hee_dispersal_location_transition")) {
    transition$transition
  } else transition
  if (is.null(dim(P)) || nrow(P) != ncol(P) ||
      is.null(rownames(P)) || !identical(rownames(P), colnames(P))) {
    stop("transition needs identical ordered row and column state names.", call. = FALSE)
  }
  states <- colnames(P)
  if (!is.null(names(land_cells))) {
    if (!all(states %in% names(land_cells))) stop("land_cells lacks states.", call. = FALSE)
    land_cells <- land_cells[states]
  }
  if (length(land_cells) != length(states) || anyNA(land_cells)) {
    stop("land_cells must have one non-missing value per state.", call. = FALSE)
  }
  land <- as.logical(land_cells)
  hubs <- unique(as.character(hub_cells))
  if (!length(hubs) || any(!hubs %in% states[land])) {
    stop("hub_cells must be non-empty terrestrial state IDs.", call. = FALSE)
  }
  rate <- as.numeric(jump_rate_per_myr)[1L]
  duration <- as.numeric(delta_t_myr)[1L]
  if (!is.finite(rate) || rate < 0 || !is.finite(duration) || duration <= 0) {
    stop("jump_rate_per_myr must be non-negative and delta_t_myr positive.", call. = FALSE)
  }
  if (!methods::is(P, "sparseMatrix")) P <- Matrix::Matrix(P, sparse = TRUE)
  if (any(!is.finite(P@x) | P@x < 0) ||
      max(abs(Matrix::colSums(P) - 1)) > 1e-8) {
    stop("transition must be non-negative and source-normalised.", call. = FALSE)
  }
  jump_probability <- 1 - exp(-rate * duration)
  if (jump_probability == 0) return(P)
  source <- which(land)
  landing <- match(hubs, states)
  n <- length(states)
  H <- Matrix::sparseMatrix(
    i = rep(landing, times = length(source)),
    j = rep(source, each = length(landing)),
    x = rep(jump_probability / length(landing), length(source) * length(landing)),
    dims = c(n, n), dimnames = list(states, states))
  stay <- rep(1, n); stay[source] <- 1 - jump_probability
  out <- P %*% Matrix::Diagonal(n = n, x = stay) + H
  dimnames(out) <- list(states, states)
  if (max(abs(Matrix::colSums(out) - 1)) > 1e-8) {
    stop("Rare jump transition lost probability mass.", call. = FALSE)
  }
  out
}

.hee_phygeo_as_rate <- function(rate, cell_order, edges) {
  n <- length(cell_order)
  if (is.null(rate)) {
    if (!"movement_rate_per_myr" %in% names(edges)) {
      stop(
        "movement_kernel needs movement_rate_per_myr when source_emigration_rate_per_myr is omitted.",
        call. = FALSE
      )
    }
    value <- suppressWarnings(as.numeric(edges$movement_rate_per_myr))
    if (any(!is.finite(value) | value < 0)) {
      stop("movement_rate_per_myr must be finite and non-negative.", call. = FALSE)
    }
    out <- numeric(n)
    index <- match(as.character(edges$from_cell_id), cell_order)
    out <- rowsum(value, index, reorder = FALSE)
    answer <- numeric(n)
    answer[as.integer(rownames(out))] <- out[, 1L]
    names(answer) <- cell_order
    return(answer)
  }
  rate <- suppressWarnings(as.numeric(rate))
  if (length(rate) == 1L) {
    answer <- rep(rate, n)
  } else if (!is.null(names(rate)) && all(cell_order %in% names(rate))) {
    answer <- rate[match(cell_order, names(rate))]
  } else if (length(rate) == n) {
    answer <- rate
  } else {
    stop(
      "source_emigration_rate_per_myr must be scalar, named by cell, or have one value per cell.",
      call. = FALSE
    )
  }
  if (any(!is.finite(answer) | answer < 0)) {
    stop("source_emigration_rate_per_myr must be finite and non-negative.",
         call. = FALSE)
  }
  stats::setNames(answer, cell_order)
}

#' Construct a Markov location transition from a landscape diffusion kernel
#'
#' This is the location-distribution counterpart of
#' [hee_dispersal_spherical_kernel()].  It converts a directional,
#' landscape-weighted spherical diffusion kernel into one *column-stochastic*
#' transition matrix.  The matrix transports a lineage-location probability
#' distribution; it is **not** a cell occupancy update and has no
#' colonisation, establishment, persistence, local-extinction, or HMSC
#' environmental-response argument.
#'
#' For source cell \eqn{i}, the probability of making a movement during
#' \eqn{\Delta t} is \eqn{1-\exp(-m_i\Delta t)}.  Conditional on movement,
#' destination probabilities are the Arias-style directional kernel
#' \eqn{K_{i\to j}}.  The remaining probability stays in the source state.
#' Dynamic land and topography therefore enter once through the documented
#' movement landscape used to construct the kernel.  HMSC-derived
#' environmental support is a separate P1 result and must not be multiplied
#' into this transition.
#'
#' Supplying `available_cells`, or separate `source_available_cells` and
#' `target_available_cells`, adds an absorbing `sink_state`. The separate
#' inputs distinguish a source that is available at the older stage from a
#' target that is available at the younger stage: a land cell can still export
#' movement before it is submerged, but its retained mass enters the sink when
#' that carrier is unavailable in the target stage.
#' This is a transparent state-space representation of unavailable geography,
#' not an inferred biological persistence or lineage-extinction process.
#'
#' @references Arias, J. S. (2024). Phylogenetic biogeography inference using
#'   dynamic paleogeography models and explicit geographic ranges. *Systematic
#'   Biology*, 73, 995-1014. https://doi.org/10.1093/sysbio/syae051.
#'
#' @param movement_kernel A [hee_dispersal_spherical_kernel()] object or an
#'   edge table with `from_cell_id`, `to_cell_id`, and
#'   `directional_weight`.
#' @param cell_order Stable state IDs. Defaults to the IDs present in the edge
#'   table. Supply it when a slice has no edges.
#' @param delta_t_myr Positive duration represented by the transition. It
#'   defaults to the kernel metadata when available.
#' @param source_emigration_rate_per_myr Optional scalar, cell-named vector, or
#'   cell vector. If omitted, the per-source sum of
#'   `movement_rate_per_myr` is used.
#' @param available_cells Optional logical/numeric vector or named vector over
#'   `cell_order`, defining land/habitat states at both ends of a static
#'   transition. Use the separate source/target arguments for a changing
#'   palaeogeographic stage.
#' @param source_available_cells,target_available_cells Optional logical/numeric
#'   vectors over `cell_order` defining geographic availability at the older
#'   and younger ends of a transition. They cannot be supplied with
#'   `available_cells`.
#' @param sink_state Name for the absorbing unavailable-geography state. It is
#'   added only when `available_cells` is supplied.
#' @param sparse Return a sparse matrix when **Matrix** is installed.
#' @return A `hee_dispersal_location_transition` list with target-by-source
#'   `transition`, source movement probabilities, and metadata.
#' @export
#'
#' @examples
#' edges <- data.frame(
#'   from_cell_id = c("a", "a", "b"),
#'   to_cell_id = c("b", "c", "a"),
#'   directional_weight = c(0.75, 0.25, 1),
#'   movement_rate_per_myr = c(0.4, 0.4, 0.2)
#' )
#' z <- hee_dispersal_location_transition(edges, c("a", "b", "c"),
#'   delta_t_myr = 1
#' )
#' colSums(as.matrix(z$transition))
hee_dispersal_location_transition <- function(
    movement_kernel,
    cell_order = NULL,
    delta_t_myr = NULL,
    source_emigration_rate_per_myr = NULL,
    available_cells = NULL,
    source_available_cells = NULL,
    target_available_cells = NULL,
    sink_state = "__outside_available_geography__",
    sparse = TRUE) {
  is_kernel <- inherits(movement_kernel, "hee_dispersal_spherical_kernel")
  is_rate_matrix <- !is_kernel && !is.data.frame(movement_kernel) &&
    !is.list(movement_kernel) && !is.null(dim(movement_kernel))
  edges <- if (is_kernel) {
    movement_kernel$edges
  } else if (is_rate_matrix) {
    if (nrow(movement_kernel) != ncol(movement_kernel) ||
        is.null(rownames(movement_kernel)) || is.null(colnames(movement_kernel)) ||
        !setequal(rownames(movement_kernel), colnames(movement_kernel))) {
      stop("A movement-rate matrix must be square with identical named source and target states.",
           call. = FALSE)
    }
    entry <- if (methods::is(movement_kernel, "Matrix")) {
      Matrix::summary(movement_kernel)
    } else {
      ij <- which(movement_kernel != 0, arr.ind = TRUE)
      data.frame(i = ij[, 1L], j = ij[, 2L], x = movement_kernel[ij])
    }
    value <- suppressWarnings(as.numeric(entry$x))
    if (any(!is.finite(value) | value < 0)) {
      stop("movement-rate matrix values must be finite and non-negative.",
           call. = FALSE)
    }
    source_rate <- rowsum(value, entry$j, reorder = FALSE)
    denom <- numeric(ncol(movement_kernel))
    denom[as.integer(rownames(source_rate))] <- source_rate[, 1L]
    data.frame(
      from_cell_id = colnames(movement_kernel)[entry$j],
      to_cell_id = rownames(movement_kernel)[entry$i],
      directional_weight = value / denom[entry$j],
      movement_rate_per_myr = value,
      stringsAsFactors = FALSE
    )
  } else {
    as.data.frame(movement_kernel, stringsAsFactors = FALSE)
  }
  .require_cols(edges, c("from_cell_id", "to_cell_id", "directional_weight"),
                "movement_kernel")
  edges$from_cell_id <- as.character(edges$from_cell_id)
  edges$to_cell_id <- as.character(edges$to_cell_id)
  direction <- suppressWarnings(as.numeric(edges$directional_weight))
  if (any(!is.finite(direction) | direction < 0)) {
    stop("directional_weight must be finite and non-negative.", call. = FALSE)
  }
  if (is.null(cell_order) && is_rate_matrix) {
    cell_order <- colnames(movement_kernel)
  } else if (is.null(cell_order)) {
    cell_order <- unique(c(edges$from_cell_id, edges$to_cell_id))
  }
  cell_order <- as.character(cell_order)
  if (!length(cell_order) || any(!nzchar(cell_order)) || anyDuplicated(cell_order)) {
    stop("cell_order must contain unique, non-blank cell IDs.", call. = FALSE)
  }
  if (any(!edges$from_cell_id %in% cell_order) ||
      any(!edges$to_cell_id %in% cell_order)) {
    stop("movement_kernel includes cells absent from cell_order.", call. = FALSE)
  }
  if (is.null(delta_t_myr) && is_kernel) {
    delta_t_myr <- movement_kernel$metadata$delta_t_myr
  }
  delta_t_myr <- suppressWarnings(as.numeric(delta_t_myr)[1L])
  if (!is.finite(delta_t_myr) || delta_t_myr <= 0) {
    stop("delta_t_myr must be a positive finite duration.", call. = FALSE)
  }
  rate <- .hee_phygeo_as_rate(source_emigration_rate_per_myr, cell_order, edges)
  source_index <- match(edges$from_cell_id, cell_order)
  target_index <- match(edges$to_cell_id, cell_order)

  if (!is.null(available_cells) &&
      (!is.null(source_available_cells) || !is.null(target_available_cells))) {
    stop("Supply available_cells or separate source_available_cells/target_available_cells, not both.",
         call. = FALSE)
  }
  align_availability <- function(x, label) {
    if (is.null(x)) return(rep(TRUE, length(cell_order)))
    if (!is.null(names(x)) && all(cell_order %in% names(x))) {
      x <- x[match(cell_order, names(x))]
    } else if (length(x) != length(cell_order)) {
      stop(label, " must be named by cell_order or have one value per cell.",
           call. = FALSE)
    }
    x <- suppressWarnings(as.numeric(x))
    is.finite(x) & x > 0
  }
  source_available <- align_availability(
    if (is.null(available_cells)) source_available_cells else available_cells,
    "source_available_cells"
  )
  target_available <- align_availability(
    if (is.null(available_cells)) target_available_cells else available_cells,
    "target_available_cells"
  )

  # A destination unavailable in the target stage cannot receive active
  # diffusion. Renormalise only over admissible destinations; the rate remains
  # a source property, preserving the Arias conditional movement semantics.
  keep_edge <- source_available[source_index] & target_available[target_index] &
    direction > 0
  source_total <- rowsum(direction[keep_edge], source_index[keep_edge],
                          reorder = FALSE)
  normalised_direction <- numeric(nrow(edges))
  if (nrow(source_total)) {
    denom <- numeric(length(cell_order))
    denom[as.integer(rownames(source_total))] <- source_total[, 1L]
    normalised_direction[keep_edge] <- direction[keep_edge] /
      denom[source_index[keep_edge]]
  }
  move_probability <- 1 - exp(-rate * delta_t_myr)
  move_probability[!source_available] <- 0
  no_destination <- source_available & !(seq_along(cell_order) %in%
    as.integer(rownames(source_total)))
  move_probability[no_destination] <- 0

  values <- normalised_direction * move_probability[source_index]
  n_state <- length(cell_order)
  use_sink <- !is.null(available_cells) || !is.null(source_available_cells) ||
    !is.null(target_available_cells)
  states <- if (use_sink) c(cell_order, sink_state) else cell_order
  if (use_sink && sink_state %in% cell_order) {
    stop("sink_state must not duplicate a cell_order ID.", call. = FALSE)
  }
  i <- target_index[values > 0]
  j <- source_index[values > 0]
  x <- values[values > 0]
  # A source retains the complementary probability only if the same carrier
  # remains geographically available at the target stage. Any unmatched mass
  # enters a sink below, rather than being mislabelled as demographic loss.
  diagonal_target <- seq_len(n_state)
  diagonal_value <- ifelse(source_available & target_available,
                           1 - move_probability, 0)
  i <- c(i, diagonal_target[diagonal_value > 0])
  j <- c(j, diagonal_target[diagonal_value > 0])
  x <- c(x, diagonal_value[diagonal_value > 0])
  if (use_sink) {
    # Complete every source column with its geographic-loss probability.
    # This handles submergence after a source has exported part of its mass.
    retained_by_source <- numeric(n_state)
    if (length(x)) retained_by_source <- rowsum(x, j, reorder = FALSE)
    retained <- numeric(n_state)
    if (length(retained_by_source)) {
      retained[as.integer(rownames(retained_by_source))] <- retained_by_source[, 1L]
    }
    sink_weight <- pmax(0, 1 - retained)
    enter_sink <- which(sink_weight > 0)
    if (length(enter_sink)) {
      i <- c(i, rep(n_state + 1L, length(enter_sink)))
      j <- c(j, enter_sink)
      x <- c(x, sink_weight[enter_sink])
    }
    i <- c(i, n_state + 1L)
    j <- c(j, n_state + 1L)
    x <- c(x, 1)
  }
  dimensions <- c(length(states), length(states))
  if (isTRUE(sparse) && requireNamespace("Matrix", quietly = TRUE)) {
    transition <- Matrix::sparseMatrix(
      i = i, j = j, x = x, dims = dimensions,
      dimnames = list(states, states)
    )
  } else {
    transition <- matrix(0, nrow = dimensions[[1L]], ncol = dimensions[[2L]],
                         dimnames = list(states, states))
    for (edge_index in seq_along(x)) {
      transition[i[[edge_index]], j[[edge_index]]] <-
        transition[i[[edge_index]], j[[edge_index]]] + x[[edge_index]]
    }
  }
  if (is.null(dim(transition)) || length(dim(transition)) != 2L) {
    stop("Internal error: location transition did not produce a two-dimensional matrix (class: ",
         paste(class(transition), collapse = ", "), ").", call. = FALSE)
  }
  transition_colsum <- if (methods::is(transition, "Matrix")) {
    Matrix::colSums(transition)
  } else base::colSums(transition)
  mass_error <- max(abs(transition_colsum - 1))
  if (!is.finite(mass_error) || mass_error > 1e-10) {
    stop("Internal error: location transition is not source-normalized.",
         call. = FALSE)
  }
  structure(
    list(
      transition = transition,
      source_move_probability = stats::setNames(move_probability, cell_order),
      source_emigration_rate_per_myr = rate,
      source_available = stats::setNames(source_available, cell_order),
      target_available = stats::setNames(target_available, cell_order),
      sink_state = if (use_sink) sink_state else NA_character_,
      diagnostics = data.frame(
        n_cell_states = n_state,
        n_total_states = length(states),
        n_source_available_cells = sum(source_available),
        n_target_available_cells = sum(target_available),
        n_active_edges = sum(values > 0),
        max_source_normalisation_error = mass_error,
        stringsAsFactors = FALSE
      ),
      interpretation = paste(
        "Arias-style location transition: directional spherical diffusion plus",
        "source retention. It has no colonisation or persistence component."
      )
    ),
    class = "hee_dispersal_location_transition"
  )
}

#' Propagate lineage location distributions through a diffusion transition
#'
#' Applies a source-normalized [hee_dispersal_location_transition()] matrix to
#' one or more lineage location distributions.  Columns are probability mass
#' over explicit geographic states, not occupancy probabilities.  Thus a
#' column sums to one before conditioning on a supplied geographic sink.
#'
#' @param location_mass Numeric state-by-lineage matrix with non-negative
#'   columns of positive mass.
#' @param transition A `hee_dispersal_location_transition` object or a
#'   target-by-source column-stochastic matrix.
#' @param condition_on_non_sink If `TRUE`, remove and renormalize an optional
#'   sink state. This is only appropriate for an explicitly conditioned
#'   surviving-lineage scenario and is reported in the return value.
#' @return A list with propagated `location_mass`, `sink_mass`, and transition
#'   metadata.
#' @export
#'
#' @examples
#' transition <- hee_dispersal_location_transition(
#'   data.frame(from_cell_id = "a", to_cell_id = "b",
#'              directional_weight = 1, movement_rate_per_myr = 1),
#'   cell_order = c("a", "b"), delta_t_myr = 1
#' )
#' hee_dispersal_location_distribution_step(
#'   matrix(c(1, 0), ncol = 1,
#'          dimnames = list(c("a", "b"), "lineage_1")), transition
#' )
hee_dispersal_location_distribution_step <- function(location_mass,
                                                      transition,
                                                      condition_on_non_sink = FALSE) {
  P <- if (inherits(transition, "hee_dispersal_location_transition")) {
    transition$transition
  } else transition
  x <- as.matrix(location_mass)
  storage.mode(x) <- "double"
  if (is.null(rownames(x)) || anyDuplicated(rownames(x))) {
    stop("location_mass must have unique named state rows.", call. = FALSE)
  }
  if (is.null(rownames(P)) || is.null(colnames(P)) ||
      !setequal(rownames(P), colnames(P)) ||
      !setequal(rownames(P), rownames(x))) {
    stop("transition and location_mass must use exactly the same named states.",
         call. = FALSE)
  }
  states <- rownames(P)
  P <- P[states, states, drop = FALSE]
  x <- x[states, , drop = FALSE]
  if (any(!is.finite(x) | x < 0) || any(colSums(x) <= 0)) {
    stop("location_mass must be finite, non-negative, and positive by lineage.",
         call. = FALSE)
  }
  P_colsum <- if (methods::is(P, "Matrix")) Matrix::colSums(P) else base::colSums(P)
  if (any(abs(P_colsum - 1) > 1e-8)) {
    stop("transition must sum to one by source state.", call. = FALSE)
  }
  x <- sweep(x, 2L, colSums(x), "/")
  next_mass <- as.matrix(P %*% x)
  sink_state <- if (inherits(transition, "hee_dispersal_location_transition")) {
    transition$sink_state
  } else NA_character_
  sink_mass <- if (!is.na(sink_state) && sink_state %in% rownames(next_mass)) {
    next_mass[sink_state, , drop = FALSE]
  } else matrix(0, nrow = 1L, ncol = ncol(next_mass),
                dimnames = list("no_sink", colnames(next_mass)))
  if (isTRUE(condition_on_non_sink) && !is.na(sink_state) &&
      sink_state %in% rownames(next_mass)) {
    keep <- setdiff(rownames(next_mass), sink_state)
    retained <- next_mass[keep, , drop = FALSE]
    total <- colSums(retained)
    if (any(total <= 0)) {
      stop("Cannot condition on non-sink geography when all mass entered the sink.",
           call. = FALSE)
    }
    next_mass <- sweep(retained, 2L, total, "/")
  }
  list(
    location_mass = next_mass,
    sink_mass = sink_mass,
    conditioned_on_non_sink = isTRUE(condition_on_non_sink),
    interpretation = paste(
      "Location-distribution propagation only. Values are not colonisation,",
      "persistence, occupancy, or extinction probabilities."
    )
  )
}

.hee_phygeo_steps <- function(branch_steps, key, states) {
  item <- branch_steps[[key]]
  if (is.null(item)) item <- branch_steps[[sub("^node_", "", key)]]
  if (is.null(item)) {
    stop("branch_steps is missing child branch: ", key, call. = FALSE)
  }
  if (is.matrix(item) || methods::is(item, "Matrix")) item <- list(item)
  if (!is.list(item) || !length(item)) {
    stop("Each branch_steps entry must be a transition or non-empty list of transitions.",
         call. = FALSE)
  }
  lapply(seq_along(item), function(i) {
    entry <- item[[i]]
    P <- if (is.list(entry) && !is.matrix(entry) && !methods::is(entry, "Matrix")) {
      entry$transition
    } else entry
    time_ma <- if (is.list(entry) && !is.null(entry$time_ma)) {
      suppressWarnings(as.numeric(entry$time_ma)[1L])
    } else NA_real_
    if (is.null(P) || is.null(dim(P)) || !identical(dim(P), c(length(states), length(states))) ||
        is.null(rownames(P)) || is.null(colnames(P)) ||
        !setequal(rownames(P), states) || !setequal(colnames(P), states)) {
      stop("Each time-ordered transition must be a square named matrix over the shared state set.",
           call. = FALSE)
    }
    P <- P[states, states, drop = FALSE]
    value <- if (methods::is(P, "Matrix")) P@x else as.numeric(P)
    P_colsum <- if (methods::is(P, "Matrix")) Matrix::colSums(P) else base::colSums(P)
    if (any(!is.finite(value) | value < 0) ||
        any(abs(P_colsum - 1) > 1e-8)) {
      stop("Each transition must be non-negative and source-normalized.",
           call. = FALSE)
    }
    list(transition = P, time_ma = time_ma, stage_index = i)
  })
}

#' Compose diffusion with a plate-carriage transition
#'
#' Combines a horizontal explicit-grid location transition with a separate,
#' source-normalized plate-carriage transition.  If `D` is the diffusion
#' transition and `A` maps a location after horizontal movement at the older
#' stage to its plate-carried position at the younger stage, the result is
#' \eqn{A D}.  Both matrices are column stochastic, so their composition is
#' also column stochastic.  This operation transports *location probability
#' mass*, not occupancy, colonisation, persistence, or extinction.
#'
#' Keeping the two operators separate makes it possible to verify that
#' continental motion has not been misinterpreted as active biological
#' dispersal.  It implements the separation of plate rotation and landscape
#' diffusion used by Arias (2024) and landscape-explicit time-ordered
#' phylogeographic analyses.
#'
#' @references Arias, J. S. (2024). Phylogenetic biogeography inference using
#'   dynamic paleogeography models and explicit geographic ranges. *Systematic
#'   Biology*, 73, 995-1014. https://doi.org/10.1093/sysbio/syae051.
#'
#' @param diffusion A [hee_dispersal_location_transition()] object or a named,
#'   target-by-source column-stochastic transition matrix.
#' @param plate_transport A named, target-by-source column-stochastic matrix
#'   over exactly the same location states. It can be a probabilistic aggregate
#'   when a coarse display cell contains more than one plate carrier.
#' @return A `hee_dispersal_location_transition` object whose `transition` is
#'   the composed plate-carried diffusion transition.
#' @export
#'
#' @examples
#' D <- matrix(c(0.8, 0.2, 0.2, 0.8), 2, 2,
#'             dimnames = list(c("a", "b"), c("a", "b")))
#' A <- matrix(c(0, 1, 1, 0), 2, 2,
#'             dimnames = list(c("a", "b"), c("a", "b")))
#' hee_dispersal_compose_plate_transport(D, A)$transition
hee_dispersal_compose_plate_transport <- function(diffusion, plate_transport) {
  D <- if (inherits(diffusion, "hee_dispersal_location_transition")) {
    diffusion$transition
  } else diffusion
  A <- plate_transport
  validate_transition <- function(P, label) {
    if (is.null(P) || is.null(dim(P)) || length(dim(P)) != 2L ||
        nrow(P) != ncol(P) || is.null(rownames(P)) || is.null(colnames(P)) ||
        !setequal(rownames(P), colnames(P))) {
      stop(label, " must be a square named target-by-source matrix.", call. = FALSE)
    }
    values <- if (methods::is(P, "Matrix")) P@x else as.numeric(P)
    sums <- if (methods::is(P, "Matrix")) Matrix::colSums(P) else base::colSums(P)
    if (any(!is.finite(values) | values < 0) || any(abs(sums - 1) > 1e-8)) {
      stop(label, " must be non-negative and source-normalized.", call. = FALSE)
    }
    P
  }
  D <- validate_transition(D, "diffusion")
  A <- validate_transition(A, "plate_transport")
  states <- rownames(D)
  if (!setequal(rownames(A), states)) {
    stop("plate_transport must use exactly the diffusion state set.", call. = FALSE)
  }
  D <- D[states, states, drop = FALSE]
  A <- A[states, states, drop = FALSE]
  P <- A %*% D
  sums <- if (methods::is(P, "Matrix")) Matrix::colSums(P) else base::colSums(P)
  error <- max(abs(sums - 1))
  if (!is.finite(error) || error > 1e-10) {
    stop("Internal error: composed transition is not source-normalized.", call. = FALSE)
  }
  structure(
    list(
      transition = P,
      plate_transport = A,
      horizontal_diffusion = D,
      sink_state = if (inherits(diffusion, "hee_dispersal_location_transition")) {
        diffusion$sink_state
      } else NA_character_,
      diagnostics = data.frame(
        n_states = length(states),
        max_source_normalisation_error = error,
        stringsAsFactors = FALSE
      ),
      interpretation = paste(
        "Plate-carried explicit location transition: horizontal diffusion is",
        "composed with a separately documented geographic transport matrix;",
        "it contains no demographic colonisation or persistence process."
      )
    ),
    class = "hee_dispersal_location_transition"
  )
}

#' Time-ordered phylogeographic pruning on an explicit dynamic grid
#'
#' Computes a Felsenstein pruning likelihood when each branch is split into a
#' time-ordered sequence of landscape-explicit location transitions.  This is
#' the dynamic-grid form needed when plate-carried geography and a movement
#' landscape change through a branch.  It integrates ancestral locations
#' conditional on the terminal geographic observations; it does not introduce
#' colonisation, establishment, persistence, local extinction, or a second
#' environmental-suitability multiplier.
#'
#' A `state_callback`, when supplied, is invoked with each branch's posterior
#' location density at every supplied stage endpoint.  This permits a runner
#' to stream map summaries to disk rather than retaining a
#' lineage-by-cell-by-time cube in memory.
#'
#' @references Arias, J. S. (2024). Phylogenetic biogeography inference using
#'   dynamic paleogeography models and explicit geographic ranges. *Systematic
#'   Biology*, 73, 995-1014. https://doi.org/10.1093/sysbio/syae051;
#'   Flannery-Sutherland, J. T. et al. (2025). Landscape-explicit
#'   phylogeography illuminates the ecographic radiation of early
#'   archosauromorph reptiles. *Nature Ecology & Evolution*, 9, 1138-1152.
#'   https://doi.org/10.1038/s41559-025-02739-y
#'
#' @param tree An `ape::phylo`-compatible rooted binary tree.
#' @param tip_likelihood Named list of non-negative named vectors over the
#'   common grid states. Each tip's vector is an explicit terminal-location or
#'   terminal-range observation density.
#' @param branch_steps Named list. Each child branch is named by its tip label
#'   or `"node_<number>"` and contains an older-to-younger list of transition
#'   matrices or `list(transition = ..., time_ma = ...)` records.
#' @param root_prior Optional named non-negative location prior.
#' @param state_callback Optional function with arguments `branch_key`,
#'   `child_node`, `time_ma`, `stage_index`, and `posterior`.
#' @return A `hee_dispersal_time_ordered_pruning` object with log likelihood,
#'   root and node location posteriors, and branch messages.
#' @export
#'
#' @examples
#' tree <- structure(list(
#'   edge = matrix(c(3, 1, 3, 2), ncol = 2, byrow = TRUE),
#'   tip.label = c("a", "b"), Nnode = 1L
#' ), class = "phylo")
#' K <- matrix(0.5, 2, 2, dimnames = list(c("x", "y"), c("x", "y")))
#' hee_dispersal_time_ordered_pruning(
#'   tree, list(a = c(x = 1, y = 0), b = c(x = 0, y = 1)),
#'   list(a = list(list(transition = K, time_ma = 0)),
#'        b = list(list(transition = K, time_ma = 0)))
#' )
hee_dispersal_time_ordered_pruning <- function(tree,
                                                tip_likelihood,
                                                branch_steps,
                                                root_prior = NULL,
                                                state_callback = NULL) {
  if (is.null(tree$edge) || is.null(tree$tip.label)) {
    stop("tree must provide edge and tip.label fields.", call. = FALSE)
  }
  edge <- as.matrix(tree$edge)
  if (ncol(edge) != 2L || nrow(edge) < 1L) {
    stop("tree$edge must be a non-empty two-column matrix.", call. = FALSE)
  }
  n_tip <- length(tree$tip.label)
  root <- setdiff(edge[, 1L], edge[, 2L])
  if (length(root) != 1L) stop("tree must have one root.", call. = FALSE)
  if (!is.list(tip_likelihood) || is.null(names(tip_likelihood)) ||
      !all(tree$tip.label %in% names(tip_likelihood))) {
    stop("tip_likelihood must be a named list covering every tree tip.",
         call. = FALSE)
  }
  states <- names(tip_likelihood[[tree$tip.label[[1L]]]])
  if (is.null(states) || !length(states) || any(!nzchar(states)) || anyDuplicated(states)) {
    stop("tip likelihoods must have unique named state values.", call. = FALSE)
  }
  normalize_named <- function(x, label) {
    if (is.null(names(x)) || !setequal(names(x), states)) {
      stop(label, " must use exactly the shared state set.", call. = FALSE)
    }
    stats::setNames(.hee_phygeo_normalise(x[states], label), states)
  }
  tips <- lapply(tree$tip.label, function(label) {
    normalize_named(tip_likelihood[[label]], paste0("tip_likelihood$", label))
  })
  names(tips) <- tree$tip.label
  if (is.null(root_prior)) root_prior <- stats::setNames(rep(1, length(states)), states)
  root_prior <- normalize_named(root_prior, "root_prior")
  if (!is.list(branch_steps) || is.null(names(branch_steps))) {
    stop("branch_steps must be a named list.", call. = FALSE)
  }
  children <- split(edge[, 2L], edge[, 1L])
  all_nodes <- sort(unique(c(edge[, 1L], edge[, 2L])))
  max_node <- max(all_nodes)
  inside <- vector("list", max_node)
  node_log_scale <- rep(NA_real_, max_node)
  messages <- list()
  parsed_steps <- list()
  steps_for <- function(child) {
    key <- .hee_phygeo_child_key(tree, child)
    if (is.null(parsed_steps[[key]])) {
      parsed_steps[[key]] <<- .hee_phygeo_steps(branch_steps, key, states)
    }
    parsed_steps[[key]]
  }
  backward_message <- function(child, value) {
    steps <- steps_for(child)
    scale <- 0
    out <- value
    for (j in rev(seq_along(steps))) {
      out <- as.vector(crossprod(steps[[j]]$transition, out))
      total <- sum(out)
      if (!is.finite(total) || total <= 0) {
        return(list(value = rep(0, length(states)), log_scale = -Inf))
      }
      out <- out / total
      scale <- scale + log(total)
    }
    list(value = stats::setNames(out, states), log_scale = scale)
  }
  visit <- function(node) {
    if (!is.null(inside[[node]])) return(invisible(NULL))
    if (node <= n_tip) {
      inside[[node]] <<- tips[[tree$tip.label[[node]]]]
      node_log_scale[[node]] <<- 0
      return(invisible(NULL))
    }
    child_nodes <- children[[as.character(node)]]
    likelihood <- rep(1, length(states)); total_scale <- 0
    for (child in child_nodes) {
      visit(child)
      message <- backward_message(child, inside[[child]])
      if (!is.finite(message$log_scale) || !is.finite(node_log_scale[[child]])) {
        inside[[node]] <<- rep(0, length(states))
        node_log_scale[[node]] <<- -Inf
        return(invisible(NULL))
      }
      key <- .hee_phygeo_child_key(tree, child)
      messages[[key]] <<- message$value
      likelihood <- likelihood * message$value
      total_scale <- total_scale + node_log_scale[[child]] + message$log_scale
    }
    scale <- sum(likelihood)
    if (!is.finite(scale) || scale <= 0) {
      inside[[node]] <<- rep(0, length(states))
      node_log_scale[[node]] <<- -Inf
    } else {
      inside[[node]] <<- stats::setNames(likelihood / scale, states)
      node_log_scale[[node]] <<- total_scale + log(scale)
    }
    invisible(NULL)
  }
  visit(root)
  root_term <- sum(root_prior * inside[[root]])
  log_likelihood <- if (!is.finite(node_log_scale[[root]]) || root_term <= 0) {
    -Inf
  } else node_log_scale[[root]] + log(root_term)
  root_posterior <- if (is.finite(log_likelihood)) {
    stats::setNames(root_prior * inside[[root]] / root_term, states)
  } else stats::setNames(rep(NA_real_, length(states)), states)

  outside <- vector("list", max_node)
  node_posterior <- vector("list", max_node)
  emit_branch_states <- function(child, context) {
    if (is.null(state_callback)) return(invisible(NULL))
    steps <- steps_for(child)
    backward <- vector("list", length(steps) + 1L)
    backward[[length(steps) + 1L]] <- inside[[child]]
    for (j in rev(seq_along(steps))) {
      z <- as.vector(crossprod(steps[[j]]$transition, backward[[j + 1L]]))
      backward[[j]] <- stats::setNames(.hee_phygeo_normalise(z, "branch backward message"), states)
    }
    forward <- context
    key <- .hee_phygeo_child_key(tree, child)
    for (j in seq_along(steps)) {
      forward <- as.vector(steps[[j]]$transition %*% forward)
      forward <- stats::setNames(.hee_phygeo_normalise(forward, "branch forward message"), states)
      posterior <- forward * backward[[j + 1L]]
      posterior <- stats::setNames(.hee_phygeo_normalise(posterior, "branch posterior"), states)
      state_callback(
        branch_key = key, child_node = child, time_ma = steps[[j]]$time_ma,
        stage_index = steps[[j]]$stage_index, posterior = posterior
      )
    }
    invisible(NULL)
  }
  forward_visit <- function(node) {
    node_posterior[[node]] <<- stats::setNames(
      .hee_phygeo_normalise(outside[[node]] * inside[[node]], "node posterior"),
      states
    )
    if (node <= n_tip) return(invisible(NULL))
    child_nodes <- children[[as.character(node)]]
    for (child in child_nodes) {
      context <- outside[[node]]
      for (sibling in setdiff(child_nodes, child)) {
        context <- context * messages[[.hee_phygeo_child_key(tree, sibling)]]
      }
      context <- stats::setNames(.hee_phygeo_normalise(context, "outside branch message"), states)
      emit_branch_states(child, context)
      forward <- context
      for (step in steps_for(child)) {
        forward <- as.vector(step$transition %*% forward)
        forward <- stats::setNames(.hee_phygeo_normalise(forward, "forward branch message"), states)
      }
      outside[[child]] <<- forward
      forward_visit(child)
    }
    invisible(NULL)
  }
  if (is.finite(log_likelihood)) {
    outside[[root]] <- root_prior
    forward_visit(root)
  }
  node_names <- vapply(all_nodes, function(node) .hee_phygeo_child_key(tree, node), character(1))
  structure(
    list(
      log_likelihood = log_likelihood,
      root_node = root,
      root_posterior = root_posterior,
      node_posterior = stats::setNames(node_posterior[all_nodes], node_names),
      node_inside = stats::setNames(inside[all_nodes], node_names),
      node_log_scale = stats::setNames(node_log_scale[all_nodes], node_names),
      branch_messages = messages,
      states = states,
      interpretation = paste(
        "Dynamic explicit-grid phylogeographic reconstruction conditional on",
        "terminal locations/ranges. It is not a demographic occupancy model;",
        "colonisation and persistence are deliberately absent."
      )
    ),
    class = "hee_dispersal_time_ordered_pruning"
  )
}
