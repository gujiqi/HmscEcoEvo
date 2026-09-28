#' Validate an area-weighted palaeogeographic cell-carriage table
#'
#' A dynamic palaeo-grid normally changes coordinates and land extent between
#' time slices. This helper validates the explicit transport table used to
#' carry occupancy from a source slice to a younger target slice before any
#' active biological dispersal is calculated. It is deliberately separate
#' from movement: continental drift is not organismal dispersal.
#'
#' The table represents target-cell mixture weights. For a target cell `i`,
#' weights must sum to at most one and its carried occupancy is
#' \deqn{q_{i,t+1}^{carry}=\sum_s w_{is,t}q_{s,t}.}
#' Weights should normally be derived from overlap area or a documented
#' plate-track approximation and area-normalised within each target cell.
#'
#' @param transport A data frame containing source and target cell IDs and a
#'   finite non-negative target-normalised weight.
#' @param source_cells,target_cells Optional source and target cell IDs used to
#'   check coverage and unknown IDs.
#' @param source_col,target_col,weight_col Column names.
#' @param tolerance Numerical tolerance for target-row weights.
#' @return A `hee_plate_grid_transport` list with the validated table and
#'   coverage diagnostics.
#' @export
hee_plate_grid_transport <- function(transport,
                                     source_cells = NULL,
                                     target_cells = NULL,
                                     source_col = "source_cell_id",
                                     target_col = "target_cell_id",
                                     weight_col = "target_weight",
                                     tolerance = 1e-8) {
  x <- as.data.frame(transport, stringsAsFactors = FALSE)
  .require_cols(x, c(source_col, target_col, weight_col), "transport")
  x[[source_col]] <- as.character(x[[source_col]])
  x[[target_col]] <- as.character(x[[target_col]])
  x[[weight_col]] <- suppressWarnings(as.numeric(x[[weight_col]]))
  if (any(!nzchar(x[[source_col]]) | !nzchar(x[[target_col]]))) {
    stop("transport contains blank source or target cell IDs.", call. = FALSE)
  }
  if (any(!is.finite(x[[weight_col]]) | x[[weight_col]] < 0)) {
    stop("transport weights must be finite and non-negative.", call. = FALSE)
  }
  key <- paste(x[[source_col]], x[[target_col]], sep = "\r")
  if (anyDuplicated(key)) {
    stop("transport must contain at most one row per source-target cell pair.",
         call. = FALSE)
  }
  target_sum <- stats::aggregate(
    x[[weight_col]], list(target_cell_id = x[[target_col]]), sum
  )
  names(target_sum)[[2L]] <- "target_weight_sum"
  if (any(target_sum$target_weight_sum > 1 + tolerance)) {
    stop("transport target weights exceed one; weights must be target-normalised.",
         call. = FALSE)
  }
  source_cells <- if (is.null(source_cells)) unique(x[[source_col]]) else
    as.character(source_cells)
  target_cells <- if (is.null(target_cells)) unique(x[[target_col]]) else
    as.character(target_cells)
  unknown_source <- setdiff(unique(x[[source_col]]), source_cells)
  unknown_target <- setdiff(unique(x[[target_col]]), target_cells)
  if (length(unknown_source) || length(unknown_target)) {
    stop("transport contains source or target IDs absent from the supplied grid.",
         call. = FALSE)
  }
  source_has_route <- source_cells %in% x[[source_col]]
  target_has_carry <- target_cells %in% x[[target_col]]
  diagnostics <- data.frame(
    n_source_cells = length(source_cells),
    n_target_cells = length(target_cells),
    n_links = nrow(x),
    source_coverage_fraction = mean(source_has_route),
    target_coverage_fraction = mean(target_has_carry),
    n_unmapped_source_cells = sum(!source_has_route),
    n_unmapped_target_cells = sum(!target_has_carry),
    max_target_weight_sum = if (nrow(target_sum)) max(target_sum$target_weight_sum) else 0,
    stringsAsFactors = FALSE
  )
  structure(
    list(
      transport = x,
      target_weight_sums = target_sum,
      source_cells = source_cells,
      target_cells = target_cells,
      source_col = source_col,
      target_col = target_col,
      weight_col = weight_col,
      diagnostics = diagnostics,
      interpretation = paste(
        "Palaeogeographic carriage before active dispersal;",
        "weights are target-area-normalised mixture weights, not a dispersal kernel."
      )
    ),
    class = "hee_plate_grid_transport"
  )
}

#' Carry grid occupancy through an explicit palaeogeographic transport table
#'
#' Applies [hee_plate_grid_transport()] to a cell-by-lineage occupancy matrix.
#' New target land without mapped source substrate starts at zero; removed land
#' has no target row. The result is a probability-density remapping, so it is
#' executed before the dispersal/colonisation process and cannot be interpreted
#' as biological range expansion.
#'
#' @param occupancy Source cell-by-lineage occupancy probabilities.
#' @param transport A `hee_plate_grid_transport` object or a compatible table.
#' @param source_cells,target_cells Optional explicit row-order identifiers.
#' @return A list containing target `occupancy`, `mapped_target`, and transport
#'   diagnostics.
#' @export
hee_plate_carry_occupancy <- function(occupancy,
                                      transport,
                                      source_cells = rownames(occupancy),
                                      target_cells = NULL) {
  q <- as.matrix(occupancy)
  storage.mode(q) <- "double"
  if (is.null(source_cells) || length(source_cells) != nrow(q)) {
    stop("source_cells must identify every occupancy row.", call. = FALSE)
  }
  if (any(!is.finite(q) | q < 0 | q > 1)) {
    stop("occupancy must be finite and within [0, 1].", call. = FALSE)
  }
  tr <- if (inherits(transport, "hee_plate_grid_transport")) transport else {
    hee_plate_grid_transport(transport, source_cells = source_cells,
                             target_cells = target_cells)
  }
  if (is.null(target_cells)) target_cells <- tr$target_cells
  target_cells <- as.character(target_cells)
  if (anyDuplicated(source_cells) || anyDuplicated(target_cells)) {
    stop("source_cells and target_cells must be unique.", call. = FALSE)
  }
  # Revalidate against the actual matrix and requested target order.
  tr <- hee_plate_grid_transport(
    tr$transport, source_cells = source_cells, target_cells = target_cells,
    source_col = tr$source_col, target_col = tr$target_col,
    weight_col = tr$weight_col
  )
  x <- tr$transport
  from <- match(x[[tr$source_col]], source_cells)
  to <- match(x[[tr$target_col]], target_cells)
  if (requireNamespace("Matrix", quietly = TRUE)) {
    A <- Matrix::sparseMatrix(
      i = to, j = from, x = x[[tr$weight_col]],
      dims = c(length(target_cells), length(source_cells)),
      dimnames = list(target_cells, source_cells)
    )
    carried <- as.matrix(A %*% q)
  } else {
    carried <- matrix(0, nrow = length(target_cells), ncol = ncol(q),
                      dimnames = list(target_cells, colnames(q)))
    for (ii in seq_len(nrow(x))) {
      carried[to[[ii]], ] <- carried[to[[ii]], ] +
        x[[tr$weight_col]][[ii]] * q[from[[ii]], ]
    }
  }
  carried <- pmin(pmax(carried, 0), 1)
  rownames(carried) <- target_cells
  colnames(carried) <- colnames(q)
  mapped_target <- rowSums(as.matrix(tr$target_weight_sums[, "target_weight_sum", drop = FALSE]) > 0)
  names(mapped_target) <- tr$target_weight_sums$target_cell_id
  list(
    occupancy = carried,
    mapped_target = target_cells %in% names(mapped_target),
    diagnostics = tr$diagnostics,
    transport = tr
  )
}
