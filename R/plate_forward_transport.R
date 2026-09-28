#' Validate a source-conservative palaeogeographic transport table
#'
#' Plate carriage is a change of geographic reference frame, not organismal
#' dispersal.  This validator accepts a source-normalised table: every retained
#' source land-cell distributes all of its area to one or more younger target
#' land-cells.  Sources that are removed from the model must be declared
#' explicitly in `lost_sources`; they are never silently dropped.
#'
#' The area-conservative mapping divides the sum of
#' `source_area * source_weight * source_occupancy` by target cell area.
#'
#' Here, `w` is a source-area fraction and therefore sums to one within each
#' retained source cell.  The function checks that mapped source area does not
#' overfill a target cell.  It intentionally has a different class from the
#' legacy target-centred [hee_plate_grid_transport()] interface.
#'
#' @param transport A data frame with source and target IDs, source-normalised
#'   weights, and source/target cell areas in square kilometres.
#' @param source_cells,target_cells Optional complete source and target IDs.
#' @param lost_sources Optional data frame with `source_cell_id` and `reason`.
#'   Valid reasons are `"land_loss"` and `"outside_model_domain"`.
#' @param source_col,target_col,weight_col Column names in `transport`.
#' @param source_area_col,target_area_col Area column names in `transport`.
#' @param tolerance Numeric tolerance for weight and area-conservation checks.
#' @return A `hee_plate_forward_transport` object.
#' @export
hee_plate_forward_transport <- function(
    transport,
    source_cells = NULL,
    target_cells = NULL,
    lost_sources = NULL,
    source_col = "source_cell_id",
    target_col = "target_cell_id",
    weight_col = "source_weight",
    source_area_col = "source_area_km2",
    target_area_col = "target_area_km2",
    tolerance = 1e-6) {
  x <- as.data.frame(transport, stringsAsFactors = FALSE)
  .require_cols(x, c(source_col, target_col, weight_col,
                     source_area_col, target_area_col), "transport")
  for (nm in c(source_col, target_col)) x[[nm]] <- as.character(x[[nm]])
  for (nm in c(weight_col, source_area_col, target_area_col)) {
    x[[nm]] <- suppressWarnings(as.numeric(x[[nm]]))
  }
  if (any(!nzchar(x[[source_col]]) | !nzchar(x[[target_col]]))) {
    stop("transport contains blank source or target cell IDs.", call. = FALSE)
  }
  if (any(!is.finite(x[[weight_col]]) | x[[weight_col]] < 0) ||
      any(!is.finite(x[[source_area_col]]) | x[[source_area_col]] <= 0) ||
      any(!is.finite(x[[target_area_col]]) | x[[target_area_col]] <= 0)) {
    stop("transport weights must be non-negative and cell areas must be positive finite values.",
         call. = FALSE)
  }
  key <- paste(x[[source_col]], x[[target_col]], sep = "\r")
  if (anyDuplicated(key)) {
    stop("transport must contain at most one row per source-target pair.",
         call. = FALSE)
  }
  source_cells <- if (is.null(source_cells)) unique(x[[source_col]]) else {
    as.character(source_cells)
  }
  target_cells <- if (is.null(target_cells)) unique(x[[target_col]]) else {
    as.character(target_cells)
  }
  if (anyDuplicated(source_cells) || anyDuplicated(target_cells)) {
    stop("source_cells and target_cells must be unique.", call. = FALSE)
  }
  if (length(setdiff(unique(x[[source_col]]), source_cells)) ||
      length(setdiff(unique(x[[target_col]]), target_cells))) {
    stop("transport contains IDs absent from the supplied source or target grid.",
         call. = FALSE)
  }

  lost <- if (is.null(lost_sources)) {
    data.frame(source_cell_id = character(), reason = character(),
               stringsAsFactors = FALSE)
  } else {
    z <- as.data.frame(lost_sources, stringsAsFactors = FALSE)
    .require_cols(z, c("source_cell_id", "reason"), "lost_sources")
    z$source_cell_id <- as.character(z$source_cell_id)
    z$reason <- as.character(z$reason)
    if (anyDuplicated(z$source_cell_id) ||
        any(!z$reason %in% c("land_loss", "outside_model_domain"))) {
      stop("lost_sources must have unique source_cell_id and a valid loss reason.",
           call. = FALSE)
    }
    z
  }
  if (length(setdiff(lost$source_cell_id, source_cells))) {
    stop("lost_sources contains IDs absent from source_cells.", call. = FALSE)
  }
  if (length(intersect(lost$source_cell_id, unique(x[[source_col]])))) {
    stop("A source cannot be both transported and declared lost.", call. = FALSE)
  }

  source_sum <- stats::aggregate(x[[weight_col]],
                                 list(source_cell_id = x[[source_col]]), sum)
  names(source_sum)[[2L]] <- "source_weight_sum"
  if (any(abs(source_sum$source_weight_sum - 1) > tolerance)) {
    stop("Forward transport weights must sum to one within every retained source cell.",
         call. = FALSE)
  }
  unmapped_source <- setdiff(source_cells, source_sum$source_cell_id)
  undeclared <- setdiff(unmapped_source, lost$source_cell_id)
  if (length(undeclared)) {
    stop("Every source cell must be mapped or explicitly declared as land_loss/outside_model_domain. Undeclared examples: ",
         paste(utils::head(undeclared, 5L), collapse = ", "), ".", call. = FALSE)
  }

  target_area <- stats::aggregate(x[[target_area_col]],
                                  list(target_cell_id = x[[target_col]]), function(v) {
                                    if (max(v) - min(v) > tolerance) NA_real_ else v[[1L]]
                                  })
  names(target_area)[[2L]] <- "target_area_km2"
  if (any(!is.finite(target_area$target_area_km2))) {
    stop("Each target cell must have one consistent target area.", call. = FALSE)
  }
  target_load <- stats::aggregate(
    x[[source_area_col]] * x[[weight_col]],
    list(target_cell_id = x[[target_col]]), sum
  )
  names(target_load)[[2L]] <- "mapped_source_area_km2"
  target_check <- merge(target_area, target_load, by = "target_cell_id", all = TRUE,
                        sort = FALSE)
  target_check$area_load_ratio <- target_check$mapped_source_area_km2 /
    target_check$target_area_km2
  if (any(target_check$area_load_ratio > 1 + tolerance)) {
    bad <- target_check$target_cell_id[target_check$area_load_ratio > 1 + tolerance]
    stop("Forward transport overfills at least one target cell; use polygon overlap or a finer source grid. Examples: ",
         paste(utils::head(bad, 5L), collapse = ", "), ".", call. = FALSE)
  }
  source_area <- stats::aggregate(x[[source_area_col]],
                                  list(source_cell_id = x[[source_col]]), function(v) {
                                    if (max(v) - min(v) > tolerance) NA_real_ else v[[1L]]
                                  })
  names(source_area)[[2L]] <- "source_area_km2"
  if (any(!is.finite(source_area$source_area_km2))) {
    stop("Each source cell must have one consistent source area.", call. = FALSE)
  }
  diagnostics <- data.frame(
    n_source_cells = length(source_cells),
    n_target_cells = length(target_cells),
    n_links = nrow(x),
    n_retained_source_cells = nrow(source_sum),
    n_declared_lost_source_cells = nrow(lost),
    source_coverage_fraction = nrow(source_sum) / length(source_cells),
    target_coverage_fraction = length(unique(x[[target_col]])) / length(target_cells),
    max_target_area_load_ratio = max(target_check$area_load_ratio),
    min_target_area_load_ratio = min(target_check$area_load_ratio),
    stringsAsFactors = FALSE
  )
  structure(list(
    transport = x,
    lost_sources = lost,
    source_weight_sums = source_sum,
    source_areas = source_area,
    target_area_load = target_check,
    source_cells = source_cells,
    target_cells = target_cells,
    source_col = source_col,
    target_col = target_col,
    weight_col = weight_col,
    source_area_col = source_area_col,
    target_area_col = target_area_col,
    diagnostics = diagnostics,
    interpretation = paste(
      "Source-normalised, area-conservative plate carriage before active dispersal;",
      "not a movement or colonisation kernel."
    )
  ), class = "hee_plate_forward_transport")
}

#' Carry occupancy forward through source-conservative plate transport
#'
#' Applies [hee_plate_forward_transport()] to a cell-by-lineage occupancy
#' matrix.  The returned `transported_area_km2` is an expected occupied-area
#' accounting diagnostic.  New target land receives zero carried occupancy and
#' can only become occupied through later colonisation.
#'
#' @param occupancy Source cell-by-lineage occupancy probabilities.
#' @param transport A `hee_plate_forward_transport` object.
#' @param source_cells,target_cells Explicit row IDs in matrix order.
#' @return A list with target occupancy, occupied-area diagnostics, new-land
#'   flags and the validated transport object.
#' @export
hee_plate_carry_forward_occupancy <- function(occupancy,
                                              transport,
                                              source_cells = rownames(occupancy),
                                              target_cells = NULL) {
  q <- as.matrix(occupancy)
  storage.mode(q) <- "double"
  if (is.null(source_cells) || length(source_cells) != nrow(q) ||
      anyDuplicated(source_cells)) {
    stop("source_cells must uniquely identify every occupancy row.", call. = FALSE)
  }
  if (any(!is.finite(q) | q < 0 | q > 1)) {
    stop("occupancy must be finite and within [0, 1].", call. = FALSE)
  }
  if (!inherits(transport, "hee_plate_forward_transport")) {
    stop("transport must be created by hee_plate_forward_transport().", call. = FALSE)
  }
  if (is.null(target_cells)) target_cells <- transport$target_cells
  target_cells <- as.character(target_cells)
  if (anyDuplicated(target_cells)) stop("target_cells must be unique.", call. = FALSE)
  if (!identical(as.character(source_cells), transport$source_cells) ||
      !identical(target_cells, transport$target_cells)) {
    stop("Matrix cell orders must exactly match the validated forward transport object.",
         call. = FALSE)
  }
  x <- transport$transport
  from <- match(x[[transport$source_col]], source_cells)
  to <- match(x[[transport$target_col]], target_cells)
  amount <- x[[transport$weight_col]] * x[[transport$source_area_col]]
  target_area <- transport$target_area_load$target_area_km2[
    match(target_cells, transport$target_area_load$target_cell_id)
  ]
  if (anyNA(target_area)) {
    target_area[is.na(target_area)] <- 1
  }
  if (requireNamespace("Matrix", quietly = TRUE)) {
    A <- Matrix::sparseMatrix(i = to, j = from, x = amount,
                              dims = c(length(target_cells), length(source_cells)),
                              dimnames = list(target_cells, source_cells))
    transported_area <- as.matrix(A %*% q)
  } else {
    transported_area <- matrix(0, nrow = length(target_cells), ncol = ncol(q),
                               dimnames = list(target_cells, colnames(q)))
    for (ii in seq_len(nrow(x))) {
      transported_area[to[[ii]], ] <- transported_area[to[[ii]], ] +
        amount[[ii]] * q[from[[ii]], ]
    }
  }
  carried <- transported_area / target_area
  if (any(carried > 1 + 1e-6)) {
    stop("Forward transport produced occupancy above one; target area conservation failed.",
         call. = FALSE)
  }
  carried <- pmin(pmax(carried, 0), 1)
  dimnames(carried) <- list(target_cells, colnames(q))
  mapped <- target_cells %in% x[[transport$target_col]]
  source_area <- transport$source_areas$source_area_km2[
    match(source_cells, transport$source_areas$source_cell_id)
  ]
  list(
    occupancy = carried,
    transported_area_km2 = transported_area,
    source_expected_area_km2 = colSums(q * source_area),
    target_expected_area_km2 = colSums(transported_area),
    mapped_target = mapped,
    new_target_land = !mapped,
    diagnostics = transport$diagnostics,
    transport = transport
  )
}
