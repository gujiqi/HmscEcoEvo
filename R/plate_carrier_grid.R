#' Match a plate-carrier grid to one palaeoenvironment time slice
#'
#' A plate reconstruction can provide stable carrier identities (for example,
#' H3 IDs) whose coordinates change through time.  This helper uses those
#' identities as the Dynamic-Earth state space and samples a complete land
#' grid only to attach time-specific environmental covariates and habitat
#' status.  It avoids pretending that a target-centred raster interpolation is
#' a biological transport operator.
#'
#' Carriers farther than `max_match_distance_km` from any valid palaeo-land
#' cell are marked inactive.  They are not assigned an arbitrary nearest land
#' environment.  Larger `time_ma` values are older.
#'
#' @param carriers Data frame with stable `track_id`, `lon` and `lat` columns
#'   at one time slice.
#' @param grid Complete palaeoenvironment grid with `cell_id`, `lon`, `lat`,
#'   and a habitat-state column.
#' @param time_ma One time in Ma.
#' @param track_id_col,lon_col,lat_col,cell_id_col,habitat_col Column names.
#' @param max_match_distance_km Maximum spherical nearest-grid distance for an
#'   active carrier.
#' @return A data frame retaining carrier coordinates, matched grid identity,
#'   match distance, `active_land`, and all non-geometry grid covariates.
#' @export
hee_plate_carrier_match_grid <- function(
    carriers,
    grid,
    time_ma,
    track_id_col = "track_id",
    lon_col = "lon",
    lat_col = "lat",
    cell_id_col = "cell_id",
    habitat_col = "H_state",
    max_match_distance_km = 250) {
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("RANN is required for plate-carrier nearest-neighbour matching.",
         call. = FALSE)
  }
  cr <- as.data.frame(carriers, stringsAsFactors = FALSE)
  gr <- as.data.frame(grid, stringsAsFactors = FALSE)
  .require_cols(cr, c(track_id_col, lon_col, lat_col), "carriers")
  .require_cols(gr, c(cell_id_col, lon_col, lat_col, habitat_col), "grid")
  max_match_distance_km <- suppressWarnings(as.numeric(max_match_distance_km)[1L])
  if (!is.finite(max_match_distance_km) || max_match_distance_km <= 0) {
    stop("max_match_distance_km must be one positive finite number.", call. = FALSE)
  }
  cr <- cr[is.finite(cr[[lon_col]]) & is.finite(cr[[lat_col]]) &
             nzchar(as.character(cr[[track_id_col]])), , drop = FALSE]
  if (anyDuplicated(cr[[track_id_col]])) {
    stop("carriers must have one row per stable track_id at a time slice.",
         call. = FALSE)
  }
  land <- gr[is.finite(gr[[lon_col]]) & is.finite(gr[[lat_col]]) &
               is.finite(suppressWarnings(as.numeric(gr[[habitat_col]]))) &
               suppressWarnings(as.numeric(gr[[habitat_col]])) > 0, , drop = FALSE]
  if (!nrow(land)) stop("grid has no valid land cells.", call. = FALSE)
  xyz <- function(lon, lat) {
    rad <- pi / 180
    lon <- as.numeric(lon) * rad
    lat <- as.numeric(lat) * rad
    cbind(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat))
  }
  hit <- RANN::nn2(xyz(land[[lon_col]], land[[lat_col]]),
                   xyz(cr[[lon_col]], cr[[lat_col]]), k = 1L)
  chord_to_km <- function(x) 2 * 6371.0088 * asin(pmin(1, x / 2))
  matched <- land[as.integer(hit$nn.idx[, 1L]), , drop = FALSE]
  distance <- chord_to_km(hit$nn.dists[, 1L])
  active <- is.finite(distance) & distance <= max_match_distance_km
  covariate_cols <- setdiff(names(matched), c(cell_id_col, lon_col, lat_col,
                                               habitat_col, "time_ma"))
  out <- data.frame(
    track_id = as.character(cr[[track_id_col]]),
    time_ma = as.numeric(time_ma),
    paleo_lon = as.numeric(cr[[lon_col]]),
    paleo_lat = as.numeric(cr[[lat_col]]),
    matched_grid_cell_id = as.character(matched[[cell_id_col]]),
    match_distance_km = distance,
    active_land = active,
    stringsAsFactors = FALSE
  )
  for (nm in covariate_cols) {
    out[[nm]] <- matched[[nm]]
    out[[nm]][!active] <- NA
  }
  out$matched_grid_cell_id[!active] <- NA_character_
  attr(out, "diagnostics") <- data.frame(
    time_ma = as.numeric(time_ma),
    n_carriers = nrow(out),
    n_active_land_carriers = sum(active),
    active_land_fraction = mean(active),
    match_distance_km_q50 = as.numeric(stats::quantile(distance, 0.5)),
    match_distance_km_q95 = as.numeric(stats::quantile(distance, 0.95)),
    match_distance_km_max = max(distance),
    max_match_distance_km = max_match_distance_km,
    stringsAsFactors = FALSE
  )
  out
}

# Internal spherical-coordinate helper shared by plate output functions.
.hee_plate_xyz <- function(lon, lat) {
  rad <- pi / 180
  lon <- as.numeric(lon) * rad
  lat <- as.numeric(lat) * rad
  cbind(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat))
}

#' Construct a cached plate-carrier-to-grid output operator
#'
#' This creates a sparse, output-only interpolation operator from active stable
#' plate-carrier IDs to a complete palaeo-land grid. It is useful when multiple
#' posterior summaries are mapped at the same time slice: calculate geographic
#' neighbours once, then apply exactly the same operator to every metric. It is
#' not a plate transport, a movement kernel, or a biological state transition.
#'
#' Larger `time_ma` values are older. Target land cells farther than
#' `max_output_distance_km` remain uncovered and later receive `NA`, never a
#' fabricated biological zero.
#'
#' @param carrier_grid Output of [hee_plate_carrier_match_grid()].
#' @param target_grid Complete palaeoenvironment grid.
#' @param cell_id_col,lon_col,lat_col,habitat_col Grid column names.
#' @param k Number of nearest active carriers.
#' @param max_output_distance_km Maximum permitted nearest-carrier distance.
#' @return A `hee_plate_carrier_output_weights` object with a sparse
#'   target-by-carrier interpolation matrix, carrier IDs and coverage table.
#' @export
hee_plate_carrier_output_weights <- function(
    carrier_grid,
    target_grid,
    cell_id_col = "cell_id",
    lon_col = "lon",
    lat_col = "lat",
    habitat_col = "H_state",
    k = 4L,
    max_output_distance_km = 300) {
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("RANN is required for plate-carrier output interpolation.", call. = FALSE)
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix is required to cache plate-carrier output weights.", call. = FALSE)
  }
  cg <- as.data.frame(carrier_grid, stringsAsFactors = FALSE)
  tg <- as.data.frame(target_grid, stringsAsFactors = FALSE)
  .require_cols(cg, c("track_id", "paleo_lon", "paleo_lat", "active_land"),
                "carrier_grid")
  .require_cols(tg, c(cell_id_col, lon_col, lat_col, habitat_col), "target_grid")
  if (anyDuplicated(cg$track_id)) {
    stop("carrier_grid must contain one row per stable track_id.", call. = FALSE)
  }
  k <- suppressWarnings(as.integer(k)[1L])
  max_output_distance_km <- suppressWarnings(as.numeric(max_output_distance_km)[1L])
  if (!is.finite(k) || k < 1L || !is.finite(max_output_distance_km) ||
      max_output_distance_km <= 0) {
    stop("k and max_output_distance_km must be positive.", call. = FALSE)
  }
  active <- cg[!is.na(cg$active_land) & cg$active_land &
                 is.finite(cg$paleo_lon) & is.finite(cg$paleo_lat), , drop = FALSE]
  if (!nrow(active)) stop("No active carriers are available for output interpolation.", call. = FALSE)
  land_index <- which(is.finite(tg[[lon_col]]) & is.finite(tg[[lat_col]]) &
                        is.finite(suppressWarnings(as.numeric(tg[[habitat_col]]))) &
                        suppressWarnings(as.numeric(tg[[habitat_col]])) > 0)
  if (!length(land_index)) stop("target_grid has no valid land cells.", call. = FALSE)
  kk <- min(k, nrow(active))
  hit <- RANN::nn2(.hee_plate_xyz(active$paleo_lon, active$paleo_lat),
                   .hee_plate_xyz(tg[[lon_col]][land_index], tg[[lat_col]][land_index]),
                   k = kk)
  chord_to_km <- function(x) 2 * 6371.0088 * asin(pmin(1, x / 2))
  hit_index <- matrix(hit$nn.idx, nrow = length(land_index), ncol = kk)
  distance <- matrix(chord_to_km(as.numeric(hit$nn.dists)),
                     nrow = length(land_index), ncol = kk)
  nearest_distance <- distance[, 1L]
  covered <- is.finite(nearest_distance) & nearest_distance <= max_output_distance_km
  i <- integer(); j <- integer(); x <- numeric()
  for (row in which(covered)) {
    take <- which(distance[row, ] <= max_output_distance_km)
    if (!length(take)) next
    weight <- 1 / pmax(distance[row, take], 1)^2
    weight <- weight / sum(weight)
    i <- c(i, rep(land_index[[row]], length(take)))
    j <- c(j, hit_index[row, take])
    x <- c(x, weight)
  }
  weights <- Matrix::sparseMatrix(
    i = i, j = j, x = x,
    dims = c(nrow(tg), nrow(active)),
    dimnames = list(as.character(tg[[cell_id_col]]), as.character(active$track_id))
  )
  coverage <- data.frame(
    cell_id = as.character(tg[[cell_id_col]][land_index]),
    nearest_carrier_distance_km = nearest_distance,
    output_covered = covered,
    stringsAsFactors = FALSE
  )
  structure(list(
    weights = weights,
    carrier_track_id = as.character(active$track_id),
    target_cell_id = as.character(tg[[cell_id_col]]),
    coverage = coverage,
    diagnostics = data.frame(
      n_land_cells = length(land_index),
      n_covered_land_cells = sum(covered),
      output_coverage_fraction = mean(covered),
      nearest_carrier_km_q95 = as.numeric(stats::quantile(nearest_distance, .95)),
      max_output_distance_km = max_output_distance_km,
      k = kk,
      stringsAsFactors = FALSE
    )
  ), class = "hee_plate_carrier_output_weights")
}

#' Apply cached plate-carrier output weights
#'
#' Applies [hee_plate_carrier_output_weights()] to one or more finite
#' carrier-indexed summaries. This is an output interpolation only and cannot
#' be used as a plate-transport or occupancy-update operator.
#'
#' @param carrier_values Numeric matrix with unique carrier IDs in row names.
#' @param output_weights Object from [hee_plate_carrier_output_weights()].
#' @return A complete target-grid cell-by-summary matrix; uncovered land and
#'   non-land cells are `NA`.
#' @export
hee_plate_carrier_apply_output_weights <- function(carrier_values,
                                                    output_weights) {
  if (!inherits(output_weights, "hee_plate_carrier_output_weights")) {
    stop("output_weights must be produced by hee_plate_carrier_output_weights().",
         call. = FALSE)
  }
  value <- as.matrix(carrier_values)
  if (is.null(rownames(value)) || anyDuplicated(rownames(value))) {
    stop("carrier_values must have unique carrier track IDs as row names.",
         call. = FALSE)
  }
  if (any(!is.finite(value))) {
    stop("carrier_values must be finite before applying output weights.",
         call. = FALSE)
  }
  required <- output_weights$carrier_track_id
  missing <- setdiff(required, rownames(value))
  if (length(missing)) {
    stop("carrier_values is missing active carrier IDs: ",
         paste(utils::head(missing, 10L), collapse = ", "), call. = FALSE)
  }
  out <- as.matrix(output_weights$weights %*% value[required, , drop = FALSE])
  uncovered <- setdiff(seq_len(nrow(out)), match(output_weights$coverage$cell_id,
                                                  rownames(out)))
  if (length(uncovered)) out[uncovered, ] <- NA_real_
  covered_cells <- output_weights$coverage$cell_id[!output_weights$coverage$output_covered]
  if (length(covered_cells)) out[match(covered_cells, rownames(out)), ] <- NA_real_
  out
}

#' Interpolate plate-carrier values to a complete palaeo-land output grid
#'
#' This is an output interpolation only: active process states remain indexed
#' by stable plate-carrier IDs. For every target land cell, the function uses
#' inverse-distance weights over nearby active carriers. An output cell with no
#' carrier within `max_output_distance_km` is `NA`, never an invented zero.
#'
#' @param carrier_values Numeric matrix with carrier IDs in row names.
#' @param carrier_grid Output of [hee_plate_carrier_match_grid()].
#' @param target_grid Complete palaeoenvironment grid.
#' @param cell_id_col,lon_col,lat_col,habitat_col Grid column names.
#' @param k Number of nearest active carriers for interpolation.
#' @param max_output_distance_km Maximum permitted nearest-carrier distance.
#' @param require_probability If `TRUE` (default), require values in `[0, 1]`.
#'   Set to `FALSE` for additive summaries such as expected lineage richness;
#'   interpolation remains an output operation and never changes process states.
#' @return A list with a cell-by-lineage output matrix and a target coverage
#'   table.
#' @export
hee_plate_carrier_to_grid <- function(
    carrier_values,
    carrier_grid,
    target_grid,
    cell_id_col = "cell_id",
    lon_col = "lon",
    lat_col = "lat",
    habitat_col = "H_state",
    k = 4L,
    max_output_distance_km = 300,
    require_probability = TRUE) {
  value <- as.matrix(carrier_values)
  if (is.null(rownames(value)) || anyDuplicated(rownames(value))) {
    stop("carrier_values must have unique carrier track IDs as row names.",
         call. = FALSE)
  }
  if (any(!is.finite(value))) {
    stop("carrier_values must be finite.", call. = FALSE)
  }
  if (isTRUE(require_probability) && any(value < 0 | value > 1)) {
    stop("carrier_values must be probabilities in [0, 1] when require_probability = TRUE.",
         call. = FALSE)
  }
  operator <- hee_plate_carrier_output_weights(
    carrier_grid = carrier_grid, target_grid = target_grid,
    cell_id_col = cell_id_col, lon_col = lon_col, lat_col = lat_col,
    habitat_col = habitat_col, k = k,
    max_output_distance_km = max_output_distance_km
  )
  out <- hee_plate_carrier_apply_output_weights(value, operator)
  list(
    values = out,
    coverage = operator$coverage,
    diagnostics = operator$diagnostics
  )
}
