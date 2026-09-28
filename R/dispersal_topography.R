.hee_movement_source_groups <- function(x) {
  group_cols <- c(intersect(c("lineage", "time_ma", "response_draw"),
                             names(x)), "from_cell_id")
  key <- do.call(paste, c(lapply(x[group_cols], as.character), sep = "\r"))
  index <- match(key, unique(key))
  list(index = index, n_groups = max(index), columns = group_cols)
}

.hee_group_sum <- function(x, group_index, n_groups) {
  # `rowsum()` is implemented in C and avoids the very expensive factor
  # interaction/ave path on global 1-degree movement graphs.
  as.numeric(rowsum(as.numeric(x), group = group_index, reorder = FALSE))
}

#' Build a continuous elevation-derived dispersal-resistance surface
#'
#' Builds the cell-level part of the HmscEE dispersal process from a continuous
#' palaeo-elevation grid.  For every land cell `c` at time `t`, local relief is
#' the median absolute elevation difference from its four or eight valid land
#' neighbours:
#' \deqn{Relief_{c,t}=\operatorname{median}_{j\in N(c)}|z_{j,t}-z_{c,t}|.}
#' The default structural resistance is then
#' \deqn{R_{c,t}=\log[1+Relief_{c,t}/z_R].}
#'
#' `R` is not suitability, occurrence probability, or a calibrated dispersal
#' rate.  It is a continuous, taxon-neutral topographic input to P2
#' dispersal.  Ocean cells are hard barriers.  Absolute elevation is not used
#' by default: a high but flat plateau is not automatically a movement barrier,
#' and absolute elevation is normally already represented by environmental
#' filtering.  An explicit absolute-elevation term is available only as a
#' documented sensitivity scenario.
#'
#' @param landscape A cell-time data frame with cell ID, time, longitude,
#'   latitude, elevation in metres, and a land mask.
#' @param elevation_col,land_col,cell_col,time_col,lon_col,lat_col Column
#'   names. `time_ma` is in Ma, with larger values older.
#' @param neighbours Either `4` or `8` neighbour movement.  The default `8`
#'   follows a queen-contiguity grid graph.
#' @param resolution_deg Grid resolution. When `NULL`, it is inferred from
#'   coordinate spacings separately for longitude and latitude.
#' @param relief_ref_m Positive relief scaling constant in metres. It changes
#'   the numerical scale only; it is not a biological threshold. Sensitivity
#'   analyses should compare plausible values such as 250, 500, and 1000 m.
#' @param absolute_elevation Logical; default `FALSE`.
#' @param elevation_ref_m,elevation_weight Optional absolute-elevation scenario
#'   settings. Ignored unless `absolute_elevation = TRUE`.
#' @param missing_elevation Either `"barrier"` (default) or `"error"`.
#'
#' @return An object of class `hee_dispersal_resistance` with `cells` and
#'   `metadata`. `cells` includes `local_relief_m`,
#'   `topographic_resistance`, `movement_allowed`, and terrain diagnostics.
#' @export
hee_dispersal_resistance_surface <- function(landscape,
                                             elevation_col = "elevation_m",
                                             land_col = "land_mask_dem",
                                             cell_col = "cell_id",
                                             time_col = "time_ma",
                                             lon_col = "lon",
                                             lat_col = "lat",
                                             neighbours = 8L,
                                             resolution_deg = NULL,
                                             relief_ref_m = 500,
                                             absolute_elevation = FALSE,
                                             elevation_ref_m = 5000,
                                             elevation_weight = 0,
                                             missing_elevation = c("barrier", "error")) {
  x <- as.data.frame(landscape)
  .require_cols(x, c(cell_col, time_col, lon_col, lat_col, elevation_col,
                     land_col), "landscape")
  .hee_validate_time_values(x[[time_col]], paste0("landscape$", time_col))
  .hee_check_unique_keys(x, c(cell_col, time_col), "landscape")
  missing_elevation <- match.arg(missing_elevation)
  neighbours <- as.integer(neighbours)[1L]
  if (!neighbours %in% c(4L, 8L)) {
    stop("neighbours must be either 4 or 8.", call. = FALSE)
  }
  if (!is.finite(relief_ref_m) || relief_ref_m <= 0 ||
      !is.finite(elevation_ref_m) || elevation_ref_m <= 0 ||
      !is.finite(elevation_weight) || elevation_weight < 0) {
    stop("Reference scales must be positive and elevation_weight non-negative.",
         call. = FALSE)
  }

  lon <- suppressWarnings(as.numeric(x[[lon_col]]))
  lat <- suppressWarnings(as.numeric(x[[lat_col]]))
  elev <- suppressWarnings(as.numeric(x[[elevation_col]]))
  land <- suppressWarnings(as.numeric(x[[land_col]]))
  if (any(!is.finite(lon) | !is.finite(lat))) {
    stop("Longitude and latitude must be finite.", call. = FALSE)
  }
  land_valid <- is.finite(land) & land > 0
  terrain_missing <- land_valid & !is.finite(elev)
  if (identical(missing_elevation, "error") && any(terrain_missing)) {
    stop("Land cells have missing elevation.", call. = FALSE)
  }

  infer_resolution <- function(v, label) {
    u <- sort(unique(v[is.finite(v)]))
    d <- diff(u)
    d <- d[d > sqrt(.Machine$double.eps)]
    if (!length(d)) stop("Cannot infer ", label, " resolution.", call. = FALSE)
    stats::median(d)
  }
  if (is.null(resolution_deg)) {
    dx <- infer_resolution(lon, "longitude")
    dy <- infer_resolution(lat, "latitude")
  } else {
    supplied <- as.numeric(resolution_deg)
    if (length(supplied) == 1L) supplied <- rep(supplied, 2L)
    if (length(supplied) != 2L || any(!is.finite(supplied) | supplied <= 0)) {
      stop("resolution_deg must be one positive value or c(longitude, latitude).",
           call. = FALSE)
    }
    dx <- supplied[[1L]]
    dy <- supplied[[2L]]
  }

  shifts <- expand.grid(dx = -1:1, dy = -1:1,
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  shifts <- shifts[!(shifts$dx == 0 & shifts$dy == 0), , drop = FALSE]
  if (neighbours == 4L) shifts <- shifts[abs(shifts$dx) + abs(shifts$dy) == 1, ]
  coordinate_key <- function(a, b) {
    paste(formatC(a, format = "f", digits = 8),
          formatC(b, format = "f", digits = 8), sep = "\r")
  }

  local_relief <- rep(NA_real_, nrow(x))
  n_neighbours <- integer(nrow(x))
  times <- unique(x[[time_col]])
  for (tm in times) {
    idx <- which(x[[time_col]] == tm & land_valid & is.finite(elev))
    if (!length(idx)) next
    key <- coordinate_key(lon[idx], lat[idx])
    differences <- matrix(NA_real_, nrow = length(idx), ncol = nrow(shifts))
    for (shift_id in seq_len(nrow(shifts))) {
      nlon <- lon[idx] + shifts$dx[[shift_id]] * dx
      nlon[nlon > 180] <- nlon[nlon > 180] - 360
      nlon[nlon < -180] <- nlon[nlon < -180] + 360
      nlat <- lat[idx] + shifts$dy[[shift_id]] * dy
      hit <- match(coordinate_key(nlon, nlat), key)
      keep <- !is.na(hit)
      if (any(keep)) {
        differences[keep, shift_id] <- abs(elev[idx[hit[keep]]] - elev[idx[keep]])
      }
    }
    n_neighbours[idx] <- rowSums(is.finite(differences))
    local_relief[idx] <- apply(differences, 1L, function(z) {
      z <- z[is.finite(z)]
      if (length(z)) stats::median(z) else 0
    })
  }
  movement_allowed <- as.integer(land_valid & is.finite(elev))
  topographic_resistance <- rep(NA_real_, nrow(x))
  valid <- movement_allowed > 0
  topographic_resistance[valid] <- log1p(local_relief[valid] / relief_ref_m)
  if (isTRUE(absolute_elevation)) {
    topographic_resistance[valid] <- topographic_resistance[valid] +
      elevation_weight * log1p(pmax(elev[valid], 0) / elevation_ref_m)
  }

  cells <- x
  cells$local_relief_m <- local_relief
  cells$n_terrain_neighbours <- n_neighbours
  cells$topographic_resistance <- topographic_resistance
  cells$topographic_permeability <- ifelse(
    is.finite(topographic_resistance), exp(-topographic_resistance), 0
  )
  cells$movement_allowed <- movement_allowed
  cells$terrain_missing <- terrain_missing
  metadata <- data.frame(
    field = c("formula", "relief_method", "neighbours", "resolution_lon_deg",
              "resolution_lat_deg", "relief_ref_m", "absolute_elevation",
              "elevation_ref_m", "elevation_weight", "ocean_rule",
              "scientific_boundary"),
    value = c(
      if (absolute_elevation) {
        "log1p(local_relief_m / relief_ref_m) + elevation_weight * log1p(max(elevation_m, 0) / elevation_ref_m)"
      } else "log1p(local_relief_m / relief_ref_m)",
      "median_abs_difference_to_valid_land_neighbours", neighbours, dx, dy,
      relief_ref_m, absolute_elevation, elevation_ref_m, elevation_weight,
      "Ocean and land cells with missing elevation have no local movement edge.",
      "Continuous structural resistance scenario; not a calibrated species-specific dispersal, gene-flow, or occurrence probability."
    ),
    stringsAsFactors = FALSE
  )
  structure(list(cells = cells, metadata = metadata),
            class = "hee_dispersal_resistance")
}

#' Build a topographic cell-to-cell connectivity graph
#'
#' Converts a [hee_dispersal_resistance_surface()] result into the sparse,
#' directed 4- or 8-neighbour graph actually used by local dispersal. For every
#' valid edge, the elevation step and barrier are
#' \deqn{S_{ij,t}=\log[1+|z_{i,t}-z_{j,t}|/z_S],}
#' \deqn{B_{ij,t}=w_R(R_{i,t}+R_{j,t})/2+w_SS_{ij,t},}
#' and the effective connection cost is
#' \deqn{C_{ij,t}=d_{ij,t}\exp(\omega_{topo}B_{ij,t}).}
#'
#' The function stores one-edge costs, not all-pairs least-cost paths. In a
#' time-stepped occupancy model repeated local transitions already represent
#' movement through the graph without assuming that an organism knows a
#' destination in advance.
#'
#' @param resistance Output from [hee_dispersal_resistance_surface()] or its
#'   `cells` table.
#' @param cell_col,time_col,lon_col,lat_col,elevation_col Column names.
#' @param neighbours,resolution_deg Grid settings. They should match the
#'   resistance-surface settings.
#' @param step_ref_m Positive elevation-step scale in metres.
#' @param weight_relief,weight_step Non-negative barrier weights.
#' @param omega_topo Non-negative topographic cost strength. With zero, the
#'   effective cost equals geographic distance exactly.
#'
#' @return An object of class `hee_dispersal_connectivity` with `edges`,
#'   `diagnostics`, and `metadata`.
#' @export
hee_dispersal_connectivity_graph <- function(resistance,
                                             cell_col = "cell_id",
                                             time_col = "time_ma",
                                             lon_col = "lon",
                                             lat_col = "lat",
                                             elevation_col = "elevation_m",
                                             neighbours = 8L,
                                             resolution_deg = NULL,
                                             step_ref_m = 500,
                                             weight_relief = 0.5,
                                             weight_step = 0.5,
                                             omega_topo = 1) {
  cells <- if (inherits(resistance, "hee_dispersal_resistance")) {
    resistance$cells
  } else as.data.frame(resistance)
  .require_cols(cells, c(cell_col, time_col, lon_col, lat_col, elevation_col,
                         "topographic_resistance", "movement_allowed"),
                "resistance")
  .hee_check_unique_keys(cells, c(cell_col, time_col), "resistance")
  .hee_validate_time_values(cells[[time_col]], paste0("resistance$", time_col))
  neighbours <- as.integer(neighbours)[1L]
  if (!neighbours %in% c(4L, 8L)) stop("neighbours must be 4 or 8.", call. = FALSE)
  settings <- c(step_ref_m, weight_relief, weight_step, omega_topo)
  if (any(!is.finite(settings)) || step_ref_m <= 0 || weight_relief < 0 ||
      weight_step < 0 || omega_topo < 0) {
    stop("Invalid connectivity-graph settings.", call. = FALSE)
  }
  lon <- as.numeric(cells[[lon_col]])
  lat <- as.numeric(cells[[lat_col]])
  elev <- as.numeric(cells[[elevation_col]])
  infer_resolution <- function(v) {
    u <- sort(unique(v[is.finite(v)])); d <- diff(u); d <- d[d > 1e-8]
    if (!length(d)) stop("Cannot infer grid resolution.", call. = FALSE)
    stats::median(d)
  }
  if (is.null(resolution_deg)) {
    dx <- infer_resolution(lon); dy <- infer_resolution(lat)
  } else {
    supplied <- as.numeric(resolution_deg)
    if (length(supplied) == 1L) supplied <- rep(supplied, 2L)
    if (length(supplied) != 2L || any(!is.finite(supplied) | supplied <= 0)) {
      stop("resolution_deg must be one positive value or c(longitude, latitude).",
           call. = FALSE)
    }
    dx <- supplied[[1L]]; dy <- supplied[[2L]]
  }
  shifts <- expand.grid(dx = -1:1, dy = -1:1,
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  shifts <- shifts[!(shifts$dx == 0 & shifts$dy == 0), , drop = FALSE]
  if (neighbours == 4L) shifts <- shifts[abs(shifts$dx) + abs(shifts$dy) == 1, ]
  key <- function(a, b) paste(formatC(a, format = "f", digits = 8),
                               formatC(b, format = "f", digits = 8), sep = "\r")
  great_circle_km <- function(lon1, lat1, lon2, lat2) {
    radians <- pi / 180
    p1 <- lat1 * radians; p2 <- lat2 * radians
    dp <- (lat2 - lat1) * radians; dl <- (lon2 - lon1) * radians
    a <- sin(dp / 2)^2 + cos(p1) * cos(p2) * sin(dl / 2)^2
    2 * 6371.0088 * atan2(sqrt(a), sqrt(pmax(0, 1 - a)))
  }
  all_edges <- list(); row_id <- 0L
  for (tm in unique(cells[[time_col]])) {
    valid <- which(cells[[time_col]] == tm & cells$movement_allowed > 0 &
                     is.finite(cells$topographic_resistance) & is.finite(elev))
    if (!length(valid)) next
    loc_key <- key(lon[valid], lat[valid])
    for (shift_id in seq_len(nrow(shifts))) {
      nlon <- lon[valid] + shifts$dx[[shift_id]] * dx
      nlon[nlon > 180] <- nlon[nlon > 180] - 360
      nlon[nlon < -180] <- nlon[nlon < -180] + 360
      nlat <- lat[valid] + shifts$dy[[shift_id]] * dy
      hit <- match(key(nlon, nlat), loc_key)
      keep <- !is.na(hit)
      if (!any(keep)) next
      from <- valid[keep]; to <- valid[hit[keep]]
      mean_resistance <- (cells$topographic_resistance[from] +
                            cells$topographic_resistance[to]) / 2
      step <- abs(elev[from] - elev[to])
      step_index <- log1p(step / step_ref_m)
      barrier <- weight_relief * mean_resistance + weight_step * step_index
      distance <- great_circle_km(lon[from], lat[from], lon[to], lat[to])
      row_id <- row_id + 1L
      all_edges[[row_id]] <- data.frame(
        from_cell_id = as.character(cells[[cell_col]][from]),
        to_cell_id = as.character(cells[[cell_col]][to]),
        time_ma = as.numeric(cells[[time_col]][from]),
        distance_km = distance,
        elev_from_m = elev[from], elev_to_m = elev[to],
        elevation_step_m = step,
        resistance_from = cells$topographic_resistance[from],
        resistance_to = cells$topographic_resistance[to],
        mean_cell_resistance = mean_resistance,
        step_index = step_index,
        barrier_index = barrier,
        effective_cost_km = distance * exp(omega_topo * barrier),
        stringsAsFactors = FALSE
      )
    }
  }
  edges <- if (length(all_edges)) do.call(rbind, all_edges) else data.frame(
    from_cell_id = character(), to_cell_id = character(), time_ma = numeric(),
    distance_km = numeric(), elev_from_m = numeric(), elev_to_m = numeric(),
    elevation_step_m = numeric(), resistance_from = numeric(),
    resistance_to = numeric(), mean_cell_resistance = numeric(),
    step_index = numeric(), barrier_index = numeric(), effective_cost_km = numeric()
  )
  .hee_check_unique_keys(edges, c("from_cell_id", "to_cell_id", "time_ma"),
                         "topographic connectivity edges")
  diagnostics <- stats::aggregate(
    edges$effective_cost_km,
    edges["time_ma"],
    function(v) c(n_edges = length(v), mean_effective_cost_km = mean(v),
                  max_effective_cost_km = max(v))
  )
  if (nrow(diagnostics)) {
    diagnostics <- cbind(diagnostics["time_ma"], as.data.frame(diagnostics$x))
  } else diagnostics <- data.frame(time_ma = numeric(), n_edges = integer(),
                                    mean_effective_cost_km = numeric(),
                                    max_effective_cost_km = numeric())
  metadata <- data.frame(
    field = c("barrier_formula", "effective_cost_formula", "neighbours",
              "resolution_lon_deg", "resolution_lat_deg", "step_ref_m",
              "weight_relief", "weight_step", "omega_topo", "ocean_rule"),
    value = c(
      "weight_relief * mean(cell resistance) + weight_step * log1p(elevation_step_m / step_ref_m)",
      "distance_km * exp(omega_topo * barrier_index)", neighbours, dx, dy,
      step_ref_m, weight_relief, weight_step, omega_topo,
      "Only valid land-to-land edges are represented; ordinary local movement cannot cross ocean."
    ), stringsAsFactors = FALSE
  )
  structure(list(edges = edges, diagnostics = diagnostics, metadata = metadata),
            class = "hee_dispersal_connectivity")
}

#' Convert topographic connection costs into a legacy dispersal kernel
#'
#' This compatibility kernel uses an exponential effective-cost allocation. New
#' landscape-explicit HmscEE analyses should prefer
#' [hee_dispersal_spherical_kernel()], which makes the spherical diffusion
#' parameterisation, target landscape measure, and separation from plate
#' carriage explicit. This function remains available for reproducibility of
#' existing Case 03--05 workflows.
#'
#' For a lineage-specific effective dispersal scale `delta`, the local kernel
#' uses \deqn{W_{\ell,ij,t}=\exp[-C_{ij,t}/\delta_{\ell,t}].} Optional
#' long-distance-dispersal (LDD) edges use the geographic distance and an
#' ocean-crossing penalty.  These components are combined before optional
#' source normalisation.  The result is a conditional allocation of a source's
#' dispersal effort, not a realised colonisation probability.
#'
#' @param connectivity Output from [hee_dispersal_connectivity_graph()] or an
#'   edge table with `effective_cost_km`.
#' @param lineage_scales Optional lineage-time table with `lineage`, `time_ma`,
#'   and `dispersal_scale_km`. Omit it to use one `dispersal_scale_km` for all
#'   edges.
#' @param dispersal_scale_km Positive fixed local scale in km.
#' @param ldd_edges Optional long-distance edges with `from_cell_id`,
#'   `to_cell_id`, `time_ma`, `distance_km`, and optional `ocean_fraction`.
#' @param ldd_fraction Fraction of movement effort allocated to LDD; it is a
#'   sensitivity parameter, not an empirically estimated probability by
#'   default.
#' @param ldd_scale_km,ocean_penalty LDD distance scale and penalty.
#' @param normalise_by_source If `TRUE` (default), outgoing combined weights
#'   sum to one within each lineage-time-source group.
#' @param max_rows Maximum permitted expanded lineage-edge rows.
#'
#' @return Edge table with `local_weight`, `ldd_weight`, `movement_weight`, and
#'   `K_movement`, suitable for [hee_dispersal_arrival_pressure()].
#' @export
hee_dispersal_edge_kernel <- function(connectivity,
                                      lineage_scales = NULL,
                                      dispersal_scale_km = 200,
                                      ldd_edges = NULL,
                                      ldd_fraction = 0,
                                      ldd_scale_km = 2000,
                                      ocean_penalty = 0,
                                      normalise_by_source = TRUE,
                                      max_rows = 5000000L) {
  local <- if (inherits(connectivity, "hee_dispersal_connectivity")) {
    connectivity$edges
  } else as.data.frame(connectivity)
  .require_cols(local, c("from_cell_id", "to_cell_id", "time_ma",
                         "effective_cost_km"), "connectivity")
  .hee_check_unique_keys(local, c("from_cell_id", "to_cell_id", "time_ma"),
                         "connectivity")
  settings <- c(dispersal_scale_km, ldd_fraction, ldd_scale_km, ocean_penalty)
  if (any(!is.finite(settings)) || dispersal_scale_km <= 0 ||
      ldd_scale_km <= 0 || ldd_fraction < 0 || ldd_fraction > 1 ||
      ocean_penalty < 0) {
    stop("Invalid dispersal-kernel settings.", call. = FALSE)
  }
  local$edge_type <- "local"
  local$local_weight <- exp(-pmax(local$effective_cost_km, 0) / dispersal_scale_km)
  local$ldd_weight <- 0
  local$movement_weight <- (1 - ldd_fraction) * local$local_weight
  edge_data <- local[, c("from_cell_id", "to_cell_id", "time_ma", "edge_type",
                         "effective_cost_km", "local_weight", "ldd_weight",
                         "movement_weight"),
                     drop = FALSE]

  if (!is.null(ldd_edges)) {
    ldd <- as.data.frame(ldd_edges)
    .require_cols(ldd, c("from_cell_id", "to_cell_id", "time_ma", "distance_km"),
                  "ldd_edges")
    .hee_check_unique_keys(ldd, c("from_cell_id", "to_cell_id", "time_ma"),
                           "ldd_edges")
    ocean_fraction <- if ("ocean_fraction" %in% names(ldd)) {
      .hee_clip01(suppressWarnings(as.numeric(ldd$ocean_fraction)))
    } else rep(0, nrow(ldd))
    ldd$edge_type <- "ldd"
    ldd$effective_cost_km <- NA_real_
    ldd$local_weight <- 0
    ldd$ldd_weight <- exp(-pmax(as.numeric(ldd$distance_km), 0) / ldd_scale_km) *
      exp(-ocean_penalty * ocean_fraction)
    ldd$movement_weight <- ldd_fraction * ldd$ldd_weight
    edge_data <- rbind(edge_data, ldd[, names(edge_data), drop = FALSE])
  }
  if (is.null(lineage_scales)) {
    edge_data$lineage <- "global"
    edge_data$dispersal_scale_km <- dispersal_scale_km
  } else {
    scales <- as.data.frame(lineage_scales)
    .require_cols(scales, c("lineage", "time_ma", "dispersal_scale_km"),
                  "lineage_scales")
    .hee_check_unique_keys(scales, c("lineage", "time_ma"), "lineage_scales")
    if (any(!is.finite(scales$dispersal_scale_km) | scales$dispersal_scale_km <= 0)) {
      stop("lineage dispersal scales must be finite and positive.", call. = FALSE)
    }
    expected <- sum(vapply(unique(edge_data$time_ma), function(tm) {
      sum(edge_data$time_ma == tm) * sum(scales$time_ma == tm)
    }, numeric(1)))
    if (expected > max_rows) {
      stop("Lineage-edge expansion exceeds max_rows; run in time/lineage blocks.",
           call. = FALSE)
    }
    expanded <- lapply(seq_len(nrow(scales)), function(i) {
      z <- edge_data[edge_data$time_ma == scales$time_ma[[i]], , drop = FALSE]
      z$lineage <- scales$lineage[[i]]
      z$dispersal_scale_km <- scales$dispersal_scale_km[[i]]
      z$local_weight <- ifelse(z$edge_type == "local",
                               exp(-pmax(z$effective_cost_km, 0) /
                                     z$dispersal_scale_km), z$local_weight)
      z$movement_weight <- ifelse(z$edge_type == "local",
                                  (1 - ldd_fraction) * z$local_weight,
                                  ldd_fraction * z$ldd_weight)
      z
    })
    edge_data <- do.call(rbind, expanded)
  }
  source_groups <- .hee_movement_source_groups(edge_data)
  denominator <- .hee_group_sum(edge_data$movement_weight,
                                source_groups$index,
                                source_groups$n_groups)[source_groups$index]
  edge_data$K_movement <- if (isTRUE(normalise_by_source)) {
    ifelse(denominator > 0, edge_data$movement_weight / denominator, 0)
  } else .hee_clip01(edge_data$movement_weight)
  edge_data$kernel_normalised_by_source <- isTRUE(normalise_by_source)
  edge_data
}

#' Convert relative edge weights to a source-conserving movement-rate kernel
#'
#' A cell graph describes *where* a lineage can move, but relative edge weights
#' alone must not determine how much total movement a source cell produces. This
#' helper separates those two quantities.  For source cell \eqn{i}, it first
#' derives a directional allocation \eqn{a_{i\to j}} that sums to one and then
#' assigns one total source rate \eqn{m_i}.  The returned edge rate is
#' \deqn{m_{i\to j}=m_i a_{i\to j}.}
#'
#' With the default `source_permeability = "mean_weight"`,
#' \eqn{m_i=m_0\bar w_i}, where \eqn{\bar w_i} is the mean unnormalised edge
#' permeability leaving source \eqn{i}.  Consequently, adding more neighbours
#' does not multiply the source's total dispersal rate, while high local
#' topographic cost can still reduce it.  This is a movement-rate construction,
#' not a realised colonisation probability or an estimate of gene flow.
#'
#' @param movement_kernel A table from [hee_dispersal_edge_kernel()]. It needs
#'   `from_cell_id`, `to_cell_id`, and `K_movement`; `movement_weight` is used
#'   when available to determine source permeability.
#' @param base_rate Positive total movement-rate ceiling per Myr before the
#'   source permeability attenuation.
#' @param source_permeability One of `"mean_weight"` (default),
#'   `"sum_capped"`, or `"none"`. The first two use unnormalised movement
#'   weights when supplied; `"none"` gives every source `base_rate`.
#' @return `movement_kernel` augmented with `source_direction_weight`,
#'   `source_permeability`, `source_total_rate`, and `movement_rate`.
#' @export
hee_dispersal_rate_kernel <- function(movement_kernel,
                                      base_rate = 1,
                                      source_permeability = c("mean_weight",
                                                              "sum_capped",
                                                              "none")) {
  edges <- as.data.frame(movement_kernel)
  .require_cols(edges, c("from_cell_id", "to_cell_id", "K_movement"),
                "movement_kernel")
  base_rate <- suppressWarnings(as.numeric(base_rate)[1L])
  if (!is.finite(base_rate) || base_rate < 0) {
    stop("base_rate must be a finite non-negative movement rate.", call. = FALSE)
  }
  source_permeability <- match.arg(source_permeability)
  edges$K_movement <- suppressWarnings(as.numeric(edges$K_movement))
  edges$K_movement[!is.finite(edges$K_movement) | edges$K_movement < 0] <- 0
  source_groups <- .hee_movement_source_groups(edges)
  .hee_check_unique_keys(edges, c(source_groups$columns, "to_cell_id"),
                         "movement_kernel")
  direction_total <- .hee_group_sum(edges$K_movement,
                                    source_groups$index,
                                    source_groups$n_groups)[source_groups$index]
  edges$source_direction_weight <- ifelse(
    direction_total > 0, edges$K_movement / direction_total, 0
  )

  has_raw_weight <- "movement_weight" %in% names(edges)
  raw_weight <- if (has_raw_weight) suppressWarnings(as.numeric(edges$movement_weight)) else {
    rep(1, nrow(edges))
  }
  raw_weight[!is.finite(raw_weight) | raw_weight < 0] <- 0
  raw_sum <- .hee_group_sum(raw_weight, source_groups$index,
                            source_groups$n_groups)[source_groups$index]
  degree <- .hee_group_sum(rep.int(1, nrow(edges)), source_groups$index,
                           source_groups$n_groups)[source_groups$index]
  permeability <- switch(
    source_permeability,
    mean_weight = if (has_raw_weight) raw_sum / pmax(degree, 1) else rep(1, nrow(edges)),
    sum_capped = if (has_raw_weight) pmin(raw_sum, 1) else rep(1, nrow(edges)),
    none = rep(1, nrow(edges))
  )
  edges$source_permeability <- pmin(pmax(permeability, 0), 1)
  edges$source_total_rate <- base_rate * edges$source_permeability
  edges$movement_rate <- edges$source_total_rate *
    edges$source_direction_weight
  edges$source_rate_definition <- paste0(source_permeability,
                                         ";base_rate_per_Myr")
  edges
}

#' Calculate raw and standardised dispersal arrival pressure
#'
#' Computes \deqn{\Lambda^{arr}_{\ell,j,t}=\sum_i q_{\ell,i,t}K_{\ell,i\to j,t}}
#' and returns the optional bounded transformation
#' \deqn{M^{arr}=1-\exp(-\Lambda^{arr}).} This is the P2 output passed to
#' colonisation; it is not an environmental-suitability score.
#'
#' @param occupancy Lineage-cell source occupancy with a probability column.
#' @param movement_kernel Output from [hee_dispersal_edge_kernel()].
#' @param probability_col Occupancy probability column.
#' @param transform Either `"one_minus_exp"` or `"none"`.
#'
#' @return A lineage-cell-time table with `raw_arrival` and
#'   `arrival_pressure`.
#' @export
hee_dispersal_arrival_pressure <- function(occupancy,
                                           movement_kernel,
                                           probability_col = "q",
                                           transform = c("one_minus_exp", "none")) {
  transform <- match.arg(transform)
  q <- as.data.frame(occupancy)
  k <- as.data.frame(movement_kernel)
  .require_cols(q, c("lineage", "cell_id", probability_col), "occupancy")
  .require_cols(k, c("lineage", "from_cell_id", "to_cell_id", "time_ma",
                     "K_movement"), "movement_kernel")
  .hee_check_unique_keys(q, c("lineage", "cell_id"), "occupancy")
  .hee_check_unique_keys(k, c("lineage", "from_cell_id", "to_cell_id", "time_ma"),
                         "movement_kernel")
  q$.q <- .hee_hmscee_prob(q[[probability_col]])
  match_key <- paste(k$lineage, k$from_cell_id, sep = "\r")
  source_key <- paste(q$lineage, q$cell_id, sep = "\r")
  k$.q <- q$.q[match(match_key, source_key)]
  k$.q[is.na(k$.q)] <- 0
  k$.contribution <- k$.q * .hee_clip01(k$K_movement)
  out <- stats::aggregate(k$.contribution,
                          k[, c("lineage", "to_cell_id", "time_ma")], sum)
  names(out) <- c("lineage", "cell_id", "time_ma", "raw_arrival")
  out$arrival_pressure <- if (identical(transform, "one_minus_exp")) {
    1 - exp(-pmax(out$raw_arrival, 0))
  } else out$raw_arrival
  out
}

#' Sample a legacy particle-based approximation to a dispersal kernel
#'
#' This retained function is named [hee_dispersal_propagule_particles()] in the
#' preferred API. Its particles route a finite number of propagules through an
#' already defined movement kernel; they are not posterior particles or an
#' ancestral-location reconstruction.
#'
#' Converts a continuous edge kernel into a finite-propagule, particle-sampled
#' edge kernel. For every source cell, `particles_per_source` propagule paths
#' are sampled with probability proportional to `K_movement`. The total
#' outgoing kernel is retained, so `particle_kernel` is an unbiased Monte
#' Carlo approximation to `K_movement`, conditional on the supplied movement
#' graph. This is a stochastic approximation to a structural dispersal
#' scenario, not a posterior particle reconstruction and not an estimate of
#' observed gene flow.
#'
#' @param movement_kernel A table produced by
#'   [hee_dispersal_edge_kernel()] with `from_cell_id`, `to_cell_id`, and
#'   `K_movement` columns. `lineage` and `time_ma`, when supplied, define
#'   independent source groups.
#' @param particles_per_source Positive integer number of propagule paths
#'   sampled for each source group.
#' @param seed Optional integer seed for reproducible particle routing.
#' @return `movement_kernel` augmented with `particle_count`,
#'   `source_kernel_total`, `particle_fraction`, and `particle_kernel`.
#'   `particle_kernel` is a non-negative movement-rate multiplier and is not
#'   individually constrained to `[0, 1]`.
#' @export
#'
#' @examples
#' edges <- data.frame(
#'   from_cell_id = c("a", "a", "b"), to_cell_id = c("b", "c", "a"),
#'   K_movement = c(0.8, 0.2, 0.5)
#' )
#' hee_dispersal_particle_kernel(edges, particles_per_source = 20, seed = 1)
hee_dispersal_particle_kernel <- function(movement_kernel,
                                          particles_per_source = 32L,
                                          seed = NULL) {
  edges <- as.data.frame(movement_kernel)
  .require_cols(edges, c("from_cell_id", "to_cell_id", "K_movement"),
                "movement_kernel")
  particles_per_source <- suppressWarnings(as.integer(particles_per_source)[1L])
  if (!is.finite(particles_per_source) || particles_per_source < 1L) {
    stop("particles_per_source must be a positive integer.", call. = FALSE)
  }
  edges$K_movement <- suppressWarnings(as.numeric(edges$K_movement))
  edges$K_movement[!is.finite(edges$K_movement) | edges$K_movement < 0] <- 0
  group_cols <- c(intersect(c("lineage", "time_ma", "response_draw"),
                             names(edges)), "from_cell_id")
  .hee_check_unique_keys(edges, c(group_cols, "to_cell_id"),
                         "movement_kernel")

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

  group_id <- do.call(paste, c(edges[group_cols], sep = "\r"))
  group_index <- match(group_id, unique(group_id))
  group_size <- tabulate(group_index)
  n_group <- length(group_size)
  order_index <- order(group_index)
  slot_ordered <- sequence(group_size)
  edge_slot <- integer(nrow(edges))
  edge_slot[order_index] <- slot_ordered
  max_degree <- max(group_size)
  weight_matrix <- matrix(0, n_group, max_degree)
  weight_matrix[cbind(group_index, edge_slot)] <- edges$K_movement
  source_total <- rowSums(weight_matrix)
  edges$particle_count <- 0L
  edges$source_kernel_total <- source_total[group_index]
  edges$particle_fraction <- 0
  edges$particle_kernel <- 0
  active <- which(is.finite(source_total) & source_total > 0)
  if (length(active)) {
    # The connectivity graph has a small local degree (normally <= 8). Sampling
    # all source cells in a matrix avoids tens of thousands of `sample()` calls
    # while retaining independent multinomial propagule paths per source.
    probability_matrix <- weight_matrix[active, , drop = FALSE] /
      source_total[active]
    cdf <- t(apply(probability_matrix, 1L, cumsum))
    count_active <- matrix(0L, nrow(cdf), ncol(cdf))
    active_row <- seq_len(nrow(cdf))
    for (draw in seq_len(particles_per_source)) {
      u <- stats::runif(nrow(cdf))
      destination <- rowSums(cdf < u) + 1L
      destination[destination > ncol(cdf)] <- ncol(cdf)
      count_active[cbind(active_row, destination)] <-
        count_active[cbind(active_row, destination)] + 1L
    }
    count_matrix <- matrix(0L, n_group, max_degree)
    count_matrix[active, ] <- count_active
    edges$particle_count <- count_matrix[cbind(group_index, edge_slot)]
    edges$particle_fraction <- edges$particle_count / particles_per_source
    # Preserves each source total exactly, and each edge in expectation.
    edges$particle_kernel <- edges$source_kernel_total * edges$particle_fraction
  }
  edges$particle_sampling <- paste0("multinomial_", particles_per_source,
                                    "_per_source")
  edges$particle_seed <- if (is.null(seed)) NA_integer_ else seed
  edges
}
