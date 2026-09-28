# Landscape-explicit dispersal ------------------------------------------------

#' Build a movement landscape from a dynamic Earth grid
#'
#' This is the Earth-driver interface for the HmscEE dispersal process. It
#' separates a hard geographical state from a continuous movement landscape
#' weight. A cell that does not exist or is declared non-traversable has zero
#' movement weight. Continuous landscape weights describe the permeability of
#' the *same* movement landscape; they are not HMSC suitability, a second
#' colonisation probability, or a residual species association.
#'
#' For a valid target cell, the numerical integration measure used by
#' [hee_dispersal_spherical_kernel()] is
#' \deqn{m_{c,k}=a_{c,k}\,\omega_{c,k},}
#' where \eqn{a} is cell area and \eqn{\omega} is the declared landscape
#' permeability. Including area prevents a latitude-dependent sampling bias on
#' a longitude-latitude grid.
#'
#' @references Arias, J. S. (2024). Phylogenetic biogeography inference using
#'   dynamic paleogeography models and explicit geographic ranges. *Systematic
#'   Biology*, 73, 995-1014. https://doi.org/10.1093/sysbio/syae051
#'
#' @param cells A cell-time data frame. It must contain unique `cell_id` values
#'   within the supplied slice.
#' @param cell_col Column identifying cells.
#' @param hard_mask_col Optional 0/1 or logical column declaring whether the
#'   ecological stage exists and is traversable. Missing values are treated as
#'   non-traversable.
#' @param movement_weight_col Optional continuous `[0, 1]` movement-landscape
#'   permeability. It must already have a single documented biological meaning
#'   (for example, a habitat-matrix permeability), rather than a product of
#'   several uncalibrated connectivity indices.
#' @param landscape_class_col Optional categorical landscape column.
#' @param landscape_class_weights Optional named numeric vector in `[0, 1]`.
#'   Values are looked up from `landscape_class_col` and multiplied with the
#'   continuous movement weight when both are supplied.
#' @param cell_area_col Optional positive cell-area column in km2. If omitted,
#'   equal cell measures are assumed and reported explicitly in metadata.
#' @return A `hee_dispersal_landscape` object with a `cells` table containing
#'   `movement_allowed`, `movement_landscape_weight`, `cell_measure_weight`, and
#'   `destination_measure`.
#' @export
#'
#' @examples
#' cells <- data.frame(
#'   cell_id = c("a", "b", "o"), land = c(1, 1, 0),
#'   habitat_permeability = c(1, 0.4, 1), cell_area_km2 = c(1, 2, 1)
#' )
#' hee_dispersal_landscape_weights(
#'   cells, hard_mask_col = "land",
#'   movement_weight_col = "habitat_permeability",
#'   cell_area_col = "cell_area_km2"
#' )
hee_dispersal_landscape_weights <- function(cells,
                                            cell_col = "cell_id",
                                            hard_mask_col = NULL,
                                            movement_weight_col = NULL,
                                            landscape_class_col = NULL,
                                            landscape_class_weights = NULL,
                                            cell_area_col = NULL) {
  x <- as.data.frame(cells, stringsAsFactors = FALSE)
  .require_cols(x, cell_col, "cells")
  x[[cell_col]] <- as.character(x[[cell_col]])
  if (any(!nzchar(x[[cell_col]])) || anyDuplicated(x[[cell_col]])) {
    stop("cells must contain unique, non-blank cell IDs within one time slice.",
         call. = FALSE)
  }
  if (!is.null(hard_mask_col)) .require_cols(x, hard_mask_col, "cells")
  if (!is.null(movement_weight_col)) {
    .require_cols(x, movement_weight_col, "cells")
  }
  if (!is.null(landscape_class_col)) {
    .require_cols(x, landscape_class_col, "cells")
    if (is.null(landscape_class_weights) || is.null(names(landscape_class_weights))) {
      stop("landscape_class_weights must be a named numeric vector when landscape_class_col is supplied.",
           call. = FALSE)
    }
  } else if (!is.null(landscape_class_weights)) {
    stop("landscape_class_col is required when landscape_class_weights is supplied.",
         call. = FALSE)
  }
  if (!is.null(cell_area_col)) .require_cols(x, cell_area_col, "cells")

  allowed <- rep(TRUE, nrow(x))
  if (!is.null(hard_mask_col)) {
    hard <- x[[hard_mask_col]]
    if (is.logical(hard)) hard <- as.numeric(hard)
    hard <- suppressWarnings(as.numeric(hard))
    allowed <- is.finite(hard) & hard > 0
  }

  movement_weight <- rep(1, nrow(x))
  if (!is.null(movement_weight_col)) {
    movement_weight <- suppressWarnings(as.numeric(x[[movement_weight_col]]))
    if (any(!is.finite(movement_weight[allowed]) |
            movement_weight[allowed] < 0 | movement_weight[allowed] > 1)) {
      stop("movement_weight_col must be finite and within [0, 1] for traversable cells.",
           call. = FALSE)
    }
    movement_weight[!is.finite(movement_weight)] <- 0
  }
  if (!is.null(landscape_class_col)) {
    class_value <- as.character(x[[landscape_class_col]])
    class_weight <- unname(landscape_class_weights[class_value])
    if (any(is.na(class_weight[allowed]))) {
      missing <- unique(class_value[allowed & is.na(class_weight)])
      stop("landscape_class_weights is missing traversable class(es): ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    class_weight <- suppressWarnings(as.numeric(class_weight))
    if (any(!is.finite(class_weight[allowed]) |
            class_weight[allowed] < 0 | class_weight[allowed] > 1)) {
      stop("landscape_class_weights must be finite and within [0, 1].",
           call. = FALSE)
    }
    class_weight[!is.finite(class_weight)] <- 0
    movement_weight <- movement_weight * class_weight
  }

  if (is.null(cell_area_col)) {
    area <- rep(1, nrow(x))
    area_source <- "implicit_equal_cell_measure"
  } else {
    area <- suppressWarnings(as.numeric(x[[cell_area_col]]))
    if (any(!is.finite(area[allowed]) | area[allowed] <= 0)) {
      stop("cell_area_col must be finite and positive for traversable cells.",
           call. = FALSE)
    }
    area[!is.finite(area) | area <= 0] <- NA_real_
    area_source <- cell_area_col
  }
  allowed <- allowed & is.finite(area) & area > 0 & movement_weight > 0
  reference_area <- if (any(allowed)) stats::median(area[allowed]) else NA_real_
  measure <- rep(0, nrow(x))
  if (any(allowed)) measure[allowed] <- area[allowed] / reference_area

  x$movement_allowed <- as.integer(allowed)
  x$movement_landscape_weight <- ifelse(allowed, movement_weight, 0)
  x$cell_measure_weight <- measure
  x$destination_measure <- x$movement_landscape_weight * x$cell_measure_weight
  diagnostics <- data.frame(
    n_cells = nrow(x),
    n_traversable_cells = sum(allowed),
    traversable_fraction = mean(allowed),
    reference_area_km2 = reference_area,
    area_source = area_source,
    stringsAsFactors = FALSE
  )
  structure(
    list(
      cells = x,
      diagnostics = diagnostics,
      metadata = list(
        hard_mask_col = hard_mask_col,
        movement_weight_col = movement_weight_col,
        landscape_class_col = landscape_class_col,
        cell_area_col = cell_area_col,
        interpretation = paste(
          "Dynamic-Earth movement landscape. Hard mask is a geographic state;",
          "continuous weights are one declared movement-permeability surface,",
          "not environmental suitability or a colonisation probability."
        )
      )
    ),
    class = "hee_dispersal_landscape"
  )
}

.hee_dispersal_slice <- function(cells, cell_col, lon_col, lat_col,
                                 time_col, time_ma) {
  x <- if (inherits(cells, "hee_dispersal_landscape")) cells$cells else {
    as.data.frame(cells, stringsAsFactors = FALSE)
  }
  .require_cols(x, c(cell_col, lon_col, lat_col), "cells")
  if (time_col %in% names(x)) {
    time_values <- sort(unique(suppressWarnings(as.numeric(x[[time_col]]))))
    time_values <- time_values[is.finite(time_values)]
    if (is.null(time_ma)) {
      if (length(time_values) != 1L) {
        stop("cells contains multiple time values; supply time_ma and one slice at a time.",
             call. = FALSE)
      }
      time_ma <- time_values[[1L]]
    }
    time_ma <- suppressWarnings(as.numeric(time_ma)[1L])
    if (!is.finite(time_ma)) stop("time_ma must be finite.", call. = FALSE)
    x <- x[is.finite(as.numeric(x[[time_col]])) &
             abs(as.numeric(x[[time_col]]) - time_ma) < 1e-10, , drop = FALSE]
  } else {
    time_ma <- suppressWarnings(as.numeric(time_ma)[1L])
    if (!is.finite(time_ma)) {
      stop("cells has no time column; supply finite time_ma.", call. = FALSE)
    }
  }
  if (!nrow(x)) stop("No cells were found for the requested time_ma.", call. = FALSE)
  x[[cell_col]] <- as.character(x[[cell_col]])
  if (any(!nzchar(x[[cell_col]])) || anyDuplicated(x[[cell_col]])) {
    stop("The selected cell-time slice must have unique non-blank cell IDs.",
         call. = FALSE)
  }
  x[[lon_col]] <- suppressWarnings(as.numeric(x[[lon_col]]))
  x[[lat_col]] <- suppressWarnings(as.numeric(x[[lat_col]]))
  if (any(!is.finite(x[[lon_col]]) | !is.finite(x[[lat_col]]) |
          abs(x[[lat_col]]) > 90)) {
    stop("cells must contain finite longitude and latitude values; latitude must be in [-90, 90].",
         call. = FALSE)
  }
  list(cells = x, time_ma = time_ma)
}

.hee_dispersal_candidate_edges <- function(cells, cell_col, lon_col, lat_col,
                                           max_distance_km, include_self) {
  if (nrow(cells) > 2000L) {
    stop("For grids with more than 2000 cells, supply sparse local edges from hee_dispersal_connectivity_graph(); all-pairs candidate generation is disabled.",
         call. = FALSE)
  }
  max_distance_km <- suppressWarnings(as.numeric(max_distance_km)[1L])
  if (!is.finite(max_distance_km) || max_distance_km <= 0) {
    stop("max_distance_km must be a positive finite value when edges is NULL.",
         call. = FALSE)
  }
  distance <- .hee_great_circle_matrix(
    cells[[lon_col]], cells[[lat_col]], cells[[lon_col]], cells[[lat_col]]
  )
  keep <- is.finite(distance) & distance <= max_distance_km
  if (!include_self) diag(keep) <- FALSE
  ij <- which(keep, arr.ind = TRUE)
  data.frame(
    from_cell_id = as.character(cells[[cell_col]][ij[, 1L]]),
    to_cell_id = as.character(cells[[cell_col]][ij[, 2L]]),
    distance_km = as.numeric(distance[keep]),
    stringsAsFactors = FALSE
  )
}

#' Construct a landscape-weighted spherical dispersal kernel
#'
#' Implements the directional component of a landscape-explicit spherical
#' normal diffusion process. For source cell \eqn{i}, target cell \eqn{j}, and
#' a lineage \eqn{l},
#' \deqn{K_{i\to j,l,k}=\frac{a_{j,k}\omega_{j,k}
#' \exp[-\tilde d_{ij,k}^2/(2\sigma^2_{l,k}\Delta t_k)]}
#' {\sum_z a_{z,k}\omega_{z,k}
#' \exp[-\tilde d_{iz,k}^2/(2\sigma^2_{l,k}\Delta t_k)]}.}
#' Here \eqn{\tilde d} is great-circle distance in radians by default.  When
#' `edges` has a documented structural path cost (by default
#' `effective_cost_km`), \eqn{\tilde d} is its great-circle-equivalent
#' effective distance. This is the normalized landscape-weighted
#' spherical-normal transition used for an explicit grid; it is equivalent to
#' the precision parameterisation
#' \eqn{\lambda=1/\sigma^2}. A separate source-emigration rate turns the
#' directional kernel into a movement-rate matrix suitable for
#' [hee_projection_global_grid_step()].
#'
#' The function intentionally operates on one time slice and sparse candidate
#' edges. For a global 1-degree grid, construct local edges first with
#' [hee_dispersal_connectivity_graph()] rather than creating an all-to-all
#' matrix. Long-distance dispersal must be added as an explicit, documented
#' edge scenario; it is never generated silently.
#' If an edge cost already includes topographic resistance, do not also put the
#' same topographic quantity into `movement_weight_col`; that would penalise the
#' same Earth feature twice.
#'
#' The spherical diffusion form follows Arias (2024). The accompanying sparse
#' time graph separates horizontal active movement from vertical plate-carriage
#' edges, following the landscape-explicit principle used by
#' Flannery-Sutherland et al. (2025). It does not convert a least-cost route
#' diagnostic into a second biological dispersal probability.
#'
#' @references Arias, J. S. (2024). Phylogenetic biogeography inference using
#'   dynamic paleogeography models and explicit geographic ranges. *Systematic
#'   Biology*, 73, 995-1014. https://doi.org/10.1093/sysbio/syae051;
#'   Flannery-Sutherland, J. T., Elsler, A., Farnsworth, A., Lunt, D. J. &
#'   Benton, M. J. (2025). Landscape-explicit phylogeography illuminates the
#'   ecographic radiation of early archosauromorph reptiles. *Nature Ecology &
#'   Evolution*, 9, 1138-1152. https://doi.org/10.1038/s41559-025-02739-y
#'
#' @param cells A `hee_dispersal_landscape` object or cell-time table with
#'   `cell_id`, `lon`, and `lat` columns. If plain cells are supplied, missing
#'   movement fields default to a traversable, unit-weight landscape.
#' @param edges Optional sparse candidate edges with `from_cell_id` and
#'   `to_cell_id`. A `distance_km` column is optional and will be calculated
#'   from coordinates when absent. If `edges` also contains `time_col`, only
#'   edges from the selected historical slice are retained.
#' @param edge_cost_col Optional non-negative edge-cost column, defaulting to
#'   `"effective_cost_km"`. When that column is present, it must be at least
#'   `distance_km` and is used as a great-circle-equivalent effective movement
#'   distance. Set to `NULL` for pure spherical diffusion.
#' @param delta_t_myr Positive duration in Myr.
#' @param diffusion_variance_rad2_per_myr Diffusion variance in radians2 per
#'   Myr. Supply exactly one of this or `diffusion_precision_myr_per_rad2`.
#' @param diffusion_precision_myr_per_rad2 Arias-style precision
#'   \eqn{\lambda=1/\sigma^2}, in Myr per radians2.
#' @param source_emigration_rate_per_myr Non-negative source-cell propagule
#'   emigration rate. It is not estimated by HMSC from one modern snapshot.
#' @param lineage Optional lineage label recorded in the output.
#' @param time_ma Time slice in Ma. Larger Ma values are older. It is inferred
#'   when `cells` contains exactly one value.
#' @param cell_col,lon_col,lat_col,time_col Column names in `cells`.
#' @param max_distance_km Maximum candidate distance when `edges = NULL`; only
#'   intended for small examples.
#' @param include_self Whether self transitions are allowed in the directional
#'   kernel. Usually `FALSE` because persistence is modelled separately.
#' @param spatial_domain `"global_explicit"` permits movement across the explicit
#'   grid. `"within_region"` retains only same-region edges and is the required
#'   setting when a BioGeoBEARS history already owns cross-region dispersal.
#' @param region_col Region column used by `within_region`.
#' @return A `hee_dispersal_spherical_kernel` object with directed `edges`.
#'   `directional_weight` sums to one by source, while `movement_rate_per_myr`
#'   sums to `source_emigration_rate_per_myr` by source.
#' @export
#'
#' @examples
#' cells <- data.frame(
#'   cell_id = c("a", "b", "c"), lon = c(179, -179, 160), lat = 0,
#'   movement_allowed = 1, movement_landscape_weight = c(1, 1, 0.5),
#'   cell_area_km2 = 1, time_ma = 10
#' )
#' k <- hee_dispersal_spherical_kernel(
#'   cells, delta_t_myr = 1, diffusion_variance_rad2_per_myr = 0.01,
#'   max_distance_km = 3000
#' )
#' k$edges
hee_dispersal_spherical_kernel <- function(
    cells,
    edges = NULL,
    delta_t_myr,
    diffusion_variance_rad2_per_myr = NULL,
    diffusion_precision_myr_per_rad2 = NULL,
    source_emigration_rate_per_myr = 1,
    lineage = "global",
    time_ma = NULL,
    cell_col = "cell_id",
    lon_col = "lon",
    lat_col = "lat",
    time_col = "time_ma",
    max_distance_km = NULL,
    include_self = FALSE,
    spatial_domain = c("global_explicit", "within_region"),
    region_col = "region",
    edge_cost_col = "effective_cost_km") {
  spatial_domain <- match.arg(spatial_domain)
  slice <- .hee_dispersal_slice(cells, cell_col, lon_col, lat_col,
                                time_col, time_ma)
  x <- slice$cells
  time_ma <- slice$time_ma
  delta_t_myr <- suppressWarnings(as.numeric(delta_t_myr)[1L])
  if (!is.finite(delta_t_myr) || delta_t_myr <= 0) {
    stop("delta_t_myr must be a positive finite duration in Myr.", call. = FALSE)
  }
  variance <- suppressWarnings(as.numeric(diffusion_variance_rad2_per_myr)[1L])
  precision <- suppressWarnings(as.numeric(diffusion_precision_myr_per_rad2)[1L])
  have_variance <- is.finite(variance)
  have_precision <- is.finite(precision)
  if (have_variance == have_precision) {
    stop("Supply exactly one of diffusion_variance_rad2_per_myr or diffusion_precision_myr_per_rad2.",
         call. = FALSE)
  }
  if (have_precision) {
    if (precision <= 0) {
      stop("diffusion_precision_myr_per_rad2 must be positive.", call. = FALSE)
    }
    variance <- 1 / precision
  }
  if (!is.finite(variance) || variance <= 0) {
    stop("diffusion variance must be positive and finite.", call. = FALSE)
  }
  source_emigration_rate_per_myr <- suppressWarnings(
    as.numeric(source_emigration_rate_per_myr)[1L]
  )
  if (!is.finite(source_emigration_rate_per_myr) ||
      source_emigration_rate_per_myr < 0) {
    stop("source_emigration_rate_per_myr must be finite and non-negative.",
         call. = FALSE)
  }

  if (!"movement_allowed" %in% names(x)) x$movement_allowed <- 1L
  if (!"movement_landscape_weight" %in% names(x)) {
    x$movement_landscape_weight <- 1
  }
  if (!"cell_measure_weight" %in% names(x)) x$cell_measure_weight <- 1
  x$movement_allowed <- as.integer(as.numeric(x$movement_allowed) > 0)
  x$movement_landscape_weight <- suppressWarnings(
    as.numeric(x$movement_landscape_weight)
  )
  x$cell_measure_weight <- suppressWarnings(as.numeric(x$cell_measure_weight))
  valid_weight <- is.finite(x$movement_landscape_weight) &
    x$movement_landscape_weight >= 0 & x$movement_landscape_weight <= 1 &
    is.finite(x$cell_measure_weight) & x$cell_measure_weight > 0
  x$movement_allowed[!valid_weight] <- 0L
  x$destination_measure <- ifelse(
    x$movement_allowed > 0,
    x$movement_landscape_weight * x$cell_measure_weight,
    0
  )
  if (spatial_domain == "within_region") .require_cols(x, region_col, "cells")

  candidate <- if (is.null(edges)) {
    .hee_dispersal_candidate_edges(x, cell_col, lon_col, lat_col,
                                   max_distance_km, include_self)
  } else {
    as.data.frame(edges, stringsAsFactors = FALSE)
  }
  .require_cols(candidate, c("from_cell_id", "to_cell_id"), "edges")
  if (!is.null(edges) && time_col %in% names(candidate)) {
    edge_time <- suppressWarnings(as.numeric(candidate[[time_col]]))
    candidate <- candidate[is.finite(edge_time) &
                             abs(edge_time - time_ma) < 1e-10, , drop = FALSE]
  }
  candidate$from_cell_id <- as.character(candidate$from_cell_id)
  candidate$to_cell_id <- as.character(candidate$to_cell_id)
  n_candidate_edges <- nrow(candidate)
  if (any(!candidate$from_cell_id %in% x[[cell_col]]) ||
      any(!candidate$to_cell_id %in% x[[cell_col]])) {
    stop("edges contains cell IDs absent from the selected cells slice.",
         call. = FALSE)
  }
  if (!include_self) {
    candidate <- candidate[candidate$from_cell_id != candidate$to_cell_id, , drop = FALSE]
  }
  from_i <- match(candidate$from_cell_id, x[[cell_col]])
  to_i <- match(candidate$to_cell_id, x[[cell_col]])
  keep <- x$movement_allowed[from_i] > 0 & x$movement_allowed[to_i] > 0
  if (spatial_domain == "within_region") {
    keep <- keep & as.character(x[[region_col]][from_i]) ==
      as.character(x[[region_col]][to_i])
  }
  candidate <- candidate[keep, , drop = FALSE]
  from_i <- from_i[keep]
  to_i <- to_i[keep]
  if (!nrow(candidate)) {
    empty_edges <- data.frame(
      lineage = character(),
      time_ma = numeric(),
      from_cell_id = character(),
      to_cell_id = character(),
      distance_km = numeric(),
      distance_rad = numeric(),
      effective_movement_distance_km = numeric(),
      effective_movement_distance_rad = numeric(),
      resistance_distance_multiplier = numeric(),
      target_landscape_weight = numeric(),
      target_cell_measure = numeric(),
      spherical_shape = numeric(),
      raw_directional_weight = numeric(),
      directional_weight = numeric(),
      movement_rate_per_myr = numeric(),
      stringsAsFactors = FALSE
    )
    return(structure(
      list(
        edges = empty_edges,
        diagnostics = data.frame(n_candidate_edges = n_candidate_edges, n_active_edges = 0L,
                                 n_active_sources = 0L, stringsAsFactors = FALSE),
        metadata = list(time_ma = time_ma, lineage = lineage,
                        spatial_domain = spatial_domain,
                        diffusion_variance_rad2_per_myr = variance,
                        source_emigration_rate_per_myr = source_emigration_rate_per_myr)
      ), class = "hee_dispersal_spherical_kernel"
    ))
  }
  if ("distance_km" %in% names(candidate)) {
    distance_km <- suppressWarnings(as.numeric(candidate$distance_km))
  } else {
    distance_km <- hee_great_circle_distance_km(
      x[[lon_col]][from_i], x[[lat_col]][from_i],
      x[[lon_col]][to_i], x[[lat_col]][to_i]
    )
  }
  if (any(!is.finite(distance_km) | distance_km < 0)) {
    stop("edges$distance_km must be finite and non-negative when supplied.",
         call. = FALSE)
  }
  use_edge_cost <- is.character(edge_cost_col) && length(edge_cost_col) == 1L &&
    !is.na(edge_cost_col) && nzchar(edge_cost_col) &&
    edge_cost_col %in% names(candidate)
  effective_distance_km <- distance_km
  if (use_edge_cost) {
    effective_distance_km <- suppressWarnings(as.numeric(candidate[[edge_cost_col]]))
    if (any(!is.finite(effective_distance_km) | effective_distance_km < 0)) {
      stop("edges$", edge_cost_col,
           " must be finite and non-negative when used as an edge cost.",
           call. = FALSE)
    }
    if (any(effective_distance_km < distance_km - 1e-8)) {
      stop("edges$", edge_cost_col,
           " must be greater than or equal to distance_km. It is a resistance-expanded movement distance, not a shortcut distance.",
           call. = FALSE)
    }
  }
  distance_rad <- distance_km / 6371.0088
  effective_distance_rad <- effective_distance_km / 6371.0088
  resistance_multiplier <- ifelse(distance_km > 0,
                                  effective_distance_km / distance_km, 1)
  interval_variance <- variance * delta_t_myr
  spherical_shape <- exp(-(effective_distance_rad^2) / (2 * interval_variance))
  raw_weight <- spherical_shape * x$destination_measure[to_i]
  group <- candidate$from_cell_id
  denominator <- stats::ave(raw_weight, group, FUN = sum)
  directional <- ifelse(denominator > 0, raw_weight / denominator, 0)
  out <- data.frame(
    lineage = as.character(lineage),
    time_ma = time_ma,
    from_cell_id = candidate$from_cell_id,
    to_cell_id = candidate$to_cell_id,
    distance_km = distance_km,
    distance_rad = distance_rad,
    effective_movement_distance_km = effective_distance_km,
    effective_movement_distance_rad = effective_distance_rad,
    resistance_distance_multiplier = resistance_multiplier,
    target_landscape_weight = x$movement_landscape_weight[to_i],
    target_cell_measure = x$cell_measure_weight[to_i],
    spherical_shape = spherical_shape,
    raw_directional_weight = raw_weight,
    directional_weight = directional,
    movement_rate_per_myr = source_emigration_rate_per_myr * directional,
    stringsAsFactors = FALSE
  )
  diagnostics <- data.frame(
    n_candidate_edges = n_candidate_edges,
    n_active_edges = nrow(out),
    n_active_sources = length(unique(out$from_cell_id)),
    maximum_directional_sum_error = max(abs(
      tapply(out$directional_weight, out$from_cell_id, sum) - 1
    )),
    stringsAsFactors = FALSE
  )
  structure(
    list(
      edges = out,
      diagnostics = diagnostics,
      metadata = list(
        time_ma = time_ma,
        lineage = as.character(lineage),
        delta_t_myr = delta_t_myr,
        diffusion_variance_rad2_per_myr = variance,
        diffusion_precision_myr_per_rad2 = 1 / variance,
        source_emigration_rate_per_myr = source_emigration_rate_per_myr,
        spatial_domain = spatial_domain,
        region_col = if (spatial_domain == "within_region") region_col else NA_character_,
        edge_cost_col = if (use_edge_cost) edge_cost_col else NA_character_,
        interpretation = paste(
          "Active P2 directional movement allocation. Plate transport is not",
          "included here; apply hee_plate_carry_occupancy() before active movement.",
          "An optional resistance-expanded edge distance is used once, not also",
          "as a duplicate destination landscape weight."
        )
      )
    ),
    class = "hee_dispersal_spherical_kernel"
  )
}

#' Convert a sparse dispersal-edge table to a target-by-source matrix
#'
#' @param movement_kernel A `hee_dispersal_spherical_kernel` object or edge
#'   table containing `from_cell_id`, `to_cell_id`, and `value_col`.
#' @param cell_order Optional stable cell ordering. Defaults to all cells present
#'   in the edge table.
#' @param value_col Edge value to place in the matrix. The default
#'   `"movement_rate_per_myr"` is compatible with
#'   [hee_projection_global_grid_step()].
#' @param sparse Return a `Matrix` sparse matrix when the optional Matrix package
#'   is installed.
#' @return A target-by-source transition or movement-rate matrix.
#' @export
#'
#' @examples
#' k <- data.frame(from_cell_id = "a", to_cell_id = "b",
#'                 movement_rate_per_myr = 0.2)
#' hee_dispersal_transition_matrix(k, c("a", "b"))
hee_dispersal_transition_matrix <- function(movement_kernel,
                                             cell_order = NULL,
                                             value_col = "movement_rate_per_myr",
                                             sparse = TRUE) {
  edges <- if (inherits(movement_kernel, "hee_dispersal_spherical_kernel")) {
    movement_kernel$edges
  } else as.data.frame(movement_kernel, stringsAsFactors = FALSE)
  required <- c("from_cell_id", "to_cell_id", value_col)
  .require_cols(edges, required, "movement_kernel")
  if (!nrow(edges)) {
    if (is.null(cell_order)) {
      stop("cell_order is required when a movement kernel has no active edges.",
           call. = FALSE)
    }
    cell_order <- as.character(cell_order)
    if (any(!nzchar(cell_order)) || anyDuplicated(cell_order)) {
      stop("cell_order must contain unique, non-blank cell IDs.", call. = FALSE)
    }
    if (isTRUE(sparse) && requireNamespace("Matrix", quietly = TRUE)) {
      return(Matrix::Matrix(
        0, nrow = length(cell_order), ncol = length(cell_order), sparse = TRUE,
        dimnames = list(cell_order, cell_order)
      ))
    }
    return(matrix(0, nrow = length(cell_order), ncol = length(cell_order),
                  dimnames = list(cell_order, cell_order)))
  }
  edges$from_cell_id <- as.character(edges$from_cell_id)
  edges$to_cell_id <- as.character(edges$to_cell_id)
  value <- suppressWarnings(as.numeric(edges[[value_col]]))
  if (any(!is.finite(value) | value < 0)) {
    stop("movement-kernel values must be finite and non-negative.", call. = FALSE)
  }
  if (is.null(cell_order)) {
    cell_order <- unique(c(edges$from_cell_id, edges$to_cell_id))
  }
  cell_order <- as.character(cell_order)
  if (any(!nzchar(cell_order)) || anyDuplicated(cell_order)) {
    stop("cell_order must contain unique, non-blank cell IDs.", call. = FALSE)
  }
  from <- match(edges$from_cell_id, cell_order)
  to <- match(edges$to_cell_id, cell_order)
  if (anyNA(from) || anyNA(to)) {
    stop("movement kernel contains cells absent from cell_order.", call. = FALSE)
  }
  if (isTRUE(sparse) && requireNamespace("Matrix", quietly = TRUE)) {
    return(Matrix::sparseMatrix(
      i = to, j = from, x = value,
      dims = c(length(cell_order), length(cell_order)),
      dimnames = list(cell_order, cell_order)
    ))
  }
  out <- matrix(0, nrow = length(cell_order), ncol = length(cell_order),
                dimnames = list(cell_order, cell_order))
  for (ii in seq_len(nrow(edges))) out[to[[ii]], from[[ii]]] <-
    out[to[[ii]], from[[ii]]] + value[[ii]]
  out
}

.hee_dispersal_node_id <- function(cell_id, time_ma) {
  paste0(as.character(cell_id), "@", formatC(as.numeric(time_ma), digits = 12,
                                                format = "fg", flag = "#"))
}

#' Build a time-ordered landscape-explicit dispersal graph
#'
#' Combines horizontal active-movement edges with vertical plate-carriage edges.
#' Horizontal edges are the P2 dispersal substrate. Vertical edges are a
#' Dynamic-Earth coordinate/land-unit transfer between two palaeogeographic
#' stages and have zero path cost; they are explicitly flagged as
#' `"plate_carriage"` and are never interpreted as organismal dispersal.
#'
#' This follows a time-ordered graph representation: local geometry can vary by
#' time slice, while plate transport links homologous land units between stages.
#' It is a sparse graph for diagnostics and constrained historical inference,
#' not an instruction to compute all-pairs least-cost paths for every occupancy
#' update.
#'
#' @param horizontal_edges A data frame with `from_cell_id`, `to_cell_id`, a
#'   time column, and a non-negative cost column. It may be output from
#'   [hee_dispersal_connectivity_graph()].
#' @param plate_transport Optional [hee_plate_grid_transport()] object or table
#'   containing `source_cell_id`, `target_cell_id`, and `target_weight`.
#' @param horizontal_time_col,cost_col Column names in `horizontal_edges`.
#' @param time_from_col,time_to_col Transport time columns when they are present.
#' @param time_from_ma,time_to_ma Required for a single transport object lacking
#'   explicit time columns. Ma decreases toward the present, so each carriage
#'   edge must have `time_from_ma > time_to_ma`.
#' @param target_cells Optional allowed target-cell IDs. Transport links ending
#'   outside this set are removed instead of inventing land.
#' @return A `hee_dispersal_spacetime_graph` object.
#' @export
#'
#' @examples
#' horizontal <- data.frame(
#'   from_cell_id = c("a", "b"), to_cell_id = c("b", "a"),
#'   time_ma = 10, effective_cost_km = 100
#' )
#' transport <- data.frame(
#'   source_cell_id = "a", target_cell_id = "a2", target_weight = 1
#' )
#' hee_dispersal_spacetime_graph(
#'   horizontal, transport, time_from_ma = 10, time_to_ma = 5
#' )
hee_dispersal_spacetime_graph <- function(
    horizontal_edges,
    plate_transport = NULL,
    horizontal_time_col = "time_ma",
    cost_col = NULL,
    time_from_col = "time_from_ma",
    time_to_col = "time_to_ma",
    time_from_ma = NULL,
    time_to_ma = NULL,
    target_cells = NULL) {
  horizontal <- as.data.frame(horizontal_edges, stringsAsFactors = FALSE)
  .require_cols(horizontal, c("from_cell_id", "to_cell_id", horizontal_time_col),
                "horizontal_edges")
  if (is.null(cost_col)) {
    candidates <- c("effective_cost_km", "path_cost", "distance_km")
    cost_col <- candidates[candidates %in% names(horizontal)][1L]
  }
  if (is.na(cost_col) || is.null(cost_col)) {
    stop("horizontal_edges needs effective_cost_km, path_cost, distance_km, or an explicit cost_col.",
         call. = FALSE)
  }
  .require_cols(horizontal, cost_col, "horizontal_edges")
  horizontal$from_cell_id <- as.character(horizontal$from_cell_id)
  horizontal$to_cell_id <- as.character(horizontal$to_cell_id)
  horizontal_time <- suppressWarnings(as.numeric(horizontal[[horizontal_time_col]]))
  horizontal_cost <- suppressWarnings(as.numeric(horizontal[[cost_col]]))
  if (any(!is.finite(horizontal_time) | !is.finite(horizontal_cost) |
          horizontal_cost < 0)) {
    stop("horizontal time and cost values must be finite; costs must be non-negative.",
         call. = FALSE)
  }
  hdist <- if ("distance_km" %in% names(horizontal)) {
    suppressWarnings(as.numeric(horizontal$distance_km))
  } else horizontal_cost
  hdist[!is.finite(hdist)] <- NA_real_
  h <- data.frame(
    from_cell_id = horizontal$from_cell_id,
    to_cell_id = horizontal$to_cell_id,
    time_from_ma = horizontal_time,
    time_to_ma = horizontal_time,
    edge_kind = "active_dispersal",
    path_cost = horizontal_cost,
    geographic_distance_km = hdist,
    transport_weight = NA_real_,
    active_dispersal = TRUE,
    stringsAsFactors = FALSE
  )
  h$from_node <- .hee_dispersal_node_id(h$from_cell_id, h$time_from_ma)
  h$to_node <- .hee_dispersal_node_id(h$to_cell_id, h$time_to_ma)

  vertical <- h[FALSE, , drop = FALSE]
  if (!is.null(plate_transport)) {
    tr_object <- inherits(plate_transport, "hee_plate_grid_transport")
    tr <- if (tr_object) plate_transport$transport else {
      as.data.frame(plate_transport, stringsAsFactors = FALSE)
    }
    source_col <- if (tr_object) plate_transport$source_col else "source_cell_id"
    target_col <- if (tr_object) plate_transport$target_col else "target_cell_id"
    weight_col <- if (tr_object) plate_transport$weight_col else "target_weight"
    .require_cols(tr, c(source_col, target_col, weight_col), "plate_transport")
    if (all(c(time_from_col, time_to_col) %in% names(tr))) {
      from_time <- suppressWarnings(as.numeric(tr[[time_from_col]]))
      to_time <- suppressWarnings(as.numeric(tr[[time_to_col]]))
    } else {
      from_time <- rep(suppressWarnings(as.numeric(time_from_ma)[1L]), nrow(tr))
      to_time <- rep(suppressWarnings(as.numeric(time_to_ma)[1L]), nrow(tr))
    }
    weight <- suppressWarnings(as.numeric(tr[[weight_col]]))
    if (any(!is.finite(from_time) | !is.finite(to_time) |
            !is.finite(weight) | weight < 0)) {
      stop("plate transport requires finite times and non-negative weights.",
           call. = FALSE)
    }
    if (any(from_time <= to_time)) {
      stop("plate transport must be directed from older to younger stages: time_from_ma > time_to_ma.",
           call. = FALSE)
    }
    vertical <- data.frame(
      from_cell_id = as.character(tr[[source_col]]),
      to_cell_id = as.character(tr[[target_col]]),
      time_from_ma = from_time,
      time_to_ma = to_time,
      edge_kind = "plate_carriage",
      path_cost = 0,
      geographic_distance_km = 0,
      transport_weight = weight,
      active_dispersal = FALSE,
      stringsAsFactors = FALSE
    )
    if (!is.null(target_cells)) {
      vertical <- vertical[vertical$to_cell_id %in% as.character(target_cells), , drop = FALSE]
    }
    vertical$from_node <- .hee_dispersal_node_id(
      vertical$from_cell_id, vertical$time_from_ma
    )
    vertical$to_node <- .hee_dispersal_node_id(
      vertical$to_cell_id, vertical$time_to_ma
    )
  }
  edges <- rbind(h, vertical)
  if (any(!is.finite(edges$path_cost) | edges$path_cost < 0)) {
    stop("All graph path costs must be finite and non-negative.", call. = FALSE)
  }
  diagnostics <- data.frame(
    n_horizontal_edges = nrow(h),
    n_plate_carriage_edges = nrow(vertical),
    n_nodes = length(unique(c(edges$from_node, edges$to_node))),
    n_time_slices = length(unique(c(edges$time_from_ma, edges$time_to_ma))),
    stringsAsFactors = FALSE
  )
  structure(
    list(
      edges = edges,
      diagnostics = diagnostics,
      metadata = list(
        cost_col = cost_col,
        interpretation = paste(
          "Horizontal edges represent active dispersal substrate; zero-cost",
          "vertical edges represent palaeogeographic carriage before dispersal."
        )
      )
    ),
    class = "hee_dispersal_spacetime_graph"
  )
}

#' Find one least-cost path on a landscape-explicit time graph
#'
#' Uses a base-R Dijkstra implementation on a sparse,
#' [hee_dispersal_spacetime_graph()] object. It is a diagnostic or
#' phylogeographic-reconstruction tool. It should not be used to construct an
#' all-pairs movement matrix for every lineage and time step.
#'
#' @param graph A `hee_dispersal_spacetime_graph` object.
#' @param origin_cell_id,destination_cell_id Origin and destination cell IDs.
#' @param origin_time_ma,destination_time_ma Corresponding times in Ma.
#' @param include_plate_carriage Whether zero-cost plate-carriage edges are
#'   available to the path search.
#' @return A list containing `found`, `path_nodes`, `path_edges`, `path_cost`,
#'   and `path_distance_km`.
#' @export
#'
#' @examples
#' graph <- hee_dispersal_spacetime_graph(
#'   data.frame(from_cell_id = c("a", "b", "a"),
#'              to_cell_id = c("b", "c", "c"), time_ma = 0,
#'              effective_cost_km = c(1, 1, 10))
#' )
#' hee_dispersal_path_diagnostic(graph, "a", "c", 0, 0)
hee_dispersal_path_diagnostic <- function(graph,
                                          origin_cell_id,
                                          destination_cell_id,
                                          origin_time_ma,
                                          destination_time_ma,
                                          include_plate_carriage = TRUE) {
  if (!inherits(graph, "hee_dispersal_spacetime_graph")) {
    stop("graph must be a hee_dispersal_spacetime_graph object.", call. = FALSE)
  }
  edges <- graph$edges
  if (!isTRUE(include_plate_carriage)) {
    edges <- edges[edges$edge_kind != "plate_carriage", , drop = FALSE]
  }
  start <- .hee_dispersal_node_id(origin_cell_id, origin_time_ma)
  goal <- .hee_dispersal_node_id(destination_cell_id, destination_time_ma)
  nodes <- unique(c(edges$from_node, edges$to_node, start, goal))
  node_index <- stats::setNames(seq_along(nodes), nodes)
  start_i <- node_index[[start]]
  goal_i <- node_index[[goal]]
  from_i <- unname(node_index[edges$from_node])
  to_i <- unname(node_index[edges$to_node])
  adjacency <- split(seq_len(nrow(edges)), from_i)
  n <- length(nodes)
  distance <- rep(Inf, n)
  previous_node <- rep(NA_integer_, n)
  previous_edge <- rep(NA_integer_, n)
  distance[[start_i]] <- 0

  heap_node <- integer(); heap_key <- numeric()
  heap_push <- function(node, key) {
    heap_node <<- c(heap_node, node); heap_key <<- c(heap_key, key)
    pos <- length(heap_key)
    while (pos > 1L) {
      parent <- pos %/% 2L
      if (heap_key[[parent]] <= heap_key[[pos]]) break
      tmp_k <- heap_key[[parent]]; heap_key[[parent]] <<- heap_key[[pos]]; heap_key[[pos]] <<- tmp_k
      tmp_n <- heap_node[[parent]]; heap_node[[parent]] <<- heap_node[[pos]]; heap_node[[pos]] <<- tmp_n
      pos <- parent
    }
  }
  heap_pop <- function() {
    node <- heap_node[[1L]]; key <- heap_key[[1L]]
    last <- length(heap_key)
    if (last == 1L) {
      heap_node <<- integer(); heap_key <<- numeric()
      return(list(node = node, key = key))
    }
    heap_node[[1L]] <<- heap_node[[last]]; heap_key[[1L]] <<- heap_key[[last]]
    heap_node <<- heap_node[-last]; heap_key <<- heap_key[-last]
    pos <- 1L
    repeat {
      left <- pos * 2L; right <- left + 1L
      if (left > length(heap_key)) break
      child <- if (right <= length(heap_key) && heap_key[[right]] < heap_key[[left]]) right else left
      if (heap_key[[pos]] <= heap_key[[child]]) break
      tmp_k <- heap_key[[pos]]; heap_key[[pos]] <<- heap_key[[child]]; heap_key[[child]] <<- tmp_k
      tmp_n <- heap_node[[pos]]; heap_node[[pos]] <<- heap_node[[child]]; heap_node[[child]] <<- tmp_n
      pos <- child
    }
    list(node = node, key = key)
  }
  heap_push(start_i, 0)
  while (length(heap_key)) {
    current <- heap_pop()
    if (current$key != distance[[current$node]]) next
    if (current$node == goal_i) break
    outgoing <- adjacency[[as.character(current$node)]]
    if (is.null(outgoing)) next
    for (edge_i in outgoing) {
      target <- to_i[[edge_i]]
      alternative <- current$key + edges$path_cost[[edge_i]]
      if (alternative + 1e-12 < distance[[target]]) {
        distance[[target]] <- alternative
        previous_node[[target]] <- current$node
        previous_edge[[target]] <- edge_i
        heap_push(target, alternative)
      }
    }
  }
  if (!is.finite(distance[[goal_i]])) {
    return(list(
      found = FALSE,
      path_nodes = data.frame(), path_edges = edges[FALSE, , drop = FALSE],
      path_cost = Inf, path_distance_km = Inf,
      interpretation = "No directed path exists under the supplied Earth graph."
    ))
  }
  reverse_nodes <- c(goal_i); reverse_edges <- integer()
  cursor <- goal_i
  while (cursor != start_i) {
    reverse_edges <- c(reverse_edges, previous_edge[[cursor]])
    cursor <- previous_node[[cursor]]
    reverse_nodes <- c(reverse_nodes, cursor)
  }
  edge_index <- rev(reverse_edges)
  node_index_path <- rev(reverse_nodes)
  path_edges <- edges[edge_index, , drop = FALSE]
  path_nodes <- data.frame(
    node_id = nodes[node_index_path],
    stringsAsFactors = FALSE
  )
  list(
    found = TRUE,
    path_nodes = path_nodes,
    path_edges = path_edges,
    path_cost = distance[[goal_i]],
    path_distance_km = sum(path_edges$geographic_distance_km, na.rm = TRUE),
    n_active_dispersal_edges = sum(path_edges$active_dispersal),
    n_plate_carriage_edges = sum(!path_edges$active_dispersal),
    interpretation = paste(
      "Least-cost diagnostic path. Plate-carriage edges change geographic",
      "reference frame and are not counted as active biological dispersal."
    )
  )
}

.hee_dispersal_node_key <- function(tree, node) {
  if (node <= length(tree$tip.label)) tree$tip.label[[node]] else paste0("node_", node)
}

#' Evaluate a grid-state phylogeographic pruning likelihood
#'
#' Computes the Felsenstein pruning likelihood for a dated-tree topology and
#' explicit cell-state transition matrices. It is an optional constraint layer
#' for narrow-range phylogeographic data, not a replacement for the HmscEE
#' occupancy model. In particular, it does not by itself represent wide ranges,
#' undetected extinction, or a community process.
#'
#' For each branch, `branch_transition[[child]]` is a target-by-source matrix.
#' The message from child to parent is \eqn{K^T L_{child}}. Child messages are
#' multiplied at internal nodes and combined with the root prior. This permits
#' an explicit landscape/plate transition to constrain a tree without treating
#' plate carriage as active dispersal.
#'
#' @param tree An `ape::phylo`-compatible tree with a valid `edge` matrix.
#' @param tip_likelihood Named list of non-negative named numeric vectors. Each
#'   vector gives the likelihood of observed terminal data conditional on each
#'   grid state and is normalised internally.
#' @param branch_transition Named list of target-by-source transition matrices.
#'   Names must be a tip label or `"node_<number>"` for every child branch.
#' @param root_prior Optional named non-negative state prior. Defaults to a
#'   uniform prior over the state set.
#' @return A `hee_dispersal_pruning` object containing log likelihood, root
#'   posterior, node likelihoods and branch messages.
#' @export
#'
#' @examples
#' tree <- structure(list(
#'   edge = matrix(c(3, 1, 3, 2), ncol = 2, byrow = TRUE),
#'   tip.label = c("a", "b"), Nnode = 1L
#' ), class = "phylo")
#' K <- matrix(0.5, 2, 2); dimnames(K) <- list(c("x", "y"), c("x", "y"))
#' hee_dispersal_pruning_likelihood(
#'   tree, list(a = c(x = 1, y = 0), b = c(x = 0, y = 1)),
#'   list(a = K, b = K), root_prior = c(x = 0.5, y = 0.5)
#' )
hee_dispersal_pruning_likelihood <- function(tree, tip_likelihood,
                                              branch_transition,
                                              root_prior = NULL) {
  if (is.null(tree$edge) || is.null(tree$tip.label)) {
    stop("tree must provide edge and tip.label fields.", call. = FALSE)
  }
  edge <- as.matrix(tree$edge)
  if (ncol(edge) != 2L || nrow(edge) < 1L) {
    stop("tree$edge must be a non-empty two-column matrix.", call. = FALSE)
  }
  n_tip <- length(tree$tip.label)
  all_nodes <- sort(unique(c(edge[, 1L], edge[, 2L])))
  root <- setdiff(edge[, 1L], edge[, 2L])
  if (length(root) != 1L) stop("tree must have one root.", call. = FALSE)
  if (!is.list(tip_likelihood) || is.null(names(tip_likelihood)) ||
      !all(tree$tip.label %in% names(tip_likelihood))) {
    stop("tip_likelihood must be a named list covering every tree tip.",
         call. = FALSE)
  }
  first_tip <- tip_likelihood[[tree$tip.label[[1L]]]]
  states <- names(first_tip)
  if (is.null(states) || any(!nzchar(states)) || anyDuplicated(states)) {
    stop("Each tip likelihood must be a named vector with unique state names.",
         call. = FALSE)
  }
  n_state <- length(states)
  validate_vector <- function(v, label) {
    v <- suppressWarnings(as.numeric(v[states]))
    if (length(v) != n_state || any(!is.finite(v) | v < 0) || sum(v) <= 0) {
      stop(label, " must be a positive, finite non-negative likelihood over the shared state set.",
           call. = FALSE)
    }
    v / sum(v)
  }
  tip_likelihood <- lapply(tree$tip.label, function(label) {
    raw <- tip_likelihood[[label]]
    if (is.null(names(raw)) || !setequal(names(raw), states)) {
      stop("Every tip likelihood must use exactly the same named state set.",
           call. = FALSE)
    }
    stats::setNames(validate_vector(raw, paste0("tip_likelihood$", label)), states)
  })
  names(tip_likelihood) <- tree$tip.label
  if (is.null(root_prior)) root_prior <- stats::setNames(rep(1 / n_state, n_state), states)
  if (is.null(names(root_prior)) || !setequal(names(root_prior), states)) {
    stop("root_prior must use exactly the same named state set.", call. = FALSE)
  }
  root_prior <- validate_vector(root_prior, "root_prior")
  if (!is.list(branch_transition) || is.null(names(branch_transition))) {
    stop("branch_transition must be a named list of target-by-source matrices.",
         call. = FALSE)
  }
  transition_for_child <- function(child) {
    key <- .hee_dispersal_node_key(tree, child)
    K <- branch_transition[[key]]
    if (is.null(K)) K <- branch_transition[[as.character(child)]]
    if (is.null(K)) {
      stop("branch_transition is missing child branch: ", key, call. = FALSE)
    }
    K <- as.matrix(K)
    storage.mode(K) <- "double"
    if (!identical(dim(K), c(n_state, n_state)) || is.null(rownames(K)) ||
        is.null(colnames(K)) || !setequal(rownames(K), states) ||
        !setequal(colnames(K), states) || any(!is.finite(K) | K < 0)) {
      stop("Each transition must be a finite non-negative target-by-source matrix over the shared state set.",
           call. = FALSE)
    }
    K <- K[states, states, drop = FALSE]
    if (any(abs(colSums(K) - 1) > 1e-8)) {
      stop("Each target-by-source transition matrix must sum to one by source.",
           call. = FALSE)
    }
    K
  }
  children <- split(edge[, 2L], edge[, 1L])
  inside <- vector("list", max(all_nodes))
  node_log_scale <- rep(NA_real_, max(all_nodes))
  messages <- list()
  visit <- function(node) {
    if (!is.null(inside[[node]])) return(invisible(NULL))
    if (node <= n_tip) {
      label <- tree$tip.label[[node]]
      inside[[node]] <<- tip_likelihood[[label]]
      node_log_scale[[node]] <<- 0
      return(invisible(NULL))
    }
    child_nodes <- children[[as.character(node)]]
    likelihood <- rep(1, n_state); total_scale <- 0
    for (child in child_nodes) {
      visit(child)
      if (!is.finite(node_log_scale[[child]])) {
        inside[[node]] <<- rep(0, n_state)
        node_log_scale[[node]] <<- -Inf
        return(invisible(NULL))
      }
      K <- transition_for_child(child)
      message <- as.vector(crossprod(K, inside[[child]]))
      message_sum <- sum(message)
      if (!is.finite(message_sum) || message_sum <= 0) {
        inside[[node]] <<- rep(0, n_state)
        node_log_scale[[node]] <<- -Inf
        return(invisible(NULL))
      }
      message <- message / message_sum
      messages[[.hee_dispersal_node_key(tree, child)]] <<- stats::setNames(message, states)
      likelihood <- likelihood * message
      total_scale <- total_scale + node_log_scale[[child]] + log(message_sum)
    }
    scale <- sum(likelihood)
    if (!is.finite(scale) || scale <= 0) {
      inside[[node]] <<- rep(0, n_state)
      node_log_scale[[node]] <<- -Inf
    } else {
      inside[[node]] <<- likelihood / scale
      node_log_scale[[node]] <<- total_scale + log(scale)
    }
    invisible(NULL)
  }
  visit(root)
  root_inside <- inside[[root]]
  root_term <- sum(root_prior * root_inside)
  log_likelihood <- if (!is.finite(node_log_scale[[root]]) || root_term <= 0) {
    -Inf
  } else node_log_scale[[root]] + log(root_term)
  root_posterior <- if (is.finite(log_likelihood)) {
    stats::setNames(root_prior * root_inside / root_term, states)
  } else stats::setNames(rep(NA_real_, n_state), states)
  node_names <- vapply(all_nodes, function(node) .hee_dispersal_node_key(tree, node), character(1))
  names(inside) <- seq_along(inside)
  structure(
    list(
      log_likelihood = log_likelihood,
      root_node = root,
      root_posterior = root_posterior,
      node_inside = stats::setNames(inside[all_nodes], node_names),
      node_log_scale = stats::setNames(node_log_scale[all_nodes], node_names),
      branch_messages = messages,
      states = states,
      interpretation = paste(
        "Exact finite-state pruning constraint for explicit grid states.",
        "It constrains a phylogeographic history but is not a complete",
        "wide-range occupancy, extinction, or community model."
      )
    ),
    class = "hee_dispersal_pruning"
  )
}

#' Canonical finite-propagule dispersal sampler
#'
#' Preferred name for [hee_dispersal_particle_kernel()]. It samples a finite
#' number of propagules from an already defined edge kernel. These are
#' propagule-routing particles, not posterior particles, reconstructed
#' ancestral locations, or observed gene-flow events.
#'
#' @param ... Arguments passed to [hee_dispersal_particle_kernel()].
#' @return An edge table from [hee_dispersal_particle_kernel()].
#' @export
hee_dispersal_propagule_particles <- function(...) {
  hee_dispersal_particle_kernel(...)
}
