#' Static dispersal accessibility from historical source cells
#'
#' Compute a static dispersal filter
#' \deqn{D_static[j,c,t] = exp(-d_km(c, source[j,t]) / scale_km[j]) * barrier[c,t]}
#' where source cells are cells in historically accessible regions. This does
#' not use the previous time slice, so it is a distance/connectivity filter
#' rather than a dynamic colonisation model.
#'
#' @param x Target cell table. It may already be species-expanded; otherwise it
#'   is expanded using `accessibility`.
#' @param accessibility Optional table with `species`, `region`, `time_ma`, and
#'   `accessibility`. Regions with accessibility greater than or equal to
#'   `source_threshold` are treated as sources. If no region passes the
#'   threshold, the best region for that species/time is used.
#' @param source_cells Optional explicit source table with `species`,
#'   `time_ma`, and either `cell_id` or coordinates.
#' @param dispersal_scale Named numeric vector or scalar in kilometres. Larger
#'   values make distance decay slower. Values must be finite and non-negative;
#'   zero means no movement across positive distance.
#' @param barrier Optional table with `barrier_passability` by cell/time or
#'   species/cell/time. If omitted, barriers are neutral. If supplied, unmatched
#'   rows or missing passability values are treated conservatively as closed
#'   (`0`), so users should provide explicit `1` values for cells known to be
#'   passable.
#' @param source_threshold Accessibility threshold for source regions.
#' @param species_col,cell_col,time_col,region_col,lon_col,lat_col Column names.
#'
#' @return A table with `static_weight` and `D_static`.
#' @export
#'
#' @examples
#' cells <- data.frame(cell_id = 1:3, time_ma = 0, lon = c(0, 1, 3),
#'   lat = 0, region = c("A", "A", "B"))
#' acc <- data.frame(species = "sp1", region = "A", time_ma = 0,
#'   accessibility = 1)
#' hee_static_dispersal(cells, acc, dispersal_scale = c(sp1 = 2))
hee_static_dispersal <- function(x,
                                 accessibility = NULL,
                                 source_cells = NULL,
                                 dispersal_scale = 1,
                                 barrier = NULL,
                                 source_threshold = 0.5,
                                 species_col = "species",
                                 cell_col = "cell_id",
                                 time_col = "time_ma",
                                 region_col = "region",
                                 lon_col = "lon",
                                 lat_col = "lat") {
  dispersal_scale <- .validate_dispersal_scale(dispersal_scale)
  cells <- as.data.frame(x)
  .require_cols(cells, c(cell_col, time_col, lon_col, lat_col), "x")
  if (!species_col %in% names(cells)) {
    if (is.null(accessibility)) {
      stop("`x` must contain `species`, or `accessibility` must be supplied.", call. = FALSE)
    }
    accessibility <- as.data.frame(accessibility)
    .require_cols(accessibility, c(species_col, region_col, time_col), "accessibility")
    sp_time <- unique(accessibility[, c(species_col, time_col), drop = FALSE])
    cells <- merge(sp_time, cells, by = time_col, all.x = TRUE, sort = FALSE)
  }
  if (!region_col %in% names(cells) && is.null(source_cells)) {
    stop("`x` must contain `region` when `source_cells` is not supplied.", call. = FALSE)
  }
  cells$.row_id <- seq_len(nrow(cells))

  if (is.null(source_cells)) {
    acc <- as.data.frame(accessibility)
    .require_cols(acc, c(species_col, region_col, time_col), "accessibility")
    if (!"accessibility" %in% names(acc)) acc$accessibility <- 1
    acc$accessibility <- pmin(pmax(suppressWarnings(as.numeric(acc$accessibility)), 0), 1)
    acc$accessibility[is.na(acc$accessibility)] <- 0
    .hee_check_unique_keys(acc, c(species_col, region_col, time_col),
                           "accessibility")
    sources <- lapply(split(acc, interaction(acc[[species_col]], acc[[time_col]], drop = TRUE)), function(z) {
      z <- z[order(z$accessibility, decreasing = TRUE), , drop = FALSE]
      if (!any(z$accessibility > 0, na.rm = TRUE)) {
        return(z[FALSE, c(species_col, region_col, time_col), drop = FALSE])
      }
      keep <- z[z$accessibility >= source_threshold, , drop = FALSE]
      if (nrow(keep) == 0L) keep <- z[1, , drop = FALSE]
      merge(unique(keep[, c(species_col, region_col, time_col), drop = FALSE]),
            unique(cells[, c(cell_col, time_col, region_col, lon_col, lat_col), drop = FALSE]),
            by = c(time_col, region_col), all.x = FALSE, sort = FALSE)
    })
    sources <- do.call(rbind, sources)
  } else {
    sources <- as.data.frame(source_cells)
    .require_cols(sources, c(species_col, time_col), "source_cells")
    if (!all(c(lon_col, lat_col) %in% names(sources))) {
      .require_cols(sources, cell_col, "source_cells")
      sources <- merge(sources,
                       unique(cells[, c(cell_col, time_col, lon_col, lat_col), drop = FALSE]),
                       by = c(cell_col, time_col), all.x = TRUE, sort = FALSE)
    }
  }
  if (is.null(sources) || nrow(sources) == 0L) {
    out <- cells
    out$distance_to_accessible_source <- NA_real_
    out$static_weight <- 0
    out$D_static <- out$static_weight
    return(out[, setdiff(names(out), ".row_id"), drop = FALSE])
  }

  scale_for <- function(sp) {
    if (length(dispersal_scale) == 1L) return(as.numeric(dispersal_scale))
    val <- dispersal_scale[as.character(sp)]
    ifelse(is.na(val), stats::median(dispersal_scale, na.rm = TRUE), as.numeric(val))
  }

  out_list <- lapply(split(cells, interaction(cells[[species_col]], cells[[time_col]], drop = TRUE)), function(target) {
    sp <- target[[species_col]][1]
    tt <- target[[time_col]][1]
    src <- sources[sources[[species_col]] == sp & sources[[time_col]] == tt, , drop = FALSE]
    if (nrow(src) == 0L) {
      src <- sources[sources[[time_col]] == tt, , drop = FALSE]
    }
    if (nrow(src) == 0L) {
      target$distance_to_accessible_source <- Inf
      target$static_weight <- 0
      return(target)
    }
    d <- .hee_great_circle_matrix(target[[lon_col]], target[[lat_col]],
                                  src[[lon_col]], src[[lat_col]])
    md <- apply(d, 1, min, na.rm = TRUE)
    sc <- max(scale_for(sp), .Machine$double.eps)
    target$distance_to_accessible_source <- md
    target$static_weight <- exp(-md / sc)
    target
  })
  out <- do.call(rbind, out_list)
  if (!is.null(barrier)) {
    b <- as.data.frame(barrier)
    if (!"barrier_passability" %in% names(b)) b$barrier_passability <- 1
    keys <- intersect(c(species_col, cell_col, time_col), names(b))
    if (!species_col %in% keys) keys <- intersect(c(cell_col, time_col), names(b))
    if (length(keys) == 0L) {
      if (nrow(b) != 1L) {
        stop("barrier must contain join keys or exactly one global row.",
             call. = FALSE)
      }
      out$barrier_passability <- .hee_clip01(b$barrier_passability[1])
      out$barrier_passability[is.na(out$barrier_passability)] <- 0
      out$static_weight <- out$static_weight * out$barrier_passability
    } else {
      .hee_check_unique_keys(b, keys, "barrier")
      before_n <- nrow(out)
      out <- merge(out, b[, unique(c(keys, "barrier_passability")), drop = FALSE],
                   by = keys, all.x = TRUE, sort = FALSE)
      if (nrow(out) != before_n) {
        stop("Join from barrier changed row count from ", before_n, " to ",
             nrow(out), ". Check keys: ", paste(keys, collapse = ", "),
             call. = FALSE)
      }
      out$barrier_passability[is.na(out$barrier_passability)] <- 0
      out$static_weight <- out$static_weight * out$barrier_passability
    }
  }
  out$static_weight <- pmin(pmax(out$static_weight, 0), 1)
  out$D_static <- out$static_weight
  out <- out[order(out$.row_id), , drop = FALSE]
  out[, setdiff(names(out), ".row_id"), drop = FALSE]
}

#' Region-level static dispersal accessibility
#'
#' Fast approximation to [hee_static_dispersal()] for deep-time grids. The
#' distance decay is computed among regions and then joined back to every cell:
#' \deqn{D_static[j,r,t] = exp(-d_km(r, accessible source[j,t]) / scale_km[j]).}
#' This is appropriate when historical accessibility is regional, as in a
#' BioGeoBEARS-style workflow.
#'
#' @param x Species-cell-time table with a `region` column.
#' @param accessibility Table with `species`, `region`, `time_ma`, and
#'   `accessibility`.
#' @param region_distance Optional square region distance matrix. If omitted,
#'   region centroids are estimated from `x`.
#' @param dispersal_scale Named numeric vector or scalar in kilometres. Values
#'   must be finite and non-negative; zero means no movement across positive
#'   distance.
#' @param source_threshold Regions with accessibility at or above this value
#'   are treated as source regions.
#' @param species_col,time_col,region_col,lon_col,lat_col Column names.
#'
#' @return `x` with `static_weight`, `D_static`, and
#'   `distance_to_accessible_source`.
#' @export
hee_static_region_dispersal <- function(x,
                                        accessibility,
                                        region_distance = NULL,
                                        dispersal_scale = 1,
                                        source_threshold = 0.5,
                                        species_col = "species",
                                        time_col = "time_ma",
                                        region_col = "region",
                                        lon_col = "lon",
                                        lat_col = "lat") {
  dispersal_scale <- .validate_dispersal_scale(dispersal_scale)
  cells <- as.data.frame(x)
  .require_cols(cells, c(species_col, time_col, region_col), "x")
  acc <- as.data.frame(accessibility)
  .require_cols(acc, c(species_col, region_col, time_col), "accessibility")
  if (!"accessibility" %in% names(acc)) acc$accessibility <- 1
  acc$accessibility <- pmin(pmax(suppressWarnings(as.numeric(acc$accessibility)), 0), 1)
  acc$accessibility[is.na(acc$accessibility)] <- 0
  .hee_check_unique_keys(acc, c(species_col, region_col, time_col),
                         "accessibility")

  regions <- sort(unique(c(as.character(cells[[region_col]]),
                           as.character(acc[[region_col]]))))
  if (is.null(region_distance)) {
    .require_cols(cells, c(lon_col, lat_col), "x")
    cent <- .hee_region_centroids(cells, region_col, lon_col, lat_col)
    cent <- cent[match(regions, cent[[region_col]]), , drop = FALSE]
    region_distance <- .hee_great_circle_matrix(
      cent[[lon_col]], cent[[lat_col]], cent[[lon_col]], cent[[lat_col]]
    )
    dimnames(region_distance) <- list(regions, regions)
  } else {
    region_distance <- as.matrix(region_distance)
  }
  if (any(!is.finite(region_distance) | region_distance < 0, na.rm = TRUE)) {
    stop("region_distance must contain finite non-negative distances in kilometres.",
         call. = FALSE)
  }

  scale_for <- function(sp) {
    if (length(dispersal_scale) == 1L) return(as.numeric(dispersal_scale))
    val <- dispersal_scale[as.character(sp)]
    ifelse(is.na(val), stats::median(dispersal_scale, na.rm = TRUE), as.numeric(val))
  }

  targets <- unique(cells[, c(species_col, time_col, region_col), drop = FALSE])
  targets$.row_id <- seq_len(nrow(targets))
  out <- lapply(split(targets, interaction(targets[[species_col]], targets[[time_col]], drop = TRUE)), function(z) {
    sp <- z[[species_col]][1]
    tt <- z[[time_col]][1]
    src <- acc[acc[[species_col]] == sp & acc[[time_col]] == tt, , drop = FALSE]
    if (nrow(src) == 0L) src <- acc[acc[[time_col]] == tt, , drop = FALSE]
    if (nrow(src) == 0L || !any(src$accessibility > 0, na.rm = TRUE)) {
      z$distance_to_accessible_source <- Inf
      z$static_weight <- 0
      return(z)
    }
    src <- src[order(src$accessibility, decreasing = TRUE), , drop = FALSE]
    keep <- src[src$accessibility >= source_threshold, , drop = FALSE]
    if (nrow(keep) == 0L) keep <- src[1, , drop = FALSE]
    source_regions <- intersect(as.character(keep[[region_col]]), colnames(region_distance))
    target_regions <- as.character(z[[region_col]])
    target_regions <- ifelse(target_regions %in% rownames(region_distance),
                             target_regions, NA_character_)
    if (length(source_regions) == 0L || anyNA(target_regions)) {
      z$distance_to_accessible_source <- Inf
      z$static_weight <- 0
      return(z)
    }
    d <- region_distance[target_regions, source_regions, drop = FALSE]
    md <- apply(d, 1, min, na.rm = TRUE)
    sc <- max(scale_for(sp), .Machine$double.eps)
    z$distance_to_accessible_source <- md
    z$static_weight <- exp(-md / sc)
    z
  })
  weights <- do.call(rbind, out)
  weights$static_weight <- pmin(pmax(weights$static_weight, 0), 1)
  weights$D_static <- weights$static_weight
  weights <- weights[order(weights$.row_id), setdiff(names(weights), ".row_id"), drop = FALSE]
  out <- merge(cells, weights, by = c(species_col, time_col, region_col),
               all.x = TRUE, sort = FALSE)
  if (nrow(out) != nrow(cells)) {
    stop("Join from static region dispersal weights changed row count from ",
         nrow(cells), " to ", nrow(out), ".", call. = FALSE)
  }
  out
}

#' Region-level dynamic dispersal filter
#'
#' Fast dynamic colonisation filter using the previous time slice summarized by
#' species and region:
#' \deqn{D_dynamic[j,r,t] = \max_{r'} P[j,r',t_{prev}] K(d_km(r',r), scale_km[j]).}
#' Use `time_direction = "forward"` for older-to-younger projection and
#' `"backward"` for envelope-style reconstructions.
#'
#' @param previous_projection Projection table with `species`, `region`,
#'   `time_ma`, and a probability column.
#' @param cell_table Optional target table. If omitted, `previous_projection` is
#'   used as the target.
#' @param region_distance Optional square region distance matrix.
#' @param dispersal_scale Named numeric vector or scalar in kilometres. Values
#'   must be finite and non-negative; zero means no movement across positive
#'   distance.
#' @param time_direction `"forward"` or `"backward"`.
#' @param kernel Distance kernel.
#' @param default Weight used only for the first time slice with no previous
#'   slice. Later slices with no source evidence get zero dynamic reachability.
#' @param probability_col Probability column in `previous_projection`.
#' @param species_col,time_col,region_col,lon_col,lat_col Column names.
#'
#' @return Target table with `dynamic_weight`, `D_dynamic`, and
#'   `dynamic_source_time_ma`.
#' @export
hee_dynamic_region_filter <- function(previous_projection,
                                      cell_table = NULL,
                                      region_distance = NULL,
                                      dispersal_scale = 1,
                                      time_direction = c("forward", "backward"),
                                      kernel = c("exponential", "gaussian"),
                                      default = 1,
                                      probability_col = "probability",
                                      species_col = "species",
                                      time_col = "time_ma",
                                      region_col = "region",
                                      lon_col = "lon",
                                      lat_col = "lat") {
  time_direction <- match.arg(time_direction)
  kernel <- match.arg(kernel)
  dispersal_scale <- .validate_dispersal_scale(dispersal_scale)
  prev <- as.data.frame(previous_projection)
  .require_cols(prev, c(species_col, time_col, region_col, probability_col), "previous_projection")
  target <- if (is.null(cell_table)) prev else as.data.frame(cell_table)
  .require_cols(target, c(species_col, time_col, region_col), "cell_table")

  regions <- sort(unique(c(as.character(prev[[region_col]]),
                           as.character(target[[region_col]]))))
  if (is.null(region_distance)) {
    coord_source <- target
    if (!all(c(lon_col, lat_col) %in% names(coord_source))) coord_source <- prev
    .require_cols(coord_source, c(lon_col, lat_col), "cell_table")
    cent <- .hee_region_centroids(coord_source, region_col, lon_col, lat_col)
    cent <- cent[match(regions, cent[[region_col]]), , drop = FALSE]
    region_distance <- .hee_great_circle_matrix(
      cent[[lon_col]], cent[[lat_col]], cent[[lon_col]], cent[[lat_col]]
    )
    dimnames(region_distance) <- list(regions, regions)
  } else {
    region_distance <- as.matrix(region_distance)
  }
  if (any(!is.finite(region_distance) | region_distance < 0, na.rm = TRUE)) {
    stop("region_distance must contain finite non-negative distances in kilometres.",
         call. = FALSE)
  }

  scale_for <- function(sp) {
    if (length(dispersal_scale) == 1L) return(as.numeric(dispersal_scale))
    val <- dispersal_scale[as.character(sp)]
    ifelse(is.na(val), stats::median(dispersal_scale, na.rm = TRUE), as.numeric(val))
  }
  kfun <- function(d, sc) {
    if (kernel == "gaussian") exp(-(d^2) / (2 * sc^2)) else exp(-d / sc)
  }

  src <- stats::aggregate(prev[[probability_col]],
                          prev[, c(species_col, time_col, region_col), drop = FALSE],
                          max, na.rm = TRUE)
  names(src)[ncol(src)] <- ".source_probability"
  times <- sort(unique(c(prev[[time_col]], target[[time_col]])), decreasing = TRUE)
  target$.row_id <- seq_len(nrow(target))
  target_key <- unique(target[, c(species_col, time_col, region_col), drop = FALSE])
  target_key$.key_id <- seq_len(nrow(target_key))
  out <- lapply(split(target_key, interaction(target_key[[species_col]], target_key[[time_col]], drop = TRUE)), function(z) {
    sp <- z[[species_col]][1]
    tt <- z[[time_col]][1]
    idx <- match(tt, times)
    src_time <- if (time_direction == "forward") times[idx - 1L] else times[idx + 1L]
    if (length(src_time) == 0L || is.na(src_time)) {
      z$dynamic_weight <- default
      z$dynamic_source_time_ma <- NA_real_
      return(z)
    }
    s <- src[src[[species_col]] == sp & src[[time_col]] == src_time, , drop = FALSE]
    if (nrow(s) == 0L) {
      z$dynamic_weight <- 0
      z$dynamic_source_time_ma <- src_time
      return(z)
    }
    tr <- as.character(z[[region_col]])
    sr <- intersect(as.character(s[[region_col]]), colnames(region_distance))
    tr <- ifelse(tr %in% rownames(region_distance), tr, NA_character_)
    if (length(sr) == 0L || anyNA(tr)) {
      z$dynamic_weight <- 0
      z$dynamic_source_time_ma <- src_time
      return(z)
    }
    s <- s[match(sr, as.character(s[[region_col]])), , drop = FALSE]
    d <- region_distance[tr, sr, drop = FALSE]
    sc <- max(scale_for(sp), .Machine$double.eps)
    score <- sweep(kfun(d, sc), 2, s$.source_probability, "*")
    z$dynamic_weight <- apply(score, 1, max, na.rm = TRUE)
    z$dynamic_source_time_ma <- src_time
    z
  })
  weights <- do.call(rbind, out)
  weights$dynamic_weight <- pmin(pmax(weights$dynamic_weight, 0), 1)
  weights$D_dynamic <- weights$dynamic_weight
  weights <- weights[order(weights$.key_id), setdiff(names(weights), ".key_id"), drop = FALSE]
  out <- merge(target, weights, by = c(species_col, time_col, region_col), all.x = TRUE, sort = FALSE)
  out <- out[order(out$.row_id), , drop = FALSE]
  out[, setdiff(names(out), ".row_id"), drop = FALSE]
}
