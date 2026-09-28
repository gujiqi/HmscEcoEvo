# Internal helper functions -------------------------------------------------

.as_matrix01 <- function(x, name = "comm") {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x)) stop(name, " must be a matrix or data.frame.", call. = FALSE)
  if (is.null(rownames(x)) || any(rownames(x) == "")) {
    stop(name, " must have site names as rownames.", call. = FALSE)
  }
  if (is.null(colnames(x)) || any(colnames(x) == "")) {
    stop(name, " must have species names as colnames.", call. = FALSE)
  }
  storage.mode(x) <- "numeric"
  ok <- is.na(x) | x %in% c(0, 1)
  if (!all(ok)) stop(name, " must contain only 0, 1, or NA.", call. = FALSE)
  x
}

.require_cols <- function(x, cols, table_name) {
  missing <- setdiff(cols, names(x))
  if (length(missing) > 0) {
    stop(table_name, " is missing required column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
}

.add_if_missing <- function(x, col, value) {
  if (!col %in% names(x)) x[[col]] <- value
  x
}

.align_rows <- function(x, row_order, name) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) {
    if (is.null(rownames(x))) stop(name, " must have rownames.", call. = FALSE)
    missing <- setdiff(row_order, rownames(x))
    if (length(missing) > 0) stop(name, " is missing rows: ", paste(missing, collapse = ", "), call. = FALSE)
    return(x[row_order, , drop = FALSE])
  }
  if (is.matrix(x)) {
    if (is.null(rownames(x))) stop(name, " must have rownames.", call. = FALSE)
    missing <- setdiff(row_order, rownames(x))
    if (length(missing) > 0) stop(name, " is missing rows: ", paste(missing, collapse = ", "), call. = FALSE)
    return(x[row_order, , drop = FALSE])
  }
  stop(name, " must be a matrix or data.frame.", call. = FALSE)
}

.align_traits <- function(x, species, name = "traits") {
  if (is.null(x)) return(NULL)
  if (!is.data.frame(x)) x <- as.data.frame(x)
  if ("species" %in% names(x)) {
    rownames(x) <- as.character(x$species)
    x$species <- NULL
  }
  if (is.null(rownames(x))) stop(name, " must have rownames or a species column.", call. = FALSE)
  missing <- setdiff(species, rownames(x))
  if (length(missing) > 0) {
    warning(name, " is missing species: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  x[intersect(species, rownames(x)), , drop = FALSE]
}

.std_prob <- function(x, default = 1) {
  if (is.null(x)) return(default)
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- default
  x
}

.validate_dispersal_scale <- function(x, name = "dispersal_scale") {
  val <- suppressWarnings(as.numeric(x))
  if (length(val) == 0L || any(!is.finite(val) | val < 0)) {
    stop(name, " must contain finite non-negative values.", call. = FALSE)
  }
  names(val) <- names(x)
  val
}

.hee_circular_mean_lon <- function(lon) {
  lon <- suppressWarnings(as.numeric(lon))
  lon <- lon[is.finite(lon)]
  if (length(lon) == 0L) return(NA_real_)
  rad <- lon * pi / 180
  atan2(mean(sin(rad)), mean(cos(rad))) * 180 / pi
}

.hee_great_circle_matrix <- function(lon1, lat1, lon2, lat2) {
  lon1 <- suppressWarnings(as.numeric(lon1))
  lat1 <- suppressWarnings(as.numeric(lat1))
  lon2 <- suppressWarnings(as.numeric(lon2))
  lat2 <- suppressWarnings(as.numeric(lat2))
  if (any(abs(lat1[is.finite(lat1)]) > 90) ||
      any(abs(lat2[is.finite(lat2)]) > 90)) {
    stop("Latitude values must be between -90 and 90 degrees.", call. = FALSE)
  }
  n1 <- length(lon1)
  n2 <- length(lon2)
  if (n1 == 0L || n2 == 0L) {
    return(matrix(numeric(), nrow = n1, ncol = n2))
  }
  bad1 <- !is.finite(lon1) | !is.finite(lat1)
  bad2 <- !is.finite(lon2) | !is.finite(lat2)
  rad <- pi / 180
  lat1r <- lat1 * rad
  lat2r <- lat2 * rad
  dlat <- outer(lat1r, lat2r, function(a, b) b - a)
  dlon_deg <- outer(lon1, lon2, function(a, b) ((b - a + 180) %% 360) - 180)
  dlon <- dlon_deg * rad
  a <- sin(dlat / 2)^2 +
    outer(cos(lat1r), cos(lat2r), `*`) * sin(dlon / 2)^2
  a <- pmin(pmax(a, 0), 1)
  out <- 2 * 6371.0088 * asin(sqrt(a))
  out[bad1, ] <- NA_real_
  out[, bad2] <- NA_real_
  out
}

.hee_region_centroids <- function(x, region_col = "region",
                                  lon_col = "lon", lat_col = "lat") {
  x <- as.data.frame(x)
  .require_cols(x, c(region_col, lon_col, lat_col), "region coordinate table")
  rows <- lapply(split(x, as.character(x[[region_col]]), drop = TRUE), function(z) {
    data.frame(
      region = as.character(z[[region_col]][1]),
      lon = .hee_circular_mean_lon(z[[lon_col]]),
      lat = mean(suppressWarnings(as.numeric(z[[lat_col]])), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  names(out)[names(out) == "region"] <- region_col
  names(out)[names(out) == "lon"] <- lon_col
  names(out)[names(out) == "lat"] <- lat_col
  out
}

.weighted_mean_safe <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

.weighted_sd_safe <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (sum(ok) <= 1) return(NA_real_)
  m <- sum(x[ok] * w[ok]) / sum(w[ok])
  sqrt(sum(w[ok] * (x[ok] - m)^2) / sum(w[ok]))
}

.weighted_mode <- function(value, weight) {
  ok <- !is.na(value) & !is.na(weight)
  if (!any(ok)) return(NA_character_)
  s <- tapply(weight[ok], value[ok], sum, na.rm = TRUE)
  names(s)[which.max(s)]
}

.safe_bind_rows <- function(x) {
  x <- x[!vapply(x, is.null, logical(1))]
  if (length(x) == 0) return(NULL)
  cols <- unique(unlist(lapply(x, names)))
  out <- lapply(x, function(d) {
    miss <- setdiff(cols, names(d))
    for (m in miss) d[[m]] <- NA
    d[cols]
  })
  do.call(rbind, out)
}

.map_event_type <- function(x) {
  y <- tolower(trimws(as.character(x)))
  out <- y
  out[grepl("^d$|disp|range.?expan|colon", y)] <- "dispersal"
  out[grepl("^e$|loss|extirp|range.?contract|local.?ext", y)] <- "loss"
  out[grepl("spec|clado|split|vicariance|sympatry", y)] <- "speciation"
  out[grepl("founder|jump", y)] <- "founder"
  out
}

.rename_first_existing <- function(d, target, candidates) {
  if (target %in% names(d)) return(d)
  hit <- intersect(candidates, names(d))
  if (length(hit) > 0) names(d)[match(hit[1], names(d))] <- target
  d
}

.find_data_frames <- function(x, max_depth = 4) {
  out <- list()
  walk <- function(obj, path, depth) {
    if (depth > max_depth) return()
    if (is.data.frame(obj)) {
      out[[paste(path, collapse = "/")]] <<- obj
    } else if (is.list(obj)) {
      nms <- names(obj)
      if (is.null(nms)) nms <- as.character(seq_along(obj))
      for (i in seq_along(obj)) walk(obj[[i]], c(path, nms[i]), depth + 1)
    }
  }
  walk(x, character(0), 0)
  out
}

.event_like_score <- function(d) {
  nms <- tolower(names(d))
  sum(grepl("event|type|time|age|from|to|region|area|species|taxon|node|branch|prob", nms))
}

.cbind_keep <- function(a, b) {
  if (is.null(b) || ncol(as.data.frame(b)) == 0) return(a)
  b <- as.data.frame(b)
  if (ncol(a) == 0) return(b[rownames(a), , drop = FALSE])
  if (!identical(rownames(a), rownames(b))) {
    b <- b[rownames(a), , drop = FALSE]
  }
  cbind(a, b)
}

.drop_all_na_cols <- function(x) {
  x <- as.data.frame(x)
  if (ncol(x) == 0) return(x)
  keep <- vapply(x, function(z) any(!is.na(z)), logical(1))
  x[, keep, drop = FALSE]
}

.z <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(ifelse(is.na(x), NA_real_, 0))
  (x - m) / s
}

.numeric_cols <- function(x) {
  if (is.null(x) || ncol(as.data.frame(x)) == 0) return(character())
  nms <- names(x)
  nms[vapply(x, is.numeric, logical(1))]
}
