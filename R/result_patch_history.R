# Derived spatial results. These are not additional biological processes.

#' Identify connected environmentally supported land patches
#'
#' `support` is an environmental-response index, not occupancy. Patches are
#' connected components of cells exceeding a declared threshold on a supplied
#' undirected neighbourhood graph. The result never links cells across sea
#' merely because they are close in longitude and latitude.
#'
#' @param cells One-time, one-lineage table with `cell_id`, `support`,
#'   `land_area_km2`, and optionally `location_density`.
#' @param edges Table with `from_cell_id`, `to_cell_id` for terrestrial
#'   adjacency at this time.
#' @param threshold Absolute support-index threshold in `[0,1]`.
#' @param min_area_km2 Minimum patch area retained.
#' @param lineage_id,time_ma Identifiers written to both outputs.
#' @return List with `membership` and `patches` tables.
#' @export
hee_result_patches <- function(cells, edges, threshold = 0.6,
                               min_area_km2 = 0, lineage_id, time_ma) {
  x <- as.data.frame(cells, stringsAsFactors = FALSE)
  e <- as.data.frame(edges, stringsAsFactors = FALSE)
  needed <- c("cell_id", "support", "land_area_km2")
  if (!all(needed %in% names(x)) ||
      !all(c("from_cell_id", "to_cell_id") %in% names(e))) {
    stop("cells or edges lack required columns.", call. = FALSE)
  }
  if (anyDuplicated(x$cell_id) || anyNA(x$cell_id) ||
      any(!is.finite(x$land_area_km2) | x$land_area_km2 < 0) ||
      any(!is.na(x$support) & (!is.finite(x$support) | x$support < 0 | x$support > 1)) ||
      length(threshold) != 1L || !is.finite(threshold) || threshold < 0 || threshold > 1 ||
      length(min_area_km2) != 1L || !is.finite(min_area_km2) || min_area_km2 < 0) {
    stop("Invalid cells, threshold, or minimum patch area.", call. = FALSE)
  }
  if (!"location_density" %in% names(x)) x$location_density <- NA_real_
  if (any(!is.na(x$location_density) &
          (!is.finite(x$location_density) | x$location_density < 0))) {
    stop("location_density must be non-negative when supplied.", call. = FALSE)
  }
  active <- which(x$land_area_km2 > 0 & !is.na(x$support) & x$support >= threshold)
  empty_membership <- data.frame(lineage_id = character(), time_ma = numeric(),
    patch_id = character(), cell_id = character(), land_area_km2 = numeric(),
    support = numeric(), location_density = numeric(), stringsAsFactors = FALSE)
  empty_patches <- data.frame(lineage_id = character(), time_ma = numeric(),
    patch_id = character(), n_cells = integer(), area_km2 = numeric(),
    mean_support = numeric(), location_mass = numeric(),
    stringsAsFactors = FALSE)
  if (!length(active)) return(list(membership = empty_membership,
                                   patches = empty_patches))
  id <- match(as.character(e$from_cell_id), x$cell_id)
  jd <- match(as.character(e$to_cell_id), x$cell_id)
  keep <- !is.na(id) & !is.na(jd) & id != jd
  is_active <- rep(FALSE, nrow(x)); is_active[active] <- TRUE
  keep <- keep & is_active[id] & is_active[jd]
  keep[is.na(keep)] <- FALSE
  parent <- seq_len(nrow(x))
  find_root <- function(i) {
    while (parent[[i]] != i) {
      parent[[i]] <<- parent[[parent[[i]]]]
      i <- parent[[i]]
    }
    i
  }
  for (k in which(keep)) {
    a <- find_root(id[[k]]); b <- find_root(jd[[k]])
    if (a != b) parent[[max(a, b)]] <- min(a, b)
  }
  roots <- vapply(active, find_root, integer(1))
  ordered_roots <- unique(roots)
  patch_number <- match(roots, ordered_roots)
  patch_id <- sprintf("%s_%s_p%04d", as.character(lineage_id),
                      format(as.numeric(time_ma), trim = TRUE), patch_number)
  m <- data.frame(lineage_id = as.character(lineage_id), time_ma = as.numeric(time_ma),
    patch_id = patch_id, cell_id = as.character(x$cell_id[active]),
    land_area_km2 = x$land_area_km2[active], support = x$support[active],
    location_density = x$location_density[active], stringsAsFactors = FALSE)
  area <- rowsum(m$land_area_km2, m$patch_id, reorder = FALSE)
  retained <- rownames(area)[area[, 1L] >= min_area_km2]
  m <- m[m$patch_id %in% retained, , drop = FALSE]
  if (!nrow(m)) return(list(membership = empty_membership,
                            patches = empty_patches))
  parts <- split(seq_len(nrow(m)), m$patch_id)
  p <- do.call(rbind, lapply(parts, function(ii) {
    z <- m[ii, , drop = FALSE]
    original <- x[match(z$cell_id, x$cell_id), , drop = FALSE]
    environment_mean <- function(field) {
      if (!field %in% names(original)) return(NA_real_)
      value <- suppressWarnings(as.numeric(original[[field]]))
      ok <- is.finite(value)
      if (!any(ok)) return(NA_real_)
      stats::weighted.mean(value[ok], z$land_area_km2[ok])
    }
    noanalog <- if ("no_analog" %in% names(original)) {
      value <- as.logical(original$no_analog)
      if (anyNA(value)) NA_real_ else stats::weighted.mean(as.numeric(value),
                                                     z$land_area_km2)
    } else NA_real_
    data.frame(lineage_id = z$lineage_id[[1L]], time_ma = z$time_ma[[1L]],
      patch_id = z$patch_id[[1L]], n_cells = nrow(z),
      area_km2 = sum(z$land_area_km2),
      mean_support = stats::weighted.mean(z$support, z$land_area_km2),
      mean_temperature = environment_mean("MAT_pohl_C"),
      mean_precipitation = environment_mean("MAP_pohl_mm_yr"),
      noanalog_area_fraction = noanalog,
      location_mass = if (all(is.na(z$location_density))) NA_real_
                      else sum(z$location_density, na.rm = TRUE),
      stringsAsFactors = FALSE)
  }))
  rownames(p) <- NULL
  list(membership = m, patches = p)
}

#' Link patches through observed land-unit carriage
#'
#' A link is created only where the supplied source-to-target land transfer
#' carries cells from a source patch into a target patch. One-to-many and
#' many-to-one links are retained; no arbitrary nearest-centroid matching is
#' performed. `carriage_weight` is the fraction of a source cell's land area
#' delivered to a target cell, not dispersal probability.
#'
#' @param older,younger Patch memberships returned by [hee_result_patches()].
#' @param carriage Table with `source_cell_id`, `target_cell_id`, and
#'   `carriage_weight`.
#' @param min_source_share Minimum transferred-area fraction of source patch.
#' @return Directed patch-link table.
#' @export
hee_result_patch_links <- function(older, younger, carriage,
                                   min_source_share = 0.05) {
  a <- as.data.frame(older); b <- as.data.frame(younger)
  tr <- as.data.frame(carriage)
  required <- c("source_cell_id", "target_cell_id", "carriage_weight")
  if (!all(c("patch_id", "cell_id", "land_area_km2", "time_ma") %in% names(a)) ||
      !all(c("patch_id", "cell_id", "time_ma") %in% names(b)) ||
      !all(required %in% names(tr))) {
    stop("Patch memberships or carriage lack required columns.", call. = FALSE)
  }
  empty <- data.frame(from_patch_id = character(), to_patch_id = character(),
    time_from_ma = numeric(), time_to_ma = numeric(),
    carried_area_km2 = numeric(), source_share = numeric(),
    stringsAsFactors = FALSE)
  if (!nrow(a) || !nrow(b) || !nrow(tr)) return(empty)
  if (length(unique(a$time_ma)) != 1L || length(unique(b$time_ma)) != 1L ||
      unique(a$time_ma) <= unique(b$time_ma) ||
      any(!is.finite(tr$carriage_weight) | tr$carriage_weight < 0 |
          tr$carriage_weight > 1) ||
      !is.finite(min_source_share) || min_source_share < 0 || min_source_share > 1) {
    stop("Carriage weights, time direction, or threshold are invalid.", call. = FALSE)
  }
  source_weight <- rowsum(tr$carriage_weight, tr$source_cell_id,
                          reorder = FALSE)
  if (any(source_weight[, 1L] > 1 + 1e-8)) {
    stop("Carriage weights must sum to at most one per source cell.",
         call. = FALSE)
  }
  source <- a[, c("cell_id", "patch_id", "land_area_km2")]
  names(source) <- c("source_cell_id", "from_patch_id", "source_area")
  target <- b[, c("cell_id", "patch_id")]
  names(target) <- c("target_cell_id", "to_patch_id")
  joined <- merge(merge(tr[, required], source, by = "source_cell_id"),
                  target, by = "target_cell_id")
  if (!nrow(joined)) return(empty)
  joined$carried <- joined$source_area * joined$carriage_weight
  sums <- stats::aggregate(joined$carried,
    list(from_patch_id = joined$from_patch_id,
         to_patch_id = joined$to_patch_id), sum)
  names(sums)[[3L]] <- "carried_area_km2"
  source_area <- stats::aggregate(a$land_area_km2,
    list(from_patch_id = a$patch_id), sum)
  names(source_area)[[2L]] <- "source_patch_area"
  out <- merge(sums, source_area, by = "from_patch_id")
  out$source_share <- out$carried_area_km2 / out$source_patch_area
  out <- out[out$source_share >= min_source_share & out$carried_area_km2 > 0,
             c("from_patch_id", "to_patch_id", "carried_area_km2", "source_share"),
             drop = FALSE]
  out$time_from_ma <- unique(a$time_ma)
  out$time_to_ma <- unique(b$time_ma)
  out[, names(empty), drop = FALSE]
}

#' Summarise connected patch-history networks
#'
#' Splits and mergers remain visible in `links`. The component identifier is
#' a history-network label, not a claim that one unchanged patch persisted.
#' @param patches Patch summary table from multiple times.
#' @param links Directed links from [hee_result_patch_links()].
#' @return Tables `patches` (with network IDs), `networks`, and link-level
#'   environmental changes. Network-wide standard deviations also contain
#'   variation among simultaneous patches, and are not purely temporal rates.
#' @export
hee_result_patch_stability <- function(patches, links) {
  p <- as.data.frame(patches); e <- as.data.frame(links)
  if (!all(c("patch_id", "time_ma", "area_km2", "mean_support") %in% names(p)) ||
      !all(c("from_patch_id", "to_patch_id") %in% names(e)) ||
      anyDuplicated(p$patch_id)) stop("Invalid patch-history inputs.", call. = FALSE)
  if (!nrow(p)) return(list(patches = p, networks = data.frame(),
                             link_changes = data.frame()))
  parent <- seq_len(nrow(p))
  root <- function(i) {
    while (parent[[i]] != i) {
      parent[[i]] <<- parent[[parent[[i]]]]
      i <- parent[[i]]
    }
    i
  }
  from <- match(e$from_patch_id, p$patch_id)
  to <- match(e$to_patch_id, p$patch_id)
  if (anyNA(from) || anyNA(to)) stop("Links reference missing patches.", call. = FALSE)
  for (k in seq_along(from)) {
    a <- root(from[[k]]); b <- root(to[[k]])
    if (a != b) parent[[max(a, b)]] <- min(a, b)
  }
  roots <- vapply(seq_len(nrow(p)), root, integer(1))
  p$network_id <- paste0("history_", sprintf("%04d", match(roots, unique(roots))))
  parts <- split(seq_len(nrow(p)), p$network_id)
  physical_sd <- function(z, field) {
    if (!field %in% names(z)) return(NA_real_)
    value <- suppressWarnings(as.numeric(z[[field]]))
    if (sum(is.finite(value)) < 2L) return(NA_real_)
    stats::sd(value, na.rm = TRUE)
  }
  net <- do.call(rbind, lapply(parts, function(ii) {
    z <- p[ii, , drop = FALSE]
    incident <- e$from_patch_id %in% z$patch_id & e$to_patch_id %in% z$patch_id
    data.frame(network_id = z$network_id[[1L]],
      oldest_ma = max(z$time_ma), youngest_ma = min(z$time_ma),
      span_ma = max(z$time_ma) - min(z$time_ma),
      n_times = length(unique(z$time_ma)), n_patch_nodes = nrow(z),
      n_links = sum(incident),
      mean_patch_area_km2 = mean(z$area_km2),
      area_weighted_support = stats::weighted.mean(z$mean_support, z$area_km2),
      support_sd = if (nrow(z) > 1L) stats::sd(z$mean_support) else NA_real_,
      temperature_sd = physical_sd(z, "mean_temperature"),
      precipitation_sd = physical_sd(z, "mean_precipitation"),
      stringsAsFactors = FALSE)
  }))
  rownames(net) <- NULL
  change <- e
  if (nrow(e)) {
    old <- p[match(e$from_patch_id, p$patch_id), , drop = FALSE]
    young <- p[match(e$to_patch_id, p$patch_id), , drop = FALSE]
    duration <- old$time_ma - young$time_ma
    if (any(!is.finite(duration) | duration <= 0)) {
      stop("Patch links must point from older to younger time.", call. = FALSE)
    }
    change$duration_myr <- duration
    change$support_change_per_myr <-
      (young$mean_support - old$mean_support) / duration
    change$temperature_change_per_myr <- if (all(c("mean_temperature") %in% names(p)))
      (young$mean_temperature - old$mean_temperature) / duration else NA_real_
    change$precipitation_change_per_myr <- if (all(c("mean_precipitation") %in% names(p)))
      (young$mean_precipitation - old$mean_precipitation) / duration else NA_real_
  }
  list(patches = p, networks = net, link_changes = change)
}

#' Identify contraction episodes and candidate ecological refugia
#'
#' An adverse interval begins when total supported land area falls by at least
#' `contraction_fraction` between two successive time slices. A patch in the
#' younger slice is a candidate only if it has a carriage link from a supported
#' patch in the older slice. This is an ecological opportunity diagnostic,
#' not evidence of population survival. Geographic location mass is reported
#' separately and never multiplied by environmental support.
#'
#' @param patches Table with `patch_id`, `time_ma`, `area_km2`,
#'   `mean_support`, and optionally `location_mass`.
#' @param links Directed temporal patch links.
#' @param contraction_fraction Relative area loss threshold `(0,1)`.
#' @param min_location_mass Optional threshold for a separate, descriptive
#'   geography-supported evidence flag. Does not turn mass into occupancy.
#' @param min_patch_area_km2 Minimum candidate patch area.
#' @param require_later_link Require an outgoing carriage link from the patch
#'   to the next available slice after the contraction. This distinguishes
#'   one-step survivors from patches with at least two-step continuity.
#' @return List with `episodes` and `candidates`.
#' @export
hee_result_refugia_candidates <- function(patches, links,
                                           contraction_fraction = 0.25,
                                           min_location_mass = 0.1,
                                           min_patch_area_km2 = 0,
                                           require_later_link = FALSE) {
  p <- as.data.frame(patches); e <- as.data.frame(links)
  if (!all(c("patch_id", "time_ma", "area_km2", "mean_support") %in% names(p)) ||
      !all(c("from_patch_id", "to_patch_id", "time_from_ma", "time_to_ma") %in% names(e)) ||
      length(contraction_fraction) != 1L || !is.finite(contraction_fraction) ||
      contraction_fraction <= 0 || contraction_fraction >= 1 ||
      length(min_location_mass) != 1L || !is.finite(min_location_mass) ||
      min_location_mass < 0 || min_location_mass > 1) {
    stop("Invalid refugia inputs or thresholds.", call. = FALSE)
  }
  if (length(min_patch_area_km2) != 1L || !is.finite(min_patch_area_km2) ||
      min_patch_area_km2 < 0 || length(require_later_link) != 1L ||
      is.na(require_later_link)) {
    stop("Invalid patch-area or persistence requirement.", call. = FALSE)
  }
  if (!"location_mass" %in% names(p)) p$location_mass <- NA_real_
  times <- sort(unique(p$time_ma), decreasing = TRUE)
  episodes <- data.frame(time_from_ma = numeric(), time_to_ma = numeric(),
    area_from_km2 = numeric(), area_to_km2 = numeric(),
    fractional_loss = numeric())
  if (length(times) > 1L) for (k in seq_len(length(times) - 1L)) {
    before <- sum(p$area_km2[p$time_ma == times[[k]]])
    after <- sum(p$area_km2[p$time_ma == times[[k + 1L]]])
    if (before > 0 && (before - after) / before >= contraction_fraction) {
      episodes[nrow(episodes) + 1L, ] <- list(times[[k]], times[[k + 1L]],
        before, after, (before - after) / before)
    }
  }
  candidates <- data.frame()
  if (nrow(episodes) && nrow(e)) {
    selected <- merge(e, episodes[, c("time_from_ma", "time_to_ma")],
      by = c("time_from_ma", "time_to_ma"))
    if (nrow(selected)) {
      if ("carried_area_km2" %in% names(selected)) {
        selected <- selected[order(selected$carried_area_km2, decreasing = TRUE),
                             , drop = FALSE]
      }
      selected <- selected[!duplicated(selected$to_patch_id), , drop = FALSE]
      candidates <- merge(selected, p, by.x = "to_patch_id", by.y = "patch_id")
      candidates$later_carriage_link <- candidates$to_patch_id %in% e$from_patch_id
      candidates <- candidates[candidates$area_km2 >= min_patch_area_km2 &
        (!require_later_link | candidates$later_carriage_link), , drop = FALSE]
      candidates$geographic_evidence <- ifelse(is.na(candidates$location_mass),
        "not_assessed", ifelse(candidates$location_mass >= min_location_mass,
                                "location_mass_above_threshold", "low_location_mass"))
      candidates$interpretation <- "candidate_ecological_refugium_not_confirmed_occupancy"
    }
  }
  list(episodes = episodes, candidates = candidates)
}
