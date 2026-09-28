#' Define legacy geological-process scenario models
#'
#' Returns the earlier geological-process scenario family. The table is kept for
#' comparison, reporting, and backward-compatible case studies. It is not the
#' HmscEE 1.0 core process list: HmscEE exposes eight explicit biological
#' processes (environmental filtering, dispersal, colonisation, biotic
#' filtering, persistence, evolution, speciation, and extinction). Geological
#' and climatic variables are external drivers that condition those processes.
#'
#' @return A data.frame with model ids, process layers, and included mechanisms.
#' @export
#'
#' @examples
#' hee_define_geoprocess_models()
hee_define_geoprocess_models <- function() {
  data.frame(
    model_id = paste0("M", 0:8),
    model_name = c(
      "climate_only",
      "climate_arena",
      "climate_arena_network",
      "dynamic_colonisation",
      "colonisation_extinction_rescue",
      "lineage_time_extinction",
      "speciation_opportunity",
      "trait_niche_evolution",
      "earth_biota_feedback"
    ),
    mechanisms = c(
      "Climate suitability only",
      "Climate suitability plus geographic or habitat arena availability",
      "M1 plus isolation and connectivity networks",
      "M2 plus dynamic colonisation from source populations",
      "M3 plus local extinction and rescue effects",
      "M4 plus lineage origin and lineage extinction constraints",
      "M5 plus allopatric, founder-event, and in-situ speciation opportunity",
      "M6 plus trait and niche evolution diagnostics",
      "M7 plus optional community-to-evolution feedback terms"
    ),
    framework_role = rep("legacy_scenario_diagnostic_not_core_process", 9),
    core_process_mapping = c(
      "environmental_filtering",
      "Drivers constraining environmental_filtering/dispersal state",
      "dispersal proxy",
      "colonisation state transition",
      "environmental_filtering/dispersal/biotic_filtering/persistence local state transition",
      "speciation/extinction lineage-time constraint",
      "speciation opportunity proxy, not speciation-rate estimation",
      "evolution feeding environmental_filtering",
      "External feedback hypothesis only"
    ),
    interpretation_boundary = c(
      "Climate suitability is not a complete past distribution.",
      "Geographic/habitat arena is a driver or hard state, not a biological probability.",
      "Connectivity requires calibration before causal dispersal claims.",
      "Colonisation is scenario-based unless calibrated with independent data.",
      "Local extinction/rescue proxies are not global lineage extinction.",
      "Extant-only trees do not recover unsampled extinct lineages.",
      "Opportunity indices cannot be reported as true speciation rates.",
      "Trait-niche diagnostics depend on response-evolution assumptions.",
      "Feedback requires independent temporal or experimental evidence."
    ),
    stringsAsFactors = FALSE
  )
}

#' Detect geological landscape events through time
#'
#' Identifies simple arena creation/loss and transformation events from a
#' region- or cell-level time table. This is a diagnostic layer: emergence,
#' submergence, expansion, contraction, uplift, subsidence, habitat gain/loss,
#' and component changes are inferred from user-supplied palaeogeographic
#' columns rather than estimated de novo.
#'
#' @param landscape Data.frame with a spatial unit, `time_ma`, and one or more
#'   state columns such as land, habitat, area, elevation, or component.
#' @param unit_col Spatial unit column. If `NULL`, the function uses `region`
#'   when present, otherwise `cell_id`.
#' @param time_col Time column in Ma. Larger values are older.
#' @param land_col Optional 0/1 geographic existence column.
#' @param habitat_col Optional 0/1 habitat availability column.
#' @param area_col Optional area column.
#' @param elevation_col Optional elevation column.
#' @param component_col Optional network/component id column.
#' @param area_rel_threshold Relative area-change threshold for expansion or
#'   contraction.
#' @param elevation_threshold Elevation-change threshold for uplift/subsidence.
#' @param uncertainty Default uncertainty value attached to inferred events.
#'
#' @return A landscape event table with event id, interval, event type, unit,
#'   area change, elevation change, and evidence.
#' @export
#'
#' @examples
#' x <- data.frame(region = "A", time_ma = c(10, 5, 0),
#'   land = c(0, 1, 1), area_km2 = c(0, 100, 180), elev_m = c(0, 20, 260))
#' hee_landscape_events(x, area_col = "area_km2", elevation_col = "elev_m")
hee_landscape_events <- function(landscape,
                                 unit_col = NULL,
                                 time_col = "time_ma",
                                 land_col = "land",
                                 habitat_col = NULL,
                                 area_col = NULL,
                                 elevation_col = NULL,
                                 component_col = NULL,
                                 area_rel_threshold = 0.1,
                                 elevation_threshold = 100,
                                 uncertainty = NA_real_) {
  x <- as.data.frame(landscape)
  .require_cols(x, time_col, "landscape")
  x[[time_col]] <- .hee_numeric_values(x[[time_col]],
                                       paste0("landscape$", time_col))
  .hee_validate_time_values(x[[time_col]], paste0("landscape$", time_col))
  if (is.null(unit_col)) {
    unit_col <- if ("region" %in% names(x)) "region" else "cell_id"
  }
  .require_cols(x, unit_col, "landscape")
  .hee_check_unique_keys(x, c(unit_col, time_col), "landscape")

  optional <- c(land_col, habitat_col, area_col, elevation_col, component_col)
  optional <- optional[!is.null(optional) & nzchar(optional)]
  missing_optional <- setdiff(optional, names(x))
  if (length(missing_optional) > 0L) {
    stop("landscape is missing requested column(s): ",
         paste(missing_optional, collapse = ", "), call. = FALSE)
  }
  if (length(intersect(optional, names(x))) == 0L) {
    stop("Provide at least one state column such as land, habitat, area, elevation, or component.",
         call. = FALSE)
  }

  x <- x[order(x[[unit_col]], x[[time_col]], decreasing = TRUE), , drop = FALSE]
  out <- list()
  event_n <- 0L

  add_event <- function(prev, curr, type, area_change = NA_real_,
                        elev_change = NA_real_, barrier_change = NA_real_,
                        evidence = NA_character_) {
    event_n <<- event_n + 1L
    out[[event_n]] <<- data.frame(
      event_id = paste0("landscape_event_", event_n),
      time_start_ma = prev[[time_col]][1],
      time_end_ma = curr[[time_col]][1],
      event_time_ma = curr[[time_col]][1],
      event_type = type,
      unit_id = as.character(curr[[unit_col]][1]),
      region_from = as.character(prev[[unit_col]][1]),
      region_to = as.character(curr[[unit_col]][1]),
      area_change_km2 = area_change,
      elevation_change_m = elev_change,
      barrier_change = barrier_change,
      uncertainty = uncertainty,
      evidence = evidence,
      stringsAsFactors = FALSE
    )
  }

  for (unit in unique(as.character(x[[unit_col]]))) {
    z <- x[as.character(x[[unit_col]]) == unit, , drop = FALSE]
    z <- z[order(z[[time_col]], decreasing = TRUE), , drop = FALSE]
    if (nrow(z) < 2L) next
    for (i in 2:nrow(z)) {
      prev <- z[i - 1L, , drop = FALSE]
      curr <- z[i, , drop = FALSE]

      if (!is.null(land_col) && land_col %in% names(x)) {
        lp <- .hee_binary(prev[[land_col]][1])
        lc <- .hee_binary(curr[[land_col]][1])
        if (isTRUE(lp == 0 && lc == 1)) {
          add_event(prev, curr, "emergence", evidence = land_col)
        } else if (isTRUE(lp == 1 && lc == 0)) {
          add_event(prev, curr, "submergence", evidence = land_col)
        }
      }

      if (!is.null(habitat_col) && habitat_col %in% names(x)) {
        hp <- .hee_binary(prev[[habitat_col]][1])
        hc <- .hee_binary(curr[[habitat_col]][1])
        if (isTRUE(hp == 0 && hc == 1)) {
          add_event(prev, curr, "habitat_gain", evidence = habitat_col)
        } else if (isTRUE(hp == 1 && hc == 0)) {
          add_event(prev, curr, "habitat_loss", evidence = habitat_col)
        }
      }

      if (!is.null(area_col) && area_col %in% names(x)) {
        ap <- suppressWarnings(as.numeric(prev[[area_col]][1]))
        ac <- suppressWarnings(as.numeric(curr[[area_col]][1]))
        if (!is.na(ap) && !is.na(ac)) {
          delta <- ac - ap
          denom <- max(abs(ap), .Machine$double.eps)
          rel <- delta / denom
          if (rel >= area_rel_threshold) {
            add_event(prev, curr, "area_expansion", area_change = delta,
                      evidence = area_col)
          } else if (rel <= -area_rel_threshold) {
            add_event(prev, curr, "area_contraction", area_change = delta,
                      evidence = area_col)
          }
        }
      }

      if (!is.null(elevation_col) && elevation_col %in% names(x)) {
        ep <- suppressWarnings(as.numeric(prev[[elevation_col]][1]))
        ec <- suppressWarnings(as.numeric(curr[[elevation_col]][1]))
        if (!is.na(ep) && !is.na(ec)) {
          delta <- ec - ep
          if (delta >= elevation_threshold) {
            add_event(prev, curr, "uplift", elev_change = delta,
                      evidence = elevation_col)
          } else if (delta <= -elevation_threshold) {
            add_event(prev, curr, "subsidence_or_erosion", elev_change = delta,
                      evidence = elevation_col)
          }
        }
      }

      if (!is.null(component_col) && component_col %in% names(x)) {
        cp <- as.character(prev[[component_col]][1])
        cc <- as.character(curr[[component_col]][1])
        if (!is.na(cp) && !is.na(cc) && !identical(cp, cc)) {
          add_event(prev, curr, "component_change", barrier_change = NA_real_,
                    evidence = component_col)
        }
      }
    }
  }

  if (length(out) == 0L) {
    return(data.frame(
      event_id = character(), time_start_ma = numeric(),
      time_end_ma = numeric(), event_time_ma = numeric(),
      event_type = character(), unit_id = character(),
      region_from = character(), region_to = character(),
      area_change_km2 = numeric(), elevation_change_m = numeric(),
      barrier_change = numeric(), uncertainty = numeric(),
      evidence = character(), stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, out)
}

#' Compute land or habitat age through deep time
#'
#' Calculates how long each cell, region, island, lake, or other palaeogeographic
#' arena has persisted as usable habitat. Time must be in Ma, with larger values
#' older. Newly emerged habitat starts at age zero; habitat that remains present
#' accumulates the true interval length between time slices; absent habitat gets
#' `NA`. This is the "land age" term in the Dynamic Earth-Biota Assembly
#' framework, not a fitted biological parameter.
#'
#' @param landscape Data.frame with a spatial unit, `time_ma`, and land/habitat
#'   columns.
#' @param unit_col Spatial unit column. Defaults to `region` when present, then
#'   `cell_id`.
#' @param time_col Time column in Ma.
#' @param land_col Geographic existence column.
#' @param habitat_col Optional habitat availability column. When supplied, the
#'   arena exists only when both land and habitat are positive.
#' @param initial_age_myr Age assigned to the oldest available present arena.
#' @param age_col Output age column.
#' @return `landscape` with `geographic_existence` and land-age columns.
#' @export
#'
#' @examples
#' x <- data.frame(cell_id = "c1", time_ma = c(10, 5, 0), land = c(0, 1, 1))
#' hee_land_age(x)
hee_land_age <- function(landscape,
                         unit_col = NULL,
                         time_col = "time_ma",
                         land_col = "land",
                         habitat_col = NULL,
                         initial_age_myr = 0,
                         age_col = "land_age_myr") {
  x <- as.data.frame(landscape)
  .require_cols(x, time_col, "landscape")
  x[[time_col]] <- .hee_numeric_values(x[[time_col]],
                                       paste0("landscape$", time_col))
  .hee_validate_time_values(x[[time_col]], paste0("landscape$", time_col))
  if (is.null(unit_col)) {
    unit_col <- if ("region" %in% names(x)) "region" else "cell_id"
  }
  .require_cols(x, c(unit_col, land_col), "landscape")
  if (!is.null(habitat_col)) .require_cols(x, habitat_col, "landscape")
  .hee_check_unique_keys(x, c(unit_col, time_col), "landscape")

  x$.row_id <- seq_len(nrow(x))
  x$geographic_existence <- .hee_binary(x[[land_col]])
  if (!is.null(habitat_col)) {
    x$geographic_existence <- x$geographic_existence * .hee_binary(x[[habitat_col]])
  }
  x[[age_col]] <- rep(NA_real_, nrow(x))
  if (nrow(x) == 0L) {
    if (!"land_age_ma" %in% names(x)) x$land_age_ma <- x[[age_col]]
    if (!"land_age_myr" %in% names(x)) x$land_age_myr <- x[[age_col]]
    return(x)
  }

  parts <- split(x, as.character(x[[unit_col]]), drop = TRUE)
  out <- lapply(parts, function(z) {
    z <- z[order(z[[time_col]], decreasing = TRUE), , drop = FALSE]
    age <- rep(NA_real_, nrow(z))
    for (i in seq_len(nrow(z))) {
      exists_now <- isTRUE(z$geographic_existence[i] > 0)
      if (!exists_now) {
        age[i] <- NA_real_
      } else if (i == 1L) {
        age[i] <- initial_age_myr
      } else if (!isTRUE(z$geographic_existence[i - 1L] > 0)) {
        age[i] <- 0
      } else {
        age[i] <- age[i - 1L] + abs(z[[time_col]][i - 1L] - z[[time_col]][i])
      }
    }
    z[[age_col]] <- age
    z
  })
  out <- do.call(rbind, out)
  if (!"land_age_ma" %in% names(out)) out$land_age_ma <- out[[age_col]]
  if (!"land_age_myr" %in% names(out)) out$land_age_myr <- out[[age_col]]
  out <- out[order(out$.row_id), setdiff(names(out), ".row_id"), drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Compute ecological opportunity from arena change and novelty
#'
#' Builds the diagnostic opportunity term described in the geological-process
#' workflow. Opportunity combines positive area gain, environmental
#' heterogeneity, no-analog or novel climate distance, and young-land opportunity
#' decay. When explicit geography or habitat-availability columns are present,
#' rows with absent or unknown geography are masked to zero. The result is a
#' transparent index used for scenario comparison, not a direct estimate of true
#' speciation or colonisation rates.
#'
#' @param cell_table Cell/region-time table.
#' @param env_cols Optional environmental columns used for heterogeneity and
#'   novelty when explicit columns are absent.
#' @param region_col Region/grouping column.
#' @param time_col Time column in Ma.
#' @param area_col Area column.
#' @param land_age_col Land-age column. Defaults to `land_age_myr` or
#'   `land_age_ma` when present.
#' @param heterogeneity_col Optional precomputed heterogeneity column.
#' @param novelty_col Optional precomputed novelty/no-analog column.
#' @param weights Named weights for `area`, `heterogeneity`, `novelty`, and
#'   `young_land`.
#' @param tau Young-land opportunity decay time in Myr.
#' @return `cell_table` with opportunity components and
#'   `ecological_opportunity`.
#' @export
#'
#' @examples
#' x <- data.frame(region = "A", cell_id = 1:3, time_ma = c(10, 5, 0),
#'   cell_area_km2 = c(0, 10, 20), land = c(0, 1, 1), bio1 = c(1, 2, 3))
#' hee_ecological_opportunity(hee_land_age(x), env_cols = "bio1")
hee_ecological_opportunity <- function(cell_table,
                                       env_cols = NULL,
                                       region_col = "region",
                                       time_col = "time_ma",
                                       area_col = "cell_area_km2",
                                       land_age_col = NULL,
                                       heterogeneity_col = NULL,
                                       novelty_col = NULL,
                                       weights = c(area = 1, heterogeneity = 1,
                                                   novelty = 1, young_land = 1),
                                       tau = 5) {
  x <- as.data.frame(cell_table)
  .require_cols(x, time_col, "cell_table")
  x[[time_col]] <- .hee_numeric_values(x[[time_col]],
                                       paste0("cell_table$", time_col))
  .hee_validate_time_values(x[[time_col]], paste0("cell_table$", time_col))
  if ("cell_id" %in% names(x)) {
    .hee_check_unique_keys(x, c("cell_id", time_col), "cell_table")
  }
  if (!region_col %in% names(x)) region_col <- if ("cell_id" %in% names(x)) "cell_id" else time_col
  if (is.null(land_age_col)) {
    land_age_col <- intersect(c("land_age_myr", "land_age_ma"), names(x))[1]
  }
  env_cols <- intersect(env_cols %||% character(), names(x))
  x$.row_id <- seq_len(nrow(x))
  if (nrow(x) == 0L) {
    x$area_gain_positive <- numeric()
    x$area_gain_index <- numeric()
    x$heterogeneity_index <- numeric()
    x$novelty_index <- numeric()
    x$young_land_index <- numeric()
    x$ecological_opportunity <- numeric()
    return(x[, setdiff(names(x), ".row_id"), drop = FALSE])
  }

  # Area-gain component, computed at region x time and mapped back to cells.
  if (area_col %in% names(x)) {
    area_tab <- stats::aggregate(x[[area_col]], x[, c(region_col, time_col), drop = FALSE],
                                 sum, na.rm = TRUE)
    names(area_tab)[ncol(area_tab)] <- ".area"
    parts <- split(area_tab, as.character(area_tab[[region_col]]), drop = TRUE)
    area_tab <- do.call(rbind, lapply(parts, function(z) {
      z <- z[order(z[[time_col]], decreasing = TRUE), , drop = FALSE]
      z$area_gain_positive <- c(0, pmax(diff(z$.area), 0))
      z
    }))
    x <- merge(x, area_tab[, c(region_col, time_col, "area_gain_positive"), drop = FALSE],
               by = c(region_col, time_col), all.x = TRUE, sort = FALSE)
  } else {
    x$area_gain_positive <- 0
  }
  x$area_gain_index <- .hee_rescale01(x$area_gain_positive)

  if (!is.null(heterogeneity_col) && heterogeneity_col %in% names(x)) {
    x$heterogeneity_index <- .hee_rescale01(x[[heterogeneity_col]])
  } else if (length(env_cols) > 0L) {
    het <- do.call(rbind, lapply(split(x, interaction(x[[region_col]], x[[time_col]], drop = TRUE)),
                                 function(z) {
      vals <- vapply(env_cols, function(nm) stats::sd(z[[nm]], na.rm = TRUE), numeric(1))
      data.frame(region = z[[region_col]][1], time_ma = z[[time_col]][1],
                 .hee_habitat_heterogeneity = mean(vals, na.rm = TRUE),
                 stringsAsFactors = FALSE)
    }))
    names(het)[1:2] <- c(region_col, time_col)
    x <- merge(x, het, by = c(region_col, time_col), all.x = TRUE, sort = FALSE)
    x$heterogeneity_index <- .hee_rescale01(x$.hee_habitat_heterogeneity)
  } else {
    x$heterogeneity_index <- 0
  }

  if (!is.null(novelty_col) && novelty_col %in% names(x)) {
    x$novelty_index <- .hee_rescale01(x[[novelty_col]])
  } else if (length(env_cols) > 0L) {
    present_time <- if (any(x[[time_col]] == 0, na.rm = TRUE)) 0 else min(x[[time_col]], na.rm = TRUE)
    ref <- x[x[[time_col]] == present_time, env_cols, drop = FALSE]
    mu <- vapply(ref, mean, numeric(1), na.rm = TRUE)
    sig <- vapply(ref, stats::sd, numeric(1), na.rm = TRUE)
    sig[!is.finite(sig) | sig == 0] <- 1
    z <- sweep(sweep(as.matrix(x[, env_cols, drop = FALSE]), 2, mu, "-"), 2, sig, "/")
    x$novel_climate_index <- .hee_rescale01(rowMeans(abs(z), na.rm = TRUE))
    x$novelty_index <- x$novel_climate_index
  } else {
    x$novelty_index <- 0
  }

  if (!is.na(land_age_col) && land_age_col %in% names(x)) {
    age <- suppressWarnings(as.numeric(x[[land_age_col]]))
    x$young_land_index <- exp(-pmax(age, 0) / max(tau, .Machine$double.eps))
    x$young_land_index[is.na(age)] <- 0
  } else {
    x$young_land_index <- 0
  }

  w <- c(area = 1, heterogeneity = 1, novelty = 1, young_land = 1)
  if (!is.null(weights)) {
    if (is.null(names(weights))) {
      w[seq_len(min(length(w), length(weights)))] <-
        weights[seq_len(min(length(w), length(weights)))]
    } else {
      hit <- intersect(names(weights), names(w))
      unknown <- setdiff(names(weights), names(w))
      if (length(unknown) > 0L) {
        warning("Ignoring unknown ecological-opportunity weight(s): ",
                paste(unknown, collapse = ", "), call. = FALSE)
      }
      w[hit] <- weights[hit]
    }
  }
  w[is.na(w) | w < 0] <- 0
  denom <- max(sum(w), .Machine$double.eps)
  x$ecological_opportunity <- .hee_clip01(
    (w["area"] * x$area_gain_index +
       w["heterogeneity"] * x$heterogeneity_index +
       w["novelty"] * x$novelty_index +
       w["young_land"] * x$young_land_index) / denom
  )
  geo_cols <- intersect(c("geographic_existence", "land", "land_mask_dem"),
                        names(x))
  if (length(geo_cols) > 0L) {
    geo_mask <- rep(1, nrow(x))
    for (nm in geo_cols) {
      v <- suppressWarnings(as.numeric(x[[nm]]))
      v[is.na(v)] <- 0
      geo_mask <- geo_mask * as.numeric(.hee_clip01(v) > 0)
    }
    x$geographic_opportunity_mask <- .hee_clip01(geo_mask)
    x$ecological_opportunity <- .hee_clip01(x$ecological_opportunity *
                                              x$geographic_opportunity_mask)
  }
  if ("habitat_availability" %in% names(x)) {
    hab <- suppressWarnings(as.numeric(x$habitat_availability))
    hab[is.na(hab)] <- 0
    x$habitat_opportunity_weight <- .hee_clip01(hab)
    x$ecological_opportunity <- .hee_clip01(x$ecological_opportunity *
                                              x$habitat_opportunity_weight)
  }
  x <- x[order(x$.row_id),
         setdiff(names(x), c(".row_id", ".hee_habitat_heterogeneity")),
         drop = FALSE]
  rownames(x) <- NULL
  x
}

#' Compute structural geological connectivity
#'
#' Structural connectivity measures the physical connection between regions,
#' cells, islands, lakes, or habitat patches. It combines route availability,
#' palaeodistance decay, and barrier strength. Species traits are deliberately
#' not used here; add them with `hee_functional_connectivity()`.
#'
#' @param pairs Region-pair/time table.
#' @param distance_col Palaeodistance column. Values must be finite
#'   non-negative distances or path costs.
#' @param barrier_col Barrier-strength column from 0 to 1. If this column is
#'   omitted, barriers are neutral; if supplied, explicit `NA` values are treated
#'   conservatively as strong barriers.
#' @param route_col Optional route-open column from 0 to 1. If this column is
#'   omitted, routes are neutral; if supplied, explicit `NA` values are treated
#'   conservatively as closed routes.
#' @param alpha Distance-decay coefficient.
#' @param barrier_weight Barrier-decay coefficient.
#' @return `pairs` with `structural_connectivity`.
#' @export
#'
#' @examples
#' x <- data.frame(from_region = "A", to_region = "B", time_ma = 0,
#'   paleodistance = 10, barrier_strength = 0.5)
#' hee_structural_connectivity(x, alpha = 0.01)
hee_structural_connectivity <- function(pairs,
                                        distance_col = "paleodistance",
                                        barrier_col = "barrier_strength",
                                        route_col = "route_open",
                                        alpha = 0.01,
                                        barrier_weight = 1) {
  x <- as.data.frame(pairs)
  .require_cols(x, distance_col, "pairs")
  if (all(c("from_region", "to_region", "time_ma") %in% names(x))) {
    x$time_ma <- .hee_numeric_values(x$time_ma, "pairs$time_ma")
    .hee_validate_time_values(x$time_ma, "pairs$time_ma")
    .hee_check_unique_keys(x, c("from_region", "to_region", "time_ma"),
                           "pairs")
  }
  d <- suppressWarnings(as.numeric(x[[distance_col]]))
  if (any(is.finite(d) & d < 0, na.rm = TRUE)) {
    stop("`", distance_col, "` must be non-negative palaeodistance values.",
         call. = FALSE)
  }
  d[is.na(d)] <- Inf
  has_barrier_col <- barrier_col %in% names(x)
  has_route_col <- route_col %in% names(x)
  barrier <- if (has_barrier_col) .hee_clip01(x[[barrier_col]]) else
    rep(0, nrow(x))
  route <- if (has_route_col) .hee_clip01(x[[route_col]]) else rep(1, nrow(x))
  barrier[is.na(barrier)] <- if (has_barrier_col) 1 else 0
  route[is.na(route)] <- if (has_route_col) 0 else 1
  x$structural_connectivity <- .hee_clip01(route * exp(-alpha * pmax(d, 0)) *
                                             exp(-barrier_weight * barrier))
  x
}

#' Compute species-specific functional connectivity
#'
#' Functional connectivity translates the same geological route into
#' species-specific permeability using dispersal traits, palaeodistance, and
#' barrier strength. This is a scenario/diagnostic calculation unless its
#' coefficients are calibrated externally.
#'
#' @param connectivity Pair/time table, usually from
#'   `hee_structural_connectivity()`.
#' @param traits Optional trait table. When omitted, one non-species-specific
#'   functional layer is returned.
#' @param dispersal_trait Trait column. Values must be finite non-negative
#'   dispersal abilities or distances.
#' @param kappa Named coefficients `k0`, `trait`, `distance`, and `barrier`.
#' @param expand_species Logical. If `connectivity` has no species column and
#'   `traits` has multiple species, `TRUE` explicitly expands each route to all
#'   trait species. Set `FALSE` to require pre-species-expanded connectivity.
#' @param species_col Species column.
#' @param distance_col,barrier_col Column names. `distance_col` values must be
#'   finite non-negative palaeodistances or path costs. If `barrier_col` is
#'   omitted, barriers are neutral; if supplied, explicit `NA` values are treated
#'   conservatively as strong barriers.
#' @return Connectivity table with `functional_connectivity`; when structural
#'   connectivity is present, a combined `connectivity` column is also returned.
#'   If `traits = NULL`, the functional layer is neutral
#'   (`functional_connectivity = 1`) so distance and barriers are not counted a
#'   second time after [hee_structural_connectivity()].
#' @export
hee_functional_connectivity <- function(connectivity,
                                        traits = NULL,
                                        dispersal_trait = "dispersal",
                                        kappa = c(k0 = 0, trait = 1,
                                                  distance = 1, barrier = 1),
                                        expand_species = TRUE,
                                        species_col = "species",
                                        distance_col = "paleodistance",
                                        barrier_col = "barrier_strength") {
  x <- as.data.frame(connectivity)
  has_barrier_col <- barrier_col %in% names(x)
  if (!distance_col %in% names(x)) x[[distance_col]] <- 0
  if (!barrier_col %in% names(x)) x[[barrier_col]] <- 0
  x$.pair_id <- seq_len(nrow(x))
  if (!is.null(traits)) {
    tr <- as.data.frame(traits)
    .require_cols(tr, c(species_col, dispersal_trait), "traits")
    tr <- tr[, c(species_col, dispersal_trait), drop = FALSE]
    .hee_check_unique_keys(tr, species_col, "traits")
    tr_disp <- suppressWarnings(as.numeric(tr[[dispersal_trait]]))
    if (any(!is.finite(tr_disp) | tr_disp < 0, na.rm = TRUE)) {
      stop("traits$", dispersal_trait,
           " must contain finite non-negative dispersal values.",
           call. = FALSE)
    }
    tr[[dispersal_trait]] <- tr_disp
    if (species_col %in% names(x)) {
      x <- merge(x, tr, by = species_col, all.x = TRUE, sort = FALSE)
      missing_trait_species <- unique(as.character(x[[species_col]])[is.na(x[[dispersal_trait]])])
      missing_trait_species <- missing_trait_species[!is.na(missing_trait_species)]
      if (length(missing_trait_species) > 0L) {
        stop("traits are missing dispersal values for species: ",
             paste(missing_trait_species, collapse = ", "), call. = FALSE)
      }
    } else {
      if (!isTRUE(expand_species) && nrow(tr) > 1L) {
        stop("connectivity has no `", species_col,
             "` column. Set expand_species = TRUE to explicitly expand routes ",
             "to all trait species, or provide species-specific connectivity.",
             call. = FALSE)
      }
      x <- merge(x, tr, by = NULL, all = TRUE, sort = FALSE)
    }
  } else {
    x[[dispersal_trait]] <- NA_real_
    x$functional_connectivity <- 1
    if ("structural_connectivity" %in% names(x)) {
      structural <- .hee_clip01(x$structural_connectivity)
      structural[is.na(structural)] <- 0
      x$connectivity <- structural
    } else {
      x$connectivity <- 1
    }
    return(x[order(x$.pair_id), setdiff(names(x), ".pair_id"), drop = FALSE])
  }
  cf <- .hee_named_coefficients(kappa, c("k0", "trait", "distance", "barrier"))
  raw_distance <- suppressWarnings(as.numeric(x[[distance_col]]))
  if (any(is.finite(raw_distance) & raw_distance < 0, na.rm = TRUE)) {
    stop("`", distance_col, "` must be non-negative palaeodistance values.",
         call. = FALSE)
  }
  d <- .hee_positive_index(raw_distance)
  b <- .hee_clip01(x[[barrier_col]])
  b[is.na(b)] <- if (has_barrier_col) 1 else 0
  disp <- suppressWarnings(as.numeric(x[[dispersal_trait]]))
  fallback <- stats::median(disp[is.finite(disp)], na.rm = TRUE)
  if (!is.finite(fallback)) fallback <- 0
  disp[!is.finite(disp)] <- fallback
  trait <- .hee_positive_index(disp)
  eta <- cf["k0"] + cf["trait"] * trait - cf["distance"] * d - cf["barrier"] * b
  x$functional_connectivity <- .hee_clip01(stats::plogis(eta))
  if ("structural_connectivity" %in% names(x)) {
    structural <- .hee_clip01(x$structural_connectivity)
    structural[is.na(structural)] <- 0
    x$connectivity <- .hee_clip01(structural * x$functional_connectivity)
  }
  x[order(x$.pair_id), setdiff(names(x), ".pair_id"), drop = FALSE]
}

#' Compute climatic corridor connectivity
#'
#' Climatic connectivity represents whether a physical route is environmentally
#' traversable. Supply a corridor/environmental distance column, a corridor
#' suitability column, or a precomputed climate-connectivity column. The function
#' does not infer least-cost paths unless the user supplies those path summaries.
#'
#' @param connectivity Pair/species/time table.
#' @param environment_distance_col Environmental distance or path cost column.
#'   Values must be finite non-negative distances when the column is supplied.
#' @param corridor_suitability_col Corridor suitability column from 0 to 1.
#' @param climate_col Output/precomputed climate-connectivity column. If a
#'   climate/corridor column is explicitly supplied, missing values are treated
#'   as no climate-corridor evidence (`0`). If no climate information is
#'   supplied at all, the climate component is neutral (`1`).
#' @param phi Environmental-distance decay.
#' @return Connectivity table with `climate_connectivity`; when structural or
#'   functional components are present, `connectivity` is updated.
#' @export
hee_climatic_connectivity <- function(connectivity,
                                      environment_distance_col = "environment_distance",
                                      corridor_suitability_col = "corridor_suitability",
                                      climate_col = "climate_connectivity",
                                      phi = 1) {
  x <- as.data.frame(connectivity)
  if (!is.null(environment_distance_col) && environment_distance_col %in% names(x)) {
    ed <- suppressWarnings(as.numeric(x[[environment_distance_col]]))
    if (any(is.finite(ed) & ed < 0, na.rm = TRUE)) {
      stop("`", environment_distance_col,
           "` must be non-negative environmental-distance values.",
           call. = FALSE)
    }
    ed[is.na(ed)] <- Inf
    x[[climate_col]] <- .hee_clip01(exp(-phi * pmax(ed, 0)))
    x[[climate_col]][is.na(x[[climate_col]])] <- 0
  } else if (!is.null(corridor_suitability_col) && corridor_suitability_col %in% names(x)) {
    x[[climate_col]] <- .hee_clip01(x[[corridor_suitability_col]])
    x[[climate_col]][is.na(x[[climate_col]])] <- 0
  } else if (climate_col %in% names(x)) {
    x[[climate_col]] <- .hee_clip01(x[[climate_col]])
    x[[climate_col]][is.na(x[[climate_col]])] <- 0
  } else {
    x[[climate_col]] <- 1
  }
  base <- rep(1, nrow(x))
  for (nm in intersect(c("structural_connectivity", "functional_connectivity"), names(x))) {
    v <- .hee_clip01(x[[nm]])
    v[is.na(v)] <- 1
    base <- base * v
  }
  x$connectivity <- .hee_clip01(base * x[[climate_col]])
  x
}

#' Build a time-stratified geological connectivity cube
#'
#' Builds structural, functional, climate, and total connectivity among regions
#' for each time slice. Structural connectivity is based on palaeodistance and
#' optional barrier passability. Functional connectivity may be species-specific
#' when dispersal traits are supplied. Climate connectivity is taken from a
#' supplied corridor table or defaults to one.
#'
#' @param distances Region-pair table or square distance matrix. A table should
#'   contain `from_region`, `to_region`, `time_ma`, and a distance column.
#' @param landscape Optional region/time table with `lon` and `lat`; used to
#'   compute pairwise distances when `distances` is omitted.
#' @param traits Optional species trait table with `species` and a dispersal
#'   distance column.
#' @param barrier Optional table with pair/time barrier strength or passability.
#'   If omitted, barriers are neutral. If only barrier strength is supplied,
#'   passability is derived as `1 - barrier_strength`. If passability is
#'   explicitly supplied, unmatched or missing passability is conservative
#'   (`0`), and unmatched or missing barrier strength is treated as a strong
#'   barrier.
#' @param climate Optional table with pair/time climate connectivity.
#' @param times Optional time vector used with matrix distances or landscape.
#' @param distance_col Distance column name. Values must be finite
#'   non-negative palaeodistances or path costs.
#' @param distance_scale Scale parameter for structural distance decay.
#' @param dispersal_col Species dispersal-distance trait column. Values must be
#'   finite non-negative dispersal distances. Zero means no movement across
#'   positive palaeodistance.
#' @param barrier_col Barrier strength column, from 0 to 1.
#' @param passability_col Optional passability column, from 0 to 1.
#' @param climate_col Climate connectivity column.
#' @param expand_species Logical. If traits are supplied, `TRUE` explicitly
#'   expands each region pair to each trait species.
#' @param species_col,from_col,to_col,time_col Column names.
#'
#' @return A connectivity cube table.
#' @export
#'
#' @examples
#' d <- data.frame(from_region = "A", to_region = "B", time_ma = 0,
#'   paleodistance = 2, barrier_strength = 0.2)
#' tr <- data.frame(species = "sp1", dispersal_distance = 4)
#' hee_build_connectivity_cube(d, traits = tr, distance_scale = 5)
hee_build_connectivity_cube <- function(distances = NULL,
                                        landscape = NULL,
                                        traits = NULL,
                                        barrier = NULL,
                                        climate = NULL,
                                        times = NULL,
                                        distance_col = "paleodistance",
                                        distance_scale = 1,
                                        dispersal_col = "dispersal_distance",
                                        barrier_col = "barrier_strength",
                                        passability_col = "barrier_passability",
                                        climate_col = "climate_connectivity",
                                        expand_species = TRUE,
                                        species_col = "species",
                                        from_col = "from_region",
                                        to_col = "to_region",
                                        time_col = "time_ma") {
  pair <- .hee_connectivity_pairs(distances, landscape, times, distance_col,
                                  from_col, to_col, time_col)

  has_pair_barrier <- barrier_col %in% names(pair)
  has_pair_passability <- passability_col %in% names(pair)
  if (!barrier_col %in% names(pair)) pair[[barrier_col]] <- 0
  if (!passability_col %in% names(pair)) pair[[passability_col]] <- NA_real_

  has_external_passability <- FALSE
  if (!is.null(barrier)) {
    b <- as.data.frame(barrier)
    has_external_passability <- passability_col %in% names(b)
    keys <- intersect(c(from_col, to_col, time_col), names(b))
    .require_cols(b, c(from_col, to_col, time_col), "barrier")
    .hee_check_unique_keys(b, c(from_col, to_col, time_col), "barrier")
    keep <- unique(c(keys, intersect(c(barrier_col, passability_col), names(b))))
    pair <- merge(pair, b[, keep, drop = FALSE], by = keys, all.x = TRUE,
                  sort = FALSE, suffixes = c("", ".barrier"))
    for (nm in c(barrier_col, passability_col)) {
      alt <- paste0(nm, ".barrier")
      if (alt %in% names(pair)) {
        pair[[nm]] <- ifelse(is.na(pair[[alt]]), pair[[nm]], pair[[alt]])
        pair[[alt]] <- NULL
      }
    }
  }

  barrier_default <- if (has_pair_barrier || !is.null(barrier)) 1 else 0
  pair[[barrier_col]][is.na(pair[[barrier_col]])] <- barrier_default
  pair[[barrier_col]] <- .hee_clip01(pair[[barrier_col]])
  pass <- pair[[passability_col]]
  pass_missing <- is.na(pass)
  pass[pass_missing] <- if (has_pair_passability || has_external_passability) {
    0
  } else {
    1 - pair[[barrier_col]][pass_missing]
  }
  pair[[passability_col]] <- .hee_clip01(pass)

  has_pair_climate <- climate_col %in% names(pair)
  if (!has_pair_climate) pair[[climate_col]] <- if (is.null(climate)) 1 else NA_real_
  if (!is.null(climate)) {
    clim <- as.data.frame(climate)
    .require_cols(clim, c(from_col, to_col, time_col, climate_col), "climate")
    .hee_check_unique_keys(clim, c(from_col, to_col, time_col), "climate")
    pair <- merge(pair, clim[, c(from_col, to_col, time_col, climate_col), drop = FALSE],
                  by = c(from_col, to_col, time_col), all.x = TRUE,
                  sort = FALSE, suffixes = c("", ".climate"))
    alt <- paste0(climate_col, ".climate")
    if (alt %in% names(pair)) {
      pair[[climate_col]] <- ifelse(is.na(pair[[alt]]), pair[[climate_col]],
                                    pair[[alt]])
      pair[[alt]] <- NULL
    }
  }
  if (is.null(climate) && !has_pair_climate) {
    pair[[climate_col]][is.na(pair[[climate_col]])] <- 1
  } else {
    pair[[climate_col]][is.na(pair[[climate_col]])] <- 0
  }
  pair[[climate_col]] <- .hee_clip01(pair[[climate_col]])
  dist_raw <- suppressWarnings(as.numeric(pair[[distance_col]]))
  if (any(is.finite(dist_raw) & dist_raw < 0, na.rm = TRUE)) {
    stop("`", distance_col, "` must be non-negative palaeodistance values.",
         call. = FALSE)
  }
  pair[[distance_col]] <- dist_raw

  pair$structural_connectivity <- .hee_clip01(
    exp(-pair[[distance_col]] / max(distance_scale, .Machine$double.eps)) *
      pair[[passability_col]]
  )

  if (is.null(traits)) {
    pair$functional_connectivity <- 1
    pair$connectivity <- .hee_clip01(pair$structural_connectivity *
                                       pair[[climate_col]])
    return(pair)
  }

  tr <- as.data.frame(traits)
  .require_cols(tr, c(species_col, dispersal_col), "traits")
  .hee_check_unique_keys(tr, species_col, "traits")
  tr_disp <- suppressWarnings(as.numeric(tr[[dispersal_col]]))
  if (any(!is.finite(tr_disp) | tr_disp < 0, na.rm = TRUE)) {
    stop("traits$", dispersal_col,
         " must contain finite non-negative dispersal distances.",
         call. = FALSE)
  }
  tr[[dispersal_col]] <- tr_disp
  pair$.pair_id <- seq_len(nrow(pair))
  if (!isTRUE(expand_species)) {
    stop("hee_build_connectivity_cube() requires expand_species = TRUE when ",
         "`traits` are supplied and `distances` are region-pair rows.",
         call. = FALSE)
  }
  z <- merge(pair, tr[, unique(c(species_col, dispersal_col)), drop = FALSE],
             all = TRUE, sort = FALSE)
  disp <- suppressWarnings(as.numeric(z[[dispersal_col]]))
  z$functional_connectivity <- ifelse(
    z[[distance_col]] <= 0,
    1,
    exp(-z[[distance_col]] / pmax(disp, .Machine$double.eps))
  )
  z$functional_connectivity <- .hee_clip01(z$functional_connectivity)
  z$connectivity <- .hee_clip01(z$structural_connectivity *
                                  z$functional_connectivity *
                                  z[[climate_col]])
  z[order(z$.pair_id), setdiff(names(z), ".pair_id"), drop = FALSE]
}

#' Build a connectivity cube using workflow-style argument names
#'
#' Convenience wrapper around [hee_build_connectivity_cube()] using the argument
#' names from the Dynamic Earth-Biota Assembly workflow.
#'
#' @param regions Optional region/time coordinate table.
#' @param times Optional ages in Ma.
#' @param cell_table Optional cell/region/time coordinate table used to compute
#'   pairwise distances when `distances` is absent.
#' @param barriers Optional barrier table.
#' @param traits Optional species trait table.
#' @param suitability Optional table with precomputed climate connectivity.
#' @param alpha Distance-decay coefficient; internally converted to
#'   `distance_scale = 1 / alpha`.
#' @param distances Optional distance table or matrix.
#' @param ... Additional arguments passed to [hee_build_connectivity_cube()].
#' @return A connectivity cube table.
#' @export
hee_connectivity_cube <- function(regions = NULL,
                                  times = NULL,
                                  cell_table = NULL,
                                  barriers = NULL,
                                  traits = NULL,
                                  suitability = NULL,
                                  alpha = 0.01,
                                  distances = NULL,
                                  ...) {
  climate <- NULL
  if (!is.null(suitability)) {
    s <- as.data.frame(suitability)
    if (all(c("from_region", "to_region", "time_ma", "climate_connectivity") %in% names(s))) {
      climate <- s
    }
  }
  hee_build_connectivity_cube(
    distances = distances,
    landscape = cell_table %||% regions,
    traits = traits,
    barrier = barriers,
    climate = climate,
    times = times,
    distance_scale = 1 / max(alpha, .Machine$double.eps),
    ...
  )
}

#' Compute isolation history from a connectivity cube
#'
#' Calculates instantaneous isolation (`1 - connectivity`) and the duration of
#' continuous low-connectivity intervals for each region pair, optionally by
#' species.
#'
#' @param connectivity_cube Output from `hee_build_connectivity_cube()` or a
#'   compatible table.
#' @param threshold Connectivity values below this threshold are treated as
#'   isolated.
#' @param connectivity_col Connectivity column.
#' @param species_col,from_col,to_col,time_col Column names.
#'
#' @return Connectivity table with isolation and isolation duration in Ma.
#' @export
#'
#' @examples
#' cc <- data.frame(from_region = "A", to_region = "B",
#'   time_ma = c(10, 5, 0), connectivity = c(0.8, 0.2, 0.1))
#' hee_isolation_history(cc, threshold = 0.5)
hee_isolation_history <- function(connectivity_cube,
                                  threshold = 0.5,
                                  connectivity_col = "connectivity",
                                  species_col = "species",
                                  from_col = "from_region",
                                  to_col = "to_region",
                                  time_col = "time_ma") {
  x <- as.data.frame(connectivity_cube)
  .require_cols(x, c(from_col, to_col, time_col, connectivity_col),
                "connectivity_cube")
  group_cols <- c(intersect(species_col, names(x)), from_col, to_col)
  x$.row_id <- seq_len(nrow(x))
  parts <- split(x, interaction(x[, group_cols, drop = FALSE], drop = TRUE))
  out <- lapply(parts, function(z) {
    z <- z[order(z[[time_col]], decreasing = TRUE), , drop = FALSE]
    z$isolation <- .hee_clip01(1 - z[[connectivity_col]])
    z$isolated <- z[[connectivity_col]] < threshold
    dt <- c(0, abs(diff(z[[time_col]])))
    dur <- numeric(nrow(z))
    for (i in seq_len(nrow(z))) {
      dur[i] <- if (isTRUE(z$isolated[i])) {
        if (i == 1L) 0 else dur[i - 1L] + dt[i]
      } else {
        0
      }
    }
    z$isolation_duration_ma <- dur
    z
  })
  out <- do.call(rbind, out)
  out[order(out$.row_id), setdiff(names(out), ".row_id"), drop = FALSE]
}

#' Compute colonisation source pressure
#'
#' Computes propagule or colonisation pressure from previously occupied source
#' regions through a geological connectivity network.
#'
#' @param previous_state Species-region-time table with an occupancy or
#'   probability column.
#' @param connectivity_cube Region-pair-time connectivity table.
#' @param target Optional target species-region-time table. If omitted, all
#'   region pairs in the connectivity cube are used.
#' @param probability_col Source occupancy probability column.
#' @param connectivity_col Connectivity column.
#' @param distance_col Optional distance column for an extra dispersal kernel.
#' @param dispersal_scale Numeric scalar or named vector by species.
#' @param kernel Distance kernel.
#' @param source_threshold Ignore source probabilities below this value.
#' @param include_self Logical; allow `from_region == to_region` links to
#'   contribute to source pressure. The default `FALSE` treats local persistence
#'   as the survival term in the occupancy recursion rather than as external
#'   colonisation/rescue.
#' @param species_col,region_col,from_col,to_col,time_col Column names.
#'
#' @return Target table with `source_pressure`.
#' @export
#'
#' @examples
#' prev <- data.frame(species = "sp1", region = "A", time_ma = 0,
#'   occupancy_probability = 0.8)
#' cc <- data.frame(from_region = "A", to_region = "B", time_ma = 0,
#'   connectivity = 0.5, paleodistance = 2)
#' hee_colonisation_pressure(prev, cc, dispersal_scale = c(sp1 = 5))
hee_colonisation_pressure <- function(previous_state,
                                      connectivity_cube,
                                      target = NULL,
                                      probability_col = "occupancy_probability",
                                      connectivity_col = "connectivity",
                                      distance_col = "paleodistance",
                                      dispersal_scale = 1,
                                      kernel = c("exponential", "gaussian"),
                                      source_threshold = 0,
                                      include_self = FALSE,
                                      species_col = "species",
                                      region_col = "region",
                                      from_col = "from_region",
                                      to_col = "to_region",
                                      time_col = "time_ma") {
  kernel <- match.arg(kernel)
  prev <- as.data.frame(previous_state)
  .require_cols(prev, c(species_col, region_col, time_col, probability_col),
                "previous_state")
  .hee_validate_time_values(prev[[time_col]], paste0("previous_state$", time_col))
  .hee_check_unique_keys(prev, c(species_col, region_col, time_col),
                         "previous_state")
  conn <- as.data.frame(connectivity_cube)
  .require_cols(conn, c(from_col, to_col, time_col, connectivity_col),
                "connectivity_cube")
  .hee_validate_time_values(conn[[time_col]], paste0("connectivity_cube$", time_col))
  conn_key <- c(from_col, to_col, time_col,
                if (species_col %in% names(conn)) species_col else character())
  .hee_check_unique_keys(conn, conn_key, "connectivity_cube")
  if (is.null(target)) {
    target <- unique(conn[, c(to_col, time_col), drop = FALSE])
    names(target)[names(target) == to_col] <- region_col
    target <- merge(unique(prev[, species_col, drop = FALSE]), target,
                    all = TRUE, sort = FALSE)
  } else {
    target <- as.data.frame(target)
    .require_cols(target, c(species_col, region_col, time_col), "target")
  }
  .hee_validate_time_values(target[[time_col]], paste0("target$", time_col))
  .hee_check_unique_keys(target, c(species_col, region_col, time_col), "target")

  src_prob <- suppressWarnings(as.numeric(prev[[probability_col]]))
  src <- prev[!is.na(src_prob) & src_prob >= source_threshold, , drop = FALSE]
  if (nrow(src) == 0L) {
    target$source_pressure <- 0
    target$n_source_regions <- 0L
    target$max_source_probability <- 0
    return(target)
  }
  src <- src[, c(species_col, region_col, time_col, probability_col), drop = FALSE]
  names(src)[names(src) == region_col] <- from_col
  names(src)[names(src) == probability_col] <- ".source_probability"

  conn_has_species <- species_col %in% names(conn)
  by_cols <- c(from_col, time_col, if (conn_has_species) species_col else character())
  z <- merge(src, conn, by = by_cols, all.x = FALSE, sort = FALSE)
  if (nrow(z) == 0L) {
    target$source_pressure <- 0
    target$n_source_regions <- 0L
    target$max_source_probability <- 0
    return(target)
  }
  if (!conn_has_species) {
    z[[species_col]] <- z[[species_col]]
  }
  z$target_region <- z[[to_col]]
  if (!isTRUE(include_self)) {
    z <- z[as.character(z[[from_col]]) != as.character(z$target_region),
           , drop = FALSE]
    if (nrow(z) == 0L) {
      target$source_pressure <- 0
      target$n_source_regions <- 0L
      target$max_source_probability <- 0
      return(target)
    }
  }
  z$kernel_weight <- 1
  if (distance_col %in% names(z)) {
    sc <- .hee_scale_for(z[[species_col]], dispersal_scale)
    d <- suppressWarnings(as.numeric(z[[distance_col]]))
    if (any(is.finite(d) & d < 0, na.rm = TRUE)) {
      stop("`", distance_col, "` must be non-negative when supplied.",
           call. = FALSE)
    }
    d[is.na(d)] <- Inf
    z$kernel_weight <- if (kernel == "gaussian") {
      exp(-(d^2) / (2 * pmax(sc, .Machine$double.eps)^2))
    } else {
      exp(-d / pmax(sc, .Machine$double.eps))
    }
  }
  z$contribution <- .hee_clip01(z$.source_probability) *
    .hee_clip01(z[[connectivity_col]]) * .hee_clip01(z$kernel_weight)
  z <- z[is.finite(z$contribution) & z$contribution > 0, , drop = FALSE]
  if (nrow(z) == 0L) {
    target$source_pressure <- 0
    target$n_source_regions <- 0L
    target$max_source_probability <- 0
    return(target)
  }

  agg <- stats::aggregate(z$contribution,
                          z[, c(species_col, "target_region", time_col), drop = FALSE],
                          function(v) 1 - prod(1 - .hee_clip01(v), na.rm = TRUE))
  names(agg)[ncol(agg)] <- "source_pressure"
  nsrc <- stats::aggregate(z$contribution,
                           z[, c(species_col, "target_region", time_col), drop = FALSE],
                           length)
  names(nsrc)[ncol(nsrc)] <- "n_source_regions"
  maxp <- stats::aggregate(z$.source_probability,
                           z[, c(species_col, "target_region", time_col), drop = FALSE],
                           max, na.rm = TRUE)
  names(maxp)[ncol(maxp)] <- "max_source_probability"
  agg <- merge(agg, nsrc, by = c(species_col, "target_region", time_col),
               all.x = TRUE, sort = FALSE)
  agg <- merge(agg, maxp, by = c(species_col, "target_region", time_col),
               all.x = TRUE, sort = FALSE)
  names(agg)[names(agg) == "target_region"] <- region_col
  out <- merge(target, agg, by = c(species_col, region_col, time_col),
               all.x = TRUE, sort = FALSE)
  out$source_pressure[is.na(out$source_pressure)] <- 0
  out$n_source_regions[is.na(out$n_source_regions)] <- 0L
  out$max_source_probability[is.na(out$max_source_probability)] <- 0
  out$source_pressure <- .hee_clip01(out$source_pressure)
  out
}

#' Compute rescue effect from neighbouring occupied sources
#'
#' The rescue effect uses the same source-pressure machinery as colonisation,
#' but it is interpreted as a reduction in local extinction risk for already
#' occupied populations. Keeping it as a separate function prevents conflating
#' colonisation of empty sites with demographic rescue of occupied sites.
#'
#' @inheritParams hee_colonisation_pressure
#' @return Target table with `rescue_effect`.
#' @export
#'
#' @examples
#' prev <- data.frame(species = "sp1", region = "A", time_ma = 0,
#'   occupancy_probability = 0.8)
#' cc <- data.frame(from_region = "A", to_region = "B", time_ma = 0,
#'   connectivity = 0.5)
#' hee_rescue_effect(prev, cc)
hee_rescue_effect <- function(previous_state,
                              connectivity_cube,
                              target = NULL,
                              probability_col = "occupancy_probability",
                              connectivity_col = "connectivity",
                              distance_col = "paleodistance",
                              dispersal_scale = 1,
                              kernel = c("exponential", "gaussian"),
                              source_threshold = 0,
                              include_self = FALSE,
                              species_col = "species",
                              region_col = "region",
                              from_col = "from_region",
                              to_col = "to_region",
                              time_col = "time_ma") {
  out <- hee_colonisation_pressure(
    previous_state = previous_state,
    connectivity_cube = connectivity_cube,
    target = target,
    probability_col = probability_col,
    connectivity_col = connectivity_col,
    distance_col = distance_col,
    dispersal_scale = dispersal_scale,
    kernel = match.arg(kernel),
    source_threshold = source_threshold,
    include_self = include_self,
    species_col = species_col,
    region_col = region_col,
    from_col = from_col,
    to_col = to_col,
    time_col = time_col
  )
  names(out)[names(out) == "source_pressure"] <- "rescue_effect"
  out
}

#' Convert source pressure into colonisation probability
#'
#' Computes the interval colonisation probability
#' `gamma = 1 - exp(-lambda_col * delta_t)`, where the colonisation rate is a
#' log-linear function of HMSC suitability, source pressure, historical
#' accessibility, and ecological opportunity. Coefficients are scenario
#' parameters unless calibrated externally. If `source_pressure` is omitted,
#' there is no propagule evidence and colonisation probability is zero. If a
#' source-pressure table is supplied, unmatched rows or explicit zero values
#' also mean no propagule evidence and force colonisation probability to zero.
#' If `accessibility` is omitted, accessibility is neutral (`1`). If an
#' accessibility table is supplied, unmatched rows or `NA` values are treated as
#' not historically accessible (`0`).
#'
#' @param suitability Table or numeric vector containing suitability values.
#' @param source_pressure Optional table/vector with source pressure.
#' @param accessibility Optional table/vector with historical accessibility.
#' @param ecological_opportunity Optional table/vector with ecological
#'   opportunity.
#' @param delta_t Time interval in Myr. May be a scalar or a column name in
#'   `suitability`.
#' @param coefficients Named coefficients `intercept`, `suitability`,
#'   `source_pressure`, `accessibility`, and `ecological_opportunity`.
#' @param keys Keys used to join optional tables.
#' @param suitability_col,source_pressure_col,accessibility_col Column names.
#' @param opportunity_col Output/input ecological-opportunity column.
#' @param eps Small value used before logs/logits.
#' @return A table with `colonisation_rate_lambda` and
#'   `colonisation_probability`.
#' @export
#'
#' @examples
#' x <- data.frame(species = "sp1", cell_id = "c1", time_ma = 0,
#'   suitability = 0.7)
#' hee_colonisation_probability(x, source_pressure = data.frame(
#'   species = "sp1", cell_id = "c1", time_ma = 0, source_pressure = 0.4))
hee_colonisation_probability <- function(suitability,
                                         source_pressure = NULL,
                                         accessibility = NULL,
                                         ecological_opportunity = NULL,
                                         delta_t = 1,
                                         coefficients = c(
                                           intercept = -3,
                                           suitability = 1,
                                           source_pressure = 1,
                                           accessibility = 1,
                                           ecological_opportunity = 0.5
                                         ),
                                         keys = c("species", "cell_id", "region", "time_ma"),
                                         suitability_col = "suitability",
                                         source_pressure_col = "source_pressure",
                                         accessibility_col = "accessibility",
                                         opportunity_col = "ecological_opportunity",
                                         eps = 1e-9) {
  x <- if (is.data.frame(suitability)) as.data.frame(suitability) else
    data.frame(suitability = suitability)
  if (!suitability_col %in% names(x)) {
    stop("suitability must contain column `", suitability_col, "`.", call. = FALSE)
  }
  x$.row_id <- seq_len(nrow(x))
  S <- .hee_clip01(x[[suitability_col]])
  S[is.na(S)] <- 0
  source_default <- 0
  M <- .hee_metric_value(x, source_pressure, source_pressure_col,
                         default = source_default, keys = keys,
                         aliases = c("colonisation_pressure",
                                     "colonization_pressure"),
                         source_name = "source_pressure")
  accessibility_default <- if (is.null(accessibility)) 1 else 0
  A <- .hee_metric_value(x, accessibility, accessibility_col,
                         default = accessibility_default,
                         keys = keys, aliases = c("A_hist", "A_BGB", "weight"),
                         source_name = "accessibility")
  O <- .hee_metric_value(x, ecological_opportunity, opportunity_col, default = 0,
                         keys = keys,
                         aliases = c("opportunity", "opportunity_index",
                                     "ecological_opportunity_index"),
                         source_name = "ecological_opportunity")
  M <- .hee_clip01(M)
  A <- .hee_clip01(A)
  O <- .hee_clip01(O)
  M[is.na(M)] <- source_default
  A[is.na(A)] <- accessibility_default
  O[is.na(O)] <- 0
  dt <- if (is.character(delta_t) && length(delta_t) == 1L && delta_t %in% names(x)) {
    suppressWarnings(as.numeric(x[[delta_t]]))
  } else {
    rep(as.numeric(delta_t)[1], nrow(x))
  }
  dt <- .hee_validate_delta_t(dt)
  cf <- .hee_named_coefficients(coefficients, c("intercept", "suitability",
                                                "source_pressure",
                                                "accessibility",
                                                "ecological_opportunity"))
  eta <- cf["intercept"] +
    cf["suitability"] * .hee_logit_bounded(S, eps) +
    cf["source_pressure"] * log(pmax(M, eps)) +
    cf["accessibility"] * log(pmax(A, eps)) +
    cf["ecological_opportunity"] * .hee_clip01(O)
  lambda <- exp(pmin(pmax(as.numeric(eta), -700), 700))
  lambda[M <= eps | A <= eps | S <= eps] <- 0
  x$source_pressure <- .hee_clip01(M)
  x$accessibility <- .hee_clip01(A)
  x$ecological_opportunity <- .hee_clip01(O)
  x$delta_t_myr <- dt
  x$colonisation_linear_predictor <- as.numeric(eta)
  x$colonisation_rate_lambda <- lambda
  x$colonisation_probability <- .hee_clip01(1 - exp(-lambda * dt))
  x[order(x$.row_id), setdiff(names(x), ".row_id"), drop = FALSE]
}

#' Compute local extinction risk
#'
#' Calculates local extinction probability from suitability decline, land or
#' habitat loss, isolation, small area, rescue effect, disturbance, and optional
#' vulnerability. If an explicit area column is supplied and `area <= 0`, the
#' local arena is treated as geographically absent and extinction is forced.
#' If a `geographic_existence` column is supplied, explicit `NA` values are
#' treated conservatively as absence (`0`); omitting the whole geography column
#' remains neutral.
#' The default coefficients are conservative diagnostics, not fitted causal
#' parameters.
#'
#' @param x Data.frame with species-cell/region-time rows.
#' @param coefficients Named numeric coefficients. Recognised names are
#'   `intercept`, `unsuitability`, `land_loss`, `isolation`, `area_inverse`,
#'   `rescue`, `disturbance`, and `vulnerability`.
#' @param suitability_col,land_loss_col,isolation_col,area_col,rescue_col
#'   Column names.
#' @param geography_col,disturbance_col,vulnerability_col Optional column names.
#'
#' @return `x` with linear predictor and `extinction_probability`.
#' @export
#'
#' @examples
#' x <- data.frame(suitability = c(0.8, 0.1), land_loss = c(0, 1),
#'   isolation = c(0.2, 0.9), area_km2 = c(100, 5), rescue_effect = c(0.5, 0))
#' hee_extinction_risk(x)
hee_extinction_risk <- function(x,
                                coefficients = c(
                                  intercept = -2,
                                  unsuitability = 2,
                                  land_loss = 3,
                                  isolation = 1,
                                  area_inverse = 0.5,
                                  rescue = -2,
                                  disturbance = 1,
                                  vulnerability = 1
                                ),
                                suitability_col = "suitability",
                                land_loss_col = "land_loss",
                                isolation_col = "isolation",
                                area_col = "area_km2",
                                rescue_col = "rescue_effect",
                                geography_col = "geographic_existence",
                                disturbance_col = "disturbance",
                                vulnerability_col = "vulnerability") {
  d <- as.data.frame(x)
  val <- function(col, default) {
    if (is.null(col) || !col %in% names(d)) return(rep(default, nrow(d)))
    out <- suppressWarnings(as.numeric(d[[col]]))
    out[is.na(out)] <- default
    out
  }
  geo_val <- function(col) {
    if (is.null(col) || !col %in% names(d)) return(rep(1, nrow(d)))
    out <- suppressWarnings(as.numeric(d[[col]]))
    out[is.na(out)] <- 0
    .hee_clip01(out)
  }
  suit <- .hee_clip01(val(suitability_col, 1))
  land_loss <- .hee_clip01(val(land_loss_col, 0))
  isolation <- .hee_clip01(val(isolation_col, 0))
  area <- val(area_col, NA_real_)
  area_inv <- .hee_area_inverse_index(area)
  rescue <- .hee_clip01(val(rescue_col, 0))
  disturbance <- .hee_clip01(val(disturbance_col, 0))
  vulnerability <- .hee_clip01(val(vulnerability_col, 0))
  geography <- geo_val(geography_col)

  cf <- .hee_named_coefficients(coefficients, c("intercept", "unsuitability",
                                                "land_loss", "isolation",
                                                "area_inverse", "rescue",
                                                "disturbance",
                                                "vulnerability"))
  eta <- cf["intercept"] +
    cf["unsuitability"] * (1 - suit) +
    cf["land_loss"] * land_loss +
    cf["isolation"] * isolation +
    cf["area_inverse"] * area_inv +
    cf["rescue"] * rescue +
    cf["disturbance"] * disturbance +
    cf["vulnerability"] * vulnerability
  d$extinction_linear_predictor <- as.numeric(eta)
  d$extinction_probability <- stats::plogis(d$extinction_linear_predictor)
  d$forced_geographic_extinction <- land_loss >= 1 | geography <= 0 |
    (!is.na(area) & area <= 0)
  d$extinction_probability[d$forced_geographic_extinction] <- 1
  d
}

#' Convert extinction drivers into interval extinction probability
#'
#' Computes `epsilon = 1 - exp(-lambda_ext * delta_t)` from unsuitability,
#' land/habitat loss, isolation, inverse area, rescue effect, disturbance, and
#' species vulnerability. Rows with `land_loss >= 1` are treated as forced
#' geographic local extinctions (`epsilon = 1`). Rows with
#' `geographic_existence <= 0`, explicit `NA` geographic existence, or an
#' explicitly supplied `area <= 0` are also forced to one, because the local
#' arena does not exist. Coefficients are scenario parameters unless calibrated
#' externally.
#'
#' @param x Species-cell/region-time table.
#' @param delta_t Time interval in Myr, scalar or column name.
#' @param coefficients Named coefficients `intercept`, `unsuitability`,
#'   `land_loss`, `isolation`, `area_inverse`, `rescue`, `disturbance`, and
#'   `vulnerability`.
#' @inheritParams hee_extinction_risk
#' @return `x` with `extinction_rate_lambda`, `extinction_probability`, and
#'   `forced_geographic_extinction`.
#' @export
#'
#' @examples
#' x <- data.frame(suitability = c(0.8, 0.1), land_loss = c(0, 1),
#'   isolation = c(0.1, 0.8), area_km2 = c(100, 2), rescue_effect = c(0.5, 0))
#' hee_extinction_probability(x, delta_t = 5)
hee_extinction_probability <- function(x,
                                       delta_t = 1,
                                       coefficients = c(
                                         intercept = -3,
                                         unsuitability = 1,
                                         land_loss = 3,
                                         isolation = 1,
                                         area_inverse = 0.5,
                                         rescue = -1,
                                         disturbance = 1,
                                         vulnerability = 1
                                       ),
                                       suitability_col = "suitability",
                                       land_loss_col = "land_loss",
                                       isolation_col = "isolation",
                                       area_col = "area_km2",
                                       rescue_col = "rescue_effect",
                                       geography_col = "geographic_existence",
                                       disturbance_col = "disturbance",
                                       vulnerability_col = "vulnerability") {
  d <- as.data.frame(x)
  val <- function(col, default) {
    if (is.null(col) || !col %in% names(d)) return(rep(default, nrow(d)))
    out <- suppressWarnings(as.numeric(d[[col]]))
    out[is.na(out)] <- default
    out
  }
  geo_val <- function(col) {
    if (is.null(col) || !col %in% names(d)) return(rep(1, nrow(d)))
    out <- suppressWarnings(as.numeric(d[[col]]))
    out[is.na(out)] <- 0
    .hee_clip01(out)
  }
  suit <- .hee_clip01(val(suitability_col, 1))
  land_loss <- .hee_clip01(val(land_loss_col, 0))
  isolation <- .hee_clip01(val(isolation_col, 0))
  area <- val(area_col, NA_real_)
  area_inv <- .hee_area_inverse_index(area)
  rescue <- .hee_clip01(val(rescue_col, 0))
  disturbance <- .hee_clip01(val(disturbance_col, 0))
  vulnerability <- .hee_clip01(val(vulnerability_col, 0))
  geography <- geo_val(geography_col)
  dt <- if (is.character(delta_t) && length(delta_t) == 1L && delta_t %in% names(d)) {
    suppressWarnings(as.numeric(d[[delta_t]]))
  } else {
    rep(as.numeric(delta_t)[1], nrow(d))
  }
  dt <- .hee_validate_delta_t(dt)

  cf <- .hee_named_coefficients(coefficients, c("intercept", "unsuitability",
                                                "land_loss", "isolation",
                                                "area_inverse", "rescue",
                                                "disturbance",
                                                "vulnerability"))
  eta <- cf["intercept"] +
    cf["unsuitability"] * (1 - suit) +
    cf["land_loss"] * land_loss +
    cf["isolation"] * isolation +
    cf["area_inverse"] * area_inv +
    cf["rescue"] * rescue +
    cf["disturbance"] * disturbance +
    cf["vulnerability"] * vulnerability
  d$delta_t_myr <- dt
  d$extinction_rate_lambda <- exp(pmin(pmax(as.numeric(eta), -700), 700))
  d$forced_geographic_extinction <- land_loss >= 1 | geography <= 0 |
    (!is.na(area) & area <= 0)
  d$extinction_probability <- .hee_clip01(1 - exp(-d$extinction_rate_lambda * dt))
  d$extinction_probability[d$forced_geographic_extinction] <- 1
  d
}

#' Compute speciation opportunity indices
#'
#' Constructs diagnostic opportunity indices for allopatric, founder-event,
#' in-situ, and radiation mechanisms. Missing process evidence is conservative:
#' unknown ecological opportunity, connectivity, land age, or persistence does
#' not create a positive opportunity by default. This function does not estimate
#' true speciation rates unless the user supplies calibrated coefficients.
#'
#' @param x Region/species/time table.
#' @param coefficients Coefficients for a log-link opportunity score.
#' @param isolation_col,niche_divergence_col,area_col,heterogeneity_col Column
#'   names.
#' @param persistence_col,source_pressure_col,connectivity_col Column names.
#' @param opportunity_col Optional ecological opportunity column.
#' @param land_age_col Optional land age column.
#' @param land_age_tau Decay scale in Myr for young-arena radiation
#'   opportunity.
#'
#' @return `x` with opportunity indices, `speciation_opportunity_lambda`, and
#'   the backwards-compatible alias `speciation_rate_lambda`.
#' @export
#'
#' @examples
#' x <- data.frame(isolation_duration_ma = 5, niche_divergence = 0.7,
#'   area_km2 = 100, habitat_heterogeneity = 0.6, population_persistence = 0.8,
#'   source_pressure = 0.4, connectivity = 0.3, ecological_opportunity = 0.9)
#' hee_speciation_opportunity(x)
hee_speciation_opportunity <- function(x,
                                       coefficients = c(
                                         intercept = -3,
                                         isolation_duration = 0.02,
                                         niche_divergence = 1,
                                         area = 0.2,
                                         habitat_heterogeneity = 1,
                                         population_persistence = 1,
                                         ecological_opportunity = 1
                                       ),
                                       isolation_col = "isolation_duration_ma",
                                       niche_divergence_col = "niche_divergence",
                                       area_col = "area_km2",
                                       heterogeneity_col = "habitat_heterogeneity",
                                       persistence_col = "population_persistence",
                                       source_pressure_col = "source_pressure",
                                       connectivity_col = "connectivity",
                                       opportunity_col = "ecological_opportunity",
                                       land_age_col = "land_age_ma",
                                       land_age_tau = 10) {
  d <- as.data.frame(x)
  val <- function(col, default) {
    if (is.null(col) || !col %in% names(d)) return(rep(default, nrow(d)))
    out <- suppressWarnings(as.numeric(d[[col]]))
    out[is.na(out)] <- default
    out
  }
  iso <- .hee_positive_index(val(isolation_col, 0))
  niche <- .hee_clip01(val(niche_divergence_col, 0))
  area <- .hee_positive_index(val(area_col, 0))
  het <- .hee_clip01(val(heterogeneity_col, 0))
  persist <- .hee_clip01(val(persistence_col, 0))
  src <- .hee_clip01(val(source_pressure_col, 0))
  conn <- .hee_clip01(val(connectivity_col, 1))
  opp <- .hee_clip01(val(opportunity_col, 0))
  raw_land_age <- val(land_age_col, NA_real_)
  young_land <- if (all(is.na(raw_land_age))) {
    rep(0, nrow(d))
  } else {
    exp(-pmax(raw_land_age, 0) / max(land_age_tau, .Machine$double.eps))
  }
  young_land[is.na(young_land)] <- 0

  d$vicariance_opportunity <- .hee_clip01(iso * niche * persist)
  d$founder_event_opportunity <- .hee_clip01(src * (1 - conn) * opp)
  d$in_situ_speciation_index <- .hee_clip01(opp * het * persist)
  d$young_land_index <- .hee_clip01(young_land)
  d$radiation_opportunity_index <- .hee_clip01(d$young_land_index * opp * het * src)
  d$speciation_opportunity <- rowMeans(d[, c("vicariance_opportunity",
                                             "founder_event_opportunity",
                                             "in_situ_speciation_index",
                                             "radiation_opportunity_index")],
                                       na.rm = TRUE)

  cf <- .hee_named_coefficients(coefficients, c("intercept",
                                                "isolation_duration",
                                                "niche_divergence", "area",
                                                "habitat_heterogeneity",
                                                "population_persistence",
                                                "ecological_opportunity"))
  eta <- cf["intercept"] +
    cf["isolation_duration"] * val(isolation_col, 0) +
    cf["niche_divergence"] * niche +
    cf["area"] * area +
    cf["habitat_heterogeneity"] * het +
    cf["population_persistence"] * persist +
    cf["ecological_opportunity"] * opp
  d$speciation_linear_predictor <- as.numeric(eta)
  d$speciation_opportunity_lambda <- exp(pmin(pmax(d$speciation_linear_predictor, -700), 700)) *
    d$speciation_opportunity
  d$speciation_rate_lambda <- d$speciation_opportunity_lambda
  d
}

#' Allopatric speciation opportunity
#'
#' Diagnostic opportunity for vicariance/allopatric divergence driven by range
#' fragmentation, persistent isolation, and population persistence.
#'
#' @param x Table with input columns.
#' @param range_fragmentation_col,isolation_col,persistence_col Column names.
#' @return `x` with `allopatric_opportunity`.
#' @export
hee_allopatric_speciation_opportunity <- function(x,
                                                  range_fragmentation_col = "range_fragmentation",
                                                  isolation_col = "isolation_duration_ma",
                                                  persistence_col = "population_persistence") {
  d <- as.data.frame(x)
  frag <- .hee_clip01(.hee_col_or_default(d, range_fragmentation_col, 0))
  iso <- .hee_positive_index(.hee_col_or_default(d, isolation_col, 0))
  persist <- .hee_clip01(.hee_col_or_default(d, persistence_col, 0))
  d$allopatric_opportunity <- .hee_clip01(frag * iso * persist)
  d
}

#' Founder-event speciation opportunity
#'
#' Diagnostic opportunity for rare dispersal followed by isolation in a newly
#' available ecological arena.
#'
#' @param x Table with input columns.
#' @param rare_dispersal_col,isolation_col,opportunity_col Column names.
#' @return `x` with `founder_event_opportunity`.
#' @export
hee_founder_speciation_opportunity <- function(x,
                                               rare_dispersal_col = "rare_dispersal",
                                               isolation_col = "isolation_after_arrival",
                                               opportunity_col = "ecological_opportunity") {
  d <- as.data.frame(x)
  rare <- .hee_clip01(.hee_col_or_default(d, rare_dispersal_col, 0))
  iso <- .hee_clip01(.hee_col_or_default(d, isolation_col, 0))
  opp <- .hee_clip01(.hee_col_or_default(d, opportunity_col, 0))
  d$founder_event_opportunity <- .hee_clip01(rare * iso * opp)
  d
}

#' In-situ ecological speciation opportunity
#'
#' Diagnostic opportunity for within-region ecological divergence. High
#' ecological opportunity, heterogeneity, and area increase the index, while
#' high connectivity can reduce isolation-mediated divergence.
#'
#' @param x Table with input columns.
#' @param opportunity_col,heterogeneity_col,area_col,connectivity_col Column
#'   names.
#' @return `x` with `in_situ_speciation_index`.
#' @export
hee_insitu_speciation_opportunity <- function(x,
                                              opportunity_col = "ecological_opportunity",
                                              heterogeneity_col = "habitat_heterogeneity",
                                              area_col = "area_km2",
                                              connectivity_col = "connectivity") {
  d <- as.data.frame(x)
  opp <- .hee_clip01(.hee_col_or_default(d, opportunity_col, 0))
  het <- .hee_clip01(.hee_col_or_default(d, heterogeneity_col, 0))
  area <- .hee_positive_index(.hee_col_or_default(d, area_col, 0))
  conn <- .hee_clip01(.hee_col_or_default(d, connectivity_col, 1))
  d$in_situ_speciation_index <- .hee_clip01(opp * het * (0.5 + 0.5 * area) * (1 - conn))
  d
}

#' Post-arena radiation opportunity
#'
#' Diagnostic opportunity for adaptive radiation after new ecological arenas
#' form. This index intentionally remains an opportunity score; it is not a
#' speciation-rate estimator without external calibration.
#'
#' @param x Table with input columns.
#' @param opportunity_col,colonisation_col,richness_col,persistence_col Column
#'   names.
#' @return `x` with `radiation_opportunity`.
#' @export
hee_radiation_opportunity <- function(x,
                                      opportunity_col = "ecological_opportunity",
                                      colonisation_col = "colonisation_input",
                                      richness_col = "richness",
                                      persistence_col = "persistence") {
  d <- as.data.frame(x)
  opp <- .hee_clip01(.hee_col_or_default(d, opportunity_col, 0))
  col <- .hee_clip01(.hee_col_or_default(d, colonisation_col, 0))
  rich <- .hee_positive_index(.hee_col_or_default(d, richness_col, 0))
  low_comp <- 1 - rich
  persist <- .hee_clip01(.hee_col_or_default(d, persistence_col, 1))
  d$radiation_opportunity <- .hee_clip01(opp * col * low_comp * persist)
  d
}

#' Update region-level species pool through time
#'
#' Updates regional species-pool size from event counts or probabilities:
#' `pool[t+1] = pool[t] + speciation + immigration - emigration - extinction`.
#'
#' @param initial_pool Data.frame with region and starting pool size. If a time
#'   column is supplied, the oldest row per region is used as the initial state.
#' @param events Optional event table with `time_ma`, event type, region columns,
#'   and optional probability/weight.
#' @param times Time axis in Ma. Larger values are older.
#' @param pool_col Initial pool size column.
#' @param event_type_col Event type column.
#' @param probability_col Event weight column.
#' @param region_col,from_col,to_col,time_col Column names.
#'
#' @return Region-time species-pool table.
#' @export
#'
#' @examples
#' init <- data.frame(region = "A", species_pool_size = 2)
#' ev <- data.frame(time_ma = c(5, 0), event_type = c("speciation", "extinction"),
#'   region = "A", probability = c(1, 0.5))
#' hee_update_species_pool(init, ev, times = c(10, 5, 0))
hee_update_species_pool <- function(initial_pool,
                                    events = NULL,
                                    times,
                                    pool_col = "species_pool_size",
                                    event_type_col = "event_type",
                                    probability_col = "probability",
                                    region_col = "region",
                                    from_col = "from_region",
                                    to_col = "to_region",
                                    time_col = "time_ma") {
  pool <- as.data.frame(initial_pool)
  .require_cols(pool, c(region_col, pool_col), "initial_pool")
  times <- sort(as.numeric(times), decreasing = TRUE)
  if (length(times) < 1L || anyNA(times)) {
    stop("times must contain at least one numeric age.", call. = FALSE)
  }
  regions <- sort(unique(as.character(pool[[region_col]])))
  if (time_col %in% names(pool)) {
    .hee_validate_time_values(pool[[time_col]], paste0("initial_pool$", time_col))
    pool <- pool[order(pool[[region_col]], pool[[time_col]], decreasing = TRUE),
                 , drop = FALSE]
    init <- do.call(rbind, lapply(split(pool, pool[[region_col]], drop = TRUE),
                                  function(z) {
      z <- z[order(z[[time_col]], decreasing = TRUE), , drop = FALSE]
      z <- z[!is.na(z[[pool_col]]), , drop = FALSE]
      if (nrow(z) == 0L) return(NULL)
      z[1, c(region_col, pool_col), drop = FALSE]
    }))
    if (is.null(init) || nrow(init) == 0L) {
      stop("initial_pool has no non-missing pool values.", call. = FALSE)
    }
    names(init)[names(init) == pool_col] <- ".pool"
  } else {
    init <- stats::aggregate(pool[[pool_col]], pool[, region_col, drop = FALSE],
                             function(v) v[which(!is.na(v))[1]])
    names(init)[2] <- ".pool"
  }
  if (is.null(events)) {
    events <- data.frame()
  } else {
    events <- as.data.frame(events)
    .require_cols(events, c(time_col, event_type_col), "events")
    if (event_type_col != "event_type") events$event_type <- events[[event_type_col]]
    if (!probability_col %in% names(events)) events[[probability_col]] <- 1
  }

  out <- list()
  n <- 0L
  for (reg in regions) {
    current <- as.numeric(init$.pool[match(reg, init[[region_col]])])
    for (tt in times) {
      ev <- events[events[[time_col]] == tt, , drop = FALSE]
      gain_spec <- .hee_event_sum(ev, reg, c("speciation",
                                             "within_area_speciation"),
                                  probability_col, region_col, from_col, to_col,
                                  direction = "in")
      immigration <- .hee_event_sum(ev, reg, c("colonisation", "immigration",
                                               "range_expansion", "founder_event"),
                                    probability_col, region_col, from_col, to_col,
                                    direction = "in")
      emigration <- .hee_event_sum(ev, reg, c("emigration"),
                                   probability_col, region_col, from_col, to_col,
                                   direction = "out")
      extinction <- .hee_event_sum(ev, reg, c("local_extinction",
                                              "regional_extinction",
                                              "lineage_extinction",
                                              "range_contraction",
                                              "extinction"),
                                   probability_col, region_col, from_col, to_col,
                                   direction = "in")
      current <- max(0, current + gain_spec + immigration - emigration - extinction)
      n <- n + 1L
      out[[n]] <- data.frame(
        region = reg,
        time_ma = tt,
        species_pool_size = current,
        origination = gain_spec,
        immigration = immigration,
        emigration = emigration,
        regional_extinction = extinction,
        source_sink_balance = immigration - emigration,
        stringsAsFactors = FALSE
      )
    }
  }
  res <- do.call(rbind, out)
  names(res)[names(res) == "region"] <- region_col
  names(res)[names(res) == "time_ma"] <- time_col
  res
}

#' Update local occupancy probability for one dynamic interval
#'
#' Applies the core dynamic-state equation
#' `p_next = G * E * (p_prev * (1 - epsilon) + (1 - p_prev) * gamma)`.
#' Inputs can be numeric vectors or a table. This function is intentionally
#' small so the transition used inside [hee_dynamic_assembly()] can be tested
#' independently.
#'
#' @param p_prev Numeric previous occupancy probability, or a data.frame.
#' @param gamma Colonisation probability, or a table joined to `p_prev`.
#' @param epsilon Local extinction probability, or a table joined to `p_prev`.
#' @param phylo_mask Lineage-existence mask. Omitted masks are neutral (`1`);
#'   supplied table/vector masks with unmatched or `NA` values are conservative
#'   (`0`).
#' @param geo_mask Geographic/habitat-existence mask. Omitted masks are neutral
#'   (`1`); supplied table/vector masks with unmatched or `NA` values are
#'   conservative (`0`).
#' @param keys Keys used when joining table inputs.
#' @param previous_col,gamma_col,epsilon_col,phylo_col,geo_col Column names.
#' @return Numeric vector or table with `occupancy_probability`.
#' @export
#'
#' @examples
#' hee_update_occupancy(p_prev = 0.8, gamma = 0.2, epsilon = 0.5,
#'   phylo_mask = 1, geo_mask = 1)
hee_update_occupancy <- function(p_prev,
                                 gamma,
                                 epsilon,
                                 phylo_mask = 1,
                                 geo_mask = 1,
                                 keys = c("species", "cell_id", "region", "time_ma"),
                                 previous_col = "previous_probability",
                                 gamma_col = "colonisation_probability",
                                 epsilon_col = "extinction_probability",
                                 phylo_col = "lineage_exists",
                                 geo_col = "geographic_existence") {
  if (!is.data.frame(p_prev)) {
    p0 <- .hee_clip01(p_prev)
    g <- .hee_clip01(gamma)
    e <- .hee_clip01(epsilon)
    phy <- .hee_clip01(phylo_mask)
    geo <- .hee_clip01(geo_mask)
    p0[is.na(p0)] <- 0
    g[is.na(g)] <- 0
    e[is.na(e)] <- 0
    phy[is.na(phy)] <- if (missing(phylo_mask)) 1 else 0
    geo[is.na(geo)] <- if (missing(geo_mask)) 1 else 0
    return(.hee_clip01(geo * phy * (p0 * (1 - e) + (1 - p0) * g)))
  }
  d <- as.data.frame(p_prev)
  d$.row_id <- seq_len(nrow(d))
  p0 <- .hee_clip01(.hee_col_or_default(d, previous_col, 0))
  g <- .hee_metric_value(d, gamma, gamma_col, default = 0, keys = keys,
                         aliases = c("gamma", "colonization_probability"),
                         source_name = "colonisation probability")
  e <- .hee_metric_value(d, epsilon, epsilon_col, default = 0, keys = keys,
                         aliases = c("epsilon", "local_extinction_probability"),
                         source_name = "extinction probability")
  phy_default <- if (missing(phylo_mask) || is.null(phylo_mask)) 1 else 0
  geo_default <- if (missing(geo_mask) || is.null(geo_mask)) 1 else 0
  phy <- .hee_metric_value(d, phylo_mask, phylo_col, default = phy_default,
                           keys = keys,
                           aliases = c("E_phylo", "phylo_existence",
                                       "phylo_present"),
                           source_name = "phylo_mask")
  geo <- .hee_metric_value(d, geo_mask, geo_col, default = geo_default,
                           keys = keys,
                           aliases = c("G", "geo_mask", "land",
                                       "land_weight", "land_mask_dem",
                                       "habitat_availability"),
                           source_name = "geo_mask")
  d$occupancy_probability <- hee_update_occupancy(p0, g, e, phy, geo)
  d[order(d$.row_id), setdiff(names(d), ".row_id"), drop = FALSE]
}

#' Dynamic local assembly state transition
#'
#' Applies the dynamic occupancy update:
#' `p_next = G * E * (p_prev * (1 - epsilon) + (1 - p_prev) * gamma)`.
#' `G` is geographic/habitat existence and `E` is lineage existence. Time is in
#' Ma, so larger ages are older. Missing colonisation and extinction tables are
#' neutral (`gamma = 0`, `epsilon = 0`); compute those process tables
#' explicitly when they should affect the state transition.
#'
#' @param previous_state Previous or full state table with occupancy
#'   probabilities. If `NULL`, `initial_probability` is used.
#' @param suitability Target species-cell-time suitability table.
#' @param colonisation Optional table with colonisation probability or source
#'   pressure. If omitted, colonisation probability is zero. If only source
#'   pressure is supplied, it is converted with [hee_colonisation_probability()]
#'   using the internally calculated interval length between adjacent time
#'   slices.
#' @param extinction Optional table with extinction probability. If omitted,
#'   extinction probability is zero.
#' @param geo_mask Optional table with geographic existence or land mask. If
#'   the whole table is omitted, geography is neutral (`1`); if the table is
#'   supplied but a row is unmatched or `NA`, the row is treated as
#'   geographically absent (`0`).
#' @param phylo_mask Optional species-time lineage mask, either as a long table
#'   with one of `E_phylo`, `phylo_present`, `lineage_exists`, or
#'   `phylo_existence`, or as a species-by-time matrix such as the output of
#'   [hee_phylo_time_mask()]. If the whole mask is omitted, lineage existence is
#'   neutral (`1`); if the mask is supplied but a row is unmatched or `NA`, the
#'   lineage is treated as absent (`0`).
#' @param time_direction `"forward"` means older-to-younger updates.
#' @param initial_probability Probability used when no previous state exists.
#' @param species_col,cell_col,time_col Column names.
#' @param suitability_col,probability_col Column names.
#'
#' @return Dynamic state table.
#' @export
#'
#' @examples
#' suit <- data.frame(species = "sp1", cell_id = "c1", time_ma = c(5, 0),
#'   suitability = c(0.8, 0.7))
#' prev <- data.frame(species = "sp1", cell_id = "c1", time_ma = 5,
#'   occupancy_probability = 0.6)
#' hee_dynamic_assembly(prev, suit)
hee_dynamic_assembly <- function(previous_state = NULL,
                                 suitability,
                                 colonisation = NULL,
                                 extinction = NULL,
                                 geo_mask = NULL,
                                 phylo_mask = NULL,
                                 time_direction = c("forward", "backward"),
                                 initial_probability = 0,
                                 species_col = "species",
                                 cell_col = "cell_id",
                                 time_col = "time_ma",
                                 suitability_col = "suitability",
                                 probability_col = "occupancy_probability") {
  time_direction <- match.arg(time_direction)
  target <- as.data.frame(suitability)
  .require_cols(target, c(species_col, cell_col, time_col, suitability_col),
                "suitability")
  .hee_validate_time_values(target[[time_col]], paste0("suitability$", time_col))
  dup_dynamic_time <- duplicated(target[, c(species_col, cell_col, time_col), drop = FALSE])
  if (any(dup_dynamic_time)) {
    bad <- unique(target[dup_dynamic_time, c(species_col, cell_col, time_col),
                         drop = FALSE])
    stop("Duplicate time values are not allowed within a species-cell dynamic ",
         "assembly series. Example duplicate: ",
         paste(paste(names(bad), bad[1, ], sep = "="), collapse = ", "),
         call. = FALSE)
  }
  .hee_check_unique_keys(target, c(species_col, cell_col, time_col),
                         "suitability")
  target$.row_id <- seq_len(nrow(target))
  target$suitability <- .hee_clip01(target[[suitability_col]])
  prev_for_interval <- if (is.null(previous_state)) NULL else
    as.data.frame(previous_state)
  if (!is.null(prev_for_interval)) {
    .require_cols(prev_for_interval, c(species_col, cell_col, time_col,
                                       probability_col), "previous_state")
    .hee_validate_time_values(prev_for_interval[[time_col]],
                              paste0("previous_state$", time_col))
    .hee_check_unique_keys(prev_for_interval,
                           c(species_col, cell_col, time_col),
                           "previous_state")
  }
  target$interval_myr <- .hee_dynamic_intervals(target, prev_for_interval,
                                                species_col, cell_col,
                                                time_col, time_direction)
  target$previous_probability <- NA_real_

  col <- .hee_join_optional(target, colonisation,
                            c(species_col, cell_col, time_col),
                            c("colonisation_probability", "source_pressure"))
  if ("colonisation_probability" %in% names(col)) {
    gamma <- col$colonisation_probability
  } else if ("source_pressure" %in% names(col)) {
    gamma_calc <- hee_colonisation_probability(
      target,
      source_pressure = col$source_pressure,
      delta_t = "interval_myr",
      keys = c(species_col, cell_col, time_col),
      suitability_col = "suitability"
    )
    gamma <- gamma_calc$colonisation_probability
  } else {
    gamma <- 0
  }
  gamma[is.na(gamma)] <- 0
  target$colonisation_probability <- .hee_clip01(gamma)
  if ("source_pressure" %in% names(col)) {
    col$source_pressure[is.na(col$source_pressure)] <- 0
    target$source_pressure <- .hee_clip01(col$source_pressure)
  }

  ext <- .hee_join_optional(target, extinction,
                            c(species_col, cell_col, time_col),
                            "extinction_probability")
  eps <- if ("extinction_probability" %in% names(ext)) {
    ext$extinction_probability
  } else {
    0
  }
  eps[is.na(eps)] <- 0
  target$extinction_probability <- .hee_clip01(eps)

  gm <- .hee_join_optional(target, geo_mask, c(cell_col, time_col),
                           c("geographic_existence", "habitat_availability",
                             "land", "L_land", "G_arena"))
  geo_cols <- intersect(c("geographic_existence", "habitat_availability",
                          "land", "L_land", "G_arena"), names(gm))
  if (!is.null(geo_mask) && length(geo_cols) == 0L) {
    stop("geo_mask must contain one of: geographic_existence, ",
         "habitat_availability, land, L_land, G_arena.", call. = FALSE)
  }
  G <- rep(1, nrow(target))
  for (nm in geo_cols) {
    v <- gm[[nm]]
    v[is.na(v)] <- if (is.null(geo_mask)) 1 else 0
    G <- G * v
  }
  target$geographic_existence <- .hee_clip01(G)

  if (is.null(phylo_mask)) {
    E <- rep(1, nrow(target))
  } else if (is.matrix(phylo_mask)) {
    E <- .hee_phylo_value(target, phylo_mask, species_col, time_col)
    E[is.na(E)] <- 0
  } else {
    pm <- .hee_join_optional(target, phylo_mask, c(species_col, time_col),
                             c("E_phylo", "phylo_present", "lineage_exists",
                               "phylo_existence"))
    phylo_cols <- intersect(c("E_phylo", "phylo_present", "lineage_exists",
                              "phylo_existence"), names(pm))
    if (length(phylo_cols) == 0L) {
      stop("phylo_mask must contain one of: E_phylo, phylo_present, ",
           "lineage_exists, phylo_existence.", call. = FALSE)
    }
    E <- rep(1, nrow(target))
    for (nm in phylo_cols) {
      v <- pm[[nm]]
      v[is.na(v)] <- 0
      E <- E * v
    }
  }
  target$lineage_exists <- .hee_clip01(E)

  p_prev <- .hee_clip01(target$previous_probability)
  gamma <- .hee_clip01(target$colonisation_probability)
  eps <- .hee_clip01(target$extinction_probability)
  G <- .hee_clip01(target$geographic_existence)
  E <- .hee_clip01(target$lineage_exists)

  prev <- prev_for_interval
  target$occupancy_probability <- NA_real_
  group_key <- interaction(target[[species_col]], target[[cell_col]], drop = TRUE)
  rows_by_group <- split(seq_len(nrow(target)), group_key)
  prev_key <- if (is.null(prev)) NULL else interaction(prev[[species_col]], prev[[cell_col]], drop = TRUE)

  for (key in names(rows_by_group)) {
    idx <- rows_by_group[[key]]
    if (any(duplicated(target[[time_col]][idx]))) {
      stop("Duplicate time values are not allowed within a species-cell dynamic ",
           "assembly series. Duplicate time: ",
           paste(unique(target[[time_col]][idx][duplicated(target[[time_col]][idx])]),
                 collapse = ", "),
           call. = FALSE)
    }
    ord <- order(target[[time_col]][idx], decreasing = identical(time_direction, "forward"))
    idx <- idx[ord]
    src <- if (is.null(prev)) data.frame() else prev[prev_key == key, , drop = FALSE]
    current <- as.numeric(initial_probability)[1]
    has_state <- FALSE
    for (jj in seq_along(idx)) {
      i <- idx[jj]
      tt <- target[[time_col]][i]
      use_initial_only <- FALSE
      if (jj == 1L) {
        if (nrow(src) > 0L) {
          exact <- src[src[[time_col]] == tt, , drop = FALSE]
          if (nrow(exact) > 0L) {
            current <- as.numeric(exact[[probability_col]][1])
            use_initial_only <- TRUE
          } else {
            candidates <- if (identical(time_direction, "forward")) {
              src[src[[time_col]] > tt, , drop = FALSE]
            } else {
              src[src[[time_col]] < tt, , drop = FALSE]
            }
            if (nrow(candidates) > 0L) {
              nearest <- which.min(abs(candidates[[time_col]] - tt))
              current <- as.numeric(candidates[[probability_col]][nearest])
              use_initial_only <- FALSE
            } else {
              use_initial_only <- TRUE
            }
          }
        } else {
          use_initial_only <- TRUE
        }
        has_state <- TRUE
      }
      if (!has_state || is.na(current)) current <- as.numeric(initial_probability)[1]
      target$previous_probability[i] <- .hee_clip01(current)
      if (isTRUE(use_initial_only)) {
        next_p <- .hee_clip01(G[i] * E[i] * target$previous_probability[i])
      } else {
        next_p <- hee_update_occupancy(target$previous_probability[i], gamma[i], eps[i],
                                       phylo_mask = E[i], geo_mask = G[i])
      }
      target$occupancy_probability[i] <- next_p
      current <- next_p
      has_state <- TRUE
    }
  }
  target[order(target$.row_id), setdiff(names(target), ".row_id"), drop = FALSE]
}

#' Attribute the limiting mechanism in dynamic deep-time predictions
#'
#' Classifies each species-cell-time prediction into the most likely limiting
#' or dominant mechanism. This is a rule-based diagnostic map, not causal proof.
#'
#' @param dynamic_state Output from `hee_dynamic_assembly()` or a compatible
#'   table.
#' @param thresholds Named thresholds for geography, lineage, suitability,
#'   source pressure, connectivity, extinction, speciation, and no-analog risk.
#' @param columns Named vector/list mapping standard mechanism inputs to column
#'   names.
#'
#' @return `dynamic_state` with `dominant_mechanism` and `mechanism_class`.
#' @export
#'
#' @examples
#' x <- data.frame(geographic_existence = c(0, 1), lineage_exists = 1,
#'   suitability = c(0.8, 0.1), source_pressure = c(0.5, 0.5),
#'   extinction_probability = c(0.1, 0.2))
#' hee_mechanism_attribution(x)
hee_mechanism_attribution <- function(dynamic_state,
                                      thresholds = c(
                                        geography = 0.5,
                                        lineage = 0.5,
                                        suitability = 0.3,
                                        source_pressure = 0.1,
                                        connectivity = 0.2,
                                        extinction = 0.7,
                                        speciation = 0.7,
                                        no_analog = 0.7
                                      ),
                                      columns = c(
                                        geography = "geographic_existence",
                                        lineage = "lineage_exists",
                                        suitability = "suitability",
                                        source_pressure = "source_pressure",
                                        connectivity = "connectivity",
                                        extinction = "extinction_probability",
                                        speciation = "speciation_opportunity",
                                        no_analog = "no_analog_score"
                                      )) {
  x <- as.data.frame(dynamic_state)
  get <- function(key, default) {
    col <- unname(columns[key])
    if (is.na(col) || !col %in% names(x)) return(rep(default, nrow(x)))
    out <- suppressWarnings(as.numeric(x[[col]]))
    out[is.na(out)] <- default
    out
  }
  g <- get("geography", 1)
  e <- get("lineage", 1)
  s <- get("suitability", 1)
  m <- get("source_pressure", 1)
  cnet <- get("connectivity", 1)
  eps <- get("extinction", 0)
  lambda <- get("speciation", 0)
  q <- get("no_analog", 0)

  reason <- rep("suitable_persistent_or_colonisable", nrow(x))
  reason[q >= thresholds["no_analog"]] <- "low_reliability_no_analog"
  reason[lambda >= thresholds["speciation"]] <- "high_speciation_opportunity"
  reason[eps >= thresholds["extinction"]] <- "high_extinction_risk"
  reason[cnet < thresholds["connectivity"] & m >= thresholds["source_pressure"]] <-
    "isolation_or_barrier_limited"
  reason[m < thresholds["source_pressure"]] <- "source_limited"
  reason[s < thresholds["suitability"]] <- "niche_limited"
  reason[e < thresholds["lineage"]] <- "lineage_absent"
  reason[g < thresholds["geography"]] <- "arena_absent_or_lost"

  x$dominant_mechanism <- reason
  x$mechanism_class <- ifelse(g < thresholds["geography"], "earth_surface",
    ifelse(e < thresholds["lineage"], "lineage_time",
      ifelse(s < thresholds["suitability"], "ecological_filter",
        ifelse(m < thresholds["source_pressure"] |
                 cnet < thresholds["connectivity"], "connectivity_dispersal",
          ifelse(eps >= thresholds["extinction"], "persistence_extinction",
            ifelse(lambda >= thresholds["speciation"], "diversification",
              ifelse(q >= thresholds["no_analog"], "uncertainty",
                     "assembly")))))))
  x
}

#' Map regional historical accessibility to cells
#'
#' Converts BioGeoBEARS/hmscHist region-level accessibility
#' `A[j, region, time]` to cell-level accessibility `A[j, cell, time]` using a
#' palaeogeographic region cube or cell table. This is the safe direction for
#' deep-time prediction: region-level exogenous history is mapped to cells,
#' while modern site-level back-calculated history should not be projected to
#' new palaeo cells.
#'
#' @param accessibility Species-region-time accessibility table, or a table that
#'   already contains `cell_id`.
#' @param cell_table Cell-time table with region ids.
#' @param region_col,cell_col,time_col,species_col Column names.
#' @param accessibility_col Accessibility column.
#' @return Species-cell-time accessibility table.
#' @export
#'
#' @examples
#' A <- data.frame(species = "sp1", region = "R1", time_ma = 0,
#'   accessibility = 0.8)
#' cells <- data.frame(cell_id = 1:2, region = "R1", time_ma = 0)
#' hee_accessibility_to_cells(A, cells)
hee_accessibility_to_cells <- function(accessibility,
                                       cell_table,
                                       region_col = "region",
                                       cell_col = "cell_id",
                                       time_col = "time_ma",
                                       species_col = "species",
                                       accessibility_col = "accessibility") {
  A <- as.data.frame(accessibility)
  if (cell_col %in% names(A)) {
    .require_cols(A, c(species_col, cell_col, time_col, accessibility_col),
                  "cell-level accessibility")
    .hee_validate_time_values(A[[time_col]],
                              paste0("accessibility$", time_col))
    .hee_check_unique_keys(A, c(species_col, cell_col, time_col),
                           "cell-level accessibility")
    A[[accessibility_col]] <- .hee_clip01(A[[accessibility_col]])
    A[[accessibility_col]][is.na(A[[accessibility_col]])] <- 0
    return(A)
  }
  cells <- as.data.frame(cell_table)
  .require_cols(A, c(species_col, region_col, time_col, accessibility_col),
                "accessibility")
  .require_cols(cells, c(cell_col, region_col, time_col), "cell_table")
  .hee_validate_time_values(A[[time_col]], paste0("accessibility$", time_col))
  .hee_validate_time_values(cells[[time_col]], paste0("cell_table$", time_col))
  .hee_check_unique_keys(A, c(species_col, region_col, time_col),
                         "accessibility")
  .hee_check_unique_keys(cells, c(cell_col, time_col),
                         "cell_table")
  keep <- unique(c(cell_col, region_col, time_col, intersect(c("lon", "lat", "paleo_lon",
                                                               "paleo_lat"), names(cells))))
  out <- merge(A, cells[, keep, drop = FALSE], by = c(region_col, time_col),
               all.x = FALSE, sort = FALSE)
  out[[accessibility_col]] <- .hee_clip01(out[[accessibility_col]])
  out
}

#' Combine static deep-time probability components
#'
#' Static snapshot integration for
#' `P_static = S * A * E * G * D`. This is a named wrapper around
#' [hee_combine()] for the geological-process workflow. Dynamic colonisation and
#' extinction should use [hee_dynamic_assembly()] instead.
#'
#' @param suitability HMSC suitability table.
#' @param accessibility Optional accessibility table.
#' @param phylo_mask Optional phylogenetic/time mask.
#' @param geo_mask Optional geographic/habitat mask.
#' @param static_filter Optional static dispersal filter.
#' @param dynamic_filter Backwards-compatible dispersal filter argument. In this
#'   static wrapper it is treated as the single static dispersal term `D`, not as
#'   a second process to multiply.
#' @param dispersal_filter Optional generic dispersal filter. Recognised columns
#'   include `static_weight`, `D_static`, `dynamic_weight`, `D_dynamic`,
#'   `dispersal_weight`, or `weight`.
#' @param extrapolation_weight Optional reliability downweight. Keep `NULL`
#'   when extrapolation is reported as a flag.
#' @param method Combination method. Only `"multiply"` is supported.
#' @return Projection table with `probability`.
#' @export
hee_combine_static <- function(suitability,
                               accessibility = NULL,
                               phylo_mask = NULL,
                               geo_mask = NULL,
                               static_filter = NULL,
                               dynamic_filter = NULL,
                               dispersal_filter = NULL,
                               extrapolation_weight = NULL,
                               method = "multiply") {
  landmask <- geo_mask
  if (!is.null(landmask)) {
    landmask <- as.data.frame(landmask)
    land_cols <- intersect(c("geographic_existence", "habitat_availability",
                             "land", "L_land", "G_arena"), names(landmask))
    if (length(land_cols) == 0L) {
      stop("geo_mask must contain one of: geographic_existence, ",
           "habitat_availability, land, L_land, G_arena.", call. = FALSE)
    }
    vals <- rep(1, nrow(landmask))
    for (nm in land_cols) {
      v <- .hee_clip01(landmask[[nm]])
      v[is.na(v)] <- 0
      vals <- vals * v
    }
    landmask$land_weight <- .hee_clip01(vals)
  }
  D <- dispersal_filter %||% static_filter %||% dynamic_filter
  if (!is.null(D)) {
    D <- as.data.frame(D)
    candidate <- intersect(c("static_weight", "D_static", "dynamic_weight",
                             "D_dynamic", "dispersal_weight", "weight"),
                           names(D))
    if (length(candidate) > 0L) {
      D$static_weight <- .hee_clip01(D[[candidate[1]]])
    }
  }
  hee_combine(suitability = suitability,
              accessibility = accessibility,
              phylo_mask = phylo_mask,
              landmask = landmask,
              static_filter = D,
              dynamic_filter = NULL,
              extrapolation_weight = extrapolation_weight,
              method = method)
}

#' Summarise limiting mechanisms to map level
#'
#' Aggregates species-level mechanism attribution to cell-time or region-time
#' maps. The dominant community mechanism is the class with the largest weighted
#' count. This is a diagnostic map and should be interpreted as "most limiting
#' in the model components", not causal proof.
#'
#' @param attribution Output from [hee_mechanism_attribution()].
#' @param level `"community"` summarises over species; `"species"` returns the
#'   input with proportions by species/cell/time.
#' @param group_cols Grouping columns for map summaries.
#' @param mechanism_col Mechanism column.
#' @param weight_col Optional weight column.
#' @return Mechanism summary table.
#' @export
hee_limitation_map <- function(attribution,
                               level = c("community", "species"),
                               group_cols = c("cell_id", "time_ma"),
                               mechanism_col = "dominant_mechanism",
                               weight_col = NULL) {
  level <- match.arg(level)
  x <- as.data.frame(attribution)
  .require_cols(x, c(group_cols, mechanism_col), "attribution")
  x$.weight <- if (!is.null(weight_col) && weight_col %in% names(x)) {
    suppressWarnings(as.numeric(x[[weight_col]]))
  } else {
    rep(1, nrow(x))
  }
  x$.weight[is.na(x$.weight)] <- 0
  if (identical(level, "species")) return(x)
  form <- stats::as.formula(paste(".weight ~", paste(c(group_cols, mechanism_col), collapse = " + ")))
  tab <- stats::aggregate(form, x, sum, na.rm = TRUE)
  names(tab)[names(tab) == ".weight"] <- "mechanism_weight"
  total <- stats::aggregate(mechanism_weight ~ ., tab[, c(group_cols, "mechanism_weight"), drop = FALSE],
                            sum, na.rm = TRUE)
  names(total)[ncol(total)] <- "total_weight"
  tab <- merge(tab, total, by = group_cols, all.x = TRUE, sort = FALSE)
  tab$mechanism_fraction <- ifelse(tab$total_weight > 0, tab$mechanism_weight / tab$total_weight, NA_real_)
  pick <- do.call(rbind, lapply(split(tab, interaction(tab[, group_cols, drop = FALSE], drop = TRUE)),
                                function(z) z[which.max(z$mechanism_weight), , drop = FALSE]))
  names(pick)[names(pick) == mechanism_col] <- "dominant_community_mechanism"
  rownames(pick) <- NULL
  pick
}

#' Compute prediction reliability index
#'
#' Combines epistemic reliability components without treating them as ecological
#' processes. The default index is
#' `PRI = (1 - extrapolation_risk) * model_agreement * phylo_validity *
#' calibration`. If an extrapolation score is already in `[0, 1]` it is kept on
#' that scale; larger non-negative scores use the saturating transform
#' `risk = score / (score + 1)`. Missing optional components are neutral
#' (`1`) and component weights exponentiate each factor, so a weight of zero
#' removes that component.
#'
#' @param extrapolation Optional table with extrapolation/no-analog risk.
#' @param model_agreement Optional table/vector with model agreement.
#' @param phylo_validity Optional table/vector with lineage/species validity.
#' @param calibration Optional table/vector with validation/calibration score.
#' @param keys Join keys.
#' @param risk_col,agreement_col,validity_col,calibration_col Column names.
#' @param weights Named component weights for `extrapolation`,
#'   `model_agreement`, `phylo_validity`, and `calibration`.
#' @return Reliability table with `prediction_reliability`.
#' @export
hee_prediction_reliability <- function(extrapolation = NULL,
                                       model_agreement = NULL,
                                       phylo_validity = NULL,
                                       calibration = NULL,
                                       keys = c("cell_id", "time_ma"),
                                       risk_col = "extrapolation_score",
                                       agreement_col = "model_agreement",
                                       validity_col = "phylo_validity",
                                       calibration_col = "calibration",
                                       weights = c(
                                         extrapolation = 1,
                                         model_agreement = 1,
                                         phylo_validity = 1,
                                         calibration = 1
                                       )) {
  base <- .hee_first_table(extrapolation, model_agreement, phylo_validity, calibration)
  if (is.null(base)) {
    n <- max(vapply(list(extrapolation, model_agreement, phylo_validity, calibration),
                    function(z) {
                      if (is.null(z)) return(0L)
                      length(z)
                    }, integer(1)))
    if (n <= 0L) stop("Provide at least one reliability component.", call. = FALSE)
    lens <- vapply(list(extrapolation, model_agreement, phylo_validity, calibration),
                   function(z) {
                     if (is.null(z)) return(1L)
                     length(z)
                   }, integer(1))
    bad <- lens[!lens %in% c(1L, n)]
    if (length(bad) > 0L) {
      stop("Numeric reliability components must have length 1 or a common length.",
           call. = FALSE)
    }
    base <- data.frame(component_id = seq_len(n))
  }
  base <- as.data.frame(base)
  base$.row_id <- seq_len(nrow(base))
  risk <- .hee_metric_value(base, extrapolation, risk_col, default = 0,
                            keys = keys,
                            aliases = c("no_analog_score",
                                        "extrapolation_risk",
                                        "no_analog_fraction"),
                            source_name = "extrapolation")
  agreement <- .hee_metric_value(base, model_agreement, agreement_col,
                                 default = 1, keys = keys,
                                 aliases = c("agreement", "agreement_score"),
                                 source_name = "model_agreement")
  validity <- .hee_metric_value(base, phylo_validity, validity_col,
                                default = 1, keys = keys,
                                aliases = c("phylo_existence", "E_phylo",
                                            "lineage_exists",
                                            "phylo_present"),
                                source_name = "phylo_validity")
  calib <- .hee_metric_value(base, calibration, calibration_col, default = 1,
                             keys = keys,
                             aliases = c("calibration_score",
                                         "validation_score"),
                             source_name = "calibration")
  wt <- c(extrapolation = 1, model_agreement = 1,
          phylo_validity = 1, calibration = 1)
  if (!is.null(weights)) {
    if (is.null(names(weights))) {
      wt[seq_len(min(length(wt), length(weights)))] <-
        weights[seq_len(min(length(wt), length(weights)))]
    } else {
      hit <- intersect(names(weights), names(wt))
      wt[hit] <- weights[hit]
    }
  }
  wt[!is.finite(wt) | wt < 0] <- 0
  base$extrapolation_risk_raw <- suppressWarnings(as.numeric(risk))
  base$extrapolation_risk_raw[is.na(base$extrapolation_risk_raw)] <- 0
  base$extrapolation_risk <- .hee_normalize_risk(base$extrapolation_risk_raw)
  base$model_agreement <- .hee_clip01(agreement)
  base$phylo_validity <- .hee_clip01(validity)
  base$calibration <- .hee_clip01(calib)
  base$model_agreement[is.na(base$model_agreement)] <- 1
  base$phylo_validity[is.na(base$phylo_validity)] <- 1
  base$calibration[is.na(base$calibration)] <- 1
  base$prediction_reliability <- .hee_clip01(
    (1 - base$extrapolation_risk)^wt["extrapolation"] *
      base$model_agreement^wt["model_agreement"] *
      base$phylo_validity^wt["phylo_validity"] *
      base$calibration^wt["calibration"]
  )
  base[order(base$.row_id), setdiff(names(base), ".row_id"), drop = FALSE]
}

#' Analyse geological processes behind deep-time projections
#'
#' Builds a complete, transparent set of diagnostic tables for the geological
#' process layer of a deep-time projection. The function summarises landscape
#' arena dynamics, ecological opportunity, source pressure, rescue effects,
#' colonisation probability, extinction probability, speciation opportunity,
#' dynamic occupancy, mechanism attribution, reliability, and species-pool
#' change. These are diagnostic calculations and scenario indices; they are not
#' automatic causal estimates.
#'
#' @param landscape_state Cell/region-time table with palaeogeographic state.
#'   Recommended columns include `cell_id`, `region`, `time_ma`, `lon`, `lat`,
#'   `land_mask_dem` or `geographic_existence`, `land_area_km2`, and
#'   environmental variables.
#' @param projection Species-cell-time projection table, usually the M5
#'   projection, with `species`, `cell_id`, `time_ma`, `probability`, and
#'   `suitability`.
#' @param connectivity_cube Optional species- or region-specific connectivity
#'   table with `from_region`, `to_region`, `time_ma`, and `connectivity`.
#' @param accessibility Optional BioGeoBEARS/history accessibility table.
#' @param phylo_mask Optional species-by-time matrix or long table.
#' @param extrapolation Optional extrapolation/no-analog risk table.
#' @param tip_ranges Optional present-day species-by-region range matrix/table.
#' @param species_origin Optional named vector/data.frame of species origin
#'   ages.
#' @param bsm_events Optional BioGeoBEARS stochastic-map event table.
#' @param traits Optional trait table used only for reporting metadata.
#' @param times Time axis in Ma. Defaults to times in `projection`.
#' @param env_cols Optional environmental columns used to compute heterogeneity
#'   and novelty.
#' @param land_col Geographic existence column.
#' @param area_col Area column.
#' @param region_col,cell_col,time_col,species_col Column names.
#' @param initial_probability Starting occupancy probability for the oldest
#'   interval when no previous regional state is available.
#' @param connectivity_threshold Connectivity threshold used to accumulate
#'   isolation duration.
#'
#' @return A named list of data.frames covering process inputs, process
#'   probabilities, dynamic state, mechanism attribution, and reliability.
#' @export
#'
#' @examples
#' land <- data.frame(cell_id = c("c1", "c1", "c2", "c2"),
#'   region = c("A", "A", "B", "B"), time_ma = c(5, 0, 5, 0),
#'   land_mask_dem = 1, land_area_km2 = c(10, 12, 8, 9),
#'   bio1 = c(1, 2, 3, 4))
#' proj <- expand.grid(species = "sp1", cell_id = c("c1", "c2"),
#'   time_ma = c(5, 0), stringsAsFactors = FALSE)
#' proj <- merge(proj, land[, c("cell_id", "region", "time_ma")])
#' proj$suitability <- 0.7
#' proj$probability <- 0.4
#' hee_geoprocess_diagnostics(land, proj, env_cols = "bio1")$process_summary
hee_geoprocess_diagnostics <- function(landscape_state,
                                       projection,
                                       connectivity_cube = NULL,
                                       accessibility = NULL,
                                       phylo_mask = NULL,
                                       extrapolation = NULL,
                                       tip_ranges = NULL,
                                       species_origin = NULL,
                                       bsm_events = NULL,
                                       traits = NULL,
                                       times = NULL,
                                       env_cols = NULL,
                                       land_col = NULL,
                                       area_col = "land_area_km2",
                                       region_col = "region",
                                       cell_col = "cell_id",
                                       time_col = "time_ma",
                                       species_col = "species",
                                       initial_probability = 0,
                                       connectivity_threshold = 0.2) {
  land <- as.data.frame(landscape_state)
  proj <- as.data.frame(projection)
  .require_cols(land, c(cell_col, time_col), "landscape_state")
  .require_cols(proj, c(species_col, cell_col, time_col, "probability", "suitability"),
                "projection")
  .hee_validate_time_values(land[[time_col]], paste0("landscape_state$", time_col))
  .hee_validate_time_values(proj[[time_col]], paste0("projection$", time_col))
  .hee_check_unique_keys(land, c(cell_col, time_col), "landscape_state")
  .hee_check_unique_keys(proj, c(species_col, cell_col, time_col), "projection")

  if (!region_col %in% names(land)) land[[region_col]] <- land[[cell_col]]
  if (!region_col %in% names(proj)) {
    region_lookup <- unique(land[, c(cell_col, region_col, time_col), drop = FALSE])
    .hee_check_unique_keys(region_lookup, c(cell_col, time_col),
                           "landscape_state region lookup")
    proj <- merge(proj, region_lookup, by = c(cell_col, time_col),
                  all.x = TRUE, sort = FALSE)
  }
  if (!region_col %in% names(proj)) proj[[region_col]] <- proj[[cell_col]]

  times <- sort(unique(as.numeric(times %||% proj[[time_col]])), decreasing = TRUE)
  if (is.null(land_col)) {
    land_col <- intersect(c("geographic_existence", "land_mask_dem", "land",
                            "land_weight", "L_land", "G_arena"), names(land))[1]
  }
  if (is.na(land_col) || !land_col %in% names(land)) {
    land$geographic_existence <- 1
    land_col <- "geographic_existence"
  }
  if (!"geographic_existence" %in% names(land)) {
    land <- hee_land_age(land, unit_col = cell_col, time_col = time_col,
                         land_col = land_col)
  } else if (!any(c("land_age_myr", "land_age_ma") %in% names(land))) {
    land <- hee_land_age(land, unit_col = cell_col, time_col = time_col,
                         land_col = "geographic_existence")
  }

  eco <- hee_ecological_opportunity(land, env_cols = env_cols,
                                    region_col = region_col,
                                    time_col = time_col,
                                    area_col = area_col)
  eco_cols <- intersect(c("geographic_existence", "land_age_myr",
                         "land_age_ma", "ecological_opportunity",
                         "area_gain_index", "heterogeneity_index",
                         "novelty_index", "young_land_index"),
                       names(eco))
  land_region <- .hee_aggregate_numeric(eco, c(region_col, time_col), eco_cols)
  if (area_col %in% names(eco)) {
    area_region <- stats::aggregate(eco[[area_col]],
                                    eco[, c(region_col, time_col), drop = FALSE],
                                    sum, na.rm = TRUE)
    names(area_region)[ncol(area_region)] <- area_col
    land_region <- merge(land_region, area_region, by = c(region_col, time_col),
                         all = TRUE, sort = FALSE)
  }
  coord_cols <- intersect(c("lon", "lat", "paleo_lon", "paleo_lat"), names(eco))
  if (length(coord_cols) > 0L) {
    coords <- .hee_aggregate_numeric(eco, c(region_col, time_col), coord_cols)
    land_region <- merge(land_region, coords, by = c(region_col, time_col),
                         all.x = TRUE, sort = FALSE)
  }
  land_region <- .hee_add_interval_and_land_loss(land_region, region_col,
                                                 time_col, area_col)

  proj_cols <- intersect(c("probability", "suitability", "accessibility",
                           "phylo_existence", "land_weight", "land_mask_dem",
                           "D_static", "D_dynamic", "static_weight",
                           "dynamic_weight", "Q_extrap",
                           "extrapolation_score"),
                         names(proj))
  region_state <- .hee_aggregate_numeric(proj,
                                         c(species_col, region_col, time_col),
                                         proj_cols)
  region_state <- merge(region_state, land_region, by = c(region_col, time_col),
                        all.x = TRUE, sort = FALSE)
  region_state <- .hee_add_interval_and_previous(region_state, species_col,
                                                 region_col, time_col,
                                                 "probability",
                                                 initial_probability)
  if (!"accessibility" %in% names(region_state)) {
    region_state$accessibility <- .hee_metric_value(
      region_state, accessibility, "accessibility",
      default = if (is.null(accessibility)) 1 else 0,
      keys = c(species_col, region_col, time_col),
      aliases = c("A_hist", "A_BGB", "weight"),
      source_name = "accessibility"
    )
  }
  if (!"phylo_existence" %in% names(region_state)) {
    region_state$phylo_existence <- .hee_phylo_value(region_state, phylo_mask,
                                                     species_col, time_col)
  }
  if (!"land_mask_dem" %in% names(region_state)) {
    region_state$land_mask_dem <- .hee_col_or_default(region_state,
                                                      "geographic_existence", 1)
  }
  if (!"D_static" %in% names(region_state)) {
    region_state$D_static <- .hee_col_or_default(region_state, "static_weight", 1)
  }
  if (!"D_dynamic" %in% names(region_state)) {
    region_state$D_dynamic <- .hee_col_or_default(region_state, "dynamic_weight", 1)
  }
  region_state$geographic_existence[is.na(region_state$geographic_existence)] <-
    region_state$land_mask_dem[is.na(region_state$geographic_existence)]

  conn <- if (is.null(connectivity_cube)) NULL else
    hee_isolation_history(connectivity_cube,
                          threshold = connectivity_threshold,
                          time_col = time_col)
  conn_to <- .hee_connectivity_to_region(conn, species_col, region_col, time_col)
  if (!is.null(conn_to)) {
    region_state <- .hee_merge_species_or_region(region_state, conn_to,
                                                 species_col, region_col,
                                                 time_col)
  }
  if (!"connectivity" %in% names(region_state)) {
    region_state$connectivity <- .hee_clip01(region_state$D_static *
                                               region_state$D_dynamic)
  }
  if (!"isolation" %in% names(region_state)) {
    region_state$isolation <- .hee_clip01(1 - region_state$connectivity)
  }
  if (!"isolation_duration_ma" %in% names(region_state)) {
    region_state$isolation_duration_ma <- 0
  }

  previous_state <- region_state[, c(species_col, region_col, time_col,
                                     "previous_probability"), drop = FALSE]
  names(previous_state)[names(previous_state) == "previous_probability"] <-
    "occupancy_probability"
  target_region <- unique(region_state[, c(species_col, region_col, time_col),
                                       drop = FALSE])
  if (is.null(conn)) {
    source_pressure <- target_region
    source_pressure$source_pressure <- 0
    source_pressure$n_source_regions <- 0L
    source_pressure$max_source_probability <- 0
    rescue <- source_pressure
    names(rescue)[names(rescue) == "source_pressure"] <- "rescue_effect"
  } else {
    source_pressure <- hee_colonisation_pressure(
      previous_state = previous_state,
      connectivity_cube = conn,
      target = target_region,
      species_col = species_col,
      region_col = region_col,
      from_col = "from_region",
      to_col = "to_region",
      time_col = time_col
    )
    rescue <- hee_rescue_effect(
      previous_state = previous_state,
      connectivity_cube = conn,
      target = target_region,
      species_col = species_col,
      region_col = region_col,
      from_col = "from_region",
      to_col = "to_region",
      time_col = time_col
    )
  }
  region_state <- merge(region_state, source_pressure,
                        by = c(species_col, region_col, time_col),
                        all.x = TRUE, sort = FALSE)
  region_state <- merge(region_state,
                        rescue[, unique(c(species_col, region_col, time_col,
                                          "rescue_effect")), drop = FALSE],
                        by = c(species_col, region_col, time_col),
                        all.x = TRUE, sort = FALSE)
  region_state$source_pressure[is.na(region_state$source_pressure)] <- 0
  region_state$rescue_effect[is.na(region_state$rescue_effect)] <- 0

  colonisation <- hee_colonisation_probability(
    region_state,
    source_pressure = region_state,
    accessibility = region_state,
    ecological_opportunity = region_state,
    delta_t = "interval_myr",
    keys = c(species_col, region_col, time_col),
    suitability_col = "suitability"
  )
  extinction_input <- region_state
  if (!"land_loss" %in% names(extinction_input)) extinction_input$land_loss <- 0
  extinction <- hee_extinction_probability(
    extinction_input,
    delta_t = "interval_myr",
    suitability_col = "suitability",
    land_loss_col = "land_loss",
    isolation_col = "isolation",
    area_col = area_col,
    rescue_col = "rescue_effect"
  )
  speciation_input <- region_state
  speciation_input$niche_divergence <- .hee_clip01(abs(
    speciation_input$suitability - speciation_input$previous_probability
  ))
  speciation_input$population_persistence <- .hee_clip01(
    pmax(speciation_input$probability, speciation_input$previous_probability,
         na.rm = TRUE)
  )
  speciation_input$range_fragmentation <- .hee_clip01(speciation_input$isolation)
  speciation_input$rare_dispersal <- .hee_clip01(speciation_input$source_pressure *
                                                  (1 - speciation_input$connectivity))
  speciation_input$isolation_after_arrival <- .hee_clip01(speciation_input$isolation)
  speciation_input$colonisation_input <- .hee_clip01(speciation_input$source_pressure)
  if (!"richness" %in% names(speciation_input)) {
    rich_region <- stats::aggregate(probability ~ region + time_ma,
                                    setNames(region_state[, c(region_col, time_col,
                                                              "probability")],
                                             c("region", "time_ma", "probability")),
                                    sum, na.rm = TRUE)
    names(rich_region)[1:2] <- c(region_col, time_col)
    names(rich_region)[3] <- "richness"
    speciation_input <- merge(speciation_input, rich_region,
                              by = c(region_col, time_col),
                              all.x = TRUE, sort = FALSE)
  }
  speciation <- hee_speciation_opportunity(speciation_input,
                                           area_col = area_col,
                                           land_age_col = "land_age_ma")
  speciation <- hee_allopatric_speciation_opportunity(speciation)
  speciation <- hee_founder_speciation_opportunity(speciation)
  speciation <- hee_insitu_speciation_opportunity(speciation)
  speciation <- hee_radiation_opportunity(speciation)
  speciation_component_cols <- intersect(
    c("allopatric_opportunity", "founder_event_opportunity",
      "in_situ_speciation_index", "radiation_opportunity"),
    names(speciation)
  )
  if (length(speciation_component_cols) > 0L) {
    speciation$speciation_opportunity <- .hee_clip01(rowMeans(
      speciation[, speciation_component_cols, drop = FALSE],
      na.rm = TRUE
    ))
  }

  dyn_suit <- region_state[, unique(c(species_col, region_col, time_col,
                                      "suitability")), drop = FALSE]
  names(dyn_suit)[names(dyn_suit) == region_col] <- cell_col
  dyn_col <- colonisation[, unique(c(species_col, region_col, time_col,
                                     "colonisation_probability")), drop = FALSE]
  names(dyn_col)[names(dyn_col) == region_col] <- cell_col
  dyn_ext <- extinction[, unique(c(species_col, region_col, time_col,
                                   "extinction_probability")), drop = FALSE]
  names(dyn_ext)[names(dyn_ext) == region_col] <- cell_col
  dyn_geo <- unique(region_state[, c(region_col, time_col,
                                     "geographic_existence"), drop = FALSE])
  names(dyn_geo)[names(dyn_geo) == region_col] <- cell_col
  dynamic_occupancy <- hee_dynamic_assembly(
    suitability = dyn_suit,
    colonisation = dyn_col,
    extinction = dyn_ext,
    geo_mask = dyn_geo,
    phylo_mask = phylo_mask,
    initial_probability = initial_probability,
    species_col = species_col,
    cell_col = cell_col,
    time_col = time_col
  )
  names(dynamic_occupancy)[names(dynamic_occupancy) == cell_col] <- region_col

  attr_input <- merge(dynamic_occupancy, region_state,
                      by = c(species_col, region_col, time_col),
                      all.x = TRUE, sort = FALSE,
                      suffixes = c("_dynamic", ""))
  attr_input <- merge(attr_input,
                      extinction[, unique(c(species_col, region_col, time_col,
                                            "extinction_probability")), drop = FALSE],
                      by = c(species_col, region_col, time_col),
                      all.x = TRUE, sort = FALSE, suffixes = c("", ".ext"))
  if ("extinction_probability.ext" %in% names(attr_input)) {
    attr_input$extinction_probability <- attr_input$extinction_probability.ext
  }
  attr_input <- merge(attr_input,
                      speciation[, unique(c(species_col, region_col, time_col,
                                            "speciation_opportunity")), drop = FALSE],
                      by = c(species_col, region_col, time_col),
                      all.x = TRUE, sort = FALSE)
  if (!"no_analog_score" %in% names(attr_input)) {
    attr_input$no_analog_score <- .hee_col_or_default(attr_input,
                                                      "extrapolation_score", 0)
  }
  mechanism_attribution <- hee_mechanism_attribution(
    attr_input,
    columns = c(geography = "geographic_existence",
                lineage = "phylo_existence",
                suitability = "suitability",
                source_pressure = "source_pressure",
                connectivity = "connectivity",
                extinction = "extinction_probability",
                speciation = "speciation_opportunity",
                no_analog = "no_analog_score")
  )
  limitation_summary <- hee_limitation_map(
    mechanism_attribution,
    group_cols = c(region_col, time_col),
    weight_col = "occupancy_probability"
  )

  pool <- .hee_species_pool_from_processes(region_state, speciation, colonisation,
                                           extinction, bsm_events, tip_ranges,
                                           times, region_col, time_col)

  extrap_region <- NULL
  if (!is.null(extrapolation)) {
    ex <- as.data.frame(extrapolation)
    if (!region_col %in% names(ex) && cell_col %in% names(ex)) {
      lookup <- unique(land[, c(cell_col, region_col, time_col), drop = FALSE])
      ex <- merge(ex, lookup, by = c(cell_col, time_col), all.x = TRUE,
                  sort = FALSE)
    }
    ex_cols <- intersect(c("extrapolation_score", "extrapolation_flag",
                           "Q_extrap", "no_analog_fraction"),
                         names(ex))
    if (length(ex_cols) > 0L && all(c(region_col, time_col) %in% names(ex))) {
      extrap_region <- .hee_aggregate_numeric(ex, c(region_col, time_col), ex_cols)
    }
  }
  reliability <- NULL
  if (!is.null(extrap_region)) {
    model_agreement_tbl <- extrap_region[, c(region_col, time_col), drop = FALSE]
    model_agreement_tbl$model_agreement <- 1 - .hee_clip01(
      .hee_col_or_default(extrap_region, "extrapolation_score", 0) * 0.25
    )
    phylo_validity_tbl <- extrap_region[, c(region_col, time_col), drop = FALSE]
    phylo_validity_tbl$phylo_validity <- 1
    reliability <- hee_prediction_reliability(
      extrapolation = extrap_region,
      model_agreement = model_agreement_tbl,
      phylo_validity = phylo_validity_tbl,
      keys = c(region_col, time_col)
    )
  }

  source_sink <- .hee_source_sink_summary(conn, source_pressure, rescue,
                                          region_col, time_col, species_col)
  process_cols <- intersect(c("probability", "suitability", "accessibility",
                              "phylo_existence", "geographic_existence",
                              "connectivity", "isolation",
                              "isolation_duration_ma", "source_pressure",
                              "rescue_effect", "land_loss",
                              "ecological_opportunity"),
                            names(region_state))
  process_summary <- .hee_aggregate_numeric(region_state, time_col, process_cols)
  process_summary <- merge(process_summary,
                           .hee_aggregate_numeric(colonisation, time_col,
                                                  "colonisation_probability"),
                           by = time_col, all = TRUE, sort = FALSE)
  process_summary <- merge(process_summary,
                           .hee_aggregate_numeric(extinction, time_col,
                                                  "extinction_probability"),
                           by = time_col, all = TRUE, sort = FALSE)
  process_summary <- merge(process_summary,
                           .hee_aggregate_numeric(speciation, time_col,
                                                  c("speciation_opportunity",
                                                    "vicariance_opportunity",
                                                    "founder_event_opportunity",
                                                    "in_situ_speciation_index",
                                                    "radiation_opportunity")),
                           by = time_col, all = TRUE, sort = FALSE)
  process_summary <- merge(process_summary,
                           .hee_aggregate_numeric(dynamic_occupancy, time_col,
                                                  "occupancy_probability"),
                           by = time_col, all = TRUE, sort = FALSE)
  if (!is.null(reliability)) {
    process_summary <- merge(process_summary,
                             .hee_aggregate_numeric(reliability, time_col,
                                                    "prediction_reliability"),
                             by = time_col, all = TRUE, sort = FALSE)
  }

  list(
    landscape_state = land,
    ecological_opportunity = eco,
    landscape_region = land_region,
    region_state = region_state,
    connectivity_region = conn_to %||% data.frame(),
    source_pressure = source_pressure,
    rescue_effect = rescue,
    colonisation_probability = colonisation,
    extinction_probability = extinction,
    speciation_opportunity = speciation,
    dynamic_occupancy = dynamic_occupancy,
    mechanism_attribution = mechanism_attribution,
    limitation_summary = limitation_summary,
    species_pool = pool,
    source_sink_summary = source_sink,
    reliability = reliability %||% data.frame(),
    process_summary = process_summary,
    metadata = data.frame(
      item = c("n_species", "n_regions", "n_times", "n_traits",
               "has_connectivity", "has_extrapolation", "interpretation"),
      value = c(length(unique(proj[[species_col]])),
                length(unique(proj[[region_col]])),
                length(times),
                if (is.null(traits)) 0 else ncol(as.data.frame(traits)),
                !is.null(connectivity_cube),
                !is.null(extrapolation),
                "diagnostic_indices_not_causal_estimates"),
      stringsAsFactors = FALSE
    )
  )
}

#' Convert HMSC Beta to niche-trait diagnostics
#'
#' Workflow alias for [calc_niche_metrics()].
#'
#' @param fit A `hmsc_ecoevo` object or object accepted by
#'   [calc_niche_metrics()].
#' @param ... Arguments passed to [calc_niche_metrics()].
#' @return A `hmsc_niche_metrics` object.
#' @export
hee_beta_to_niche_traits <- function(fit, ...) {
  calc_niche_metrics(fit, ...)
}

#' Compute phylogenetic niche signal
#'
#' Workflow alias for [calc_phylo_signal_metrics()].
#'
#' @param evo A `hmsc_ecoevo` object.
#' @param ... Arguments passed to [calc_phylo_signal_metrics()].
#' @return A `hmsc_phylo_signal_metrics` object.
#' @export
hee_phylo_niche_signal <- function(evo, ...) {
  calc_phylo_signal_metrics(evo, ...)
}

#' Compute trait-mediated niche signal
#'
#' Workflow alias for [calc_trait_mediation_metrics()].
#'
#' @param evo A `hmsc_ecoevo` object.
#' @param ... Arguments passed to [calc_trait_mediation_metrics()].
#' @return A `hmsc_trait_mediation_metrics` object.
#' @export
hee_trait_mediation_signal <- function(evo, ...) {
  calc_trait_mediation_metrics(evo, ...)
}

#' Compute Gamma shift diagnostics
#'
#' Workflow alias for [calc_gamma_evolution_metrics()].
#'
#' @param evo A `hmsc_ecoevo` object.
#' @param ... Arguments passed to [calc_gamma_evolution_metrics()].
#' @return A `hmsc_gamma_evolution_metrics` object.
#' @export
hee_gamma_shift <- function(evo, ...) {
  calc_gamma_evolution_metrics(evo, ...)
}

# Internal helpers ---------------------------------------------------------

.hee_aggregate_numeric <- function(x, group_cols, value_cols) {
  x <- as.data.frame(x)
  group_cols <- intersect(group_cols, names(x))
  value_cols <- intersect(value_cols, names(x))
  if (length(group_cols) == 0L || length(value_cols) == 0L) {
    return(data.frame())
  }
  for (nm in value_cols) x[[nm]] <- suppressWarnings(as.numeric(x[[nm]]))
  stats::aggregate(x[, value_cols, drop = FALSE],
                   x[, group_cols, drop = FALSE],
                   mean, na.rm = TRUE)
}

.hee_add_interval_and_land_loss <- function(x, region_col, time_col,
                                            area_col = "land_area_km2") {
  if (nrow(x) == 0L) return(x)
  x$.row_id <- seq_len(nrow(x))
  parts <- split(x, as.character(x[[region_col]]), drop = TRUE)
  out <- lapply(parts, function(z) {
    z <- z[order(z[[time_col]], decreasing = TRUE), , drop = FALSE]
    z$interval_myr <- c(0, abs(diff(z[[time_col]])))
    if (area_col %in% names(z)) {
      area <- suppressWarnings(as.numeric(z[[area_col]]))
      prev_area <- c(area[1], utils::head(area, -1))
      z$land_area_previous <- prev_area
      z$land_loss <- ifelse(prev_area > 0, pmax(prev_area - area, 0) / prev_area, 0)
    } else {
      z$land_loss <- 0
    }
    if ("geographic_existence" %in% names(z)) {
      z$land_loss[.hee_clip01(z$geographic_existence) <= 0] <- 1
    }
    z$land_loss <- .hee_clip01(z$land_loss)
    z
  })
  out <- do.call(rbind, out)
  out <- out[order(out$.row_id), setdiff(names(out), ".row_id"), drop = FALSE]
  rownames(out) <- NULL
  out
}

.hee_dynamic_intervals <- function(target, previous_state, species_col, cell_col,
                                   time_col, time_direction) {
  if (nrow(target) == 0L) return(numeric(0))
  out <- rep(0, nrow(target))
  group_key <- interaction(target[[species_col]], target[[cell_col]],
                           drop = TRUE)
  rows_by_group <- split(seq_len(nrow(target)), group_key)
  prev_key <- if (is.null(previous_state)) NULL else
    interaction(previous_state[[species_col]], previous_state[[cell_col]],
                drop = TRUE)
  for (key in names(rows_by_group)) {
    idx <- rows_by_group[[key]]
    ord <- order(target[[time_col]][idx],
                 decreasing = identical(time_direction, "forward"))
    idx_ord <- idx[ord]
    tt <- as.numeric(target[[time_col]][idx_ord])
    dt <- c(0, abs(diff(tt)))
    if (!is.null(previous_state) && key %in% as.character(prev_key)) {
      src <- previous_state[prev_key == key, , drop = FALSE]
      first_time <- tt[1]
      src_time <- as.numeric(src[[time_col]])
      if (!any(src_time == first_time)) {
        candidates <- if (identical(time_direction, "forward")) {
          src_time[src_time > first_time]
        } else {
          src_time[src_time < first_time]
        }
        if (length(candidates) > 0L) {
          dt[1] <- min(abs(candidates - first_time), na.rm = TRUE)
        }
      }
    }
    out[idx_ord] <- .hee_validate_delta_t(dt, "interval_myr")
  }
  out
}

.hee_add_interval_and_previous <- function(x, species_col, region_col, time_col,
                                           value_col,
                                           initial_probability = 0) {
  if (nrow(x) == 0L) return(x)
  x$.row_id <- seq_len(nrow(x))
  parts <- split(x, interaction(x[[species_col]], x[[region_col]], drop = TRUE),
                 drop = TRUE)
  out <- lapply(parts, function(z) {
    z <- z[order(z[[time_col]], decreasing = TRUE), , drop = FALSE]
    z$interval_myr <- c(0, abs(diff(z[[time_col]])))
    val <- .hee_clip01(z[[value_col]])
    z$previous_probability <- c(as.numeric(initial_probability)[1],
                                utils::head(val, -1))
    z$previous_probability[is.na(z$previous_probability)] <-
      as.numeric(initial_probability)[1]
    z
  })
  out <- do.call(rbind, out)
  out <- out[order(out$.row_id), setdiff(names(out), ".row_id"), drop = FALSE]
  rownames(out) <- NULL
  out
}

.hee_phylo_value <- function(target, phylo_mask, species_col, time_col) {
  if (is.null(phylo_mask)) return(rep(1, nrow(target)))
  if (is.matrix(phylo_mask)) {
    times <- attr(phylo_mask, "times")
    if (is.null(times)) times <- as.numeric(gsub("Ma$", "", colnames(phylo_mask)))
    out <- rep(0, nrow(target))
    sp <- match(as.character(target[[species_col]]), rownames(phylo_mask))
    tt <- match(target[[time_col]], times)
    ok <- !is.na(sp) & !is.na(tt)
    out[ok] <- phylo_mask[cbind(sp[ok], tt[ok])]
    return(out)
  }
  .hee_metric_value(target, phylo_mask, "phylo_existence", default = 0,
                    keys = c(species_col, time_col),
                    aliases = c("E_phylo", "phylo_present",
                                "lineage_exists"),
                    source_name = "phylo_mask")
}

.hee_connectivity_to_region <- function(conn, species_col, region_col, time_col) {
  if (is.null(conn) || nrow(conn) == 0L) return(NULL)
  x <- as.data.frame(conn)
  if (!all(c("to_region", time_col) %in% names(x))) return(NULL)
  vals <- intersect(c("structural_connectivity", "functional_connectivity",
                      "climate_connectivity", "connectivity", "isolation",
                      "isolation_duration_ma"),
                    names(x))
  if (length(vals) == 0L) return(NULL)
  group_cols <- c(intersect(species_col, names(x)), "to_region", time_col)
  out <- .hee_aggregate_numeric(x, group_cols, vals)
  names(out)[names(out) == "to_region"] <- region_col
  out
}

.hee_merge_species_or_region <- function(target, source, species_col, region_col,
                                         time_col) {
  if (is.null(source) || nrow(source) == 0L) return(target)
  src <- as.data.frame(source)
  if (species_col %in% names(src)) {
    merge(target, src, by = c(species_col, region_col, time_col),
          all.x = TRUE, sort = FALSE, suffixes = c("", ".conn"))
  } else {
    merge(target, src, by = c(region_col, time_col),
          all.x = TRUE, sort = FALSE, suffixes = c("", ".conn"))
  }
}

.hee_species_pool_from_processes <- function(region_state, speciation,
                                             colonisation, extinction,
                                             bsm_events, tip_ranges, times,
                                             region_col, time_col) {
  regions <- sort(unique(as.character(region_state[[region_col]])))
  if (!is.null(tip_ranges)) {
    if (is.matrix(tip_ranges) || is.data.frame(tip_ranges)) {
      m <- as.matrix(tip_ranges)
      suppressWarnings(storage.mode(m) <- "numeric")
      common <- intersect(regions, colnames(m))
      init <- data.frame(region = regions, species_pool_size = 0)
      if (length(common) > 0L) {
        init$species_pool_size[match(common, init$region)] <-
          colSums(m[, common, drop = FALSE] > 0, na.rm = TRUE)
      }
    } else {
      init <- data.frame(region = regions, species_pool_size = 1)
    }
  } else {
    present_time <- min(region_state[[time_col]], na.rm = TRUE)
    cur <- region_state[region_state[[time_col]] == present_time, , drop = FALSE]
    init <- stats::aggregate(probability ~ ., cur[, c(region_col, "probability"),
                                                  drop = FALSE],
                             sum, na.rm = TRUE)
    names(init) <- c("region", "species_pool_size")
  }
  names(init)[names(init) == "region"] <- region_col

  events <- data.frame()
  make_events <- function(tab, value, event_type) {
    if (is.null(tab) || !value %in% names(tab)) return(data.frame())
    z <- .hee_aggregate_numeric(tab, c(region_col, time_col), value)
    if (nrow(z) == 0L) return(z)
    out <- data.frame(
      region = z[[region_col]],
      time_ma = z[[time_col]],
      event_type = event_type,
      probability = .hee_clip01(z[[value]]),
      from_region = NA_character_,
      to_region = NA_character_,
      stringsAsFactors = FALSE
    )
    names(out)[1:2] <- c(region_col, time_col)
    out
  }
  events <- rbind(
    make_events(speciation, "speciation_opportunity", "speciation"),
    make_events(colonisation, "colonisation_probability", "colonisation"),
    make_events(extinction, "extinction_probability", "local_extinction")
  )
  if (!is.null(bsm_events) && nrow(as.data.frame(bsm_events)) > 0L) {
    bsm <- as.data.frame(bsm_events)
    if (all(c(time_col, "event_type") %in% names(bsm))) {
      if (!"probability" %in% names(bsm)) {
        bsm$probability <- if ("event_prob" %in% names(bsm)) bsm$event_prob else 1
      }
      if (region_col %in% names(bsm)) {
        keep <- c(region_col, time_col, "event_type", "probability",
                  "from_region", "to_region")
        for (nm in setdiff(keep, names(bsm))) bsm[[nm]] <- NA
        events <- rbind(events, bsm[, keep, drop = FALSE])
      }
    }
  }
  names(events)[names(events) == "region"] <- region_col
  names(events)[names(events) == "time_ma"] <- time_col
  hee_update_species_pool(init, events = events, times = times,
                          region_col = region_col, time_col = time_col)
}

.hee_source_sink_summary <- function(conn, source_pressure, rescue, region_col,
                                     time_col, species_col) {
  ss <- .hee_aggregate_numeric(source_pressure, c(region_col, time_col),
                               c("source_pressure", "n_source_regions",
                                 "max_source_probability"))
  rr <- .hee_aggregate_numeric(rescue, c(region_col, time_col), "rescue_effect")
  out <- if (nrow(ss) > 0L && nrow(rr) > 0L) {
    merge(ss, rr, by = c(region_col, time_col), all = TRUE, sort = FALSE)
  } else if (nrow(ss) > 0L) ss else rr
  if (!is.null(conn) && nrow(as.data.frame(conn)) > 0L &&
      all(c("from_region", "to_region", time_col) %in% names(conn))) {
    flow <- .hee_aggregate_numeric(conn, c("from_region", "to_region", time_col),
                                   intersect(c("connectivity", "isolation"),
                                             names(conn)))
    names(flow)[names(flow) == "from_region"] <- "source_region"
    names(flow)[names(flow) == "to_region"] <- "sink_region"
    attr(out, "pair_flow") <- flow
  }
  out
}

.hee_binary <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), NA_real_, ifelse(x > 0, 1, 0))
}

.hee_clip01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- NA_real_
  pmin(pmax(x, 0), 1)
}

.hee_rescale01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(rep(0, length(x)))
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (!any(ok)) return(rep(0, length(x)))
  if (all(x[ok] >= 0 & x[ok] <= 1)) {
    out[ok] <- x[ok]
    out[!ok] <- 0
    return(.hee_clip01(out))
  }
  if (all(x[ok] >= 0)) {
    z <- pmax(x[ok], 0)
    if (all(z == 0)) {
      out[ok] <- 0
      out[!ok] <- 0
      return(out)
    }
    if (length(unique(z)) == 1L) {
      scale <- stats::median(z[z > 0], na.rm = TRUE)
      if (!is.finite(scale) || scale <= 0) scale <- 1
      out[ok] <- z / (z + scale)
      out[!ok] <- 0
      return(.hee_clip01(out))
    }
  }
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(diff(rng)) || diff(rng) == 0) return(rep(0, length(x)))
  out <- (x - rng[1]) / diff(rng)
  out[!is.finite(out) | is.na(out)] <- 0
  .hee_clip01(out)
}

.hee_normalize_risk <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x) | x < 0 | is.nan(x)] <- 0
  x[is.infinite(x) & x > 0] <- .Machine$double.xmax
  out <- ifelse(x <= 1, x, x / (x + 1))
  .hee_clip01(out)
}

.hee_standardize_finite <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep(0, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)
  s <- stats::sd(x[ok], na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(out)
  out[ok] <- (x[ok] - mean(x[ok], na.rm = TRUE)) / s
  out
}

.hee_positive_index <- function(x, scale = NULL) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep(0, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)
  z <- pmax(x[ok], 0)
  if (all(z <= 1, na.rm = TRUE)) {
    out[ok] <- z
    return(.hee_clip01(out))
  }
  if (is.null(scale)) {
    scale <- stats::median(z[z > 0], na.rm = TRUE)
  }
  if (!is.finite(scale) || scale <= 0) scale <- 1
  out[ok] <- z / (z + scale)
  .hee_clip01(out)
}

.hee_area_inverse_index <- function(area, scale = NULL) {
  area <- suppressWarnings(as.numeric(area))
  out <- rep(0, length(area))
  finite <- is.finite(area)
  if (!any(finite)) return(out)
  positive <- finite & area > 0
  if (is.null(scale)) {
    scale <- stats::median(area[positive], na.rm = TRUE)
  }
  if (!is.finite(scale) || scale <= 0) scale <- 1
  out[positive] <- scale / (area[positive] + scale)
  out[finite & area <= 0] <- 1
  .hee_clip01(out)
}

.hee_named_coefficients <- function(x, names_required) {
  out <- setNames(rep(0, length(names_required)), names_required)
  if (is.null(names(x))) {
    out[seq_len(min(length(out), length(x)))] <- x[seq_len(min(length(out), length(x)))]
  } else {
    hit <- intersect(names(x), names_required)
    out[hit] <- x[hit]
  }
  out
}

.hee_scale_for <- function(species, scale) {
  if (length(scale) == 1L) return(rep(as.numeric(scale), length(species)))
  val <- scale[as.character(species)]
  fallback <- stats::median(as.numeric(scale), na.rm = TRUE)
  val[is.na(val)] <- fallback
  as.numeric(val)
}

.hee_numeric_values <- function(x, label = "value") {
  raw <- x
  if (is.factor(raw)) raw <- as.character(raw)
  v <- suppressWarnings(as.numeric(raw))
  bad <- is.na(v) & !is.na(raw)
  if (any(bad)) {
    stop(label, " must be numeric or coercible to numeric.", call. = FALSE)
  }
  v
}

.hee_validate_time_values <- function(x, label = "time_ma") {
  v <- .hee_numeric_values(x, label)
  if (length(v) == 0L) return(invisible(TRUE))
  if (any(!is.finite(v))) {
    stop(label, " must contain finite ages in Ma.", call. = FALSE)
  }
  if (any(v < 0)) {
    stop(label, " must not contain negative ages in Ma.", call. = FALSE)
  }
  invisible(TRUE)
}

.hee_validate_delta_t <- function(dt, label = "delta_t") {
  dt <- suppressWarnings(as.numeric(dt))
  if (length(dt) == 0L) return(dt)
  if (any(!is.finite(dt) | is.na(dt))) {
    stop(label, " must contain finite non-negative time intervals.", call. = FALSE)
  }
  if (any(dt < 0)) {
    stop(label, " must be non-negative. Use absolute adjacent time differences ",
         "after sorting time_ma from older to younger.", call. = FALSE)
  }
  dt
}

.hee_duplicate_key_table <- function(x, keys) {
  if (length(keys) == 0L || nrow(x) == 0L) return(x[FALSE, , drop = FALSE])
  dup <- duplicated(x[, keys, drop = FALSE]) |
    duplicated(x[, keys, drop = FALSE], fromLast = TRUE)
  unique(x[dup, keys, drop = FALSE])
}

.hee_format_duplicate_keys <- function(x, keys, max_rows = 5L) {
  d <- .hee_duplicate_key_table(x, keys)
  if (nrow(d) == 0L) return("")
  d <- utils::head(d, max_rows)
  apply(d, 1, function(r) paste(paste(keys, r, sep = "="), collapse = ", "))
}

.hee_check_unique_keys <- function(x, keys, table_name = "table") {
  keys <- intersect(keys, names(x))
  if (length(keys) == 0L || nrow(x) == 0L) return(invisible(TRUE))
  bad <- .hee_duplicate_key_table(x, keys)
  if (nrow(bad) > 0L) {
    ex <- .hee_format_duplicate_keys(x, keys)
    stop("Duplicate key(s) in ", table_name, " for ",
         paste(keys, collapse = ", "), ": ", paste(ex, collapse = " | "),
         call. = FALSE)
  }
  invisible(TRUE)
}

.hee_safe_left_join <- function(target, source, keys, value_cols,
                                target_name = "target",
                                source_name = "source") {
  if (is.null(source)) return(data.frame(.row_id = target$.row_id))
  x <- as.data.frame(source)
  keys <- intersect(keys, names(x))
  if (length(keys) == 0L) return(data.frame(.row_id = target$.row_id))
  values <- intersect(value_cols, names(x))
  if (length(values) == 0L) return(data.frame(.row_id = target$.row_id))
  .hee_check_unique_keys(x, keys, source_name)
  .hee_check_unique_keys(target[, unique(c(".row_id", keys)), drop = FALSE],
                         ".row_id", target_name)
  src <- x[, unique(c(keys, values)), drop = FALSE]
  out <- merge(target[, c(".row_id", keys), drop = FALSE], src, by = keys,
               all.x = TRUE, sort = FALSE)
  if (nrow(out) != nrow(target)) {
    stop("Join from ", source_name, " to ", target_name,
         " changed row count from ", nrow(target), " to ", nrow(out),
         ". Check keys: ", paste(keys, collapse = ", "), call. = FALSE)
  }
  out[match(target$.row_id, out$.row_id), , drop = FALSE]
}

.hee_connectivity_pairs <- function(distances, landscape, times, distance_col,
                                    from_col, to_col, time_col) {
  if (!is.null(distances)) {
    if (is.matrix(distances)) {
      regs <- rownames(distances)
      if (is.null(regs)) regs <- as.character(seq_len(nrow(distances)))
      if (is.null(times)) times <- 0
      base <- expand.grid(from_region = regs, to_region = colnames(distances) %||% regs,
                          time_ma = times, stringsAsFactors = FALSE)
      base[[distance_col]] <- as.numeric(distances[cbind(match(base$from_region, regs),
                                                         match(base$to_region, colnames(distances) %||% regs))])
      names(base)[names(base) == "from_region"] <- from_col
      names(base)[names(base) == "to_region"] <- to_col
      names(base)[names(base) == "time_ma"] <- time_col
      return(base)
    }
    d <- as.data.frame(distances)
    if (!distance_col %in% names(d)) {
      alt <- intersect(c("distance", "distance_km", "paleodistance_km"), names(d))
      if (length(alt) == 0L) {
        stop("distances must contain a distance column.", call. = FALSE)
      }
      names(d)[match(alt[1], names(d))] <- distance_col
    }
    .require_cols(d, c(from_col, to_col, time_col, distance_col), "distances")
    .hee_validate_time_values(d[[time_col]], paste0("distances$", time_col))
    .hee_check_unique_keys(d, c(from_col, to_col, time_col), "distances")
    dist_raw <- suppressWarnings(as.numeric(d[[distance_col]]))
    if (any(!is.finite(dist_raw) | dist_raw < 0, na.rm = TRUE)) {
      stop("`", distance_col, "` in distances must contain finite non-negative values.",
           call. = FALSE)
    }
    d[[distance_col]] <- dist_raw
    return(d)
  }

  if (is.null(landscape)) {
    stop("Provide either distances or landscape with region/time coordinates.",
         call. = FALSE)
  }
  l <- as.data.frame(landscape)
  .require_cols(l, c("region", time_col, "lon", "lat"), "landscape")
  if (is.null(times)) times <- sort(unique(l[[time_col]]), decreasing = TRUE)
  out <- list()
  n <- 0L
  for (tt in times) {
    z <- l[l[[time_col]] == tt, , drop = FALSE]
    cent <- .hee_region_centroids(z, "region", "lon", "lat")
    regs <- as.character(cent$region)
    if (length(regs) == 0L) next
    pair <- expand.grid(from_region = regs, to_region = regs,
                        stringsAsFactors = FALSE)
    a <- cent[match(pair$from_region, cent$region), ]
    b <- cent[match(pair$to_region, cent$region), ]
    pair[[time_col]] <- tt
    pair[[distance_col]] <- hee_great_circle_distance_km(
      a$lon, a$lat, b$lon, b$lat
    )
    names(pair)[names(pair) == "from_region"] <- from_col
    names(pair)[names(pair) == "to_region"] <- to_col
    n <- n + 1L
    out[[n]] <- pair
  }
  do.call(rbind, out)
}

.hee_event_sum <- function(events, region, types, probability_col, region_col,
                           from_col, to_col, direction = c("in", "out")) {
  direction <- match.arg(direction)
  if (nrow(events) == 0L) return(0)
  ev_type <- tolower(as.character(events$event_type))
  keep <- ev_type %in% tolower(types)
  if (region_col %in% names(events)) keep <- keep & as.character(events[[region_col]]) == region
  if (to_col %in% names(events) && direction == "in") {
    keep <- keep | (ev_type %in% tolower(types) & as.character(events[[to_col]]) == region)
  }
  if (from_col %in% names(events) && direction == "out") {
    keep <- keep | (ev_type %in% tolower(types) & as.character(events[[from_col]]) == region)
  }
  sum(suppressWarnings(as.numeric(events[[probability_col]][keep])), na.rm = TRUE)
}

.hee_join_optional <- function(target, table, keys, value_cols) {
  .hee_safe_left_join(target, table, keys, value_cols,
                      target_name = "target", source_name = "optional table")
}

.hee_previous_probability <- function(target, previous_state, initial_probability,
                                      species_col, cell_col, time_col,
                                      probability_col, time_direction) {
  if (is.null(previous_state)) return(rep(initial_probability, nrow(target)))
  prev <- as.data.frame(previous_state)
  .require_cols(prev, c(species_col, cell_col, time_col, probability_col),
                "previous_state")
  out <- rep(initial_probability, nrow(target))
  prev_key <- interaction(prev[[species_col]], prev[[cell_col]], drop = TRUE)
  target_key <- interaction(target[[species_col]], target[[cell_col]], drop = TRUE)
  split_target <- split(seq_len(nrow(target)), target_key)
  for (key in names(split_target)) {
    idx <- split_target[[key]]
    src <- prev[prev_key == key, , drop = FALSE]
    if (nrow(src) == 0L) next
    for (i in idx) {
      tt <- target[[time_col]][i]
      candidates <- if (time_direction == "forward") {
        src[src[[time_col]] > tt, , drop = FALSE]
      } else {
        src[src[[time_col]] < tt, , drop = FALSE]
      }
      if (nrow(candidates) == 0L) {
        same <- src[src[[time_col]] == tt, , drop = FALSE]
        if (nrow(same) > 0L) candidates <- same
      }
      if (nrow(candidates) == 0L) next
      j <- which.min(abs(candidates[[time_col]] - tt))
      out[i] <- candidates[[probability_col]][j]
    }
  }
  out
}

.hee_col_or_default <- function(x, col, default) {
  if (is.null(col) || !col %in% names(x)) return(rep(default, nrow(x)))
  out <- suppressWarnings(as.numeric(x[[col]]))
  out[is.na(out)] <- default
  out
}

.hee_logit_bounded <- function(x, eps = 1e-9) {
  x <- .hee_clip01(x)
  x <- pmin(pmax(x, eps), 1 - eps)
  stats::qlogis(x)
}

.hee_metric_value <- function(target, source, value_col, default = 1,
                              keys = c("species", "cell_id", "region", "time_ma"),
                              aliases = character(),
                              source_name = "metric table") {
  if (is.null(source)) return(rep(default, nrow(target)))
  target_local <- as.data.frame(target)
  if (!".row_id" %in% names(target_local)) {
    target_local$.row_id <- seq_len(nrow(target_local))
  }
  if (!is.data.frame(source)) {
    val <- suppressWarnings(as.numeric(source))
    if (length(val) == 1L) val <- rep(val, nrow(target_local))
    if (length(val) != nrow(target_local)) {
      stop("Vector input for `", value_col, "` must have length 1 or match target rows.",
           call. = FALSE)
    }
    val[is.na(val)] <- default
    return(val)
  }
  src <- as.data.frame(source)
  candidates <- unique(c(value_col, aliases))
  hit <- intersect(candidates, names(src))
  if (length(hit) == 0L) {
    stop(source_name, " must contain column `", value_col, "`",
         if (length(aliases) > 0L) {
           paste0(" or one of: ", paste(aliases, collapse = ", "))
         } else {
           ""
         },
         ".", call. = FALSE)
  }
  if (!value_col %in% names(src)) src[[value_col]] <- src[[hit[1]]]
  join_keys <- intersect(keys, intersect(names(target_local), names(src)))
  if (length(join_keys) == 0L) {
    val <- src[[value_col]]
    if (length(val) == 1L) return(rep(val, nrow(target_local)))
    stop(source_name, " has multiple rows but no shared join keys with the target. ",
         "Supply shared keys (for example species, cell_id, region, time_ma) ",
         "or pass a scalar/vector instead.", call. = FALSE)
  }
  tmp <- target_local[, unique(c(".row_id", join_keys)), drop = FALSE]
  src <- src[, unique(c(join_keys, value_col)), drop = FALSE]
  if (any(duplicated(src[, join_keys, drop = FALSE]))) {
    ex <- .hee_format_duplicate_keys(src, join_keys)
    stop("Duplicate key(s) in metric table for ", value_col, " by ",
         paste(join_keys, collapse = ", "), ": ", paste(ex, collapse = " | "),
         call. = FALSE)
  }
  merged <- merge(tmp, src, by = join_keys, all.x = TRUE, sort = FALSE)
  merged <- merged[match(target_local$.row_id, merged$.row_id), , drop = FALSE]
  val <- suppressWarnings(as.numeric(merged[[value_col]]))
  val[is.na(val)] <- default
  val
}

.hee_first_table <- function(...) {
  xs <- list(...)
  for (x in xs) {
    if (is.data.frame(x)) return(x)
  }
  NULL
}

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) {
    if (is.null(x)) y else x
  }
}
