# Paleoenvironment adapter --------------------------------------------------

#' Detect roles of palaeoenvironmental variables
#'
#' Builds a role table that maps variables in a palaeoenvironmental time cube to
#' HmscEcoEvo process components. The defaults are tuned for the Phanerozoic
#' Environment 540-0 Ma land-harmonized products, but the function also works
#' with smaller cubes when variables use common names such as `bio1`, `bio12`,
#' `elev`, `land_mask_dem`, and `land_area_km2`.
#'
#' A variable may appear in more than one role. For example, temperature can be
#' an HMSC predictor, an extrapolation variable, and a climatic-connectivity
#' driver. The table records how the package should use a layer; it does not
#' claim that the layer is a direct causal estimate.
#'
#' @param env_cube Optional time cube or path accepted by [hee_load_timecube()].
#' @param variables Optional character vector of variable names.
#' @param metadata Optional metadata table with at least `variable`; when
#'   `env_cube` is supplied, cube metadata are used automatically.
#' @param extra_roles Optional data.frame with columns matching the output
#'   columns. User rows are appended and can document project-specific layers.
#'
#' @return A data.frame with variable-role mappings, availability, units,
#'   target framework component, recommended transform, and interpretation
#'   boundary.
#' @export
#'
#' @examples
#' hee_paleoenv_variable_roles(variables = c("bio1", "bio12", "land_mask_dem"))
hee_paleoenv_variable_roles <- function(env_cube = NULL,
                                        variables = NULL,
                                        metadata = NULL,
                                        extra_roles = NULL) {
  if (!is.null(env_cube)) {
    cube <- .hee_as_timecube(env_cube)
    variables <- variables %||% cube$variables
    if (is.null(metadata) && is.data.frame(cube$metadata)) {
      metadata <- cube$metadata
    }
  }
  variables <- unique(as.character(variables %||% character()))
  if (length(variables) == 0L) {
    stop("Provide `env_cube` or a non-empty `variables` vector.", call. = FALSE)
  }

  role_defs <- .hee_paleoenv_role_definitions()
  rows <- lapply(seq_len(nrow(role_defs)), function(i) {
    def <- role_defs[i, , drop = FALSE]
    matched <- .hee_match_role_variables(variables, def$exact[[1]], def$pattern)
    if (length(matched) == 0L) return(NULL)
    data.frame(
      variable = matched,
      role = def$role,
      target_component = def$target_component,
      target_function = def$target_function,
      recommended_use = def$recommended_use,
      default_transform = def$default_transform,
      interpretation = def$interpretation,
      estimate_type = def$estimate_type,
      priority = def$priority,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    out <- data.frame(
      variable = variables,
      role = "unmapped",
      target_component = NA_character_,
      target_function = NA_character_,
      recommended_use = "available but not automatically mapped",
      default_transform = NA_character_,
      interpretation = "Inspect manually before using in a model.",
      estimate_type = "unclassified",
      priority = 99L,
      stringsAsFactors = FALSE
    )
  }

  missing_vars <- setdiff(variables, unique(out$variable))
  if (length(missing_vars) > 0L) {
    out <- rbind(out, data.frame(
      variable = missing_vars,
      role = "unmapped",
      target_component = NA_character_,
      target_function = NA_character_,
      recommended_use = "available but not automatically mapped",
      default_transform = NA_character_,
      interpretation = "Inspect manually before using in a model.",
      estimate_type = "unclassified",
      priority = 99L,
      stringsAsFactors = FALSE
    ))
  }

  if (!is.null(extra_roles)) {
    er <- as.data.frame(extra_roles, stringsAsFactors = FALSE)
    for (nm in setdiff(names(out), names(er))) er[[nm]] <- NA
    out <- rbind(out, er[, names(out), drop = FALSE])
  }

  out$available <- out$variable %in% variables
  if (is.data.frame(metadata) && nrow(metadata) > 0L &&
      "variable" %in% names(metadata)) {
    keep <- intersect(c("variable", "long_name", "unit", "units", "source",
                        "method", "aggregation_method"), names(metadata))
    meta <- metadata[, keep, drop = FALSE]
    meta <- meta[!duplicated(meta$variable), , drop = FALSE]
    out <- merge(out, meta, by = "variable", all.x = TRUE, sort = FALSE)
    if ("units" %in% names(out) && !"unit" %in% names(out)) {
      names(out)[names(out) == "units"] <- "unit"
    }
  }
  out <- out[order(out$priority, out$variable, out$role), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Select default palaeoenvironment variables for process modules
#'
#' @param env_cube Optional time cube or path accepted by [hee_load_timecube()].
#' @param variables Optional character vector of variable names.
#' @param metadata Optional metadata table.
#'
#' @return A data.frame with one preferred variable per process role, plus
#'   alternatives and status.
#' @export
#'
#' @examples
#' hee_paleoenv_default_variable_map(
#'   variables = c("MAT_pohl_C", "MAP_pohl_mm_yr", "land_mask_dem"))
hee_paleoenv_default_variable_map <- function(env_cube = NULL,
                                              variables = NULL,
                                              metadata = NULL) {
  roles <- hee_paleoenv_variable_roles(env_cube = env_cube,
                                       variables = variables,
                                       metadata = metadata)
  preferred_roles <- c(
    "land_mask", "land_area", "elevation", "slope", "relief",
    "coastal_distance", "temperature", "precipitation", "moisture",
    "wetland", "bryophyte_habitat", "temperature_uncertainty",
    "precipitation_uncertainty", "climate_class", "hydro_class"
  )
  rows <- lapply(preferred_roles, function(role) {
    z <- roles[roles$role == role & roles$available, , drop = FALSE]
    pref <- .hee_paleoenv_preferred_variables(role)
    z$.preferred_rank <- match(z$variable, pref)
    z$.preferred_rank[is.na(z$.preferred_rank)] <- Inf
    z <- z[order(z$priority, z$.preferred_rank, z$variable), , drop = FALSE]
    if (nrow(z) == 0L) {
      return(data.frame(
        role = role, variable = NA_character_, status = "missing",
        alternatives = NA_character_, target_component = NA_character_,
        target_function = NA_character_, stringsAsFactors = FALSE
      ))
    }
    data.frame(
      role = role,
      variable = z$variable[1],
      status = "available",
      alternatives = paste(unique(z$variable), collapse = ";"),
      target_component = z$target_component[1],
      target_function = z$target_function[1],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Prepare palaeoenvironment-derived inputs for geoprocess functions
#'
#' Extracts a full palaeoenvironmental grid and derives standard input tables
#' for the Dynamic Earth-Biota Assembly functions. This adapter is intentionally
#' conservative: by default, `G_arena` equals land availability only, so climate
#' and moisture layers are not double-counted when HMSC also uses them as
#' environmental predictors. If a study wants a habitat-gated scenario, set
#' `habitat_strategy` explicitly.
#'
#' @param env_cube Time cube or path accepted by [hee_load_timecube()].
#' @param times Ages in Ma. Defaults to the full cube axis.
#' @param hmsc_variables Optional HMSC predictor variables to include in the
#'   extracted grid.
#' @param role_map Optional output from `hee_paleoenv_default_variable_map()`.
#'   When omitted it is detected.
#' @param habitat_strategy One of `"land_only"`, `"moisture"`, `"wetland"`,
#'   `"bryophyte"`, or `"none"`. The default avoids double-counting habitat
#'   proxies already used by HMSC.
#' @param land_only Logical. If `FALSE` the full grid is retained and ocean or
#'   unknown geography receives `L_land = 0`.
#' @param include_opportunity Logical; compute land age, landscape events, and
#'   ecological opportunity when possible.
#'
#' @return A list with `role_table`, `variable_map`, `paleo_grid`, `arena`,
#'   `landscape_state`, and optional `land_age`, `landscape_events`, and
#'   `ecological_opportunity` tables.
#' @export
#'
#' @examples
#' cube <- hee_load_timecube(object = local({
#'   e <- new.env(parent = emptyenv())
#'   e$age_ma <- c(10, 0)
#'   e$lon <- c(0, 1)
#'   e$lat <- c(0, 1)
#'   e$variable_names <- c("bio1", "bio12", "land_mask_dem", "land_area_km2")
#'   grid <- expand.grid(lon = e$lon, lat = e$lat)
#'   e$get_slice <- function(age, variables = NULL, as_data_table = TRUE) {
#'     d <- data.frame(time_ma = age, lon = grid$lon, lat = grid$lat)
#'     d$bio1 <- 10; d$bio12 <- 100; d$land_mask_dem <- 1; d$land_area_km2 <- 1
#'     d[, c("time_ma", "lon", "lat", variables), drop = FALSE]
#'   }
#'   e$get_metadata <- function(variable = NULL) data.frame(variable = e$variable_names)
#'   e
#' }))
#' hee_prepare_paleoenv_geoprocess_inputs(cube, times = 0)
hee_prepare_paleoenv_geoprocess_inputs <- function(env_cube,
                                                   times = NULL,
                                                   hmsc_variables = NULL,
                                                   role_map = NULL,
                                                   habitat_strategy = c("land_only", "moisture", "wetland", "bryophyte", "none"),
                                                   land_only = FALSE,
                                                   include_opportunity = TRUE) {
  habitat_strategy <- match.arg(habitat_strategy)
  cube <- .hee_as_timecube(env_cube)
  times <- times %||% cube$times
  times <- .hee_match_times(times, cube$times)
  role_table <- hee_paleoenv_variable_roles(cube)
  variable_map <- role_map %||% hee_paleoenv_default_variable_map(cube)

  pick <- function(role) {
    z <- variable_map[variable_map$role == role &
                        variable_map$status == "available", , drop = FALSE]
    if (nrow(z) == 0L) return(NA_character_)
    z$variable[1]
  }
  land_var <- pick("land_mask")
  if (is.na(land_var)) {
    stop("env_cube must contain a land/geographic availability layer such as ",
         "`land_mask_dem`, `LANDFRAC_li_fraction`, or `land`.", call. = FALSE)
  }
  area_var <- pick("land_area")
  elevation_var <- pick("elevation")
  slope_var <- pick("slope")
  relief_var <- pick("relief")
  coast_var <- pick("coastal_distance")
  temp_var <- pick("temperature")
  precip_var <- pick("precipitation")
  moisture_var <- pick("moisture")
  wetland_var <- pick("wetland")
  bryophyte_var <- pick("bryophyte_habitat")
  temp_uncert_var <- pick("temperature_uncertainty")
  precip_uncert_var <- pick("precipitation_uncertainty")

  habitat_var <- switch(
    habitat_strategy,
    land_only = NA_character_,
    none = NA_character_,
    moisture = moisture_var,
    wetland = wetland_var,
    bryophyte = bryophyte_var
  )
  if (habitat_strategy %in% c("moisture", "wetland", "bryophyte") &&
      is.na(habitat_var)) {
    stop("Requested habitat_strategy = `", habitat_strategy,
         "`, but the matching palaeoenvironment layer is not available.",
         call. = FALSE)
  }

  selected <- unique(stats::na.omit(c(
    hmsc_variables, land_var, area_var, elevation_var, slope_var, relief_var,
    coast_var, temp_var, precip_var, moisture_var, wetland_var, bryophyte_var,
    temp_uncert_var, precip_uncert_var, habitat_var
  )))
  selected <- intersect(selected, cube$variables)
  grid <- hee_make_paleo_grid(cube, times = times, variables = selected,
                              landmask = land_var, land_only = land_only)

  grid$L_land <- .hee_clip01(grid[[land_var]])
  grid$land_weight <- grid$L_land
  grid$geographic_existence <- ifelse(is.finite(grid$L_land) & grid$L_land > 0, 1, 0)
  if (!is.na(area_var) && area_var %in% names(grid)) {
    grid$cell_area_km2 <- suppressWarnings(as.numeric(grid[[area_var]]))
  }
  if (!is.na(elevation_var) && elevation_var %in% names(grid)) {
    grid$elev_m <- suppressWarnings(as.numeric(grid[[elevation_var]]))
  }
  if (!is.na(temp_var) && temp_var %in% names(grid)) {
    grid$bio1 <- suppressWarnings(as.numeric(grid[[temp_var]]))
  }
  if (!is.na(precip_var) && precip_var %in% names(grid)) {
    grid$bio12 <- suppressWarnings(as.numeric(grid[[precip_var]]))
  }
  if (!is.na(moisture_var) && moisture_var %in% names(grid)) {
    grid$moisture_index <- suppressWarnings(as.numeric(grid[[moisture_var]]))
  }
  if (!is.na(habitat_var) && habitat_var %in% names(grid)) {
    grid$habitat_availability <- .hee_z_to_unit_interval(grid[[habitat_var]])
  } else {
    grid$habitat_availability <- ifelse(grid$geographic_existence > 0, 1, 0)
  }

  arena <- hee_geographic_stage(
    grid,
    land_col = "L_land",
    habitat_col = if (habitat_strategy %in% c("moisture", "wetland", "bryophyte")) {
      "habitat_availability"
    } else {
      NULL
    },
    output_col = "G_arena"
  )
  arena$L_land <- grid$L_land
  arena$geographic_existence <- grid$geographic_existence
  arena_cols <- intersect(c("cell_id", "time_ma", "lon", "lat", "L_land",
                            "G_arena", "geographic_existence",
                            "habitat_availability", "land_area_km2",
                            "cell_area_km2"), names(arena))
  arena <- arena[, arena_cols, drop = FALSE]

  landscape_cols <- intersect(c(
    "cell_id", "time_ma", "lon", "lat", "L_land", "G_arena",
    "geographic_existence", "habitat_availability", "cell_area_km2",
    "land_area_km2", "elev_m", "bio1", "bio12", "moisture_index",
    land_var, area_var, elevation_var, slope_var, relief_var, coast_var,
    temp_var, precip_var, moisture_var, wetland_var, bryophyte_var,
    temp_uncert_var, precip_uncert_var
  ), names(grid))
  landscape <- grid[, unique(landscape_cols), drop = FALSE]

  out <- list(
    role_table = role_table,
    variable_map = variable_map,
    paleo_grid = grid,
    arena = arena,
    landscape_state = landscape,
    hmsc_variables = intersect(hmsc_variables %||% character(), names(grid)),
    habitat_strategy = habitat_strategy,
    land_variable = land_var
  )
  if (isTRUE(include_opportunity)) {
    aged <- hee_land_age(landscape, unit_col = "cell_id",
                         land_col = "geographic_existence",
                         habitat_col = NULL)
    out$land_age <- aged
    event_area <- if ("cell_area_km2" %in% names(aged)) "cell_area_km2" else NULL
    event_elev <- if ("elev_m" %in% names(aged)) "elev_m" else NULL
    out$landscape_events <- hee_landscape_events(
      aged,
      unit_col = "cell_id",
      land_col = "geographic_existence",
      habitat_col = "habitat_availability",
      area_col = event_area,
      elevation_col = event_elev
    )
    env_cols <- intersect(c("bio1", "bio12", "moisture_index", slope_var,
                            relief_var, coast_var), names(aged))
    out$ecological_opportunity <- hee_ecological_opportunity(
      aged,
      env_cols = env_cols,
      region_col = "cell_id",
      area_col = if ("cell_area_km2" %in% names(aged)) "cell_area_km2" else "cell_area_km2"
    )
  }
  class(out) <- c("hee_paleoenv_geoprocess_inputs", "list")
  out
}

.hee_paleoenv_role_definitions <- function() {
  data.frame(
    role = c(
      "land_mask", "land_mask", "land_area", "elevation", "elevation",
      "bathymetry", "slope", "relief", "relief", "coastal_distance",
      "shelf", "deep_ocean", "temperature", "temperature", "temperature",
      "temperature_uncertainty", "temperature_uncertainty", "precipitation",
      "precipitation", "precipitation", "precipitation_uncertainty",
      "precipitation_uncertainty", "water_balance", "water_balance",
      "aridity", "thermal_months", "wet_months", "dry_season",
      "moisture", "wetland", "bryophyte_habitat", "climate_class",
      "climate_class", "hydro_class", "wind", "wind", "albedo"
    ),
    exact = I(rep(list(character()), 37L)),
    pattern = c(
      "^land_mask|^land$|^land_weight$|^L_land$|^G_arena$", "LANDFRAC", "land.*area",
      "^elevation_m$|^elev_m$|^elev$", "paleodem|PHIS", "bathym",
      "slope", "relief", "tpi", "coast", "shelf", "deep.*ocean",
      "^MAT_pohl_C$|^MAT_gridded_model_mean_C$|^bio1$", "^MAT_li_C$", "^T_|seasonality",
      "^MAT_gridded_model_sd_C$", "delta_MAT", "^MAP_pohl_mm_yr$|^MAP_gridded_model_mean_mm_yr$|^bio12$",
      "^MAP_li_mm_yr$", "^P_.*month|P_seasonality", "^MAP_gridded_model_sd_mm_yr$",
      "delta_MAP|ratio_MAP", "P_minus_E", "runoff", "aridity",
      "months_above", "months_wet", "dry_season", "^moisture_|moisture_availability",
      "wetland", "bryophyte", "koppen_recomputed", "koppen_pohl",
      "hydro", "U1000", "V1000", "SALB"
    ),
    target_component = c(
      "L_land/G_arena", "L_land/G_arena", "ecological_opportunity",
      "S_HMSC/topographic_resistance", "S_HMSC/topographic_resistance",
      "marine_boundary", "topographic_resistance", "topographic_heterogeneity",
      "topographic_heterogeneity", "coastal_connectivity",
      "marine_or_coastal_mask", "marine_mask", "S_HMSC/climatic_connectivity",
      "S_HMSC/climatic_connectivity", "S_HMSC/climate_extremes",
      "Q_extrap/climate_uncertainty", "Q_extrap/climate_uncertainty",
      "S_HMSC/climatic_connectivity", "S_HMSC/climatic_connectivity",
      "S_HMSC/climate_extremes", "Q_extrap/climate_uncertainty",
      "Q_extrap/climate_uncertainty", "habitat_water_balance",
      "habitat_water_balance", "habitat_dryness", "thermal_window",
      "wetness_window", "disturbance_or_drought_proxy", "habitat_availability",
      "habitat_availability", "habitat_availability", "habitat_class",
      "habitat_class", "habitat_class", "climate_dynamics", "climate_dynamics",
      "climate_surface"
    ),
    target_function = c(
      "hee_geographic_stage; hee_combine", "hee_geographic_stage; hee_combine",
      "hee_ecological_opportunity", "hee_project_hmsc_timecube; hee_build_connectivity_cube",
      "hee_project_hmsc_timecube; hee_build_connectivity_cube",
      "map/report only for terrestrial workflow", "hee_build_connectivity_cube",
      "hee_ecological_opportunity; hee_build_connectivity_cube",
      "hee_ecological_opportunity; hee_build_connectivity_cube",
      "hee_build_connectivity_cube", "optional marine/coastal workflow",
      "optional marine workflow", "hee_project_hmsc_timecube; hee_extrapolation_risk",
      "hee_project_hmsc_timecube; hee_extrapolation_risk",
      "hee_project_hmsc_timecube; hee_extrapolation_risk",
      "hee_prediction_reliability", "hee_prediction_reliability",
      "hee_project_hmsc_timecube; hee_extrapolation_risk",
      "hee_project_hmsc_timecube; hee_extrapolation_risk",
      "hee_project_hmsc_timecube; hee_extrapolation_risk",
      "hee_prediction_reliability", "hee_prediction_reliability",
      "hee_ecological_opportunity", "hee_ecological_opportunity",
      "hee_ecological_opportunity", "hee_ecological_opportunity",
      "hee_ecological_opportunity", "hee_ecological_opportunity",
      "hee_geographic_stage optional habitat gate",
      "hee_geographic_stage optional habitat gate",
      "hee_geographic_stage optional habitat gate",
      "reporting/stratified validation", "reporting/stratified validation",
      "reporting/stratified validation", "optional dispersal/climate dynamics",
      "optional dispersal/climate dynamics", "optional climate diagnostics"
    ),
    recommended_use = c(
      "Hard land/geography multiplier; fraction is allowed.",
      "Secondary land fraction if DEM land mask is absent.",
      "Cell area, area gain and land-stage size.",
      "HMSC predictor and topographic resistance; do not call absolute elevation slope.",
      "Alternative elevation source; inspect before mixing with elevation_m.",
      "Ocean context only unless running marine analyses.",
      "Terrain resistance and heterogeneity.",
      "Terrain heterogeneity and resistance.",
      "Terrain position and heterogeneity.",
      "Coastal access or isolation proxy.",
      "Optional shallow-marine/coastal context.",
      "Ocean mask for terrestrial exclusion.",
      "Primary temperature predictor.",
      "Alternative model temperature; useful for sensitivity.",
      "Temperature extremes and seasonality.",
      "Climate-model uncertainty diagnostic.",
      "Climate-model disagreement diagnostic.",
      "Primary precipitation predictor.",
      "Alternative model precipitation; useful for sensitivity.",
      "Precipitation extremes and seasonality.",
      "Climate-model uncertainty diagnostic.",
      "Climate-model disagreement diagnostic.",
      "Water-balance habitat/opportunity proxy.",
      "Hydrological opportunity proxy.",
      "Dryness/aridity proxy.",
      "Thermal physiological-window proxy.",
      "Wet-season availability proxy.",
      "Drought/seasonality proxy.",
      "Generic moisture opportunity; optional habitat gate.",
      "Wetland opportunity; optional habitat gate.",
      "Bryophyte-specific habitat proxy; optional habitat gate.",
      "Categorical climate class for summaries.",
      "Source Koppen class for sensitivity summaries.",
      "Hydrogeomorphic class for summaries.",
      "Optional atmospheric-flow context.",
      "Optional atmospheric-flow context.",
      "Optional surface-energy context."
    ),
    default_transform = c(
      "clip01; supplied NA -> 0", "clip01; supplied NA -> 0",
      "non-negative area", "raw or recipe scaling", "raw or recipe scaling",
      "report separately", "non-negative; optionally rescale", "positive index",
      "signed or absolute index", "non-negative km", "clip01", "clip01",
      "0 Ma recipe scaling", "0 Ma recipe scaling", "0 Ma recipe scaling",
      "positive reliability-risk component", "absolute difference or positive risk",
      "0 Ma recipe scaling", "0 Ma recipe scaling", "0 Ma recipe scaling",
      "positive reliability-risk component", "positive reliability-risk component",
      "recipe or opportunity proxy", "positive hydrological index",
      "positive or inverse dryness index", "positive month count", "positive month count",
      "positive drought index", "logistic z-to-[0,1] if used as habitat",
      "logistic z-to-[0,1] if used as habitat",
      "logistic z-to-[0,1] if used as habitat", "categorical",
      "categorical", "categorical", "optional", "optional", "optional"
    ),
    interpretation = c(
      rep("External palaeoenvironment-derived layer; not inferred by HmscEcoEvo.", 2),
      "Arena size and area change, not habitat quality by itself.",
      "Environmental predictor or resistance proxy depending on use.",
      "Alternative elevation source; mixing sources should be documented.",
      "Marine context; not terrestrial land existence.",
      "Topographic resistance/heterogeneity proxy.",
      "Topographic heterogeneity proxy.",
      "Topographic position proxy.",
      "Coastal isolation/access proxy.",
      "Marine/coastal context.",
      "Marine/ocean exclusion context.",
      rep("Climate predictor; suitability comes from fitted HMSC, not this layer alone.", 3),
      rep("Climate uncertainty diagnostic, not biological uncertainty.", 2),
      rep("Climate predictor; suitability comes from fitted HMSC, not this layer alone.", 3),
      rep("Climate uncertainty diagnostic, not biological uncertainty.", 2),
      rep("Hydrological proxy; use cautiously and avoid double-counting with HMSC predictors.", 9),
      rep("Categorical summary/stratification layer, not continuous suitability.", 3),
      rep("Optional climate-dynamics context.", 3)
    ),
    estimate_type = c(
      rep("external_geography", 12),
      rep("external_climate_predictor", 3),
      rep("external_uncertainty_proxy", 2),
      rep("external_climate_predictor", 3),
      rep("external_uncertainty_proxy", 2),
      rep("proxy_or_predictor", 9),
      rep("categorical_context", 3),
      rep("optional_context", 3)
    ),
    priority = c(
      1, 3, 1, 1, 4, 5, 1, 1, 2, 1, 4, 5,
      1, 3, 4, 1, 2, 1, 3, 4, 1, 2,
      2, 2, 2, 3, 3, 3, 1, 1, 1, 1, 2, 1, 4, 4, 4
    ),
    stringsAsFactors = FALSE
  )
}

.hee_match_role_variables <- function(variables, exact, pattern) {
  exact <- exact %||% character()
  hit_exact <- intersect(exact, variables)
  hit_pattern <- if (!is.na(pattern) && nzchar(pattern)) {
    variables[grepl(pattern, variables, ignore.case = TRUE)]
  } else {
    character()
  }
  unique(c(hit_exact, hit_pattern))
}

.hee_paleoenv_preferred_variables <- function(role) {
  switch(
    role,
    land_mask = c("land_mask_dem", "land", "land_weight", "L_land",
                  "G_arena", "LANDFRAC_li_fraction"),
    land_area = c("land_area_km2", "cell_area_km2", "area_km2"),
    elevation = c("elevation_m", "elev_m", "elev", "paleodem_m",
                  "PHIS_li_m"),
    temperature = c("MAT_pohl_C", "MAT_gridded_model_mean_C", "bio1",
                    "MAT_li_C"),
    precipitation = c("MAP_pohl_mm_yr", "MAP_gridded_model_mean_mm_yr",
                      "bio12", "MAP_li_mm_yr"),
    moisture = c("moisture_availability_index_z", "P_minus_E_pohl_mm_yr",
                 "P_minus_E_recomputed_mm_yr"),
    wetland = c("wetland_potential_index_z"),
    bryophyte_habitat = c("bryophyte_moisture_score_z"),
    temperature_uncertainty = c("MAT_gridded_model_sd_C",
                                "delta_MAT_li_minus_pohl_C"),
    precipitation_uncertainty = c("MAP_gridded_model_sd_mm_yr",
                                  "delta_MAP_li_minus_pohl_mm_yr",
                                  "ratio_MAP_li_over_pohl"),
    climate_class = c("koppen_recomputed_major", "koppen_pohl_source_nearest"),
    hydro_class = c("hydro_class_simple"),
    character()
  )
}

.hee_z_to_unit_interval <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep(0, length(x))
  ok <- is.finite(x)
  out[ok] <- stats::plogis(x[ok])
  .hee_clip01(out)
}
