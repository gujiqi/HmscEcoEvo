# Case 03 real-data infrastructure -----------------------------------------

#' Validate Case 03 real-data input paths
#'
#' Separates inputs required for the core modern-HMSC/M1 workflow from inputs
#' required only by phylogenetic, historical-biogeographic, plate-track, or
#' validation modules. Missing optional inputs never trigger simulated
#' replacements.
#'
#' @param paths Named list of input paths. Recognised names are `comm`, `sites`,
#'   `traits`, `tree`, `species_origin`, `tip_ranges`, `env_rds`, `region_cube`,
#'   `landmask_cube`, `plate_points`, `bgb_accessibility`, `connectivity_paths`,
#'   `bsm_events`, and `fossils`.
#' @return A list with `manifest`, `missing_core`, `missing_optional`, and `ok`.
#' @export
hee_case03_validate_real_inputs <- function(paths) {
  paths <- as.list(paths)
  spec <- data.frame(
    input = c("comm", "sites", "env_rds", "traits", "tree",
              "species_origin", "tip_ranges", "region_cube", "landmask_cube",
              "plate_points", "bgb_accessibility", "connectivity_paths",
              "bsm_events", "fossils"),
    argument = c("--comm", "--sites", "--env_rds", "--traits", "--tree",
                 "--species_origin", "--tip_ranges", "--region_cube",
                 "--landmask_cube", "--plate_points", "--bgb_accessibility",
                 "--connectivity_paths", "--bsm_events", "--fossils"),
    requirement_level = c(
      rep("required_core", 3),
      "optional_hmsc_traits", "required_M2_M5", "required_M2_M5",
      "required_BioGeoBEARS", "required_regional_M3_M5",
      "optional_if_env_has_landmask", "optional_plate_track",
      "required_M3_M5", "optional_external_least_cost_paths",
      "optional_BSM_events", "optional_validation"
    ),
    affected_modules = c(
      "HMSC,M1-M5", "HMSC,M1-M5", "HMSC,M1-M5",
      "HMSC trait mediation", "M2-M5", "M2-M5",
      "BioGeoBEARS ancestral ranges", "regional connectivity,M3-M5",
      "land masking", "plate-track validation", "M3-M5",
      "M4-M5 connectivity", "BSM event diagnostics", "fossil validation"
    ),
    stringsAsFactors = FALSE
  )
  supplied <- vapply(spec$input, function(nm) {
    value <- paths[[nm]]
    !is.null(value) && length(value) == 1L && !is.na(value) && nzchar(value)
  }, logical(1))
  path_value <- vapply(spec$input, function(nm) {
    value <- paths[[nm]]
    if (is.null(value) || length(value) == 0L || is.na(value)) "" else as.character(value[[1]])
  }, character(1))
  exists <- supplied & file.exists(path_value)
  spec$path <- path_value
  spec$status <- ifelse(exists, "AVAILABLE", "MISSING")
  spec$required_for_core <- spec$requirement_level == "required_core"
  missing <- spec[spec$status == "MISSING", , drop = FALSE]
  list(
    manifest = spec,
    missing_core = missing[missing$required_for_core, , drop = FALSE],
    missing_optional = missing[!missing$required_for_core, , drop = FALSE],
    ok = !any(missing$required_for_core)
  )
}

#' Great-circle distance in kilometres
#'
#' Uses the haversine formula with a mean Earth radius of 6371.0088 km. Input
#' longitudes are wrapped automatically, so 179 and -179 degrees are two
#' degrees apart rather than 358 degrees apart.
#'
#' @param lon1,lat1,lon2,lat2 Coordinates in decimal degrees. Arguments are
#'   recycled using base-R vector rules.
#' @return Numeric great-circle distances in kilometres.
#' @export
hee_great_circle_distance_km <- function(lon1, lat1, lon2, lat2) {
  vals <- list(lon1 = lon1, lat1 = lat1, lon2 = lon2, lat2 = lat2)
  vals <- lapply(vals, function(z) suppressWarnings(as.numeric(z)))
  n <- max(lengths(vals))
  vals <- lapply(vals, rep_len, length.out = n)
  bad <- Reduce(`|`, lapply(vals, function(z) !is.finite(z)))
  if (any(abs(vals$lat1[!bad]) > 90 | abs(vals$lat2[!bad]) > 90)) {
    stop("Latitude values must be between -90 and 90 degrees.", call. = FALSE)
  }
  rad <- pi / 180
  dlat <- (vals$lat2 - vals$lat1) * rad
  dlon_deg <- ((vals$lon2 - vals$lon1 + 180) %% 360) - 180
  dlon <- dlon_deg * rad
  a <- sin(dlat / 2)^2 + cos(vals$lat1 * rad) *
    cos(vals$lat2 * rad) * sin(dlon / 2)^2
  a <- pmin(pmax(a, 0), 1)
  out <- 2 * 6371.0088 * asin(sqrt(a))
  out[bad] <- NA_real_
  out
}

#' Match modern survey sites to the 0 Ma environmental grid
#'
#' @param sites Table with `site_id`, `lon`, and `lat`.
#' @param grid0 A 0 Ma grid with `cell_id`, `lon`, and `lat`.
#' @param max_distance_km Maximum acceptable nearest-cell distance.
#' @param action Whether excessive distances are flagged, warned about, or
#'   treated as an error.
#' @return Site rows plus matched grid columns, `matched_distance_km`, and
#'   `match_status`.
#' @export
hee_case03_match_sites_to_grid <- function(sites, grid0,
                                           max_distance_km = Inf,
                                           action = c("flag", "warn", "error")) {
  action <- match.arg(action)
  sites <- as.data.frame(sites)
  grid0 <- as.data.frame(grid0)
  .require_cols(sites, c("site_id", "lon", "lat"), "sites")
  .require_cols(grid0, c("cell_id", "lon", "lat"), "grid0")
  .hee_check_unique_keys(sites, "site_id", "sites")
  .hee_check_unique_keys(grid0, "cell_id", "grid0")
  if (nrow(grid0) == 0L) stop("grid0 has no candidate cells.", call. = FALSE)
  if (!is.finite(max_distance_km) && !is.infinite(max_distance_km)) {
    stop("max_distance_km must be non-negative or Inf.", call. = FALSE)
  }
  if (max_distance_km < 0) stop("max_distance_km must be non-negative.", call. = FALSE)
  matched <- lapply(seq_len(nrow(sites)), function(i) {
    d <- hee_great_circle_distance_km(
      sites$lon[[i]], sites$lat[[i]], grid0$lon, grid0$lat
    )
    if (all(!is.finite(d))) {
      stop("No finite grid distance is available for site `", sites$site_id[[i]], "`.",
           call. = FALSE)
    }
    j <- which.min(d)
    extra <- grid0[j, setdiff(names(grid0), c("lon", "lat")), drop = FALSE]
    cbind(sites[i, , drop = FALSE], extra,
          matched_grid_lon = grid0$lon[[j]],
          matched_grid_lat = grid0$lat[[j]],
          matched_distance_km = d[[j]],
          stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, matched)
  out$match_status <- ifelse(out$matched_distance_km <= max_distance_km,
                             "MATCHED", "DISTANCE_EXCEEDED")
  bad <- out$match_status == "DISTANCE_EXCEEDED"
  if (any(bad)) {
    msg <- paste0(sum(bad), " site(s) exceeded max_distance_km = ",
                  max_distance_km, ". Worst distance: ",
                  round(max(out$matched_distance_km[bad]), 3), " km.")
    if (action == "error") stop(msg, call. = FALSE)
    if (action == "warn") warning(msg, call. = FALSE)
  }
  rownames(out) <- NULL
  out
}

#' Draw spatially stratified demonstration sites from a modern grid
#'
#' Selects unique modern grid cells while retaining at least one site from every
#' stratum. This helper is intended only for simulated Case 03 demonstrations
#' and experimental-design sensitivity checks. Real-data workflows must read
#' observed survey coordinates from `sites.csv`; this function must not be used
#' to manufacture empirical sampling locations.
#'
#' @param grid Modern grid containing `cell_id`, `lon`, `lat`, and a stratum
#'   column such as `region`.
#' @param n_sites Number of unique cells to select. It must be at least the
#'   number of non-empty strata and no greater than the number of grid cells.
#' @param strata Name of the categorical spatial-stratum column. The default
#'   `"region"` gives broad geographical coverage.
#' @return A randomly ordered subset of `grid`, with `sampling_stratum` and a
#'   sequential `site_id`. Selection is reproducible after [set.seed()].
#' @details The function samples without replacement. It prevents the ordered
#'   grid-index sampling that can accidentally place demonstration sites along a
#'   visually misleading transect. It does not correct survey bias or infer
#'   real-world sampling design.
#' @export
hee_case03_stratified_site_sample <- function(grid, n_sites,
                                              strata = "region") {
  grid <- as.data.frame(grid)
  .require_cols(grid, c("cell_id", "lon", "lat", strata), "grid")
  .hee_check_unique_keys(grid, "cell_id", "grid")
  if (length(n_sites) != 1L || !is.finite(n_sites) ||
      n_sites < 1 || n_sites != as.integer(n_sites)) {
    stop("n_sites must be one positive integer.", call. = FALSE)
  }
  n_sites <- as.integer(n_sites)
  if (n_sites > nrow(grid)) {
    stop("n_sites cannot exceed the number of available modern grid cells.",
         call. = FALSE)
  }
  labels <- as.character(grid[[strata]])
  if (anyNA(labels) || any(!nzchar(labels))) {
    stop("grid spatial strata contain missing or empty labels.", call. = FALSE)
  }
  groups <- split(seq_len(nrow(grid)), labels, drop = TRUE)
  n_groups <- length(groups)
  if (n_sites < n_groups) {
    stop("n_sites must be at least the number of non-empty spatial strata (",
         n_groups, ") to retain geographical coverage.", call. = FALSE)
  }

  # Seed every stratum once, then allocate remaining cells by unused capacity.
  allocation <- rep.int(1L, n_groups)
  capacity <- lengths(groups)
  remaining <- n_sites - n_groups
  while (remaining > 0L) {
    eligible <- which(allocation < capacity)
    if (length(eligible) == 0L) break
    weight <- capacity[eligible] - allocation[eligible]
    # `sample(6, 1)` means sample from 1:6 rather than from the one-element
    # vector c(6); sample an index to avoid that base-R special case.
    chosen <- eligible[[sample.int(length(eligible), size = 1L, prob = weight)]]
    allocation[[chosen]] <- allocation[[chosen]] + 1L
    remaining <- remaining - 1L
  }
  selected <- unlist(Map(function(idx, n) sample(idx, size = n, replace = FALSE),
                          groups, allocation), use.names = FALSE)
  selected <- sample(selected, size = length(selected), replace = FALSE)
  out <- grid[selected, , drop = FALSE]
  out$sampling_stratum <- as.character(out[[strata]])
  out$site_id <- sprintf("site_%03d", seq_len(nrow(out)))
  rownames(out) <- NULL
  out
}

#' Map modern spatial sampling gaps
#'
#' Computes the great-circle distance from every 0 Ma prediction cell to the
#' nearest observed survey site. This diagnoses spatial coverage of the modern
#' HMSC training data; it is not a correction for sampling bias by itself.
#'
#' @param grid0 Modern grid with `cell_id`, `lon`, and `lat`.
#' @param sites Survey table with `site_id`, `lon`, and `lat`.
#' @param coverage_distance_km Distance threshold defining survey coverage.
#' @param chunk_size Number of grid cells processed together.
#' @return Grid rows plus nearest survey distance/site and binary coverage.
#' @export
hee_case03_sampling_gap <- function(grid0, sites,
                                    coverage_distance_km = 250,
                                    chunk_size = 5000L) {
  grid0 <- as.data.frame(grid0)
  sites <- as.data.frame(sites)
  .require_cols(grid0, c("cell_id", "lon", "lat"), "grid0")
  .require_cols(sites, c("site_id", "lon", "lat"), "sites")
  if (nrow(grid0) == 0L || nrow(sites) == 0L) {
    stop("grid0 and sites must both contain at least one row.", call. = FALSE)
  }
  if (!is.finite(coverage_distance_km) || coverage_distance_km < 0) {
    stop("coverage_distance_km must be finite and non-negative.", call. = FALSE)
  }
  chunk_size <- max(1L, as.integer(chunk_size))
  nearest_distance <- rep(NA_real_, nrow(grid0))
  nearest_site <- rep(NA_character_, nrow(grid0))
  chunks <- split(seq_len(nrow(grid0)),
                  ceiling(seq_len(nrow(grid0)) / chunk_size))
  for (idx in chunks) {
    for (i in idx) {
      d <- hee_great_circle_distance_km(
        grid0$lon[[i]], grid0$lat[[i]], sites$lon, sites$lat
      )
      j <- which.min(d)
      nearest_distance[[i]] <- d[[j]]
      nearest_site[[i]] <- as.character(sites$site_id[[j]])
    }
  }
  grid0$nearest_survey_site_id <- nearest_site
  grid0$nearest_survey_distance_km <- nearest_distance
  grid0$survey_covered <- as.integer(
    is.finite(nearest_distance) & nearest_distance <= coverage_distance_km
  )
  grid0
}

#' Summarise actual grid coverage for every palaeo time slice
#'
#' @param grid Palaeoenvironmental cell-time table.
#' @param predictors Predictor columns required for complete environmental rows.
#' @param landmask_col Land-mask column.
#' @return One row per `time_ma` with total, land, and complete-cell counts.
#' @export
hee_case03_grid_inventory <- function(grid, predictors,
                                      landmask_col = "land_mask_dem") {
  grid <- as.data.frame(grid)
  .require_cols(grid, c("cell_id", "time_ma", landmask_col), "grid")
  .require_cols(grid, predictors, "grid")
  .hee_check_unique_keys(grid, c("cell_id", "time_ma"), "grid")
  .hee_validate_time_values(grid$time_ma, "grid$time_ma")
  times <- sort(unique(as.numeric(grid$time_ma)), decreasing = TRUE)
  rows <- lapply(times, function(tm) {
    z <- grid[as.numeric(grid$time_ma) == tm, , drop = FALSE]
    land <- is.finite(as.numeric(z[[landmask_col]])) &
      as.numeric(z[[landmask_col]]) > 0
    complete <- land & stats::complete.cases(z[, predictors, drop = FALSE])
    data.frame(
      time_ma = tm,
      n_total_cells = nrow(z),
      n_land_cells = sum(land),
      n_complete_environment_cells = sum(complete),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Validate time-specific region assignments against a palaeo grid
#'
#' @param env_grid Cell-time palaeoenvironment grid.
#' @param region_grid Cell-time region assignments.
#' @param require_complete Require every environmental row to have a region.
#' @param coordinate_tolerance Maximum lon/lat mismatch in degrees when both
#'   tables contain coordinates.
#' @return Region assignments with a unique `region_time_id`.
#' @export
hee_case03_validate_region_grid <- function(env_grid, region_grid,
                                            require_complete = TRUE,
                                            coordinate_tolerance = 1e-8) {
  env_grid <- as.data.frame(env_grid)
  region_grid <- as.data.frame(region_grid)
  .require_cols(env_grid, c("cell_id", "time_ma"), "env_grid")
  .require_cols(region_grid, c("cell_id", "time_ma", "region"), "region_grid")
  .hee_check_unique_keys(env_grid, c("cell_id", "time_ma"), "env_grid")
  .hee_check_unique_keys(region_grid, c("cell_id", "time_ma"), "region_grid")
  .hee_validate_time_values(env_grid$time_ma, "env_grid$time_ma")
  .hee_validate_time_values(region_grid$time_ma, "region_grid$time_ma")
  env_key <- paste(env_grid$cell_id, env_grid$time_ma, sep = "\r")
  reg_key <- paste(region_grid$cell_id, region_grid$time_ma, sep = "\r")
  missing <- setdiff(env_key, reg_key)
  if (isTRUE(require_complete) && length(missing) > 0L) {
    stop("region_grid has missing region assignments for ", length(missing),
         " cell_id x time_ma key(s).", call. = FALSE)
  }
  extra <- setdiff(reg_key, env_key)
  if (length(extra) > 0L) {
    warning("region_grid contains ", length(extra),
            " cell-time key(s) absent from env_grid; they were removed.",
            call. = FALSE)
    region_grid <- region_grid[reg_key %in% env_key, , drop = FALSE]
    reg_key <- paste(region_grid$cell_id, region_grid$time_ma, sep = "\r")
  }
  if (all(c("lon", "lat") %in% names(env_grid)) &&
      all(c("lon", "lat") %in% names(region_grid))) {
    j <- match(reg_key, env_key)
    coord_bad <- abs(region_grid$lon - env_grid$lon[j]) > coordinate_tolerance |
      abs(region_grid$lat - env_grid$lat[j]) > coordinate_tolerance
    coord_bad[is.na(coord_bad)] <- TRUE
    if (any(coord_bad)) {
      stop("region_grid coordinates do not align with env_grid for ",
           sum(coord_bad), " cell-time row(s).", call. = FALSE)
    }
  }
  region_grid$region <- as.character(region_grid$region)
  if (anyNA(region_grid$region) || any(!nzchar(region_grid$region))) {
    stop("region_grid contains missing or empty region labels.", call. = FALSE)
  }
  region_grid$region_id <- if ("region_id" %in% names(region_grid)) {
    as.character(region_grid$region_id)
  } else {
    region_grid$region
  }
  region_grid$region_time_id <- paste0(region_grid$time_ma, "Ma::",
                                        region_grid$region_id)
  region_grid
}

#' Determine which Case 03 projection models can legally run
#'
#' @param has_tree,has_regions,has_accessibility,has_static,has_dynamic Logical
#'   indicators for real input availability.
#' @return A model-status table. Missing historical accessibility always blocks
#'   M3-M5; it is never replaced by a mock table.
#' @export
hee_case03_model_availability <- function(has_tree,
                                          has_regions,
                                          has_accessibility,
                                          has_static,
                                          has_dynamic) {
  flags <- vapply(list(has_tree, has_regions, has_accessibility,
                       has_static, has_dynamic), function(z) {
    isTRUE(z) && length(z) == 1L && !is.na(z)
  }, logical(1))
  names(flags) <- c("tree", "regions", "accessibility", "static", "dynamic")
  models <- c("M1_env_only", "M2_env_phylo",
              "M3_env_phylo_bgb_no_dispersal",
              "M4_env_phylo_bgb_static_dispersal",
              "M5_env_phylo_bgb_dynamic_dispersal")
  status <- rep("READY", length(models))
  reason <- rep("all required components are available", length(models))
  if (!flags[["tree"]]) {
    status[2:5] <- "NOT_RUN_MISSING_PHYLOGENY"
    reason[2:5] <- "dated tree/species origin evidence is missing"
  }
  if (!flags[["accessibility"]]) {
    status[3:5] <- "NOT_RUN_MISSING_ACCESSIBILITY"
    reason[3:5] <- "BioGeoBEARS or equivalent accessibility is missing"
  } else if (!flags[["regions"]]) {
    status[3:5] <- "NOT_RUN_MISSING_REGIONS"
    reason[3:5] <- "time-specific region cube is missing"
  }
  if (status[[4]] == "READY" && !flags[["static"]]) {
    status[[4]] <- "NOT_RUN_MISSING_STATIC_CONNECTIVITY"
    reason[[4]] <- "static connectivity is missing"
  }
  if (status[[5]] == "READY" && !flags[["dynamic"]]) {
    status[[5]] <- "NOT_RUN_MISSING_DYNAMIC_CONNECTIVITY"
    reason[[5]] <- "dynamic connectivity is missing"
  }
  data.frame(model_id = models, status = status, reason = reason,
             stringsAsFactors = FALSE)
}

#' Read and align modern Case 03 community data
#'
#' @param comm_path,sites_path CSV paths.
#' @param response_type Community response type passed to [hee_check_comm()].
#' @return A list with numeric site-by-species matrix `Y`, aligned `sites`, and
#'   `species`. Site identifiers must match exactly and are never generated.
#' @export
hee_case03_read_community <- function(comm_path, sites_path,
                                      response_type = "presence") {
  if (!file.exists(comm_path)) stop("Missing required file: ", comm_path,
                                    call. = FALSE)
  if (!file.exists(sites_path)) stop("Missing required file: ", sites_path,
                                     call. = FALSE)
  comm <- utils::read.csv(comm_path, check.names = FALSE,
                          stringsAsFactors = FALSE)
  sites <- utils::read.csv(sites_path, check.names = FALSE,
                           stringsAsFactors = FALSE)
  .require_cols(comm, "site_id", "comm.csv")
  .require_cols(sites, c("site_id", "lon", "lat"), "sites.csv")
  .hee_check_unique_keys(comm, "site_id", "comm.csv")
  .hee_check_unique_keys(sites, "site_id", "sites.csv")
  if (!setequal(as.character(comm$site_id), as.character(sites$site_id))) {
    missing_sites <- setdiff(as.character(comm$site_id), as.character(sites$site_id))
    extra_sites <- setdiff(as.character(sites$site_id), as.character(comm$site_id))
    stop("comm.csv and sites.csv site_id values do not match exactly. ",
         "Missing in sites: ", paste(missing_sites, collapse = ", "),
         "; extra in sites: ", paste(extra_sites, collapse = ", "),
         call. = FALSE)
  }
  species <- setdiff(names(comm), "site_id")
  if (length(species) == 0L) stop("comm.csv has no species columns.",
                                  call. = FALSE)
  Y <- as.matrix(comm[, species, drop = FALSE])
  suppressWarnings(storage.mode(Y) <- "numeric")
  if (any(!is.finite(Y) & !is.na(Y))) {
    stop("comm.csv contains non-finite community values.", call. = FALSE)
  }
  rownames(Y) <- as.character(comm$site_id)
  chk <- hee_check_comm(Y, response_type = response_type)
  if (!isTRUE(chk$ok)) {
    stop("Invalid comm.csv: ", paste(chk$warnings, collapse = "; "),
         call. = FALSE)
  }
  sites <- sites[match(rownames(Y), as.character(sites$site_id)), , drop = FALSE]
  if (any(!is.finite(as.numeric(sites$lon))) ||
      any(!is.finite(as.numeric(sites$lat)))) {
    stop("sites.csv lon/lat must be finite numeric coordinates.", call. = FALSE)
  }
  sites$lon <- as.numeric(sites$lon)
  sites$lat <- as.numeric(sites$lat)
  list(Y = Y, sites = sites, species = species)
}

#' Read a Case 03 cell-time layer
#'
#' Supports CSV tables, RDS data.frames, and compact HmscEcoEvo time cubes.
#' Raster products that are not compact time cubes should first be indexed and
#' converted to a `cell_id` x `time_ma` table so alignment remains auditable.
#'
#' @param path Input CSV/RDS path.
#' @param times Optional required ages.
#' @param variables Optional compact-cube variables to extract.
#' @return A cell-time data.frame.
#' @export
hee_case03_read_cell_time_layer <- function(path, times = NULL,
                                             variables = NULL) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stop("Missing cell-time layer: ", path %||% "", call. = FALSE)
  }
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    out <- utils::read.csv(path, check.names = FALSE,
                           stringsAsFactors = FALSE)
  } else if (ext == "rds") {
    obj <- readRDS(path)
    if (is.data.frame(obj)) {
      out <- obj
    } else if (inherits(obj, "hee_timecube") || is.environment(obj)) {
      cube <- .hee_as_timecube(obj)
      variables <- variables %||% cube$variables
      times <- times %||% hee_time_axis_from_env_cube(cube)
      out <- hee_make_paleo_grid(cube, times = times,
                                 variables = variables,
                                 land_only = FALSE)
    } else {
      stop("RDS cell-time layer must contain a data.frame or compact time cube.",
           call. = FALSE)
    }
  } else {
    stop("Cell-time layers currently require .csv or .rds input.",
         call. = FALSE)
  }
  out <- as.data.frame(out)
  if (!"time_ma" %in% names(out) && "age_ma" %in% names(out)) {
    names(out)[names(out) == "age_ma"] <- "time_ma"
  }
  if (!"cell_id" %in% names(out) && all(c("lon", "lat") %in% names(out))) {
    out$cell_id <- .hee_cell_id(out$lon, out$lat)
  }
  .require_cols(out, c("cell_id", "time_ma"), "cell-time layer")
  .hee_validate_time_values(out$time_ma, "cell-time layer$time_ma")
  if (!is.null(times)) {
    missing_times <- setdiff(as.numeric(times), unique(as.numeric(out$time_ma)))
    if (length(missing_times) > 0L) {
      stop("Cell-time layer is missing age(s): ",
           paste(missing_times, collapse = ", "), " Ma.", call. = FALSE)
    }
    out <- out[as.numeric(out$time_ma) %in% as.numeric(times), , drop = FALSE]
  }
  .hee_check_unique_keys(out, c("cell_id", "time_ma"), "cell-time layer")
  out
}

#' Check real Case 03 species-name alignment
#'
#' @param species Community-matrix species names.
#' @param traits Optional trait table with a species column or species row names.
#' @param phy Optional `phylo` object.
#' @param tip_ranges Optional species-region table/matrix.
#' @param accessibility Optional BioGeoBEARS accessibility table.
#' @return An alignment table and invisibly validated inputs.
#' @export
hee_case03_check_species_alignment <- function(species,
                                               traits = NULL,
                                               phy = NULL,
                                               tip_ranges = NULL,
                                               accessibility = NULL) {
  species <- as.character(species)
  if (anyNA(species) || any(!nzchar(species)) || anyDuplicated(species)) {
    stop("Community species names must be unique non-empty strings.",
         call. = FALSE)
  }
  sets <- list(community = species)
  if (!is.null(traits)) {
    tr <- as.data.frame(traits)
    tr_species <- if ("species" %in% names(tr)) as.character(tr$species) else rownames(tr)
    sets$traits <- tr_species
  }
  if (!is.null(phy)) {
    if (!inherits(phy, "phylo")) stop("phy must inherit from `phylo`.",
                                      call. = FALSE)
    sets$tree <- as.character(phy$tip.label)
  }
  if (!is.null(tip_ranges)) {
    trng <- as.data.frame(tip_ranges)
    sets$tip_ranges <- if ("species" %in% names(trng))
      as.character(trng$species) else rownames(trng)
  }
  if (!is.null(accessibility)) {
    acc <- as.data.frame(accessibility)
    .require_cols(acc, "species", "accessibility")
    sets$accessibility <- unique(as.character(acc$species))
  }
  rows <- do.call(rbind, lapply(names(sets), function(nm) {
    values <- sets[[nm]]
    if (is.null(values) || anyNA(values) || any(!nzchar(values)) ||
        anyDuplicated(values)) {
      stop(nm, " species names must be unique non-empty strings.",
           call. = FALSE)
    }
    data.frame(source = nm,
               n_species = length(values),
               missing_from_source = paste(setdiff(species, values), collapse = ";"),
               extra_in_source = paste(setdiff(values, species), collapse = ";"),
               exact_match = setequal(species, values),
               stringsAsFactors = FALSE)
  }))
  bad <- rows$source != "community" & !rows$exact_match
  if (any(bad)) {
    stop("Species names do not align exactly across: ",
         paste(rows$source[bad], collapse = ", "), ".", call. = FALSE)
  }
  rows
}

#' Project HMSC Beta posterior draws through a palaeoenvironment time cube
#'
#' Calculates environmental fixed-effect predictions only. Modern site-level
#' random effects are deliberately excluded because they are not transferable
#' to palaeo grids. Computation is chunked by cell and species so all posterior
#' draws need not be materialised for all cells and taxa simultaneously.
#'
#' @param beta_draws Numeric array in draw x species x axis order.
#' @param env_cube A time cube accepted by [hee_load_timecube()].
#' @param recipe Optional locked 0 Ma recipe from [hee_lock_recipe()].
#' @param times Ages in Ma; larger values are older.
#' @param variables Environmental predictor names.
#' @param link Inverse link, normally `"probit"` or `"logit"`.
#' @param intercept_axis Intercept axis name.
#' @param landmask Land-mask variable.
#' @param land_only Keep only valid land cells.
#' @param cell_chunk_size Number of cells processed together.
#' @param max_draws Optional maximum number of evenly spaced posterior draws.
#' @return Long table with posterior mean, SD, 2.5%, median, 97.5%, and the
#'   number of posterior draws. `suitability` aliases `suitability_mean` for
#'   compatibility with [hee_combine()].
#' @export
hee_project_hmsc_posterior_timecube <- function(beta_draws,
                                                 env_cube,
                                                 recipe = NULL,
                                                 times,
                                                 variables = NULL,
                                                 link = "probit",
                                                 intercept_axis = "(Intercept)",
                                                 landmask = "land_mask_dem",
                                                 land_only = TRUE,
                                                 cell_chunk_size = 250L,
                                                 max_draws = Inf) {
  if (!is.array(beta_draws) || length(dim(beta_draws)) != 3L) {
    stop("beta_draws must be a draw x species x axis array.", call. = FALSE)
  }
  dn <- dimnames(beta_draws)
  if (is.null(dn[[2]]) || is.null(dn[[3]]) ||
      any(!nzchar(dn[[2]])) || any(!nzchar(dn[[3]]))) {
    stop("beta_draws must have species and axis dimnames.", call. = FALSE)
  }
  cell_chunk_size <- as.integer(cell_chunk_size)
  if (!is.finite(cell_chunk_size) || cell_chunk_size < 1L) {
    stop("cell_chunk_size must be a positive integer.", call. = FALSE)
  }
  n_draw <- dim(beta_draws)[1]
  if (n_draw < 2L) stop("At least two posterior draws are required.", call. = FALSE)
  if (is.finite(max_draws)) {
    max_draws <- as.integer(max_draws)
    if (max_draws < 2L) stop("max_draws must be at least 2.", call. = FALSE)
    keep <- unique(round(seq(1, n_draw, length.out = min(max_draws, n_draw))))
    beta_draws <- beta_draws[keep, , , drop = FALSE]
    n_draw <- length(keep)
  }
  species <- dn[[2]]
  axes <- dn[[3]]
  variables <- variables %||% setdiff(axes, intercept_axis)
  missing_axes <- setdiff(setdiff(axes, intercept_axis), variables)
  if (length(missing_axes) > 0L) {
    stop("variables are missing Beta axis/axes: ",
         paste(missing_axes, collapse = ", "), call. = FALSE)
  }
  times <- sort(unique(as.numeric(times)), decreasing = TRUE)
  .hee_validate_time_values(times, "times")
  cube <- .hee_as_timecube(env_cube)
  axes_no_intercept <- setdiff(axes, intercept_axis)
  out <- list()
  k <- 1L
  for (tm in times) {
    grid <- hee_make_paleo_grid(cube, times = tm, variables = variables,
                                landmask = landmask, land_only = land_only)
    X <- if (is.null(recipe)) grid[, variables, drop = FALSE] else
      hee_apply_recipe(recipe, grid)
    missing_predictors <- setdiff(axes_no_intercept, names(X))
    if (length(missing_predictors) > 0L) {
      stop("Palaeo predictors are missing Beta axis/axes: ",
           paste(missing_predictors, collapse = ", "), call. = FALSE)
    }
    chunks <- split(seq_len(nrow(grid)),
                    ceiling(seq_len(nrow(grid)) / cell_chunk_size))
    for (idx in chunks) {
      Xmat <- as.matrix(X[idx, axes_no_intercept, drop = FALSE])
      finite_rows <- if (ncol(Xmat) == 0L) rep(TRUE, nrow(Xmat)) else
        apply(Xmat, 1, function(z) all(is.finite(z)))
      ids <- grid[idx, intersect(c("time_ma", "lon", "lat", "cell_id",
                                   landmask, "land_area_km2"), names(grid)),
                  drop = FALSE]
      for (sp in species) {
        B <- matrix(beta_draws[, sp, axes_no_intercept, drop = FALSE],
                    nrow = n_draw, ncol = length(axes_no_intercept))
        eta <- if (length(axes_no_intercept) == 0L) {
          matrix(0, nrow = nrow(Xmat), ncol = n_draw)
        } else {
          Xmat %*% t(B)
        }
        if (intercept_axis %in% axes) {
          eta <- sweep(eta, 2, beta_draws[, sp, intercept_axis], `+`)
        }
        pred <- .hee_link_inverse(eta, link)
        pred[!finite_rows, ] <- NA_real_
        q <- t(apply(pred, 1, stats::quantile,
                     probs = c(0.025, 0.5, 0.975), na.rm = TRUE,
                     names = FALSE))
        row <- cbind(
          ids,
          species = sp,
          suitability_mean = rowMeans(pred, na.rm = TRUE),
          suitability_sd = apply(pred, 1, stats::sd, na.rm = TRUE),
          suitability_q025 = q[, 1],
          suitability_q50 = q[, 2],
          suitability_q975 = q[, 3],
          n_posterior_draws = n_draw,
          stringsAsFactors = FALSE
        )
        row$suitability_mean[!finite_rows] <- NA_real_
        row$suitability_sd[!finite_rows] <- NA_real_
        row$suitability <- row$suitability_mean
        out[[k]] <- row
        k <- k + 1L
      }
    }
  }
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  class(ans) <- c("hee_suitability_posterior", "hee_suitability", class(ans))
  ans
}

#' Compute palaeotopographic neighbourhood metrics
#'
#' Uses `terra` on each complete time slice. Raw elevation is never relabelled
#' as slope or disturbance. The returned topographic heterogeneity index is a
#' bounded descriptive proxy built from local SD, relief, TRI, and roughness;
#' `disturbance` is not generated.
#'
#' @param grid Full cell-time palaeo grid with longitude, latitude, elevation,
#'   and land mask.
#' @param elevation_col,landmask_col Column names.
#' @param window Odd focal-window width in cells.
#' @param crs Coordinate reference system for the palaeo grid.
#' @return Input rows plus local elevation mean, SD, range, slope (degrees),
#'   TRI, roughness, local relief, topographic heterogeneity, and metadata.
#' @export
hee_case03_topography <- function(grid,
                                  elevation_col = "elev",
                                  landmask_col = "land_mask_dem",
                                  window = 3L,
                                  crs = "EPSG:4326") {
  .require_pkg("terra", "palaeotopographic metrics")
  grid <- as.data.frame(grid)
  .require_cols(grid, c("cell_id", "time_ma", "lon", "lat",
                        elevation_col, landmask_col), "grid")
  .hee_check_unique_keys(grid, c("cell_id", "time_ma"), "grid")
  window <- as.integer(window)
  if (!is.finite(window) || window < 3L || window %% 2L == 0L) {
    stop("window must be an odd integer >= 3.", call. = FALSE)
  }
  one_time <- function(z) {
    land <- is.finite(as.numeric(z[[landmask_col]])) &
      as.numeric(z[[landmask_col]]) > 0
    xyz <- data.frame(x = z$lon, y = z$lat,
                      elevation = ifelse(land, as.numeric(z[[elevation_col]]), NA_real_))
    lon_values <- sort(unique(xyz$x[is.finite(xyz$x)]))
    lon_step <- if (length(lon_values) > 1L) {
      stats::median(diff(lon_values), na.rm = TRUE)
    } else {
      NA_real_
    }
    global_longitude <- is.finite(lon_step) && lon_step > 0 &&
      diff(range(lon_values)) + lon_step >= 359
    if (global_longitude) {
      n_pad <- (window - 1L) / 2L
      left_values <- utils::tail(lon_values, n_pad)
      right_values <- utils::head(lon_values, n_pad)
      left <- xyz[xyz$x %in% left_values, , drop = FALSE]
      right <- xyz[xyz$x %in% right_values, , drop = FALSE]
      left$x <- left$x - 360
      right$x <- right$x + 360
      xyz_for_raster <- rbind(left, xyz, right)
    } else {
      xyz_for_raster <- xyz
    }
    r <- terra::rast(xyz_for_raster, type = "xyz", crs = crs)
    w <- matrix(1, nrow = window, ncol = window)
    local_mean <- terra::focal(r, w = w, fun = mean, na.rm = TRUE,
                               na.policy = "all", expand = TRUE)
    local_sd <- terra::focal(r, w = w, fun = stats::sd, na.rm = TRUE,
                             na.policy = "all", expand = TRUE)
    local_min <- terra::focal(r, w = w, fun = min, na.rm = TRUE,
                              na.policy = "all", expand = TRUE)
    local_max <- terra::focal(r, w = w, fun = max, na.rm = TRUE,
                              na.policy = "all", expand = TRUE)
    slope <- terra::terrain(r, v = "slope", unit = "degrees", neighbors = 8)
    tri <- terra::terrain(r, v = "TRI", neighbors = 8)
    rough <- terra::terrain(r, v = "roughness", neighbors = 8)
    extract_values <- function(x) {
      pts <- terra::vect(z[, c("lon", "lat")], geom = c("lon", "lat"), crs = crs)
      terra::extract(x, pts)[, 2]
    }
    z$mean_elevation <- extract_values(local_mean)
    z$elevation_sd <- extract_values(local_sd)
    z$elevation_range <- extract_values(local_max) - extract_values(local_min)
    z$slope_deg <- extract_values(slope)
    z$TRI <- extract_values(tri)
    z$roughness <- extract_values(rough)
    z$local_relief <- z$elevation_range
    components <- cbind(z$elevation_sd, z$local_relief, z$TRI, z$roughness)
    scale_component <- function(x) {
      pos <- x[is.finite(x) & x > 0]
      scale <- stats::median(pos, na.rm = TRUE)
      if (!is.finite(scale) || scale <= 0) scale <- 1
      out <- pmax(x, 0) / (pmax(x, 0) + scale)
      out[!is.finite(out)] <- NA_real_
      out
    }
    scaled <- apply(components, 2, scale_component)
    z$topographic_heterogeneity <- rowMeans(scaled, na.rm = TRUE)
    z$topographic_heterogeneity[!land | !is.finite(z$topographic_heterogeneity)] <- NA_real_
    z$topography_window_cells <- window
    z$topography_crs <- crs
    z$topography_longitude_boundary <- if (global_longitude) {
      "cyclic_longitude_padding"
    } else {
      "non_global_grid_no_longitude_wrap"
    }
    z
  }
  out <- lapply(split(grid, grid$time_ma), one_time)
  ans <- do.call(rbind, out)
  ans <- ans[order(-as.numeric(ans$time_ma), ans$cell_id), , drop = FALSE]
  rownames(ans) <- NULL
  ans
}

#' Build time-specific region distances in kilometres
#'
#' Region centroids use a circular mean for longitude, avoiding the 180-degree
#' seam. The default is an explicitly labelled centroid great-circle
#' approximation, not a least-cost path. Barrier and corridor summaries can be
#' supplied from the full palaeo grid and retained for later connectivity
#' calculations.
#'
#' @param region_grid Cell-time grid with region, lon, and lat.
#' @param distance_scale_km Structural distance-decay scale in kilometres.
#' @return Ordered region-pair rows with `distance_km`, `paleodistance`, method,
#'   units, and structural distance decay.
#' @export
hee_case03_region_distances <- function(region_grid,
                                        distance_scale_km = 1000) {
  x <- as.data.frame(region_grid)
  .require_cols(x, c("cell_id", "time_ma", "region", "lon", "lat"),
                "region_grid")
  if (!is.finite(distance_scale_km) || distance_scale_km <= 0) {
    stop("distance_scale_km must be a positive distance in kilometres.",
         call. = FALSE)
  }
  circular_mean <- function(deg) {
    rad <- deg[is.finite(deg)] * pi / 180
    if (length(rad) == 0L) return(NA_real_)
    atan2(mean(sin(rad)), mean(cos(rad))) * 180 / pi
  }
  centres <- do.call(rbind, lapply(split(x, list(x$time_ma, x$region), drop = TRUE),
                                   function(z) {
    data.frame(time_ma = z$time_ma[[1]], region = z$region[[1]],
               centre_lon = circular_mean(z$lon),
               centre_lat = mean(z$lat, na.rm = TRUE),
               n_region_cells = nrow(z), stringsAsFactors = FALSE)
  }))
  rows <- lapply(split(centres, centres$time_ma), function(z) {
    if (nrow(z) < 2L) return(data.frame())
    ij <- expand.grid(i = seq_len(nrow(z)), j = seq_len(nrow(z)))
    ij <- ij[ij$i != ij$j, , drop = FALSE]
    data.frame(
      from_region = z$region[ij$i],
      to_region = z$region[ij$j],
      time_ma = z$time_ma[ij$i],
      from_lon = z$centre_lon[ij$i],
      from_lat = z$centre_lat[ij$i],
      to_lon = z$centre_lon[ij$j],
      to_lat = z$centre_lat[ij$j],
      distance_km = hee_great_circle_distance_km(
        z$centre_lon[ij$i], z$centre_lat[ij$i],
        z$centre_lon[ij$j], z$centre_lat[ij$j]
      ),
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[vapply(rows, nrow, integer(1)) > 0L]
  if (length(rows) == 0L) return(data.frame())
  out <- do.call(rbind, rows)
  out$paleodistance <- out$distance_km
  out$least_cost_distance <- NA_real_
  out$distance_method <- "region_centroid_great_circle_approximation"
  out$distance_unit <- "km"
  out$distance_scale_km <- distance_scale_km
  out$structural_distance_decay <- exp(-out$distance_km / distance_scale_km)
  rownames(out) <- NULL
  out
}

#' Add spherical cell area and inferred grid resolution
#'
#' Cell area is calculated with [terra::cellSize()] in square kilometres from
#' the actual longitude-latitude grid at each time slice. It is not copied from
#' a modern sample or assumed constant across latitude.
#'
#' @param grid Regular cell-time longitude-latitude grid.
#' @return `grid` plus `cell_area_km2`, `grid_resolution_lon_deg`, and
#'   `grid_resolution_lat_deg`.
#' @export
hee_case03_add_cell_area <- function(grid) {
  .require_pkg("terra", "palaeogeographic cell area")
  grid <- as.data.frame(grid)
  .require_cols(grid, c("cell_id", "time_ma", "lon", "lat"), "grid")
  .hee_check_unique_keys(grid, c("cell_id", "time_ma"), "grid")
  one_time <- function(z) {
    dummy <- data.frame(lon = z$lon, lat = z$lat, value = 1)
    r <- terra::rast(dummy, type = "xyz", crs = "EPSG:4326")
    area <- terra::cellSize(r, unit = "km")
    pts <- terra::vect(z[, c("lon", "lat")], geom = c("lon", "lat"),
                       crs = "EPSG:4326")
    z$cell_area_km2 <- terra::extract(area, pts)[, 2]
    lon_values <- sort(unique(z$lon[is.finite(z$lon)]))
    lat_values <- sort(unique(z$lat[is.finite(z$lat)]))
    z$grid_resolution_lon_deg <- if (length(lon_values) > 1L) {
      stats::median(diff(lon_values), na.rm = TRUE)
    } else {
      NA_real_
    }
    z$grid_resolution_lat_deg <- if (length(lat_values) > 1L) {
      stats::median(diff(lat_values), na.rm = TRUE)
    } else {
      NA_real_
    }
    z
  }
  out <- lapply(split(grid, grid$time_ma), one_time)
  ans <- do.call(rbind, out)
  ans <- ans[order(-as.numeric(ans$time_ma), ans$cell_id), , drop = FALSE]
  rownames(ans) <- NULL
  ans
}

#' Import and model-average external BioGeoBEARS accessibility
#'
#' Multiple candidate models must remain separate in the input. For each model,
#' `species x region x time_ma` keys are required to be unique. The supplied
#' model weights are normalised once and accessibility is averaged as
#' \deqn{A_{jrt}=\sum_m w_m A_{jrtm}.}
#' This function imports external historical-biogeographic evidence; it does not
#' fit BioGeoBEARS or infer ancestral ranges internally.
#'
#' @param x External accessibility table with `species`, `region`, `time_ma`,
#'   `model`, `model_weight`, and `accessibility`.
#' @return A list containing model-averaged `accessibility`, normalised
#'   `model_weights`, a model-level `model_summary`, and `provenance`.
#' @export
hee_case03_import_bgb_accessibility <- function(x) {
  x <- as.data.frame(x)
  required <- c("species", "region", "time_ma", "model", "model_weight",
                "accessibility")
  .require_cols(x, required, "bgb_accessibility")
  if (nrow(x) == 0L) {
    stop("bgb_accessibility must contain at least one row.", call. = FALSE)
  }
  if (anyNA(x$species) || anyNA(x$region) || anyNA(x$model) ||
      any(!nzchar(as.character(x$species))) ||
      any(!nzchar(as.character(x$region))) ||
      any(!nzchar(as.character(x$model)))) {
    stop("bgb_accessibility species, region, and model must be non-missing.",
         call. = FALSE)
  }
  .hee_validate_time_values(x$time_ma, "bgb_accessibility$time_ma")
  .hee_check_unique_keys(x, c("model", "species", "region", "time_ma"),
                         "bgb_accessibility")
  x$accessibility <- .hee_clip01(x$accessibility)
  if (anyNA(x$accessibility)) {
    stop("bgb_accessibility$accessibility must be finite; encode absence as 0.",
         call. = FALSE)
  }
  weight_by_model <- lapply(split(x$model_weight, x$model), function(z) {
    z <- unique(suppressWarnings(as.numeric(z)))
    z <- z[is.finite(z)]
    if (length(z) != 1L || z < 0) {
      stop("Each BioGeoBEARS model requires one finite non-negative model_weight.",
           call. = FALSE)
    }
    z
  })
  weights <- unlist(weight_by_model, use.names = TRUE)
  if (!is.finite(sum(weights)) || sum(weights) <= 0) {
    stop("BioGeoBEARS model weights must sum to a positive value.", call. = FALSE)
  }
  weights <- weights / sum(weights)
  model_tables <- split(x[, c("species", "region", "time_ma", "accessibility"),
                          drop = FALSE], x$model)
  averaged <- hee_bgb_accessibility(model_tables, weights = weights)
  model_weights <- data.frame(
    model = names(weights), model_weight = as.numeric(weights),
    stringsAsFactors = FALSE
  )
  summary_cols <- intersect(c("model", "logLik", "n_parameters", "d", "e",
                              "j", "AIC", "AICc", "delta_AICc"), names(x))
  model_summary <- unique(x[, summary_cols, drop = FALSE])
  model_summary <- merge(model_weights, model_summary, by = "model",
                         all.x = TRUE, sort = FALSE)
  provenance <- data.frame(
    source = "external_biogeographic_model_output",
    execution_status = "IMPORTED_EXTERNAL_MODEL_AVERAGED_ACCESSIBILITY",
    fitted_by_case03 = FALSE,
    n_models = nrow(model_weights),
    interpretation = paste(
      "Historical accessibility constraint; not HMSC suitability,",
      "not a causal dispersal estimate, and not fitted by this importer."
    ),
    stringsAsFactors = FALSE
  )
  list(accessibility = averaged, model_weights = model_weights,
       model_summary = model_summary, provenance = provenance)
}

#' Static geological dispersal filter from external source evidence
#'
#' Source regions are defined by external accessibility at or above
#' `source_threshold`. For each target region the filter is
#' \deqn{D_{static,jrt}=\max_{r':A_{jr't}\geq\tau} C_{jr'rt}.}
#' Accessibility is used to identify sources and is not multiplied into this
#' filter a second time; target accessibility remains the separate `A_BGB`
#' factor in M4.
#'
#' @param target Unique species-region-time target table.
#' @param accessibility Model-averaged species-region-time accessibility.
#' @param connectivity Species-specific from-region/to-region/time table.
#' @param source_threshold Minimum accessibility defining a source region.
#' @param connectivity_col Connectivity column in `[0,1]`.
#' @return `target` plus `D_static`, selected source region, source evidence,
#'   and the threshold used.
#' @export
hee_case03_static_filter_from_connectivity <- function(
    target, accessibility, connectivity, source_threshold = 0.5,
    connectivity_col = "connectivity") {
  target <- as.data.frame(target)
  accessibility <- as.data.frame(accessibility)
  connectivity <- as.data.frame(connectivity)
  .require_cols(target, c("species", "region", "time_ma"), "target")
  .require_cols(accessibility,
                c("species", "region", "time_ma", "accessibility"),
                "accessibility")
  .require_cols(connectivity,
                c("species", "from_region", "to_region", "time_ma",
                  connectivity_col), "connectivity")
  .hee_check_unique_keys(target, c("species", "region", "time_ma"), "target")
  .hee_check_unique_keys(accessibility, c("species", "region", "time_ma"),
                         "accessibility")
  .hee_check_unique_keys(connectivity,
                         c("species", "from_region", "to_region", "time_ma"),
                         "connectivity")
  source_threshold <- as.numeric(source_threshold)[1]
  if (!is.finite(source_threshold) || source_threshold < 0 ||
      source_threshold > 1) {
    stop("source_threshold must be within [0,1].", call. = FALSE)
  }
  accessibility$accessibility <- .hee_clip01(accessibility$accessibility)
  accessibility$accessibility[is.na(accessibility$accessibility)] <- 0
  connectivity[[connectivity_col]] <- .hee_clip01(connectivity[[connectivity_col]])
  connectivity[[connectivity_col]][is.na(connectivity[[connectivity_col]])] <- 0
  rows <- lapply(seq_len(nrow(target)), function(i) {
    q <- target[i, , drop = FALSE]
    src <- accessibility[
      accessibility$species == q$species &
        accessibility$time_ma == q$time_ma &
        accessibility$accessibility > 0 &
        accessibility$accessibility >= source_threshold,
      , drop = FALSE
    ]
    edge <- connectivity[
      connectivity$species == q$species &
        connectivity$time_ma == q$time_ma &
        connectivity$to_region == q$region &
        connectivity$from_region %in% src$region,
      , drop = FALSE
    ]
    if (nrow(edge) == 0L) {
      return(cbind(q, D_static = 0, source_region = NA_character_,
                   source_accessibility = 0,
                   source_threshold = source_threshold))
    }
    j <- which.max(edge[[connectivity_col]])
    source_region <- as.character(edge$from_region[[j]])
    source_accessibility <- src$accessibility[match(source_region, src$region)]
    cbind(q, D_static = edge[[connectivity_col]][[j]],
          source_region = source_region,
          source_accessibility = source_accessibility,
          source_threshold = source_threshold)
  })
  out <- do.call(rbind, rows)
  out$D_static <- .hee_clip01(out$D_static)
  rownames(out) <- NULL
  out
}

#' Initialise the legacy Case 03 M5 dynamic filter from external range support
#'
#' This helper belongs to the older product-form M1-M5 scenario workflow that
#' is retained for backward-compatible sensitivity contrasts. It is not used by
#' the HmscEE 1.0 nested Dynamic Earth-Biota Assembly core. In the legacy M5
#' model, continuous historical accessibility is already present as the
#' separate factor `A_BGB`; at the oldest processed time, and when a lineage
#' first crosses its externally supplied origin age, this helper therefore
#' converts accessibility to a support indicator rather than multiplying the
#' same probability twice:
#' \deqn{D_{init,jrt}=I(A_{BGB,jrt} \geq \tau).}
#' The legacy M5 scenario then remains
#' `S_HMSC * E_phylo * A_BGB * D_dynamic * L_land`. For the main HmscEE 1.0
#' model, use BioGeoBEARS/BSM `R` with [hee_bsm_region_history()] and
#' within-region occupancy with [hee_nested_region_cell_occupancy()].
#' The threshold is a sensitivity parameter, not an estimated dispersal rate.
#' Missing species-region-time matches are conservative zero.
#'
#' @param target Unique species-region-time target table.
#' @param accessibility External model-averaged accessibility with columns
#'   `species`, `region`, `time_ma`, and `accessibility`.
#' @param source_threshold Minimum accessibility interpreted as supported
#'   initial range membership.
#' @return `target` plus `source_accessibility`, binary `D_dynamic`, the
#'   threshold, and an explicit initialisation label.
#' @export
hee_case03_initial_dynamic_filter <- function(
    target, accessibility, source_threshold = 0.5) {
  target <- as.data.frame(target)
  accessibility <- as.data.frame(accessibility)
  keys <- c("species", "region", "time_ma")
  .require_cols(target, keys, "target")
  .require_cols(accessibility, c(keys, "accessibility"), "accessibility")
  .hee_check_unique_keys(target, keys, "target")
  .hee_check_unique_keys(accessibility, keys, "accessibility")
  source_threshold <- as.numeric(source_threshold)[1]
  if (!is.finite(source_threshold) || source_threshold < 0 ||
      source_threshold > 1) {
    stop("source_threshold must be within [0,1].", call. = FALSE)
  }
  target$.row_id <- seq_len(nrow(target))
  joined <- .hee_safe_left_join(
    target,
    accessibility[, c(keys, "accessibility"), drop = FALSE],
    keys = keys, value_cols = "accessibility",
    target_name = "target", source_name = "accessibility"
  )
  source_accessibility <- .hee_clip01(joined$accessibility)
  source_accessibility[!is.finite(source_accessibility)] <- 0
  target$source_accessibility <- source_accessibility
  target$D_dynamic <- as.numeric(
    source_accessibility > 0 & source_accessibility >= source_threshold
  )
  target$dynamic_source_threshold <- source_threshold
  target$dynamic_initialisation <-
    "external_accessibility_support_threshold"
  target$.row_id <- NULL
  target
}

#' Dynamic geological filter from previous occupancy and current paths
#'
#' For each target region, this calculates
#' \deqn{D_{dynamic,jrt}=\max_{r'} P_{jr't_{prev}} C_{jr'rt}.}
#' The supplied connectivity can already combine structural, functional, and
#' climatic path components. Missing source/path matches are conservative zero.
#'
#' @param previous_projection Previous older-time projection with species,
#'   region, and a probability column.
#' @param target Current species-region-time target table.
#' @param connectivity Current species-specific path table.
#' @param probability_col,connectivity_col Input value columns.
#' @return Target table plus `D_dynamic`, source region, source probability,
#'   path connectivity, and source time.
#' @export
hee_case03_dynamic_filter_from_connectivity <- function(
    previous_projection, target, connectivity,
    probability_col = "probability", connectivity_col = "connectivity") {
  previous_projection <- as.data.frame(previous_projection)
  target <- as.data.frame(target)
  connectivity <- as.data.frame(connectivity)
  .require_cols(previous_projection,
                c("species", "region", "time_ma", probability_col),
                "previous_projection")
  .require_cols(target, c("species", "region", "time_ma"), "target")
  .require_cols(connectivity,
                c("species", "from_region", "to_region", "time_ma",
                  connectivity_col), "connectivity")
  .hee_check_unique_keys(target, c("species", "region", "time_ma"), "target")
  .hee_check_unique_keys(connectivity,
                         c("species", "from_region", "to_region", "time_ma"),
                         "connectivity")
  previous <- stats::aggregate(
    previous_projection[[probability_col]],
    previous_projection[, c("species", "region", "time_ma"), drop = FALSE],
    max, na.rm = TRUE
  )
  names(previous)[ncol(previous)] <- ".source_probability"
  previous$.source_probability <- .hee_clip01(previous$.source_probability)
  previous$.source_probability[is.na(previous$.source_probability)] <- 0
  connectivity[[connectivity_col]] <- .hee_clip01(connectivity[[connectivity_col]])
  connectivity[[connectivity_col]][is.na(connectivity[[connectivity_col]])] <- 0
  rows <- lapply(seq_len(nrow(target)), function(i) {
    q <- target[i, , drop = FALSE]
    src <- previous[previous$species == q$species, , drop = FALSE]
    edge <- connectivity[
      connectivity$species == q$species &
        connectivity$time_ma == q$time_ma &
        connectivity$to_region == q$region &
        connectivity$from_region %in% src$region,
      , drop = FALSE
    ]
    if (nrow(edge) == 0L) {
      return(cbind(q, D_dynamic = 0, source_region = NA_character_,
                   source_probability = 0, path_connectivity = 0,
                   dynamic_source_time_ma = if (nrow(src)) max(src$time_ma) else NA_real_))
    }
    src_prob <- src$.source_probability[match(edge$from_region, src$region)]
    score <- src_prob * edge[[connectivity_col]]
    score[!is.finite(score)] <- 0
    j <- which.max(score)
    cbind(q, D_dynamic = score[[j]],
          source_region = as.character(edge$from_region[[j]]),
          source_probability = src_prob[[j]],
          path_connectivity = edge[[connectivity_col]][[j]],
          dynamic_source_time_ma = src$time_ma[match(edge$from_region[[j]],
                                                     src$region)])
  })
  out <- do.call(rbind, rows)
  out$D_dynamic <- .hee_clip01(out$D_dynamic)
  rownames(out) <- NULL
  out
}

#' Propagate HMSC suitability uncertainty through deterministic filters
#'
#' For multiplicative deterministic filters with combined weight `W`, posterior
#' summaries transform as `P_mean = W * S_mean`, `P_sd = W * S_sd`, and
#' `P_q = W * S_q`. This does not represent uncertainty in BioGeoBEARS,
#' phylogeny, land masks, dispersal coefficients, or palaeoenvironment layers.
#'
#' @param projection Output from [hee_combine()] retaining HMSC suitability
#'   posterior summaries.
#' @return Projection plus `probability_mean`, `probability_sd`, and posterior
#'   probability quantiles.
#' @export
hee_case03_propagate_projection_uncertainty <- function(projection) {
  x <- as.data.frame(projection)
  .require_cols(x, c("probability", "suitability_mean", "suitability_sd",
                     "suitability_q025", "suitability_q50",
                     "suitability_q975"), "projection")
  component_cols <- c("accessibility", "phylo_existence", "land_weight",
                      "static_weight", "dynamic_weight", "Q_extrap")
  .require_cols(x, component_cols, "projection")
  parts <- lapply(component_cols, function(nm) {
    z <- .hee_clip01(x[[nm]])
    z[is.na(z)] <- 0
    z
  })
  w <- Reduce(`*`, parts)
  x$probability_mean <- .hee_clip01(x$suitability_mean * w)
  x$probability_sd <- pmax(as.numeric(x$suitability_sd), 0) * w
  x$probability_q025 <- .hee_clip01(x$suitability_q025 * w)
  x$probability_q50 <- .hee_clip01(x$suitability_q50 * w)
  x$probability_q975 <- .hee_clip01(x$suitability_q975 * w)
  hard_zero <- !is.finite(w) | w <= 0
  for (nm in c("probability_mean", "probability_sd", "probability_q025",
               "probability_q50", "probability_q975")) {
    x[[nm]][hard_zero] <- 0
  }
  x$probability <- x$probability_mean
  x
}

#' Estimate species-specific binary suitability thresholds for Case 03 maps
#'
#' Thresholds are derived from modern fitted or cross-validated probabilities
#' against the observed community matrix by maximizing TSS
#' (`sensitivity + specificity - 1`). Species with only one observed class are
#' assigned the neutral threshold 0.5 and flagged as not validated.
#'
#' @param Y Site by species occurrence matrix.
#' @param predicted Site by species predicted occurrence probabilities.
#' @return Data frame with one row per species and threshold diagnostics.
#' @export
hee_case03_species_thresholds <- function(Y, predicted) {
  y <- as.matrix(Y)
  p <- as.matrix(predicted)
  if (!all(dim(y) == dim(p))) {
    stop("Y and predicted must have identical site x species dimensions.",
         call. = FALSE)
  }
  if (is.null(colnames(y))) colnames(y) <- paste0("sp", seq_len(ncol(y)))
  if (is.null(colnames(p))) colnames(p) <- colnames(y)
  rows <- lapply(seq_len(ncol(y)), function(j) {
    obs <- as.numeric(y[, j])
    pr <- .hee_clip01(p[, j])
    ok <- is.finite(obs) & is.finite(pr)
    obs <- obs[ok] > 0
    pr <- pr[ok]
    prevalence <- if (length(obs)) mean(obs) else NA_real_
    if (length(unique(obs)) < 2L || length(unique(pr)) < 2L) {
      return(data.frame(
        species = colnames(y)[[j]], threshold = 0.5,
        sensitivity = NA_real_, specificity = NA_real_, TSS = NA_real_,
        prevalence = prevalence, n_sites = length(obs),
        threshold_source = "neutral_0.5_insufficient_modern_classes",
        stringsAsFactors = FALSE
      ))
    }
    cand <- sort(unique(pr))
    score <- vapply(cand, function(th) {
      pred_bin <- pr >= th
      sens <- sum(pred_bin & obs) / max(sum(obs), 1)
      spec <- sum(!pred_bin & !obs) / max(sum(!obs), 1)
      sens + spec - 1
    }, numeric(1))
    best <- which(score == max(score, na.rm = TRUE))
    if (length(best) > 1L) {
      best <- best[which.min(abs(cand[best] - stats::median(pr)))]
    }
    th <- cand[[best]]
    pred_bin <- pr >= th
    sens <- sum(pred_bin & obs) / max(sum(obs), 1)
    spec <- sum(!pred_bin & !obs) / max(sum(!obs), 1)
    data.frame(
      species = colnames(y)[[j]], threshold = th,
      sensitivity = sens, specificity = spec, TSS = sens + spec - 1,
      prevalence = prevalence, n_sites = length(obs),
      threshold_source = "modern_HMSC_probability_TSS_maximisation",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Build single-species map layers for Case 03
#'
#' Converts posterior probability summaries into species-level map attributes:
#' binary suitability, historical accessibility, accessible suitability, and
#' environment-suitable-but-inaccessible areas. These layers are map products;
#' they are not observed ranges.
#'
#' @param projection One model/time projection table after
#'   [hee_case03_propagate_projection_uncertainty()].
#' @param thresholds Species threshold table from
#'   `hee_case03_species_thresholds()`.
#' @param risk Optional extrapolation-risk table keyed by `cell_id,time_ma`.
#' @return Species by cell map table.
#' @export
hee_case03_species_map_layers <- function(projection, thresholds, risk = NULL) {
  x <- as.data.frame(projection)
  .require_cols(x, c("species", "cell_id", "time_ma", "lon", "lat",
                     "probability_mean", "probability_sd",
                     "probability_q025", "probability_q50",
                     "probability_q975",
                     "suitability_mean"), "projection")
  th <- as.data.frame(thresholds)
  .require_cols(th, c("species", "threshold"), "thresholds")
  if (anyDuplicated(th$species)) {
    stop("thresholds must be unique by species.", call. = FALSE)
  }
  x$threshold <- th$threshold[match(x$species, th$species)]
  if (any(!is.finite(x$threshold))) {
    missing <- unique(x$species[!is.finite(x$threshold)])
    stop("Missing suitability threshold for species: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  x$posterior_mean <- .hee_clip01(x$probability_mean)
  x$posterior_sd <- pmax(as.numeric(x$probability_sd), 0)
  x$posterior_q025 <- .hee_clip01(x$probability_q025)
  x$posterior_q50 <- .hee_clip01(x$probability_q50)
  x$posterior_q975 <- .hee_clip01(x$probability_q975)
  x$environment_suitable <- as.integer(.hee_clip01(x$suitability_mean) >=
                                         x$threshold)
  x$binary_suitability <- as.integer(x$posterior_mean >= x$threshold)
  has_history <- any(c("accessibility", "D_static", "D_dynamic") %in% names(x))
  if (has_history) {
    hist_parts <- list()
    if ("accessibility" %in% names(x)) hist_parts$accessibility <- x$accessibility
    if ("D_static" %in% names(x)) hist_parts$D_static <- x$D_static
    if ("D_dynamic" %in% names(x)) hist_parts$D_dynamic <- x$D_dynamic
    hist_weight <- Reduce(`*`, lapply(hist_parts, function(z) {
      z <- .hee_clip01(z)
      z[is.na(z)] <- 0
      z
    }))
    x$historically_accessible <- as.integer(hist_weight > 0)
  } else {
    x$historically_accessible <- NA_integer_
  }
  x$accessible_suitability <- ifelse(
    is.na(x$historically_accessible), NA_integer_,
    as.integer(x$binary_suitability > 0 & x$historically_accessible > 0)
  )
  x$environment_suitable_but_inaccessible <- ifelse(
    is.na(x$historically_accessible), NA_integer_,
    as.integer(x$environment_suitable > 0 & x$historically_accessible <= 0)
  )
  if (!is.null(risk)) {
    r <- as.data.frame(risk)
    .require_cols(r, c("cell_id", "time_ma"), "risk")
    risk_col <- intersect(c("extrapolation_risk", "extrapolation_score",
                            "no_analog"), names(r))
    if (length(risk_col) > 0L) {
      rr <- r[, c("cell_id", "time_ma", risk_col[[1]]), drop = FALSE]
      names(rr)[[3]] <- "extrapolation_risk"
      x <- merge(x, rr, by = c("cell_id", "time_ma"), all.x = TRUE,
                 sort = FALSE)
    } else {
      x$extrapolation_risk <- NA_real_
    }
  } else {
    x$extrapolation_risk <- NA_real_
  }
  keep <- intersect(
    c("model_id", "species", "cell_id", "time_ma", "lon", "lat",
      "threshold", "posterior_mean", "posterior_sd", "posterior_q025",
      "posterior_q50", "posterior_q975", "binary_suitability",
      "environment_suitable",
      "historically_accessible", "accessible_suitability",
      "environment_suitable_but_inaccessible", "extrapolation_risk",
      "origin_ma", "region", "paleo_cell_id", "coordinate_frame",
      "grid_resolution_lon_deg", "grid_resolution_lat_deg"),
    names(x)
  )
  x[, keep, drop = FALSE]
}

.case03_diversity_from_prob <- function(p, range_weight = NULL) {
  p <- .hee_clip01(p)
  p[!is.finite(p)] <- 0
  richness <- sum(p)
  total <- sum(p)
  rel <- if (total > 0) p / total else p
  shannon <- if (total > 0) -sum(rel[rel > 0] * log(rel[rel > 0])) else 0
  simpson <- if (total > 0) 1 - sum(rel^2) else 0
  we <- if (!is.null(range_weight)) sum(p * range_weight, na.rm = TRUE) else NA_real_
  c(expected_species_richness = richness,
    shannon_diversity = shannon,
    simpson_diversity = simpson,
    weighted_endemism = we)
}

.case03_phylo_helpers <- function(phy, species) {
  if (is.null(phy) || !inherits(phy, "phylo")) return(NULL)
  species <- intersect(species, phy$tip.label)
  if (!length(species)) return(NULL)
  phy <- ape::keep.tip(phy, species)
  ntip <- length(phy$tip.label)
  descendants <- vector("list", nrow(phy$edge))
  children <- split(phy$edge[, 2], phy$edge[, 1])
  get_tips <- function(node) {
    if (node <= ntip) return(node)
    kids <- children[[as.character(node)]]
    unique(unlist(lapply(kids, get_tips), use.names = FALSE))
  }
  for (i in seq_len(nrow(phy$edge))) descendants[[i]] <- get_tips(phy$edge[i, 2])
  branch_species <- matrix(FALSE, nrow = nrow(phy$edge), ncol = ntip,
                           dimnames = list(NULL, phy$tip.label))
  for (i in seq_along(descendants)) {
    branch_species[i, descendants[[i]]] <- TRUE
  }
  list(phy = phy, branch_species = branch_species,
       branch_length = if (is.null(phy$edge.length)) {
         rep(1, nrow(phy$edge))
       } else {
         phy$edge.length
       },
       cophenetic = ape::cophenetic.phylo(phy))
}

.case03_faith_pd <- function(present, phy_helper) {
  if (is.null(phy_helper)) return(NA_real_)
  spp <- intersect(names(present)[present], colnames(phy_helper$branch_species))
  if (!length(spp)) return(0)
  used <- rowSums(phy_helper$branch_species[, spp, drop = FALSE]) > 0
  sum(phy_helper$branch_length[used], na.rm = TRUE)
}

.case03_mpd <- function(present, phy_helper) {
  if (is.null(phy_helper)) return(NA_real_)
  spp <- intersect(names(present)[present], rownames(phy_helper$cophenetic))
  if (length(spp) < 2L) return(0)
  d <- phy_helper$cophenetic[spp, spp, drop = FALSE]
  mean(d[upper.tri(d)], na.rm = TRUE)
}

.case03_trait_helpers <- function(traits, species) {
  if (is.null(traits)) return(NULL)
  tr <- as.data.frame(traits)
  if ("species" %in% names(tr)) rownames(tr) <- tr$species
  numeric_cols <- names(tr)[vapply(tr, is.numeric, logical(1))]
  numeric_cols <- setdiff(numeric_cols, "species")
  species <- intersect(species, rownames(tr))
  if (!length(species) || !length(numeric_cols)) return(NULL)
  mat <- as.matrix(tr[species, numeric_cols, drop = FALSE])
  scale <- apply(mat, 2, function(z) diff(range(z, na.rm = TRUE)))
  scale[!is.finite(scale) | scale == 0] <- 1
  list(mat = sweep(mat, 2, scale, "/"))
}

.case03_fric <- function(present, trait_helper) {
  if (is.null(trait_helper)) return(NA_real_)
  spp <- intersect(names(present)[present], rownames(trait_helper$mat))
  if (length(spp) < 2L) return(0)
  mat <- trait_helper$mat[spp, , drop = FALSE]
  sqrt(sum(apply(mat, 2, function(z) diff(range(z, na.rm = TRUE)))^2))
}

.case03_fdisp <- function(prob, trait_helper) {
  if (is.null(trait_helper)) return(NA_real_)
  spp <- intersect(names(prob), rownames(trait_helper$mat))
  if (!length(spp)) return(0)
  p <- .hee_clip01(prob[spp])
  if (sum(p, na.rm = TRUE) <= 0) return(0)
  mat <- trait_helper$mat[spp, , drop = FALSE]
  center <- colSums(mat * p, na.rm = TRUE) / sum(p, na.rm = TRUE)
  d <- sqrt(rowSums((sweep(mat, 2, center, "-"))^2))
  sum(d * p, na.rm = TRUE) / sum(p, na.rm = TRUE)
}

#' Compute Case 03 diversity map layers from species probabilities
#'
#' Metrics are calculated on each cell in a single time slice. Richness is the
#' sum of occurrence probabilities; Shannon and Simpson use normalized
#' occurrence probabilities and should be interpreted as probability-weighted
#' diversity proxies, not abundance diversity. Thresholded layers are calculated
#' from species-specific modern validation thresholds. Accessible richness and
#' unsuitable-but-accessible/suitable-but-inaccessible counts use the geographic
#' filters already present in the M1--M5 projection table, so geography can
#' visibly change the diversity maps instead of appearing only in separate
#' diagnostic tables. Functional and phylogenetic metrics require observed
#' traits and a dated tree.
#'
#' @param projection One model/time projection table.
#' @param thresholds Species threshold table.
#' @param traits Optional observed trait table.
#' @param phy Optional phylogeny.
#' @param hotspot_quantile Quantile used to flag diversity hotspots.
#' @return Cell-level diversity map table.
#' @export
hee_case03_diversity_map_layers <- function(projection, thresholds,
                                            traits = NULL, phy = NULL,
                                            hotspot_quantile = 0.9) {
  x <- hee_case03_species_map_layers(projection, thresholds)
  spp <- sort(unique(x$species))
  range_size <- stats::aggregate(posterior_mean ~ species, x, sum,
                                 na.rm = TRUE)
  range_weight <- stats::setNames(1 / pmax(range_size$posterior_mean, 1e-12),
                                  range_size$species)
  phy_helper <- .case03_phylo_helpers(phy, spp)
  trait_helper <- .case03_trait_helpers(traits, spp)
  cells <- split(x, interaction(x$cell_id, x$time_ma, drop = TRUE))
  rows <- lapply(cells, function(z) {
    z <- z[match(spp, z$species), , drop = FALSE]
    p_mean <- stats::setNames(z$posterior_mean, z$species)
    p_q025 <- stats::setNames(z$posterior_q025, z$species)
    p_q975 <- stats::setNames(z$posterior_q975, z$species)
    bin_mean <- stats::setNames(z$binary_suitability > 0, z$species)
    bin_low <- stats::setNames(z$posterior_q025 >= z$threshold, z$species)
    bin_high <- stats::setNames(z$posterior_q975 >= z$threshold, z$species)
    accessible <- if ("accessible_suitability" %in% names(z)) {
      stats::setNames(z$accessible_suitability > 0, z$species)
    } else {
      stats::setNames(rep(NA, nrow(z)), z$species)
    }
    env_suitable <- if ("environment_suitable" %in% names(z)) {
      stats::setNames(z$environment_suitable > 0, z$species)
    } else {
      bin_mean
    }
    hist_accessible <- if ("historically_accessible" %in% names(z)) {
      stats::setNames(z$historically_accessible > 0, z$species)
    } else {
      stats::setNames(rep(NA, nrow(z)), z$species)
    }
    inaccessible_env <- if ("environment_suitable_but_inaccessible" %in%
        names(z)) {
      stats::setNames(z$environment_suitable_but_inaccessible > 0, z$species)
    } else {
      stats::setNames(rep(NA, nrow(z)), z$species)
    }
    mean_metrics <- .case03_diversity_from_prob(p_mean, range_weight[z$species])
    low_metrics <- .case03_diversity_from_prob(p_q025, range_weight[z$species])
    high_metrics <- .case03_diversity_from_prob(p_q975, range_weight[z$species])
    sum_binary <- function(v) {
      if (all(is.na(v))) return(NA_real_)
      sum(v, na.rm = TRUE)
    }
    out <- data.frame(
      model_id = z$model_id[[1]], cell_id = z$cell_id[[1]],
      time_ma = z$time_ma[[1]], lon = z$lon[[1]], lat = z$lat[[1]],
      expected_species_richness_mean = mean_metrics[["expected_species_richness"]],
      expected_species_richness_sd = sqrt(sum(z$posterior_sd^2, na.rm = TRUE)),
      expected_species_richness_q025 = low_metrics[["expected_species_richness"]],
      expected_species_richness_q975 = high_metrics[["expected_species_richness"]],
      binary_species_richness = sum(bin_mean, na.rm = TRUE),
      environment_suitable_species_richness = sum(env_suitable, na.rm = TRUE),
      accessible_species_richness = sum_binary(accessible),
      historically_accessible_species_richness = sum_binary(hist_accessible),
      environment_suitable_but_inaccessible_richness =
        sum_binary(inaccessible_env),
      shannon_diversity_mean = mean_metrics[["shannon_diversity"]],
      shannon_diversity_q025 = low_metrics[["shannon_diversity"]],
      shannon_diversity_q975 = high_metrics[["shannon_diversity"]],
      simpson_diversity_mean = mean_metrics[["simpson_diversity"]],
      simpson_diversity_q025 = low_metrics[["simpson_diversity"]],
      simpson_diversity_q975 = high_metrics[["simpson_diversity"]],
      weighted_endemism_mean = mean_metrics[["weighted_endemism"]],
      weighted_endemism_q025 = low_metrics[["weighted_endemism"]],
      weighted_endemism_q975 = high_metrics[["weighted_endemism"]],
      phylogenetic_diversity_mean = .case03_faith_pd(bin_mean, phy_helper),
      phylogenetic_diversity_q025 = .case03_faith_pd(bin_low, phy_helper),
      phylogenetic_diversity_q975 = .case03_faith_pd(bin_high, phy_helper),
      mean_pairwise_phylogenetic_distance_mean = .case03_mpd(bin_mean, phy_helper),
      mean_pairwise_phylogenetic_distance_q025 = .case03_mpd(bin_low, phy_helper),
      mean_pairwise_phylogenetic_distance_q975 = .case03_mpd(bin_high, phy_helper),
      functional_richness_mean = .case03_fric(bin_mean, trait_helper),
      functional_richness_q025 = .case03_fric(bin_low, trait_helper),
      functional_richness_q975 = .case03_fric(bin_high, trait_helper),
      functional_dispersion_mean = .case03_fdisp(p_mean, trait_helper),
      functional_dispersion_q025 = .case03_fdisp(p_q025, trait_helper),
      functional_dispersion_q975 = .case03_fdisp(p_q975, trait_helper),
      stringsAsFactors = FALSE
    )
    out$shannon_diversity_sd <- pmax(out$shannon_diversity_q975 -
                                       out$shannon_diversity_q025, 0) / 3.92
    out$simpson_diversity_sd <- pmax(out$simpson_diversity_q975 -
                                       out$simpson_diversity_q025, 0) / 3.92
    out$weighted_endemism_sd <- pmax(out$weighted_endemism_q975 -
                                      out$weighted_endemism_q025, 0) / 3.92
    out$phylogenetic_diversity_sd <- pmax(out$phylogenetic_diversity_q975 -
                                           out$phylogenetic_diversity_q025, 0) / 3.92
    out$functional_richness_sd <- pmax(out$functional_richness_q975 -
                                        out$functional_richness_q025, 0) / 3.92
    out$functional_dispersion_sd <- pmax(out$functional_dispersion_q975 -
                                          out$functional_dispersion_q025, 0) / 3.92
    out$mean_pairwise_phylogenetic_distance_sd <-
      pmax(out$mean_pairwise_phylogenetic_distance_q975 -
             out$mean_pairwise_phylogenetic_distance_q025, 0) / 3.92
    out
  })
  out <- do.call(rbind, rows)
  mean_prob <- stats::aggregate(posterior_mean ~ species, x, mean,
                                na.rm = TRUE)
  names(mean_prob)[[2]] <- "mean_time_probability"
  wide <- stats::reshape(x[, c("cell_id", "species", "posterior_mean")],
                         idvar = "cell_id", timevar = "species",
                         direction = "wide")
  probs <- as.matrix(wide[, grep("^posterior_mean[.]", names(wide)), drop = FALSE])
  colnames(probs) <- sub("^posterior_mean[.]", "", colnames(probs))
  center <- mean_prob$mean_time_probability[match(colnames(probs),
                                                  mean_prob$species)]
  denom_center <- sqrt(sum(center^2, na.rm = TRUE))
  spatial_beta <- apply(probs, 1, function(v) {
    denom <- sqrt(sum(v^2, na.rm = TRUE)) * denom_center
    if (!is.finite(denom) || denom <= 0) return(0)
    .hee_clip01(1 - sum(v * center, na.rm = TRUE) / denom)
  })
  out$spatial_beta_diversity <- spatial_beta[match(out$cell_id, wide$cell_id)]
  cutoff <- stats::quantile(out$expected_species_richness_mean,
                            probs = hotspot_quantile, na.rm = TRUE,
                            names = FALSE, type = 7)
  out$diversity_hotspot <- as.integer(out$expected_species_richness_mean >= cutoff)
  out
}

#' Compute temporal transition maps between adjacent Case 03 time slices
#'
#' @param current Current species map layer table.
#' @param previous Previous older species map layer table. If `NULL`, transition
#'   metrics are returned as `NA`.
#' @return Cell-level temporal transition map table.
#' @export
hee_case03_temporal_transition_layers <- function(current, previous = NULL) {
  cur <- as.data.frame(current)
  .require_cols(cur, c("species", "cell_id", "time_ma", "lon", "lat",
                       "binary_suitability"), "current")
  base <- unique(cur[, c("model_id", "cell_id", "time_ma", "lon", "lat"),
                     drop = FALSE])
  if (is.null(previous)) {
    base$colonization <- NA_real_
    base$local_extinction <- NA_real_
    base$persistence <- NA_real_
    base$temporal_beta_diversity <- NA_real_
    return(base)
  }
  prev <- as.data.frame(previous)
  .require_cols(prev, c("species", "cell_id", "binary_suitability"), "previous")
  joined <- merge(
    cur[, c("species", "cell_id", "time_ma", "lon", "lat", "model_id",
            "binary_suitability"), drop = FALSE],
    prev[, c("species", "cell_id", "binary_suitability"), drop = FALSE],
    by = c("species", "cell_id"), all.x = TRUE, sort = FALSE,
    suffixes = c("_current", "_previous")
  )
  joined$binary_suitability_previous[is.na(joined$binary_suitability_previous)] <- 0
  rows <- split(joined, joined$cell_id)
  out <- lapply(rows, function(z) {
    cur_bin <- z$binary_suitability_current > 0
    prev_bin <- z$binary_suitability_previous > 0
    union_n <- sum(cur_bin | prev_bin, na.rm = TRUE)
    data.frame(
      model_id = z$model_id[[1]], cell_id = z$cell_id[[1]],
      time_ma = z$time_ma[[1]], lon = z$lon[[1]], lat = z$lat[[1]],
      colonization = sum(cur_bin & !prev_bin, na.rm = TRUE),
      local_extinction = sum(!cur_bin & prev_bin, na.rm = TRUE),
      persistence = sum(cur_bin & prev_bin, na.rm = TRUE),
      temporal_beta_diversity = if (union_n > 0) {
        1 - sum(cur_bin & prev_bin, na.rm = TRUE) / union_n
      } else 0,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}
