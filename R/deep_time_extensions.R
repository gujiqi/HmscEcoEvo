# Additional deep-time workflow helpers ------------------------------------

#' Deep-time workflow helper functions
#'
#' These functions fill out the full `hee_*` workflow described in the
#' deep-time manuscripts. Lightweight data checks, summaries, folds, tiling,
#' HMSC input construction, prediction summaries, and file outputs are computed
#' directly. Specialist analyses such as BioGeoBEARS fitting, plate
#' reconstruction, and long HMSC MCMC runs are exposed as explicit wrappers or
#' import/precomputed interfaces; HmscEcoEvo does not fabricate those results.
#'
#' @param prepared Prepared BioGeoBEARS/HMSC input object.
#' @param range Longitude range convention.
#' @param group_col Column used for grouping rows.
#' @param phylo_mask Optional species-by-time existence mask.
#' @param landmask Optional land/ocean mask table.
#' @param env_cube Deep-time environmental cube or wrapper.
#' @param grid Grid table with longitude and latitude.
#' @param value Optional value column to summarize.
#' @param area Cell area column.
#' @param group_cols Grouping columns for summaries.
#' @param tile_size Number of rows per tile.
#' @param id_col Optional id column.
#' @param tile A tile table to project.
#' @param tasks Character vector or list of tasks.
#' @param status Task status label.
#' @param inputs Files or objects to hash.
#' @param method Method name.
#' @param cor_threshold Correlation threshold.
#' @param vif_threshold Variance inflation factor threshold.
#' @param traits Trait table or logical switch, depending on function.
#' @param base_formula Base HMSC formula.
#' @param phylogeny Logical; include phylogenetic structure.
#' @param historical_exogenous Optional historical exogenous predictor names.
#' @param hmsc_input HMSC input list.
#' @param model_set Model definition table.
#' @param fitter Fitting function.
#' @param fit Fitted HMSC-like object.
#' @param XData Predictor table.
#' @param evo Optional `hmsc_ecoevo` object.
#' @param beta Species-by-predictor Beta matrix.
#' @param recipe Locked preprocessing recipe.
#' @param link Link function.
#' @param predicted Predicted values.
#' @param scale Prediction scale.
#' @param probs Quantile probabilities.
#' @param old Previous summary table.
#' @param new New summary table.
#' @param summary Summary object.
#' @param .grid Logical; expand scenario combinations.
#' @param scenarios Scenario table.
#' @param results Result list.
#' @param metric Metric column to summarize.
#' @param n Number of observations or simulations.
#' @param k Number of folds or blocks.
#' @param seed Optional random seed.
#' @param env Environmental table.
#' @param variable Variable name.
#' @param clade Clade labels.
#' @param true True simulated values.
#' @param estimated Estimated values.
#' @param noise_sd Simulation noise standard deviation.
#' @param n_species Number of species.
#' @param n_site Number of sites.
#' @param variables Predictor names.
#' @param species_origin Named species origin ages in Ma.
#' @param times Ages in Ma.
#' @param breaks Age breakpoints.
#' @param points Point table.
#' @param backend Plate reconstruction backend.
#' @param raw_points Raw coordinate table.
#' @param corrected_points Plate-corrected coordinate table.
#' @param time_col Time column.
#' @param distance Region distance matrix.
#' @param previous_projection Previous time-step projection.
#' @param cell_table Target species-cell-time table with coordinates for dynamic
#'   dispersal filtering.
#' @param dispersal_matrix Region dispersal matrix.
#' @param dispersal_scale Species-specific or scalar dispersal distance scale
#'   in kilometres when coordinates are used.
#' @param barrier Optional table or scalar passability multiplier for dynamic
#'   dispersal filtering.
#' @param time_direction Direction used to connect adjacent time slices; use
#'   `"forward"` for older-to-younger projections and `"backward"` for envelope
#'   reconstructions.
#' @param kernel Distance kernel for dynamic dispersal.
#' @param max_distance Optional maximum source-target distance in kilometres.
#' @param default Default weight.
#' @param comm Community matrix.
#' @param threshold Threshold object or numeric cutoff.
#' @param start_ma Oldest age in Ma.
#' @param end_ma Youngest age in Ma.
#' @param by Time-step size in Ma.
#' @param result Workflow result object.
#' @param file Optional output file.
#' @param dir Output directory.
#' @param prefix Output file prefix.
#' @name hee_deep_time_helpers
NULL

#' @rdname hee_deep_time_helpers
#' @param message Error message.
#' @param class Optional error class.
#' @export
hee_abort <- function(message, class = "hee_error") {
  e <- simpleError(message)
  class(e) <- c(class, class(e))
  stop(e)
}

#' @rdname hee_deep_time_helpers
#' @param expected Expected object names or ids.
#' @param observed Observed object names or ids.
#' @param context Context label.
#' @export
hee_error_data_mismatch <- function(expected, observed, context = "data") {
  miss <- setdiff(expected, observed)
  extra <- setdiff(observed, expected)
  msg <- paste0(context, " mismatch.",
                if (length(miss) > 0) paste0(" Missing: ", paste(miss, collapse = ", "), ".") else "",
                if (length(extra) > 0) paste0(" Extra: ", paste(extra, collapse = ", "), ".") else "")
  hee_abort(msg, class = "hee_error_data_mismatch")
}

#' @rdname hee_deep_time_helpers
#' @param ... Arguments passed to [hee_prepare_data()].
#' @export
hee_data <- function(...) {
  hee_prepare_data(...)
}

#' @rdname hee_deep_time_helpers
#' @param x A `hee_data` object or list of diagnostic objects.
#' @export
hee_audit_data <- function(x) {
  if (inherits(x, "hee_data")) x <- x$diagnostics
  tabs <- lapply(names(x), function(nm) {
    z <- x[[nm]]
    if (is.null(z)) return(NULL)
    data.frame(component = nm,
               ok = if (!is.null(z$ok)) isTRUE(z$ok) else NA,
               n_warning = length(z$warnings %||% character()),
               n_message = length(z$messages %||% character()),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, tabs[!vapply(tabs, is.null, logical(1))])
  class(out) <- c("hee_audit", class(out))
  out
}

#' @rdname hee_deep_time_helpers
#' @param dat A `hee_data` object.
#' @param XFormula HMSC environmental formula.
#' @param TrFormula HMSC trait formula.
#' @param spatial Logical; include spatial coordinates if present.
#' @param distr HMSC distribution.
#' @export
hee_as_hmsc_input <- function(dat, XFormula, TrFormula = NULL,
                              spatial = FALSE, distr = "probit") {
  if (!inherits(dat, "hee_data")) stop("dat must be a hee_data object.", call. = FALSE)
  out <- list(Y = dat$comm, XData = dat$env_now, XFormula = XFormula,
              TrData = dat$traits, TrFormula = TrFormula,
              phyloTree = dat$tree, distr = distr)
  if (spatial && all(c("lon", "lat") %in% names(dat$env_now))) {
    out$studyDesign <- data.frame(site = factor(rownames(dat$env_now)))
    out$sData <- as.matrix(dat$env_now[, c("lon", "lat")])
  }
  class(out) <- c("hee_hmsc_input", class(out))
  out
}

#' @rdname hee_deep_time_helpers
#' @param package Package/backend name.
#' @export
hee_backend_hmsc <- function(package = "Hmsc") {
  list(package = package, available = requireNamespace(package, quietly = TRUE),
       role = "HMSC model fitting and posterior prediction")
}

#' @rdname hee_deep_time_helpers
#' @export
hee_backend_bgb <- function(package = "BioGeoBEARS") {
  list(package = package, available = requireNamespace(package, quietly = TRUE),
       role = "historical biogeography model fitting; use import/precomputed if unavailable")
}

#' @rdname hee_deep_time_helpers
#' @export
hee_backend_plate <- function(package = "rgplates") {
  list(package = package, available = requireNamespace(package, quietly = TRUE),
       role = "plate reconstruction; HmscEcoEvo expects plate-corrected coordinates")
}

#' @rdname hee_deep_time_helpers
#' @param tip_ranges Species by region range table.
#' @param tree Optional tree.
#' @param regions Optional region names.
#' @export
hee_bgb_prepare <- function(tip_ranges, tree = NULL, regions = NULL) {
  x <- as.matrix(tip_ranges)
  regions <- regions %||% colnames(x)
  list(tip_ranges = x[, regions, drop = FALSE],
       tree = tree, regions = regions,
       range_size = rowSums(x[, regions, drop = FALSE] > 0, na.rm = TRUE))
}

#' @rdname hee_deep_time_helpers
#' @param models BioGeoBEARS-like model table/list.
#' @param accessibility Optional accessibility table.
#' @param events Optional BSM/event table.
#' @export
hee_bgb_import <- function(models = NULL, accessibility = NULL, events = NULL) {
  list(models = if (is.null(models)) data.frame() else as.data.frame(models),
       accessibility = if (is.null(accessibility)) data.frame() else as.data.frame(accessibility),
       events = if (is.null(events)) data.frame() else as.data.frame(events))
}

#' @rdname hee_deep_time_helpers
#' @param runner Optional user-supplied function that runs BioGeoBEARS models.
#' @param precomputed Optional precomputed result.
#' @export
hee_bgb_run_models <- function(prepared = NULL, models = c("DEC", "DEC+J", "DIVALIKE",
                                                           "DIVALIKE+J", "BAYAREALIKE",
                                                           "BAYAREALIKE+J"),
                               runner = NULL, precomputed = NULL, ...) {
  if (!is.null(precomputed)) return(precomputed)
  if (!is.null(runner) && is.function(runner)) return(runner(prepared = prepared, models = models, ...))
  stop("BioGeoBEARS model fitting is not run internally. Supply `runner` or `precomputed`.",
       call. = FALSE)
}

#' @rdname hee_deep_time_helpers
#' @param model_accessibility List of model accessibility tables.
#' @param weights Named model weights.
#' @details For `hee_bgb_accessibility()`, a direct data-frame input is treated
#'   as a pre-averaged accessibility layer. It must contain `accessibility`,
#'   `A_hist`, `A_BGB`, or `weight`; duplicate shared keys are rejected;
#'   `time_ma` must be finite and non-negative; and values are clipped to
#'   `[0, 1]`. These
#'   accessibility values are model-derived or imported historical
#'   biogeographic constraints, not HMSC-estimated environmental suitability.
#' @export
hee_bgb_accessibility <- function(model_accessibility, weights = NULL) {
  if (is.data.frame(model_accessibility)) {
    d <- as.data.frame(model_accessibility)
    if (!"accessibility" %in% names(d)) {
      hit <- intersect(c("A_hist", "A_BGB", "weight"), names(d))
      if (length(hit) == 0L) {
        stop("BioGeoBEARS accessibility data.frame must contain `accessibility`, ",
             "`A_hist`, `A_BGB`, or `weight`.", call. = FALSE)
      }
      d$accessibility <- d[[hit[1]]]
    }
    keys <- intersect(c("species", "region", "time_ma", "cell_id"), names(d))
    if (length(keys) == 0L && nrow(d) > 1L) {
      stop("BioGeoBEARS accessibility data.frame needs at least one key column ",
           "among species, region, time_ma, cell_id.", call. = FALSE)
    }
    if ("time_ma" %in% names(d)) {
      .hee_validate_time_values(d$time_ma, "model_accessibility$time_ma")
    }
    .hee_check_unique_keys(d, keys, "model_accessibility")
    d$accessibility <- .hee_clip01(d$accessibility)
    d$accessibility[is.na(d$accessibility)] <- 0
    return(d)
  }
  if (!is.list(model_accessibility) || length(model_accessibility) == 0) {
    return(data.frame())
  }
  nms <- names(model_accessibility)
  if (is.null(nms) || any(nms == "")) nms <- paste0("model", seq_along(model_accessibility))
  if (is.null(weights)) {
    weights <- stats::setNames(rep(1 / length(model_accessibility), length(model_accessibility)), nms)
  } else {
    weight_names <- names(weights)
    weights <- suppressWarnings(as.numeric(weights))
    names(weights) <- weight_names
    if (is.null(weight_names)) {
      if (length(weights) != length(nms)) {
        stop("Unnamed BioGeoBEARS accessibility weights must match the number of models.",
             call. = FALSE)
      }
      names(weights) <- nms
    }
    missing_weights <- setdiff(nms, names(weights))
    if (length(missing_weights) > 0L) {
      stop("BioGeoBEARS accessibility weights are missing model(s): ",
           paste(missing_weights, collapse = ", "), call. = FALSE)
    }
    weights <- weights[nms]
    if (any(!is.finite(weights) | weights < 0)) {
      stop("BioGeoBEARS accessibility weights must be finite non-negative values.",
           call. = FALSE)
    }
    sw <- sum(weights)
    if (!is.finite(sw) || sw <= 0) {
      stop("BioGeoBEARS accessibility weights must sum to a positive value.",
           call. = FALSE)
    }
    weights <- weights / sw
  }
  rows <- list()
  for (nm in nms) {
    d <- as.data.frame(model_accessibility[[nm]])
    if (!"accessibility" %in% names(d)) {
      val <- names(d)[vapply(d, is.numeric, logical(1))]
      val <- setdiff(val, c("time_ma", "lon", "lat"))
      d$accessibility <- if (length(val) > 0) d[[val[1]]] else 1
    }
    keys <- intersect(c("species", "region", "time_ma", "cell_id"), names(d))
    if (length(keys) == 0L && nrow(d) > 1L) {
      stop("BioGeoBEARS accessibility table for model ", nm,
           " needs at least one key column among species, region, time_ma, cell_id.",
           call. = FALSE)
    }
    .hee_check_unique_keys(d, keys, paste0("model_accessibility$", nm))
    d$accessibility <- .hee_clip01(d$accessibility)
    d$.weighted <- d$accessibility * as.numeric(weights[nm])
    rows[[nm]] <- d
  }
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(d) {
    miss <- setdiff(all_names, names(d))
    for (nm in miss) d[[nm]] <- NA
    d[, all_names, drop = FALSE]
  })
  all <- do.call(rbind, rows)
  keys <- intersect(c("species", "region", "time_ma", "cell_id"), names(all))
  out <- stats::aggregate(.weighted ~ ., all[, c(keys, ".weighted"), drop = FALSE], sum, na.rm = TRUE)
  names(out)[names(out) == ".weighted"] <- "accessibility"
  out$accessibility <- .hee_clip01(out$accessibility)
  out
}

#' @rdname hee_deep_time_helpers
#' @param bsm_events BSM/event table.
#' @export
hee_bgb_bsm <- function(bsm_events) {
  x <- as.data.frame(bsm_events)
  if ("time" %in% names(x) && !"time_ma" %in% names(x)) names(x)[names(x) == "time"] <- "time_ma"
  x
}

#' @rdname hee_deep_time_helpers
#' @param events Standardized event table.
#' @export
hee_bgb_event_summary <- function(events) {
  ev <- hee_bgb_bsm(events)
  if (nrow(ev) == 0) return(list(region = data.frame(), dispersal = data.frame()))
  if (!"event_prob" %in% names(ev)) ev$event_prob <- 1
  region <- data.frame()
  if ("region" %in% names(ev)) {
    region <- stats::aggregate(event_prob ~ region + time_ma + event_type, ev, sum, na.rm = TRUE)
    names(region)[4] <- "event_weight"
  }
  dispersal <- data.frame()
  if (all(c("from_region", "to_region") %in% names(ev))) {
    d <- ev[ev$event_type %in% c("dispersal", "range expansion", "founder-event"), , drop = FALSE]
    if (nrow(d) > 0) {
      dispersal <- stats::aggregate(event_prob ~ from_region + to_region + time_ma, d, sum, na.rm = TRUE)
      names(dispersal)[4] <- "dispersal_rate"
    }
  }
  list(region = region, dispersal = dispersal)
}

#' @rdname hee_deep_time_helpers
#' @param projection Projection table.
#' @export
hee_can_compute_richness <- function(projection) {
  all(c("cell_id", "time_ma", "probability") %in% names(as.data.frame(projection)))
}

#' @rdname hee_deep_time_helpers
#' @param lon Longitude centers.
#' @param lat Latitude centers.
#' @param resolution Resolution in degrees.
#' @export
hee_cell_area <- function(lon, lat, resolution = 1) {
  R <- 6371.0088
  dlon <- resolution * pi / 180
  lat1 <- (lat - resolution / 2) * pi / 180
  lat2 <- (lat + resolution / 2) * pi / 180
  abs(R^2 * dlon * (sin(lat2) - sin(lat1)))
}

#' @rdname hee_deep_time_helpers
#' @param coords Coordinate table.
#' @param lon_col Longitude column.
#' @param lat_col Latitude column.
#' @export
hee_normalize_coordinates <- function(coords, lon_col = "lon", lat_col = "lat",
                                      range = c("-180_180", "0_360")) {
  range <- match.arg(range)
  x <- as.data.frame(coords)
  .require_cols(x, c(lon_col, lat_col), "coords")
  x[[lon_col]] <- if (range == "-180_180") ((x[[lon_col]] + 180) %% 360) - 180 else x[[lon_col]] %% 360
  x[[lat_col]] <- pmax(pmin(x[[lat_col]], 90), -90)
  x
}

#' @rdname hee_deep_time_helpers
#' @param lon Numeric longitudes.
#' @export
hee_convert_lon_range <- function(lon, range = c("-180_180", "0_360")) {
  range <- match.arg(range)
  if (range == "-180_180") ((lon + 180) %% 360) - 180 else lon %% 360
}

#' @rdname hee_deep_time_helpers
#' @export
hee_check_antimeridian <- function(coords, lon_col = "lon", group_col = NULL) {
  x <- as.data.frame(coords)
  groups <- if (is.null(group_col)) list(all = x) else split(x, x[[group_col]])
  out <- lapply(names(groups), function(nm) {
    z <- groups[[nm]]
    jump <- if (nrow(z) < 2) FALSE else any(abs(diff(z[[lon_col]])) > 180, na.rm = TRUE)
    data.frame(group = nm, crosses_antimeridian = jump, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

#' @rdname hee_deep_time_helpers
#' @param obj Object with CRS metadata or a CRS string.
#' @param expected Expected CRS string.
#' @export
hee_check_crs <- function(obj, expected = "EPSG:4326") {
  crs <- if (is.character(obj)) obj else attr(obj, "crs") %||% NA_character_
  list(ok = is.na(crs) || identical(crs, expected), crs = crs, expected = expected)
}

#' @rdname hee_deep_time_helpers
#' @param train Training data.
#' @param newdata New data.
#' @export
hee_check_factor_levels <- function(train, newdata) {
  train <- as.data.frame(train); newdata <- as.data.frame(newdata)
  fac <- names(train)[vapply(train, function(z) is.factor(z) || is.character(z), logical(1))]
  out <- lapply(fac, function(v) {
    unseen <- setdiff(unique(as.character(newdata[[v]])), unique(as.character(train[[v]])))
    data.frame(variable = v, n_unseen = length(unseen),
               unseen = paste(unseen, collapse = ";"), stringsAsFactors = FALSE)
  })
  if (length(out) == 0) data.frame() else do.call(rbind, out)
}

#' @rdname hee_deep_time_helpers
#' @param suitability Suitability table.
#' @export
hee_check_combine_inputs <- function(suitability, accessibility = NULL,
                                     phylo_mask = NULL, landmask = NULL) {
  P <- as.data.frame(suitability)
  missing <- setdiff(c("species", "cell_id", "time_ma", "suitability"), names(P))
  list(ok = length(missing) == 0,
       missing_suitability_columns = missing,
       has_accessibility = !is.null(accessibility),
       has_phylo_mask = !is.null(phylo_mask),
       has_landmask = !is.null(landmask))
}

#' @rdname hee_deep_time_helpers
#' @param a First table.
#' @param b Second table.
#' @param keys Alignment keys.
#' @export
hee_align_by_cell_id <- function(a, b, keys = c("cell_id", "time_ma")) {
  merge(as.data.frame(a), as.data.frame(b), by = keys, all = FALSE, sort = FALSE)
}

#' @rdname hee_deep_time_helpers
#' @export
hee_check_projection_alignment <- function(a, b, keys = c("cell_id", "time_ma")) {
  A <- unique(as.data.frame(a)[, keys, drop = FALSE])
  B <- unique(as.data.frame(b)[, keys, drop = FALSE])
  both <- merge(A, B, by = keys)
  list(ok = nrow(both) == nrow(A) && nrow(both) == nrow(B),
       n_a = nrow(A), n_b = nrow(B), n_overlap = nrow(both))
}

#' @rdname hee_deep_time_helpers
#' @export
hee_make_cell_table <- function(env_cube = NULL, grid = NULL, resolution = NULL) {
  if (!is.null(grid)) {
    d <- as.data.frame(grid)
  } else {
    cube <- .hee_as_timecube(env_cube)
    d <- expand.grid(lon = cube$cube$lon, lat = cube$cube$lat)
  }
  d$cell_id <- .hee_cell_id(d$lon, d$lat)
  if (!is.null(resolution)) d$cell_area_km2 <- hee_cell_area(d$lon, d$lat, resolution)
  d
}

#' @rdname hee_deep_time_helpers
#' @export
hee_area_summary <- function(x, value = NULL, area = "land_area_km2",
                             group_cols = "time_ma") {
  d <- as.data.frame(x)
  if (!area %in% names(d)) d[[area]] <- 1
  if (is.null(value)) {
    stats::aggregate(d[[area]], d[group_cols], sum, na.rm = TRUE)
  } else {
    .require_cols(d, c(value, area, group_cols), "x")
    d$.weighted <- d[[value]] * d[[area]]
    num <- stats::aggregate(.weighted ~ ., d[, c(group_cols, ".weighted"), drop = FALSE], sum, na.rm = TRUE)
    den <- stats::aggregate(d[[area]], d[group_cols], sum, na.rm = TRUE)
    names(den)[ncol(den)] <- ".area"
    out <- merge(num, den, by = group_cols)
    out$weighted_mean <- out$.weighted / out$.area
    out
  }
}

#' @rdname hee_deep_time_helpers
#' @export
hee_make_tiles <- function(x, tile_size = 5000, id_col = NULL) {
  n <- if (is.data.frame(x)) nrow(x) else length(x)
  tile <- ceiling(seq_len(n) / tile_size)
  data.frame(row_id = seq_len(n), tile_id = tile,
             id = if (!is.null(id_col) && is.data.frame(x)) x[[id_col]] else seq_len(n))
}

#' @rdname hee_deep_time_helpers
#' @param tiles List of tile result tables.
#' @export
hee_merge_tiles <- function(tiles) {
  if (length(tiles) == 0) data.frame() else do.call(rbind, lapply(tiles, as.data.frame))
}

#' @rdname hee_deep_time_helpers
#' @export
hee_project_hmsc_tile <- function(tile, ...) {
  hee_project_hmsc_table(X_past = tile, ...)
}

#' @rdname hee_deep_time_helpers
#' @export
hee_task_table <- function(tasks, status = "pending") {
  data.frame(task_id = seq_along(tasks), task = as.character(tasks),
             status = status, stringsAsFactors = FALSE)
}

#' @rdname hee_deep_time_helpers
#' @param task A function or expression-like object.
#' @export
hee_run_task <- function(task, ...) {
  if (is.function(task)) return(task(...))
  stop("task must be a function.", call. = FALSE)
}

#' @rdname hee_deep_time_helpers
#' @param task_table Task table.
#' @export
hee_resume <- function(task_table) {
  x <- as.data.frame(task_table)
  x[x$status != "done", , drop = FALSE]
}

#' @rdname hee_deep_time_helpers
#' @param path Project path.
#' @param overwrite Overwrite existing directory.
#' @export
hee_init_project <- function(path, overwrite = FALSE) {
  if (dir.exists(path) && !overwrite) stop("Project path already exists: ", path, call. = FALSE)
  dirs <- c("data_raw", "data_processed", "env_cube", "plate_reconstruction",
            "biogeobears", "hmsc_models", "projections", "rasters",
            "summaries", "diagnostics", "reports", "logs", "config")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  for (d in dirs) dir.create(file.path(path, d), recursive = TRUE, showWarnings = FALSE)
  writeLines(c("parameter,value", "created_at,auto"), file.path(path, "config", "params.csv"))
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(path, "logs", "session_info.txt"))
  invisible(normalizePath(path))
}

#' @rdname hee_deep_time_helpers
#' @export
hee_hash_inputs <- function(inputs) {
  if (is.character(inputs) && all(file.exists(inputs))) {
    return(data.frame(input = inputs, md5 = unname(tools::md5sum(inputs)), stringsAsFactors = FALSE))
  }
  raw <- serialize(inputs, NULL)
  data.frame(input = "object", md5 = paste(as.character(sum(as.integer(raw))), nchar(raw), sep = "-"),
             stringsAsFactors = FALSE)
}

#' @rdname hee_deep_time_helpers
#' @export
hee_filter_collinearity <- function(x, method = c("correlation", "vif"),
                                    cor_threshold = 0.7, vif_threshold = 5) {
  method <- match.arg(method)
  d <- as.data.frame(x)
  num <- d[vapply(d, is.numeric, logical(1))]
  if (ncol(num) < 2) return(list(keep = names(num), drop = character(), diagnostics = data.frame()))
  if (method == "correlation") {
    C <- abs(stats::cor(num, use = "pairwise.complete.obs"))
    drop <- character()
    vars <- colnames(C)
    for (i in seq_along(vars)) {
      if (vars[i] %in% drop) next
      high <- vars[which(C[vars[i], ] > cor_threshold)]
      high <- setdiff(high, c(vars[i], drop))
      drop <- c(drop, high)
    }
    return(list(keep = setdiff(names(num), drop), drop = unique(drop), correlation = C))
  }
  vif <- vapply(names(num), function(v) {
    others <- setdiff(names(num), v)
    r2 <- .safe_lm_r2(num[[v]], num[, others, drop = FALSE])
    1 / max(1 - r2, .Machine$double.eps)
  }, numeric(1))
  list(keep = names(vif)[vif <= vif_threshold], drop = names(vif)[vif > vif_threshold], vif = vif)
}

#' @rdname hee_deep_time_helpers
#' @export
hee_impute_traits <- function(traits, method = c("mean", "mode")) {
  method <- match.arg(method)
  x <- as.data.frame(traits)
  flags <- data.frame(row.names = rownames(x))
  for (v in names(x)) {
    miss <- is.na(x[[v]])
    flags[[paste0(v, "_missing")]] <- miss
    if (!any(miss)) next
    if (is.numeric(x[[v]])) {
      x[[v]][miss] <- mean(x[[v]], na.rm = TRUE)
    } else {
      tab <- sort(table(x[[v]]), decreasing = TRUE)
      x[[v]][miss] <- names(tab)[1]
    }
  }
  list(traits = x, missing_flags = flags)
}

#' @rdname hee_deep_time_helpers
#' @export
hee_define_hmsc_models <- function(base_formula, traits = TRUE, phylogeny = TRUE,
                                   spatial = TRUE, historical_exogenous = NULL) {
  data.frame(
    model_id = c("M0", "M1", "M2", "M3", "M4", "M5", "M6"),
    environment = TRUE,
    traits = c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE) & traits,
    phylogeny = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE) & phylogeny,
    latent_factors = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE),
    spatial = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE) & spatial,
    historical_exogenous = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE) &
      !is.null(historical_exogenous),
    formula = paste(deparse(base_formula), collapse = " "),
    stringsAsFactors = FALSE
  )
}

#' @rdname hee_deep_time_helpers
#' @export
hee_fit_hmsc <- function(hmsc_input = NULL, precomputed = NULL, ...) {
  if (!is.null(precomputed)) return(precomputed)
  if (!requireNamespace("Hmsc", quietly = TRUE)) {
    stop("Package 'Hmsc' is required for model fitting. Supply `precomputed` to avoid long MCMC.",
         call. = FALSE)
  }
  do.call(Hmsc::Hmsc, c(hmsc_input, list(...)))
}

#' @rdname hee_deep_time_helpers
#' @export
hee_fit_hmsc_set <- function(hmsc_input, model_set, fitter = hee_fit_hmsc, ...) {
  ids <- model_set$model_id %||% paste0("M", seq_len(nrow(model_set)))
  out <- lapply(ids, function(id) fitter(hmsc_input = hmsc_input, ...))
  names(out) <- ids
  out
}

#' @rdname hee_deep_time_helpers
#' @export
hee_hmsc_predict_adapter <- function(fit, XData = NULL, expected = NULL, ...) {
  if (!is.null(expected)) return(expected)
  if (inherits(fit, "hmsc_ecoevo")) return(NULL)
  if (requireNamespace("Hmsc", quietly = TRUE)) {
    pred <- tryCatch(Hmsc::computePredictedValues(fit, ...), error = function(e) NULL)
    if (!is.null(pred)) return(pred)
  }
  stop("Could not obtain HMSC predictions. Supply `expected` or a supported Hmsc fit.",
       call. = FALSE)
}

#' @rdname hee_deep_time_helpers
#' @export
hee_predict_hmsc <- function(...) hee_hmsc_predict_adapter(...)

#' @rdname hee_deep_time_helpers
#' @export
hee_predict_modern <- function(evo = NULL, beta = NULL, XData, recipe = NULL,
                               link = "probit") {
  hee_project_hmsc_table(evo = evo, beta = beta, X_past = XData,
                         recipe = recipe, link = link)
}

#' @rdname hee_deep_time_helpers
#' @export
hee_validate_modern <- function(observed, predicted) {
  hee_evaluate_hmsc(observed, predicted)
}

#' @rdname hee_deep_time_helpers
#' @export
hee_evaluate_hmsc <- function(observed, predicted) {
  obs <- as.numeric(observed); pred <- as.numeric(predicted)
  ok <- !is.na(obs) & !is.na(pred)
  obs <- obs[ok]; pred <- pred[ok]
  auc <- .hee_auc(obs, pred)
  data.frame(AUC = auc,
             Tjur_R2 = mean(pred[obs > 0], na.rm = TRUE) - mean(pred[obs <= 0], na.rm = TRUE),
             RMSE = sqrt(mean((obs - pred)^2, na.rm = TRUE)),
             calibration_slope = .hee_calibration_slope(obs, pred),
             stringsAsFactors = FALSE)
}

#' @rdname hee_deep_time_helpers
#' @export
hee_mcmc_diagnostics <- function(evo = NULL, precomputed = NULL) {
  if (!is.null(precomputed)) return(as.data.frame(precomputed))
  if (is.null(evo)) return(data.frame())
  calc_validation_metrics(evo)$mcmc_diagnostics
}

#' @rdname hee_deep_time_helpers
#' @export
hee_prediction_scale <- function(x, scale = c("response", "linear", "logit")) {
  scale <- match.arg(scale)
  if (scale == "response") return(pmin(pmax(x, 0), 1))
  if (scale == "logit") return(stats::qlogis(pmin(pmax(x, 1e-6), 1 - 1e-6)))
  x
}

#' @rdname hee_deep_time_helpers
#' @export
hee_summarise_prediction_array <- function(x, probs = c(0.025, 0.5, 0.975)) {
  if (is.data.frame(x)) {
    val <- names(x)[vapply(x, is.numeric, logical(1))]
    val <- setdiff(val, c("time_ma", "lon", "lat"))
    return(summary(x[[val[1]]]))
  }
  apply(x, seq_len(length(dim(x)))[-1], function(z) {
    c(mean = mean(z, na.rm = TRUE), stats::quantile(z, probs, na.rm = TRUE))
  })
}

#' @rdname hee_deep_time_helpers
#' @export
hee_online_summary <- function(x, value = "probability", group_cols = c("species", "time_ma")) {
  d <- as.data.frame(x)
  stats::aggregate(d[[value]], d[group_cols], function(z) {
    c(mean = mean(z, na.rm = TRUE), sd = stats::sd(z, na.rm = TRUE), n = sum(!is.na(z)))
  })
}

#' @rdname hee_deep_time_helpers
#' @export
hee_update_summary <- function(old, new) {
  if (is.null(old) || nrow(as.data.frame(old)) == 0) return(new)
  rbind(as.data.frame(old), as.data.frame(new))
}

#' @rdname hee_deep_time_helpers
#' @export
hee_finalize_summary <- function(summary) {
  as.data.frame(summary)
}

#' @rdname hee_deep_time_helpers
#' @export
hee_define_scenarios <- function(..., .grid = TRUE) {
  vals <- list(...)
  if (.grid) do.call(expand.grid, c(vals, stringsAsFactors = FALSE)) else as.data.frame(vals)
}

#' @rdname hee_deep_time_helpers
#' @export
hee_run_scenarios <- function(scenarios, runner, ...) {
  sc <- as.data.frame(scenarios)
  out <- lapply(seq_len(nrow(sc)), function(i) runner(sc[i, , drop = FALSE], ...))
  names(out) <- paste0("scenario", seq_along(out))
  out
}

#' @rdname hee_deep_time_helpers
#' @export
hee_summarise_scenarios <- function(results, metric = "expected_richness") {
  rows <- lapply(names(results), function(nm) {
    z <- results[[nm]]
    if (is.list(z) && !is.null(z$richness)) z <- z$richness
    z <- as.data.frame(z)
    val <- if (metric %in% names(z)) z[[metric]] else z[[names(z)[vapply(z, is.numeric, logical(1))][1]]]
    data.frame(scenario = nm, mean = mean(val, na.rm = TRUE), sd = stats::sd(val, na.rm = TRUE))
  })
  do.call(rbind, rows)
}

#' @rdname hee_deep_time_helpers
#' @export
hee_uncertainty_budget <- function(results, metric = "expected_richness") {
  s <- hee_summarise_scenarios(results, metric = metric)
  total <- stats::var(s$mean, na.rm = TRUE)
  data.frame(metric = metric, scenario_variance = total,
             mean_within_scenario_sd = mean(s$sd, na.rm = TRUE))
}

#' @rdname hee_deep_time_helpers
#' @export
hee_cv_random <- function(n, k = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sample(rep(seq_len(k), length.out = n))
}

#' @rdname hee_deep_time_helpers
#' @export
hee_cv_spatial_block <- function(coords, k = 5) {
  x <- as.data.frame(coords)
  .require_cols(x, c("lon", "lat"), "coords")
  cut(x$lon, breaks = k, labels = FALSE, include.lowest = TRUE)
}

#' @rdname hee_deep_time_helpers
#' @export
hee_cv_environment_block <- function(env, variable = NULL, k = 5) {
  x <- as.data.frame(env)
  variable <- variable %||% names(x)[vapply(x, is.numeric, logical(1))][1]
  cut(x[[variable]], breaks = k, labels = FALSE, include.lowest = TRUE)
}

#' @rdname hee_deep_time_helpers
#' @export
hee_cv_leave_clade_out <- function(clade) {
  as.integer(factor(clade))
}

#' @rdname hee_deep_time_helpers
#' @export
hee_evaluate_recovery <- function(true, estimated) {
  true <- as.numeric(true); estimated <- as.numeric(estimated)
  ok <- !is.na(true) & !is.na(estimated)
  data.frame(R2 = if (sum(ok) > 2) stats::cor(true[ok], estimated[ok])^2 else NA_real_,
             bias = mean(estimated[ok] - true[ok]),
             RMSE = sqrt(mean((estimated[ok] - true[ok])^2)))
}

#' @rdname hee_deep_time_helpers
#' @export
hee_simulate_recovery <- function(n = 100, noise_sd = 0.1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  true <- stats::rnorm(n)
  estimated <- true + stats::rnorm(n, 0, noise_sd)
  cbind(data.frame(true = true, estimated = estimated), hee_evaluate_recovery(true, estimated))
}

#' @rdname hee_deep_time_helpers
#' @export
hee_simulate_example_data <- function(n_species = 20, n_site = 50,
                                      variables = c("temp", "precip")) {
  species <- paste0("sp", seq_len(n_species))
  sites <- paste0("site", seq_len(n_site))
  env <- as.data.frame(matrix(stats::rnorm(n_site * length(variables)), n_site,
                              dimnames = list(sites, variables)))
  beta <- matrix(stats::rnorm(n_species * (length(variables) + 1), 0, 0.5),
                 n_species, dimnames = list(species, c("(Intercept)", variables)))
  eta <- as.matrix(cbind("(Intercept)" = 1, env[, variables, drop = FALSE])) %*% t(beta)
  prob <- stats::pnorm(eta)
  comm <- matrix(stats::rbinom(length(prob), 1, prob), n_site,
                 dimnames = list(sites, species))
  list(comm = comm, env_now = env, beta = beta,
       species_origin = stats::setNames(seq(5, 100, length.out = n_species), species))
}

#' @rdname hee_deep_time_helpers
#' @export
hee_estimate_species_origin <- function(tree = NULL, species_origin = NULL) {
  if (!is.null(species_origin)) {
    if (is.data.frame(species_origin)) return(stats::setNames(species_origin$origin_ma, species_origin$species))
    return(species_origin)
  }
  .require_pkg("ape", "species origin estimation")
  if (is.null(tree) || !inherits(tree, "phylo")) stop("Supply a dated ape::phylo tree.", call. = FALSE)
  depth <- ape::node.depth.edgelength(tree)
  root_age <- max(depth)
  parent <- tree$edge[match(seq_along(tree$tip.label), tree$edge[, 2]), 1]
  origin <- root_age - depth[parent]
  stats::setNames(origin, tree$tip.label)
}

#' @rdname hee_deep_time_helpers
#' @export
hee_lineage_mode <- function(times, breaks = c(2, 20, 50)) {
  cut(times, breaks = c(-Inf, breaks, Inf),
      labels = c("species", "species_or_lineage", "lineage_or_clade", "clade"),
      right = TRUE)
}

#' @rdname hee_deep_time_helpers
#' @export
hee_reconstruct_points <- function(points, precomputed = NULL,
                                   lon_col = "lon", lat_col = "lat",
                                   backend = c("user_table", "rgplates", "pygplates"),
                                   ...) {
  backend <- match.arg(backend)
  if (!is.null(precomputed)) return(as.data.frame(precomputed))
  if (backend != "user_table") {
    stop("Plate reconstruction backend '", backend,
         "' is not run internally in this version. Supply `precomputed` or use backend = 'user_table'.",
         call. = FALSE)
  }
  x <- as.data.frame(points)
  if (all(c("paleo_lon", "paleo_lat", "time_ma") %in% names(x))) return(x)
  stop("Plate reconstruction is not performed internally. Supply plate-corrected `precomputed` points.",
       call. = FALSE)
}

#' @rdname hee_deep_time_helpers
#' @export
hee_import_plate_points <- function(points) {
  x <- as.data.frame(points)
  need <- c("time_ma", "paleo_lon", "paleo_lat")
  if (all(need %in% names(x))) {
    if (!"coord_status" %in% names(x)) x$coord_status <- "plate_reconstructed"
    return(x)
  }
  if (all(c("lon", "lat") %in% names(x))) {
    if (!"time_ma" %in% names(x) || all(x$time_ma == 0, na.rm = TRUE)) {
      x$time_ma <- x$time_ma %||% 0
      x$paleo_lon <- x$lon
      x$paleo_lat <- x$lat
      x$coord_status <- "modern_coordinates_0Ma"
      return(x)
    }
    stop("Plate points for time_ma > 0 require plate-reconstructed `paleo_lon` ",
         "and `paleo_lat`. HmscEcoEvo will not silently use modern lon/lat as ",
         "deep-time coordinates.", call. = FALSE)
  }
  missing <- setdiff(need, names(x))
  stop("plate points are missing required column(s): ",
       paste(missing, collapse = ", "), call. = FALSE)
}

#' @rdname hee_deep_time_helpers
#' @export
hee_compare_plate_corrected <- function(raw_points, corrected_points) {
  a <- as.data.frame(raw_points); b <- as.data.frame(corrected_points)
  if (!"point_id" %in% names(a)) a$point_id <- seq_len(nrow(a))
  if (!"point_id" %in% names(b)) b$point_id <- seq_len(nrow(b))
  .hee_check_unique_keys(a, "point_id", "raw_points")
  .hee_check_unique_keys(b, "point_id", "corrected_points")
  m <- merge(a, b, by = "point_id", suffixes = c("_raw", "_corrected"))
  if (nrow(m) != nrow(a) || nrow(m) != nrow(b)) {
    stop("Plate-correction comparison must be one-to-one by point_id.",
         call. = FALSE)
  }
  lon_raw <- if ("lon_raw" %in% names(m)) m$lon_raw else m$lon
  lat_raw <- if ("lat_raw" %in% names(m)) m$lat_raw else m$lat
  paleo_lon <- if ("paleo_lon" %in% names(m)) m$paleo_lon else m$paleo_lon_corrected
  paleo_lat <- if ("paleo_lat" %in% names(m)) m$paleo_lat else m$paleo_lat_corrected
  m$shift_deg <- sqrt((lon_raw - paleo_lon)^2 + (lat_raw - paleo_lat)^2)
  m$shift_km <- hee_great_circle_distance_km(lon_raw, lat_raw,
                                              paleo_lon, paleo_lat)
  m
}

#' @rdname hee_deep_time_helpers
#' @export
hee_centroid_shift <- function(points, group_col = "track_id",
                               lon_col = "paleo_lon", lat_col = "paleo_lat",
                               time_col = "time_ma") {
  x <- as.data.frame(points)
  x <- x[order(x[[group_col]], x[[time_col]], decreasing = TRUE), ]
  out <- lapply(split(x, x[[group_col]]), function(z) {
    data.frame(group = z[[group_col]][1],
               total_shift_deg = sum(sqrt(diff(z[[lon_col]])^2 + diff(z[[lat_col]])^2), na.rm = TRUE),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

#' @rdname hee_deep_time_helpers
#' @export
hee_make_dispersal_matrices <- function(regions, distance = NULL, scale = 1) {
  scale <- .validate_dispersal_scale(scale, name = "scale")
  if (length(scale) != 1L) {
    stop("scale must be a single finite non-negative value.", call. = FALSE)
  }
  regions <- as.character(regions)
  if (is.null(distance)) {
    distance <- as.matrix(stats::dist(seq_along(regions)))
    dimnames(distance) <- list(regions, regions)
  }
  distance <- as.matrix(distance)
  if (any(is.finite(distance) & distance < 0, na.rm = TRUE)) {
    stop("distance must contain non-negative values.", call. = FALSE)
  }
  distance[is.na(distance)] <- Inf
  if (scale <= 0) {
    W <- ifelse(distance <= 0, 1, 0)
  } else {
    W <- exp(-distance / scale)
  }
  W[!is.finite(W)] <- 0
  diag(W) <- 1
  W
}

#' @rdname hee_deep_time_helpers
#' @export
hee_dynamic_filter <- function(previous_projection = NULL,
                               cell_table = NULL,
                               dispersal_matrix = NULL,
                               dispersal_scale = 1,
                               barrier = NULL,
                               time_direction = c("forward", "backward"),
                               kernel = c("exponential", "gaussian"),
                               max_distance = NULL,
                               default = 1) {
  time_direction <- match.arg(time_direction)
  kernel <- match.arg(kernel)
  dispersal_scale <- .validate_dispersal_scale(dispersal_scale)
  if (is.null(previous_projection)) return(data.frame(dynamic_weight = numeric()))
  P <- as.data.frame(previous_projection)
  if (!"probability" %in% names(P)) {
    out <- P[, intersect(c("species", "cell_id", "time_ma", "region"), names(P)), drop = FALSE]
    out$dynamic_weight <- default
    out$D_dynamic <- out$dynamic_weight
    out$dynamic_direction <- time_direction
    return(out)
  }
  .require_cols(P, c("species", "cell_id", "time_ma", "probability"), "previous_projection")
  if (is.null(cell_table)) {
    out <- P[, intersect(c("species", "cell_id", "time_ma", "region", "lon", "lat"), names(P)), drop = FALSE]
    out$dynamic_weight <- pmin(pmax(P$probability, 0), 1)
    out$D_dynamic <- out$dynamic_weight
    out$dynamic_direction <- time_direction
    return(out)
  }
  targets <- as.data.frame(cell_table)
  .require_cols(targets, c("cell_id", "time_ma", "lon", "lat"), "cell_table")
  species <- unique(P$species)
  if (!"species" %in% names(targets)) {
    targets <- merge(data.frame(species = species, stringsAsFactors = FALSE),
                     targets, all = TRUE, sort = FALSE)
  }
  targets$.row_id <- seq_len(nrow(targets))
  times <- sort(unique(c(P$time_ma, targets$time_ma)), decreasing = TRUE)
  scale_for <- function(sp) {
    if (length(dispersal_scale) == 1L) return(as.numeric(dispersal_scale))
    val <- dispersal_scale[as.character(sp)]
    ifelse(is.na(val), stats::median(dispersal_scale, na.rm = TRUE), as.numeric(val))
  }
  kernel_fun <- function(distance, scale) {
    if (kernel == "gaussian") exp(-(distance / scale)^2) else exp(-distance / scale)
  }
  parts <- lapply(split(targets, interaction(targets$species, targets$time_ma, drop = TRUE)), function(tg) {
    sp <- tg$species[1]
    tt <- tg$time_ma[1]
    prev_time <- if (time_direction == "forward") {
      older <- times[times > tt]
      if (length(older) == 0L) NA_real_ else min(older)
    } else {
      younger <- times[times < tt]
      if (length(younger) == 0L) NA_real_ else max(younger)
    }
    src <- if (is.na(prev_time)) data.frame() else P[P$species == sp & P$time_ma == prev_time, , drop = FALSE]
    if (is.na(prev_time)) {
      tg$dynamic_weight <- default
      tg$dynamic_source_time_ma <- prev_time
      tg$dynamic_min_distance <- NA_real_
      return(tg)
    }
    if (nrow(src) == 0L) {
      tg$dynamic_weight <- 0
      tg$dynamic_source_time_ma <- prev_time
      tg$dynamic_min_distance <- NA_real_
      return(tg)
    }
    if (!all(c("lon", "lat") %in% names(src))) {
      stop("previous_projection must contain lon/lat for dynamic distance ",
           "filtering when a previous time slice exists.", call. = FALSE)
    }
    sc <- max(scale_for(sp), .Machine$double.eps)
    dist <- .hee_great_circle_matrix(tg$lon, tg$lat, src$lon, src$lat)
    if (!is.null(max_distance)) dist[dist > max_distance] <- Inf
    K <- kernel_fun(dist, sc)
    value <- sweep(K, 2, pmin(pmax(src$probability, 0), 1), `*`)
    dyn <- apply(value, 1, max, na.rm = TRUE)
    mind <- apply(dist, 1, min, na.rm = TRUE)
    dyn[!is.finite(dyn)] <- 0
    mind[!is.finite(mind)] <- NA_real_
    tg$dynamic_weight <- dyn
    tg$dynamic_source_time_ma <- prev_time
    tg$dynamic_min_distance <- mind
    tg
  })
  out <- do.call(rbind, parts)
  if (!is.null(barrier)) {
    b <- as.data.frame(barrier)
    if (!"barrier_passability" %in% names(b)) b$barrier_passability <- 1
    keys <- intersect(c("species", "cell_id", "time_ma"), names(b))
    if (!"species" %in% keys) keys <- intersect(c("cell_id", "time_ma"), names(b))
    if (length(keys) == 0L) {
      if (nrow(b) != 1L) {
        stop("barrier must contain join keys or exactly one global row.",
             call. = FALSE)
      }
      out$barrier_passability <- .hee_clip01(b$barrier_passability[1])
      out$barrier_passability[is.na(out$barrier_passability)] <- 0
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
    }
    out$dynamic_weight <- out$dynamic_weight * out$barrier_passability
  }
  out$dynamic_weight <- pmin(pmax(out$dynamic_weight, 0), 1)
  out$D_dynamic <- out$dynamic_weight
  out$dynamic_direction <- time_direction
  out <- out[order(out$.row_id), , drop = FALSE]
  out[, setdiff(names(out), ".row_id"), drop = FALSE]
}

#' @rdname hee_deep_time_helpers
#' @export
hee_cwm <- function(comm, traits) {
  Y <- as.matrix(comm)
  T <- as.matrix(traits[colnames(Y), , drop = FALSE])
  rs <- rowSums(Y, na.rm = TRUE)
  rs[rs == 0] <- NA_real_
  out <- Y %*% T / rs
  as.data.frame(out)
}

#' @rdname hee_deep_time_helpers
#' @export
hee_functional_diversity <- function(comm, traits) {
  CWM <- hee_cwm(comm, traits)
  data.frame(site = rownames(CWM),
             functional_dispersion = apply(CWM, 1, stats::sd, na.rm = TRUE),
             stringsAsFactors = FALSE)
}

#' @rdname hee_deep_time_helpers
#' @export
hee_phylo_diversity <- function(comm, tree = NULL) {
  Y <- as.matrix(comm)
  if (is.null(tree) || !inherits(tree, "phylo")) {
    return(data.frame(site = rownames(Y), phylo_diversity = rowSums(Y > 0, na.rm = TRUE)))
  }
  .require_pkg("ape", "phylogenetic diversity")
  pd <- vapply(seq_len(nrow(Y)), function(i) {
    spp <- colnames(Y)[Y[i, ] > 0]
    if (length(spp) < 2) return(0)
    tr <- ape::keep.tip(tree, intersect(spp, tree$tip.label))
    sum(tr$edge.length, na.rm = TRUE)
  }, numeric(1))
  data.frame(site = rownames(Y), phylo_diversity = pd)
}

#' @rdname hee_deep_time_helpers
#' @export
hee_threshold <- function(observed = NULL, predicted = NULL, threshold = NULL, ...) {
  if (!is.null(threshold) && is.null(observed)) return(hee_threshold_apply(threshold, predicted))
  hee_threshold_train(observed, predicted, ...)
}

#' @rdname hee_deep_time_helpers
#' @export
hee_time_axis <- function(start_ma = 540, end_ma = 0, by = 5) {
  seq(start_ma, end_ma, by = -abs(by))
}

#' @rdname hee_deep_time_helpers
#' @export
hee_pseudo_hindcast <- function(beta, env_cube, recipe, times, variables, species_origin = NULL) {
  hee_run_deep_time_pipeline(beta = beta, env_cube = env_cube, recipe = recipe,
                             times = times, variables = variables,
                             species_origin = species_origin)
}

#' @rdname hee_deep_time_helpers
#' @export
hee_run_full_pipeline <- function(...) {
  hee_run_deep_time_pipeline(...)
}

#' @rdname hee_deep_time_helpers
#' @export
hee_report <- function(result, file = NULL) {
  txt <- c("# HmscEcoEvo Deep-Time Report", "",
           paste("Generated:", Sys.time()), "",
           paste("Objects:", paste(names(result), collapse = ", ")))
  if (!is.null(file)) writeLines(txt, file)
  txt
}

#' @rdname hee_deep_time_helpers
#' @export
hee_write_outputs <- function(result, dir = "outputs") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  if (is.list(result)) {
    for (nm in names(result)) {
      z <- result[[nm]]
      if (is.data.frame(z)) utils::write.csv(z, file.path(dir, paste0(nm, ".csv")), row.names = FALSE)
    }
  } else if (is.data.frame(result)) {
    utils::write.csv(result, file.path(dir, "result.csv"), row.names = FALSE)
  }
  invisible(normalizePath(dir))
}

#' @rdname hee_deep_time_helpers
#' @export
hee_write_rasters <- function(x, dir = "rasters", prefix = "layer") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  d <- as.data.frame(x)
  utils::write.csv(d, file.path(dir, paste0(prefix, ".csv")), row.names = FALSE)
  invisible(file.path(dir, paste0(prefix, ".csv")))
}

.hee_auc <- function(obs, pred) {
  if (length(unique(obs)) < 2) return(NA_real_)
  r <- rank(pred)
  n1 <- sum(obs > 0); n0 <- sum(obs <= 0)
  (sum(r[obs > 0]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

.hee_calibration_slope <- function(obs, pred) {
  if (length(obs) < 3 || stats::sd(pred, na.rm = TRUE) == 0) return(NA_real_)
  stats::coef(stats::lm(obs ~ pred))[2]
}
