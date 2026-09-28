# Deep-time HMSC hindcasting workflow --------------------------------------

#' Load a Phanerozoic/deep-time environmental time cube
#'
#' Loads the compact 4-degree/1-degree Phanerozoic RDS product, or wraps an
#' already loaded compact environment. The returned object is deliberately thin:
#' it keeps the original RDS accessors and records the available ages and
#' variables without copying the full cube.
#'
#' @param path Character path to an `.rds` file or product directory. If a
#'   directory is supplied, the function searches `08_rds_export` first.
#' @param object Optional already loaded compact RDS environment.
#' @return A `hee_timecube` object.
#' @export
#' @examples
#' \dontrun{
#' cube <- hee_load_timecube("Phanerozoic_Environment_540_0Ma_5Myr_4deg_landharmonized_v2")
#' hee_index_timecube(cube)
#' }
hee_load_timecube <- function(path = NULL, object = NULL) {
  if (!is.null(object)) {
    cube <- object
    src <- NA_character_
  } else {
    if (is.null(path)) stop("Supply `path` or `object`.", call. = FALSE)
    src <- .hee_find_timecube_rds(path)
    cube <- readRDS(src)
  }
  .hee_validate_cube(cube)
  out <- list(
    cube = cube,
    path = src,
    times = as.numeric(cube$age_ma),
    variables = as.character(cube$variable_names),
    metadata = tryCatch(cube$get_metadata(), error = function(e) data.frame()),
    product = if (!is.null(cube$product)) cube$product else list()
  )
  class(out) <- "hee_timecube"
  out
}

#' Index a deep-time environmental time cube
#'
#' @param env_cube A `hee_timecube`, compact RDS environment, data.frame, or
#'   path accepted by [hee_load_timecube()].
#' @param variables Optional variables to index.
#' @param times Optional ages in Ma.
#' @param compute_summary Logical; compute min/max/mean by reading slices.
#'   Defaults to `FALSE` so indexing large products remains quick.
#' @return A data.frame with one row per time-variable pair.
#' @export
hee_index_timecube <- function(env_cube,
                               variables = NULL,
                               times = NULL,
                               compute_summary = FALSE) {
  cube <- .hee_as_timecube(env_cube)
  variables <- variables %||% cube$variables
  times <- times %||% cube$times
  miss <- setdiff(variables, cube$variables)
  if (length(miss) > 0) {
    stop("Time cube is missing variable(s): ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
  times <- .hee_match_times(times, cube$times)
  grid_n <- length(cube$cube$lon) * length(cube$cube$lat)
  idx <- expand.grid(time_ma = times, variable = variables,
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  idx$n_cell <- grid_n
  idx$resolution <- if (!is.null(cube$product$grid)) cube$product$grid else NA_character_
  idx$source_path <- cube$path
  if (nrow(cube$metadata) > 0) {
    idx <- merge(idx, cube$metadata, by = "variable", all.x = TRUE, sort = FALSE)
  }
  if (compute_summary) {
    idx$min <- idx$max <- idx$mean <- NA_real_
    for (tm in unique(idx$time_ma)) {
      sl <- .hee_get_slice(cube, tm, variables)
      for (v in variables) {
        ii <- idx$time_ma == tm & idx$variable == v
        idx$min[ii] <- suppressWarnings(min(sl[[v]], na.rm = TRUE))
        idx$max[ii] <- suppressWarnings(max(sl[[v]], na.rm = TRUE))
        idx$mean[ii] <- mean(sl[[v]], na.rm = TRUE)
      }
    }
  }
  idx
}

#' Check that a time cube contains required variables and ages
#'
#' @param env_cube A time cube.
#' @param variables Required variables.
#' @param times Required ages in Ma.
#' @param landmask Optional land mask variable to require.
#' @return A diagnostic list.
#' @export
hee_check_timecube <- function(env_cube,
                               variables = NULL,
                               times = NULL,
                               landmask = "land_mask_dem") {
  cube <- .hee_as_timecube(env_cube)
  variables <- variables %||% character()
  required <- unique(c(variables, landmask))
  missing_variables <- setdiff(required, cube$variables)
  missing_times <- numeric()
  if (!is.null(times)) {
    missing_times <- setdiff(as.numeric(times), cube$times)
  }
  ok <- length(missing_variables) == 0 && length(missing_times) == 0
  out <- list(
    ok = ok,
    n_time = length(cube$times),
    n_variable = length(cube$variables),
    missing_variables = missing_variables,
    missing_times = missing_times,
    messages = character(),
    warnings = character()
  )
  if (!ok) {
    out$warnings <- c(
      if (length(missing_variables) > 0) paste("Missing variables:", paste(missing_variables, collapse = ", ")),
      if (length(missing_times) > 0) paste("Missing ages:", paste(missing_times, collapse = ", "))
    )
  }
  class(out) <- "hee_check"
  out
}

#' Create file-safe taxon identifiers
#'
#' @param species_names Character vector of species names.
#' @param taxon_id_prefix Prefix for generated taxon ids.
#' @return A data.frame mapping original names to safe ids.
#' @export
hee_make_taxon_table <- function(species_names, taxon_id_prefix = "sp") {
  species_names <- as.character(species_names)
  if (anyNA(species_names) || any(species_names == "")) {
    stop("species_names must be non-empty strings.", call. = FALSE)
  }
  clean <- gsub("\\s+", "_", trimws(species_names))
  safe <- gsub("[^A-Za-z0-9_]+", "_", clean)
  safe <- gsub("_+", "_", safe)
  safe <- gsub("^_|_$", "", safe)
  safe[safe == ""] <- paste0(taxon_id_prefix, seq_len(sum(safe == "")))
  data.frame(
    taxon_id = sprintf("%s%03d", taxon_id_prefix, seq_along(species_names)),
    species_original = species_names,
    species_clean = clean,
    file_safe_name = make.unique(safe, sep = "_"),
    tree_tip_label = make.unique(clean, sep = "_"),
    stringsAsFactors = FALSE
  )
}

#' Check a community matrix
#'
#' @param comm Site by species community matrix/data.frame.
#' @param response_type Response type.
#' @return A diagnostic list.
#' @export
hee_check_comm <- function(comm,
                           response_type = c("presence", "count", "abundance", "cover")) {
  response_type <- match.arg(response_type)
  ok <- TRUE
  warnings <- messages <- character()
  if (is.data.frame(comm)) comm <- as.matrix(comm)
  if (!is.matrix(comm)) stop("comm must be a matrix or data.frame.", call. = FALSE)
  suppressWarnings(storage.mode(comm) <- "numeric")
  if (is.null(rownames(comm)) || any(rownames(comm) == "")) {
    ok <- FALSE; warnings <- c(warnings, "comm needs site row names.")
  }
  if (is.null(colnames(comm)) || any(colnames(comm) == "")) {
    ok <- FALSE; warnings <- c(warnings, "comm needs species column names.")
  }
  if (anyNA(comm)) warnings <- c(warnings, "comm contains NA values.")
  if (any(duplicated(colnames(comm)))) warnings <- c(warnings, "comm contains duplicated species names.")
  if (response_type == "presence" && any(!comm %in% c(0, 1, NA))) {
    ok <- FALSE; warnings <- c(warnings, "presence data must contain only 0/1/NA.")
  }
  if (response_type == "count" && any(comm < 0 | abs(comm - round(comm)) > .Machine$double.eps^0.5, na.rm = TRUE)) {
    ok <- FALSE; warnings <- c(warnings, "count data must be non-negative integers.")
  }
  rare <- colSums(comm > 0, na.rm = TRUE)
  zero_sites <- rownames(comm)[rowSums(comm > 0, na.rm = TRUE) == 0]
  zero_species <- colnames(comm)[rare == 0]
  list(ok = ok, n_site = nrow(comm), n_species = ncol(comm),
       zero_sites = zero_sites, zero_species = zero_species,
       rare_species = names(rare)[rare <= 1],
       zero_fraction = mean(comm == 0, na.rm = TRUE),
       messages = messages, warnings = unique(warnings))
}

#' Check modern environmental predictors against sites
#'
#' @param env_now Site-level environmental data.
#' @param sites Optional site ids.
#' @return A diagnostic list.
#' @export
hee_check_env_now <- function(env_now, sites = NULL) {
  env_now <- as.data.frame(env_now)
  sites <- sites %||% rownames(env_now)
  missing_sites <- setdiff(sites, rownames(env_now))
  extra_sites <- setdiff(rownames(env_now), sites)
  numeric_vars <- names(env_now)[vapply(env_now, is.numeric, logical(1))]
  na_vars <- names(env_now)[vapply(env_now, function(z) anyNA(z), logical(1))]
  cor_pairs <- data.frame()
  if (length(numeric_vars) > 1) {
    C <- stats::cor(env_now[, numeric_vars, drop = FALSE], use = "pairwise.complete.obs")
    wh <- which(abs(C) > 0.7 & upper.tri(C), arr.ind = TRUE)
    if (nrow(wh) > 0) {
      cor_pairs <- data.frame(var1 = rownames(C)[wh[, 1]], var2 = colnames(C)[wh[, 2]],
                              correlation = C[wh], stringsAsFactors = FALSE)
    }
  }
  list(ok = length(missing_sites) == 0,
       missing_sites = missing_sites, extra_sites = extra_sites,
       numeric_vars = numeric_vars, variables_with_na = na_vars,
       high_correlation_pairs = cor_pairs,
       warnings = c(if (length(missing_sites) > 0) "env_now is missing sites.",
                    if (length(na_vars) > 0) "env_now contains NA values."))
}

#' Check and match a phylogeny to species
#'
#' @param tree An `ape::phylo` tree.
#' @param species Species names.
#' @param drop_unmatched Drop tips not present in `species`.
#' @return A list with matched tree and diagnostics.
#' @export
hee_match_tree_species <- function(tree, species, drop_unmatched = TRUE) {
  .require_pkg("ape", "deep-time tree matching")
  if (!inherits(tree, "phylo")) stop("tree must be an ape::phylo object.", call. = FALSE)
  species <- as.character(species)
  missing_in_tree <- setdiff(species, tree$tip.label)
  extra_in_tree <- setdiff(tree$tip.label, species)
  matched <- tree
  if (drop_unmatched && length(extra_in_tree) > 0) {
    matched <- ape::drop.tip(matched, extra_in_tree)
  }
  list(tree = matched,
       missing_in_tree = missing_in_tree,
       extra_in_tree = extra_in_tree,
       ok = length(missing_in_tree) == 0)
}

#' @rdname hee_match_tree_species
#' @export
hee_check_tree <- function(tree, species) {
  out <- hee_match_tree_species(tree, species, drop_unmatched = FALSE)
  tr <- out$tree
  out$has_branch_lengths <- !is.null(tr$edge.length)
  out$negative_branch_lengths <- if (is.null(tr$edge.length)) FALSE else any(tr$edge.length < 0, na.rm = TRUE)
  out$is_ultrametric <- tryCatch(ape::is.ultrametric(tr), error = function(e) NA)
  out
}

#' Check a trait table
#'
#' @param traits Species by trait data.frame.
#' @param species Species names.
#' @return A diagnostic list.
#' @export
hee_check_traits <- function(traits, species) {
  traits <- as.data.frame(traits)
  missing_species <- setdiff(species, rownames(traits))
  extra_species <- setdiff(rownames(traits), species)
  miss_prop <- vapply(traits, function(z) mean(is.na(z)), numeric(1))
  list(ok = length(missing_species) == 0,
       missing_species = missing_species,
       extra_species = extra_species,
       missing_fraction = miss_prop,
       numeric_traits = names(traits)[vapply(traits, is.numeric, logical(1))])
}

#' Check a species-by-region range table
#'
#' @param tip_ranges Matrix/data.frame of species by region occupancy or probability.
#' @param species Species names.
#' @param max_range_size Optional maximum allowed occupied regions.
#' @return A diagnostic list.
#' @export
hee_check_tip_ranges <- function(tip_ranges, species, max_range_size = Inf) {
  x <- as.matrix(tip_ranges)
  suppressWarnings(storage.mode(x) <- "numeric")
  missing_species <- setdiff(species, rownames(x))
  range_size <- rowSums(x > 0, na.rm = TRUE)
  list(ok = length(missing_species) == 0 && all(range_size <= max_range_size),
       missing_species = missing_species,
       species_without_region = names(range_size)[range_size == 0],
       species_over_max_range = names(range_size)[range_size > max_range_size],
       probability_like = all(x >= 0 & x <= 1, na.rm = TRUE))
}

#' Lock and apply a projection recipe
#'
#' Stores modern centering/scaling and factor levels so paleo-environmental
#' projection uses exactly the same transformations as model calibration.
#'
#' @param XData Modern predictor table.
#' @param formula Optional model formula; variables in the formula are locked.
#' @param transform Named list. Numeric variables default to `"scale"`;
#'   categorical variables default to `"factor"`.
#' @return A `hee_recipe` object.
#' @export
hee_lock_recipe <- function(XData, formula = NULL, transform = NULL) {
  XData <- as.data.frame(XData)
  vars <- if (is.null(formula)) names(XData) else all.vars(formula)
  missing <- setdiff(vars, names(XData))
  if (length(missing) > 0) stop("XData is missing formula variable(s): ",
                                paste(missing, collapse = ", "), call. = FALSE)
  spec <- lapply(vars, function(v) {
    z <- XData[[v]]
    tr <- transform[[v]] %||% if (is.numeric(z)) "scale" else "factor"
    if (is.numeric(z)) {
      zz <- suppressWarnings(as.numeric(z))
      if (any(!is.na(zz) & !is.finite(zz))) {
        stop("Numeric recipe variable `", v, "` contains non-finite values.",
             call. = FALSE)
      }
      finite <- zz[is.finite(zz)]
      if (length(finite) == 0L) {
        stop("Numeric recipe variable `", v, "` has no finite values.",
             call. = FALSE)
      }
      mu <- mean(finite)
      sdv <- stats::sd(finite)
      if (!is.finite(sdv) || sdv == 0) sdv <- 1
    } else {
      mu <- NA_real_
      sdv <- NA_real_
    }
    list(variable = v, transform = tr,
         mean = mu,
         sd = sdv,
         levels = if (is.factor(z) || is.character(z)) sort(unique(as.character(z))) else NULL)
  })
  names(spec) <- vars
  out <- list(variables = vars, formula = formula, spec = spec,
              created_at = as.character(Sys.time()))
  class(out) <- "hee_recipe"
  out
}

#' @rdname hee_lock_recipe
#' @param recipe A `hee_recipe`.
#' @param newdata New predictor data.
#' @export
hee_apply_recipe <- function(recipe, newdata) {
  if (!inherits(recipe, "hee_recipe")) stop("recipe must be a hee_recipe.", call. = FALSE)
  newdata <- as.data.frame(newdata)
  missing <- setdiff(recipe$variables, names(newdata))
  if (length(missing) > 0) stop("newdata is missing variable(s): ",
                                paste(missing, collapse = ", "), call. = FALSE)
  rn <- rownames(newdata)
  if (is.null(rn) || any(!nzchar(rn))) rn <- as.character(seq_len(nrow(newdata)))
  out <- data.frame(row.names = rn)
  for (v in recipe$variables) {
    sp <- recipe$spec[[v]]
    z <- newdata[[v]]
    if (identical(sp$transform, "scale")) {
      if (!is.numeric(z)) stop("Variable `", v, "` must be numeric for scaling.", call. = FALSE)
      s <- sp$sd
      if (is.na(s) || s == 0) s <- 1
      out[[v]] <- (z - sp$mean) / s
    } else if (identical(sp$transform, "center")) {
      out[[v]] <- z - sp$mean
    } else if (identical(sp$transform, "identity")) {
      out[[v]] <- z
    } else if (identical(sp$transform, "factor")) {
      unseen <- setdiff(unique(as.character(z)), sp$levels)
      if (length(unseen) > 0) stop("Variable `", v, "` has unseen factor level(s): ",
                                   paste(unseen, collapse = ", "), call. = FALSE)
      out[[v]] <- factor(as.character(z), levels = sp$levels)
    } else {
      stop("Unknown transform for `", v, "`: ", sp$transform, call. = FALSE)
    }
  }
  out
}

#' Tag predictors and detect deep-time information leakage
#'
#' @param predictors Predictor names.
#' @param roles Named character vector/list of roles.
#' @return A data.frame of predictor roles.
#' @export
hee_tag_predictors <- function(predictors, roles = NULL) {
  predictors <- as.character(predictors)
  role <- rep("unknown", length(predictors))
  names(role) <- predictors
  if (!is.null(roles)) {
    role[names(roles)] <- as.character(roles)
  }
  data.frame(predictor = predictors, role = unname(role), stringsAsFactors = FALSE)
}

#' @rdname hee_tag_predictors
#' @param XFormula Formula to inspect.
#' @param predictor_roles Output from [hee_tag_predictors()].
#' @param task Task name; `"deep_time_projection"` blocks composition-derived predictors.
#' @export
hee_detect_leakage <- function(XFormula, predictor_roles,
                               task = c("deep_time_projection", "modern_explanation")) {
  task <- match.arg(task)
  vars <- all.vars(XFormula)
  pr <- as.data.frame(predictor_roles)
  bad <- pr$predictor[pr$predictor %in% vars &
                        pr$role %in% c("composition_derived", "current_community_derived")]
  if (task == "deep_time_projection" && length(bad) > 0) {
    stop("Information leakage: predictor(s) derived from current community composition ",
         "cannot be used for deep-time projection: ", paste(bad, collapse = ", "),
         call. = FALSE)
  }
  list(ok = TRUE, used_predictors = vars, blocked_predictors = bad)
}

#' Build a phylogenetic existence mask through time
#'
#' @param tree Optional dated phylogeny.
#' @param times Ages in Ma; larger values are older.
#' @param species_origin Named vector/data.frame of species origin ages in Ma.
#' @param species Species names. Defaults to names in `species_origin` or tree tips.
#' @param level Interpretation level.
#' @param missing_origin How to handle missing `species_origin`. The default
#'   `"error"` is conservative for species-level deep-time projection because
#'   modern species should not be allowed to exist before their origin age.
#'   `"oldest"` restores the historical neutral/demo behaviour and assumes each
#'   supplied species or tip existed from the oldest requested time. Use that
#'   option only for lineage/clade summaries, sensitivity analysis, or explicit
#'   demonstrations where this assumption is scientifically stated.
#' @return A species by time matrix with 0/1 existence weights.
#' @export
hee_phylo_time_mask <- function(tree = NULL,
                                times,
                                species_origin = NULL,
                                species = NULL,
                                level = c("species", "lineage", "clade"),
                                missing_origin = c("error", "oldest")) {
  level <- match.arg(level)
  missing_origin <- match.arg(missing_origin)
  times <- as.numeric(times)
  if (length(times) == 0L || any(!is.finite(times) | times < 0)) {
    stop("times must be finite non-negative ages in Ma.", call. = FALSE)
  }
  if (is.null(species_origin)) {
    if (!is.null(tree) && inherits(tree, "phylo")) {
      species <- species %||% tree$tip.label
    }
    if (is.null(species)) stop("Supply `species_origin`, `species`, or a tree.", call. = FALSE)
    if (identical(missing_origin, "error")) {
      stop("species_origin is required by default for species-level deep-time ",
           "projection. Assuming presence from the oldest time slice can ",
           "overestimate past distributions; set missing_origin = \"oldest\" ",
           "only for an explicit neutral lineage/clade or sensitivity scenario.",
           call. = FALSE)
    }
    origin <- stats::setNames(rep(max(times, na.rm = TRUE), length(species)), species)
  } else {
    if (is.data.frame(species_origin)) {
      .require_cols(species_origin, c("species", "origin_ma"), "species_origin")
      origin <- stats::setNames(species_origin$origin_ma, species_origin$species)
    } else {
      origin <- as.numeric(species_origin)
      origin_names <- names(species_origin)
      if (is.null(origin_names) || any(origin_names == "")) {
        if (is.null(species) || length(species) != length(origin)) {
          stop("Unnamed species_origin vectors require a `species` vector of the same length.",
               call. = FALSE)
        }
        origin_names <- as.character(species)
      }
      names(origin) <- origin_names
    }
    species <- species %||% names(origin)
    species <- as.character(species)
    origin <- origin[species]
  }
  if (any(is.na(origin) | !is.finite(origin) | origin < 0)) {
    bad <- names(origin)[is.na(origin) | !is.finite(origin) | origin < 0]
    stop("species_origin must contain finite non-negative origin_ma values for all species. Problem species: ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
  M <- outer(origin, times, FUN = function(o, t) as.numeric(t <= o))
  dimnames(M) <- list(species, paste0(times, "Ma"))
  attr(M, "times") <- times
  attr(M, "level") <- level
  M
}

#' Extract a paleoenvironmental grid for one or more ages
#'
#' @param env_cube Time cube.
#' @param times Ages in Ma.
#' @param variables Variables to extract.
#' @param landmask Land mask variable.
#' @param land_only Logical; keep only land cells with landmask > 0.
#' @return A data.frame.
#' @export
hee_make_paleo_grid <- function(env_cube,
                                times,
                                variables,
                                landmask = "land_mask_dem",
                                land_only = TRUE) {
  cube <- .hee_as_timecube(env_cube)
  missing_vars <- setdiff(variables, cube$variables)
  if (length(missing_vars) > 0L) {
    stop("env_cube is missing requested paleo variable(s): ",
         paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  if (isTRUE(land_only) && !landmask %in% cube$variables) {
    stop("land_only = TRUE requires landmask variable `", landmask,
         "` in env_cube.", call. = FALSE)
  }
  vars <- unique(c(variables, intersect(c(landmask, "land_area_km2"), cube$variables)))
  times <- .hee_match_times(times, cube$times)
  out <- lapply(times, function(tm) {
    sl <- .hee_get_slice(cube, tm, vars)
    sl$cell_id <- .hee_cell_id(sl$lon, sl$lat)
    if (land_only && landmask %in% names(sl)) sl <- sl[!is.na(sl[[landmask]]) & sl[[landmask]] > 0, , drop = FALSE]
    sl
  })
  do.call(rbind, out)
}

#' Estimate extrapolation risk for paleo predictors
#'
#' This function supports two call patterns. The original lightweight pattern is
#' `hee_extrapolation_risk(recipe, X_past_raw)`. The full deep-time pattern is
#' `hee_extrapolation_risk(env_now, env_cube, variables = ..., times = ...)`,
#' where the modern 0 Ma predictors define the training envelope and all
#' paleo-grid rows are scored against that envelope.
#'
#' @param recipe A `hee_recipe`, or a modern predictor data.frame for the full
#'   deep-time interface.
#' @param X_past_raw Paleo predictor table, or a time cube for the full
#'   deep-time interface.
#' @param threshold Maximum standardized/risk value before down-weighting.
#' @param method Extrapolation method: `"range"`, `"mahalanobis"`, `"pca"`, or
#'   `"mess"`. `"pca"` and `"mess"` are conservative range-based summaries in
#'   this lightweight implementation.
#' @param variables Predictor variables to score when using the time-cube
#'   interface.
#' @param times Ages in Ma; defaults to all ages in the time cube.
#' @param landmask Land mask variable.
#' @param land_only Keep only land cells when extracting the paleo grid.
#' @return Input rows with risk, flag, weight, and optional variable-wise
#'   out-of-range columns. Palaeo rows with missing predictor values are treated
#'   conservatively as extrapolation-risk rows: they are flagged and receive a
#'   reduced `extrapolation_weight` rather than being interpreted as reliable
#'   analog environments.
#' @export
hee_extrapolation_risk <- function(recipe, X_past_raw, threshold = 3,
                                   method = c("range", "mahalanobis", "pca", "mess"),
                                   variables = NULL, times = NULL,
                                   landmask = "land_mask_dem", land_only = TRUE) {
  method <- match.arg(method)
  if (!inherits(recipe, "hee_recipe")) {
    env_now <- as.data.frame(recipe)
    env_cube <- X_past_raw
    variables <- variables %||% names(env_now)[vapply(env_now, is.numeric, logical(1))]
    if (length(variables) == 0) stop("No numeric variables supplied for extrapolation scoring.", call. = FALSE)
    cube <- .hee_as_timecube(env_cube)
    missing_now <- setdiff(variables, names(env_now))
    if (length(missing_now) > 0L) {
      stop("Modern reference data are missing extrapolation variable(s): ",
           paste(missing_now, collapse = ", "), call. = FALSE)
    }
    missing_cube <- setdiff(variables, cube$variables)
    if (length(missing_cube) > 0L) {
      stop("Paleoenvironment cube is missing extrapolation variable(s): ",
           paste(missing_cube, collapse = ", "), call. = FALSE)
    }
    bad_reference <- variables[!vapply(variables, function(v) {
      z <- suppressWarnings(as.numeric(env_now[[v]]))
      any(is.finite(z))
    }, logical(1))]
    if (length(bad_reference) > 0L) {
      stop("Modern reference variable(s) have no finite values: ",
           paste(bad_reference, collapse = ", "), call. = FALSE)
    }
    times <- times %||% sort(cube$times, decreasing = TRUE)
    X_past_raw <- hee_make_paleo_grid(cube, times = times, variables = variables,
                                      landmask = landmask, land_only = land_only)
    out <- X_past_raw
    lo <- vapply(variables, function(v) min(env_now[[v]], na.rm = TRUE), numeric(1))
    hi <- vapply(variables, function(v) max(env_now[[v]], na.rm = TRUE), numeric(1))
    mu <- vapply(variables, function(v) mean(env_now[[v]], na.rm = TRUE), numeric(1))
    sdv <- vapply(variables, function(v) stats::sd(env_now[[v]], na.rm = TRUE), numeric(1))
    sdv[!is.finite(sdv) | sdv == 0] <- 1
    X <- as.matrix(out[, variables, drop = FALSE])
    if (method == "mahalanobis") {
      S <- stats::cov(env_now[, variables, drop = FALSE], use = "pairwise.complete.obs")
      if (any(!is.finite(S)) || qr(S)$rank < ncol(S)) S <- diag(sdv^2, length(sdv))
      score <- sqrt(stats::mahalanobis(X, center = mu, cov = S))
      flag_mat <- sweep(X, 2, lo, "<") | sweep(X, 2, hi, ">")
    } else {
      below <- sweep(X, 2, lo, "<")
      above <- sweep(X, 2, hi, ">")
      flag_mat <- below | above
      z <- sweep(sweep(X, 2, mu, "-"), 2, sdv, "/")
      excess_low <- sweep(X, 2, lo, function(x, y) y - x)
      excess_high <- sweep(X, 2, hi, "-")
      excess <- pmax(excess_low, excess_high, 0)
      finite_n <- rowSums(is.finite(X))
      all_missing <- finite_n == 0
      score <- if (method == "mess") {
        apply(excess / matrix(sdv, nrow = nrow(X), ncol = length(sdv), byrow = TRUE), 1, function(v) {
          if (all(!is.finite(v))) NA_real_ else max(v, na.rm = TRUE)
        })
      } else {
        max_abs_z <- apply(abs(z), 1, function(v) {
          if (all(!is.finite(v))) NA_real_ else max(v, na.rm = TRUE)
        })
        n_out <- rowSums(flag_mat, na.rm = TRUE)
        pmax(max_abs_z, ifelse(n_out > 0, threshold + n_out, 0))
      }
      score[all_missing] <- threshold + length(variables)
      missing_n <- length(variables) - finite_n
      score[missing_n > 0] <- pmax(
        score[missing_n > 0],
        threshold + missing_n[missing_n > 0] / length(variables),
        na.rm = TRUE
      )
    }
    score[!is.finite(score)] <- threshold + length(variables)
    out$extrapolation_score <- score
    out$extrapolation_risk <- score
    missing_env <- rowSums(is.finite(X)) < length(variables)
    out$extrapolation_flag <- score > threshold | rowSums(flag_mat, na.rm = TRUE) > 0 |
      missing_env
    out$extrapolation_flag[is.na(out$extrapolation_flag)] <- TRUE
    out$no_analog_fraction <- stats::ave(as.numeric(out$extrapolation_flag), out$time_ma, FUN = mean)
    score_for_weight <- score
    score_for_weight[is.na(score_for_weight)] <- 0
    out$extrapolation_weight <- exp(-pmax(score_for_weight - threshold, 0))
    for (v in variables) out[[paste0("out_of_range_", v)]] <- flag_mat[, v]
    out$variable_out_of_range <- apply(flag_mat, 1, function(z) {
      hit <- variables[as.logical(z)]
      if (length(hit) == 0) "" else paste(hit, collapse = ";")
    })
    return(out)
  }
  X <- hee_apply_recipe(recipe, X_past_raw)
  num <- X[vapply(X, is.numeric, logical(1))]
  risk <- if (ncol(num) == 0) rep(0, nrow(X)) else apply(abs(as.matrix(num)), 1, function(v) {
    if (all(!is.finite(v))) threshold + ncol(num) else max(v, na.rm = TRUE)
  })
  risk[!is.finite(risk)] <- threshold + max(ncol(num), 1L)
  if (ncol(num) == 0) {
    missing_env <- rep(FALSE, nrow(X))
  } else {
    finite_n <- rowSums(is.finite(as.matrix(num)))
    missing_n <- ncol(num) - finite_n
    missing_env <- missing_n > 0
    risk[missing_env] <- pmax(
      risk[missing_env],
      threshold + missing_n[missing_env] / ncol(num),
      na.rm = TRUE
    )
  }
  risk_for_weight <- risk
  weight <- exp(-pmax(risk_for_weight - threshold, 0))
  cbind(as.data.frame(X_past_raw), extrapolation_score = risk,
        extrapolation_risk = risk,
        extrapolation_flag = risk > threshold | missing_env,
        extrapolation_weight = weight)
}

#' Project HMSC Beta responses through a deep-time environmental cube
#'
#' This lightweight projector uses a species by axis Beta matrix and an optional
#' recipe. It is intended for fast environmental-only hindcasts and examples.
#' Full Hmsc posterior prediction can be supplied externally and then combined
#' with [hee_combine()].
#'
#' @param evo Optional `hmsc_ecoevo` object supplying `beta` and `distr`.
#' @param beta Species by predictor matrix.
#' @param env_cube Time cube.
#' @param recipe Optional `hee_recipe`.
#' @param times Ages in Ma.
#' @param variables Paleoenvironmental variables to extract.
#' @param species Optional species subset.
#' @param link Link/response scale: `"probit"`, `"logit"`, `"poisson"`, or `"identity"`.
#' @param intercept_axis Name of intercept column in Beta.
#' @param landmask Land mask variable.
#' @param land_only Keep only land cells.
#' @return A `hee_suitability` data.frame.
#' @export
hee_project_hmsc_timecube <- function(evo = NULL,
                                      beta = NULL,
                                      env_cube,
                                      recipe = NULL,
                                      times,
                                      variables = NULL,
                                      species = NULL,
                                      link = NULL,
                                      intercept_axis = "(Intercept)",
                                      landmask = "land_mask_dem",
                                      land_only = TRUE) {
  if (!is.null(evo)) {
    evo <- as_hmsc_ecoevo(evo)
    beta <- beta %||% evo$beta
    link <- link %||% evo$distr
  }
  if (is.null(beta)) stop("Supply `beta` or `evo`.", call. = FALSE)
  B <- .ensure_named_matrix(beta, "beta")
  if (!is.null(species)) {
    missing_species <- setdiff(species, rownames(B))
    if (length(missing_species) > 0L) {
      stop("beta is missing requested species: ",
           paste(missing_species, collapse = ", "), call. = FALSE)
    }
    B <- B[species, , drop = FALSE]
  }
  variables <- variables %||% setdiff(colnames(B), intercept_axis)
  grid <- hee_make_paleo_grid(env_cube, times = times, variables = variables,
                              landmask = landmask, land_only = land_only)
  X <- if (is.null(recipe)) {
    missing <- setdiff(variables, names(grid))
    if (length(missing) > 0) stop("Paleo grid missing variable(s): ",
                                  paste(missing, collapse = ", "), call. = FALSE)
    grid[, variables, drop = FALSE]
  } else {
    hee_apply_recipe(recipe, grid)
  }
  axes <- setdiff(colnames(B), intercept_axis)
  missing_axes <- setdiff(axes, names(X))
  if (length(missing_axes) > 0L) {
    stop("Paleo predictors are missing Beta axis/axes: ",
         paste(missing_axes, collapse = ", "), call. = FALSE)
  }
  if (length(axes) == 0) stop("No non-intercept Beta axis is available for projection.", call. = FALSE)
  Xmat <- as.matrix(X[, axes, drop = FALSE])
  Bsub <- B[, axes, drop = FALSE]
  eta <- Xmat %*% t(Bsub)
  if (intercept_axis %in% colnames(B)) {
    eta <- sweep(eta, 2, B[, intercept_axis], `+`)
  }
  pred <- .hee_link_inverse(eta, link %||% "probit")
  out <- cbind(grid[, intersect(c("age_ma", "time_ma", "lon", "lat", "cell_id",
                                  landmask, "land_area_km2"), names(grid)), drop = FALSE],
               as.data.frame(pred, check.names = FALSE))
  if (!"time_ma" %in% names(out)) names(out)[names(out) == "age_ma"] <- "time_ma"
  long <- stats::reshape(out,
    varying = rownames(B),
    v.names = "suitability",
    timevar = "species",
    times = rownames(B),
    direction = "long"
  )
  rownames(long) <- NULL
  class(long) <- c("hee_suitability", class(long))
  long
}

#' Extract paleoenvironment along plate-corrected tracks
#'
#' The function assumes the supplied coordinates are already in the paleo frame.
#' It does not perform plate reconstruction internally; use GPlates or another
#' plate model first, then pass the reconstructed coordinates here.
#'
#' @param points Data.frame with paleo coordinates and ages.
#' @param env_cube Time cube.
#' @param variables Variables to extract.
#' @param time_col Age column.
#' @param lon_col Paleo-longitude column.
#' @param lat_col Paleo-latitude column.
#' @param landmask Land mask variable.
#' @param land_only Logical; if `TRUE`, match to the nearest land cell.
#' @return Points with nearest-cell paleoenvironmental values, matched grid
#'   coordinates, legacy angular distance (`matched_distance_deg`), and
#'   auditable great-circle distance in kilometres (`matched_distance_km`).
#' @export
hee_extract_paleoenv_track <- function(points,
                                       env_cube,
                                       variables,
                                       time_col = "time_ma",
                                       lon_col = "paleo_lon",
                                       lat_col = "paleo_lat",
                                       landmask = "land_mask_dem",
                                       land_only = FALSE) {
  points <- as.data.frame(points)
  .require_cols(points, c(time_col, lon_col, lat_col), "points")
  cube <- .hee_as_timecube(env_cube)
  missing_vars <- setdiff(variables, cube$variables)
  if (length(missing_vars) > 0L) {
    stop("env_cube is missing requested paleo track variable(s): ",
         paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  if (isTRUE(land_only) && !landmask %in% cube$variables) {
    stop("land_only = TRUE requires landmask variable `", landmask,
         "` in env_cube.", call. = FALSE)
  }
  vars <- unique(c(variables, if (land_only) landmask else character()))
  out <- lapply(seq_len(nrow(points)), function(i) {
    pt <- points[i, , drop = FALSE]
    tm <- .hee_match_times(pt[[time_col]], cube$times)
    sl <- .hee_get_slice(cube, tm, vars)
    if (land_only && landmask %in% names(sl)) {
      sl <- sl[!is.na(sl[[landmask]]) & sl[[landmask]] > 0, , drop = FALSE]
    }
    if (nrow(sl) == 0) stop("No candidate cells available at time ", tm, " Ma.", call. = FALSE)
    distance_km <- hee_great_circle_distance_km(
      pt[[lon_col]], pt[[lat_col]], sl$lon, sl$lat
    )
    j <- which.min(distance_km)
    selected <- sl[j, setdiff(names(sl), "time_ma"), drop = FALSE]
    matched_lon <- selected$lon
    matched_lat <- selected$lat
    if ("cell_id" %in% names(selected)) {
      names(selected)[names(selected) == "cell_id"] <- "matched_grid_cell_id"
    }
    selected$lon <- NULL
    selected$lat <- NULL
    angular_distance <- sqrt((matched_lon - pt[[lon_col]])^2 +
                               (matched_lat - pt[[lat_col]])^2)
    cbind(pt, selected,
          matched_grid_lon = matched_lon,
          matched_grid_lat = matched_lat,
          matched_time_ma = tm,
          matched_distance_deg = angular_distance,
          matched_distance_km = distance_km[[j]])
  })
  do.call(rbind, out)
}

#' Project HMSC Beta responses for an arbitrary paleo predictor table
#'
#' @param evo Optional `hmsc_ecoevo` object.
#' @param beta Optional species by predictor matrix.
#' @param X_past Paleo predictor table.
#' @param recipe Optional projection recipe.
#' @param id_cols Columns to keep in the long output.
#' @param species Optional species subset.
#' @param link Link/response scale.
#' @param intercept_axis Intercept column in Beta.
#' @return A long `hee_suitability` table.
#' @export
hee_project_hmsc_table <- function(evo = NULL,
                                   beta = NULL,
                                   X_past,
                                   recipe = NULL,
                                   id_cols = c("site_id", "track_id", "time_ma",
                                               "matched_time_ma", "lon", "lat",
                                               "paleo_lon", "paleo_lat", "cell_id"),
                                   species = NULL,
                                   link = NULL,
                                   intercept_axis = "(Intercept)") {
  if (!is.null(evo)) {
    evo <- as_hmsc_ecoevo(evo)
    beta <- beta %||% evo$beta
    link <- link %||% evo$distr
  }
  if (is.null(beta)) stop("Supply `beta` or `evo`.", call. = FALSE)
  B <- .ensure_named_matrix(beta, "beta")
  if (!is.null(species)) {
    missing_species <- setdiff(species, rownames(B))
    if (length(missing_species) > 0L) {
      stop("beta is missing requested species: ",
           paste(missing_species, collapse = ", "), call. = FALSE)
    }
    B <- B[species, , drop = FALSE]
  }
  X_past <- as.data.frame(X_past)
  if (!"cell_id" %in% names(X_past) && all(c("lon", "lat") %in% names(X_past))) {
    X_past$cell_id <- .hee_cell_id(X_past$lon, X_past$lat)
  }
  X <- if (is.null(recipe)) X_past else hee_apply_recipe(recipe, X_past)
  axes <- setdiff(colnames(B), intercept_axis)
  missing_axes <- setdiff(axes, names(X))
  if (length(missing_axes) > 0L) {
    stop("X_past predictors are missing Beta axis/axes: ",
         paste(missing_axes, collapse = ", "), call. = FALSE)
  }
  if (length(axes) == 0) stop("No non-intercept Beta axis is available for projection.", call. = FALSE)
  eta <- as.matrix(X[, axes, drop = FALSE]) %*% t(B[, axes, drop = FALSE])
  if (intercept_axis %in% colnames(B)) eta <- sweep(eta, 2, B[, intercept_axis], `+`)
  pred <- .hee_link_inverse(eta, link %||% "probit")
  out <- cbind(X_past[, intersect(id_cols, names(X_past)), drop = FALSE],
               as.data.frame(pred, check.names = FALSE))
  long <- stats::reshape(out, varying = rownames(B), v.names = "suitability",
                         timevar = "species", times = rownames(B),
                         direction = "long")
  rownames(long) <- NULL
  if (!"time_ma" %in% names(long) && "matched_time_ma" %in% names(long)) {
    long$time_ma <- long$matched_time_ma
  }
  class(long) <- c("hee_suitability", class(long))
  long
}

#' Combine environmental suitability with historical and phylogenetic filters
#'
#' @param suitability Output from [hee_project_hmsc_timecube()] or equivalent
#'   long table with `species`, `cell_id`, `time_ma`, and `suitability` (or
#'   the framework alias `S_HMSC`). Missing suitability is recorded in
#'   `missing_suitability` and treated conservatively as zero probability.
#' @param accessibility Optional table or matrix of accessibility weights. If
#'   the whole component is `NULL` it is neutral (`1`); if a component table is
#'   supplied, unmatched projection rows are treated as no evidence for
#'   accessibility (`0`). Recognised aliases include `A_hist` and `A_BGB`.
#' @param phylo_mask Optional species by time matrix from [hee_phylo_time_mask()].
#' @param landmask Optional cell/time land mask table. Recognised columns are
#'   `land_weight`, `geographic_existence`, `land_mask_dem`, `land`, and
#'   `habitat_availability`, plus framework aliases `L_land` and `G_arena`;
#'   when several are supplied they are multiplied. If a landmask table is
#'   supplied, unmatched rows are treated as geographically absent (`0`).
#' @param static_filter Optional static dispersal filter table. Recognised
#'   composite columns are `static_weight`, `D_static`, and `weight`. If no
#'   composite column is supplied, `structural_connectivity`,
#'   `functional_connectivity`, and `climatic_connectivity` are multiplied.
#'   When a static-filter table is supplied, unmatched rows are treated as not
#'   statically reachable (`0`).
#' @param dynamic_filter Optional dynamic dispersal filter table. Recognised
#'   columns are `dynamic_weight`, `D_dynamic`, and `weight`. When supplied,
#'   unmatched rows are treated as not dynamically reachable (`0`).
#' @param extrapolation_weight Optional cell/time `Q_extrap` table. Recognised
#'   columns are `Q_extrap`, `extrapolation_weight`, and `weight`. Omitting the
#'   whole table is neutral (`Q_extrap = 1`). If a table is supplied, explicit
#'   `NA` or unmatched projection rows are conservative (`Q_extrap = 0`) so
#'   unassessed no-analog risk is not treated as full reliability. Probabilities
#'   are clipped to `[0, 1]`.
#' @param method Currently `"multiply"`.
#' @param allow_species_broadcast Logical; default `FALSE`. Species-level
#'   components such as accessibility and phylogenetic existence must normally
#'   contain a `species` key when more than one species is projected. Set this
#'   to `TRUE` only for an explicitly shared species-neutral scenario.
#' @return A `hee_projection` data.frame with final `probability` and explicit
#'   framework components `S_HMSC`, `A_hist`, `E_phylo`, `L_land`,
#'   `D_static`, `D_dynamic`, and `Q_extrap`. Projection identity may be a
#'   full-grid `cell_id`, reconstructed-point `track_id`, `region`, or
#'   `site_id`; if none is present, `species x time_ma` must be unique.
#' @export
hee_combine <- function(suitability,
                        accessibility = NULL,
                        phylo_mask = NULL,
                        landmask = NULL,
                        static_filter = NULL,
                        dynamic_filter = NULL,
                        extrapolation_weight = NULL,
                        method = "multiply",
                        allow_species_broadcast = FALSE) {
  if (!identical(method, "multiply")) stop("Only method = 'multiply' is currently supported.", call. = FALSE)
  P <- as.data.frame(suitability)
  if (!"suitability" %in% names(P) && "S_HMSC" %in% names(P)) {
    P$suitability <- P$S_HMSC
  }
  .require_cols(P, c("species", "time_ma", "suitability"), "suitability")
  identity_col <- intersect(c("cell_id", "track_id", "region", "site_id"),
                            names(P))
  if (length(identity_col) > 1L) identity_col <- identity_col[[1]]
  .hee_check_unique_keys(P, c("species", identity_col, "time_ma"),
                         "suitability")
  .hee_require_any_columns(accessibility, c("accessibility", "A_hist",
                                            "A_BGB", "weight"),
                           "accessibility")
  .hee_require_any_columns(landmask, c("land_weight", "geographic_existence",
                                       "land_mask_dem", "land",
                                       "habitat_availability", "L_land",
                                       "G_arena"),
                           "landmask")
  .hee_require_any_columns(static_filter, c("static_weight", "D_static", "weight",
                                            "structural_connectivity",
                                            "functional_connectivity",
                                            "climatic_connectivity",
                                            "climate_connectivity"),
                           "static_filter")
  .hee_require_any_columns(dynamic_filter, c("dynamic_weight", "D_dynamic", "weight"),
                           "dynamic_filter")
  .hee_require_any_columns(extrapolation_weight, c("Q_extrap", "extrapolation_weight",
                                                   "weight"),
                           "extrapolation_weight")
  .hee_reject_species_broadcast(P, accessibility, "accessibility",
                                allow_species_broadcast)
  .hee_reject_species_broadcast(P, phylo_mask, "phylo_mask",
                                allow_species_broadcast)
  P$accessibility <- .hee_lookup_weight(P, accessibility, "accessibility",
                                        default = if (is.null(accessibility)) 1 else 0,
                                        aliases = c("A_hist", "A_BGB", "weight"),
                                        table_name = "accessibility")
  P$phylo_existence <- .hee_lookup_phylo(P, phylo_mask)
  P$land_weight <- .hee_lookup_product(
    P, landmask,
    cols = c("land_weight", "geographic_existence", "land_mask_dem",
             "land", "habitat_availability", "L_land", "G_arena"),
    default = if (is.null(landmask)) 1 else 0,
    table_name = "landmask"
  )
  if (is.null(landmask)) {
    local_land <- intersect(c("land_weight", "geographic_existence",
                              "land_mask_dem", "land",
                              "habitat_availability", "L_land",
                              "G_arena"), names(P))
    if (length(local_land) > 0L) {
      P$land_weight <- .hee_clip01(Reduce(`*`, lapply(local_land, function(nm) {
        v <- .hee_clip01(P[[nm]])
        v[is.na(v)] <- 0
        v
      })))
    }
  } else {
    P$land_weight[!.hee_rows_matched(P, landmask, "landmask")] <- 0
  }
  static_tab <- if (is.null(static_filter)) NULL else as.data.frame(static_filter)
  static_component_default <- function(cols) {
    if (is.null(static_tab) || !any(cols %in% names(static_tab))) 1 else 0
  }
  P$structural_connectivity <- .hee_lookup_weight(
    P, static_filter, "structural_connectivity",
    default = static_component_default("structural_connectivity"),
    table_name = "static_filter"
  )
  P$functional_connectivity <- .hee_lookup_weight(
    P, static_filter, "functional_connectivity",
    default = static_component_default("functional_connectivity"),
    table_name = "static_filter"
  )
  P$climatic_connectivity <- .hee_lookup_weight(
    P, static_filter, "climatic_connectivity",
    default = static_component_default(c("climatic_connectivity",
                                         "climate_connectivity")),
    aliases = "climate_connectivity",
    table_name = "static_filter"
  )
  P$static_weight <- .hee_lookup_weight(
    P, static_filter, "static_weight", default = NA_real_,
    aliases = c("D_static", "weight"),
    table_name = "static_filter"
  )
  static_has_composite <- !is.null(static_filter) &&
    any(c("static_weight", "D_static", "weight") %in% names(as.data.frame(static_filter)))
  missing_static <- is.na(P$static_weight)
  if (any(missing_static) && !static_has_composite) {
    P$static_weight[missing_static] <- P$structural_connectivity[missing_static] *
      P$functional_connectivity[missing_static] *
      P$climatic_connectivity[missing_static]
  } else if (any(missing_static)) {
    P$static_weight[missing_static] <- 0
  }
  if (!is.null(static_filter)) {
    P$static_weight[!.hee_rows_matched(P, static_filter, "static_filter")] <- 0
  }
  P$dynamic_weight <- .hee_lookup_weight(
    P, dynamic_filter, "dynamic_weight",
    default = if (is.null(dynamic_filter)) 1 else 0,
    aliases = c("D_dynamic", "weight"),
    table_name = "dynamic_filter"
  )
  if (!is.null(dynamic_filter)) {
    P$dynamic_weight[!.hee_rows_matched(P, dynamic_filter, "dynamic_filter")] <- 0
  }
  P$Q_extrap <- .hee_lookup_weight(
    P, extrapolation_weight, "Q_extrap",
    default = if (is.null(extrapolation_weight)) 1 else 0,
    aliases = c("extrapolation_weight", "weight"),
    table_name = "extrapolation_weight"
  )
  if (!is.null(extrapolation_weight)) {
    P$Q_extrap[!.hee_rows_matched(P, extrapolation_weight, "extrapolation_weight")] <- 0
  }
  S <- .hee_clip01(P$suitability)
  P$missing_suitability <- is.na(S)
  S[is.na(S)] <- 0
  P$S_HMSC <- S
  P$A_hist <- P$accessibility
  P$E_phylo <- P$phylo_existence
  P$L_land <- P$land_weight
  P$probability <- S * P$accessibility * P$phylo_existence *
    P$land_weight * P$static_weight * P$dynamic_weight * P$Q_extrap
  hard_zero <- (P$accessibility <= 0 | P$phylo_existence <= 0 |
                  P$land_weight <= 0 | P$static_weight <= 0 |
                  P$dynamic_weight <= 0 | P$Q_extrap <= 0)
  hard_zero[is.na(hard_zero)] <- FALSE
  P$probability[hard_zero] <- 0
  P$D_static <- P$static_weight
  P$D_dynamic <- P$dynamic_weight
  P$probability <- pmin(pmax(P$probability, 0), 1)
  class(P) <- c("hee_projection", class(P))
  P
}

#' Calculate expected richness through time
#'
#' @param projection Output from [hee_combine()] or table with probability.
#' @param threshold Optional threshold for binary richness.
#' @param group_cols Grouping columns. Defaults to `cell_id` and `time_ma`;
#'   use `c("track_id", "time_ma")` for tectonic-track summaries.
#' @return A data.frame with expected richness per cell/time.
#' @export
hee_richness <- function(projection, threshold = NULL,
                         group_cols = c("cell_id", "time_ma")) {
  P <- as.data.frame(projection)
  .require_cols(P, c(group_cols, "probability"), "projection")
  P$.value <- if (is.null(threshold)) P$probability else as.numeric(P$probability >= threshold)
  f <- stats::as.formula(paste(".value ~", paste(group_cols, collapse = " + ")))
  out <- stats::aggregate(f, P, sum, na.rm = TRUE)
  names(out)[names(out) == ".value"] <- if (is.null(threshold)) "expected_richness" else "binary_richness"
  coords <- unique(P[, intersect(c(group_cols, "lon", "lat", "paleo_lon", "paleo_lat"), names(P)), drop = FALSE])
  coord_cols <- setdiff(names(coords), group_cols)
  if (length(coord_cols) > 0L) {
    .hee_check_unique_keys(coords, group_cols, "projection coordinates")
    out <- merge(out, coords, by = group_cols, all.x = TRUE, sort = FALSE)
  }
  class(out) <- c("hee_richness", class(out))
  out
}

#' Calculate community turnover between adjacent time slices
#'
#' @param projection Projection table.
#' @return A data.frame with Bray-Curtis-like turnover.
#' @export
hee_turnover <- function(projection) {
  P <- as.data.frame(projection)
  .require_cols(P, c("cell_id", "time_ma", "species", "probability"), "projection")
  times <- sort(unique(P$time_ma), decreasing = TRUE)
  rows <- list()
  k <- 1L
  for (i in seq_len(length(times) - 1L)) {
    a <- times[i]; b <- times[i + 1L]
    Pa <- P[P$time_ma == a, c("cell_id", "species", "probability")]
    Pb <- P[P$time_ma == b, c("cell_id", "species", "probability")]
    M <- merge(Pa, Pb, by = c("cell_id", "species"), all = TRUE, suffixes = c("_older", "_younger"))
    M$probability_older[is.na(M$probability_older)] <- 0
    M$probability_younger[is.na(M$probability_younger)] <- 0
    sp <- split(M, M$cell_id)
    rows[[k]] <- do.call(rbind, lapply(sp, function(z) {
      denom <- sum(pmax(z$probability_older, z$probability_younger), na.rm = TRUE)
      data.frame(cell_id = z$cell_id[1], time_from_ma = a, time_to_ma = b,
                 turnover = if (denom > 0) 1 - sum(pmin(z$probability_older, z$probability_younger), na.rm = TRUE) / denom else NA_real_,
                 stringsAsFactors = FALSE)
    }))
    k <- k + 1L
  }
  if (length(rows) == 0) data.frame() else do.call(rbind, rows)
}

#' Summarize long-term refugia potential
#'
#' @param richness Output from [hee_richness()].
#' @param group_col Spatial or track grouping column.
#' @return A cell-level refugia summary.
#' @export
hee_refugia <- function(richness, group_col = NULL) {
  R <- as.data.frame(richness)
  group_col <- group_col %||% if ("cell_id" %in% names(R)) "cell_id" else names(R)[1]
  val <- if ("expected_richness" %in% names(R)) "expected_richness" else "binary_richness"
  mn <- stats::aggregate(R[[val]], list(group = R[[group_col]]), mean, na.rm = TRUE)
  sdv <- stats::aggregate(R[[val]], list(group = R[[group_col]]), stats::sd, na.rm = TRUE)
  names(mn)[1] <- group_col
  names(sdv)[1] <- group_col
  names(mn)[2] <- "mean_richness"
  names(sdv)[2] <- "richness_sd"
  out <- merge(mn, sdv, by = group_col, all = TRUE, sort = FALSE)
  ntime <- stats::aggregate(R[[val]], list(group = R[[group_col]]), function(z) sum(!is.na(z)))
  names(ntime)[1] <- group_col
  names(ntime)[2] <- "n_time"
  out <- merge(out, ntime, by = group_col, all.x = TRUE, sort = FALSE)
  out$richness_sd[is.na(out$richness_sd) & out$n_time <= 1] <- 0
  out$stability <- 1 / (1 + out$richness_sd)
  z <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    ans <- rep(0, length(x))
    ok <- is.finite(x)
    if (!any(ok)) return(ans)
    s <- stats::sd(x[ok], na.rm = TRUE)
    if (!is.finite(s) || s <= 0) return(ans)
    ans[ok] <- (x[ok] - mean(x[ok], na.rm = TRUE)) / s
    ans
  }
  out$refugia_score <- z(out$mean_richness) + z(out$stability)
  coord_cols <- intersect(c("lon", "lat", "paleo_lon", "paleo_lat"), names(R))
  if (length(coord_cols) > 0) {
    coords <- stats::aggregate(R[, coord_cols, drop = FALSE],
                               list(group = R[[group_col]]), mean, na.rm = TRUE)
    names(coords)[1] <- group_col
    out <- merge(out, coords, by = group_col, all.x = TRUE, sort = FALSE)
  }
  out
}

#' Train and apply probability thresholds
#'
#' @param observed Observed 0/1 values.
#' @param predicted Predicted probabilities.
#' @param method Threshold criterion.
#' @return Threshold object.
#' @export
hee_threshold_train <- function(observed, predicted, method = c("max_tss", "prevalence")) {
  method <- match.arg(method)
  observed <- as.numeric(observed)
  predicted <- as.numeric(predicted)
  ok <- !is.na(observed) & !is.na(predicted)
  observed <- observed[ok]; predicted <- predicted[ok]
  if (method == "prevalence") {
    thr <- mean(observed > 0, na.rm = TRUE)
  } else {
    cand <- sort(unique(predicted))
    score <- vapply(cand, function(th) {
      pred <- predicted >= th
      sens <- sum(pred & observed == 1) / max(sum(observed == 1), 1)
      spec <- sum(!pred & observed == 0) / max(sum(observed == 0), 1)
      sens + spec - 1
    }, numeric(1))
    thr <- cand[which.max(score)]
  }
  structure(list(threshold = thr, method = method), class = "hee_threshold")
}

#' @rdname hee_threshold_train
#' @param threshold A threshold object or numeric threshold.
#' @param values Predicted probabilities.
#' @export
hee_threshold_apply <- function(threshold, values) {
  th <- if (inherits(threshold, "hee_threshold")) threshold$threshold else as.numeric(threshold)[1]
  as.numeric(values >= th)
}

#' Validate fossil or paleo-record evidence against predictions
#'
#' @param fossils Table with `lon`, `lat`, `age_min`, `age_max`; species may be
#'   supplied as `species` or `taxon_id`.
#' @param projection Projection table.
#' @param threshold Optional hit threshold.
#' @return A validation table. `matched_time_ma` is the available projection
#'   time nearest the midpoint of the fossil age interval. Spatial matching is
#'   then performed within that time slice using great-circle distance. The
#'   returned `temporal_distance_myr` and `spatial_distance_km` make both
#'   approximations explicit. Coordinates must already be expressed in the
#'   same palaeogeographic frame as the projection; this function does not
#'   perform plate reconstruction.
#' @export
hee_validate_fossils <- function(fossils, projection, threshold = 0.5) {
  fossils <- as.data.frame(fossils)
  P <- as.data.frame(projection)
  .require_cols(fossils, c("lon", "lat", "age_min", "age_max"), "fossils")
  .require_cols(P, c("lon", "lat", "time_ma", "species", "probability"), "projection")
  sp_col <- if ("species" %in% names(fossils)) "species" else if ("taxon_id" %in% names(fossils)) "taxon_id" else NA_character_
  out <- lapply(seq_len(nrow(fossils)), function(i) {
    f <- fossils[i, , drop = FALSE]
    ages <- sort(as.numeric(c(f$age_min, f$age_max)))
    if (any(!is.finite(ages)) || any(ages < 0)) {
      stop("fossils age_min and age_max must be finite non-negative Ma.",
           call. = FALSE)
    }
    age_midpoint <- mean(ages)
    species_pool <- if (!is.na(sp_col)) {
      P[P$species == f[[sp_col]], , drop = FALSE]
    } else {
      P
    }
    available_times <- sort(unique(as.numeric(species_pool$time_ma)),
                            decreasing = TRUE)
    if (length(available_times) == 0L) {
      cand <- species_pool
      matched_time <- NA_real_
    } else {
      matched_time <- available_times[[which.min(abs(available_times -
                                                       age_midpoint))]]
      cand <- species_pool[species_pool$time_ma == matched_time, , drop = FALSE]
    }
    if (nrow(cand) == 0) {
      return(cbind(f, predicted_probability = NA_real_, matched_time_ma = NA_real_,
                   matched_cell_id = NA_character_, temporal_distance_myr = NA_real_,
                   matched_within_age_interval = NA,
                   spatial_distance_km = NA_real_, hit = NA))
    }
    distance_km <- hee_great_circle_distance_km(
      f$lon, f$lat, cand$lon, cand$lat
    )
    if (all(!is.finite(distance_km))) {
      return(cbind(f, predicted_probability = NA_real_,
                   matched_time_ma = matched_time,
                   matched_cell_id = NA_character_,
                   temporal_distance_myr = abs(matched_time - age_midpoint),
                   matched_within_age_interval = matched_time >= ages[[1]] &&
                     matched_time <= ages[[2]], spatial_distance_km = NA_real_,
                   hit = NA))
    }
    j <- which.min(distance_km)
    cbind(f, predicted_probability = cand$probability[j],
          matched_time_ma = matched_time,
          matched_cell_id = cand$cell_id[j],
          temporal_distance_myr = abs(matched_time - age_midpoint),
          matched_within_age_interval = matched_time >= ages[[1]] &&
            matched_time <= ages[[2]],
          spatial_distance_km = distance_km[j],
          hit = cand$probability[j] >= threshold)
  })
  do.call(rbind, out)
}

#' Build standardized output paths for deep-time products
#'
#' @param type Product type.
#' @param time_ma Age in Ma.
#' @param species_id Optional species id.
#' @param stat Statistic label.
#' @param ext File extension.
#' @param root Output root.
#' @return A file path.
#' @export
hee_output_path <- function(type,
                            time_ma,
                            species_id = NULL,
                            stat = "mean",
                            ext = "csv",
                            root = "outputs") {
  tstr <- sprintf("%05.1fMa", as.numeric(time_ma))
  safe <- function(x) gsub("[^A-Za-z0-9_\\-]+", "_", as.character(x))
  if (!is.null(species_id)) {
    file.path(root, type, paste0("species=", safe(species_id)),
              paste0("time=", tstr),
              paste0("stat=", safe(stat), ".", ext))
  } else {
    file.path(root, type, paste0(type, "_time=", tstr, "_stat=", safe(stat), ".", ext))
  }
}

#' Prepare a lightweight deep-time data object
#'
#' @param comm Community matrix.
#' @param env_now Modern environment.
#' @param env_cube Optional time cube.
#' @param tree Optional phylogeny.
#' @param traits Optional traits.
#' @param tip_ranges Optional range table.
#' @param fossils Optional fossil table.
#' @param response_type Response type.
#' @param strict Logical; stop on critical alignment failures.
#' @return A `hee_data` object.
#' @export
hee_prepare_data <- function(comm,
                             env_now,
                             env_cube = NULL,
                             tree = NULL,
                             traits = NULL,
                             tip_ranges = NULL,
                             fossils = NULL,
                             response_type = "presence",
                             strict = TRUE) {
  comm_check <- hee_check_comm(comm, response_type = response_type)
  env_check <- hee_check_env_now(env_now, rownames(comm))
  if (isTRUE(strict) && !isTRUE(env_check$ok)) {
    stop("env_now row names must match community site names. Missing site(s): ",
         paste(env_check$missing_sites, collapse = ", "), call. = FALSE)
  }
  tree_check <- if (!is.null(tree)) hee_check_tree(tree, colnames(comm)) else NULL
  trait_check <- if (!is.null(traits)) hee_check_traits(traits, colnames(comm)) else NULL
  range_check <- if (!is.null(tip_ranges)) hee_check_tip_ranges(tip_ranges, colnames(comm)) else NULL
  out <- list(comm = as.matrix(comm), env_now = as.data.frame(env_now),
              env_cube = if (!is.null(env_cube)) .hee_as_timecube(env_cube) else NULL,
              tree = tree, traits = traits, tip_ranges = tip_ranges,
              fossils = fossils,
              diagnostics = list(comm = comm_check, env_now = env_check,
                                 tree = tree_check, traits = trait_check,
                                 tip_ranges = range_check))
  class(out) <- "hee_data"
  out
}

#' Run a compact deep-time hindcasting pipeline
#'
#' @param evo Optional `hmsc_ecoevo` object.
#' @param beta Optional Beta matrix.
#' @param env_cube Time cube.
#' @param recipe Projection recipe.
#' @param times Ages in Ma.
#' @param variables Variables to project.
#' @param species_origin Optional species origin ages.
#' @param accessibility Optional accessibility weights.
#' @param link Link/response scale passed to [hee_project_hmsc_timecube()].
#' @return A list with suitability, projection, richness, turnover, and refugia.
#' @export
hee_run_deep_time_pipeline <- function(evo = NULL,
                                       beta = NULL,
                                       env_cube,
                                       recipe,
                                       times,
                                       variables,
                                       species_origin = NULL,
                                       accessibility = NULL,
                                       link = "probit") {
  suitability <- hee_project_hmsc_timecube(
    evo = evo, beta = beta, env_cube = env_cube, recipe = recipe,
    times = times, variables = variables, link = link
  )
  species <- unique(suitability$species)
  phylo_mask <- hee_phylo_time_mask(times = times, species_origin = species_origin,
                                    species = species)
  projection <- hee_combine(suitability, accessibility = accessibility,
                            phylo_mask = phylo_mask)
  richness <- hee_richness(projection)
  turnover <- hee_turnover(projection)
  refugia <- hee_refugia(richness)
  out <- list(suitability = suitability, projection = projection,
              richness = richness, turnover = turnover, refugia = refugia)
  class(out) <- "hee_deep_time_result"
  out
}

#' Compare BioGeoBEARS-style model table and model-average accessibility
#'
#' @param models Table with unique `model` names, `logLik`, and
#'   `n_parameters`. Duplicate candidate model names are rejected because they
#'   would duplicate AICc weights in model averaging.
#' @param n Number of observations for AICc. AICc is undefined when
#'   `n <= k + 1` for any candidate model; in that case the function stops
#'   instead of manufacturing model weights.
#' @return Model comparison table.
#' @export
hee_bgb_compare <- function(models, n) {
  models <- as.data.frame(models)
  .require_cols(models, c("model", "logLik", "n_parameters"), "models")
  if (any(is.na(models$model) | !nzchar(as.character(models$model)))) {
    stop("models$model must contain non-empty model names.", call. = FALSE)
  }
  dup_model <- duplicated(as.character(models$model))
  if (any(dup_model)) {
    stop("Duplicate BioGeoBEARS model name(s): ",
         paste(unique(as.character(models$model)[dup_model]), collapse = ", "),
         call. = FALSE)
  }
  n <- suppressWarnings(as.numeric(n)[1])
  if (!is.finite(n) || n <= 0) {
    stop("`n` must be a finite positive number of observations for AICc.",
         call. = FALSE)
  }
  k <- suppressWarnings(as.numeric(models$n_parameters))
  ll <- suppressWarnings(as.numeric(models$logLik))
  if (any(!is.finite(k) | k < 0 | !is.finite(ll))) {
    stop("models must contain finite logLik and non-negative n_parameters.",
         call. = FALSE)
  }
  bad <- models$model[n <= k + 1]
  if (length(bad) > 0L) {
    stop("AICc is undefined when n <= k + 1. Increase `n` or remove ",
         "over-parameterized model(s): ", paste(bad, collapse = ", "),
         call. = FALSE)
  }
  models$AIC <- -2 * models$logLik + 2 * k
  models$AICc <- models$AIC + (2 * k * (k + 1)) / (n - k - 1)
  models$delta_AICc <- models$AICc - min(models$AICc, na.rm = TRUE)
  w <- exp(-0.5 * models$delta_AICc)
  models$AICc_weight <- w / sum(w, na.rm = TRUE)
  models[order(models$AICc), , drop = FALSE]
}

#' @rdname hee_bgb_compare
#' @param null_model Null model row.
#' @param alternative_model Alternative model row.
#' @export
hee_bgb_lrt <- function(null_model, alternative_model) {
  n0 <- as.data.frame(null_model)
  n1 <- as.data.frame(alternative_model)
  .require_cols(n0, c("model", "logLik", "n_parameters"), "null_model")
  .require_cols(n1, c("model", "logLik", "n_parameters"), "alternative_model")
  valid_pairs <- c(
    "DEC" = "DEC+J",
    "DIVALIKE" = "DIVALIKE+J",
    "BAYAREALIKE" = "BAYAREALIKE+J"
  )
  m0 <- as.character(n0$model[1])
  m1 <- as.character(n1$model[1])
  if (m0 %in% names(valid_pairs) || m1 %in% unname(valid_pairs)) {
    if (!identical(unname(valid_pairs[m0]), m1)) {
      stop("BioGeoBEARS LRT is only valid for recognised nested pairs: ",
           "DEC vs DEC+J, DIVALIKE vs DIVALIKE+J, and BAYAREALIKE vs BAYAREALIKE+J.",
           call. = FALSE)
    }
  } else {
    warning("Unrecognised BioGeoBEARS model names; LRT assumes the supplied ",
            "models are genuinely nested.", call. = FALSE)
  }
  LR <- 2 * (n1$logLik[1] - n0$logLik[1])
  df <- n1$n_parameters[1] - n0$n_parameters[1]
  if (!is.finite(LR) || !is.finite(df) || df <= 0) {
    stop("BioGeoBEARS LRT requires finite log-likelihoods and a positive ",
         "parameter-count difference.", call. = FALSE)
  }
  if (LR < 0) {
    stop("BioGeoBEARS LRT alternative model has lower logLik than the null model.",
         call. = FALSE)
  }
  data.frame(null_model = n0$model[1], alternative_model = n1$model[1],
             LR = LR, df = df, p_value = stats::pchisq(LR, df = df, lower.tail = FALSE),
             stringsAsFactors = FALSE)
}

#' Extract one time slice from a deep-time environmental cube
#'
#' Public wrapper for the compact Phanerozoic time-cube accessor. Use this in
#' scripts and templates instead of private namespace helpers.
#'
#' @param env_cube A `hee_timecube`, compact RDS environment, or path accepted
#'   by [hee_load_timecube()].
#' @param time Age in Ma.
#' @param variables Variables to extract.
#' @return A data.frame for the requested time slice.
#' @export
#'
#' @examples
#' \dontrun{
#' cube <- hee_load_timecube("Phanerozoic_Environment_540_0Ma_5Myr_4deg_landharmonized_v2")
#' hee_get_time_slice(cube, 0, c("bio1", "bio12"))
#' }
hee_get_time_slice <- function(env_cube, time, variables = NULL) {
  cube <- .hee_as_timecube(env_cube)
  if (is.null(variables)) variables <- cube$variables
  .hee_get_slice(cube, time, variables)
}

# Internal helpers ----------------------------------------------------------

.hee_find_timecube_rds <- function(path) {
  if (file.exists(path) && !dir.exists(path)) return(path)
  if (!dir.exists(path)) stop("Cannot find path: ", path, call. = FALSE)
  roots <- c(file.path(path, "08_rds_export"), path)
  hits <- unique(unlist(lapply(roots, function(d) {
    if (dir.exists(d)) list.files(d, pattern = "\\.rds$", full.names = TRUE) else character()
  })))
  preferred <- hits[grepl("compact.*\\.rds$|float32.*\\.rds$", basename(hits), ignore.case = TRUE)]
  hit <- if (length(preferred) > 0) preferred[1] else hits[1]
  if (is.na(hit) || length(hit) == 0) stop("No RDS time cube found under: ", path, call. = FALSE)
  hit
}

.hee_validate_cube <- function(cube) {
  ok <- is.environment(cube) &&
    all(c("age_ma", "variable_names", "get_slice", "get_metadata") %in% ls(cube))
  if (!ok) stop("env_cube must be a compact Phanerozoic RDS environment with get_slice/get_metadata methods.",
                call. = FALSE)
  invisible(TRUE)
}

.hee_as_timecube <- function(x) {
  if (inherits(x, "hee_timecube")) return(x)
  if (is.character(x) && length(x) == 1) return(hee_load_timecube(x))
  if (is.environment(x)) return(hee_load_timecube(object = x))
  stop("Expected a hee_timecube, compact RDS environment, or path.", call. = FALSE)
}

.hee_match_times <- function(times, available) {
  times <- as.numeric(times)
  missing <- setdiff(times, available)
  if (length(missing) > 0) stop("Time cube missing age(s): ", paste(missing, collapse = ", "),
                                call. = FALSE)
  times
}

.hee_get_slice <- function(cube, time, variables) {
  sl <- cube$cube$get_slice(age = time, variables = variables)
  sl <- as.data.frame(sl)
  if ("age_ma" %in% names(sl) && !"time_ma" %in% names(sl)) {
    names(sl)[names(sl) == "age_ma"] <- "time_ma"
  }
  sl
}

.hee_cell_id <- function(lon, lat) {
  paste0("lon", sprintf("%+.3f", lon), "_lat", sprintf("%+.3f", lat))
}

.hee_link_inverse <- function(eta, link) {
  link <- tolower(as.character(link)[1])
  if (grepl("probit|bernoulli|binomial", link)) return(stats::pnorm(eta))
  if (grepl("logit", link)) return(stats::plogis(eta))
  if (grepl("poisson|count", link)) return(exp(eta))
  eta
}

.hee_reject_species_broadcast <- function(P, tab, table_name,
                                          allow_species_broadcast = FALSE) {
  if (isTRUE(allow_species_broadcast) || is.null(tab) || is.matrix(tab)) {
    return(invisible(TRUE))
  }
  tab <- as.data.frame(tab)
  if ("species" %in% names(tab) || !"species" %in% names(P)) {
    return(invisible(TRUE))
  }
  spp <- unique(as.character(P$species))
  spp <- spp[!is.na(spp)]
  if (length(spp) <= 1L) return(invisible(TRUE))
  stop(table_name, " was supplied without a `species` column while the ",
       "projection contains ", length(spp), " species. This would broadcast a ",
       "species-level process to every species. Add a `species` key, aggregate ",
       "the projection to a species-neutral level, or set ",
       "`allow_species_broadcast = TRUE` for an explicit neutral scenario.",
       call. = FALSE)
}

.hee_lookup_weight <- function(P, tab, value_col, default = 1,
                               aliases = character(), table_name = "table") {
  if (is.null(tab)) return(rep(default, nrow(P)))
  if (is.matrix(tab)) {
    return(.hee_lookup_matrix_weight(P, tab, default = default,
                                     table_name = table_name))
  }
  tab <- as.data.frame(tab)
  candidates <- unique(c(value_col, aliases))
  hit <- intersect(candidates, names(tab))
  if (length(hit) == 0L) {
    return(rep(default, nrow(P)))
  }
  if (!value_col %in% names(tab)) {
    tab[[value_col]] <- tab[[hit[1]]]
  }
  keys <- intersect(c("species", "cell_id", "time_ma", "region"), names(tab))
  keys <- keys[keys %in% names(P)]
  if (length(keys) == 0) {
    if (nrow(tab) == 1L) {
      out <- rep(tab[[value_col]][1], nrow(P))
      out[is.na(out)] <- default
      return(.hee_clip01(out))
    }
    return(rep(default, nrow(P)))
  }
  .hee_check_unique_keys(tab, keys, table_name)
  P$.lookup_row_id <- seq_len(nrow(P))
  z <- merge(P[, c(".lookup_row_id", keys), drop = FALSE],
             tab[, c(keys, value_col), drop = FALSE],
             by = keys, all.x = TRUE, sort = FALSE)
  if (nrow(z) != nrow(P)) {
    stop("Join from ", table_name, " changed row count from ", nrow(P),
         " to ", nrow(z), ". Check keys: ", paste(keys, collapse = ", "),
         call. = FALSE)
  }
  out <- rep(default, nrow(P))
  out[z$.lookup_row_id] <- z[[value_col]]
  out[is.na(out)] <- default
  .hee_clip01(out)
}

.hee_lookup_matrix_weight <- function(P, mat, default = 1,
                                      table_name = "matrix table") {
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop(table_name, " supplied as a matrix must have row and column names ",
         "so it can be matched to species/cell/region and time.", call. = FALSE)
  }
  out <- rep(default, nrow(P))
  row_keys <- intersect(c("species", "cell_id", "region"), names(P))
  row_key <- row_keys[vapply(row_keys, function(k) {
    all(as.character(P[[k]]) %in% rownames(mat))
  }, logical(1))]
  if (length(row_key) == 0L) {
    stop(table_name, " matrix row names do not match any available key in ",
         "the projection table (`species`, `cell_id`, or `region`).",
         call. = FALSE)
  }
  row_key <- row_key[1]

  col_key <- NULL
  col_index <- rep(NA_integer_, nrow(P))
  if ("time_ma" %in% names(P)) {
    mat_time <- suppressWarnings(as.numeric(gsub("Ma$", "", colnames(mat))))
    if (all(is.finite(mat_time)) && any(as.numeric(P$time_ma) %in% mat_time)) {
      col_key <- "time_ma"
      col_index <- match(as.numeric(P$time_ma), mat_time)
    }
  }
  if (is.null(col_key)) {
    other_cols <- setdiff(intersect(c("region", "cell_id", "species"), names(P)),
                          row_key)
    hit <- other_cols[vapply(other_cols, function(k) {
      all(as.character(P[[k]]) %in% colnames(mat))
    }, logical(1))]
    if (length(hit) > 0L) {
      col_key <- hit[1]
      col_index <- match(as.character(P[[col_key]]), colnames(mat))
    }
  }
  if (is.null(col_key)) {
    stop(table_name, " matrix column names do not match `time_ma` or another ",
         "available key in the projection table.", call. = FALSE)
  }

  row_index <- match(as.character(P[[row_key]]), rownames(mat))
  ok <- !is.na(row_index) & !is.na(col_index)
  out[ok] <- mat[cbind(row_index[ok], col_index[ok])]
  out[is.na(out)] <- default
  .hee_clip01(out)
}

.hee_require_any_columns <- function(tab, cols, table_name = "table") {
  if (is.null(tab) || is.matrix(tab)) return(invisible(TRUE))
  tab <- as.data.frame(tab)
  if (length(intersect(cols, names(tab))) == 0L) {
    stop(table_name, " was supplied but contains none of the recognised value columns: ",
         paste(cols, collapse = ", "), ".", call. = FALSE)
  }
  invisible(TRUE)
}

.hee_lookup_product <- function(P, tab, cols, default = 1, table_name = "table") {
  if (is.null(tab)) return(rep(default, nrow(P)))
  tab <- as.data.frame(tab)
  present <- intersect(cols, names(tab))
  if (length(present) == 0L) return(rep(default, nrow(P)))
  weights <- lapply(present, function(nm) {
    .hee_lookup_weight(P, tab, nm, default = default, table_name = table_name)
  })
  .hee_clip01(Reduce(`*`, weights))
}

.hee_lookup_phylo <- function(P, phylo_mask) {
  if (is.null(phylo_mask)) return(rep(1, nrow(P)))
  if (is.matrix(phylo_mask)) {
    times <- attr(phylo_mask, "times")
    if (is.null(times)) times <- as.numeric(gsub("Ma$", "", colnames(phylo_mask)))
    out <- rep(0, nrow(P))
    sp_i <- match(as.character(P$species), rownames(phylo_mask))
    tm_i <- match(P$time_ma, times)
    ok <- !is.na(sp_i) & !is.na(tm_i)
    out[ok] <- phylo_mask[cbind(sp_i[ok], tm_i[ok])]
    return(out)
  }
  tab <- as.data.frame(phylo_mask)
  .hee_lookup_weight(P, tab, "phylo_existence", default = 0,
                     aliases = c("E_phylo", "phylo_present", "lineage_exists"),
                     table_name = "phylo_mask")
}

.hee_rows_matched <- function(P, tab, table_name = "table") {
  if (is.null(tab)) return(rep(FALSE, nrow(P)))
  if (is.matrix(tab)) {
    if (is.null(rownames(tab)) || is.null(colnames(tab))) return(rep(FALSE, nrow(P)))
    row_keys <- intersect(c("species", "cell_id", "region"), names(P))
    row_key <- row_keys[vapply(row_keys, function(k) {
      all(as.character(P[[k]]) %in% rownames(tab))
    }, logical(1))]
    if (length(row_key) == 0L) return(rep(FALSE, nrow(P)))
    row_key <- row_key[1]
    col_index <- rep(NA_integer_, nrow(P))
    if ("time_ma" %in% names(P)) {
      mat_time <- suppressWarnings(as.numeric(gsub("Ma$", "", colnames(tab))))
      if (all(is.finite(mat_time))) col_index <- match(as.numeric(P$time_ma), mat_time)
    }
    if (all(is.na(col_index))) {
      other_cols <- setdiff(intersect(c("region", "cell_id", "species"), names(P)),
                            row_key)
      hit <- other_cols[vapply(other_cols, function(k) {
        all(as.character(P[[k]]) %in% colnames(tab))
      }, logical(1))]
      if (length(hit) > 0L) col_index <- match(as.character(P[[hit[1]]]), colnames(tab))
    }
    row_index <- match(as.character(P[[row_key]]), rownames(tab))
    return(!is.na(row_index) & !is.na(col_index))
  }
  tab <- as.data.frame(tab)
  keys <- intersect(c("species", "cell_id", "time_ma", "region"), names(tab))
  keys <- keys[keys %in% names(P)]
  if (length(keys) == 0L) return(rep(nrow(tab) == 1L, nrow(P)))
  .hee_check_unique_keys(tab, keys, table_name)
  P$.lookup_row_id <- seq_len(nrow(P))
  key_tab <- unique(tab[, keys, drop = FALSE])
  key_tab$.matched <- TRUE
  z <- merge(P[, c(".lookup_row_id", keys), drop = FALSE],
             key_tab,
             by = keys, all.x = TRUE, sort = FALSE)
  if (nrow(z) != nrow(P)) {
    stop("Join from ", table_name, " changed row count from ", nrow(P),
         " to ", nrow(z), ". Check keys: ", paste(keys, collapse = ", "),
         call. = FALSE)
  }
  matched <- rep(FALSE, nrow(P))
  matched[z$.lookup_row_id] <- z$.matched %in% TRUE
  matched
}
