# HmscEE nested Dynamic Earth-Biota Assembly core ---------------------------

.hee_stop_obsolete_product_core <- function(fun) {
  stop(
    fun, "() belonged to the obsolete multiplicative HmscEE core ",
    "`S_HMSC * A_hist * E_phylo * L_land * D_dynamic * Q_extrap`. ",
    "That equation is no longer the main model because BioGeoBEARS/BSM ",
    "regional history and HmscEE within-region movement can otherwise explain ",
    "the same dispersal process twice. Use `hee_nested_region_cell_occupancy()` ",
    "with BioGeoBEARS/BSM `region_history`, `hee_within_region_movement_kernel()`, ",
    "`hee_arrival_pressure()`, `hee_establishment_probability()`, and ",
    "`hee_persistence_probability()` instead.",
    call. = FALSE
  )
}

.hee_hmscee_key <- function(...) paste(..., sep = "\r")

.hee_hmscee_first_col <- function(x, candidates, label) {
  hit <- intersect(candidates, names(x))
  if (length(hit) == 0L) {
    stop(label, " must contain one of: ", paste(candidates, collapse = ", "),
         call. = FALSE)
  }
  hit[[1L]]
}

.hee_hmscee_prob <- function(x) {
  out <- .hee_clip01(x)
  out[!is.finite(out) | is.na(out)] <- 0
  out
}

.hee_hmscee_normalize_weights <- function(w) {
  w <- suppressWarnings(as.numeric(w))
  w[!is.finite(w) | is.na(w) | w < 0] <- 0
  s <- sum(w)
  if (!is.finite(s) || s <= 0) rep(1 / length(w), length(w)) else w / s
}

.hee_hmscee_weighted_quantile <- function(x, w, probs) {
  x <- suppressWarnings(as.numeric(x))
  w <- suppressWarnings(as.numeric(w))
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(rep(NA_real_, length(probs)))
  x <- x[ok]
  w <- .hee_hmscee_normalize_weights(w[ok])
  ord <- order(x)
  x <- x[ord]
  cw <- cumsum(w[ord])
  vapply(probs, function(p) x[which(cw >= p)[1L]], numeric(1))
}

.hee_hmscee_pick_scenario_col <- function(x, canonical, legacy = character(),
                                          explicit = NULL) {
  if (!is.null(explicit) && explicit %in% names(x)) return(explicit)
  hit <- intersect(c(canonical, legacy), names(x))
  if (length(hit)) hit[[1L]] else NULL
}

.hee_hmscee_normalize_group_weight <- function(keys, weight) {
  w <- suppressWarnings(as.numeric(weight))
  w[!is.finite(w) | is.na(w) | w < 0] <- 0
  stats::ave(w, keys, FUN = .hee_hmscee_normalize_weights)
}

#' Define dynamic palaeo-Earth state variables for HmscEE
#'
#' HmscEE treats geology as state variables, not as independent biological
#' probability multipliers. This helper standardizes the cell-time palaeo-Earth
#' table into the three objects used by the nested model:
#' \deqn{E_k = \{X_k, H_k, W_k\}.}
#' `X` is the supplied environmental table, `H_state` is a hard cell existence
#' and usable-habitat state, and movement resistance or least-cost distances are
#' handled separately by `hee_within_region_movement_kernel()`. Supplied `NA`
#' values in hard-state columns are treated as absent geography (`0`).
#'
#' @param earth Cell-time table with `cell_id`, `time_ma`, and `region`.
#' @param habitat_col Optional usable habitat column.
#' @param land_col Land or geographic existence column. If absent but
#'   `geographic_existence` or `H_state` exists, that column is used.
#' @param region_col,cell_col,time_col Column names.
#' @param env_cols Optional environmental columns to keep in a stable order.
#' @param earth_scenario_col Legacy scenario column. Prefer
#'   `geography_scenario_col`; the legacy value is interpreted as a geography
#'   scenario for backward compatibility.
#' @param geography_scenario_col Optional palaeogeographic/plate scenario
#'   column `g`. If absent, all rows are assigned to `"g1"`.
#' @param climate_scenario_col Optional palaeoclimate scenario column `c`. If
#'   absent, all rows are assigned to `"c1"`.
#' @param climate_weight_col Optional conditional climate-scenario weight
#'   column. Weights are normalized within each `geography_scenario`, i.e.
#'   `w[c | g]`.
#'
#' @return A cell-time table with `H_state`, `land_state`, `habitat_state`, and
#'   optional environmental columns. `H_state` is a hard constraint: if it is
#'   zero, occupancy is forced to zero.
#' @export
#'
#' @examples
#' x <- data.frame(cell_id = "c1", time_ma = 0, region = "R1",
#'   land = 1, bio1 = 20)
#' hee_paleo_earth_state(x, env_cols = "bio1")
hee_paleo_earth_state <- function(earth,
                                  habitat_col = NULL,
                                  land_col = NULL,
                                  region_col = "region",
                                  cell_col = "cell_id",
                                  time_col = "time_ma",
                                  earth_scenario_col = NULL,
                                  env_cols = NULL,
                                  geography_scenario_col = NULL,
                                  climate_scenario_col = NULL,
                                  climate_weight_col = NULL) {
  x <- as.data.frame(earth)
  .require_cols(x, c(cell_col, time_col, region_col), "earth")
  .hee_validate_time_values(x[[time_col]], paste0("earth$", time_col))
  geography_scenario_col <- .hee_hmscee_pick_scenario_col(
    x, "geography_scenario", legacy = "earth_scenario",
    explicit = geography_scenario_col
  )
  if (is.null(geography_scenario_col) && !is.null(earth_scenario_col) &&
      earth_scenario_col %in% names(x)) {
    geography_scenario_col <- earth_scenario_col
  }
  if (is.null(geography_scenario_col) ||
      !geography_scenario_col %in% names(x)) {
    x$geography_scenario <- "g1"
    geography_scenario_col <- "geography_scenario"
  }
  climate_scenario_col <- .hee_hmscee_pick_scenario_col(
    x, "climate_scenario", explicit = climate_scenario_col
  )
  if (is.null(climate_scenario_col) ||
      !climate_scenario_col %in% names(x)) {
    x$climate_scenario <- "c1"
    climate_scenario_col <- "climate_scenario"
  }
  if (is.null(climate_weight_col)) {
    climate_weight_col <- intersect(c("climate_weight_given_geography",
                                      "climate_weight", "climate_scenario_weight"),
                                    names(x))[1L]
  }
  if (is.na(climate_weight_col) || !climate_weight_col %in% names(x)) {
    x$climate_weight_given_geography <- 1
    climate_weight_col <- "climate_weight_given_geography"
  }
  .hee_check_unique_keys(x, c(geography_scenario_col, climate_scenario_col,
                              cell_col, time_col), "earth")
  if (is.null(land_col)) {
    land_col <- intersect(c("H_state", "land", "land_mask", "L_land",
                            "geographic_existence", "habitat_available"),
                          names(x))[1L]
  }
  if (is.na(land_col) || !land_col %in% names(x)) {
    stop("earth must contain a land/geographic hard-state column such as ",
         "`land`, `land_mask`, `geographic_existence`, or `H_state`.",
         call. = FALSE)
  }
  if (!is.null(habitat_col)) .require_cols(x, habitat_col, "earth")
  if (!is.null(env_cols)) .require_cols(x, env_cols, "earth")
  land <- .hee_hmscee_prob(x[[land_col]])
  habitat <- if (!is.null(habitat_col)) .hee_hmscee_prob(x[[habitat_col]]) else
    rep(1, nrow(x))
  x$land_state <- as.integer(land > 0)
  x$habitat_state <- as.integer(habitat > 0)
  x$H_state <- as.integer(x$land_state > 0 & x$habitat_state > 0)
  x$region <- as.character(x[[region_col]])
  x$cell_id <- as.character(x[[cell_col]])
  x$time_ma <- as.numeric(x[[time_col]])
  x$geography_scenario <- as.character(x[[geography_scenario_col]])
  x$climate_scenario <- as.character(x[[climate_scenario_col]])
  cw <- unique(x[, c("geography_scenario", "climate_scenario",
                    climate_weight_col), drop = FALSE])
  names(cw)[names(cw) == climate_weight_col] <-
    "climate_weight_given_geography"
  .hee_check_unique_keys(cw, c("geography_scenario", "climate_scenario"),
                         "climate scenario weights")
  cw$climate_weight_given_geography <- .hee_hmscee_normalize_group_weight(
    cw$geography_scenario, cw$climate_weight_given_geography
  )
  x$climate_weight_given_geography <- cw$climate_weight_given_geography[
    match(.hee_hmscee_key(x$geography_scenario, x$climate_scenario),
          .hee_hmscee_key(cw$geography_scenario, cw$climate_scenario))
  ]
  x$earth_scenario <- x$geography_scenario
  keep <- unique(c(geography_scenario_col, climate_scenario_col,
                   cell_col, time_col, region_col,
                   "geography_scenario", "climate_scenario", "earth_scenario",
                   "climate_weight_given_geography", "cell_id", "time_ma",
                   "region", "H_state", "land_state", "habitat_state",
                   env_cols, setdiff(names(x), c(cell_col, time_col,
                                                 region_col,
                                                 geography_scenario_col,
                                                 climate_scenario_col))))
  x[, intersect(keep, names(x)), drop = FALSE]
}

#' Validate BioGeoBEARS/BSM regional history draws
#'
#' BioGeoBEARS or an equivalent historical biogeographic model provides the
#' macro-scale state
#' \deqn{R_{\ell,r,k}^{(h)} \in \{0,1\},}
#' indicating whether lineage `lineage` occupies biogeographic region `region`
#' at time `time_ma` in BSM/history draw `history_draw`, conditional on
#' palaeogeographic/plate scenario `geography_scenario` (`g`). This avoids
#' using `e` for Earth scenario, because BioGeoBEARS already uses parameter
#' `e` for range contraction/local extinction. HmscEE does not turn this into
#' a continuous accessibility multiplier; it conditions cell-scale occupancy
#' on this regional state.
#'
#' @param region_history Table with lineage, region, time, and occupancy state.
#' @param state_col Regional occupancy column. If absent, common aliases
#'   `R_region`, `region_occupied`, `occupied`, and `accessibility` are tried.
#' @param history_col Optional BSM/history draw id column. If absent, a single
#'   draw `"h1"` is used.
#' @param earth_scenario_col Legacy scenario column. Prefer
#'   `geography_scenario_col`; the legacy value is interpreted as geography
#'   scenario `g` for backward compatibility.
#' @param geography_scenario_col Optional palaeogeographic/plate scenario
#'   column `g`. If absent, all histories are assigned to `"g1"`.
#' @param history_weight_col Optional conditional history weight column. The
#'   weights are normalized within each `geography_scenario`, i.e.
#'   `w[h | g]`, not globally.
#' @param earth_weight_col Legacy scenario-weight column.
#' @param geography_weight_col Optional geography-scenario weight column
#'   `w[g]`. Weights are normalized across unique `geography_scenario`.
#' @param lineage_col,region_col,time_col Column names.
#'
#' @return A validated table with `lineage`, `region`, `time_ma`,
#'   `geography_scenario`, `history_draw`, conditional weights, and binary
#'   `R_region`. Legacy `earth_scenario` aliases are retained in the output for
#'   old scripts.
#' @export
#'
#' @examples
#' r <- data.frame(lineage = "sp1", region = "R1", time_ma = 0,
#'   R_region = 1)
#' hee_bsm_region_history(r)
hee_bsm_region_history <- function(region_history,
                                   state_col = NULL,
                                   history_col = NULL,
                                   earth_scenario_col = NULL,
                                   history_weight_col = NULL,
                                   earth_weight_col = NULL,
                                   geography_scenario_col = NULL,
                                   geography_weight_col = NULL,
                                   lineage_col = "lineage",
                                   region_col = "region",
                                   time_col = "time_ma") {
  x <- as.data.frame(region_history)
  .require_cols(x, c(lineage_col, region_col, time_col), "region_history")
  .hee_validate_time_values(x[[time_col]], paste0("region_history$", time_col))
  if (is.null(state_col)) {
    state_col <- .hee_hmscee_first_col(
      x, c("R_region", "region_occupied", "occupied", "accessibility"),
      "region_history"
    )
  } else {
    .require_cols(x, state_col, "region_history")
  }
  if (is.null(history_col) && "history_draw" %in% names(x)) {
    history_col <- "history_draw"
  }
  geography_scenario_col <- .hee_hmscee_pick_scenario_col(
    x, "geography_scenario", legacy = "earth_scenario",
    explicit = geography_scenario_col
  )
  if (is.null(geography_scenario_col) && !is.null(earth_scenario_col) &&
      earth_scenario_col %in% names(x)) {
    geography_scenario_col <- earth_scenario_col
  }
  if (is.null(history_weight_col) &&
      "history_weight_given_geography" %in% names(x)) {
    history_weight_col <- "history_weight_given_geography"
  }
  if (is.null(history_weight_col) &&
      "history_weight_given_earth" %in% names(x)) {
    history_weight_col <- "history_weight_given_earth"
  }
  if (is.null(history_weight_col) && "history_weight" %in% names(x)) {
    history_weight_col <- "history_weight"
  }
  if (is.null(geography_weight_col) && "geography_weight" %in% names(x)) {
    geography_weight_col <- "geography_weight"
  }
  if (is.null(earth_weight_col) && "earth_weight" %in% names(x)) {
    earth_weight_col <- "earth_weight"
  }
  if (is.null(history_col) || !history_col %in% names(x)) {
    x$history_draw <- "h1"
    history_col <- "history_draw"
  }
  if (is.null(geography_scenario_col) ||
      !geography_scenario_col %in% names(x)) {
    x$geography_scenario <- "g1"
    geography_scenario_col <- "geography_scenario"
  }
  if (is.null(history_weight_col) || !history_weight_col %in% names(x)) {
    x$history_weight_given_geography <- 1
    history_weight_col <- "history_weight_given_geography"
  }
  if ((is.null(geography_weight_col) ||
       !geography_weight_col %in% names(x)) &&
      !is.null(earth_weight_col) && earth_weight_col %in% names(x)) {
    geography_weight_col <- earth_weight_col
  }
  if (is.null(geography_weight_col) || !geography_weight_col %in% names(x)) {
    x$geography_weight <- 1
    geography_weight_col <- "geography_weight"
  }
  out <- data.frame(
    geography_scenario = as.character(x[[geography_scenario_col]]),
    lineage = as.character(x[[lineage_col]]),
    region = as.character(x[[region_col]]),
    time_ma = as.numeric(x[[time_col]]),
    history_draw = as.character(x[[history_col]]),
    geography_weight = suppressWarnings(as.numeric(x[[geography_weight_col]])),
    history_weight_given_geography = suppressWarnings(as.numeric(x[[history_weight_col]])),
    R_region = as.integer(.hee_hmscee_prob(x[[state_col]]) > 0),
    stringsAsFactors = FALSE
  )
  .hee_check_unique_keys(out, c("geography_scenario", "history_draw", "lineage",
                                "region", "time_ma"), "region_history")
  gw <- unique(out[, c("geography_scenario", "geography_weight"),
                   drop = FALSE])
  .hee_check_unique_keys(gw, "geography_scenario", "geography scenario weights")
  gw$geography_weight <- .hee_hmscee_normalize_weights(gw$geography_weight)
  hw <- unique(out[, c("geography_scenario", "history_draw",
                       "history_weight_given_geography"), drop = FALSE])
  .hee_check_unique_keys(hw, c("geography_scenario", "history_draw"),
                         "history weights conditional on geography")
  hw$history_weight_given_geography <- .hee_hmscee_normalize_group_weight(
    hw$geography_scenario, hw$history_weight_given_geography
  )
  out$geography_weight <- gw$geography_weight[
    match(out$geography_scenario, gw$geography_scenario)
  ]
  out$history_weight_given_geography <- hw$history_weight_given_geography[
    match(.hee_hmscee_key(out$geography_scenario, out$history_draw),
          .hee_hmscee_key(hw$geography_scenario, hw$history_draw))
  ]
  out$earth_scenario <- out$geography_scenario
  out$earth_weight <- out$geography_weight
  out$history_weight_given_earth <- out$history_weight_given_geography
  out[order(out$geography_scenario, out$history_draw, out$lineage, out$region,
            -out$time_ma),
      , drop = FALSE]
}

.hee_response_tree_geometry <- function(tree) {
  .require_pkg("ape", "ancestral environmental-response reconstruction")
  if (!inherits(tree, "phylo")) {
    stop("tree must be an ape::phylo object.", call. = FALSE)
  }
  if (is.null(tree$edge.length)) {
    stop("tree must be dated and contain branch lengths.", call. = FALSE)
  }
  if (any(!is.finite(tree$edge.length) | tree$edge.length < 0)) {
    stop("tree branch lengths must be finite and non-negative.", call. = FALSE)
  }
  tree <- ape::reorder.phylo(tree, "cladewise")
  n_tip <- length(tree$tip.label)
  node_ids <- seq_len(n_tip + tree$Nnode)
  depth <- ape::node.depth.edgelength(tree)
  root <- setdiff(tree$edge[, 1], tree$edge[, 2])
  root <- root[1L]
  root_age <- max(depth[seq_len(n_tip)], na.rm = TRUE)
  node_age <- root_age - depth
  labels <- rep(NA_character_, length(node_ids))
  labels[seq_len(n_tip)] <- tree$tip.label
  internal <- (n_tip + 1L):(n_tip + tree$Nnode)
  node_labels <- tree$node.label
  if (is.null(node_labels)) node_labels <- rep("", tree$Nnode)
  labels[internal] <- ifelse(nzchar(node_labels),
                             node_labels,
                             paste0("node_", internal))
  list(tree = tree, n_tip = n_tip, root = root, root_age = root_age,
       node_age = node_age, labels = labels)
}

.hee_tip_beta_draw_table <- function(tip_beta_draws,
                                     basis_cols = NULL,
                                     species_col = "species",
                                     response_draw_col = "response_draw",
                                     intercept_axis = "(Intercept)",
                                     intercept_col = NULL) {
  if (is.array(tip_beta_draws)) {
    if (length(dim(tip_beta_draws)) != 3L) {
      stop("tip_beta_draws must be a 3D array with dimensions draw x species x axis, ",
           "or a data.frame with species, draw and beta-axis columns.",
           call. = FALSE)
    }
    dn <- dimnames(tip_beta_draws)
    d2 <- dn[[2]] %||% character()
    d3 <- dn[[3]] %||% character()
    if (length(d2) == 0L || length(d3) == 0L ||
        any(!nzchar(d2)) || any(!nzchar(d3))) {
      stop("tip_beta_draws array must have species and beta-axis dimnames.",
           call. = FALSE)
    }
    if (is.null(basis_cols)) {
      basis_cols <- setdiff(d3, intercept_axis)
    }
    if (!all(basis_cols %in% d3) && all(basis_cols %in% d2)) {
      tip_beta_draws <- aperm(tip_beta_draws, c(1, 3, 2))
      dn <- dimnames(tip_beta_draws)
      d2 <- dn[[2]]
      d3 <- dn[[3]]
    }
    missing_axes <- setdiff(basis_cols, d3)
    if (length(missing_axes) > 0L) {
      stop("tip_beta_draws is missing beta response axis/axes: ",
           paste(missing_axes, collapse = ", "), call. = FALSE)
    }
    draw_ids <- dn[[1]] %||% paste0("s", seq_len(dim(tip_beta_draws)[1]))
    if (length(draw_ids) != dim(tip_beta_draws)[1]) {
      draw_ids <- paste0("s", seq_len(dim(tip_beta_draws)[1]))
    }
    rows <- vector("list", dim(tip_beta_draws)[1] * dim(tip_beta_draws)[2])
    k <- 1L
    for (i in seq_len(dim(tip_beta_draws)[1])) {
      for (j in seq_len(dim(tip_beta_draws)[2])) {
        row <- data.frame(
          species = d2[[j]],
          response_draw = as.character(draw_ids[[i]]),
          stringsAsFactors = FALSE
        )
        for (ax in basis_cols) row[[ax]] <- tip_beta_draws[i, j, ax]
        if (intercept_axis %in% d3) {
          row$.tip_intercept <- tip_beta_draws[i, j, intercept_axis]
        }
        rows[[k]] <- row
        k <- k + 1L
      }
    }
    out <- do.call(rbind, rows)
    attr(out, "basis_cols") <- basis_cols
    return(out)
  }

  x <- as.data.frame(tip_beta_draws)
  if (!species_col %in% names(x)) {
    stop("tip_beta_draws data.frame must contain a `", species_col,
         "` column with modern tip species names. Modern suitability maps ",
         "are not valid input for ancestral response reconstruction.",
         call. = FALSE)
  }
  if (is.null(basis_cols)) {
    ids <- c(species_col, response_draw_col, "cell_id", "time_ma",
             "lon", "lat", "suitability", "probability")
    basis_cols <- setdiff(names(x)[vapply(x, is.numeric, logical(1))], ids)
  }
  .require_cols(x, c(species_col, basis_cols), "tip_beta_draws")
  if (!response_draw_col %in% names(x)) x[[response_draw_col]] <- "s1"
  out <- data.frame(
    species = as.character(x[[species_col]]),
    response_draw = as.character(x[[response_draw_col]]),
    stringsAsFactors = FALSE
  )
  for (ax in basis_cols) out[[ax]] <- suppressWarnings(as.numeric(x[[ax]]))
  if (is.null(intercept_col)) {
    hit <- intersect(c("intercept", intercept_axis, "(Intercept)"), names(x))
    intercept_col <- hit[1L]
  }
  if (!is.na(intercept_col) && intercept_col %in% names(x)) {
    out$.tip_intercept <- suppressWarnings(as.numeric(x[[intercept_col]]))
  }
  .hee_check_unique_keys(out, c("species", "response_draw"),
                         "tip_beta_draws")
  attr(out, "basis_cols") <- basis_cols
  out
}

#' Reconstruct ancestral environmental-response coefficients from HMSC Beta
#'
#' This is the direct-response evolution path for HmscEE. It reconstructs
#' environmental-response coefficients, not modern prediction maps:
#' \deqn{\beta_{\ell,t}^{(s)} \sim p(\beta_{\ell,t}\mid
#' \beta_{tips}^{(s)}, Tree).}
#' For each HMSC posterior draw `s`, the function treats each Beta axis as a
#' continuous response trait on a dated phylogeny, estimates internal-node
#' conditional means under BM, OU or early-burst (EB) evolution (including
#' polytomies), and predicts active branch responses at requested `time_ma`
#' values. The optional AICc mode compares these models separately for each
#' environmental-response axis and HMSC posterior draw.
#' Larger `time_ma` values are older.
#'
#' By default `reconstruct_intercept = FALSE`, so HMSC prevalence/sampling
#' intercepts are not carried into ancestral niche interpretation. The output
#' can be passed to [hee_lineage_suitability()] or the convenience wrapper
#' [hee_environmental_filtering_ancestral_suitability()]. BioGeoBEARS/BSM
#' ranges are deliberately not used here; ancestral geographic range and
#' ancestral environmental response are different objects.
#'
#' @param tip_beta_draws HMSC Beta posterior draws as a 3D array
#'   `draw x species x axis` (or `draw x axis x species` when axis names make
#'   orientation clear), or a data.frame with modern tip `species`,
#'   `response_draw`, and Beta-axis columns. Do not pass modern prediction
#'   maps.
#' @param tree Dated `ape::phylo` tree with branch lengths in Ma.
#' @param times Ages in Ma where active-lineage responses are needed. Larger
#'   values are older. Times older than the tree root are omitted with a warning.
#' @param basis_cols Beta/environmental basis columns to reconstruct. If
#'   `NULL`, all non-intercept axes are used.
#' @param species_col,response_draw_col Column names for data.frame input.
#' @param intercept_axis,intercept_col Intercept names. Intercepts are ignored
#'   unless `reconstruct_intercept = TRUE`.
#' @param reconstruct_intercept Logical. Reconstruct the HMSC intercept as an
#'   ancestral response axis. Default `FALSE` is recommended because HMSC
#'   intercepts include prevalence, sampling and accessibility signals.
#' @param evolution_model `"BM"` (default), `"OU"`, `"EB"`, or `"AICc"`.
#'   The last fits all three and selects the lowest AICc per response axis and
#'   posterior draw. OU has one shared optimum-like mean and attraction rate;
#'   EB has an exponentially declining rate from the tree root. Model support
#'   does not establish stabilizing selection or an adaptive burst. Fits at a
#'   parameter bound are flagged in `model_comparison`; AICc on posterior Beta
#'   draws is a sensitivity diagnostic, not a joint model probability.
#' @param branch_uncertainty `"none"` returns conditional mean branch
#'   trajectories. `"process_bridge"` draws a time-coherent Gaussian bridge
#'   under the selected BM, OU or EB process. `"brownian_bridge"` remains an
#'   alias for older scripts. Node states are conditional means rather than
#'   sampled node posteriors, so this is not a full joint model.
#' @param seed Optional random seed used only when
#'   `branch_uncertainty` requests a stochastic bridge.
#' @param drop_tree_extra Logical; drop tree tips not present in Beta draws.
#' @param time_tolerance Numeric tolerance in Ma for matching branch endpoints.
#'   The default (`1e-5`) treats tips in near-ultrametric dated trees as present
#'   at 0 Ma despite harmless floating-point or file-rounding differences.
#'
#' @return A data.frame with `lineage`, `lineage_type`, `node`, `parent_node`,
#'   `time_ma`, `response_draw`, `intercept`, and one column per Beta axis.
#'   Attribute `model_comparison` contains per-axis fit diagnostics and AICc.
#' @export
#'
#' @examples
#' \dontrun{
#' evo <- as_hmsc_ecoevo(hmsc_model)
#' response <- hee_evolution_ancestral_response_direct(
#'   evo$beta_draws, tree = evo$phylo, times = c(100, 50, 0),
#'   basis_cols = c("bio1", "bio12")
#' )
#' }
hee_evolution_ancestral_response_direct <- function(tip_beta_draws,
                                                    tree,
                                                    times,
                                                    basis_cols = NULL,
                                                    species_col = "species",
                                                    response_draw_col = "response_draw",
                                                    intercept_axis = "(Intercept)",
                                                    intercept_col = NULL,
                                                    reconstruct_intercept = FALSE,
                                                    branch_uncertainty = c("none", "process_bridge", "brownian_bridge"),
                                                    seed = NULL,
                                                    drop_tree_extra = TRUE,
                                                    time_tolerance = 1e-5,
                                                    evolution_model = c("BM", "OU", "EB", "AICc")) {
  branch_uncertainty <- match.arg(branch_uncertainty)
  evolution_model <- match.arg(evolution_model)
  time_tolerance <- suppressWarnings(as.numeric(time_tolerance)[1L])
  if (!is.finite(time_tolerance) || time_tolerance < 0) {
    stop("time_tolerance must be a non-negative finite number.", call. = FALSE)
  }
  geom <- .hee_response_tree_geometry(tree)
  tip_tab <- .hee_tip_beta_draw_table(
    tip_beta_draws,
    basis_cols = basis_cols,
    species_col = species_col,
    response_draw_col = response_draw_col,
    intercept_axis = intercept_axis,
    intercept_col = intercept_col
  )
  basis_cols <- attr(tip_tab, "basis_cols")
  species <- unique(tip_tab$species)
  missing_in_tree <- setdiff(species, geom$tree$tip.label)
  if (length(missing_in_tree) > 0L) {
    stop("tree is missing tip species from tip_beta_draws: ",
         paste(missing_in_tree, collapse = ", "), call. = FALSE)
  }
  extra_tips <- setdiff(geom$tree$tip.label, species)
  if (length(extra_tips) > 0L) {
    if (!isTRUE(drop_tree_extra)) {
      stop("tree contains tips not present in tip_beta_draws: ",
           paste(extra_tips, collapse = ", "), call. = FALSE)
    }
    geom <- .hee_response_tree_geometry(ape::drop.tip(geom$tree, extra_tips))
  }
  times <- sort(unique(suppressWarnings(as.numeric(times))), decreasing = TRUE)
  .hee_validate_time_values(times, "times")
  old <- times > geom$root_age + time_tolerance
  if (any(old)) {
    warning("Omitting ", sum(old), " time value(s) older than the tree root (",
            round(geom$root_age, 6), " Ma).", call. = FALSE)
    times <- times[!old]
  }
  if (length(times) == 0L) {
    stop("No requested times fall inside the dated tree temporal domain.",
         call. = FALSE)
  }
  if (!is.null(seed) && branch_uncertainty != "none") {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
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

  axes <- basis_cols
  if (isTRUE(reconstruct_intercept) && ".tip_intercept" %in% names(tip_tab)) {
    tip_tab$.reconstruct_intercept <- tip_tab$.tip_intercept
    axes_to_fit <- c(".reconstruct_intercept", axes)
  } else {
    axes_to_fit <- axes
  }
  draw_ids <- unique(tip_tab$response_draw)
  rows <- list()
  fit_rows <- list()
  rr <- 0L
  edge <- geom$tree$edge
  n_node_total <- length(geom$labels)
  if (geom$n_tip < 2L) {
    stop("Ancestral response evolution needs at least two tree tips.",
         call. = FALSE)
  }
  node_ids <- seq.int(geom$n_tip + 1L, n_node_total)
  branch_time_index <- lapply(seq_len(nrow(edge)), function(ei) {
    which(times >= geom$node_age[edge[ei, 2L]] - time_tolerance &
          times < geom$node_age[edge[ei, 1L]] - time_tolerance)
  })
  for (draw in draw_ids) {
    dd <- tip_tab[tip_tab$response_draw == draw, , drop = FALSE]
    .hee_check_unique_keys(dd, "species", "tip_beta_draws within response_draw")
    dd <- dd[match(geom$tree$tip.label, dd$species), , drop = FALSE]
    values <- matrix(NA_real_, nrow = n_node_total, ncol = length(axes_to_fit),
                     dimnames = list(as.character(seq_len(n_node_total)),
                                     axes_to_fit))
    axis_fits <- vector("list", length(axes_to_fit))
    names(axis_fits) <- axes_to_fit
    for (ax in axes_to_fit) {
      y <- suppressWarnings(as.numeric(dd[[ax]]))
      if (any(!is.finite(y))) {
        stop("tip_beta_draws contains non-finite beta values for axis `", ax,
             "` in draw `", draw, "`.", call. = FALSE)
      }
      names(y) <- dd$species
      values[as.character(seq_len(geom$n_tip)), ax] <-
        y[geom$tree$tip.label]
      candidates <- if (evolution_model == "AICc") c("BM", "OU", "EB") else
        evolution_model
      fits <- lapply(candidates, function(model) {
        .hee_fit_response_process(y, geom, model)
      })
      names(fits) <- candidates
      selected <- if (evolution_model == "AICc") {
        candidates[which.min(vapply(fits, `[[`, numeric(1), "aicc"))]
      } else evolution_model
      fit <- fits[[selected]]
      axis_fits[[ax]] <- fit
      for (candidate in fits) {
        fit_rows[[length(fit_rows) + 1L]] <- data.frame(
          response_draw = draw, axis = ax,
          model = candidate$model,
          parameter_per_ma = candidate$parameter,
          boundary_status = candidate$boundary_status,
          sigma2_ml = candidate$sigma2_ml,
          log_likelihood = candidate$log_likelihood,
          aicc = candidate$aicc,
          selected = candidate$model == selected,
          stringsAsFactors = FALSE)
      }
      ancestor_tip_cov <- fit$covariance[node_ids, seq_len(geom$n_tip),
                                         drop = FALSE]
      values[as.character(node_ids), ax] <- as.numeric(
        fit$mu + ancestor_tip_cov %*% fit$inverse_residual
      )
    }

    branch_values <- vector("list", nrow(edge))
    for (ei in seq_len(nrow(edge))) {
      index <- branch_time_index[[ei]]
      if (!length(index)) next
      parent <- edge[ei, 1L]
      child <- edge[ei, 2L]
      depths <- geom$root_age - times[index]
      projected <- matrix(NA_real_, nrow = length(index),
                          ncol = length(axes_to_fit),
                          dimnames = list(NULL, axes_to_fit))
      for (ax in axes_to_fit) {
        projected[, ax] <- .hee_response_branch_bridge(
          parent_value = values[parent, ax],
          child_value = values[child, ax],
          parent_depth = geom$root_age - geom$node_age[parent],
          child_depth = geom$root_age - geom$node_age[child],
          sample_depths = depths,
          fit = axis_fits[[ax]],
          stochastic = branch_uncertainty != "none")
      }
      branch_values[[ei]] <- list(time_index = index, values = projected)
    }

    for (tm in times) {
      time_index <- match(tm, times)
      if (abs(tm - geom$root_age) <= time_tolerance) {
        beta <- values[as.character(geom$root), axes_to_fit, drop = TRUE]
        if (is.null(names(beta))) names(beta) <- axes_to_fit
        rr <- rr + 1L
        row <- data.frame(
          lineage = geom$labels[[geom$root]],
          lineage_type = "root",
          node = geom$root,
          parent_node = NA_integer_,
          branch_start_ma = geom$root_age,
          branch_end_ma = geom$root_age,
          branch_position = 0,
          time_ma = tm,
          response_draw = draw,
          intercept = if (".reconstruct_intercept" %in% names(beta)) {
            unname(beta[[".reconstruct_intercept"]])
          } else 0,
          stringsAsFactors = FALSE
        )
        for (ax in axes) row[[ax]] <- unname(beta[[ax]])
        rows[[rr]] <- row
        next
      }
      active <- which(tm >= geom$node_age[edge[, 2]] - time_tolerance &
                        tm < geom$node_age[edge[, 1]] - time_tolerance)
      for (ei in active) {
        parent <- edge[ei, 1]
        child <- edge[ei, 2]
        start <- geom$node_age[[parent]]
        end <- geom$node_age[[child]]
        L <- max(start - end, 0)
        u <- if (L <= 0) 1 else (start - tm) / L
        branch <- branch_values[[ei]]
        position <- match(time_index, branch$time_index)
        beta <- branch$values[position, ]
        rr <- rr + 1L
        row <- data.frame(
          lineage = geom$labels[[child]],
          lineage_type = if (child <= geom$n_tip) "tip_branch" else
            "internal_branch",
          node = child,
          parent_node = parent,
          branch_start_ma = start,
          branch_end_ma = end,
          branch_position = u,
          time_ma = tm,
          response_draw = draw,
          intercept = if (".reconstruct_intercept" %in% names(beta)) {
            unname(beta[[".reconstruct_intercept"]])
          } else 0,
          stringsAsFactors = FALSE
        )
        for (ax in axes) row[[ax]] <- unname(beta[[ax]])
        rows[[rr]] <- row
      }
    }
  }
  out <- do.call(rbind, rows)
  out$response_weight <- 1
  out$response_source <- paste0(
    "direct_", evolution_model, "_ancestral_environmental_response;",
    "tip_beta_posterior_to_tree_to_palaeo_environment;",
    if (isTRUE(reconstruct_intercept)) "intercept_reconstructed" else
      "intercept_set_to_zero"
  )
  out <- out[order(out$response_draw, -out$time_ma, out$lineage), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "model_comparison") <- do.call(rbind, fit_rows)
  out
}

.hee_summarise_ancestral_suitability <- function(draws) {
  x <- as.data.frame(draws)
  group_cols <- intersect(c("geography_scenario", "climate_scenario",
                            "earth_scenario", "lineage", "cell_id",
                            "time_ma", "region", "lon", "lat", "H_state"),
                          names(x))
  split_key <- interaction(x[, group_cols, drop = FALSE], drop = TRUE,
                           lex.order = TRUE)
  rows <- lapply(split(x, split_key, drop = TRUE), function(z) {
    base <- z[1, group_cols, drop = FALSE]
    w <- if ("response_weight" %in% names(z)) z$response_weight else
      rep(1, nrow(z))
    q <- .hee_hmscee_weighted_quantile(z$suitability, w,
                                       c(0.025, 0.5, 0.975))
    ww <- .hee_hmscee_normalize_weights(w)
    mu <- sum(z$suitability * ww, na.rm = TRUE)
    sdv <- sqrt(sum(ww * (z$suitability - mu)^2, na.rm = TRUE))
    data.frame(
      base,
      eta_mean = sum(z$eta * ww, na.rm = TRUE),
      suitability_mean = mu,
      suitability_sd = sdv,
      suitability_q025 = q[[1]],
      suitability_q50 = q[[2]],
      suitability_q975 = q[[3]],
      n_response_draws = length(unique(z$response_draw)),
      suitability_interpretation =
        "ancestral_environmental_support_not_final_occupancy",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Compute ancestral environmental support on palaeo-Earth cells
#'
#' Convenience wrapper for the scientifically preferred path:
#' `HMSC tip Beta posterior -> dated-tree ancestral Beta ->
#' ancestral Beta x palaeoenvironment -> S_env`. It does not use modern
#' prediction maps and it does not use BioGeoBEARS ranges to define ecological
#' response. The returned suitability is
#' \deqn{S_{\ell,c,t}^{(s)} = g^{-1}[b(X_{c,t})^\top
#' \beta_{\ell,t}^{(s)}],}
#' and should be interpreted as ancestral environmental support or potential
#' suitability, not final historical occupancy probability.
#'
#' @param earth_state Cell-time palaeo-Earth table. It may be raw environmental
#'   data or output from [hee_paleo_earth_state()]. If `recipe` is supplied,
#'   environmental columns are transformed by [hee_apply_recipe()] before
#'   projection.
#' @param tip_beta_draws HMSC Beta posterior draws for modern tips; see
#'   [hee_evolution_ancestral_response_direct()].
#' @param tree Dated `ape::phylo`.
#' @param times Optional ages in Ma. Defaults to unique `earth_state$time_ma`.
#' @param basis_cols Environmental basis columns.
#' @param recipe Optional locked modern recipe. Do not standardize each paleo
#'   time slice independently.
#' @param link `"probit"`, `"logit"`, or `"identity"`.
#' @param ... Additional arguments passed to
#'   [hee_evolution_ancestral_response_direct()].
#'
#' @return A `hee_ancestral_suitability` list with `responses`, `draws`, and
#'   posterior `summary` tables.
#' @export
#'
#' @examples
#' \dontrun{
#' ans <- hee_environmental_filtering_ancestral_suitability(
#'   earth_state = paleo_grid, tip_beta_draws = evo$beta_draws,
#'   tree = dated_tree, basis_cols = c("bio1", "bio12"), link = "probit"
#' )
#' }
hee_environmental_filtering_ancestral_suitability <- function(earth_state,
                                                              tip_beta_draws,
                                                              tree,
                                                              times = NULL,
                                                              basis_cols,
                                                              recipe = NULL,
                                                              link = c("probit", "logit", "identity"),
                                                              ...) {
  link <- match.arg(link)
  e <- as.data.frame(earth_state)
  .require_cols(e, c("cell_id", "time_ma", basis_cols), "earth_state")
  times <- times %||% sort(unique(as.numeric(e$time_ma)), decreasing = TRUE)
  if (!is.null(recipe)) {
    X <- hee_apply_recipe(recipe, e)
    missing_basis <- setdiff(basis_cols, names(X))
    if (length(missing_basis) > 0L) {
      stop("recipe output is missing basis column(s): ",
           paste(missing_basis, collapse = ", "), call. = FALSE)
    }
    for (ax in basis_cols) e[[ax]] <- X[[ax]]
  }
  responses <- hee_evolution_ancestral_response_direct(
    tip_beta_draws = tip_beta_draws,
    tree = tree,
    times = times,
    basis_cols = basis_cols,
    ...
  )
  draws <- hee_lineage_suitability(
    earth_state = e,
    responses = responses,
    basis_cols = basis_cols,
    link = link
  )
  draws$environmental_support <- draws$suitability
  draws$suitability_interpretation <-
    "ancestral_environmental_support_not_final_occupancy"
  summary <- .hee_summarise_ancestral_suitability(draws)
  out <- list(responses = responses, draws = draws, summary = summary)
  class(out) <- c("hee_ancestral_suitability", "list")
  out
}

#' Build trait-mediated ancestral environmental responses
#'
#' This helper keeps three ancestral objects separate. `ancestor_traits`
#' supplies functional traits
#' \deqn{T_{\ell,t}.}
#' `gamma` supplies the HMSC trait-to-response matrix
#' \deqn{\Gamma,}
#' and optional `residual_response` supplies the unexplained lineage response
#' \deqn{u_{\ell,t}.}
#' The returned environmental response is
#' \deqn{\beta_{\ell,t} = \Gamma T_{\ell,t} + u_{\ell,t}.}
#' BioGeoBEARS ancestral ranges are deliberately not used here; geographic
#' history belongs in `hee_bsm_region_history()`.
#'
#' @param ancestor_traits Table with `lineage`, optional `time_ma`, optional
#'   `response_draw`, and trait columns.
#' @param gamma Trait-by-environment response matrix, or long table with
#'   columns `trait`, `axis`, and an effect column such as `gamma` or
#'   `estimate`.
#' @param residual_response Optional table with lineage/time/draw and
#'   environment-axis residual response columns.
#' @param trait_cols Trait columns in `ancestor_traits`. If `NULL`, numeric
#'   columns other than identifiers are used.
#' @param basis_cols Environment-response axes. If `NULL`, axes are inferred
#'   from `gamma`.
#' @param lineage_col,time_col,response_draw_col Column names.
#' @param intercept Baseline establishment/environment linear predictor
#'   intercept for output. This is not interpreted as an ancestral HMSC
#'   prevalence parameter.
#'
#' @return Wide response table suitable for `hee_lineage_suitability()`, with
#'   columns `lineage`, `time_ma`, `response_draw`, `intercept`, and one column
#'   per environmental basis axis.
#' @export
#'
#' @examples
#' traits <- data.frame(lineage = "anc1", time_ma = 10, leaf = 1)
#' G <- matrix(0.5, nrow = 1, dimnames = list("leaf", "bio1"))
#' hee_trait_mediated_ancestral_response(traits, G)
hee_trait_mediated_ancestral_response <- function(ancestor_traits,
                                                  gamma,
                                                  residual_response = NULL,
                                                  trait_cols = NULL,
                                                  basis_cols = NULL,
                                                  lineage_col = "lineage",
                                                  time_col = "time_ma",
                                                  response_draw_col = "response_draw",
                                                  intercept = 0) {
  tr <- as.data.frame(ancestor_traits)
  .require_cols(tr, lineage_col, "ancestor_traits")
  if (!time_col %in% names(tr)) tr[[time_col]] <- NA_real_
  if (!response_draw_col %in% names(tr)) tr[[response_draw_col]] <- "d1"
  if (is.null(trait_cols)) {
    ids <- c(lineage_col, time_col, response_draw_col, "region",
             "R_region", "history_draw", "geography_scenario",
             "climate_scenario", "earth_scenario")
    trait_cols <- setdiff(names(tr)[vapply(tr, is.numeric, logical(1))], ids)
  }
  .require_cols(tr, trait_cols, "ancestor_traits")
  if (length(intersect(c("region", "R_region", "history_draw"), trait_cols))) {
    stop("BioGeoBEARS range columns cannot be used as ancestral niche traits.",
         call. = FALSE)
  }
  if (is.matrix(gamma) || is.data.frame(gamma) &&
      !all(c("trait", "axis") %in% names(gamma))) {
    G <- as.matrix(gamma)
    if (is.null(rownames(G)) || is.null(colnames(G))) {
      stop("matrix gamma must have trait rownames and axis colnames.",
           call. = FALSE)
    }
  } else {
    g <- as.data.frame(gamma)
    .require_cols(g, c("trait", "axis"), "gamma")
    val_col <- .hee_hmscee_first_col(g, c("gamma", "estimate", "mean",
                                          "value"), "gamma")
    G <- stats::xtabs(stats::as.formula(paste(val_col, "~ trait + axis")),
                      data = g)
  }
  if (is.null(basis_cols)) basis_cols <- colnames(G)
  missing_traits <- setdiff(rownames(G), trait_cols)
  if (length(missing_traits)) {
    stop("ancestor_traits is missing gamma trait(s): ",
         paste(missing_traits, collapse = ", "), call. = FALSE)
  }
  G <- G[trait_cols, basis_cols, drop = FALSE]
  key_cols <- c(lineage_col, time_col, response_draw_col)
  .hee_check_unique_keys(tr, key_cols, "ancestor_traits")
  Tmat <- as.matrix(tr[, trait_cols, drop = FALSE])
  beta <- Tmat %*% G
  out <- data.frame(
    lineage = as.character(tr[[lineage_col]]),
    time_ma = suppressWarnings(as.numeric(tr[[time_col]])),
    response_draw = as.character(tr[[response_draw_col]]),
    intercept = intercept,
    stringsAsFactors = FALSE
  )
  out <- cbind(out, as.data.frame(beta, stringsAsFactors = FALSE))
  if (!is.null(residual_response)) {
    u <- as.data.frame(residual_response)
    .require_cols(u, c(lineage_col, time_col, response_draw_col),
                  "residual_response")
    .require_cols(u, basis_cols, "residual_response")
    .hee_check_unique_keys(u, key_cols, "residual_response")
    u2 <- u[, c(key_cols, basis_cols), drop = FALSE]
    names(u2)[seq_along(key_cols)] <- key_cols
    out <- merge(out, u2, by.x = c("lineage", "time_ma", "response_draw"),
                 by.y = key_cols, all.x = TRUE, sort = FALSE,
                 suffixes = c("", ".u"))
    for (ax in basis_cols) {
      ucol <- paste0(ax, ".u")
      if (ucol %in% names(out)) {
        val <- suppressWarnings(as.numeric(out[[ucol]]))
        val[!is.finite(val)] <- 0
        out[[ax]] <- suppressWarnings(as.numeric(out[[ax]])) + val
        out[[ucol]] <- NULL
      }
    }
  }
  out$response_source <- "trait_mediated_beta_equals_Gamma_T_plus_u"
  out
}

#' Compute lineage environmental suitability on palaeo-Earth cells
#'
#' Computes
#' \deqn{S_{\ell,c,k}=g^{-1}[\alpha_{\ell,k}+
#' b(X_{c,k})^\top\theta_{\ell,k}],}
#' from a cell-time environmental design matrix and lineage response draws. The
#' response table may contain posterior or ancestral-response draws; HmscEE does
#' not fit a second phylogenetic model here.
#'
#' @param earth_state Cell-time table with environmental columns.
#' @param responses Lineage response table with `lineage` and coefficients.
#' @param basis_cols Environmental basis columns in `earth_state`.
#' @param intercept_col Intercept column in `responses`.
#' @param response_draw_col Optional posterior/ancestor response draw column.
#' @param earth_scenario_col Legacy scenario column. Prefer
#'   `geography_scenario_col`.
#' @param geography_scenario_col Optional palaeogeographic scenario column `g`
#'   in `earth_state` and, if present, in `responses`. If responses omit it,
#'   response draws are reused across geography scenarios.
#' @param climate_scenario_col Optional palaeoclimate scenario column `c` in
#'   `earth_state` and, if present, in `responses`. If responses omit it,
#'   response draws are reused across climate scenarios.
#' @param response_weight_col Optional response-draw weight column `w[s]`.
#' @param lineage_col,cell_col,time_col Column names.
#' @param link `"logit"`, `"probit"`, or `"identity"`.
#'
#' @return Lineage-cell-time table with `eta`, `suitability`, and draw ids.
#' @export
hee_lineage_suitability <- function(earth_state,
                                    responses,
                                    basis_cols,
                                    intercept_col = "intercept",
                                    response_draw_col = NULL,
                                    earth_scenario_col = NULL,
                                    geography_scenario_col = NULL,
                                    climate_scenario_col = NULL,
                                    response_weight_col = NULL,
                                    lineage_col = "lineage",
                                    cell_col = "cell_id",
                                    time_col = "time_ma",
                                    link = c("logit", "probit", "identity")) {
  link <- match.arg(link)
  e <- as.data.frame(earth_state)
  r <- as.data.frame(responses)
  .require_cols(e, c(cell_col, time_col, basis_cols), "earth_state")
  .require_cols(r, c(lineage_col, intercept_col, basis_cols), "responses")
  geography_scenario_col <- .hee_hmscee_pick_scenario_col(
    e, "geography_scenario", legacy = "earth_scenario",
    explicit = geography_scenario_col
  )
  if (is.null(geography_scenario_col) && !is.null(earth_scenario_col) &&
      earth_scenario_col %in% names(e)) {
    geography_scenario_col <- earth_scenario_col
  }
  if (is.null(geography_scenario_col) ||
      !geography_scenario_col %in% names(e)) {
    e$geography_scenario <- "g1"
    geography_scenario_col <- "geography_scenario"
  }
  climate_scenario_col <- .hee_hmscee_pick_scenario_col(
    e, "climate_scenario", explicit = climate_scenario_col
  )
  if (is.null(climate_scenario_col) ||
      !climate_scenario_col %in% names(e)) {
    e$climate_scenario <- "c1"
    climate_scenario_col <- "climate_scenario"
  }
  if (!"climate_weight_given_geography" %in% names(e)) {
    e$climate_weight_given_geography <- 1
  }
  if (is.null(response_draw_col) && "response_draw" %in% names(r)) {
    response_draw_col <- "response_draw"
  }
  if (is.null(response_weight_col) && "response_weight" %in% names(r)) {
    response_weight_col <- "response_weight"
  }
  if (is.null(response_draw_col) || !response_draw_col %in% names(r)) {
    r$response_draw <- "d1"
    response_draw_col <- "response_draw"
  }
  if (is.null(response_weight_col) || !response_weight_col %in% names(r)) {
    r$response_weight <- 1
    response_weight_col <- "response_weight"
  }
  response_geography_col <- .hee_hmscee_pick_scenario_col(
    r, "geography_scenario", legacy = "earth_scenario"
  )
  response_climate_col <- .hee_hmscee_pick_scenario_col(
    r, "climate_scenario"
  )
  response_has_geography <- !is.null(response_geography_col)
  response_has_climate <- !is.null(response_climate_col)
  response_has_time <- time_col %in% names(r)
  response_keys <- c(
    if (response_has_geography) response_geography_col else character(),
    if (response_has_climate) response_climate_col else character(),
    lineage_col, response_draw_col
  )
  if (response_has_time) response_keys <- c(response_keys, time_col)
  .hee_check_unique_keys(r, response_keys, "responses")
  bad <- !stats::complete.cases(e[, basis_cols, drop = FALSE])
  if (any(bad)) {
    warning(sum(bad), " earth_state row(s) have incomplete basis values; ",
            "suitability is set to NA before hard-state masking.",
            call. = FALSE)
  }
  rows <- vector("list", nrow(r))
  for (i in seq_len(nrow(r))) {
    e_i <- if (response_has_geography) {
      e[as.character(e[[geography_scenario_col]]) ==
          as.character(r[[response_geography_col]][i]), , drop = FALSE]
    } else {
      e
    }
    if (response_has_climate) {
      e_i <- e_i[as.character(e_i[[climate_scenario_col]]) ==
                   as.character(r[[response_climate_col]][i]), , drop = FALSE]
    }
    if (response_has_time) {
      e_i <- e_i[as.numeric(e_i[[time_col]]) ==
                   as.numeric(r[[time_col]][i]), , drop = FALSE]
    }
    if (!nrow(e_i)) next
    bad_i <- !stats::complete.cases(e_i[, basis_cols, drop = FALSE])
    coefs <- suppressWarnings(as.numeric(r[i, basis_cols, drop = TRUE]))
    eta <- suppressWarnings(as.numeric(r[[intercept_col]][i])) +
      as.matrix(e_i[, basis_cols, drop = FALSE]) %*% coefs
    eta <- as.numeric(eta)
    eta[bad_i] <- NA_real_
    s <- switch(link,
                logit = stats::plogis(eta),
                probit = stats::pnorm(eta),
                identity = eta)
    rows[[i]] <- data.frame(
      geography_scenario = as.character(e_i[[geography_scenario_col]]),
      climate_scenario = as.character(e_i[[climate_scenario_col]]),
      earth_scenario = as.character(e_i[[geography_scenario_col]]),
      lineage = as.character(r[[lineage_col]][i]),
      response_draw = as.character(r[[response_draw_col]][i]),
      response_weight = suppressWarnings(as.numeric(r[[response_weight_col]][i])),
      climate_weight_given_geography = suppressWarnings(as.numeric(
        e_i$climate_weight_given_geography
      )),
      cell_id = as.character(e_i[[cell_col]]),
      time_ma = as.numeric(e_i[[time_col]]),
      eta = eta,
      suitability = .hee_clip01(s),
      stringsAsFactors = FALSE
    )
    keep <- intersect(c("region", "lon", "lat", "H_state"), names(e_i))
    if (length(keep)) rows[[i]] <- cbind(rows[[i]], e_i[, keep, drop = FALSE])
  }
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  out$response_weight[!is.finite(out$response_weight) |
                        out$response_weight < 0] <- 0
  out
}

#' Build the within-region movement kernel
#'
#' Converts least-cost distances and optional dispersal traits into the only
#' movement term used by the nested HmscEE core:
#' \deqn{K_{\ell,c'\rightarrow c,k}=1-\exp[-\lambda^{mov}_{\ell,c'\rightarrow c,k}\Delta t_k].}
#' Cross-region movement is set to zero by default because BioGeoBEARS/BSM
#' already handles range expansion, contraction, and founder events among
#' regions.
#'
#' @param paths Edge table with source/target cells or regions, time, and a
#'   least-cost distance column.
#' @param traits Optional lineage traits with `dispersal_trait_col`.
#' @param delta_t Optional scalar, vector, or edge column with interval length
#'   in Myr.
#' @param movement_intercept,movement_trait_coef,distance_decay Movement-rate
#'   parameters. Distances must use the same unit as `distance_decay`.
#' @param allow_cross_region Logical. If `FALSE`, cross-region edges get
#'   `K_movement = 0`.
#' @param lineage_col,from_col,to_col,from_region_col,to_region_col,time_col
#'   Column names.
#' @param distance_col Least-cost distance column.
#' @param dispersal_trait_col Optional trait column.
#'
#' @return Edge table with `lambda_movement`, `K_movement`, `delta_t`, and
#'   `movement_scale`.
#' @export
hee_within_region_movement_kernel <- function(paths,
                                              traits = NULL,
                                              delta_t = 1,
                                              movement_intercept = -2,
                                              movement_trait_coef = 0,
                                              distance_decay = 0.01,
                                              allow_cross_region = FALSE,
                                              lineage_col = "lineage",
                                              from_col = "from_cell_id",
                                              to_col = "to_cell_id",
                                              from_region_col = "from_region",
                                              to_region_col = "to_region",
                                              time_col = "time_ma",
                                              distance_col = "least_cost_distance_km",
                                              dispersal_trait_col = NULL) {
  x <- as.data.frame(paths)
  .require_cols(x, c(lineage_col, from_col, to_col, time_col, distance_col),
                "paths")
  if (!isTRUE(allow_cross_region)) {
    .require_cols(x, c(from_region_col, to_region_col), "paths")
  }
  .hee_validate_time_values(x[[time_col]], paste0("paths$", time_col))
  path_keys <- c(lineage_col, from_col, to_col, time_col)
  if ("geography_scenario" %in% names(x)) {
    path_keys <- c("geography_scenario", path_keys)
  } else if ("earth_scenario" %in% names(x)) {
    x$geography_scenario <- as.character(x$earth_scenario)
    path_keys <- c("geography_scenario", path_keys)
  }
  if ("climate_scenario" %in% names(x)) {
    path_keys <- c("climate_scenario", path_keys)
  }
  .hee_check_unique_keys(x, path_keys, "paths")
  d <- suppressWarnings(as.numeric(x[[distance_col]]))
  if (any(!is.finite(d) | d < 0, na.rm = TRUE)) {
    stop("least-cost distances must be finite and non-negative.", call. = FALSE)
  }
  dt <- if (length(delta_t) == 1L && is.character(delta_t) &&
            delta_t %in% names(x)) {
    suppressWarnings(as.numeric(x[[delta_t]]))
  } else {
    rep_len(suppressWarnings(as.numeric(delta_t)), nrow(x))
  }
  if (any(!is.finite(dt) | dt < 0, na.rm = TRUE)) {
    stop("delta_t must be finite and non-negative.", call. = FALSE)
  }
  trait <- rep(0, nrow(x))
  if (!is.null(traits)) {
    if (is.null(dispersal_trait_col)) {
      dispersal_trait_col <- intersect(c("dispersal_trait",
                                         "dispersal_distance_km",
                                         "dispersal_ability"),
                                       names(traits))[1L]
    }
    if (is.na(dispersal_trait_col) || !dispersal_trait_col %in% names(traits)) {
      stop("traits must contain a dispersal trait column.", call. = FALSE)
    }
    tr <- as.data.frame(traits)
    .require_cols(tr, c(lineage_col, dispersal_trait_col), "traits")
    .hee_check_unique_keys(tr, lineage_col, "traits")
    trait <- suppressWarnings(as.numeric(tr[[dispersal_trait_col]][
      match(as.character(x[[lineage_col]]), as.character(tr[[lineage_col]]))
    ]))
    trait[!is.finite(trait)] <- 0
  }
  eta <- movement_intercept + movement_trait_coef * trait -
    distance_decay * pmax(d, 0)
  lambda <- exp(eta)
  K <- 1 - exp(-lambda * dt)
  if (!isTRUE(allow_cross_region)) {
    cross <- as.character(x[[from_region_col]]) != as.character(x[[to_region_col]])
    K[cross | is.na(cross)] <- 0
  }
  x$delta_t <- dt
  x$lambda_movement <- lambda
  x$K_movement <- .hee_clip01(K)
  x$movement_scale <- ifelse(isTRUE(allow_cross_region),
                             "user_allowed_cross_region",
                             "within_BioGeoBEARS_region_only")
  x
}

#' Compute within-region arrival pressure
#'
#' For target cell `c`, arrival pressure is
#' \deqn{A^{arr}_{\ell,c,k}=1-\prod_{c'\in r(c,k)}
#' [1-q_{\ell,c',k}K_{\ell,c'\rightarrow c,k}],}
#' the probability that at least one occupied source cell in the same
#' BioGeoBEARS region sends propagules to the target.
#'
#' @param previous_state Lineage-cell previous occupancy table.
#' @param movement_kernel Output from `hee_within_region_movement_kernel()`.
#' @param target_cells Optional target lineage-cell-time table.
#' @param lineage_col,cell_col,time_col,probability_col Column names.
#'
#' @return Target table with `arrival_pressure`.
#' @export
hee_arrival_pressure <- function(previous_state,
                                 movement_kernel,
                                 target_cells = NULL,
                                 lineage_col = "lineage",
                                 cell_col = "cell_id",
                                 time_col = "time_ma",
                                 probability_col = "q") {
  prev <- as.data.frame(previous_state)
  ker <- as.data.frame(movement_kernel)
  .require_cols(prev, c(lineage_col, cell_col, probability_col),
                "previous_state")
  .require_cols(ker, c(lineage_col, "from_cell_id", "to_cell_id", time_col,
                       "K_movement"), "movement_kernel")
  .hee_check_unique_keys(prev, c(lineage_col, cell_col), "previous_state")
  .hee_check_unique_keys(ker, c(lineage_col, "from_cell_id", "to_cell_id",
                                time_col), "movement_kernel")
  prev$.source_q <- .hee_hmscee_prob(prev[[probability_col]])
  src <- prev[, c(lineage_col, cell_col, ".source_q"), drop = FALSE]
  names(src)[names(src) == cell_col] <- "from_cell_id"
  z <- merge(ker, src, by = c(lineage_col, "from_cell_id"), all.x = TRUE,
             sort = FALSE)
  z$.source_q[is.na(z$.source_q)] <- 0
  z$contribution <- .hee_clip01(z$.source_q * .hee_hmscee_prob(z$K_movement))
  out <- stats::aggregate(
    z$contribution,
    z[, c(lineage_col, "to_cell_id", time_col), drop = FALSE],
    function(v) 1 - prod(1 - .hee_clip01(v), na.rm = TRUE)
  )
  names(out) <- c(lineage_col, cell_col, time_col, "arrival_pressure")
  out$arrival_pressure <- .hee_clip01(out$arrival_pressure)
  if (!is.null(target_cells)) {
    tgt <- as.data.frame(target_cells)
    .require_cols(tgt, c(lineage_col, cell_col, time_col), "target_cells")
    keys <- c(lineage_col, cell_col, time_col)
    .hee_check_unique_keys(tgt, keys, "target_cells")
    .hee_check_unique_keys(out, keys, "arrival_pressure")
    n_before <- nrow(tgt)
    out <- merge(tgt, out, by = keys, all.x = TRUE, sort = FALSE)
    if (nrow(out) != n_before) {
      stop("Joining arrival_pressure changed the number of target rows.",
           call. = FALSE)
    }
    out$arrival_pressure[is.na(out$arrival_pressure)] <- 0
  }
  out
}

#' Convert environmental response into establishment probability
#'
#' Establishment is conditional on arrival:
#' \deqn{Est_{\ell,c,k}=logit^{-1}(a_{\ell}+b_S\eta_{\ell,c,k}).}
#'
#' @param x Table or numeric vector containing the environmental linear
#'   predictor.
#' @param eta_col Column name when `x` is a data.frame.
#' @param intercept,slope Logistic establishment parameters.
#' @return `x` with `establishment_probability`, or a numeric vector.
#' @export
hee_establishment_probability <- function(x,
                                          eta_col = "eta",
                                          intercept = 0,
                                          slope = 1) {
  validate_parameter <- function(value, n, label) {
    value <- suppressWarnings(as.numeric(value))
    if (!length(value) || any(!is.finite(value)) ||
        !(length(value) %in% c(1L, n))) {
      stop(label, " must be one finite value or one finite value per eta.",
           call. = FALSE)
    }
    value
  }
  if (is.data.frame(x)) {
    .require_cols(x, eta_col, "x")
    eta <- suppressWarnings(as.numeric(x[[eta_col]]))
    intercept <- validate_parameter(intercept, length(eta), "intercept")
    slope <- validate_parameter(slope, length(eta), "slope")
    x$establishment_probability <- .hee_clip01(stats::plogis(intercept + slope * eta))
    x$establishment_probability[!is.finite(eta)] <- 0
    return(x)
  }
  eta <- suppressWarnings(as.numeric(x))
  intercept <- validate_parameter(intercept, length(eta), "intercept")
  slope <- validate_parameter(slope, length(eta), "slope")
  out <- .hee_clip01(stats::plogis(intercept + slope * eta))
  out[!is.finite(eta)] <- 0
  out
}

#' Combine arrival and establishment into colonisation probability
#'
#' \deqn{\gamma_{\ell,c,k}=A^{arr}_{\ell,c,k}Est_{\ell,c,k}.}
#'
#' @param x Table with arrival and establishment columns.
#' @param arrival_col,establishment_col Column names.
#' @return `x` with `colonisation_probability`.
#' @export
hee_colonisation_from_arrival <- function(x,
                                          arrival_col = "arrival_pressure",
                                          establishment_col = "establishment_probability") {
  d <- as.data.frame(x)
  .require_cols(d, c(arrival_col, establishment_col), "x")
  d$colonisation_probability <- .hee_clip01(
    .hee_hmscee_prob(d[[arrival_col]]) *
      .hee_hmscee_prob(d[[establishment_col]])
  )
  d
}

#' Convert environmental response into persistence probability
#'
#' Core persistence is
#' \deqn{\phi_{\ell,c,k}=logit^{-1}(u_{\ell}+b_P\eta_{\ell,c,k}).}
#' Optional habitat, area, rescue, and disturbance terms are accepted only when
#' supplied explicitly; omitted optional terms are neutral. Local extinction is
#' `1 - phi`.
#'
#' @param x Table or numeric vector containing `eta`.
#' @param eta_col Column name when `x` is a data.frame.
#' @param intercept,slope Logistic persistence parameters.
#' @param habitat_col,area_col,rescue_col,disturbance_col Optional columns.
#' @param habitat_coef,area_coef,rescue_coef,disturbance_coef Optional
#'   coefficients. Disturbance is subtracted.
#' @return `x` with `persistence_probability` and
#'   `local_extinction_probability`, or a numeric vector of persistence values.
#' @export
hee_persistence_probability <- function(x,
                                        eta_col = "eta",
                                        intercept = 1,
                                        slope = 1,
                                        habitat_col = NULL,
                                        area_col = NULL,
                                        rescue_col = NULL,
                                        disturbance_col = NULL,
                                        habitat_coef = 0,
                                        area_coef = 0,
                                        rescue_coef = 0,
                                        disturbance_coef = 0) {
  calc <- function(d, eta) {
    validate_parameter <- function(value, n, label) {
      value <- suppressWarnings(as.numeric(value))
      if (!length(value) || any(!is.finite(value)) ||
          !(length(value) %in% c(1L, n))) {
        stop(label, " must be one finite value or one finite value per eta.",
             call. = FALSE)
      }
      value
    }
    intercept <- validate_parameter(intercept, length(eta), "intercept")
    slope <- validate_parameter(slope, length(eta), "slope")
    lp <- intercept + slope * eta
    add <- function(col, coef, sign = 1) {
      if (!is.null(col)) {
        .require_cols(d, col, "x")
        v <- suppressWarnings(as.numeric(d[[col]]))
        v[!is.finite(v)] <- 0
        lp <<- lp + sign * coef * v
      }
    }
    add(habitat_col, habitat_coef)
    add(area_col, area_coef)
    add(rescue_col, rescue_coef)
    add(disturbance_col, disturbance_coef, sign = -1)
    out <- .hee_clip01(stats::plogis(lp))
    out[!is.finite(eta)] <- 0
    out
  }
  if (is.data.frame(x)) {
    .require_cols(x, eta_col, "x")
    eta <- suppressWarnings(as.numeric(x[[eta_col]]))
    x$persistence_probability <- calc(x, eta)
    x$local_extinction_probability <- .hee_clip01(1 - x$persistence_probability)
    return(x)
  }
  eta <- suppressWarnings(as.numeric(x))
  calc(data.frame(.dummy = seq_along(eta)), eta)
}

#' Run global-grid dynamic occupancy without BioGeoBEARS regions
#'
#' This is the no-BioGeoBEARS dynamic core used by Case04 when discrete
#' regions are intentionally removed. It updates lineage-cell occupancy
#' directly on a global cell graph:
#' \deqn{\psi_{\ell,c,t+\delta t}=(1-\psi_{\ell,c,t})P_{01}+
#' \psi_{\ell,c,t}P_{11}.}
#' Arrival is computed from occupied source cells and the supplied movement
#' kernel, establishment is conditional on arrival, and persistence controls
#' local loss. `H_state = 0` is a hard palaeogeographic constraint and forces
#' occupancy to zero. This function does not use `R_region`, `A_hist`,
#' BioGeoBEARS/BSM histories, or M1-M5 multipliers.
#'
#' @param suitability Lineage-cell-time environmental support table. Required
#'   columns are `lineage`, `cell_id`, `time_ma`, `eta`, and `suitability`.
#'   `H_state` is optional and defaults to one.
#' @param movement_kernel Optional global movement edge table with `lineage`,
#'   `from_cell_id`, `to_cell_id`, `time_ma`, and `K_movement`. Unlike
#'   [hee_nested_region_cell_occupancy()], no region boundary is used here.
#' @param root_state Optional initial occupancy table with `lineage`,
#'   `cell_id`, and `q_root`. If omitted, the oldest time slice is initialized
#'   from environmental support and `initial_rho`.
#' @param initial_rho Initial occupancy mass used when `root_state` is absent.
#' @param establishment_intercept,establishment_slope Parameters for
#'   [hee_establishment_probability()].
#' @param persistence_intercept,persistence_slope Parameters for
#'   [hee_persistence_probability()].
#' @param internal_dt Maximum internal integration step in Myr. Smaller values
#'   repeatedly recompute arrival pressure within each palaeoenvironmental
#'   interval.
#' @param summarise_draws If `TRUE`, summarize response draws using
#'   `response_weight`.
#'
#' @return A list with `draws`, `summary`, and `formula_catalog`.
#' @export
hee_global_grid_dynamic_occupancy <- function(suitability,
                                              movement_kernel = NULL,
                                              root_state = NULL,
                                              initial_rho = 0.01,
                                              establishment_intercept = 0,
                                              establishment_slope = 1,
                                              persistence_intercept = 1,
                                              persistence_slope = 1,
                                              internal_dt = 0.5,
                                              summarise_draws = TRUE) {
  s <- as.data.frame(suitability)
  .require_cols(s, c("lineage", "cell_id", "time_ma", "eta", "suitability"),
                "suitability")
  .hee_validate_time_values(s$time_ma, "suitability$time_ma")
  .hee_check_unique_keys(
    s,
    intersect(c("geography_scenario", "climate_scenario", "response_draw",
                "lineage", "cell_id", "time_ma"), names(s)),
    "suitability"
  )
  if (!"geography_scenario" %in% names(s)) s$geography_scenario <- "g1"
  if (!"climate_scenario" %in% names(s)) s$climate_scenario <- "c1"
  if (!"response_draw" %in% names(s)) s$response_draw <- "d1"
  if (!"response_weight" %in% names(s)) s$response_weight <- 1
  if (!"H_state" %in% names(s)) s$H_state <- 1
  s$lineage <- as.character(s$lineage)
  s$cell_id <- as.character(s$cell_id)
  s$time_ma <- as.numeric(s$time_ma)
  s$eta <- suppressWarnings(as.numeric(s$eta))
  s$suitability <- .hee_hmscee_prob(s$suitability)
  s$H_state <- as.integer(.hee_hmscee_prob(s$H_state) > 0)
  s$response_weight <- suppressWarnings(as.numeric(s$response_weight))
  s$response_weight[!is.finite(s$response_weight) | s$response_weight < 0] <- 0
  id_cols <- c("geography_scenario", "climate_scenario", "response_draw",
               "lineage")
  times <- sort(unique(s$time_ma), decreasing = TRUE)
  dt_max <- suppressWarnings(as.numeric(internal_dt))[1L]
  if (!is.finite(dt_max) || dt_max <= 0) {
    stop("internal_dt must be a finite positive number of Myr.", call. = FALSE)
  }

  mk <- if (!is.null(movement_kernel)) as.data.frame(movement_kernel) else NULL
  if (!is.null(mk)) {
    .require_cols(mk, c("lineage", "from_cell_id", "to_cell_id", "time_ma",
                        "K_movement"), "movement_kernel")
    .hee_check_unique_keys(
      mk,
      intersect(c("geography_scenario", "climate_scenario", "response_draw",
                  "lineage", "from_cell_id", "to_cell_id", "time_ma"),
                names(mk)),
      "movement_kernel"
    )
    if (!"geography_scenario" %in% names(mk)) mk$geography_scenario <- "g1"
    if (!"climate_scenario" %in% names(mk)) mk$climate_scenario <- "c1"
    if (!"response_draw" %in% names(mk)) mk$response_draw <- NA_character_
    mk$lineage <- as.character(mk$lineage)
    mk$from_cell_id <- as.character(mk$from_cell_id)
    mk$to_cell_id <- as.character(mk$to_cell_id)
    mk$time_ma <- as.numeric(mk$time_ma)
    mk$K_movement <- .hee_hmscee_prob(mk$K_movement)
  }

  root <- if (!is.null(root_state)) as.data.frame(root_state) else NULL
  if (!is.null(root)) {
    .require_cols(root, c("lineage", "cell_id", "q_root"), "root_state")
    if (!"geography_scenario" %in% names(root)) root$geography_scenario <- "g1"
    if (!"climate_scenario" %in% names(root)) root$climate_scenario <- "c1"
    if (!"response_draw" %in% names(root)) root$response_draw <- NA_character_
    root$lineage <- as.character(root$lineage)
    root$cell_id <- as.character(root$cell_id)
    root$q_root <- .hee_hmscee_prob(root$q_root)
    .hee_check_unique_keys(
      root,
      intersect(c("geography_scenario", "climate_scenario", "response_draw",
                  "lineage", "cell_id"), names(root)),
      "root_state"
    )
  }

  s$q <- 0
  s$arrival_pressure <- 0
  s$establishment_probability <- 0
  s$colonisation_probability <- 0
  s$persistence_probability <- 0
  s$local_extinction_probability <- 1
  s$colonisation_hazard <- 0
  s$local_loss_hazard <- 0
  s$delta_t <- NA_real_
  s$n_internal_steps <- 0L
  s$transition_case <- "not_updated"
  draw_key <- do.call(.hee_hmscee_key, s[id_cols])

  prob_to_rate <- function(p, dt) {
    p <- .hee_clip01(p)
    p[p >= 1] <- 1 - 1e-12
    out <- -log1p(-p) / dt
    out[!is.finite(out) | is.na(out) | out < 0] <- 0
    out
  }
  persistence_to_loss_rate <- function(phi, dt) {
    phi <- .hee_clip01(phi)
    phi[phi <= 0] <- 1e-12
    out <- -log(phi) / dt
    out[!is.finite(out) | is.na(out) | out < 0] <- 0
    out
  }
  ctmc_update <- function(q0, lambda_c, lambda_l, h) {
    ss <- lambda_c + lambda_l
    p01 <- ifelse(ss > 0, lambda_c / ss * (1 - exp(-ss * h)), 0)
    p11 <- ifelse(ss > 0, lambda_c / ss + lambda_l / ss * exp(-ss * h), 1)
    .hee_clip01((1 - q0) * p01 + q0 * p11)
  }

  for (grp in unique(draw_key)) {
    idx_grp <- which(draw_key == grp)
    ztimes <- times[times %in% s$time_ma[idx_grp]]
    prev_state <- NULL
    prev_time <- NA_real_
    for (tt in ztimes) {
      idx <- idx_grp[s$time_ma[idx_grp] == tt]
      cur <- s[idx, , drop = FALSE]
      absent <- cur$H_state <= 0
      cur$establishment_probability <- hee_establishment_probability(
        cur, eta_col = "eta", intercept = establishment_intercept,
        slope = establishment_slope
      )$establishment_probability
      cur$persistence_probability <- hee_persistence_probability(
        cur, eta_col = "eta", intercept = persistence_intercept,
        slope = persistence_slope
      )$persistence_probability
      cur$local_extinction_probability <- .hee_clip01(1 - cur$persistence_probability)

      if (is.null(prev_state)) {
        cur$transition_case <- ifelse(absent, "hard_zero", "root_initialisation")
        if (!is.null(root)) {
          rcur <- root[
            root$geography_scenario %in% unique(cur$geography_scenario) &
              root$climate_scenario %in% unique(cur$climate_scenario) &
              (is.na(root$response_draw) |
                 root$response_draw %in% unique(cur$response_draw)) &
              root$lineage %in% unique(cur$lineage),
            , drop = FALSE
          ]
          if (nrow(rcur)) {
            rkey <- .hee_hmscee_key(rcur$lineage, rcur$cell_id)
            ckey <- .hee_hmscee_key(cur$lineage, cur$cell_id)
            cur$q <- .hee_hmscee_prob(rcur$q_root[match(ckey, rkey)])
            cur$q[is.na(cur$q)] <- 0
          }
        } else if (any(!absent)) {
          w <- .hee_hmscee_prob(cur$suitability) * as.numeric(cur$H_state > 0)
          total <- sum(w, na.rm = TRUE)
          if (is.finite(total) && total > 0) {
            cur$q <- .hee_clip01(initial_rho * w / total)
          }
        }
        cur$q[absent] <- 0
      } else {
        dt <- abs(prev_time - tt)
        if (!is.finite(dt) || dt < 0) {
          stop("time_ma intervals must be finite.", call. = FALSE)
        }
        n_steps <- max(1L, ceiling(dt / dt_max))
        h <- if (n_steps > 0) dt / n_steps else 0
        q_lookup <- stats::setNames(prev_state$q, prev_state$cell_id)
        qcur <- .hee_hmscee_prob(q_lookup[cur$cell_id])
        qcur[is.na(qcur)] <- 0
        cur$delta_t <- dt
        cur$n_internal_steps <- n_steps
        cur$transition_case <- ifelse(absent, "hard_zero",
                                      "global_grid_ctmc_update")

        for (step in seq_len(n_steps)) {
          step_prev <- cur[, c("lineage", "cell_id"), drop = FALSE]
          step_prev$q <- qcur
          if (is.null(mk)) {
            arr <- cur[, c("lineage", "cell_id", "time_ma"), drop = FALSE]
            arr$arrival_pressure <- qcur
          } else {
            kcur <- mk[mk$time_ma == tt &
                         mk$lineage %in% unique(cur$lineage) &
                         mk$geography_scenario %in%
                         unique(cur$geography_scenario) &
                         mk$climate_scenario %in%
                         unique(cur$climate_scenario), , drop = FALSE]
            if (any(!is.na(kcur$response_draw))) {
              kcur <- kcur[is.na(kcur$response_draw) |
                             kcur$response_draw %in%
                             unique(cur$response_draw), , drop = FALSE]
            }
            arr <- hee_arrival_pressure(
              step_prev,
              kcur,
              cur[, c("lineage", "cell_id", "time_ma"), drop = FALSE]
            )
          }
          akey <- .hee_hmscee_key(arr$lineage, arr$cell_id, arr$time_ma)
          ckey <- .hee_hmscee_key(cur$lineage, cur$cell_id, cur$time_ma)
          arrival <- .hee_hmscee_prob(arr$arrival_pressure[match(ckey, akey)])
          arrival[is.na(arrival)] <- 0
          gamma_prob <- .hee_clip01(arrival * cur$establishment_probability)
          lambda_c <- prob_to_rate(gamma_prob, max(dt, .Machine$double.eps))
          lambda_l <- persistence_to_loss_rate(cur$persistence_probability,
                                               max(dt, .Machine$double.eps))
          qcur <- ctmc_update(qcur, lambda_c, lambda_l, h)
          qcur[absent] <- 0
          if (step == n_steps) {
            cur$arrival_pressure <- arrival
            cur$colonisation_probability <- gamma_prob
            cur$colonisation_hazard <- lambda_c
            cur$local_loss_hazard <- lambda_l
          }
        }
        cur$q <- .hee_clip01(qcur)
      }
      cur$q[absent] <- 0
      s[idx, names(cur)] <- cur
      prev_state <- cur[, c("lineage", "cell_id", "q"), drop = FALSE]
      prev_time <- tt
    }
  }
  s$probability_draw <- .hee_clip01(s$q)
  if (!isTRUE(summarise_draws)) {
    return(list(draws = s, summary = NULL,
                formula_catalog = hee_global_grid_dynamic_formula_catalog()))
  }
  s$.summary_key <- .hee_hmscee_key(s$lineage, s$cell_id, s$time_ma)
  parts <- split(seq_len(nrow(s)), s$.summary_key)
  summary <- do.call(rbind, lapply(parts, function(ii) {
    z <- s[ii, , drop = FALSE]
    w <- .hee_hmscee_normalize_weights(z$response_weight)
    p <- .hee_hmscee_prob(z$probability_draw)
    qs <- .hee_hmscee_weighted_quantile(p, w, c(0.025, 0.5, 0.975))
    mu <- sum(p * w)
    data.frame(
      lineage = z$lineage[[1]],
      cell_id = z$cell_id[[1]],
      time_ma = z$time_ma[[1]],
      probability_mean = .hee_clip01(mu),
      probability_q025 = .hee_clip01(qs[[1]]),
      probability_q50 = .hee_clip01(qs[[2]]),
      probability_q975 = .hee_clip01(qs[[3]]),
      probability_sd = sqrt(sum(w * (p - mu)^2)),
      n_response_draws = length(unique(z$response_draw)),
      stringsAsFactors = FALSE
    )
  }))
  coords <- unique(s[, intersect(c("cell_id", "time_ma", "lon", "lat",
                                  "H_state"), names(s)), drop = FALSE])
  if (all(c("cell_id", "time_ma") %in% names(coords))) {
    coords <- stats::aggregate(
      coords[, intersect(c("lon", "lat", "H_state"), names(coords)),
             drop = FALSE],
      coords[, c("cell_id", "time_ma"), drop = FALSE],
      function(v) if (is.numeric(v) || is.integer(v)) mean(v, na.rm = TRUE) else
        paste(unique(as.character(v)), collapse = ";")
    )
    summary <- merge(summary, coords, by = c("cell_id", "time_ma"),
                     all.x = TRUE, sort = FALSE)
  }
  list(draws = s, summary = summary,
       formula_catalog = hee_global_grid_dynamic_formula_catalog())
}

#' Formula catalog for no-BioGeoBEARS global grid dynamics
#'
#' @return Data frame describing the formulas used by
#'   [hee_global_grid_dynamic_occupancy()].
#' @export
hee_global_grid_dynamic_formula_catalog <- function() {
  data.frame(
    formula_id = c("global_grid_ctmc_update", "arrival_pressure",
                   "colonisation_hazard", "local_loss_hazard"),
    formula = c(
      "psi_next = (1 - psi) * P01 + psi * P11",
      "A_arr = 1 - prod(1 - q_source * K_movement)",
      "lambda_C = -log(1 - A_arr * Establishment) / delta_t",
      "lambda_L = -log(Persistence) / delta_t"
    ),
    interpretation = c(
      "No-BioGeoBEARS global cell-level occupancy update.",
      "Propagule arrival from occupied source cells on the global cell graph.",
      "Successful colonisation rate; zero arrival gives zero colonisation.",
      "Local loss rate derived once from persistence."
    ),
    scientific_boundary = c(
      "Endpoint/fossil conditioning must be added outside this transition for reconstruction.",
      "Movement is not a BioGeoBEARS regional history and must use explicit edges.",
      "Scenario or prior parameter unless calibrated by independent dynamic data.",
      "Scenario or prior parameter unless calibrated by independent dynamic data."
    ),
    stringsAsFactors = FALSE
  )
}

#' Initialize cells when BioGeoBEARS records entry into a new region
#'
#' When BSM histories switch from `R=0` to `R=1`, HmscEE does not re-estimate
#' cross-region dispersal. It allocates the already-accepted regional entry
#' event to cells inside the target region using
#' \deqn{w^{entry}_{\ell,c,k}=H_{c,k}S_{\ell,c,k}B_{c|source\rightarrow r,k}.}
#'
#' @param target_cells Lineage-cell table for one or more target regions.
#' @param rho Initial within-region occupancy mass. May be scalar or column.
#' @param suitability_col,habitat_col,entry_bias_col Column names.
#' @param group_cols Grouping columns, usually lineage, region, time, and draw.
#' @return `target_cells` with `entry_weight` and `q_entry`.
#' @export
hee_entry_initialization <- function(target_cells,
                                     rho = 0.25,
                                     suitability_col = "suitability",
                                     habitat_col = "H_state",
                                     entry_bias_col = NULL,
                                     group_cols = c("history_draw", "response_draw",
                                                    "lineage", "region", "time_ma")) {
  x <- as.data.frame(target_cells)
  group_cols <- intersect(group_cols, names(x))
  .require_cols(x, c(group_cols, suitability_col), "target_cells")
  if (!habitat_col %in% names(x)) x[[habitat_col]] <- 1
  s <- .hee_hmscee_prob(x[[suitability_col]])
  h <- as.numeric(.hee_hmscee_prob(x[[habitat_col]]) > 0)
  b <- if (!is.null(entry_bias_col)) {
    .require_cols(x, entry_bias_col, "target_cells")
    .hee_hmscee_prob(x[[entry_bias_col]])
  } else {
    rep(1, nrow(x))
  }
  x$entry_weight <- h * s * b
  rho_vec <- if (length(rho) == 1L && is.character(rho) && rho %in% names(x)) {
    .hee_hmscee_prob(x[[rho]])
  } else {
    rep_len(.hee_hmscee_prob(rho), nrow(x))
  }
  key <- if (length(group_cols)) do.call(.hee_hmscee_key, x[group_cols]) else
    rep("all", nrow(x))
  x$q_entry <- 0
  for (kk in unique(key)) {
    idx <- which(key == kk)
    total <- sum(x$entry_weight[idx], na.rm = TRUE)
    if (is.finite(total) && total > 0) {
      x$q_entry[idx] <- .hee_clip01(rho_vec[idx] * x$entry_weight[idx] / total)
    }
  }
  x
}

#' Run the nested BioGeoBEARS-region/HmscEE-cell occupancy model
#'
#' This is the replacement for the obsolete product model. It implements
#' \deqn{P(cell occupied)=P(region occupied)\times
#' P(cell occupied \mid region occupied)}
#' by conditioning cell-scale dynamics on BSM regional histories. Cross-region
#' transitions are allowed only through changes in `R_region`; normal movement
#' kernels are restricted to cells within the same BioGeoBEARS region.
#'
#' @param suitability Lineage-cell-time suitability table from
#'   `hee_lineage_suitability()` or an equivalent table. Required columns are
#'   `lineage`, `cell_id`, `time_ma`, `region`, `suitability`, and `eta`.
#' @param region_history BioGeoBEARS/BSM region history, validated by
#'   `hee_bsm_region_history()`.
#' @param earth_state Cell-time hard-state table from `hee_paleo_earth_state()`.
#' @param movement_kernel Optional within-region movement kernel. If omitted,
#'   only same-cell carryover is possible between time slices; no silent
#'   long-distance expansion is generated.
#' @param initial_rho Initial occupancy mass for root or region-entry events.
#' @param establishment_intercept,establishment_slope Establishment parameters.
#' @param persistence_intercept,persistence_slope Persistence parameters.
#' @param min_region_q Minimum support enforced when BSM says a region is
#'   occupied but numerical dynamics leave every cell at zero.
#' @param summarise_draws If `TRUE`, summarize `R_region * q` using
#'   `geography_weight * history_weight_given_geography *
#'   climate_weight_given_geography * response_weight`. This is
#'   `sum_g w[g] sum_h w[h|g] sum_c w[c|g] sum_s w[s]`, not an independent
#'   Cartesian average of palaeogeography, palaeoclimate, BioGeoBEARS histories
#'   and HMSC/ancestral-response draws.
#' @return A list with `draws`, `summary`, and `formula_catalog`.
#' @export
hee_nested_region_cell_occupancy <- function(suitability,
                                             region_history,
                                             earth_state,
                                             movement_kernel = NULL,
                                             initial_rho = 0.25,
                                             establishment_intercept = 0,
                                             establishment_slope = 1,
                                             persistence_intercept = 1,
                                             persistence_slope = 1,
                                             min_region_q = 1e-6,
                                             summarise_draws = TRUE) {
  s <- as.data.frame(suitability)
  e <- as.data.frame(earth_state)
  r <- hee_bsm_region_history(region_history)
  .require_cols(s, c("lineage", "cell_id", "time_ma", "suitability", "eta"),
                "suitability")
  .require_cols(e, c("cell_id", "time_ma", "region", "H_state"), "earth_state")
  if (!"geography_scenario" %in% names(s)) {
    s$geography_scenario <- if ("earth_scenario" %in% names(s)) {
      as.character(s$earth_scenario)
    } else "g1"
  }
  if (!"geography_scenario" %in% names(e)) {
    e$geography_scenario <- if ("earth_scenario" %in% names(e)) {
      as.character(e$earth_scenario)
    } else "g1"
  }
  if (!"climate_scenario" %in% names(s)) s$climate_scenario <- "c1"
  if (!"climate_scenario" %in% names(e)) e$climate_scenario <- "c1"
  if (!"earth_scenario" %in% names(s)) s$earth_scenario <- s$geography_scenario
  if (!"earth_scenario" %in% names(e)) e$earth_scenario <- e$geography_scenario
  if (!"response_draw" %in% names(s)) s$response_draw <- "d1"
  if (!"response_weight" %in% names(s)) s$response_weight <- 1
  if (!"climate_weight_given_geography" %in% names(s)) {
    s$climate_weight_given_geography <- 1
  }
  if (!"climate_weight_given_geography" %in% names(e)) {
    e$climate_weight_given_geography <- 1
  }
  .hee_check_unique_keys(s, c("geography_scenario", "climate_scenario",
                              "response_draw", "lineage", "cell_id",
                              "time_ma"), "suitability")
  .hee_check_unique_keys(e, c("geography_scenario", "climate_scenario",
                              "cell_id", "time_ma"),
                         "earth_state")
  base <- merge(s, e[, c("geography_scenario", "climate_scenario",
                         "cell_id", "time_ma", "region", "H_state",
                         "climate_weight_given_geography"),
                     drop = FALSE],
                by = c("geography_scenario", "climate_scenario",
                       "cell_id", "time_ma"),
                all.x = TRUE, sort = FALSE,
                suffixes = c("", ".earth"))
  if ("climate_weight_given_geography.earth" %in% names(base)) {
    base$climate_weight_given_geography <-
      base$climate_weight_given_geography.earth
    base$climate_weight_given_geography.earth <- NULL
  }
  if ("H_state.earth" %in% names(base)) {
    base$H_state <- base$H_state.earth
    base$H_state.earth <- NULL
  }
  for (nm in c("earth_scenario.x", "earth_scenario.y")) {
    if (nm %in% names(base)) base[[nm]] <- NULL
  }
  if ("region.earth" %in% names(base)) {
    if ("region" %in% names(base)) {
      bad_region <- !is.na(base$region) & !is.na(base$region.earth) &
        as.character(base$region) != as.character(base$region.earth)
      if (any(bad_region)) {
        stop("suitability and earth_state disagree on region for ",
             sum(bad_region), " row(s).", call. = FALSE)
      }
    }
    base$region <- as.character(base$region.earth)
    base$region.earth <- NULL
  }
  base$H_state <- as.integer(.hee_hmscee_prob(base$H_state) > 0)
  base$suitability <- .hee_hmscee_prob(base$suitability)
  times <- sort(unique(base$time_ma), decreasing = TRUE)
  .hee_validate_time_values(times, "suitability$time_ma")
  combo <- merge(
    base,
    unique(r[, c("geography_scenario", "history_draw", "lineage",
                 "geography_weight", "history_weight_given_geography"),
             drop = FALSE]),
    by = c("geography_scenario", "lineage"), all.x = TRUE, sort = FALSE
  )
  combo$history_draw[is.na(combo$history_draw)] <- "h1"
  combo$geography_weight[is.na(combo$geography_weight)] <- 0
  combo$history_weight_given_geography[is.na(combo$history_weight_given_geography)] <- 0
  combo <- merge(combo, r, by = c("geography_scenario", "history_draw",
                                  "lineage", "region", "time_ma"),
                 all.x = TRUE, sort = FALSE)
  for (nm in c("geography_weight", "history_weight_given_geography")) {
    y <- paste0(nm, ".y")
    x <- paste0(nm, ".x")
    if (y %in% names(combo)) {
      combo[[nm]] <- combo[[y]]
      combo[[y]] <- NULL
    }
    if (x %in% names(combo)) combo[[x]] <- NULL
  }
  for (nm in c("earth_scenario.x", "earth_scenario.y",
               "earth_weight.x", "earth_weight.y",
               "history_weight_given_earth.x",
               "history_weight_given_earth.y")) {
    if (nm %in% names(combo)) combo[[nm]] <- NULL
  }
  combo$earth_scenario <- combo$geography_scenario
  combo$earth_weight <- combo$geography_weight
  combo$history_weight_given_earth <- combo$history_weight_given_geography
  combo$R_region[is.na(combo$R_region)] <- 0
  combo$R_region <- as.integer(combo$R_region > 0)
  combo$response_weight <- suppressWarnings(as.numeric(combo$response_weight))
  combo$response_weight[!is.finite(combo$response_weight) |
                          combo$response_weight < 0] <- 0
  combo$q <- 0
  combo$arrival_pressure <- 0
  combo$establishment_probability <- 0
  combo$colonisation_probability <- 0
  combo$persistence_probability <- 0
  combo$local_extinction_probability <- 1
  combo$transition_case <- "not_updated"
  combo$regional_support_enforced <- FALSE
  mk <- if (!is.null(movement_kernel)) as.data.frame(movement_kernel) else NULL
  if (!is.null(mk)) {
    .require_cols(mk, c("lineage", "from_cell_id", "to_cell_id", "time_ma",
                        "K_movement"), "movement_kernel")
    mk$K_movement <- .hee_hmscee_prob(mk$K_movement)
  }
  draw_key <- .hee_hmscee_key(combo$geography_scenario,
                              combo$climate_scenario,
                              combo$history_draw,
                              combo$response_draw, combo$lineage)
  for (grp in unique(draw_key)) {
    idx_grp <- which(draw_key == grp)
    ztimes <- times[times %in% combo$time_ma[idx_grp]]
    prev_state <- NULL
    prev_R <- NULL
    for (tt in ztimes) {
      idx <- idx_grp[combo$time_ma[idx_grp] == tt]
      cur <- combo[idx, , drop = FALSE]
      cur$establishment_probability <- hee_establishment_probability(
        cur, intercept = establishment_intercept,
        slope = establishment_slope
      )$establishment_probability
      cur$persistence_probability <- hee_persistence_probability(
        cur, intercept = persistence_intercept,
        slope = persistence_slope
      )$persistence_probability
      cur$local_extinction_probability <- .hee_clip01(1 - cur$persistence_probability)
      absent <- cur$H_state <= 0 | cur$R_region <= 0
      if (is.null(prev_state)) {
        cur$transition_case <- ifelse(absent, "hard_zero",
                                      "initial_region_state")
        if (any(!absent)) {
          entry <- hee_entry_initialization(cur[!absent, , drop = FALSE],
                                            rho = initial_rho)
          cur$q[!absent] <- entry$q_entry
        }
      } else {
        prev_lookup <- stats::setNames(prev_state$q, prev_state$cell_id)
        cur$q_prev_same_cell <- .hee_hmscee_prob(prev_lookup[cur$cell_id])
        prev_R_lookup <- stats::setNames(prev_R$R_region, prev_R$region)
        cur$R_prev <- .hee_hmscee_prob(prev_R_lookup[cur$region])
        entered <- !absent & cur$R_prev <= 0 & cur$R_region > 0
        persisted_region <- !absent & cur$R_prev > 0 & cur$R_region > 0
        cur$transition_case <- ifelse(absent, "hard_zero",
                               ifelse(entered, "BioGeoBEARS_region_entry",
                               ifelse(persisted_region, "within_region_update",
                                      "hard_zero")))
        if (any(entered)) {
          entry <- hee_entry_initialization(cur[entered, , drop = FALSE],
                                            rho = initial_rho)
          cur$q[entered] <- entry$q_entry
        }
        if (any(persisted_region)) {
          if (is.null(mk)) {
            arr <- cur[persisted_region,
                       c("lineage", "cell_id", "time_ma"), drop = FALSE]
            arr$arrival_pressure <- cur$q_prev_same_cell[persisted_region]
          } else {
            kcur <- mk[mk$time_ma == tt &
                         mk$lineage %in% unique(cur$lineage), , drop = FALSE]
            if ("geography_scenario" %in% names(kcur)) {
              kcur <- kcur[kcur$geography_scenario %in%
                             unique(cur$geography_scenario), , drop = FALSE]
            } else if ("earth_scenario" %in% names(kcur)) {
              kcur <- kcur[kcur$earth_scenario %in%
                             unique(cur$geography_scenario), , drop = FALSE]
            }
            if ("climate_scenario" %in% names(kcur)) {
              kcur <- kcur[kcur$climate_scenario %in%
                             unique(cur$climate_scenario), , drop = FALSE]
            }
            arr <- hee_arrival_pressure(
              prev_state[, c("lineage", "cell_id", "q"), drop = FALSE],
              kcur,
              cur[persisted_region,
                  c("lineage", "cell_id", "time_ma"), drop = FALSE]
            )
          }
          arr_key <- .hee_hmscee_key(arr$lineage, arr$cell_id, arr$time_ma)
          cur_key <- .hee_hmscee_key(cur$lineage, cur$cell_id, cur$time_ma)
          cur$arrival_pressure[persisted_region] <-
            arr$arrival_pressure[match(cur_key[persisted_region], arr_key)]
          cur$arrival_pressure[is.na(cur$arrival_pressure)] <- 0
          cur$colonisation_probability[persisted_region] <-
            .hee_clip01(cur$arrival_pressure[persisted_region] *
                          cur$establishment_probability[persisted_region])
          q_new <- cur$q_prev_same_cell[persisted_region] *
            cur$persistence_probability[persisted_region] +
            (1 - cur$q_prev_same_cell[persisted_region]) *
            cur$colonisation_probability[persisted_region]
          cur$q[persisted_region] <- .hee_clip01(q_new)
        }
      }
      cur$q[absent] <- 0
      region_key <- .hee_hmscee_key(cur$geography_scenario,
                                    cur$climate_scenario, cur$history_draw,
                                    cur$response_draw, cur$lineage,
                                    cur$region, cur$time_ma)
      for (rk in unique(region_key[cur$R_region > 0 & cur$H_state > 0])) {
        ridx <- which(region_key == rk & cur$R_region > 0 & cur$H_state > 0)
        if (length(ridx) > 0L && sum(cur$q[ridx], na.rm = TRUE) <= 0 &&
            min_region_q > 0) {
          entry <- hee_entry_initialization(cur[ridx, , drop = FALSE],
                                            rho = min_region_q)
          cur$q[ridx] <- pmax(cur$q[ridx], entry$q_entry)
          cur$regional_support_enforced[ridx] <- TRUE
        }
      }
      combo[idx, names(cur)] <- cur
      prev_state <- cur[, c("lineage", "cell_id", "q"), drop = FALSE]
      prev_R <- unique(cur[, c("region", "R_region"), drop = FALSE])
    }
  }
  combo$cell_probability_given_region <- .hee_clip01(combo$q)
  combo$probability_draw <- .hee_clip01(combo$R_region * combo$q)
  combo$climate_weight_given_geography <- suppressWarnings(as.numeric(
    combo$climate_weight_given_geography
  ))
  combo$climate_weight_given_geography[
    !is.finite(combo$climate_weight_given_geography) |
      combo$climate_weight_given_geography < 0
  ] <- 0
  combo$draw_weight <- combo$geography_weight *
    combo$history_weight_given_geography *
    combo$climate_weight_given_geography * combo$response_weight
  if (!isTRUE(summarise_draws)) {
    return(list(draws = combo, summary = NULL,
                formula_catalog = hee_hmscee_formula_catalog()))
  }
  combo$.summary_key <- .hee_hmscee_key(combo$lineage, combo$cell_id,
                                        combo$time_ma)
  parts <- split(seq_len(nrow(combo)), combo$.summary_key)
  summary <- do.call(rbind, lapply(parts, function(ii) {
    z <- combo[ii, , drop = FALSE]
    w <- .hee_hmscee_normalize_weights(z$draw_weight)
    p <- .hee_hmscee_prob(z$probability_draw)
    qs <- .hee_hmscee_weighted_quantile(p, w, c(0.025, 0.5, 0.975))
    mu <- sum(p * w)
    data.frame(
      lineage = z$lineage[[1]],
      cell_id = z$cell_id[[1]],
      time_ma = z$time_ma[[1]],
      probability_mean = .hee_clip01(mu),
      probability_q025 = .hee_clip01(qs[[1]]),
      probability_q50 = .hee_clip01(qs[[2]]),
      probability_q975 = .hee_clip01(qs[[3]]),
      probability_sd = sqrt(sum(w * (p - mu)^2)),
      n_geography_scenarios = length(unique(z$geography_scenario)),
      n_climate_scenarios = length(unique(z$climate_scenario)),
      n_earth_scenarios = length(unique(z$geography_scenario)),
      n_history_draws = length(unique(.hee_hmscee_key(z$geography_scenario,
                                                      z$history_draw))),
      n_response_draws = length(unique(z$response_draw)),
      stringsAsFactors = FALSE
    )
  }))
  coords <- unique(combo[, intersect(c("geography_scenario",
                                      "climate_scenario", "earth_scenario",
                                      "cell_id", "time_ma",
                                      "region", "lon", "lat", "H_state"),
                                    names(combo)),
                         drop = FALSE])
  if (nrow(coords)) {
    coords <- stats::aggregate(
      coords[, intersect(c("lon", "lat", "H_state"), names(coords)),
             drop = FALSE],
      coords[, c("cell_id", "time_ma"), drop = FALSE],
      function(v) if (is.numeric(v) || is.integer(v)) {
        mean(v, na.rm = TRUE)
      } else {
        paste(unique(as.character(v)), collapse = ";")
      }
    )
    reg <- unique(combo[, intersect(c("cell_id", "time_ma", "region"),
                                   names(combo)), drop = FALSE])
    if (all(c("cell_id", "time_ma", "region") %in% names(reg))) {
      reg <- stats::aggregate(reg$region, reg[, c("cell_id", "time_ma"),
                                              drop = FALSE],
                              function(v) paste(unique(v), collapse = ";"))
      names(reg)[names(reg) == "x"] <- "region"
      coords <- merge(coords, reg, by = c("cell_id", "time_ma"), all.x = TRUE,
                      sort = FALSE)
    }
  }
  if (all(c("cell_id", "time_ma") %in% names(coords))) {
    summary <- merge(summary, coords, by = c("cell_id", "time_ma"),
                     all.x = TRUE, sort = FALSE)
  }
  list(draws = combo, summary = summary,
       formula_catalog = hee_hmscee_formula_catalog())
}

#' Aggregate cell occupancy to BioGeoBEARS-region occupancy
#'
#' Region-level presence is derived from cell probabilities:
#' \deqn{P_{\ell,r,k}=1-\prod_{c\in r}(1-p_{\ell,c,k}).}
#' This replaces independent regional species-pool dynamics in the HmscEE core.
#'
#' @param cell_occupancy Lineage-cell-time table.
#' @param probability_col Probability column.
#' @return Lineage-region-time table.
#' @export
hee_region_pool_from_cell_occupancy <- function(cell_occupancy,
                                                probability_col = "probability_mean") {
  x <- as.data.frame(cell_occupancy)
  .require_cols(x, c("lineage", "region", "time_ma", probability_col),
                "cell_occupancy")
  out <- stats::aggregate(
    .hee_hmscee_prob(x[[probability_col]]),
    x[, c("lineage", "region", "time_ma"), drop = FALSE],
    function(v) .hee_clip01(1 - prod(1 - .hee_hmscee_prob(v), na.rm = TRUE))
  )
  names(out) <- c("lineage", "region", "time_ma", "regional_probability")
  out
}

#' Report extrapolation and uncertainty separately from biology
#'
#' Extrapolation risk is epistemic uncertainty, not absence. This helper joins
#' no-analog or extrapolation diagnostics to a biological probability table
#' without multiplying them into probability.
#'
#' @param probability Biological probability table. It may contain additional
#'   dimensions such as `lineage`, posterior draw, or model id; extrapolation is
#'   joined as a many-to-one cell-time diagnostic.
#' @param extrapolation Extrapolation/no-analog table keyed by cell and time.
#' @param score_col Extrapolation score column.
#' @return Probability table with `extrapolation_score` and `no_analog_flag`;
#'   `probability_mean` or `probability` is unchanged.
#' @export
hee_extrapolation_uncertainty_report <- function(probability,
                                                 extrapolation,
                                                 score_col = NULL) {
  p <- as.data.frame(probability)
  e <- as.data.frame(extrapolation)
  .require_cols(p, c("cell_id", "time_ma"), "probability")
  .require_cols(e, c("cell_id", "time_ma"), "extrapolation")
  if (is.null(score_col)) {
    score_col <- .hee_hmscee_first_col(
      e, c("extrapolation_score", "no_analog_score", "risk", "Q_extrap"),
      "extrapolation"
    )
  }
  .hee_check_unique_keys(e, c("cell_id", "time_ma"), "extrapolation")
  n_before <- nrow(p)
  out <- merge(p, e[, c("cell_id", "time_ma", score_col), drop = FALSE],
               by = c("cell_id", "time_ma"), all.x = TRUE, sort = FALSE)
  if (nrow(out) != n_before) {
    stop("Joining extrapolation diagnostics changed the number of probability rows.",
         call. = FALSE)
  }
  names(out)[names(out) == score_col] <- "extrapolation_score"
  out$extrapolation_score <- .hee_clip01(out$extrapolation_score)
  out$no_analog_flag <- as.integer(is.na(out$extrapolation_score) |
                                     out$extrapolation_score > 0)
  out$probability_interpretation <- "biological_probability_not_downweighted_by_extrapolation"
  out
}

#' HmscEE core process catalog
#'
#' Compatibility entry point for [hee_core_process_catalog()]. The current
#' HmscEE user-facing API names six biological processes: environmental
#' filtering, dispersal, biotic filtering, evolution, speciation, and
#' extinction. Colonisation and persistence remain terms of optional legacy
#' occupancy scenarios, not separate processes in the current framework.
#' Climate, geology,
#' geography, habitat history, BioGeoBEARS/BSM, dated trees, fossils,
#' phyloregion analyses and extrapolation diagnostics remain drivers,
#' constraints, observations, uncertainty layers, or result summaries.
#'
#' @return A data.frame from [hee_core_process_catalog()].
#' @export
#'
#' @examples
#' hee_hmscee_process_catalog()
hee_hmscee_process_catalog <- function() {
  hee_core_process_catalog()
}

#' HmscEE nested framework formula catalog
#'
#' Returns the equations and interpretation boundaries of the current
#' six-process HmscEE core. The obsolete `A_hist * D_dynamic * Q_extrap`
#' product and uncalibrated occupancy transitions are deliberately absent.
#'
#' @return A data.frame with formula ids, equations, scale, and interpretation.
#' @export
#'
#' @examples
#' hee_hmscee_formula_catalog()
hee_hmscee_formula_catalog <- function() {
  row <- function(id, formula, scale, interpretation, estimate_type,
                  framework_layer, core_process) {
    data.frame(formula_id = id, formula = formula, scale = scale,
               interpretation = interpretation, estimate_type = estimate_type,
               framework_layer = framework_layer, core_process = core_process,
               stringsAsFactors = FALSE)
  }
  do.call(rbind, list(
    row("paleo_earth_state", "E[k] = {X[k], H[k], W[k]}", "cell-time",
        "Palaeoclimate, land, habitat and resistance are upstream conditions, not extra biological probabilities.",
        "external_state", "external_drivers", "environmental_filtering/dispersal"),
    row("ancestral_environmental_response", "beta[lineage,time] is reconstructed from HMSC tip Beta draws on the dated tree",
        "lineage-time", "Ancestral response is not an average of modern maps or a geographic range.",
        "posterior_draw", "core_process", "evolution"),
    row("trait_mediated_environmental_response", "beta[lineage,time] = Gamma T[lineage,time] + u[lineage,time]",
        "lineage-time", "Optional trait-mediated model retains an unexplained phylogenetic response term.",
        "posterior_or_sensitivity", "core_process", "evolution"),
    row("environmental_support", "S_env[lineage,cell,time] = link^{-1}(b(X[cell,time])^T beta[lineage,time])",
        "lineage-cell-time", "S_env measures environmental support, not historical occupancy.",
        "posterior_prediction", "core_process", "environmental_filtering"),
    row("movement_kernel", "K[from,to,time] is a land-constrained directional active-movement kernel",
        "cell-transition", "Plate carriage changes the spatial reference frame; it is not active dispersal.",
        "model_or_sensitivity_kernel", "core_process", "dispersal"),
    row("biotic_effect", "B[lineage,cell,time] = 0 unless independent interaction evidence supports a term",
        "lineage-cell-time", "HMSC residual association alone is not a causal interaction coefficient.",
        "external_or_unparameterised", "core_process", "biotic_filtering"),
    row("speciation_identity", "active_lineages[time] follow dated-tree branch and node times",
        "lineage-time", "Descendant lineages cannot occur before their split; this does not estimate speciation rates.",
        "dated_tree_constraint", "core_process", "speciation"),
    row("extinction_scope", "geographic_loss != local_population_extirpation != lineage_extinction",
        "cell-or-lineage", "Marine or vanished-land sink mass is not a measured biological extinction rate.",
        "scope_diagnostic", "core_process", "extinction"),
    row("regional_history", "R[lineage,region,time]^(history,geography) is optional BSM information conditional on geography",
        "lineage-region-time", "BioGeoBEARS is an optional historical constraint, not a seventh process or an accessibility multiplier.",
        "external_conditional_history", "historical_inference_or_constraint", "dispersal/speciation"),
    row("final_probability", "p(location at time | tip locations, dated tree, dynamic Earth, K) is inferred by time-ordered pruning",
        "lineage-cell-time", "Current global-grid Case05 yields conditional lineage-location mass, not demographic occupancy or historical species richness.",
        "conditional_location_posterior", "integrated_result", "dispersal/speciation"),
    row("extrapolation_uncertainty", "No-analog diagnostics are reported beside S_env and p, never multiplied into them",
        "epistemic_uncertainty", "Reliability evidence guides interpretation; it is not biological absence.",
        "diagnostic_not_probability_multiplier", "observation_uncertainty", "not_core_process")
  ))
}

#' Dynamic Earth-Biota formula catalog
#'
#' Compatibility alias for `hee_hmscee_formula_catalog()`. The returned catalog
#' describes the current nested HmscEE framework, not the obsolete product
#' model.
#'
#' @export
hee_earth_biota_formula_catalog <- function() hee_hmscee_formula_catalog()

# Deprecated/removed old product-core entry points -------------------------

#' Obsolete product-form Dynamic Earth-Biota core
#'
#' These functions belonged to the obsolete multiplicative HmscEE core
#' `S_HMSC * A_hist * E_phylo * L_land * D_dynamic * Q_extrap`. They now stop
#' with an explanatory error because the current HmscEE core uses nested
#' BioGeoBEARS/BSM regional history plus within-region dynamic occupancy.
#'
#' @param ... Obsolete compatibility arguments. They are not evaluated.
#' @return These functions do not return a value; they stop with guidance to
#'   use `hee_nested_region_cell_occupancy()`.
#' @name hee_obsolete_product_core
#' @export
hee_dynamic_earth_biota_probability <- function(...) {
  .hee_stop_obsolete_product_core("hee_dynamic_earth_biota_probability")
}

#' @rdname hee_obsolete_product_core
#' @export
hee_static_snapshot_approximation <- function(...) {
  .hee_stop_obsolete_product_core("hee_static_snapshot_approximation")
}

#' @rdname hee_obsolete_product_core
#' @export
hee_dynamic_state_transition <- function(...) {
  .hee_stop_obsolete_product_core("hee_dynamic_state_transition")
}

#' @rdname hee_obsolete_product_core
#' @export
hee_regional_species_pool_transition <- function(...) {
  .hee_stop_obsolete_product_core("hee_regional_species_pool_transition")
}

#' @rdname hee_paleo_earth_state
#' @param arena Compatibility input table for `hee_geographic_stage()`.
#' @param output_col Compatibility output column name, usually `G_arena`.
#' @param ... Additional arguments passed to `hee_paleo_earth_state()`.
#' @export
hee_geographic_stage <- function(arena, output_col = "G_arena", ...) {
  warning(
    "`hee_geographic_stage()` is a compatibility alias. The HmscEE core now ",
    "uses hard-state `H_state` from `hee_paleo_earth_state()` rather than a ",
    "geographic probability multiplier.",
    call. = FALSE
  )
  x <- as.data.frame(arena)
  if (!"region" %in% names(x)) x$region <- "__arena__"
  out <- hee_paleo_earth_state(x, ...)
  out[[output_col]] <- out$H_state
  if (!"L_land" %in% names(out)) out$L_land <- out$land_state
  out
}
