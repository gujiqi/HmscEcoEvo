#' Legacy M1-M5 static projection scenario formulas
#'
#' Return the five older HmscEcoEvo product-form projection scenarios in a
#' reader-facing table. These formulas are retained for backward-compatible
#' Case 1/Case 2 sensitivity contrasts and reports. They are not the HmscEE
#' 1.0 Dynamic Earth-Biota Assembly core, where BioGeoBEARS/BSM regional
#' histories are conditional on palaeogeographic/plate scenarios `g` and
#' HmscEE updates
#' cell occupancy inside occupied regions with [hee_nested_region_cell_occupancy()].
#'
#' @param include_q Logical. Include `Q_extrap` in the written formula.
#' @param extrapolation_mode How `Q_extrap` is interpreted. `"flag"` means it
#'   is reported as reliability information but not used to reduce probability;
#'   `"downweight"` means it enters as a multiplicative weight.
#'
#' @return A data.frame with legacy scenario model definitions.
#' @export
#'
#' @examples
#' hee_projection_formula()
hee_projection_formula <- function(include_q = TRUE,
                                   extrapolation_mode = c("flag", "downweight")) {
  extrapolation_mode <- match.arg(extrapolation_mode)
  q <- if (include_q && extrapolation_mode == "downweight") " x Q_extrap" else ""
  q_note <- if (include_q && extrapolation_mode == "flag") {
    "Q_extrap is reported as an extrapolation/reliability flag, not multiplied into P."
  } else if (include_q) {
    "Q_extrap is multiplied into P as an explicit extrapolation downweight."
  } else {
    "No extrapolation term is included."
  }
  data.frame(
    model_id = c(
      "M1_env_only",
      "M2_env_phylo",
      "M3_env_phylo_bgb_no_dispersal",
      "M4_env_phylo_bgb_static_dispersal",
      "M5_env_phylo_bgb_dynamic_dispersal"
    ),
    formula = c(
      paste0("P = S_HMSC x L_land", q),
      paste0("P = S_HMSC x E_phylo x L_land", q),
      paste0("P = S_HMSC x E_phylo x A_BGB x L_land", q),
      paste0("P = S_HMSC x E_phylo x A_BGB x D_static x L_land", q),
      paste0("P = S_HMSC x E_phylo x A_BGB x D_dynamic x L_land", q)
    ),
    added_constraint = c(
      "environment plus land mask only",
      "adds lineage/species existence through time",
      "adds historical regional accessibility",
      "adds static distance/connectivity dispersal",
      "adds dynamic previous-time dispersal continuity"
    ),
    interpretation = c(
      "Where would the environment be suitable on land?",
      "Where is the environment suitable and the species/lineage already exists?",
      "Where is it suitable, extant in time, and historically reachable at regional scale?",
      "Where is it suitable, reachable, and near/connectable enough under a static dispersal kernel?",
      "Where is it suitable, reachable, and dynamically reachable from the previous time slice?"
    ),
    caveat = c(
      "Over-broad because it ignores lineage age and accessibility.",
      "Still ignores historical range accessibility and dispersal limits.",
      "Good conservative baseline, but region-level accessibility is not cell-level dispersal.",
      "Depends on dispersal scale, distances, and barrier assumptions.",
      "Most mechanistic, but most assumption-dependent; report direction and time step clearly."
    ),
    extrapolation_policy = q_note,
    stringsAsFactors = FALSE
  )
}

#' Convert extrapolation risk to `Q_extrap`
#'
#' HmscEcoEvo treats extrapolation risk as a reliability diagnostic by default.
#' Use `mode = "downweight"` only when the user explicitly wants risk to reduce
#' final probabilities.
#'
#' @param extrapolation Extrapolation table with `cell_id`, `time_ma`, and a
#'   score column.
#' @param mode `"flag"` keeps `Q_extrap = 1` and reports risk; `"downweight"`
#'   sets `Q_extrap = exp(-lambda * max(score - threshold, 0))`. Missing,
#'   `NaN`, and negative scores are treated as neutral zero risk; `Inf` is
#'   treated as extreme extrapolation risk and downweights to zero.
#' @param lambda Downweight strength for `mode = "downweight"`.
#' @param threshold Risk score below which no downweight is applied.
#' @param score_col Score column.
#'
#' @return A data.frame with `Q_extrap`, `extrapolation_score`, and flags.
#' @export
#'
#' @examples
#' x <- data.frame(cell_id = 1:3, time_ma = 0, extrapolation_score = c(0, 1, 2))
#' hee_extrapolation_weight(x)
hee_extrapolation_weight <- function(extrapolation,
                                     mode = c("flag", "downweight"),
                                     lambda = 1,
                                     threshold = 0,
                                     score_col = "extrapolation_score") {
  mode <- match.arg(mode)
  x <- as.data.frame(extrapolation)
  .require_cols(x, c("cell_id", "time_ma", score_col), "extrapolation")
  score <- suppressWarnings(as.numeric(x[[score_col]]))
  score[is.na(score) | is.nan(score)] <- 0
  score[is.infinite(score) & score < 0] <- 0
  score[is.infinite(score) & score > 0] <- .Machine$double.xmax
  score[score < 0] <- 0
  q <- if (mode == "flag") {
    rep(1, nrow(x))
  } else {
    exp(-lambda * pmax(score - threshold, 0))
  }
  out <- x[, intersect(c("cell_id", "time_ma", "lon", "lat", score_col,
                         "extrapolation_flag", "no_analog_fraction",
                         "variable_out_of_range"), names(x)), drop = FALSE]
  names(out)[names(out) == score_col] <- "extrapolation_score"
  out$extrapolation_score <- score
  out$Q_extrap <- pmin(pmax(q, 0), 1)
  out$extrapolation_mode <- mode
  out$extrapolation_used_as_downweight <- mode == "downweight"
  out
}

#' Work through the five projection models for one numeric example
#'
#' This helper is deliberately simple so the result can be placed in reports.
#' It shows how a high environmental suitability can become low final
#' probability once historical accessibility or dispersal are considered.
#'
#' @param S_HMSC Environmental suitability.
#' @param A_BGB Historical accessibility.
#' @param E_phylo Phylogenetic/time existence mask.
#' @param L_land Land mask.
#' @param D_static Static dispersal weight.
#' @param D_dynamic Dynamic dispersal weight.
#' @param Q_extrap Optional extrapolation weight.
#' @param extrapolation_mode Whether `Q_extrap` is a flag or downweight.
#'
#' @return A data.frame with M1-M5 probabilities.
#' @export
#'
#' @examples
#' hee_projection_component_example(S_HMSC = 0.8, A_BGB = 0.25,
#'   E_phylo = 1, L_land = 1, D_static = 0.6, D_dynamic = 0.1)
hee_projection_component_example <- function(S_HMSC = 0.8,
                                             A_BGB = 0.25,
                                             E_phylo = 1,
                                             L_land = 1,
                                             D_static = 0.6,
                                             D_dynamic = 0.1,
                                             Q_extrap = 1,
                                             extrapolation_mode = c("flag", "downweight")) {
  extrapolation_mode <- match.arg(extrapolation_mode)
  q <- if (extrapolation_mode == "downweight") Q_extrap else 1
  formula <- hee_projection_formula(include_q = TRUE, extrapolation_mode = extrapolation_mode)
  probability <- c(
    S_HMSC * L_land * q,
    S_HMSC * E_phylo * L_land * q,
    S_HMSC * E_phylo * A_BGB * L_land * q,
    S_HMSC * E_phylo * A_BGB * D_static * L_land * q,
    S_HMSC * E_phylo * A_BGB * D_dynamic * L_land * q
  )
  data.frame(
    model_id = formula$model_id,
    S_HMSC = S_HMSC,
    E_phylo = c(1, E_phylo, E_phylo, E_phylo, E_phylo),
    A_BGB = c(1, 1, A_BGB, A_BGB, A_BGB),
    D_static = c(1, 1, 1, D_static, 1),
    D_dynamic = c(1, 1, 1, 1, D_dynamic),
    L_land = L_land,
    Q_extrap = q,
    probability = pmin(pmax(probability, 0), 1),
    formula = formula$formula,
    interpretation = formula$interpretation,
    stringsAsFactors = FALSE
  )
}

#' Combine S, A, E, L, D, and Q into legacy M1-M5 projections
#'
#' This is a convenience wrapper around [hee_combine()] that returns a single
#' long table for the older product-form M1-M5 scenario family. Use it for
#' transparent backward-compatible sensitivity contrasts, not as the HmscEE 1.0
#' nested core. For the main Dynamic Earth-Biota Assembly model use
#' [hee_paleo_earth_state()], [hee_bsm_region_history()],
#' [hee_lineage_suitability()], [hee_within_region_movement_kernel()], and
#' [hee_nested_region_cell_occupancy()].
#'
#' @param suitability HMSC environmental suitability table.
#' @param accessibility Optional historical accessibility table.
#' @param phylo_mask Optional species-by-time mask.
#' @param landmask Optional land mask table.
#' @param static_filter Optional static dispersal filter.
#' @param dynamic_filter Optional dynamic dispersal filter.
#' @param extrapolation_weight Optional `Q_extrap` table from
#'   `hee_extrapolation_weight()`.
#' @param strict Logical. When `TRUE` (default), the function refuses to label
#'   M1-M5 as complete unless the formula-required land, phylogenetic,
#'   accessibility, static-dispersal, and dynamic-dispersal components are
#'   supplied. Use `strict = FALSE` only for formula demonstrations or deliberate
#'   sensitivity scenarios, and report the omitted components as neutral.
#'
#' @return A long data.frame with `model_id` and component columns.
#' @export
hee_combine_projection_models <- function(suitability,
                                          accessibility = NULL,
                                          phylo_mask = NULL,
                                          landmask = NULL,
                                          static_filter = NULL,
                                          dynamic_filter = NULL,
                                          extrapolation_weight = NULL,
                                          strict = TRUE) {
  if (isTRUE(strict)) {
    missing_components <- c(
      if (is.null(landmask)) "landmask/L_land for M1-M5" else character(),
      if (is.null(phylo_mask)) "phylo_mask/E_phylo for M2-M5" else character(),
      if (is.null(accessibility)) "accessibility/A_BGB for M3-M5" else character(),
      if (is.null(static_filter)) "static_filter/D_static for M4" else character(),
      if (is.null(dynamic_filter)) "dynamic_filter/D_dynamic for M5" else character()
    )
    if (length(missing_components) > 0L) {
      stop("Strict M1-M5 projection requires all named process components. ",
           "Missing: ", paste(missing_components, collapse = "; "),
           ". Use strict = FALSE only for explicit neutral-component ",
           "demonstrations, not as complete deep-time results.",
           call. = FALSE)
    }
  }
  extrap_mode <- "flag"
  if (!is.null(extrapolation_weight)) {
    ew <- as.data.frame(extrapolation_weight)
    if ("extrapolation_used_as_downweight" %in% names(ew) &&
        any(ew$extrapolation_used_as_downweight %in% TRUE, na.rm = TRUE)) {
      extrap_mode <- "downweight"
    }
  }
  defs <- hee_projection_formula(
    include_q = !is.null(extrapolation_weight),
    extrapolation_mode = extrap_mode
  )
  make <- function(model_id, access = NULL, phylo = NULL, static = NULL, dynamic = NULL) {
    p <- hee_combine(suitability, accessibility = access, phylo_mask = phylo,
                     landmask = landmask, static_filter = static,
                     dynamic_filter = dynamic,
                     extrapolation_weight = extrapolation_weight)
    p$model_id <- model_id
    p
  }
  out <- rbind(
    make(defs$model_id[1]),
    make(defs$model_id[2], phylo = phylo_mask),
    make(defs$model_id[3], access = accessibility, phylo = phylo_mask),
    make(defs$model_id[4], access = accessibility, phylo = phylo_mask, static = static_filter),
    make(defs$model_id[5], access = accessibility, phylo = phylo_mask, dynamic = dynamic_filter)
  )
  out <- merge(out, defs[, c("model_id", "formula", "interpretation")],
               by = "model_id", all.x = TRUE, sort = FALSE)
  out
}
