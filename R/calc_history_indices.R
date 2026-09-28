#' Calculate HMSC-ready historical indices
#'
#' Calculates site-level and species-level historical indices from standardized
#' history tables stored in a `HmscEcoEvo_data` object. Event tables are optional:
#' a single `species_region_history` table is enough for the basic workflow.
#'
#' @param proj A `HmscEcoEvo_data` object with standardized history tables.
#' @param compute Character vector of groups to compute or "all". Supported
#'   groups are "occupancy_dynamics", "speciation", "dispersal_sources", and
#'   "trait_history".
#' @param probability_weighted Logical; use probability columns when available.
#' @param shannon_base Base for Shannon diversity; default `exp(1)`.
#' @param x A historical-index object.
#' @param object A historical-index object.
#' @param ... Additional arguments, currently ignored.
#' @return An object of class `c("HmscEcoEvo_history_indices", "hmsc_history_indices")`.
#' @export
calc_history_indices <- function(proj,
                                 compute = "all",
                                 probability_weighted = TRUE,
                                 shannon_base = exp(1)) {
  if (!inherits(proj, "hmscHist_data")) {
    stop("proj must be an hmscHist_data object.", call. = FALSE)
  }
  groups <- c("occupancy_dynamics", "speciation", "dispersal_sources", "trait_history")
  if (length(compute) == 1 && compute == "all") compute <- groups
  bad <- setdiff(compute, groups)
  if (length(bad) > 0) stop("Unknown compute group(s): ", paste(bad, collapse = ", "), call. = FALSE)

  tabs <- proj$history_tables
  sites <- rownames(proj$comm)
  species <- colnames(proj$comm)

  XData_history <- data.frame(row.names = sites)
  TrData_history <- NULL
  random_history <- data.frame(row.names = sites)
  messages <- character()
  warnings <- character()

  if ("occupancy_dynamics" %in% compute) {
    occ <- calc_occupancy_dynamics_indices(proj, probability_weighted = probability_weighted)
    XData_history <- .cbind_keep(XData_history, occ$XData)
    messages <- c(messages, occ$diagnostics$messages)
    warnings <- c(warnings, occ$diagnostics$warnings)
  }

  if ("speciation" %in% compute) {
    spc <- calc_speciation_indices(proj, probability_weighted = probability_weighted)
    XData_history <- .cbind_keep(XData_history, spc$XData)
    messages <- c(messages, spc$diagnostics$messages)
    warnings <- c(warnings, spc$diagnostics$warnings)
  }

  if ("dispersal_sources" %in% compute) {
    dsp <- calc_dispersal_source_indices(proj, probability_weighted = probability_weighted,
                                         shannon_base = shannon_base)
    XData_history <- .cbind_keep(XData_history, dsp$XData)
    random_history <- .cbind_keep(random_history, dsp$random)
    messages <- c(messages, dsp$diagnostics$messages)
    warnings <- c(warnings, dsp$diagnostics$warnings)
  }

  if ("trait_history" %in% compute) {
    tr <- calc_trait_history_indices(proj, probability_weighted = probability_weighted)
    TrData_history <- tr$TrData
    messages <- c(messages, tr$diagnostics$messages)
    warnings <- c(warnings, tr$diagnostics$warnings)
  }

  groups_out <- list(
    occupancy_dynamics = intersect(c(
      "colonization_age", "colonization_age_dispersion",
      "lineage_retention_index", "range_loss_fraction",
      "historical_gain_fraction", "historical_loss_fraction",
      "historical_turnover_predictor", "historical_gain_loss_balance",
      "historical_occupancy_balance"
    ), names(XData_history)),
    speciation = intersect(c(
      "insitu_speciation_prop", "insitu_diversification_strength", "insitu_exsitu_balance"
    ), names(XData_history)),
    dispersal_sources = intersect(c(
      "prob_weighted_source_diversity", "dispersal_source_diversity",
      "dominant_source_prop", "historical_inflow_outflow_balance"
    ), names(XData_history)),
    trait_history = if (is.null(TrData_history)) character() else names(TrData_history),
    random_history = names(random_history)
  )

  diagnostics <- list(
    messages = unique(messages),
    warnings = unique(warnings),
    computed_groups = names(groups_out)[vapply(groups_out, length, integer(1)) > 0],
    skipped_groups = setdiff(groups, names(groups_out)[vapply(groups_out, length, integer(1)) > 0])
  )

  out <- list(
    XData_history = XData_history,
    TrData_history = TrData_history,
    random_history = random_history,
    groups = groups_out,
    diagnostics = diagnostics,
    dictionary = history_dictionary(),
    settings = list(compute = compute, probability_weighted = probability_weighted,
                    shannon_base = shannon_base),
    project = proj
  )
  class(out) <- c("HmscEcoEvo_history_indices", "hmsc_history_indices")
  out
}

#' @rdname calc_history_indices
#' @export
print.HmscEcoEvo_history_indices <- function(x, ...) {
  cat("HmscEcoEvo_history_indices object\n")
  cat("  site-level historical variables: ", ncol(x$XData_history), "\n", sep = "")
  cat("  species-level historical variables: ", if (is.null(x$TrData_history)) 0 else ncol(x$TrData_history), "\n", sep = "")
  cat("  random historical variables: ", ncol(x$random_history), "\n", sep = "")
  if (length(x$diagnostics$warnings) > 0) {
    cat("  warnings: ", length(x$diagnostics$warnings), "\n", sep = "")
  }
  invisible(x)
}

#' @rdname calc_history_indices
#' @export
summary.HmscEcoEvo_history_indices <- function(object, ...) {
  print(object)
  cat("\nGroups:\n")
  for (nm in names(object$groups)) {
    vals <- object$groups[[nm]]
    cat("  ", nm, ": ", if (length(vals) == 0) "none" else paste(vals, collapse = ", "), "\n", sep = "")
  }
  if (length(object$diagnostics$messages) > 0) {
    cat("\nMessages:\n")
    cat(paste0("  - ", object$diagnostics$messages, collapse = "\n"), "\n")
  }
  if (length(object$diagnostics$warnings) > 0) {
    cat("\nWarnings:\n")
    cat(paste0("  - ", object$diagnostics$warnings, collapse = "\n"), "\n")
  }
  invisible(object)
}

#' @rdname calc_history_indices
#' @export
print.hmsc_history_indices <- print.HmscEcoEvo_history_indices

#' @rdname calc_history_indices
#' @export
summary.hmsc_history_indices <- summary.HmscEcoEvo_history_indices
