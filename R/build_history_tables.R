#' Build standardized history tables
#'
#' Build the four internal history tables used by HmscEcoEvo. The function can
#' standardize user-provided tables, or it can try to extract event-like tables
#' from supplied BioGeoBEARS/BSM outputs. It does not infer ancestral ranges and
#' does not reimplement BioGeoBEARS.
#'
#' @param proj A `HmscEcoEvo_data` object.
#' @param biogeobears Optional BioGeoBEARS result object. Direct support is limited;
#'   event-like data.frames are searched recursively.
#' @param bsm_events Optional BioGeoBEARS stochastic-map/BSM events as a data.frame
#'   or a list of data.frames.
#' @param tip_ranges Optional present species by region ranges, as a matrix or
#'   data.frame. Used to infer `retained`.
#' @param node_ranges Optional node range summaries. Currently only stored/checked;
#'   raw ancestral reconstruction is not performed.
#' @param species_region_history Optional manual species-region history table.
#' @param speciation_events Optional manual speciation event table.
#' @param history_events Optional manual dispersal/loss event table.
#' @param traits Optional species trait table.
#' @param ancestral_traits Optional precomputed ancestral trait summaries.
#' @param trait_history Optional manual trait history table.
#' @param method One of "auto", "bsm", "node", or "manual".
#' @param probability_weighted Logical; retain probability columns and use them later.
#' @param time_unit Character; default "Ma".
#' @return Updated `HmscEcoEvo_data` object with `history_tables`.
#' @export
build_history_tables <- function(proj,
                                 biogeobears = NULL,
                                 bsm_events = NULL,
                                 tip_ranges = NULL,
                                 node_ranges = NULL,
                                 species_region_history = NULL,
                                 speciation_events = NULL,
                                 history_events = NULL,
                                 traits = NULL,
                                 ancestral_traits = NULL,
                                 trait_history = NULL,
                                 method = c("auto", "bsm", "node", "manual"),
                                 probability_weighted = TRUE,
                                 time_unit = "Ma") {
  if (!inherits(proj, "hmscHist_data")) stop("proj must be an hmscHist_data object.", call. = FALSE)
  method <- match.arg(method)
  messages <- character()
  warnings <- character()

  raw_events <- NULL
  if (method %in% c("auto", "bsm")) {
    if (!is.null(bsm_events) || !is.null(biogeobears)) {
      raw_events <- tryCatch(
        extract_biogeobears_events(biogeobears = biogeobears, bsm_events = bsm_events),
        error = function(e) {
          warnings <<- c(warnings, paste("Could not extract BioGeoBEARS/BSM events:", conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(raw_events)) messages <- c(messages, "Extracted event table from supplied BioGeoBEARS/BSM-like input.")
    }
  }

  if (method == "node" && is.null(node_ranges)) {
    warnings <- c(warnings, "method = 'node' requested but node_ranges is NULL. No node-based tables were built.")
  }

  species <- colnames(proj$comm)
  regions <- unique(proj$site_region$region)

  if (is.null(species_region_history) && !is.null(raw_events)) {
    species_region_history <- build_species_region_history(
      raw_events = raw_events,
      tip_ranges = tip_ranges,
      species = species,
      regions = regions
    )
    messages <- c(messages, "Built species_region_history from standardized event table.")
  }

  if (is.null(speciation_events) && !is.null(raw_events)) {
    speciation_events <- build_speciation_events(raw_events = raw_events, phy = proj$phy)
    if (!is.null(speciation_events) && nrow(speciation_events) > 0) {
      messages <- c(messages, "Built speciation_events from standardized event table.")
    }
  }

  if (is.null(history_events) && !is.null(raw_events)) {
    history_events <- build_history_events(raw_events = raw_events)
    if (!is.null(history_events) && nrow(history_events) > 0) {
      messages <- c(messages, "Built history_events from standardized event table.")
    }
  }

  if (is.null(trait_history) && !is.null(ancestral_traits)) {
    trait_history <- build_trait_history(
      traits = NULL,
      ancestral_traits = ancestral_traits,
      trait_history = NULL,
      species = species
    )
    if (!is.null(trait_history) && nrow(trait_history) > 0) {
      messages <- c(messages, "Built trait_history from supplied ancestral_trait summaries.")
    }
  }
  if (is.null(trait_history) && !is.null(traits)) {
    messages <- c(messages, "traits are stored as ordinary HMSC traits; no trait_history table was created unless trait_history or ancestral_traits is supplied.")
  }

  proj$history_tables <- list(
    species_region_history = species_region_history,
    speciation_events = speciation_events,
    history_events = history_events,
    trait_history = trait_history
  )
  proj$history_tables <- proj$history_tables[!vapply(proj$history_tables, is.null, logical(1))]

  proj <- standardize_history_tables(proj, time_unit = time_unit)
  proj$history_diagnostics$messages <- unique(c(proj$history_diagnostics$messages, messages))
  proj$history_diagnostics$warnings <- unique(c(proj$history_diagnostics$warnings, warnings))
  proj$settings <- modifyList(proj$settings %||% list(), list(
    method = method,
    probability_weighted = probability_weighted,
    time_unit = time_unit
  ))
  proj
}

`%||%` <- function(a, b) if (is.null(a)) b else a
