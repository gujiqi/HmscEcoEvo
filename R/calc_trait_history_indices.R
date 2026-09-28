#' Calculate species-level trait history indices
#'
#' @param proj A `HmscEcoEvo_data` object.
#' @param probability_weighted Logical; use probabilities where available.
#' @return A list with `TrData` and `diagnostics`.
#' @export
calc_trait_history_indices <- function(proj, probability_weighted = TRUE) {
  species <- colnames(proj$comm)
  th <- proj$history_tables$trait_history
  srh <- proj$history_tables$species_region_history
  se <- proj$history_tables$speciation_events
  messages <- warnings <- character()

  if (is.null(th) || nrow(th) == 0) {
    messages <- c(messages, "trait_history not available; TrData_history skipped.")
    return(list(TrData = NULL, diagnostics = list(messages = messages, warnings = warnings)))
  }
  if (!"species" %in% names(th)) th$species <- rownames(th)
  rownames(th) <- th$species
  Tr <- th[intersect(species, rownames(th)), setdiff(names(th), "species"), drop = FALSE]

  # Coupling: dispersal trait x number of source regions.
  if ("dispersal_trait" %in% names(Tr) && !is.null(srh) && "source_region" %in% names(srh)) {
    nsrc <- tapply(srh$source_region, srh$species, function(z) length(unique(na.omit(z[z != ""]))))
    nsrc <- nsrc[rownames(Tr)]
    Tr$trait_dispersal_coupling <- as.numeric(.z(Tr$dispersal_trait) * .z(as.numeric(nsrc)))
  } else {
    messages <- c(messages, "trait_dispersal_coupling skipped; requires dispersal_trait and source_region history.")
  }

  # Coupling: trait shift magnitude x species-level in situ speciation proportion.
  if ("trait_shift_magnitude" %in% names(Tr) && !is.null(se) && nrow(se) > 0) {
    se$event_prob <- if ("event_prob" %in% names(se) && probability_weighted) .std_prob(se$event_prob, 1) else rep(1, nrow(se))
    total <- tapply(se$event_prob, se$descendant_species, sum, na.rm = TRUE)
    insitu <- total # in the standardized table speciation_events are in situ by construction in v0.1
    prop <- insitu / total
    prop <- prop[rownames(Tr)]
    Tr$trait_insitu_coupling <- as.numeric(.z(Tr$trait_shift_magnitude) * .z(as.numeric(prop)))
  } else {
    messages <- c(messages, "trait_insitu_coupling skipped; requires trait_shift_magnitude and speciation_events.")
  }

  Tr <- .drop_all_na_cols(Tr)
  list(TrData = Tr, diagnostics = list(messages = unique(messages), warnings = unique(warnings)))
}
