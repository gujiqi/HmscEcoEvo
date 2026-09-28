#' Calculate speciation-related historical indices
#'
#' @param proj A `HmscEcoEvo_data` object.
#' @param probability_weighted Logical; use event probabilities when available.
#' @return A list with `XData` and `diagnostics`.
#' @export
calc_speciation_indices <- function(proj, probability_weighted = TRUE) {
  comm <- proj$comm
  sr <- proj$site_region
  ev <- proj$history_tables$speciation_events
  srh <- proj$history_tables$species_region_history
  sites <- rownames(comm)
  out <- data.frame(row.names = sites)
  messages <- warnings <- character()

  if (is.null(ev) || nrow(ev) == 0) {
    messages <- c(messages, "speciation_events not available; speciation indices skipped.")
    return(list(XData = out, diagnostics = list(messages = messages, warnings = warnings)))
  }

  ev$event_prob <- if ("event_prob" %in% names(ev) && probability_weighted) .std_prob(ev$event_prob, 1) else rep(1, nrow(ev))
  out$insitu_speciation_prop <- NA_real_
  out$insitu_diversification_strength <- NA_real_

  for (site in sites) {
    reg <- sr$region[match(site, sr$site)]
    spp <- names(which(comm[site, ] > 0))
    tab <- ev[ev$descendant_species %in% spp, , drop = FALSE]
    if (nrow(tab) == 0) next
    total <- sum(tab$event_prob, na.rm = TRUE)
    insitu <- sum(tab$event_prob[tab$region == reg], na.rm = TRUE)
    out[site, "insitu_speciation_prop"] <- if (total > 0) insitu / total else NA_real_
    out[site, "insitu_diversification_strength"] <- if (length(spp) > 0) insitu / length(spp) else NA_real_
  }

  # Ex situ colonization from species_region_history when available.
  if (!is.null(srh) && nrow(srh) > 0 && "source_region" %in% names(srh)) {
    out$exsitu_colonization_prop <- NA_real_
    for (site in sites) {
      reg <- sr$region[match(site, sr$site)]
      spp <- names(which(comm[site, ] > 0))
      tab <- srh[srh$species %in% spp & srh$region == reg, , drop = FALSE]
      if (nrow(tab) == 0) next
      w <- if (probability_weighted && "entry_prob" %in% names(tab)) .std_prob(tab$entry_prob, 1) else rep(1, nrow(tab))
      ex <- !is.na(tab$source_region) & tab$source_region != "" & tab$source_region != reg
      out[site, "exsitu_colonization_prop"] <- sum(w[ex], na.rm = TRUE) / sum(w, na.rm = TRUE)
    }
    out$insitu_exsitu_balance <- out$insitu_speciation_prop - out$exsitu_colonization_prop
  } else {
    messages <- c(messages, "species_region_history with source_region not available; insitu_exsitu_balance skipped.")
  }

  warnings <- c(warnings, "Speciation indices depend on upstream phylogenetic and ancestral-area resolution; interpret cautiously with heavily imputed trees.")
  out <- .drop_all_na_cols(out)
  list(XData = out, diagnostics = list(messages = unique(messages), warnings = unique(warnings)))
}
