#' Calculate historical dispersal source indices
#'
#' @param proj A `HmscEcoEvo_data` object.
#' @param probability_weighted Logical; use entry probabilities when available.
#' @param shannon_base Base for Shannon diversity.
#' @return A list with `XData`, `random`, and `diagnostics`.
#' @export
calc_dispersal_source_indices <- function(proj, probability_weighted = TRUE,
                                          shannon_base = exp(1)) {
  comm <- proj$comm
  sr <- proj$site_region
  srh <- proj$history_tables$species_region_history
  hev <- proj$history_tables$history_events
  sites <- rownames(comm)
  out <- data.frame(row.names = sites)
  ran <- data.frame(row.names = sites)
  messages <- warnings <- character()

  if (!is.null(srh) && nrow(srh) > 0 && "source_region" %in% names(srh)) {
    out$prob_weighted_source_diversity <- NA_real_
    out$dominant_source_prop <- NA_real_
    ran$dominant_source_region <- NA_character_
    for (site in sites) {
      reg <- sr$region[match(site, sr$site)]
      spp <- names(which(comm[site, ] > 0))
      tab <- srh[srh$species %in% spp & srh$region == reg & !is.na(srh$source_region) & srh$source_region != "", , drop = FALSE]
      if (nrow(tab) == 0) next
      w <- if (probability_weighted && "entry_prob" %in% names(tab)) .std_prob(tab$entry_prob, 1) else rep(1, nrow(tab))
      s <- tapply(w, tab$source_region, sum, na.rm = TRUE)
      p <- s / sum(s, na.rm = TRUE)
      p <- p[p > 0 & !is.na(p)]
      out[site, "prob_weighted_source_diversity"] <- -sum(p * (log(p) / log(shannon_base)))
      out[site, "dominant_source_prop"] <- max(p)
      ran[site, "dominant_source_region"] <- names(p)[which.max(p)]
    }
  } else {
    messages <- c(messages, "species_region_history with source_region not available; source diversity skipped.")
  }

  # Region-level inflow/outflow from history_events if available.
  if (!is.null(hev) && nrow(hev) > 0) {
    hev$event_prob <- if ("event_prob" %in% names(hev) && probability_weighted) .std_prob(hev$event_prob, 1) else rep(1, nrow(hev))
    regions <- unique(sr$region)
    bal <- data.frame(region = regions, historical_inflow_outflow_balance = NA_real_, stringsAsFactors = FALSE)
    for (i in seq_along(regions)) {
      reg <- regions[i]
      inflow <- sum(hev$event_prob[hev$event_type == "dispersal" & hev$to_region == reg], na.rm = TRUE)
      outflow <- sum(hev$event_prob[hev$event_type == "dispersal" & hev$from_region == reg], na.rm = TRUE)
      denom <- inflow + outflow
      bal$historical_inflow_outflow_balance[i] <- if (denom > 0) (inflow - outflow) / denom else NA_real_
    }
    out$historical_inflow_outflow_balance <- bal$historical_inflow_outflow_balance[match(sr$region, bal$region)]
  } else {
    messages <- c(messages, "history_events not available; historical_inflow_outflow_balance skipped.")
  }

  out <- .drop_all_na_cols(out)
  ran <- ran[, vapply(ran, function(z) any(!is.na(z)), logical(1)), drop = FALSE]
  list(XData = out, random = ran, diagnostics = list(messages = unique(messages), warnings = unique(warnings)))
}
