#' Calculate historical occupancy dynamics indices
#'
#' Calculates colonization, retention, gain/loss fractions, turnover predictors,
#' and occupancy balance predictors. A `species_region_history` table is
#' sufficient for the basic indices. A `history_events` table adds event-based
#' gain/loss fractions.
#'
#' @param proj A `HmscEcoEvo_data` object.
#' @param probability_weighted Logical; use probability columns when available.
#' @return A list with `XData` and `diagnostics`.
#' @export
calc_occupancy_dynamics_indices <- function(proj, probability_weighted = TRUE) {
  comm <- proj$comm
  sr <- proj$site_region
  tabs <- proj$history_tables
  sites <- rownames(comm)
  out <- data.frame(row.names = sites)
  messages <- warnings <- character()

  srh <- tabs$species_region_history
  if (!is.null(srh) && nrow(srh) > 0) {
    out$colonization_age <- NA_real_
    out$colonization_age_dispersion <- NA_real_
    out$lineage_retention_index <- NA_real_
    for (site in sites) {
      reg <- sr$region[match(site, sr$site)]
      spp <- names(which(comm[site, ] > 0))
      tab <- srh[srh$species %in% spp & srh$region == reg, , drop = FALSE]
      if (nrow(tab) == 0) next
      w <- if (probability_weighted && "entry_prob" %in% names(tab)) .std_prob(tab$entry_prob, 1) else rep(1, nrow(tab))
      if ("entry_time" %in% names(tab)) {
        out[site, "colonization_age"] <- .weighted_mean_safe(tab$entry_time, w)
        out[site, "colonization_age_dispersion"] <- .weighted_sd_safe(tab$entry_time, w)
      }
      if ("retained" %in% names(tab)) {
        rr <- suppressWarnings(as.numeric(tab$retained))
        out[site, "lineage_retention_index"] <- .weighted_mean_safe(rr, w)
      }
    }
  } else {
    messages <- c(messages, "species_region_history not available; colonization and retention indices skipped.")
  }

  region_fractions <- .regional_gain_loss_fractions(tabs$history_events, srh, unique(sr$region), probability_weighted)
  if (!is.null(region_fractions) && nrow(region_fractions) > 0) {
    m <- match(sr$region, region_fractions$region)
    out$historical_gain_fraction <- region_fractions$historical_gain_fraction[m]
    out$historical_loss_fraction <- region_fractions$historical_loss_fraction[m]
    out$range_loss_fraction <- region_fractions$range_loss_fraction[m]
    out$historical_turnover_predictor <- region_fractions$historical_turnover_predictor[m]
    out$historical_gain_loss_balance <- region_fractions$historical_gain_loss_balance[m]
  } else {
    messages <- c(messages, "No usable history_events or loss information; gain/loss fractions and turnover predictors skipped.")
  }

  # Composite predictor. Use only if all components exist.
  need <- c("historical_gain_fraction", "lineage_retention_index", "historical_loss_fraction")
  if (all(need %in% names(out))) {
    out$historical_occupancy_balance <- as.numeric(.z(out$historical_gain_fraction) +
                                                     .z(out$lineage_retention_index) -
                                                     .z(out$historical_loss_fraction))
  } else {
    messages <- c(messages, "historical_occupancy_balance skipped because gain fraction, retention, or loss fraction component is missing.")
  }

  out <- .drop_all_na_cols(out)
  list(XData = out, diagnostics = list(messages = unique(messages), warnings = unique(warnings)))
}

.regional_gain_loss_fractions <- function(history_events, species_region_history, regions, probability_weighted = TRUE) {
  rows <- list()

  if (!is.null(history_events) && nrow(history_events) > 0) {
    ev <- history_events
    ev$event_prob <- if ("event_prob" %in% names(ev) && probability_weighted) .std_prob(ev$event_prob, 1) else rep(1, nrow(ev))
    gain_by_region <- stats::setNames(rep(0, length(regions)), regions)
    loss_by_region <- stats::setNames(rep(0, length(regions)), regions)
    gains <- tapply(ev$event_prob[ev$event_type == "dispersal"], ev$to_region[ev$event_type == "dispersal"], sum, na.rm = TRUE)
    losses <- tapply(ev$event_prob[ev$event_type == "loss"], ev$region[ev$event_type == "loss"], sum, na.rm = TRUE)
    gain_by_region[names(gains)] <- gains
    loss_by_region[names(losses)] <- losses
    total_gain <- sum(gain_by_region, na.rm = TRUE)
    total_loss <- sum(loss_by_region, na.rm = TRUE)
    total_turnover <- total_gain + total_loss
    for (reg in regions) {
      gain <- gain_by_region[[reg]]
      loss <- loss_by_region[[reg]]
      local_total <- gain + loss
      rows[[length(rows) + 1]] <- data.frame(
        region = reg,
        historical_gain_fraction = if (total_gain > 0) gain / total_gain else NA_real_,
        historical_loss_fraction = if (total_loss > 0) loss / total_loss else NA_real_,
        range_loss_fraction = if (local_total > 0) loss / local_total else NA_real_,
        historical_turnover_predictor = if (total_turnover > 0) local_total / total_turnover else NA_real_,
        historical_gain_loss_balance = if (local_total > 0) (gain - loss) / local_total else NA_real_,
        stringsAsFactors = FALSE
      )
    }
    return(.safe_bind_rows(rows))
  }

  # Fallback from species_region_history: use loss_time or retained == 0.
  if (!is.null(species_region_history) && nrow(species_region_history) > 0) {
    srh <- species_region_history
    srh$entry_prob <- if ("entry_prob" %in% names(srh) && probability_weighted) .std_prob(srh$entry_prob, 1) else rep(1, nrow(srh))
    gain_by_region <- stats::setNames(rep(0, length(regions)), regions)
    loss_by_region <- stats::setNames(rep(0, length(regions)), regions)
    for (reg in regions) {
      tab <- srh[srh$region == reg, , drop = FALSE]
      if (nrow(tab) == 0) next
      gain_by_region[[reg]] <- sum(tab$entry_prob, na.rm = TRUE)
      loss_flag <- rep(FALSE, nrow(tab))
      if ("loss_time" %in% names(tab)) loss_flag <- loss_flag | !is.na(tab$loss_time)
      if ("retained" %in% names(tab)) loss_flag <- loss_flag | (!is.na(tab$retained) & tab$retained == 0)
      loss_by_region[[reg]] <- sum(tab$entry_prob[loss_flag], na.rm = TRUE)
    }
    total_gain <- sum(gain_by_region, na.rm = TRUE)
    total_loss <- sum(loss_by_region, na.rm = TRUE)
    total_turnover <- total_gain + total_loss
    for (reg in regions) {
      gain <- gain_by_region[[reg]]
      loss <- loss_by_region[[reg]]
      local_total <- gain + loss
      rows[[length(rows) + 1]] <- data.frame(
        region = reg,
        historical_gain_fraction = if (total_gain > 0) gain / total_gain else NA_real_,
        historical_loss_fraction = if (total_loss > 0) loss / total_loss else NA_real_,
        range_loss_fraction = if (local_total > 0) loss / local_total else NA_real_,
        historical_turnover_predictor = if (total_turnover > 0) local_total / total_turnover else NA_real_,
        historical_gain_loss_balance = if (local_total > 0) (gain - loss) / local_total else NA_real_,
        stringsAsFactors = FALSE
      )
    }
    return(.safe_bind_rows(rows))
  }
  NULL
}
