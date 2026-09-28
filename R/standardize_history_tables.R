#' Standardize internal history tables
#'
#' Standardizes already stored history tables in A `HmscEcoEvo_data` object.
#' Missing optional probability columns are filled with 1. This function does
#' not infer historical processes.
#'
#' @param proj A `HmscEcoEvo_data` object.
#' @param time_unit Character; default "Ma".
#' @return Updated `HmscEcoEvo_data` object.
#' @export
standardize_history_tables <- function(proj, time_unit = "Ma") {
  if (!inherits(proj, "hmscHist_data")) stop("proj must be an hmscHist_data object.", call. = FALSE)
  tabs <- proj$history_tables
  messages <- proj$history_diagnostics$messages %||% character()
  warnings <- proj$history_diagnostics$warnings %||% character()
  species <- colnames(proj$comm)
  regions <- unique(proj$site_region$region)

  if (!is.null(tabs$species_region_history)) {
    res <- .std_species_region_history(tabs$species_region_history, species, regions)
    tabs$species_region_history <- res$table
    messages <- c(messages, res$messages)
    warnings <- c(warnings, res$warnings)
  }

  if (!is.null(tabs$speciation_events)) {
    res <- .std_speciation_events(tabs$speciation_events, species, regions)
    tabs$speciation_events <- res$table
    messages <- c(messages, res$messages)
    warnings <- c(warnings, res$warnings)
  }

  if (!is.null(tabs$history_events)) {
    res <- .std_history_events(tabs$history_events, species, regions)
    tabs$history_events <- res$table
    messages <- c(messages, res$messages)
    warnings <- c(warnings, res$warnings)
  }

  if (!is.null(tabs$trait_history)) {
    res <- .std_trait_history(tabs$trait_history, species)
    tabs$trait_history <- res$table
    messages <- c(messages, res$messages)
    warnings <- c(warnings, res$warnings)
  }

  proj$history_tables <- tabs
  proj$available_indices <- .available_indices_from_tables(tabs)
  proj$history_diagnostics <- list(
    messages = unique(messages),
    warnings = unique(warnings),
    time_unit = time_unit,
    available_indices = proj$available_indices
  )
  proj
}

.std_species_region_history <- function(x, species, regions) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  .require_cols(x, c("species", "region"), "species_region_history")
  msg <- warn <- character()

  x$species <- as.character(x$species)
  x$region <- as.character(x$region)
  x <- .add_if_missing(x, "entry_time", NA_real_)
  x <- .add_if_missing(x, "entry_time_dispersion", NA_real_)
  x <- .add_if_missing(x, "source_region", NA_character_)
  x <- .add_if_missing(x, "entry_prob", 1)
  x <- .add_if_missing(x, "retained", NA_real_)
  x <- .add_if_missing(x, "loss_time", NA_real_)

  x$entry_time <- suppressWarnings(as.numeric(x$entry_time))
  x$entry_time_dispersion <- suppressWarnings(as.numeric(x$entry_time_dispersion))
  x$entry_prob <- .std_prob(x$entry_prob, 1)
  x$retained <- suppressWarnings(as.numeric(x$retained))
  x$loss_time <- suppressWarnings(as.numeric(x$loss_time))

  missing_sp <- setdiff(species, unique(x$species))
  extra_sp <- setdiff(unique(x$species), species)
  if (length(missing_sp) > 0) msg <- c(msg, paste0(length(missing_sp), " species in comm are missing from species_region_history."))
  if (length(extra_sp) > 0) warn <- c(warn, paste0(length(extra_sp), " species in species_region_history are not in comm and were dropped."))
  x <- x[x$species %in% species, , drop = FALSE]

  bad_reg <- setdiff(unique(na.omit(c(x$region, x$source_region))), regions)
  if (length(bad_reg) > 0) warn <- c(warn, paste("Unknown region(s) in species_region_history:", paste(bad_reg, collapse = ", ")))

  x <- x[x$region %in% regions, , drop = FALSE]
  x <- x[order(x$species, x$region), , drop = FALSE]
  rownames(x) <- NULL
  list(table = x, messages = msg, warnings = warn)
}

.std_speciation_events <- function(x, species, regions) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  .require_cols(x, c("event_id", "region", "descendant_species"), "speciation_events")
  msg <- warn <- character()

  x$event_id <- as.character(x$event_id)
  x$region <- as.character(x$region)
  x$descendant_species <- as.character(x$descendant_species)
  x <- .add_if_missing(x, "time", NA_real_)
  x <- .add_if_missing(x, "event_prob", 1)
  x <- .add_if_missing(x, "event_type", "insitu")
  x$time <- suppressWarnings(as.numeric(x$time))
  x$event_prob <- .std_prob(x$event_prob, 1)

  extra_sp <- setdiff(unique(x$descendant_species), species)
  if (length(extra_sp) > 0) warn <- c(warn, paste0(length(extra_sp), " descendant species not in comm were dropped from speciation_events."))
  x <- x[x$descendant_species %in% species, , drop = FALSE]

  bad_reg <- setdiff(unique(na.omit(x$region)), regions)
  if (length(bad_reg) > 0) warn <- c(warn, paste("Unknown region(s) in speciation_events:", paste(bad_reg, collapse = ", ")))
  x <- x[x$region %in% regions, , drop = FALSE]
  rownames(x) <- NULL
  list(table = x, messages = msg, warnings = warn)
}

.std_history_events <- function(x, species, regions) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  .require_cols(x, "event_type", "history_events")
  msg <- warn <- character()

  x <- .add_if_missing(x, "event_id", paste0("event_", seq_len(nrow(x))))
  x <- .add_if_missing(x, "species", NA_character_)
  x <- .add_if_missing(x, "from_region", NA_character_)
  x <- .add_if_missing(x, "to_region", NA_character_)
  x <- .add_if_missing(x, "region", NA_character_)
  x <- .add_if_missing(x, "time", NA_real_)
  x <- .add_if_missing(x, "event_prob", 1)

  x$event_type <- .map_event_type(x$event_type)
  x$event_id <- as.character(x$event_id)
  x$species <- as.character(x$species)
  x$from_region <- as.character(x$from_region)
  x$to_region <- as.character(x$to_region)
  x$region <- as.character(x$region)
  x$time <- suppressWarnings(as.numeric(x$time))
  x$event_prob <- .std_prob(x$event_prob, 1)

  allowed <- c("dispersal", "loss", "founder")
  bad_type <- setdiff(unique(x$event_type), allowed)
  if (length(bad_type) > 0) warn <- c(warn, paste("Unsupported event_type(s) dropped:", paste(bad_type, collapse = ", ")))
  x <- x[x$event_type %in% allowed, , drop = FALSE]
  x$event_type[x$event_type == "founder"] <- "dispersal"

  has_sp <- !is.na(x$species) & x$species != "NA" & x$species != ""
  extra_sp <- setdiff(unique(x$species[has_sp]), species)
  if (length(extra_sp) > 0) warn <- c(warn, paste0(length(extra_sp), " species in history_events are not in comm."))

  regs <- unique(na.omit(c(x$from_region, x$to_region, x$region)))
  regs <- regs[regs != "NA" & regs != ""]
  bad_reg <- setdiff(regs, regions)
  if (length(bad_reg) > 0) warn <- c(warn, paste("Unknown region(s) in history_events:", paste(bad_reg, collapse = ", ")))

  ok_disp <- x$event_type == "dispersal" & x$from_region %in% regions & x$to_region %in% regions
  ok_loss <- x$event_type == "loss" & x$region %in% regions
  dropped <- sum(!(ok_disp | ok_loss))
  if (dropped > 0) msg <- c(msg, paste0(dropped, " history_events rows lacked enough region information and were dropped."))
  x <- x[ok_disp | ok_loss, , drop = FALSE]
  rownames(x) <- NULL
  list(table = x, messages = msg, warnings = warn)
}

.std_trait_history <- function(x, species) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!"species" %in% names(x)) {
    if (is.null(rownames(x))) stop("trait_history must have a species column or species rownames.", call. = FALSE)
    x$species <- rownames(x)
  }
  msg <- warn <- character()
  x$species <- as.character(x$species)
  extra_sp <- setdiff(unique(x$species), species)
  missing_sp <- setdiff(species, unique(x$species))
  if (length(extra_sp) > 0) warn <- c(warn, paste0(length(extra_sp), " species in trait_history are not in comm and were dropped."))
  if (length(missing_sp) > 0) msg <- c(msg, paste0(length(missing_sp), " species in comm are missing from trait_history."))
  x <- x[x$species %in% species, , drop = FALSE]
  rownames(x) <- x$species
  list(table = x, messages = msg, warnings = warn)
}

.available_indices_from_tables <- function(tabs) {
  out <- character()
  srh <- tabs$species_region_history
  if (!is.null(srh)) {
    if ("entry_time" %in% names(srh) && any(!is.na(srh$entry_time))) {
      out <- c(out, "colonization_age", "colonization_age_dispersion")
    }
    if ("entry_prob" %in% names(srh) && any(!is.na(srh$entry_prob))) {
      out <- c(out, "historical_gain_fraction", "historical_turnover_predictor")
    }
    if ("source_region" %in% names(srh) && any(!is.na(srh$source_region))) {
      out <- c(out, "prob_weighted_source_diversity", "dominant_source_prop", "dominant_source_region")
    }
    if ("retained" %in% names(srh) && any(!is.na(srh$retained))) {
      out <- c(out, "lineage_retention_index")
    }
    if (("loss_time" %in% names(srh) && any(!is.na(srh$loss_time))) ||
        ("retained" %in% names(srh) && any(!is.na(srh$retained) & srh$retained == 0))) {
      out <- c(out, "range_loss_fraction", "historical_gain_fraction", "historical_loss_fraction",
               "historical_turnover_predictor", "historical_gain_loss_balance",
               "historical_occupancy_balance")
    }
  }
  if (!is.null(tabs$speciation_events)) {
    out <- c(out, "insitu_speciation_prop", "insitu_diversification_strength", "insitu_exsitu_balance")
  }
  if (!is.null(tabs$history_events)) {
    out <- c(out, "historical_inflow_outflow_balance", "range_loss_fraction",
             "historical_gain_fraction", "historical_loss_fraction",
             "historical_turnover_predictor", "historical_gain_loss_balance",
             "historical_occupancy_balance")
  }
  if (!is.null(tabs$trait_history)) {
    out <- c(out, intersect(c("trait_conservatism", "trait_transition_count", "trait_stasis_time",
                              "trait_shift_magnitude", "dispersal_trait"), names(tabs$trait_history)))
  }
  unique(out)
}
