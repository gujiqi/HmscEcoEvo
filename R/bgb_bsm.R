#' Summarise BioGeoBEARS stochastic-map events into accessibility modifiers
#'
#' BioGeoBEARS stochastic maps are not cell-level suitability maps. They are
#' historical-biogeographic event records. This helper converts an event table
#' into transparent region/time summaries that can be reported or used as
#' optional modifiers of `A_BGB`.
#'
#' @param events BSM/event table. Common columns are `time_ma`, `event_type`,
#'   `from_region`, `to_region`, `region`, `species`, and `event_prob`.
#' @param time_bin Optional bin width in Ma. If supplied, events are grouped by
#'   `floor(time_ma / time_bin) * time_bin`.
#'
#' @return A list with `region_metrics`, `dispersal_rates`,
#'   `source_sink_long`, and `events`.
#' @export
#'
#' @examples
#' ev <- data.frame(time_ma = c(50, 50), event_type = "dispersal",
#'   from_region = c("A", "A"), to_region = c("B", "C"), event_prob = 1)
#' hee_bgb_bsm_metrics(ev)$region_metrics
hee_bgb_bsm_metrics <- function(events, time_bin = NULL) {
  ev <- hee_bgb_bsm(events)
  if (nrow(ev) == 0L) {
    empty <- data.frame()
    return(list(region_metrics = empty, dispersal_rates = empty,
                source_sink_long = empty, events = ev))
  }
  if (!"event_prob" %in% names(ev)) ev$event_prob <- 1
  if (!"time_ma" %in% names(ev)) ev$time_ma <- 0
  if (!is.null(time_bin)) ev$time_ma <- floor(ev$time_ma / time_bin) * time_bin
  ev$event_type_norm <- tolower(gsub("[ _-]+", "_", ev$event_type %||% "event"))

  disp_types <- c("dispersal", "range_expansion", "founder_event", "founder-event")
  ext_types <- c("extinction", "local_extinction", "range_contraction", "extirpation")
  d <- ev[ev$event_type_norm %in% disp_types, , drop = FALSE]
  e <- ev[ev$event_type_norm %in% ext_types, , drop = FALSE]

  source <- data.frame()
  sink <- data.frame()
  dispersal <- data.frame()
  if (nrow(d) > 0L && all(c("from_region", "to_region") %in% names(d))) {
    source <- stats::aggregate(event_prob ~ from_region + time_ma, d, sum, na.rm = TRUE)
    names(source) <- c("region", "time_ma", "region_source_strength")
    sink <- stats::aggregate(event_prob ~ to_region + time_ma, d, sum, na.rm = TRUE)
    names(sink) <- c("region", "time_ma", "region_sink_strength")
    dispersal <- stats::aggregate(event_prob ~ from_region + to_region + time_ma, d, sum, na.rm = TRUE)
    names(dispersal)[4] <- "dispersal_rate"
  }

  extinction <- data.frame()
  if (nrow(e) > 0L) {
    if (!"region" %in% names(e) && "from_region" %in% names(e)) e$region <- e$from_region
    if ("region" %in% names(e)) {
      extinction <- stats::aggregate(event_prob ~ region + time_ma, e, sum, na.rm = TRUE)
      names(extinction)[3] <- "extinction_rate"
    }
  }

  regions <- unique(c(source$region, sink$region, extinction$region))
  times <- sort(unique(ev$time_ma), decreasing = TRUE)
  grid <- expand.grid(region = regions, time_ma = times, stringsAsFactors = FALSE)
  region_metrics <- Reduce(function(a, b) merge(a, b, by = c("region", "time_ma"), all.x = TRUE, sort = FALSE),
                           Filter(nrow, list(grid, source, sink, extinction)))
  if (nrow(region_metrics) == 0L) {
    region_metrics <- grid
  }
  for (nm in c("region_source_strength", "region_sink_strength", "extinction_rate")) {
    if (!nm %in% names(region_metrics)) region_metrics[[nm]] <- 0
    region_metrics[[nm]][is.na(region_metrics[[nm]])] <- 0
  }
  region_metrics$historical_turnover <- with(
    region_metrics,
    region_source_strength + region_sink_strength + extinction_rate
  )
  region_metrics <- region_metrics[order(region_metrics$region, -region_metrics$time_ma), , drop = FALSE]
  region_metrics$lineage_accumulation <- stats::ave(
    region_metrics$region_sink_strength - region_metrics$extinction_rate,
    region_metrics$region,
    FUN = cumsum
  )
  source_sink_long <- rbind(
    data.frame(region = region_metrics$region, time_ma = region_metrics$time_ma,
               metric = "region_source_strength",
               value = region_metrics$region_source_strength, stringsAsFactors = FALSE),
    data.frame(region = region_metrics$region, time_ma = region_metrics$time_ma,
               metric = "region_sink_strength",
               value = region_metrics$region_sink_strength, stringsAsFactors = FALSE),
    data.frame(region = region_metrics$region, time_ma = region_metrics$time_ma,
               metric = "extinction_rate",
               value = region_metrics$extinction_rate, stringsAsFactors = FALSE),
    data.frame(region = region_metrics$region, time_ma = region_metrics$time_ma,
               metric = "historical_turnover",
               value = region_metrics$historical_turnover, stringsAsFactors = FALSE)
  )
  list(region_metrics = region_metrics,
       dispersal_rates = dispersal,
       source_sink_long = source_sink_long,
       events = ev)
}

#' Adjust historical accessibility using BSM-derived region metrics
#'
#' This is optional. It makes the adjustment explicit instead of hiding it in
#' the final projection: sinks and sources can increase accessibility, whereas
#' extinction and high turnover can reduce it.
#'
#' @param accessibility Table with `species`, `region`, `time_ma`,
#'   and `accessibility`.
#' @param bsm_metrics Output from `hee_bgb_bsm_metrics()` or its
#'   `region_metrics` table.
#' @param sink_bonus,source_bonus Multipliers for sink/source strength.
#' @param extinction_penalty,turnover_penalty Multipliers for extinction and
#'   turnover penalties.
#'
#' @return Accessibility table with adjusted values and modifier columns.
#' @export
hee_adjust_accessibility_bsm <- function(accessibility,
                                         bsm_metrics,
                                         sink_bonus = 0.10,
                                         source_bonus = 0.05,
                                         extinction_penalty = 0.10,
                                         turnover_penalty = 0) {
  acc <- as.data.frame(accessibility)
  .require_cols(acc, c("region", "time_ma", "accessibility"), "accessibility")
  metrics <- if (is.list(bsm_metrics) && "region_metrics" %in% names(bsm_metrics)) {
    bsm_metrics$region_metrics
  } else {
    as.data.frame(bsm_metrics)
  }
  if (nrow(metrics) == 0L) {
    acc$bsm_modifier <- 1
    acc$accessibility_bsm_adjusted <- acc$accessibility
    return(acc)
  }
  acc <- merge(acc, metrics, by = c("region", "time_ma"), all.x = TRUE, sort = FALSE)
  for (nm in c("region_source_strength", "region_sink_strength", "extinction_rate", "historical_turnover")) {
    if (!nm %in% names(acc)) acc[[nm]] <- 0
    acc[[nm]][is.na(acc[[nm]])] <- 0
  }
  bonus <- 1 + source_bonus * acc$region_source_strength + sink_bonus * acc$region_sink_strength
  penalty <- exp(-extinction_penalty * acc$extinction_rate - turnover_penalty * acc$historical_turnover)
  acc$bsm_modifier <- pmin(pmax(bonus * penalty, 0), 2)
  acc$accessibility_bsm_adjusted <- pmin(pmax(acc$accessibility * acc$bsm_modifier, 0), 1)
  acc
}
