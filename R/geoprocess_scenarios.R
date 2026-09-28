#' Run an end-to-end geological-process scenario suite
#'
#' Builds compact deterministic scenarios for the Dynamic Earth-Biota Assembly
#' workflow described in the geological-process framework: arena creation and
#' loss, land age, area and heterogeneity, structural-functional-climatic
#' connectivity, isolation duration, source pressure, rescue, colonisation,
#' local extinction, regional species-pool change, speciation opportunity,
#' dynamic occupancy, richness, turnover, refugia, dominant limitation, and
#' prediction reliability.
#'
#' The suite is designed for examples, regression tests, and report generation.
#' It deliberately includes normal and pathological cases. Pathological cases
#' such as duplicate keys are expected to stop with a clear validation error;
#' the suite records those errors instead of silently repairing them.
#'
#' @param scenarios Character vector of scenario names. Supported values are
#'   `"normal"`, `"extreme"`, `"missing_optional"`, `"duplicate_key"`,
#'   `"single_time"`, `"uneven_time"`, `"shuffled_time"`, and
#'   `"land_appearance_loss"`.
#' @param species Character vector of species names.
#' @param seed Optional seed used only for deterministic row shuffling.
#' @param include_expected_errors Logical; include scenarios that are expected
#'   to error, such as `"duplicate_key"`.
#' @return A `hee_geoprocess_scenario_suite` list with `scenarios`,
#'   `input_tables`, `result_tables`, `validation_table`,
#'   `scenario_summary`, `formula_catalog`, and `interpretation_table`.
#' @export
#'
#' @examples
#' suite <- hee_geoprocess_scenario_suite(
#'   scenarios = c("normal", "single_time"),
#'   species = c("sp1", "sp2")
#' )
#' suite$scenario_summary
hee_geoprocess_scenario_suite <- function(
    scenarios = c("normal", "extreme", "missing_optional", "duplicate_key",
                  "single_time", "uneven_time", "shuffled_time",
                  "land_appearance_loss"),
    species = c("sp1", "sp2", "sp3"),
    seed = 1,
    include_expected_errors = TRUE) {
  supported <- c("normal", "extreme", "missing_optional", "duplicate_key",
                 "single_time", "uneven_time", "shuffled_time",
                 "land_appearance_loss")
  scenarios <- unique(as.character(scenarios))
  bad <- setdiff(scenarios, supported)
  if (length(bad) > 0L) {
    stop("Unknown scenario(s): ", paste(bad, collapse = ", "), call. = FALSE)
  }
  if (!include_expected_errors) {
    scenarios <- setdiff(scenarios, "duplicate_key")
  }
  if (length(species) < 1L || anyNA(species) || any(!nzchar(species))) {
    stop("species must contain at least one non-empty species name.",
         call. = FALSE)
  }
  species <- unique(as.character(species))

  if (!is.null(seed)) set.seed(seed)
  scenario_results <- vector("list", length(scenarios))
  names(scenario_results) <- scenarios
  for (scenario in scenarios) {
    scenario_results[[scenario]] <- .hee_run_one_geoprocess_scenario(
      scenario = scenario,
      species = species
    )
  }

  collect_inputs <- function(name) {
    rows <- lapply(scenario_results, function(z) {
      if (is.null(z$inputs) || is.null(z$inputs[[name]])) return(NULL)
      d <- as.data.frame(z$inputs[[name]])
      d$scenario <- if (nrow(d) == 0L) character() else z$scenario
      d
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows) == 0L) data.frame() else .hee_bind_fill(rows)
  }
  collect_results <- function(name) {
    rows <- lapply(scenario_results, function(z) {
      if (is.null(z$results) || is.null(z$results[[name]])) return(NULL)
      d <- as.data.frame(z$results[[name]])
      d$scenario <- if (nrow(d) == 0L) character() else z$scenario
      d
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows) == 0L) data.frame() else .hee_bind_fill(rows)
  }

  input_names <- unique(unlist(lapply(scenario_results, function(z) {
    names(z$inputs %||% list())
  })))
  result_names <- unique(unlist(lapply(scenario_results, function(z) {
    names(z$results %||% list())
  })))
  input_tables <- stats::setNames(lapply(input_names, collect_inputs),
                                  input_names)
  result_tables <- stats::setNames(lapply(result_names, collect_results),
                                   result_names)
  validation_table <- .hee_bind_fill(lapply(scenario_results, `[[`,
                                           "validation"))
  scenario_summary <- .hee_bind_fill(lapply(scenario_results, function(z) {
    data.frame(
      scenario = z$scenario,
      status = z$status,
      expected_error = z$expected_error,
      error_message = z$error_message %||% "",
      n_time_slices = z$n_time_slices,
      n_species = length(species),
      n_validation_failures = sum(z$validation$status == "FAIL", na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  out <- list(
    scenarios = scenario_results,
    input_tables = input_tables,
    result_tables = result_tables,
    validation_table = validation_table,
    scenario_summary = scenario_summary,
    formula_catalog = .hee_geoprocess_formula_catalog(),
    interpretation_table = .hee_geoprocess_interpretation_table()
  )
  class(out) <- c("hee_geoprocess_scenario_suite", class(out))
  out
}

#' Plot a geological-process scenario suite
#'
#' @param suite Output from `hee_geoprocess_scenario_suite()`.
#' @return A named list of `ggplot` objects.
#' @export
#'
#' @examples
#' \donttest{
#' suite <- hee_geoprocess_scenario_suite(c("normal", "single_time"))
#' plot_geoprocess_scenario_suite(suite)
#' }
plot_geoprocess_scenario_suite <- function(suite) {
  if (!inherits(suite, "hee_geoprocess_scenario_suite")) {
    stop("suite must be from hee_geoprocess_scenario_suite().", call. = FALSE)
  }
  .require_pkg("ggplot2")
  ggplot2 <- asNamespace("ggplot2")

  plots <- list()
  ps <- suite$result_tables$process_summary
  if (!is.null(ps) && nrow(ps) > 0L) {
    process_cols <- intersect(
      c("geographic_existence", "ecological_opportunity", "connectivity",
        "isolation", "source_pressure", "rescue_effect",
        "colonisation_probability", "extinction_probability",
        "speciation_opportunity", "occupancy_probability",
        "prediction_reliability"),
      names(ps)
    )
    if (length(process_cols) > 0L) {
      means <- stats::aggregate(ps[, process_cols, drop = FALSE],
                                ps[, "scenario", drop = FALSE],
                                mean, na.rm = TRUE)
      long <- stats::reshape(means, varying = process_cols,
                             v.names = "mean_value", timevar = "process",
                             times = process_cols, direction = "long")
      plots$process_heatmap <- ggplot2$ggplot(
        long,
        ggplot2$aes(x = .data$process, y = .data$scenario,
                    fill = .data$mean_value)
      ) +
        ggplot2$geom_tile(color = "white", linewidth = 0.2) +
        ggplot2$scale_fill_gradient(low = "white", high = "#225EA8",
                                    na.value = "grey90") +
        ggplot2$labs(title = "Geological-process scenario means",
                     x = "Process", y = "Scenario", fill = "Mean") +
        ggplot2$theme_bw(base_size = 9) +
        ggplot2$theme(axis.text.x = ggplot2$element_text(
          angle = 35, hjust = 1
        ))
    }
  }

  val <- suite$validation_table
  if (!is.null(val) && nrow(val) > 0L) {
    vv <- as.data.frame(table(val$scenario, val$status),
                        stringsAsFactors = FALSE)
    names(vv) <- c("scenario", "status", "n")
    plots$validation_status <- ggplot2$ggplot(
      vv,
      ggplot2$aes(x = .data$scenario, y = .data$n, fill = .data$status)
    ) +
      ggplot2$geom_col(position = "stack") +
      ggplot2$coord_flip() +
      ggplot2$scale_fill_manual(values = c(PASS = "#238B45",
                                           FAIL = "#CB181D",
                                           INFO = "#6BAED6"),
                                na.value = "grey80") +
      ggplot2$labs(title = "Automatic validation status",
                   x = NULL, y = "Check count", fill = "Status") +
      ggplot2$theme_bw(base_size = 9)
  }

  events <- suite$result_tables$landscape_events
  if (!is.null(events) && nrow(events) > 0L) {
    ee <- as.data.frame(table(events$scenario, events$event_type),
                        stringsAsFactors = FALSE)
    names(ee) <- c("scenario", "event_type", "n")
    plots$landscape_events <- ggplot2$ggplot(
      ee,
      ggplot2$aes(x = .data$event_type, y = .data$n,
                  fill = .data$scenario)
    ) +
      ggplot2$geom_col(position = "dodge") +
      ggplot2$labs(title = "Arena creation, loss, and transformation events",
                   x = "Event type", y = "Count", fill = "Scenario") +
      ggplot2$theme_bw(base_size = 9) +
      ggplot2$theme(axis.text.x = ggplot2$element_text(
        angle = 35, hjust = 1
      ))
  }

  dyn <- suite$result_tables$dynamic_occupancy
  if (!is.null(dyn) && nrow(dyn) > 0L &&
      all(c("time_ma", "occupancy_probability", "scenario") %in% names(dyn))) {
    dt <- stats::aggregate(occupancy_probability ~ scenario + time_ma,
                           dyn, mean, na.rm = TRUE)
    dt$.n_time <- stats::ave(dt$time_ma, dt$scenario,
                             FUN = function(z) length(unique(z)))
    dt_line <- dt[dt$.n_time > 1L, , drop = FALSE]
    plots$dynamic_occupancy <- ggplot2$ggplot(
      dt,
      ggplot2$aes(x = .data$time_ma,
                  y = .data$occupancy_probability,
                  color = .data$scenario)
    ) +
      ggplot2$geom_line(data = dt_line, linewidth = 0.8) +
      ggplot2$geom_point(size = 1.3) +
      ggplot2$scale_x_reverse() +
      ggplot2$labs(title = "Dynamic occupancy through scenario time axes",
                   x = "Time (Ma)", y = "Mean occupancy probability",
                   color = "Scenario") +
      ggplot2$theme_bw(base_size = 9)
  }

  net <- suite$result_tables$network_summary
  if (!is.null(net) && nrow(net) > 0L &&
      all(c("time_ma", "fragmentation_index", "largest_component_fraction",
            "scenario") %in% names(net))) {
    nm <- stats::aggregate(
      cbind(fragmentation_index, largest_component_fraction) ~ scenario + time_ma,
      net, mean, na.rm = TRUE
    )
    nl <- stats::reshape(nm,
                         varying = c("fragmentation_index",
                                     "largest_component_fraction"),
                         v.names = "value", timevar = "network_metric",
                         times = c("fragmentation_index",
                                   "largest_component_fraction"),
                         direction = "long")
    nl$.n_time <- stats::ave(
      nl$time_ma,
      interaction(nl$scenario, nl$network_metric, drop = TRUE),
      FUN = function(z) length(unique(z))
    )
    nl_line <- nl[nl$.n_time > 1L, , drop = FALSE]
    plots$network_fragmentation <- ggplot2$ggplot(
      nl,
      ggplot2$aes(x = .data$time_ma, y = .data$value,
                  color = .data$network_metric)
    ) +
      ggplot2$geom_line(data = nl_line, linewidth = 0.8) +
      ggplot2$geom_point(size = 1.1) +
      ggplot2$facet_wrap(~ scenario) +
      ggplot2$scale_x_reverse() +
      ggplot2$labs(title = "Network connectivity and fragmentation",
                   x = "Time (Ma)", y = "Index", color = "Metric") +
      ggplot2$theme_bw(base_size = 9)
  }

  lin_ext <- suite$result_tables$lineage_extinction
  reg_ext <- suite$result_tables$regional_extinction
  loc_ext <- suite$result_tables$local_extinction
  ext_rows <- list()
  if (!is.null(loc_ext) && nrow(loc_ext) > 0L) {
    x <- stats::aggregate(local_extinction_risk ~ scenario + time_ma,
                          loc_ext, mean, na.rm = TRUE)
    names(x)[names(x) == "local_extinction_risk"] <- "value"
    x$scale <- "local"
    ext_rows[[length(ext_rows) + 1L]] <- x
  }
  if (!is.null(reg_ext) && nrow(reg_ext) > 0L) {
    x <- stats::aggregate(regional_extinction_proxy ~ scenario + time_ma,
                          reg_ext, mean, na.rm = TRUE)
    names(x)[names(x) == "regional_extinction_proxy"] <- "value"
    x$scale <- "regional"
    ext_rows[[length(ext_rows) + 1L]] <- x
  }
  if (!is.null(lin_ext) && nrow(lin_ext) > 0L) {
    x <- stats::aggregate(lineage_extinction_proxy ~ scenario + time_ma,
                          lin_ext, mean, na.rm = TRUE)
    names(x)[names(x) == "lineage_extinction_proxy"] <- "value"
    x$scale <- "lineage"
    ext_rows[[length(ext_rows) + 1L]] <- x
  }
  if (length(ext_rows) > 0L) {
    ex <- .hee_bind_fill(ext_rows)
    ex$.n_time <- stats::ave(
      ex$time_ma,
      interaction(ex$scenario, ex$scale, drop = TRUE),
      FUN = function(z) length(unique(z))
    )
    ex_line <- ex[ex$.n_time > 1L, , drop = FALSE]
    plots$extinction_layers <- ggplot2$ggplot(
      ex,
      ggplot2$aes(x = .data$time_ma, y = .data$value,
                  color = .data$scale)
    ) +
      ggplot2$geom_line(data = ex_line, linewidth = 0.8) +
      ggplot2$geom_point(size = 1.1) +
      ggplot2$facet_wrap(~ scenario) +
      ggplot2$scale_x_reverse() +
      ggplot2$scale_y_continuous(limits = c(0, 1)) +
      ggplot2$labs(title = "Local, regional, and lineage extinction proxies",
                   x = "Time (Ma)", y = "Mean proxy", color = "Scale") +
      ggplot2$theme_bw(base_size = 9)
  }

  cmg <- suite$result_tables$cradle_museum_grave
  if (!is.null(cmg) && nrow(cmg) > 0L &&
      all(c("scenario", "dominant_region_status") %in% names(cmg))) {
    cc <- as.data.frame(table(cmg$scenario, cmg$dominant_region_status),
                        stringsAsFactors = FALSE)
    names(cc) <- c("scenario", "dominant_region_status", "n")
    plots$cradle_museum_grave <- ggplot2$ggplot(
      cc,
      ggplot2$aes(x = .data$dominant_region_status, y = .data$n,
                  fill = .data$scenario)
    ) +
      ggplot2$geom_col(position = "dodge") +
      ggplot2$labs(title = "Cradle, museum, grave, source, and sink proxies",
                   x = "Dominant regional status", y = "Region-time count",
                   fill = "Scenario") +
      ggplot2$theme_bw(base_size = 9) +
      ggplot2$theme(axis.text.x = ggplot2$element_text(
        angle = 30, hjust = 1
      ))
  }

  class(plots) <- c("hee_plot_list", class(plots))
  plots
}

.hee_run_one_geoprocess_scenario <- function(scenario, species) {
  expected_error <- identical(scenario, "duplicate_key")
  inputs <- .hee_geoprocess_scenario_inputs(scenario, species)
  error_message <- NULL
  results <- NULL
  status <- "PASS"
  validation <- NULL

  run <- tryCatch({
    geo <- hee_geoprocess_diagnostics(
      landscape_state = inputs$landscape_state,
      projection = inputs$projection,
      connectivity_cube = inputs$connectivity_cube,
      accessibility = inputs$accessibility,
      phylo_mask = inputs$phylo_mask,
      extrapolation = inputs$extrapolation,
      tip_ranges = inputs$tip_ranges,
      species_origin = inputs$species_origin,
      bsm_events = inputs$bsm_events,
      traits = inputs$traits,
      times = inputs$times,
      env_cols = inputs$env_cols,
      land_col = "geographic_existence",
      area_col = "land_area_km2",
      initial_probability = 0
    )
    derived <- .hee_geoprocess_scenario_derived(geo, inputs)
    c(geo, derived)
  }, error = function(e) e)

  if (inherits(run, "error")) {
    error_message <- conditionMessage(run)
    status <- if (expected_error) "EXPECTED_ERROR" else "ERROR"
    validation <- data.frame(
      scenario = scenario,
      module = "input_validation",
      check = if (expected_error) "expected_clear_error" else "unexpected_error",
      status = if (expected_error) "PASS" else "FAIL",
      severity = if (expected_error) "info" else "critical",
      message = error_message,
      interpretation_type = if (expected_error) "automatic_validation" else
        "blocking_error",
      stringsAsFactors = FALSE
    )
  } else {
    results <- run
    validation <- .hee_validate_geoprocess_scenario(scenario, inputs, results)
    status <- if (any(validation$status == "FAIL")) "FAIL" else "PASS"
  }

  list(
    scenario = scenario,
    expected_error = expected_error,
    status = status,
    error_message = error_message,
    inputs = inputs,
    results = results,
    validation = validation,
    n_time_slices = length(inputs$times)
  )
}

.hee_geoprocess_scenario_inputs <- function(scenario, species) {
  times <- switch(
    scenario,
    single_time = 0,
    uneven_time = c(540, 517, 251, 66, 0),
    c(540, 520, 500, 250, 0)
  )
  regions <- c("R1", "R2", "R3")
  cells0 <- data.frame(
    cell_id = paste0("c", seq_len(6)),
    region = rep(regions, each = 2),
    lon = c(-40, -34, 0, 6, 42, 48),
    lat = c(-10, 2, -4, 8, -2, 10),
    cell_area_base = c(80, 95, 120, 110, 70, 85),
    stringsAsFactors = FALSE
  )
  landscape <- merge(
    expand.grid(cell_id = cells0$cell_id, time_ma = times,
                KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE),
    cells0,
    by = "cell_id", all.x = TRUE, sort = FALSE
  )
  age01 <- (max(times) - landscape$time_ma) / max(max(times), 1)
  region_shift <- match(landscape$region, regions) - 2
  landscape$geographic_existence <- 1
  landscape$habitat_availability <- .hee_clip01(0.85 - 0.12 * abs(region_shift) +
                                                  0.08 * sin(age01 * pi))
  if (scenario == "extreme") {
    landscape$habitat_availability <- ifelse(landscape$region == "R2",
                                             0.08, 0.95)
  }
  if (scenario == "land_appearance_loss") {
    landscape$geographic_existence[
      landscape$cell_id %in% c("c1", "c2") & landscape$time_ma >= 520
    ] <- 0
    landscape$geographic_existence[
      landscape$cell_id %in% c("c5", "c6") & landscape$time_ma == 500
    ] <- 0
  }
  if (scenario == "extreme") {
    landscape$geographic_existence[
      landscape$region == "R3" & landscape$time_ma %in% c(500, 250)
    ] <- 0
  }
  landscape$land <- landscape$geographic_existence
  landscape$land_mask_dem <- landscape$geographic_existence
  landscape$land_area_km2 <- ifelse(
    landscape$geographic_existence > 0,
    landscape$cell_area_base * (0.55 + 0.7 * age01) *
      pmax(0.15, landscape$habitat_availability),
    0
  )
  landscape$elev_m <- ifelse(landscape$geographic_existence > 0,
                             100 + 900 * age01 + 120 * region_shift, 0)
  landscape$bio1 <- 18 - 10 * age01 - 0.04 * landscape$elev_m + region_shift
  landscape$bio12 <- 800 + 350 * sin(age01 * 2 * pi + region_shift) -
    2 * landscape$lat
  landscape$habitat_heterogeneity <- .hee_clip01(
    0.25 + 0.45 * abs(sin(age01 * pi + match(landscape$cell_id, cells0$cell_id)))
  )
  landscape$paleo_habitat <- ifelse(landscape$habitat_availability > 0.6,
                                    "forest_wet", "open_dry")
  if (scenario == "duplicate_key") {
    landscape <- rbind(landscape, landscape[1, , drop = FALSE])
  }

  traits <- data.frame(
    species = species,
    dispersal_distance = seq(4, 18, length.out = length(species)),
    dispersal = seq(0.2, 1, length.out = length(species)),
    vulnerability = rev(seq(0.15, 0.9, length.out = length(species))),
    stringsAsFactors = FALSE
  )

  projection <- merge(
    expand.grid(species = species, cell_id = cells0$cell_id, time_ma = times,
                KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE),
    landscape[, c("cell_id", "time_ma", "region", "lon", "lat",
                  "geographic_existence", "land_mask_dem",
                  "habitat_availability", "land_area_km2", "bio1", "bio12",
                  "habitat_heterogeneity"), drop = FALSE],
    by = c("cell_id", "time_ma"), all.x = TRUE, sort = FALSE
  )
  sp_index <- match(projection$species, species)
  projection$suitability <- .hee_clip01(stats::pnorm(
    scale(projection$bio1)[, 1] * (0.5 - 0.12 * sp_index) +
      scale(projection$bio12)[, 1] * (-0.1 + 0.08 * sp_index) +
      projection$habitat_heterogeneity - 0.2
  ))
  projection$probability <- .hee_clip01(projection$suitability *
                                          projection$geographic_existence)
  projection$accessibility <- .hee_clip01(0.55 + 0.12 * sp_index -
                                            0.15 * (projection$region == "R3"))
  projection$D_static <- .hee_clip01(0.85 - 0.1 * (projection$region == "R3"))
  projection$D_dynamic <- .hee_clip01(0.6 + 0.25 * projection$probability)
  projection$phylo_existence <- 1
  projection$Q_extrap <- .hee_clip01(1 - 0.25 * abs(scale(projection$bio1)[, 1]))
  projection$extrapolation_score <- .hee_clip01(1 - projection$Q_extrap)
  if (scenario == "extreme") {
    projection$suitability[projection$region == "R2"] <- 0.02
    projection$probability <- .hee_clip01(projection$suitability *
                                            projection$geographic_existence)
    projection$extrapolation_score[projection$region == "R2"] <- 0.95
  }

  origins <- stats::setNames(seq(max(times), min(times), length.out = length(species)),
                             species)
  if (length(species) > 1L) origins[length(species)] <- min(max(times) / 2, 250)
  phylo_mask <- expand.grid(species = species, time_ma = times,
                            KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  phylo_mask$phylo_existence <- as.numeric(
    phylo_mask$time_ma <= origins[phylo_mask$species]
  )
  projection$phylo_existence <- phylo_mask$phylo_existence[
    match(paste(projection$species, projection$time_ma),
          paste(phylo_mask$species, phylo_mask$time_ma))
  ]

  pair <- expand.grid(from_region = regions, to_region = regions,
                      time_ma = times, KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)
  pair$paleodistance <- ifelse(pair$from_region == pair$to_region, 0,
                               15 + 25 * abs(match(pair$from_region, regions) -
                                               match(pair$to_region, regions)) +
                                 0.04 * pair$time_ma)
  pair$barrier_strength <- ifelse(pair$from_region == pair$to_region, 0,
                                  .hee_clip01(0.15 + 0.5 *
                                                (pair$time_ma %in% c(520, 500))))
  pair$route_open <- ifelse(pair$from_region == pair$to_region, 1,
                            .hee_clip01(1 - pair$barrier_strength))
  pair$climate_connectivity <- ifelse(pair$from_region == pair$to_region, 1,
                                      .hee_clip01(0.8 - 0.001 * pair$paleodistance))
  if (scenario == "extreme") {
    pair$paleodistance[pair$from_region != pair$to_region] <-
      pair$paleodistance[pair$from_region != pair$to_region] * 4
    pair$barrier_strength[pair$from_region != pair$to_region] <- 0.95
    pair$route_open[pair$from_region != pair$to_region] <- 0.03
    pair$climate_connectivity[pair$from_region != pair$to_region] <- 0.05
  }
  conn <- hee_build_connectivity_cube(pair, traits = traits,
                                      distance_scale = 50)

  accessibility <- expand.grid(species = species, region = regions,
                               time_ma = times, KEEP.OUT.ATTRS = FALSE,
                               stringsAsFactors = FALSE)
  accessibility$accessibility <- .hee_clip01(
    0.45 + 0.15 * match(accessibility$species, species) -
      0.1 * (accessibility$region == "R3") -
      0.0003 * accessibility$time_ma
  )

  extrapolation <- unique(projection[, c("cell_id", "time_ma",
                                         "extrapolation_score"), drop = FALSE])
  extrapolation$extrapolation_flag <- extrapolation$extrapolation_score > 0.7
  extrapolation$Q_extrap <- .hee_clip01(1 - extrapolation$extrapolation_score)

  tip_ranges <- matrix(0, nrow = length(species), ncol = length(regions),
                       dimnames = list(species, regions))
  tip_ranges[, "R1"] <- 1
  if (length(regions) > 1L) tip_ranges[seq_along(species) %% 2 == 0, "R2"] <- 1

  bsm_events <- data.frame(
    time_ma = rep(stats::median(times), 4),
    event_type = c("speciation", "range_expansion",
                   "local_extinction", "founder_event"),
    region = c("R1", "R2", "R3", "R2"),
    from_region = c("R1", "R1", "R3", "R1"),
    to_region = c("R1", "R2", "R3", "R2"),
    probability = c(0.4, 0.5, 0.3, 0.2),
    stringsAsFactors = FALSE
  )

  if (scenario == "missing_optional") {
    conn <- NULL
    accessibility <- NULL
    extrapolation <- NULL
    bsm_events <- NULL
  }
  if (scenario == "shuffled_time") {
    landscape <- landscape[sample(seq_len(nrow(landscape))), , drop = FALSE]
    projection <- projection[sample(seq_len(nrow(projection))), , drop = FALSE]
    if (!is.null(conn)) conn <- conn[sample(seq_len(nrow(conn))), , drop = FALSE]
  }

  list(
    times = times,
    env_cols = c("bio1", "bio12", "habitat_heterogeneity"),
    landscape_state = landscape,
    projection = projection,
    connectivity_cube = conn,
    accessibility = accessibility,
    phylo_mask = phylo_mask,
    extrapolation = extrapolation,
    tip_ranges = tip_ranges,
    species_origin = origins,
    bsm_events = bsm_events,
    traits = traits
  )
}

.hee_geoprocess_scenario_derived <- function(geo, inputs) {
  dyn <- geo$dynamic_occupancy
  proj <- dyn
  names(proj)[names(proj) == "region"] <- "cell_id"
  proj$species <- as.character(proj$species)
  proj$probability <- proj$occupancy_probability
  richness <- hee_richness(proj, group_cols = c("cell_id", "time_ma"))
  turnover <- if (length(unique(proj$time_ma)) > 1L) hee_turnover(proj) else
    data.frame()
  refugia <- hee_refugia(richness, group_col = "cell_id")
  landscape_events <- hee_landscape_events(
    inputs$landscape_state,
    unit_col = "cell_id",
    time_col = "time_ma",
    land_col = "geographic_existence",
    habitat_col = "habitat_availability",
    area_col = "land_area_km2",
    elevation_col = "elev_m"
  )
  network <- if (!is.null(inputs$connectivity_cube)) {
    hee_network_metrics(inputs$connectivity_cube, threshold = 0.2)
  } else {
    list(node_metrics = data.frame(),
         network_summary = data.frame(),
         corridor_persistence = data.frame())
  }
  extinction_layers <- hee_extinction_layers(
    geo$dynamic_occupancy,
    geo$extinction_probability
  )
  cradle_museum_grave <- hee_cradle_museum_grave(
    geo$species_pool,
    geo$source_sink_summary
  )
  list(
    landscape_events = landscape_events,
    richness = richness,
    turnover = turnover,
    refugia = refugia,
    network_node_metrics = network$node_metrics,
    network_summary = network$network_summary,
    corridor_persistence = network$corridor_persistence,
    local_extinction = extinction_layers$local_extinction,
    regional_extinction = extinction_layers$regional_extinction,
    lineage_extinction = extinction_layers$lineage_extinction,
    cradle_museum_grave = cradle_museum_grave,
    dynamic_model_family = hee_dynamic_model_family(),
    mechanism_hypotheses = hee_geoprocess_hypotheses()
  )
}

.hee_validate_geoprocess_scenario <- function(scenario, inputs, results) {
  rows <- list()
  add <- function(module, check, ok, message, severity = "major",
                  interpretation_type = "automatic_validation") {
    rows[[length(rows) + 1L]] <<- data.frame(
      scenario = scenario,
      module = module,
      check = check,
      status = if (isTRUE(ok)) "PASS" else "FAIL",
      severity = if (isTRUE(ok)) "info" else severity,
      message = message,
      interpretation_type = interpretation_type,
      stringsAsFactors = FALSE
    )
  }
  bounded_cols <- c("geographic_existence", "habitat_availability",
                    "ecological_opportunity", "structural_connectivity",
                    "functional_connectivity", "climate_connectivity",
                    "connectivity", "isolation", "source_pressure",
                    "rescue_effect", "colonisation_probability",
                    "extinction_probability", "speciation_opportunity",
                    "occupancy_probability", "prediction_reliability",
                    "expected_richness", "turnover", "refugia_score",
                    "edge_density", "largest_component_fraction",
                    "fragmentation_index", "mean_connectivity",
                    "centrality_index", "corridor_persistence",
                    "local_extinction_risk", "regional_occupancy",
                    "regional_absence_probability",
                    "regional_extinction_proxy",
                    "lineage_extinction_proxy",
                    "lineage_persistence_proxy",
                    "cradle_score", "museum_score", "grave_score",
                    "source_score", "sink_score")
  for (nm in names(results)) {
    d <- results[[nm]]
    if (!is.data.frame(d) || nrow(d) == 0L) next
    cols <- intersect(bounded_cols, names(d))
    if (length(cols) == 0L) next
    vals <- unlist(d[, cols, drop = FALSE], use.names = FALSE)
    vals <- suppressWarnings(as.numeric(vals))
    finite_vals <- vals[!is.na(vals) & is.finite(vals)]
    finite_only <- intersect(cols, "refugia_score")
    nonnegative_only <- intersect(cols, "expected_richness")
    probability_only <- setdiff(cols, c(finite_only, nonnegative_only))
    ok <- all(is.na(vals) | is.finite(vals)) &&
      (length(finite_vals) == 0L ||
         all(is.finite(finite_vals)))
    if (length(nonnegative_only) > 0L) {
      nvals <- suppressWarnings(as.numeric(unlist(
        d[, nonnegative_only, drop = FALSE], use.names = FALSE
      )))
      ok <- ok && all(is.na(nvals) | (is.finite(nvals) & nvals >= -1e-8))
    }
    if (length(probability_only) > 0L) {
      pvals <- suppressWarnings(as.numeric(unlist(
        d[, probability_only, drop = FALSE], use.names = FALSE
      )))
      ok <- ok && all(is.na(pvals) | (is.finite(pvals) &
                                        pvals >= -1e-8 & pvals <= 1 + 1e-8))
    }
    add(nm, "bounded_or_finite_outputs", ok,
        paste("Checked columns:", paste(cols, collapse = ", ")))
  }

  dyn <- results$dynamic_occupancy
  if (is.data.frame(dyn) && nrow(dyn) > 0L) {
    add("dynamic_occupancy", "unique_species_region_time",
        !any(duplicated(dyn[, c("species", "region", "time_ma"),
                            drop = FALSE])),
        "Dynamic occupancy keeps one row per species-region-time.")
    lost <- dyn$geographic_existence <= 0
    add("dynamic_occupancy", "geographic_absence_forces_zero",
        !any(lost & dyn$occupancy_probability > 1e-9, na.rm = TRUE),
        "Rows with geographic existence equal to zero have zero occupancy.")
  }

  if (scenario == "land_appearance_loss") {
    events <- results$landscape_events
    add("landscape_events", "detects_emergence_and_submergence",
        is.data.frame(events) &&
          all(c("emergence", "submergence") %in% events$event_type),
        "Sudden land appearance/loss should generate emergence and submergence events.",
        severity = "critical")
  }
  if (is.data.frame(results$network_summary) &&
      nrow(results$network_summary) > 0L) {
    add("network_summary", "network_metrics_available",
        all(c("fragmentation_index", "largest_component_fraction",
              "network_components") %in% names(results$network_summary)),
        "Dynamic connectivity network metrics are available for each time slice.")
  }
  if (is.data.frame(results$lineage_extinction) &&
      nrow(results$lineage_extinction) > 0L) {
    add("extinction_layers", "separates_local_regional_lineage",
        all(c("expected_occupied_regions", "lineage_extinction_proxy",
              "lineage_persistence_proxy") %in%
              names(results$lineage_extinction)),
        "Extinction diagnostics are split into local, regional, and lineage proxies.")
  }
  if (is.data.frame(results$dynamic_model_family) &&
      nrow(results$dynamic_model_family) > 0L) {
    add("model_family", "contains_M0_to_M8",
        identical(results$dynamic_model_family$model_id, paste0("M", 0:8)),
        "Dynamic Earth-Biota model family contains M0-M8.")
  }
  if (is.data.frame(results$mechanism_hypotheses) &&
      nrow(results$mechanism_hypotheses) > 0L) {
    add("hypotheses", "contains_H1_to_H8",
        identical(results$mechanism_hypotheses$hypothesis_id, paste0("H", 1:8)),
        "Mechanism hypothesis table contains H1-H8.")
  }
  if (scenario == "single_time") {
    add("time_axis", "single_time_slice_allowed",
        length(inputs$times) == 1L && is.data.frame(results$process_summary),
        "Single time-slice diagnostics should return static summaries.")
  }
  if (scenario == "uneven_time") {
    ps <- results$process_summary
    add("time_axis", "uneven_intervals_are_recorded",
        "interval_myr" %in% names(results$region_state) &&
          any(results$region_state$interval_myr !=
                results$region_state$interval_myr[1], na.rm = TRUE),
        "Unequal adjacent Ma intervals should propagate into interval_myr.")
  }
  if (scenario == "shuffled_time") {
    add("time_axis", "shuffled_input_order_does_not_break_outputs",
        is.data.frame(results$dynamic_occupancy) &&
          nrow(results$dynamic_occupancy) > 0L,
        "Unsorted inputs still produce dynamic outputs after internal sorting.")
  }
  if (scenario == "missing_optional") {
    add("optional_inputs", "missing_optional_processes_are_neutral_or_explicit",
        is.data.frame(results$source_pressure) &&
          all(results$source_pressure$source_pressure == 0, na.rm = TRUE),
        "Without connectivity, source pressure is explicitly zero rather than guessed.")
  }
  .hee_bind_fill(rows)
}

.hee_geoprocess_formula_catalog <- function() {
  core <- hee_hmscee_process_catalog()
  core_rows <- data.frame(
    process = paste0("core_", core$process_slug),
    formula = core$main_hmscee_responsibility,
    interpretation_type = core$layer,
    external_data_needed = core$typical_inputs,
    causal_warning = core$interpretation_boundary,
    core_process_mapping = paste(core$process_name_zh, core$process_name_en,
                                 sep = " / "),
    stringsAsFactors = FALSE
  )
  legacy <- data.frame(
    process = c("static_projection_M1", "static_projection_M2",
                "static_projection_M3", "static_projection_M4",
                "static_projection_M5", "source_pressure",
                "colonisation_probability", "extinction_probability",
                "dynamic_occupancy", "structural_connectivity",
                "functional_connectivity", "climatic_connectivity",
                "prediction_reliability", "network_fragmentation",
                "extinction_layers", "cradle_museum_grave",
                paste0("dynamic_model_", 0:8), paste0("hypothesis_", 1:8)),
    formula = c(
      "P = S_HMSC * L_land",
      "P = S_HMSC * E_phylo * L_land",
      "P = S_HMSC * E_phylo * A_BGB * L_land",
      "P = S_HMSC * E_phylo * A_BGB * D_static * L_land",
      "P = S_HMSC * E_phylo * A_BGB * D_dynamic * L_land",
      "M[j,c,t] = 1 - product(1 - p[j,c',t] * K[d] * C[j,c',c,t])",
      "gamma = 1 - exp(-exp(eta_col) * delta_t)",
      "epsilon = 1 - exp(-exp(eta_ext) * delta_t)",
      "p[t+1] = G * E * (p[t] * (1 - epsilon) + (1 - p[t]) * gamma)",
      "C_struct = route_open * exp(-alpha * distance) * exp(-barrier_weight * barrier)",
      "C_func = plogis(k0 + k_trait * trait - k_distance * distance - k_barrier * barrier)",
      "C_clim = exp(-phi * climate_mismatch) or supplied corridor suitability",
      "PRI = (1 - extrapolation_risk)^w1 * model_agreement^w2 * phylo_validity^w3 * calibration^w4",
      "fragmentation = 1 - largest connected component / total nodes",
      "local epsilon; regional absence = 1 - mean occupancy; lineage absence = I(sum occupancy <= 0)",
      "status = argmax(cradle, museum, grave, source, sink proxy scores)",
      hee_dynamic_model_family()$formula,
      hee_geoprocess_hypotheses()$prediction
    ),
    interpretation_type = c(
      rep("multiplicative_scenario", 5),
      "proxy_index", "heuristic_transition_probability",
      "heuristic_transition_probability", "dynamic_state_recursion",
      "proxy_connectivity", "proxy_connectivity", "proxy_connectivity",
      "epistemic_diagnostic", "proxy_connectivity",
      "proxy_extinction_summary", "proxy_region_role",
      hee_dynamic_model_family()$implementation_status,
      rep("testable_mechanism_hypothesis", 8)
    ),
    external_data_needed = c(
      "HMSC suitability and land mask",
      "species origin/extinction times or dated tree",
      "BioGeoBEARS/accessibility or equivalent historical range model",
      "static palaeodistance or dispersal surface",
      "dynamic dispersal filter from prior state and connectivity",
      "previous occupancy and connectivity",
      "calibrated coefficients for inference",
      "calibrated coefficients for inference",
      "colonisation, extinction, phylo and geography masks",
      "palaeodistance/barrier reconstruction",
      "dispersal traits and palaeodistance/barriers",
      "climate corridor or environmental mismatch",
      "diagnostics and optional external validation",
      "region-pair connectivity through time",
      "dynamic occupancy and extinction probability, or observed extinction events",
      "species-pool and source-sink summaries",
      hee_dynamic_model_family()$required_components,
      hee_geoprocess_hypotheses()$proxy_outputs
    ),
    core_process_mapping = c(
      "Legacy environmental_filtering sensitivity scenario",
      "Legacy environmental_filtering plus speciation/extinction lineage-time mask sensitivity scenario",
      "Legacy dispersal/speciation regional-history sensitivity scenario",
      "Legacy dispersal static within-region filter sensitivity scenario",
      "Legacy dispersal dynamic within-region filter sensitivity scenario",
      "colonisation arrival pressure proxy",
      "colonisation transition proxy",
      "environmental_filtering/dispersal/biotic_filtering/persistence local loss proxy",
      "environmental_filtering/dispersal/colonisation/biotic_filtering/persistence local occupancy recursion",
      "dispersal movement-resistance proxy",
      "dispersal movement trait proxy",
      "environmental_filtering/dispersal climate-corridor proxy",
      "Uncertainty layer, not a biological process",
      "dispersal graph diagnostic",
      "environmental_filtering/dispersal/biotic_filtering/speciation/extinction result proxy",
      "Result diagnostic",
      hee_dynamic_model_family()$core_process_mapping,
      hee_geoprocess_hypotheses()$core_process_mapping
    ),
    causal_warning = c(
      rep("Scenario product, not causal proof.", 5),
      "Source pressure is a model component, not observed propagules.",
      "Use as a calibrated rate only with independent calibration.",
      "Use as a calibrated rate only with independent calibration.",
      "Dynamic recursion is conditional on supplied process tables.",
      "Geological proxy unless calibrated to dispersal evidence.",
      "Trait-mediated proxy unless calibrated.",
      "Corridor proxy unless derived from explicit dispersal paths.",
      "Reliability index is not a credible interval.",
      "Network fragmentation is a graph proxy, not observed dispersal.",
      "Layered extinction terms are diagnostic proxies unless externally observed.",
      "Regional labels are descriptive roles, not causal classification.",
      rep("Model-family definition; fitting each model requires supplied components.", 9),
      hee_geoprocess_hypotheses()$interpretation_limit
    ),
    stringsAsFactors = FALSE
  )
  .hee_bind_fill(list(core_rows, legacy))
}

.hee_geoprocess_interpretation_table <- function() {
  data.frame(
    module = c("external_drivers", "eight_core_processes",
               "historical_constraints", "observation_uncertainty",
               "arena", "land_age", "area_heterogeneity", "connectivity",
               "network_metrics", "isolation", "colonisation",
               "rescue_extinction", "extinction_layers", "species_pool",
               "cradle_museum_grave", "speciation_opportunity",
               "dynamic_model_family", "mechanism_hypotheses",
               "dynamic_occupancy",
               "richness_turnover_refugia", "limitation", "reliability",
               "external_processes"),
    interpretation_type = c("driver_layer",
                            "core_process_catalog",
                            "inference_constraint_layer",
                            "uncertainty_validation_layer",
                            "external_geological_input_or_proxy",
                            "proxy_index",
                            "proxy_index",
                            "proxy_connectivity",
                            "proxy_connectivity",
                            "proxy_index",
                            "heuristic_transition_probability",
                            "heuristic_transition_probability",
                            "proxy_extinction_summary",
                            "event_accounting_or_external_BGB",
                            "proxy_region_role",
                            "proxy_index",
                            "model_family_definition",
                            "testable_mechanism_hypotheses",
                            "dynamic_state_recursion",
                            "derived_community_summary",
                            "rule_based_diagnostic",
                            "epistemic_diagnostic",
                            "requires_external_data"),
    scientific_explanation = c(
      "Climate, geology, geography and habitat history define X, H, W and cell-region maps; they condition the eight biological processes but are not extra processes.",
      "The core processes are environmental filtering, dispersal, colonisation, biotic filtering, persistence, evolution, speciation, and extinction.",
      "BioGeoBEARS/BSM and dated trees constrain dispersal, speciation and extinction histories; they are inference tools rather than extra processes.",
      "Fossils, pollen, aDNA, no-analog diagnostics, posterior uncertainty and phyloregion summaries validate or interpret outputs; they are not biological mechanisms.",
      "Land, habitat, area, elevation and component histories define the ecological stage; the package detects changes in supplied maps but does not infer tectonics.",
      "Newly emerged habitat starts at zero age and persistent habitat accumulates true Ma intervals.",
      "Area and environmental heterogeneity can raise opportunity and lower extinction risk, but are not direct causes without calibration.",
      "Structural, functional and climatic connectivity separate physical routes, species traits and climatic corridors.",
      "Network metrics summarise components, edge density, corridor persistence and fragmentation in the supplied dynamic graph.",
      "Isolation is not simply one minus connectivity; duration matters for speciation and extinction risk.",
      "Colonisation requires source pressure, accessibility, suitability and a time interval; omitted source pressure is treated as no propagule evidence and yields zero colonisation.",
      "Rescue reduces extinction risk; land disappearance forces local geographic extinction.",
      "Extinction is separated into local probability, regional absence, and lineage-level persistence proxies; only external event data can make these observed extinctions.",
      "Regional species pools update through speciation, immigration, emigration and extinction events.",
      "Cradle, museum, grave, source and sink labels are transparent regional role proxies from pool change and source-sink balance.",
      "Allopatric, founder-event and in-situ opportunities are indices, not estimated speciation rates.",
      "M0-M8 define a nested mechanism family from climate-only to feedback-enabled Earth-Biota assembly; the table is a design map, not automatic model fitting.",
      "H1-H8 list mechanism hypotheses and the package outputs that can act as proxy tests.",
      "Dynamic occupancy combines previous state, colonisation, persistence, geography and lineage existence.",
      "Richness, turnover and refugia are summaries derived from probabilities.",
      "Dominant limitation labels the lowest or most restrictive supplied component and is not causal proof.",
      "Prediction reliability combines extrapolation and agreement diagnostics; it is not posterior uncertainty.",
      "Real plate motion, BioGeoBEARS fitting, fossil-calibrated extinction and true speciation rates require external data or specialist models."
    ),
    stringsAsFactors = FALSE
  )
}

.hee_bind_fill <- function(rows) {
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0L) return(data.frame())
  cols <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(x) {
    x <- as.data.frame(x)
    miss <- setdiff(cols, names(x))
    for (m in miss) x[[m]] <- rep(NA, nrow(x))
    x[, cols, drop = FALSE]
  })
  do.call(rbind, rows)
}
