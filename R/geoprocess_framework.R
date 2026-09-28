#' Summarise dynamic connectivity networks
#'
#' Computes graph-level proxy metrics from a region-pair connectivity cube:
#' weighted degree, network component count, largest-component fraction,
#' fragmentation index, and corridor persistence. The function uses only the
#' supplied edge table; it does not infer palaeogeography or plate motion.
#'
#' @param connectivity_cube Pairwise connectivity table.
#' @param threshold Connectivity threshold below which an edge is treated as
#'   closed for component and corridor-open summaries.
#' @param time_col,from_col,to_col,species_col Column names.
#' @param connectivity_col Optional connectivity column. If `NULL`, the first
#'   available column among `connectivity`, `total_connectivity`,
#'   `functional_connectivity`, and `structural_connectivity` is used.
#' @return A list with `node_metrics`, `network_summary`, and
#'   `corridor_persistence` data.frames.
#' @export
#'
#' @examples
#' x <- expand.grid(from_region = c("A", "B"), to_region = c("A", "B"),
#'   time_ma = c(10, 0), stringsAsFactors = FALSE)
#' x$connectivity <- c(1, 0.4, 0.4, 1, 1, 0.1, 0.1, 1)
#' hee_network_metrics(x)
hee_network_metrics <- function(connectivity_cube,
                                threshold = 0.2,
                                time_col = "time_ma",
                                from_col = "from_region",
                                to_col = "to_region",
                                species_col = "species",
                                connectivity_col = NULL) {
  x <- as.data.frame(connectivity_cube)
  .require_cols(x, c(from_col, to_col, time_col), "connectivity_cube")
  .hee_validate_time_values(x[[time_col]], paste0("connectivity_cube$", time_col))
  if (is.null(connectivity_col)) {
    connectivity_col <- intersect(c("connectivity", "total_connectivity",
                                    "functional_connectivity",
                                    "structural_connectivity"), names(x))[1]
  }
  if (is.na(connectivity_col) || !connectivity_col %in% names(x)) {
    stop("connectivity_cube must contain a connectivity column.", call. = FALSE)
  }
  keys <- c(if (species_col %in% names(x)) species_col else character(),
            from_col, to_col, time_col)
  .hee_check_unique_keys(x, keys, "connectivity_cube")
  x[[connectivity_col]] <- .hee_clip01(x[[connectivity_col]])
  x$.species_group <- if (species_col %in% names(x)) as.character(x[[species_col]]) else "all_species"

  groups <- split(x, interaction(x$.species_group, x[[time_col]], drop = TRUE),
                  drop = TRUE)
  node_rows <- list()
  summary_rows <- list()
  nr <- 0L
  sr <- 0L
  for (g in groups) {
    sp <- g$.species_group[1]
    tt <- as.numeric(g[[time_col]][1])
    nodes <- sort(unique(c(as.character(g[[from_col]]),
                           as.character(g[[to_col]]))))
    if (length(nodes) == 0L) next
    open <- g[g[[connectivity_col]] >= threshold &
                as.character(g[[from_col]]) != as.character(g[[to_col]]),
              , drop = FALSE]
    adj <- stats::setNames(vector("list", length(nodes)), nodes)
    if (nrow(open) > 0L) {
      for (i in seq_len(nrow(open))) {
        a <- as.character(open[[from_col]][i])
        b <- as.character(open[[to_col]][i])
        adj[[a]] <- unique(c(adj[[a]], b))
        adj[[b]] <- unique(c(adj[[b]], a))
      }
    }
    comp_id <- stats::setNames(rep(NA_integer_, length(nodes)), nodes)
    comp_n <- 0L
    for (nd in nodes) {
      if (!is.na(comp_id[[nd]])) next
      comp_n <- comp_n + 1L
      queue <- nd
      comp_id[[nd]] <- comp_n
      while (length(queue) > 0L) {
        cur <- queue[1]
        queue <- queue[-1]
        nxt <- adj[[cur]]
        for (nb in nxt) {
          if (is.na(comp_id[[nb]])) {
            comp_id[[nb]] <- comp_n
            queue <- c(queue, nb)
          }
        }
      }
    }
    weighted_degree <- stats::aggregate(
      g[[connectivity_col]],
      data.frame(node = as.character(g[[from_col]]), stringsAsFactors = FALSE),
      sum, na.rm = TRUE
    )
    names(weighted_degree) <- c("node", "weighted_degree")
    open_degree <- stats::aggregate(
      as.numeric(g[[connectivity_col]] >= threshold),
      data.frame(node = as.character(g[[from_col]]), stringsAsFactors = FALSE),
      sum, na.rm = TRUE
    )
    names(open_degree) <- c("node", "open_degree")
    node <- data.frame(
      species = sp,
      time_ma = tt,
      region = nodes,
      component_id = as.integer(comp_id[nodes]),
      stringsAsFactors = FALSE
    )
    node <- merge(node, weighted_degree, by.x = "region", by.y = "node",
                  all.x = TRUE, sort = FALSE)
    node <- merge(node, open_degree, by.x = "region", by.y = "node",
                  all.x = TRUE, sort = FALSE)
    node$weighted_degree[is.na(node$weighted_degree)] <- 0
    node$open_degree[is.na(node$open_degree)] <- 0
    node$centrality_index <- .hee_positive_index(node$weighted_degree)
    nr <- nr + 1L
    node_rows[[nr]] <- node

    comp_sizes <- table(comp_id)
    largest <- if (length(comp_sizes) == 0L) 0 else max(as.numeric(comp_sizes))
    possible_edges <- max(length(nodes) * (length(nodes) - 1), 1)
    open_edges <- nrow(open)
    sr <- sr + 1L
    summary_rows[[sr]] <- data.frame(
      species = sp,
      time_ma = tt,
      n_nodes = length(nodes),
      n_open_edges = open_edges,
      edge_density = open_edges / possible_edges,
      network_components = comp_n,
      largest_component_fraction = largest / length(nodes),
      fragmentation_index = .hee_clip01(1 - largest / length(nodes)),
      mean_connectivity = mean(g[[connectivity_col]], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  node_metrics <- if (length(node_rows)) do.call(rbind, node_rows) else data.frame()
  network_summary <- if (length(summary_rows)) do.call(rbind, summary_rows) else data.frame()

  edge_keys <- c(".species_group", from_col, to_col)
  persistence <- stats::aggregate(
    as.numeric(x[[connectivity_col]] >= threshold),
    x[, edge_keys, drop = FALSE],
    mean, na.rm = TRUE
  )
  names(persistence)[ncol(persistence)] <- "corridor_persistence"
  names(persistence)[names(persistence) == ".species_group"] <- "species"
  persistence$corridor_persistence <- .hee_clip01(persistence$corridor_persistence)
  list(
    node_metrics = node_metrics,
    network_summary = network_summary,
    corridor_persistence = persistence
  )
}

#' Separate local, regional, and lineage extinction proxies
#'
#' Summarises extinction at three ecological scales from dynamic occupancy and
#' extinction-probability tables. Regional and lineage extinction are proxy
#' summaries unless independently observed extinction events are supplied.
#'
#' @param dynamic_state Species-cell/region-time occupancy table.
#' @param extinction_probability Optional extinction probability table.
#' @param species_col,region_col,time_col,probability_col Column names.
#' @param extinction_col Extinction-probability column.
#' @param threshold Occupancy threshold used to define likely presence.
#' @return A list with `local_extinction`, `regional_extinction`, and
#'   `lineage_extinction` data.frames.
#' @export
#'
#' @examples
#' x <- expand.grid(species = "sp1", region = "R1", time_ma = c(10, 0))
#' x$occupancy_probability <- c(0.8, 0.1)
#' hee_extinction_layers(x, threshold = 0.5)
hee_extinction_layers <- function(dynamic_state,
                                  extinction_probability = NULL,
                                  species_col = "species",
                                  region_col = "region",
                                  time_col = "time_ma",
                                  probability_col = "occupancy_probability",
                                  extinction_col = "extinction_probability",
                                  threshold = 0.5) {
  x <- as.data.frame(dynamic_state)
  .require_cols(x, c(species_col, region_col, time_col, probability_col),
                "dynamic_state")
  .hee_validate_time_values(x[[time_col]], paste0("dynamic_state$", time_col))
  .hee_check_unique_keys(x, c(species_col, region_col, time_col),
                         "dynamic_state")
  x$.row_id <- seq_len(nrow(x))
  if (!is.null(extinction_probability)) {
    e <- as.data.frame(extinction_probability)
    if (extinction_col %in% names(e)) {
      joined <- .hee_safe_left_join(
        x, e, c(species_col, region_col, time_col), extinction_col,
        target_name = "dynamic_state",
        source_name = "extinction_probability"
      )
      x[[extinction_col]] <- joined[[extinction_col]]
    }
  }
  if (!extinction_col %in% names(x)) x[[extinction_col]] <- 0
  x[[probability_col]] <- .hee_clip01(x[[probability_col]])
  x[[extinction_col]] <- .hee_clip01(x[[extinction_col]])

  local <- x[, c(species_col, region_col, time_col, probability_col,
                 extinction_col), drop = FALSE]
  local$likely_present <- local[[probability_col]] >= threshold
  local$local_extinction_risk <- local[[extinction_col]]

  reg <- stats::aggregate(
    x[[probability_col]],
    x[, c(species_col, region_col, time_col), drop = FALSE],
    mean, na.rm = TRUE
  )
  names(reg)[ncol(reg)] <- "regional_occupancy"
  reg$regional_absence_probability <- .hee_clip01(1 - reg$regional_occupancy)
  reg$regional_extinction_proxy <- reg$regional_absence_probability

  lin <- stats::aggregate(
    x[[probability_col]],
    x[, c(species_col, time_col), drop = FALSE],
    sum, na.rm = TRUE
  )
  names(lin)[ncol(lin)] <- "expected_occupied_regions"
  lin$lineage_extinction_proxy <- as.numeric(lin$expected_occupied_regions <=
                                               .Machine$double.eps)
  max_occ <- max(lin$expected_occupied_regions, na.rm = TRUE)
  if (!is.finite(max_occ) || max_occ <= 0) max_occ <- 1
  lin$lineage_persistence_proxy <- .hee_clip01(lin$expected_occupied_regions /
                                                 max_occ)
  list(local_extinction = local,
       regional_extinction = reg,
       lineage_extinction = lin)
}

#' Classify regions as cradle, museum, grave, source, or sink proxies
#'
#' Builds a transparent regional status table from species-pool or source-sink
#' summaries. These labels are diagnostic proxies and should not be interpreted
#' as causal proof without independent fossil or phylogenetic calibration.
#'
#' @param species_pool Region-time species-pool table.
#' @param source_sink Optional source-sink summary table.
#' @param region_col,time_col Column names.
#' @return A data.frame with cradle/museum/grave/source/sink proxy scores and a
#'   dominant status label.
#' @export
#'
#' @examples
#' pool <- data.frame(region = "A", time_ma = c(10, 0),
#'   species_pool_size = c(2, 5), origination = c(1, 0),
#'   regional_extinction = c(0, 0.2), source_sink_balance = c(0.5, -0.2))
#' hee_cradle_museum_grave(pool)
hee_cradle_museum_grave <- function(species_pool,
                                    source_sink = NULL,
                                    region_col = "region",
                                    time_col = "time_ma") {
  x <- as.data.frame(species_pool)
  .require_cols(x, c(region_col, time_col), "species_pool")
  .hee_validate_time_values(x[[time_col]], paste0("species_pool$", time_col))
  .hee_check_unique_keys(x, c(region_col, time_col), "species_pool")
  x$.row_id <- seq_len(nrow(x))
  get <- function(nm, default = 0) {
    if (nm %in% names(x)) suppressWarnings(as.numeric(x[[nm]])) else
      rep(default, nrow(x))
  }
  x$cradle_score <- .hee_positive_index(get("origination") +
                                          get("speciation_opportunity"))
  pool_size <- get("species_pool_size")
  x$museum_score <- .hee_positive_index(pool_size)
  x$grave_score <- .hee_positive_index(get("regional_extinction") +
                                         get("lineage_extinction_proxy"))
  x$source_score <- .hee_positive_index(pmax(get("source_sink_balance"), 0))
  x$sink_score <- .hee_positive_index(pmax(-get("source_sink_balance"), 0))
  score_cols <- c("cradle_score", "museum_score", "grave_score",
                  "source_score", "sink_score")
  status <- c("cradle", "museum", "grave", "source", "sink")
  score_matrix <- as.matrix(x[, score_cols, drop = FALSE])
  score_sum <- rowSums(score_matrix, na.rm = TRUE)
  idx <- max.col(score_matrix, ties.method = "first")
  x$dominant_region_status <- ifelse(score_sum > 0, status[idx],
                                     "undetermined")
  if (!is.null(source_sink)) {
    ss <- as.data.frame(source_sink)
    keep <- intersect(c(region_col, time_col, "source_pressure",
                        "rescue_effect"), names(ss))
    if (all(c(region_col, time_col) %in% keep)) {
      joined <- .hee_safe_left_join(
        x, ss, c(region_col, time_col),
        setdiff(keep, c(region_col, time_col)),
        target_name = "species_pool",
        source_name = "source_sink"
      )
      for (nm in setdiff(keep, c(region_col, time_col))) {
        x[[nm]] <- joined[[nm]]
      }
    }
  }
  x$.row_id <- NULL
  x
}

#' Define H1-H8 geological-driver hypotheses
#'
#' Returns geological-driver hypotheses and maps each hypothesis to package
#' outputs that can be used as proxy tests. These hypotheses are not additional
#' biological process engines. They ask how climate, geology, geography and
#' habitat history condition the eight explicit HmscEE processes:
#' environmental filtering, dispersal, colonisation, biotic filtering,
#' persistence, evolution, speciation, and extinction.
#'
#' @return A data.frame with hypothesis, prediction, indicators, and limits.
#' @export
#'
#' @examples
#' hee_geoprocess_hypotheses()
hee_geoprocess_hypotheses <- function() {
  data.frame(
    hypothesis_id = paste0("H", 1:8),
    hypothesis = c(
      "arena_generation",
      "isolation_divergence",
      "connectivity_rescue",
      "intermediate_connectivity",
      "land_loss_extinction",
      "corridor_pulse",
      "trait_mediation",
      "post_crisis_radiation"
    ),
    prediction = c(
      "New land or habitat first increases empty opportunity and colonisation, then richness after a lag.",
      "Persistent isolation raises allopatric opportunity but can also raise extinction risk.",
      "High connectivity raises colonisation and rescue while lowering local extinction.",
      "Diversification opportunity peaks at intermediate connectivity, not at the extremes.",
      "Submergence, area loss, and habitat loss increase local/regional extinction and may create extinction debt.",
      "Opening and closing corridors create pulses of dispersal, turnover, and source-sink switching.",
      "Dispersal and vulnerability traits modulate connectivity, colonisation, rescue, and extinction.",
      "Crisis or reorganisation releases ecological space and can increase radiation opportunity in surviving lineages."
    ),
    proxy_outputs = c(
      "landscape_events, land_age, ecological_opportunity, colonisation_probability, richness",
      "isolation_duration, allopatric_opportunity, extinction_probability",
      "connectivity, source_pressure, rescue_effect, extinction_probability",
      "network_summary, speciation_opportunity, cradle_museum_grave",
      "submergence, land_loss, regional_extinction, lineage_extinction",
      "corridor_persistence, source_sink_summary, turnover",
      "functional_connectivity, species traits, extinction vulnerability",
      "radiation_opportunity, ecological_opportunity, richness recovery"
    ),
    core_process_mapping = c(
      "Drivers feeding environmental_filtering, dispersal and colonisation",
      "Drivers feeding dispersal, speciation and extinction diagnostics",
      "Drivers feeding dispersal, colonisation and persistence",
      "Drivers feeding dispersal, speciation and extinction diagnostics",
      "Drivers forcing habitat loss in environmental_filtering/dispersal state transitions",
      "Drivers feeding dispersal-driven range expansion and turnover",
      "evolution traits modifying environmental_filtering/dispersal/colonisation transition components",
      "Post-event diagnostics for environmental_filtering/dispersal/speciation/extinction outcomes"
    ),
    interpretation_limit = c(
      rep("Proxy test unless calibrated with independent palaeo events, fossils, or simulations.", 8)
    ),
    stringsAsFactors = FALSE
  )
}

#' Define the M0-M8 geological scenario-diagnostic family
#'
#' Extends the simpler M1-M5 static projection set with a scenario-diagnostic
#' sequence. This table is retained for backward-compatible sensitivity
#' analyses and teaching workflows. It is not the HmscEE 1.0 core process list:
#' the core biological processes are the six explicit process names in
#' [hee_core_process_catalog()]. This function defines required components
#' for scenario comparisons; it does not fit all models automatically.
#'
#' @return A data.frame with model id, mechanism, and required components.
#' @export
#'
#' @examples
#' hee_dynamic_model_family()
hee_dynamic_model_family <- function() {
  data.frame(
    model_id = paste0("M", 0:8),
    model_name = c(
      "climate_only",
      "climate_plus_arena",
      "arena_plus_network",
      "dynamic_colonisation",
      "local_extinction_rescue",
      "lineage_time_extinction",
      "speciation_opportunity",
      "trait_niche_evolution",
      "full_earth_biota_feedback"
    ),
    formula = c(
      "legacy diagnostic: S",
      "legacy diagnostic: S constrained by hard geography/habitat G",
      "legacy diagnostic: M1 plus supplied network/connectivity proxy C",
      "legacy diagnostic recursion: G * E * ((1 - p) * gamma + p)",
      "legacy diagnostic recursion: G * E * (p * (1 - epsilon) + (1 - p) * gamma)",
      "legacy diagnostic: M4 with lineage-time mask from dated tree/external extinction evidence",
      "post hoc proxy diagnostics for allopatric, founder, in-situ and radiation opportunity",
      "evolution/environmental_filtering diagnostic from beta/gamma trait-niche evolution outputs",
      "external-only sensitivity: optional feedback diagnostics, not estimated from one modern community table"
    ),
    required_components = c(
      "suitability",
      "suitability; geographic_existence",
      "suitability; geographic_existence; connectivity",
      "previous occupancy; colonisation_probability",
      "previous occupancy; colonisation_probability; extinction_probability; rescue_effect",
      "phylo_time_mask; species_origin; extinction_time if known",
      "speciation_opportunity",
      "beta_niche_traits; gamma_shift; trait mediation",
      "external feedback path or precomputed feedback diagnostics"
    ),
    implementation_status = c(
      "legacy_sensitivity_diagnostic",
      "legacy_sensitivity_diagnostic",
      "legacy_or_proxy_diagnostic",
      "legacy_dynamic_recursion",
      "legacy_dynamic_recursion",
      "lineage_constraint_diagnostic",
      "post_hoc_proxy_index",
      "evolution_to_environmental_filtering_diagnostic_from_HMSC_outputs",
      "external_or_sensitivity_only"
    ),
    framework_role = rep("scenario_diagnostic_not_core_process", 9),
    core_process_mapping = c(
      "environmental_filtering",
      "Drivers constraining environmental_filtering/dispersal state",
      "dispersal proxy",
      "colonisation state transition",
      "environmental_filtering/dispersal/biotic_filtering/persistence local state transition",
      "speciation/extinction lineage-time constraint",
      "speciation opportunity proxy, not a speciation-rate model",
      "evolution feeding environmental_filtering",
      "External feedback hypothesis, not estimated by default"
    ),
    interpretation_boundary = c(
      "Use the nested HmscEE core for final deep-time occupancy inference.",
      "Geography/habitat is a hard state or driver, not a sixth process.",
      "Connectivity proxies require calibration before causal interpretation.",
      "Gamma is a transition probability only under the supplied scenario and coefficients.",
      "Epsilon is local cell loss/persistence, not full lineage extinction.",
      "Extinct unsampled lineages require fossils/extinct tips/FBD models.",
      "Opportunity indices cannot be reported as true speciation rates.",
      "Trait-mediated response is a model assumption and uncertainty source.",
      "Feedback terms require independent time-series or mechanistic evidence."
    ),
    stringsAsFactors = FALSE
  )
}

#' Compare geological-process model outputs across multiple metrics
#'
#' Creates a compact model-comparison table for M0-M8 or user-supplied model
#' outputs. Metrics are descriptive and may include validation and ecological
#' summaries when provided.
#'
#' @param model_outputs Named list of projection tables or metric tables.
#' @param fossil_validation Optional fossil-validation table with model id.
#' @param richness Optional richness table.
#' @param turnover Optional turnover table.
#' @param source_sink Optional source/sink table.
#' @param model_col,time_col Column names.
#' @return A data.frame with model-level comparison metrics.
#' @export
#'
#' @examples
#' p <- data.frame(model_id = "M1", probability = c(0.2, 0.8))
#' hee_geoprocess_model_comparison(list(M1 = p))
hee_geoprocess_model_comparison <- function(model_outputs,
                                            fossil_validation = NULL,
                                            richness = NULL,
                                            turnover = NULL,
                                            source_sink = NULL,
                                            model_col = "model_id",
                                            time_col = "time_ma") {
  if (!is.list(model_outputs) || length(model_outputs) == 0L) {
    stop("model_outputs must be a non-empty named list.", call. = FALSE)
  }
  if (is.null(names(model_outputs)) || any(!nzchar(names(model_outputs)))) {
    stop("model_outputs must be a named list so each result has a stable model_id.",
         call. = FALSE)
  }
  rows <- lapply(names(model_outputs), function(id) {
    x <- as.data.frame(model_outputs[[id]])
    prob_col <- intersect(c("probability", "occupancy_probability",
                            "suitability"), names(x))[1]
    prob <- if (!is.na(prob_col)) .hee_clip01(x[[prob_col]]) else numeric()
    data.frame(
      model_id = id,
      n_rows = nrow(x),
      mean_probability = if (length(prob)) mean(prob, na.rm = TRUE) else NA_real_,
      impossible_fraction = if (all(c("geographic_existence", prob_col) %in%
                                    names(x))) {
        mean(x$geographic_existence <= 0 & prob > 0, na.rm = TRUE)
      } else NA_real_,
      range_size_proxy = if (length(prob)) sum(prob, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  add_metric <- function(tab, value_col, out_col, fun = mean) {
    if (is.null(tab)) return()
    d <- as.data.frame(tab)
    if (!all(c(model_col, value_col) %in% names(d))) return()
    z <- stats::aggregate(d[[value_col]], d[, model_col, drop = FALSE],
                          fun, na.rm = TRUE)
    names(z) <- c("model_id", out_col)
    out <<- merge(out, z, by = "model_id", all.x = TRUE, sort = FALSE)
  }
  add_metric(fossil_validation, "predicted_probability",
             "mean_fossil_probability")
  add_metric(richness, "expected_richness", "mean_richness")
  add_metric(turnover, "turnover", "mean_turnover")
  add_metric(source_sink, "source_sink_balance", "mean_source_sink_balance")
  out
}
