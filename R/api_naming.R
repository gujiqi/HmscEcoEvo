# Public API naming layer ----------------------------------------------------

#' HmscEE core process catalog
#'
#' HmscEE exposes six biological process names for user-facing API design and
#' reporting: environmental filtering, dispersal, biotic filtering, evolution,
#' speciation, and extinction. Colonisation and local persistence are optional
#' state-transition terms in legacy dynamic-occupancy scenarios, not additional
#' processes or outputs of the current Case05 phylogeographic reconstruction.
#' Earth history,
#' BioGeoBEARS/BSM, dated trees, historical predictors, validation data and
#' uncertainty diagnostics provide inputs, constraints or outputs around these
#' process functions.
#'
#' @return A data.frame describing the core biological processes and the
#'   surrounding driver/constraint/result layers.
#' @export
#'
#' @examples
#' hee_core_process_catalog()
hee_core_process_catalog <- function() {
  data.frame(
    layer = c(
      "external_driver",
      rep("core_process", 6),
      "historical_information",
      "historical_inference_or_constraint",
      "observation_uncertainty",
      "result_diagnostic"
    ),
    process_slug = c(
      "dynamic_earth_driver",
      "environmental_filtering",
      "dispersal",
      "biotic_filtering",
      "evolution",
      "speciation",
      "extinction",
      "history_predictors",
      "bgb_tree_constraints",
      "validation_uncertainty",
      "biodiversity_outputs"
    ),
    process_name_zh = c(
      "\u52a8\u6001\u5730\u7403\u9a71\u52a8",
      "\u73af\u5883\u8fc7\u6ee4",
      "\u6269\u6563",
      "\u751f\u7269\u8fc7\u6ee4",
      "\u6f14\u5316",
      "\u7269\u79cd\u5f62\u6210",
      "\u706d\u7edd",
      "\u5386\u53f2\u9884\u6d4b\u53d8\u91cf",
      "\u5386\u53f2\u751f\u7269\u5730\u7406\u4e0e\u7cfb\u7edf\u6811\u7ea6\u675f",
      "\u9a8c\u8bc1\u4e0e\u4e0d\u786e\u5b9a\u6027",
      "\u751f\u7269\u591a\u6837\u6027\u7ed3\u679c"
    ),
    process_name_en = c(
      "Dynamic Earth drivers",
      "Environmental filtering",
      "Dispersal",
      "Biotic filtering",
      "Evolution",
      "Speciation",
      "Extinction",
      "Historical predictors",
      "Historical biogeography and tree constraints",
      "Validation and uncertainty",
      "Biodiversity outputs"
    ),
    function_prefix = c(
      "hee_earth_*",
      "hee_environmental_filtering_*",
      "hee_dispersal_*",
      "hee_biotic_filtering_*",
      "hee_evolution_*",
      "hee_speciation_*",
      "hee_extinction_*",
      "hee_history_*",
      "hee_bgb_* / hee_phylo_*",
      "hee_validation_* / hee_uncertainty_*",
      "hee_diversity_* / hee_result_*"
    ),
    biological_question = c(
      "What did the palaeo-Earth stage, habitat, region map and resistance surface look like?",
      "How does a lineage respond to the environment at this time?",
      "How does a lineage move on an explicit dynamic landscape grid?",
      "What independently supported effects do other organisms exert?",
      "How did traits and environmental responses change through the tree?",
      "When do tree nodes produce daughter lineages?",
      "How is disappearance represented at local, regional or lineage scale?",
      "How are legacy hmscHist historical indices organised for modern HMSC?",
      "How do BioGeoBEARS/BSM and dated trees constrain regional histories and node inheritance?",
      "How reliable are model predictions and extrapolation conditions?",
      "What patterns are derived from lineage-cell-time probability?"
    ),
    main_hmscee_responsibility = c(
      "Standardise X, H, W and the cell-region map; Earth drivers are not biological processes.",
      "Use HMSC and ancestral environmental responses to compute environmental support, not historical occupancy.",
      "Infer lineage locations or construct an active movement kernel; plate carriage is upstream geography, not dispersal.",
      "Do not treat HMSC residual association as causal interaction by default; add explicit external evidence only.",
      "Reconstruct ancestral traits, beta or beta = Gamma*T + u; do not use BGB range to define niche.",
      "Use dated trees for parent-to-daughter identity; external BGB node events are optional only in regional analyses.",
      "Separate geographic loss from biological extinction; extant-only trees cannot recover full extinct diversity.",
      "Preserve hmscHist historical indices, standardisation and HMSC input conversion.",
      "BioGeoBEARS answers which regions; HmscEE refines cells inside those regions.",
      "Fossils, CV, PPC, no-analog and posterior uncertainty audit or interpret outputs.",
      "Richness, turnover, refugia, PD, CWM and phyloregion are derived summaries, not mechanisms."
    ),
    typical_inputs = c(
      "palaeogeography g, palaeoclimate c|g, land, habitat, topography, resistance",
      "earth X, ancestral beta/response draws, HMSC posterior",
      "dynamic land/habitat mask, grid geometry, plate transport, sparse movement graph, dispersal traits, delta_t",
      "external interaction networks or justified interaction covariates",
      "dated tree, tip traits, HMSC beta/Gamma posterior, optional comparative-model outputs",
      "dated tree nodes, daughter labels, optional external node geography",
      "land history, fossils or extinct tips when available",
      "species-region histories, trait histories, event tables, modern sites",
      "tip ranges, region histories, BSM draws, dated tree",
      "fossils, pollen, aDNA, CV folds, posterior draws, no-analog metrics",
      "cell probabilities, traits, tree, region map, time axis"
    ),
    main_outputs = c(
      "E_k = {X_k,H_k,W_k}",
      "S_env[lineage,cell,time]",
      "K[from_cell,to_cell,time] and conditional lineage-location distributions",
      "B[lineage,cell,time] or explicit zero/omitted effect",
      "T[lineage,time], beta[lineage,time], response draws",
      "active lineages and node inheritance events",
      "geographic-loss and separately identified extinction summaries",
      "HMSC-ready XData, TrData, random effects and leakage diagnostics",
      "R[lineage,region,time] and node/range constraints",
      "uncertainty and validation tables reported beside biology",
      "richness, turnover, refugia, PD, CWM, functional diversity"
    ),
    interpretation_boundary = c(
      "Climate, geology, geography and habitat are upstream drivers, not independent biological probabilities.",
      "S_env is potential environmental support, not historical occupancy or a guaranteed fundamental niche.",
      "Use within-region movement when BioGeoBEARS owns cross-region history; use global explicit movement only in a no-BioGeoBEARS model. Do not treat plate carriage as active dispersal.",
      "Biotic filtering is optional; HMSC Omega/residual association is not causal evidence by itself.",
      "Evolution changes traits and responses; ancestral traits, responses and geographic ranges are separate objects.",
      "Speciation timing comes from tree/node evidence; opportunity indices are not speciation rates.",
      "Extinction must state scale; without fossil/extinct tips, full lineage extinction is not identifiable.",
      "History predictors can explain modern HMSC patterns; current-composition-derived site history must not leak into palaeo-grid prediction.",
      "BioGeoBEARS is a constraint/inference tool for dispersal/speciation histories, not an extra process.",
      "Uncertainty and extrapolation are not biological absence probabilities.",
      "Derived patterns can support hypotheses but require validation before causal claims."
    ),
    stringsAsFactors = FALSE
  )
}

#' HmscEE function catalog
#'
#' Provides a machine-readable map between the six-process API names,
#' older compatibility names, and non-process layers such as Earth drivers,
#' history predictors, BioGeoBEARS constraints, validation, uncertainty, and
#' biodiversity outputs.
#'
#' @return A data.frame with canonical names, legacy names, layers, process
#'   slugs, export status, and interpretation notes.
#' @export
#'
#' @examples
#' hee_function_catalog()
hee_function_catalog <- function() {
  rows <- list(
    c("hee_core_process_catalog", "hee_hmscee_process_catalog",
      "api_catalog", "api_catalog", "canonical",
      "Names the six HmscEE processes and surrounding layers."),
    c("hee_environmental_filtering_suitability", "hee_lineage_suitability",
      "core_process", "environmental_filtering", "canonical_alias",
      "Computes S from ancestral response and palaeoenvironment X."),
    c("hee_environmental_filtering_ancestral_suitability", NA,
      "core_process", "environmental_filtering", "canonical",
      "Runs HMSC tip Beta posterior -> ancestral Beta -> palaeoenvironmental support S_env."),
    c("hee_dispersal_landscape_weights", NA,
      "external_driver", "dynamic_earth_driver", "canonical",
      "Builds one documented movement-permeability surface plus a hard land/habitat state; it is not suitability or a second colonisation probability."),
    c("hee_dispersal_spherical_kernel", NA,
      "core_process", "dispersal", "canonical",
      "Builds the area-aware, landscape-weighted spherical-normal directional movement kernel on sparse explicit-grid edges."),
    c("hee_dispersal_transition_matrix", NA,
      "core_process", "dispersal", "canonical",
      "Converts a sparse active-movement edge kernel into the target-by-source rate matrix used by global occupancy updates."),
    c("hee_dispersal_location_transition", NA,
      "historical_inference_or_constraint", "phylogeographic_diffusion", "canonical",
      "Converts an Arias-style directional kernel into a source-normalized lineage-location transition with source retention; it deliberately has no colonisation or persistence term."),
    c("hee_dispersal_terrestrial_sink_transition", NA,
      "historical_inference_or_constraint", "phylogeographic_diffusion", "canonical",
      "Sends marine or vanished-land location mass to an explicit geographic sink after active movement and plate carriage; this is not biological extinction."),
    c("hee_dispersal_rare_land_jump_transition", NA,
      "core_process", "dispersal", "canonical",
      "Adds a declared low-rate land-to-land long-distance jump tail without allowing marine ancestral location states; its rate is a sensitivity parameter unless calibrated."),
    c("hee_dispersal_location_distribution_step", NA,
      "historical_inference_or_constraint", "phylogeographic_diffusion", "canonical",
      "Propagates probability mass over explicit geographic states; it is not a grid occupancy update or demographic state model."),
    c("hee_dispersal_compose_plate_transport", NA,
      "historical_inference_or_constraint", "phylogeographic_diffusion", "canonical",
      "Composes horizontal location diffusion with a separately source-normalized plate-carriage operator; it deliberately has no colonisation, persistence, or extinction term."),
    c("hee_dispersal_spacetime_graph", NA,
      "external_driver", "dynamic_earth_driver", "canonical",
      "Builds time-ordered sparse horizontal movement and zero-cost vertical plate-carriage edges without conflating carriage with biological dispersal."),
    c("hee_dispersal_path_diagnostic", NA,
      "core_process", "dispersal", "canonical",
      "Computes one explicit least-cost diagnostic/reconstruction path; it is not an all-pairs dynamic-occupancy kernel."),
    c("hee_dispersal_pruning_likelihood", NA,
      "historical_inference_or_constraint", "historical_constraint", "canonical",
      "Evaluates finite-state phylogeographic pruning likelihoods from explicit transition matrices; broad ranges and extinction require an extended state model."),
    c("hee_dispersal_time_ordered_pruning", NA,
      "historical_inference_or_constraint", "phylogeographic_diffusion", "canonical",
      "Runs time-ordered explicit-grid Felsenstein pruning and streams ancestral location densities without introducing uncalibrated colonisation or persistence rates."),
    c("hee_dispersal_range_support_calibrate", NA,
      "result_or_diagnostic", "dispersal_conditioned_diversity", "canonical",
      "Calibrates a terminal range-support scale from observed grid support; it is explicitly not a demographic occupancy calibration."),
    c("hee_dispersal_range_support_diversity", NA,
      "result_or_diagnostic", "dispersal_conditioned_diversity", "canonical",
      "Summarises terminal-calibrated lineage-support intensity and mixture diversity from location densities without calling them historical richness or occupancy."),
    c("hee_audit_observation_design", NA,
      "observation_uncertainty", "endpoint_observation_design", "canonical",
      "Audits whether matrix zeroes are confirmed absences or occurrence-grid non-detections before endpoint calibration."),
    c("hee_endpoint_presence_only_score", NA,
      "observation_uncertainty", "endpoint_observation_design", "canonical",
      "Scores spatial allocation of an occupancy map at occurrence records without treating non-records as absences."),
    c("hee_endpoint_presence_only_particle_weights", NA,
      "historical_inference_or_constraint", "endpoint_particle_conditioning", "canonical",
      "Reweights complete history particles with independent or spatially held-out presence-only endpoint evidence; training reuse remains diagnostic only."),
    c("hee_posterior_occupancy_diversity", NA,
      "result_or_diagnostic", "posterior_occupancy_diversity", "canonical",
      "Computes sampled-lineage incidence diversity from binary posterior occupancy draws, never from a mean probability field."),
    c("hee_dispersal_propagule_particles", "hee_dispersal_particle_kernel",
      "core_process", "dispersal", "canonical_alias",
      "Samples finite propagules from a pre-defined movement kernel; these are routing particles, not posterior particles or observed gene flow."),
    c("hee_dispersal_kernel", "hee_within_region_movement_kernel",
      "legacy_compatibility", "dispersal", "legacy_compatibility",
      "Legacy within-region least-cost kernel retained for scripts. New analyses should use hee_dispersal_spherical_kernel() or a documented sparse graph."),
    c("hee_dispersal_resistance_surface", NA,
      "core_process", "dispersal", "canonical",
      "Derives a continuous land-only structural resistance surface from raw palaeo-elevation; it is not suitability or a dispersal probability."),
    c("hee_dispersal_connectivity_graph", NA,
      "core_process", "dispersal", "canonical",
      "Builds the directed land-only 4/8-neighbour graph with elevation-step and local-relief effective costs."),
    c("hee_dispersal_edge_kernel", NA,
      "legacy_compatibility", "dispersal", "legacy_compatibility",
      "Legacy exponential local edge kernel retained for reproducibility; optional LDD remains a separately labelled sensitivity component."),
    c("hee_dispersal_rate_kernel", NA,
      "core_process", "dispersal", "canonical",
      "Separates one source-cell movement rate from directional edge allocation so graph degree does not alter total movement intensity."),
    c("hee_dispersal_particle_kernel", "hee_dispersal_propagule_particles",
      "legacy_compatibility", "dispersal", "legacy_alias",
      "Legacy name for finite-propagule routing; it is not a posterior-particle reconstruction or observed gene flow."),
    c("hee_dispersal_arrival_pressure", NA,
      "core_process", "dispersal", "canonical",
      "Aggregates source-cell movement into raw and bounded arrival pressure in optional occupancy scenarios."),
    c("hee_colonisation_arrival_pressure", "hee_arrival_pressure",
      "legacy_compatibility", "legacy_occupancy_transition", "legacy_alias",
      "Legacy alias for an optional dynamic-occupancy arrival term; not a HmscEE core process."),
    c("hee_colonisation_establishment_probability",
      "hee_establishment_probability", "legacy_compatibility", "legacy_occupancy_transition",
      "legacy_alias", "Optional legacy establishment scenario requiring independent calibration; not used by the current Case05."),
    c("hee_colonisation_probability_from_arrival",
      "hee_colonisation_from_arrival", "legacy_compatibility", "legacy_occupancy_transition",
      "legacy_alias", "Optional legacy occupancy transition; not a seventh biological process."),
    c("hee_colonisation_entry_seed", "hee_entry_initialization",
      "legacy_compatibility", "legacy_occupancy_transition", "legacy_alias",
      "Allocates BSM-accepted regional entry events in the legacy nested-occupancy scenario."),
    c("hee_biotic_filtering_effect", NA, "core_process",
      "biotic_filtering", "canonical",
      "Adds an explicit optional biotic effect; default is neutral zero when no external evidence is supplied."),
    c("hee_persistence_probability", NA,
      "legacy_compatibility", "legacy_occupancy_transition", "legacy_compatibility",
      "Optional calibrated or scenario-based local state-retention term; not a core process or current Case05 output."),
    c("hee_projection_global_grid_step", NA, "projection",
      "global_grid_dynamic_occupancy", "canonical",
      "Updates a calibrated or explicitly scenario-labelled global cell occupancy model. It must not be used as the default empirical reconstruction when colonisation and persistence are not independently identifiable."),
    c("hee_evolution_ancestral_response",
      "hee_trait_mediated_ancestral_response", "core_process", "evolution",
      "canonical_alias", "Computes beta = Gamma*T + u while keeping range separate from niche."),
    c("hee_evolution_ancestral_response_direct", NA,
      "core_process", "evolution", "canonical",
      "Directly reconstructs ancestral environmental-response Beta from HMSC tip Beta posterior draws."),
    c("hee_speciation_events", "build_speciation_events",
      "core_process", "speciation", "canonical_alias",
      "Standardizes dated/BGB node or founder speciation events."),
    c("hee_extinction_probability", NA, "core_process", "extinction",
      "canonical", "Scenario/local extinction probability; scope must be reported."),
    c("hee_extinction_layers", NA, "core_process", "extinction",
      "canonical", "Separates local, regional, and lineage extinction proxies."),
    c("hee_earth_state", "hee_paleo_earth_state", "external_driver",
      "dynamic_earth_driver", "canonical_alias",
      "Builds X/H/W-ready palaeo-Earth cell-time state."),
    c("hee_earth_land_age", "hee_land_age", "external_driver",
      "dynamic_earth_driver", "canonical_alias",
      "Computes land age from hard land/habitat state."),
    c("hee_earth_landscape_events", "hee_landscape_events",
      "external_driver", "dynamic_earth_driver", "canonical_alias",
      "Detects land gain/loss and landscape changes."),
    c("hee_history_data", "HmscEcoEvo_data", "historical_information",
      "history_predictors", "canonical_alias",
      "Constructs historical-predictor project data while retaining hmscHist compatibility."),
    c("hee_history_build_tables", "build_history_tables",
      "historical_information", "history_predictors", "canonical_alias",
      "Builds species/site/trait history tables."),
    c("hee_history_standardize", "standardize_history_tables",
      "historical_information", "history_predictors", "canonical_alias",
      "Standardizes historical tables."),
    c("hee_history_indices", "calc_history_indices",
      "historical_information", "history_predictors", "canonical_alias",
      "Calculates historical predictors for modern HMSC interpretation."),
    c("hee_history_as_xdata", "as_hmsc_xdata", "historical_information",
      "history_predictors", "canonical_alias",
      "Converts history predictors to HMSC XData."),
    c("hee_history_as_trdata", "as_hmsc_trdata", "historical_information",
      "history_predictors", "canonical_alias",
      "Converts history predictors to HMSC TrData."),
    c("hee_history_as_random", "as_hmsc_random", "historical_information",
      "history_predictors", "canonical_alias",
      "Builds HMSC random-effect inputs from history tables."),
    c("hee_history_detect_leakage", "hee_detect_leakage",
      "historical_information", "history_predictors", "canonical_alias",
      "Flags current-composition-derived predictors that must not leak into palaeo-grid projection."),
    c("hee_bgb_region_history", "hee_bsm_region_history",
      "historical_inference_or_constraint", "historical_constraint",
      "canonical_alias", "Validates BSM/BioGeoBEARS region states R[h|g]."),
    c("hee_projection_nested_occupancy", "hee_nested_region_cell_occupancy",
      "projection", "integrated_occupancy", "canonical_alias",
      "Runs nested region-cell occupancy and averages g, h|g, c|g and response draws."),
    c("hee_global_grid_dynamic_occupancy", NA,
      "legacy_compatibility", "legacy_occupancy_transition", "legacy_compatibility",
      "Runs the older scenario-calibrated global CTMC occupancy model; current Case05 uses location pruning instead."),
    c("hee_process_attribution", NA,
      "result_diagnostic", "biodiversity_outputs", "canonical",
      "Attributes dynamic occupancy, range, richness, PD, FD, and turnover targets to state transitions, upstream Shapley mechanisms, and lineage events without flat double-counting."),
    c("hee_uncertainty_extrapolation_report",
      "hee_extrapolation_uncertainty_report", "observation_uncertainty",
      "validation_uncertainty", "canonical_alias",
      "Reports no-analog/extrapolation beside p without downweighting biology."),
    c("hee_diversity_region_pool", "hee_region_pool_from_cell_occupancy",
      "result_diagnostic", "biodiversity_outputs", "canonical_alias",
      "Aggregates cell occupancy into region-level pool probability."),
    c("hee_geoprocess_diagnostics", NA, "legacy_compatibility",
      "legacy_geoprocess_observer", "legacy_compatibility",
      "Legacy geological-process diagnostic bundle; use the six-process catalog in new writing."),
    c("hee_dynamic_assembly", NA, "legacy_compatibility",
      "legacy_dynamic_observer", "legacy_compatibility",
      "Legacy scenario/dynamic observer, not the preferred process-name entry point."),
    c("hee_combine_projection_models", NA, "legacy_compatibility",
      "legacy_m1_m5_sensitivity", "legacy_compatibility",
      "Older M1-M5 product-form sensitivity contrast, not the HmscEE 1.0 core.")
  )
  out <- do.call(rbind, lapply(rows, function(z) {
    data.frame(
      canonical_name = z[[1]],
      legacy_name = ifelse(is.na(z[[2]]), NA_character_, z[[2]]),
      layer = z[[3]],
      process_slug = z[[4]],
      status = z[[5]],
      note = z[[6]],
      exported = TRUE,
      stringsAsFactors = FALSE
    )
  }))
  out
}

#' Canonical environmental filtering suitability function
#'
#' Six-process name for [hee_lineage_suitability()].
#'
#' @param ... Arguments passed to [hee_lineage_suitability()].
#' @return Output from [hee_lineage_suitability()].
#' @export
hee_environmental_filtering_suitability <- function(...) {
  hee_lineage_suitability(...)
}

#' Legacy within-region dispersal kernel
#'
#' Compatibility entry point for [hee_within_region_movement_kernel()]. New
#' landscape-explicit analyses should use [hee_dispersal_spherical_kernel()] for
#' the active P2 kernel and [hee_dispersal_spacetime_graph()] to represent plate
#' carriage as a separate Dynamic-Earth operation. This legacy alias remains so
#' that existing scripts and historical results continue to run unchanged.
#'
#' @param ... Arguments passed to [hee_within_region_movement_kernel()].
#' @return Output from [hee_within_region_movement_kernel()].
#' @export
hee_dispersal_kernel <- function(...) {
  hee_within_region_movement_kernel(...)
}

#' Canonical colonisation arrival pressure
#'
#' Legacy optional-occupancy alias for [hee_arrival_pressure()].
#'
#' @param ... Arguments passed to [hee_arrival_pressure()].
#' @return Output from [hee_arrival_pressure()].
#' @export
hee_colonisation_arrival_pressure <- function(...) {
  hee_arrival_pressure(...)
}

#' Canonical establishment probability
#'
#' Legacy optional-occupancy alias for [hee_establishment_probability()]. This is the
#' establishment part of colonisation: arrived propagules still need suitable
#' local conditions to form a population.
#'
#' @param ... Arguments passed to [hee_establishment_probability()].
#' @return Output from [hee_establishment_probability()].
#' @export
hee_colonisation_establishment_probability <- function(...) {
  hee_establishment_probability(...)
}

#' Canonical colonisation probability from arrival and establishment
#'
#' Legacy optional-occupancy alias for [hee_colonisation_from_arrival()].
#'
#' @param ... Arguments passed to [hee_colonisation_from_arrival()].
#' @return Output from [hee_colonisation_from_arrival()].
#' @export
hee_colonisation_probability_from_arrival <- function(...) {
  hee_colonisation_from_arrival(...)
}

#' Canonical entry seeding after BioGeoBEARS regional entry
#'
#' Legacy nested-occupancy alias for [hee_entry_initialization()].
#'
#' @param ... Arguments passed to [hee_entry_initialization()].
#' @return Output from [hee_entry_initialization()].
#' @export
hee_colonisation_entry_seed <- function(...) {
  hee_entry_initialization(...)
}

#' Explicit biotic filtering effect
#'
#' Adds a biotic-filtering term for establishment or persistence models. When
#' no independent interaction/abundance evidence is supplied, the function
#' returns a neutral zero effect and labels the source as not parameterised.
#' This prevents HMSC residual associations from being silently interpreted as
#' causal biotic interactions.
#'
#' @param x Data.frame or numeric vector.
#' @param biotic_col Optional column containing an externally justified biotic
#'   effect. Non-finite values are treated as zero.
#' @param output_col Output column name.
#' @return A data.frame with `output_col` and `biotic_filtering_source`, or a
#'   numeric vector when `x` is numeric.
#' @export
#'
#' @examples
#' hee_biotic_filtering_effect(data.frame(cell_id = "c1"))
hee_biotic_filtering_effect <- function(x,
                                        biotic_col = NULL,
                                        output_col = "biotic_filtering_effect") {
  if (is.data.frame(x)) {
    d <- as.data.frame(x)
    if (is.null(biotic_col)) {
      d[[output_col]] <- 0
      d$biotic_filtering_source <-
        "not_parameterised_no_external_interaction_evidence"
      return(d)
    }
    .require_cols(d, biotic_col, "x")
    effect <- suppressWarnings(as.numeric(d[[biotic_col]]))
    effect[!is.finite(effect) | is.na(effect)] <- 0
    d[[output_col]] <- effect
    d$biotic_filtering_source <- paste0("external_column:", biotic_col)
    return(d)
  }
  effect <- suppressWarnings(as.numeric(x))
  effect[!is.finite(effect) | is.na(effect)] <- 0
  effect
}

#' Canonical evolution ancestral response
#'
#' Six-process name for [hee_trait_mediated_ancestral_response()].
#'
#' @param ... Arguments passed to [hee_trait_mediated_ancestral_response()].
#' @return Output from [hee_trait_mediated_ancestral_response()].
#' @export
hee_evolution_ancestral_response <- function(...) {
  hee_trait_mediated_ancestral_response(...)
}

#' Canonical speciation event standardisation
#'
#' Six-process name for [build_speciation_events()].
#'
#' @param ... Arguments passed to [build_speciation_events()].
#' @return Output from [build_speciation_events()].
#' @export
hee_speciation_events <- function(...) {
  build_speciation_events(...)
}

#' Canonical palaeo-Earth state
#'
#' Preferred driver-layer name for [hee_paleo_earth_state()].
#'
#' @inheritParams hee_paleo_earth_state
#' @return Output from [hee_paleo_earth_state()].
#' @export
hee_earth_state <- function(...) {
  hee_paleo_earth_state(...)
}

#' Canonical land-age driver
#'
#' Preferred driver-layer name for [hee_land_age()].
#'
#' @param ... Arguments passed to [hee_land_age()].
#' @return Output from [hee_land_age()].
#' @export
hee_earth_land_age <- function(...) {
  hee_land_age(...)
}

#' Canonical landscape-event driver
#'
#' Preferred driver-layer name for [hee_landscape_events()].
#'
#' @param ... Arguments passed to [hee_landscape_events()].
#' @return Output from [hee_landscape_events()].
#' @export
hee_earth_landscape_events <- function(...) {
  hee_landscape_events(...)
}

#' Canonical historical-predictor data constructor
#'
#' Preferred history-layer name for [HmscEcoEvo_data()].
#'
#' @param ... Arguments passed to [HmscEcoEvo_data()].
#' @return Output from [HmscEcoEvo_data()].
#' @export
hee_history_data <- function(...) {
  HmscEcoEvo_data(...)
}

#' @rdname hee_history_data
#' @export
hee_history_build_tables <- function(...) {
  build_history_tables(...)
}

#' @rdname hee_history_data
#' @export
hee_history_standardize <- function(...) {
  standardize_history_tables(...)
}

#' @rdname hee_history_data
#' @export
hee_history_indices <- function(...) {
  calc_history_indices(...)
}

#' @rdname hee_history_data
#' @export
hee_history_as_xdata <- function(...) {
  as_hmsc_xdata(...)
}

#' @rdname hee_history_data
#' @export
hee_history_as_trdata <- function(...) {
  as_hmsc_trdata(...)
}

#' @rdname hee_history_data
#' @export
hee_history_as_random <- function(...) {
  as_hmsc_random(...)
}

#' @rdname hee_history_data
#' @export
hee_history_detect_leakage <- function(...) {
  hee_detect_leakage(...)
}

#' Canonical BioGeoBEARS region-history constraint
#'
#' Preferred constraint-layer name for [hee_bsm_region_history()].
#'
#' @param ... Arguments passed to [hee_bsm_region_history()].
#' @return Output from [hee_bsm_region_history()].
#' @export
hee_bgb_region_history <- function(...) {
  hee_bsm_region_history(...)
}

#' Canonical nested projection/occupancy function
#'
#' Preferred projection-layer name for [hee_nested_region_cell_occupancy()].
#'
#' @param ... Arguments passed to [hee_nested_region_cell_occupancy()].
#' @return Output from [hee_nested_region_cell_occupancy()].
#' @export
hee_projection_nested_occupancy <- function(...) {
  hee_nested_region_cell_occupancy(...)
}

#' Canonical extrapolation uncertainty report
#'
#' Preferred uncertainty-layer name for
#' [hee_extrapolation_uncertainty_report()].
#'
#' @param ... Arguments passed to [hee_extrapolation_uncertainty_report()].
#' @return Output from [hee_extrapolation_uncertainty_report()].
#' @export
hee_uncertainty_extrapolation_report <- function(...) {
  hee_extrapolation_uncertainty_report(...)
}

#' Canonical regional diversity pool summary
#'
#' Preferred result-layer name for [hee_region_pool_from_cell_occupancy()].
#'
#' @param ... Arguments passed to [hee_region_pool_from_cell_occupancy()].
#' @return Output from [hee_region_pool_from_cell_occupancy()].
#' @export
hee_diversity_region_pool <- function(...) {
  hee_region_pool_from_cell_occupancy(...)
}
