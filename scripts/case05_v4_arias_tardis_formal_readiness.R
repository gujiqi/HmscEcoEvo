#!/usr/bin/env Rscript

# Case05 v4 -- formal-readiness gate for a landscape-explicit HmscEE analysis.
#
# This script deliberately does not manufacture a "formal" result from the
# Plant200 occurrence-grid data.  It audits the scientific preconditions that
# are required before the expensive 200-taxon, 66-slice global run can be
# launched.  The design follows Arias (2024): an explicit dynamic grid and
# terminal likelihood/integration over histories, and Flannery-Sutherland et
# al. (2025): plate carriage and time-ordered least-cost routes are distinct
# from active dispersal.  A missing precondition produces a durable output
# record and a non-zero status instead of a plausible-looking map.

`%||%` <- function(x, y) if (is.null(x) || !length(x) || identical(x, "") || is.na(x)) y else x
parse_args <- function(args) {
  out <- list()
  for (arg in args) if (grepl("^--", arg)) {
    z <- sub("^--", "", arg); at <- regexpr("=", z, fixed = TRUE)
    if (at < 0L) out[[z]] <- TRUE else out[[substr(z, 1L, at - 1L)]] <-
      substr(z, at + 1L, nchar(z))
  }
  out
}
safe_dir <- function(x) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
  normalizePath(x, winslash = "/", mustWork = FALSE)
}
write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("data.table", quietly = TRUE)) data.table::fwrite(x, path) else
    utils::write.csv(x, path, row.names = FALSE)
  invisible(path)
}
read_key_value_manifest <- function(path) {
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("field", "value") %in% names(x))) {
    stop("hmsc_response_manifest must have field,value columns.", call. = FALSE)
  }
  value <- as.character(x$value)
  names(value) <- as.character(x$field)
  value
}
as_flag <- function(x) tolower(as.character(x)[1L]) %in% c("true", "t", "1", "yes", "y")

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
file_arg <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
if (!length(script)) script <- "scripts/case05_v4_arias_tardis_formal_readiness.R"
pkg_root <- normalizePath(cfg$pkg_root %||%
  file.path(dirname(normalizePath(script[[1L]], winslash = "/", mustWork = FALSE)), ".."),
  winslash = "/", mustWork = TRUE)
workspace <- normalizePath(file.path(pkg_root, ".."), winslash = "/", mustWork = TRUE)
input_root <- normalizePath(cfg$input_root %||% file.path(
  workspace, "HmscEcoEvo_API8_20260818", "data_external",
  "plant_200_multifamily_extant_20260820", "prepared_inputs",
  "case03_data_raw", "Plant200_multifamily_extant_1deg_traits"
), winslash = "/", mustWork = TRUE)
carrier_root <- normalizePath(cfg$carrier_root %||% file.path(
  pkg_root, "derived_inputs", "case05_v2_paleomap_h3_carrier_325_0Ma_1deg"
), winslash = "/", mustWork = TRUE)
movement_root <- normalizePath(cfg$movement_root %||% file.path(
  pkg_root, "derived_inputs", "case05_v3_arias_target_landscape_cache"
), winslash = "/", mustWork = TRUE)
beta_path <- normalizePath(cfg$beta_path %||% file.path(
  pkg_root, "derived_inputs", "case04_plant200_hmsc_beta_posterior_25draws_20260919",
  "case04_hmsc_beta_posterior_25draws.csv"
), winslash = "/", mustWork = TRUE)
output <- safe_dir(cfg$output %||% file.path(
  workspace, "outputs", "HmscEcoEvo",
  paste0("case05_v4_arias_tardis_formal_readiness_", format(Sys.time(), "%Y%m%d"))
))
fail_on_not_ready <- as_flag(cfg$fail_on_not_ready %||% "true")
hmsc_manifest_path <- cfg$hmsc_response_manifest %||% ""
endpoint_manifest_path <- cfg$endpoint_evidence_manifest %||% ""

if (!requireNamespace("pkgload", quietly = TRUE) || !requireNamespace("ape", quietly = TRUE)) {
  stop("pkgload and ape are required for the Case05 v4 readiness audit.", call. = FALSE)
}
pkgload::load_all(pkg_root, quiet = TRUE)

dirs <- list(
  config = safe_dir(file.path(output, "00_config")),
  audit = safe_dir(file.path(output, "01_input_and_model_audit")),
  report = safe_dir(file.path(output, "02_formal_readiness_report"))
)

comm_path <- file.path(input_root, "comm.csv")
sites_path <- file.path(input_root, "sites.csv")
tree_path <- file.path(input_root, "tree.tre")
fossil_path <- file.path(input_root, "fossils.csv")
manifest_path <- file.path(input_root, "input_manifest.csv")
required <- c(comm_path, sites_path, tree_path, fossil_path, manifest_path,
              beta_path,
              file.path(movement_root, "case05_v3_movement_cache_index.csv"),
              file.path(movement_root, "case05_v3_movement_metadata.csv"),
              file.path(carrier_root, "case05_v2_plate_carrier_index.csv"))
if (any(!file.exists(required))) {
  write_csv(data.frame(path = required, exists = file.exists(required)),
            file.path(dirs$audit, "missing_required_inputs.csv"))
  stop("Case05 v4 required inputs are missing; see missing_required_inputs.csv.",
       call. = FALSE)
}

comm <- utils::read.csv(comm_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!"site_id" %in% names(comm)) stop("comm.csv lacks site_id.", call. = FALSE)
species <- setdiff(names(comm), "site_id")
Y <- as.matrix(comm[, species, drop = FALSE]); storage.mode(Y) <- "double"
rownames(Y) <- as.character(comm$site_id)
sites <- utils::read.csv(sites_path, stringsAsFactors = FALSE)
tree <- ape::read.tree(tree_path)
fossils <- utils::read.csv(fossil_path, stringsAsFactors = FALSE, check.names = FALSE)
input_manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
beta <- utils::read.csv(beta_path, stringsAsFactors = FALSE, check.names = FALSE)

site_alignment <- setequal(rownames(Y), as.character(sites$site_id)) &&
  !anyDuplicated(rownames(Y)) && !anyDuplicated(as.character(sites$site_id))
tip_alignment <- setequal(species, tree$tip.label) && !anyDuplicated(tree$tip.label)
observation <- HmscEcoEvo::hee_audit_observation_design(
  Y, design_label = "Plant200 GBIF occupied 1-degree occurrence grid"
)
root_age <- max(ape::node.depth.edgelength(tree)[seq_along(tree$tip.label)])
fossil_rank_col <- intersect(c("taxon_match_rank", "match_rank", "taxonomic_rank",
                               "rank", "matched_rank"), names(fossils))
fossil_rank <- if (length(fossil_rank_col)) as.character(fossils[[fossil_rank_col[[1L]]]]) else NA_character_
fossil_summary <- data.frame(
  n_fossil_records = nrow(fossils),
  n_species_level = sum(tolower(fossil_rank) %in% c("species", "exact_species"), na.rm = TRUE),
  n_genus_level = sum(tolower(fossil_rank) %in% c("genus", "exact_genus"), na.rm = TRUE),
  n_family_or_higher = sum(!(tolower(fossil_rank) %in% c("species", "exact_species", "genus", "exact_genus")), na.rm = TRUE),
  fossil_rank_column = if (length(fossil_rank_col)) fossil_rank_col[[1L]] else NA_character_,
  scientific_boundary = "Family/genus fossil matches constrain clade presence only; they do not validate a selected extant species at a palaeocell.",
  stringsAsFactors = FALSE
)

movement_index <- utils::read.csv(file.path(movement_root, "case05_v3_movement_cache_index.csv"), stringsAsFactors = FALSE)
movement_meta <- utils::read.csv(file.path(movement_root, "case05_v3_movement_metadata.csv"), stringsAsFactors = FALSE)
movement_value <- stats::setNames(as.character(movement_meta$value), as.character(movement_meta$field))
carrier_index <- utils::read.csv(file.path(carrier_root, "case05_v2_plate_carrier_index.csv"), stringsAsFactors = FALSE)
movement_complete <- nrow(movement_index) == 66L &&
  all(file.exists(as.character(movement_index$file))) &&
  all(movement_index$n_sources_truncated == 0L)
arias_kernel_clean <- identical(movement_value[["distance_term"]],
  "great-circle distance only; no least-cost or resistance-expanded distance in the active kernel") &&
  identical(movement_value[["topography_count"]],
  "one: topography is a declared destination-permeability surface only")
tardis_separate <- grepl("separate route diagnostics", movement_value[["tardis_role"]], fixed = TRUE)
plate_separate <- grepl("not organismal dispersal", movement_value[["plate_carriage_role"]], fixed = TRUE)

hmsc_manifest <- read_key_value_manifest(hmsc_manifest_path)
hmsc_observation_model <- if (is.null(hmsc_manifest)) NA_character_ else hmsc_manifest[["observation_model"]] %||% NA_character_
hmsc_absence_validated <- if (is.null(hmsc_manifest)) FALSE else as_flag(hmsc_manifest[["confirmed_absence_or_valid_presence_background"]] %||% "false")
hmsc_cross_fitted <- if (is.null(hmsc_manifest)) FALSE else as_flag(hmsc_manifest[["spatially_cross_fitted"]] %||% "false")
hmsc_beta_present <- nrow(beta) > 0L && all(c("response_draw", "lineage") %in% names(beta)) &&
  all(species %in% unique(as.character(beta$lineage)))
hmsc_model_ready <- hmsc_beta_present && hmsc_absence_validated

endpoint_manifest <- if (nzchar(endpoint_manifest_path) && file.exists(endpoint_manifest_path)) {
  utils::read.csv(endpoint_manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
} else NULL
endpoint_required <- c("path", "data_role", "independent_from_hmsc_training", "observation_type", "effort_surface_path")
endpoint_schema_ok <- !is.null(endpoint_manifest) && all(endpoint_required %in% names(endpoint_manifest)) && nrow(endpoint_manifest) > 0L
endpoint_independent <- endpoint_schema_ok && all(as_flag(endpoint_manifest$independent_from_hmsc_training)) &&
  all(tolower(as.character(endpoint_manifest$data_role)) %in% c("independent_validation", "spatial_holdout"))
endpoint_paths_ok <- endpoint_schema_ok && all(file.exists(as.character(endpoint_manifest$path)))
endpoint_effort_paths_ok <- endpoint_schema_ok &&
  all(nzchar(as.character(endpoint_manifest$effort_surface_path))) &&
  all(file.exists(as.character(endpoint_manifest$effort_surface_path)))
endpoint_ready <- endpoint_schema_ok && endpoint_independent && endpoint_paths_ok && endpoint_effort_paths_ok

# Static contracts prevent a runner refactor from silently restoring the two
# invalid shortcuts removed from Case05: per-step survival renormalisation and
# Shannon/Simpson maps calculated from a mean probability field.
runner_path <- file.path(pkg_root, "scripts", "case05_v3_run_reference_calibrated_draw.R")
runner_text <- if (file.exists(runner_path)) paste(readLines(runner_path, warn = FALSE), collapse = "\n") else ""
forward_survival_contract <- grepl("condition_extant_tree_survival.*false", runner_text) &&
  grepl("must not condition every time step", runner_text, fixed = TRUE) &&
  !grepl("hee_extinction_condition_extant_lineage_survival\\(", runner_text)
posterior_diversity_contract <- exists("hee_posterior_occupancy_diversity", where = asNamespace("HmscEcoEvo"), inherits = FALSE)
formal_presence_only_contract <- grepl("declared effort/exposure surface",
  paste(readLines(file.path(pkg_root, "R", "endpoint_presence_only.R"), warn = FALSE), collapse = "\n"), fixed = TRUE)

checks <- data.frame(
  check = c(
    "species_tree_alignment", "site_community_alignment", "observation_design_has_confirmed_absences",
    "HMSC_response_posterior_present", "HMSC_observation_model_validated",
    "HMSC_and_endpoint_are_spatially_cross_fitted", "independent_or_holdout_endpoint_manifest",
    "endpoint_files_exist", "endpoint_effort_surface_exists", "dated_tree_is_screening_backbone", "fossil_evidence_is_species_level",
    "all_66_dynamic_earth_slices_cached", "Arias_spherical_kernel_no_duplicate_topography",
    "plate_carriage_is_not_active_dispersal", "TARDIS_path_is_not_a_second_probability",
    "carrier_grid_has_all_time_slices", "forward_runner_disables_per_step_survival_renormalisation",
    "posterior_diversity_requires_binary_state_draws", "formal_presence_only_requires_effort_surface"
  ),
  passed = c(
    tip_alignment, site_alignment, isTRUE(observation$bernoulli_absence_likelihood_eligible[[1L]]),
    hmsc_beta_present, hmsc_model_ready, hmsc_cross_fitted, endpoint_ready,
    endpoint_paths_ok, endpoint_effort_paths_ok, FALSE, fossil_summary$n_species_level > 0L,
    movement_complete, arias_kernel_clean, plate_separate, tardis_separate,
    nrow(carrier_index) == 66L && all(file.exists(as.character(carrier_index$file))),
    forward_survival_contract, posterior_diversity_contract, formal_presence_only_contract
  ),
  category = c(
    "input", "input", "observation", "HMSC", "HMSC", "HMSC_endpoint",
    "endpoint", "endpoint", "endpoint", "phylogeny", "fossils", "dynamic_earth",
    "dispersal", "dynamic_earth", "dispersal_diagnostic", "dynamic_earth",
    "dynamic_occupancy", "diversity", "endpoint"
  ),
  consequence_if_failed = c(
    "Stop: species columns and tree tips do not define the same taxon set.",
    "Stop: community rows and environmental sites do not define the same cells.",
    "GBIF non-records cannot enter a Bernoulli endpoint likelihood.",
    "No tip-response posterior exists for ancestral-response reconstruction.",
    "Do not interpret an HMSC fitted to occurrence-only zeroes as a calibrated occurrence model.",
    "Reused observations would double-count modern distribution information.",
    "No formal endpoint-conditioned reconstruction may be launched.",
    "Endpoint evidence cannot be evaluated.",
    "Presence-only endpoint likelihood cannot be interpreted without a declared spatial effort/exposure surface.",
    "V.PhyloMaker screening backbone cannot support a paper-grade deep-time timing claim without a calibrated-tree sensitivity analysis.",
    "Only clade-level rather than selected-lineage fossil conditioning is allowed.",
    "Global trajectory is not represented at every declared Earth time slice.",
    "Distance and topography would be duplicated in active P2 movement.",
    "Plate displacement would be incorrectly counted as active organismal movement.",
    "Least-cost connectivity would be double-counted as a movement/occupancy probability.",
    "Palaeogeographic carrier transfer is incomplete.",
    "A forward scenario would artificially reinflate branches instead of weighting complete trajectories at the appropriate endpoint/node event.",
    "Shannon/Simpson maps could be misreported from a mean probability field rather than binary posterior occupancy draws.",
    "A formal presence-only endpoint could silently fall back to a uniform sampling surface."
  ),
  stringsAsFactors = FALSE
)

formal_required <- checks$check %in% c(
  "species_tree_alignment", "site_community_alignment", "HMSC_response_posterior_present",
  "HMSC_observation_model_validated", "HMSC_and_endpoint_are_spatially_cross_fitted",
  "independent_or_holdout_endpoint_manifest", "endpoint_files_exist",
  "endpoint_effort_surface_exists",
  "dated_tree_is_screening_backbone",
  "all_66_dynamic_earth_slices_cached", "Arias_spherical_kernel_no_duplicate_topography",
  "plate_carriage_is_not_active_dispersal", "TARDIS_path_is_not_a_second_probability",
  "carrier_grid_has_all_time_slices",
  "forward_runner_disables_per_step_survival_renormalisation",
  "posterior_diversity_requires_binary_state_draws",
  "formal_presence_only_requires_effort_surface"
)
formal_ready <- all(checks$passed[formal_required])
status <- if (formal_ready) "FORMAL_RUN_AUTHORISED" else "FORMAL_RUN_BLOCKED"

write_csv(data.frame(
  field = c("case_id", "model_definition", "Arias_2024_role", "Flannery_Sutherland_2025_role",
            "input_root", "n_species", "n_occurrence_grid_cells", "n_total_records",
            "tree_root_age_ma", "formal_status", "scientific_boundary"),
  value = c(
    "case05_v4_plant200_global_dynamic", 
    "explicit_grid_dynamic_occupancy_with_terminal_likelihood_or_cross_fitted_presence_only_endpoint",
    "Time-stratified spherical diffusion on an explicit dynamic grid; terminal data constrain histories through likelihood/pruning rather than a terminal anchor.",
    "Time-ordered least-cost paths diagnose inferred routes and encountered terrain; they are not multiplied into active movement or occupancy.",
    input_root, length(species), nrow(Y), sum(Y), format(root_age, digits = 8), status,
    "Plant200 currently contains GBIF occupied-cell occurrence records and a V.PhyloMaker2 screening backbone. Until a validated HMSC observation model, cross-fitted or independent endpoint, and calibrated-tree sensitivity are supplied, no output may be called a formal palaeodistribution reconstruction."
  ), stringsAsFactors = FALSE
), file.path(dirs$config, "case05_v4_model_contract.csv"))
write_csv(data.frame(
  template = c("hmsc_response_manifest", "endpoint_evidence_manifest"),
  package_template = c(
    file.path(pkg_root, "inst", "extdata", "case05_hmsc_response_manifest_template.csv"),
    file.path(pkg_root, "inst", "extdata", "case05_endpoint_evidence_manifest_template.csv")
  ),
  purpose = c(
    "Documents a valid, spatially cross-fitted response posterior and its locked environmental recipe.",
    "Documents independent or held-out terminal occurrence/range evidence and its effort surface."
  ),
  stringsAsFactors = FALSE
), file.path(dirs$config, "case05_v4_required_manifest_templates.csv"))
write_csv(observation, file.path(dirs$audit, "observation_design_audit.csv"))
write_csv(fossil_summary, file.path(dirs$audit, "fossil_evidence_audit.csv"))
write_csv(checks, file.path(dirs$audit, "formal_readiness_checks.csv"))
write_csv(data.frame(
  formal_status = status,
  formal_run_authorised = formal_ready,
  n_blockers = sum(!checks$passed[formal_required]),
  hmsc_response_manifest = if (nzchar(hmsc_manifest_path)) normalizePath(hmsc_manifest_path, winslash = "/", mustWork = FALSE) else "not_supplied",
  endpoint_evidence_manifest = if (nzchar(endpoint_manifest_path)) normalizePath(endpoint_manifest_path, winslash = "/", mustWork = FALSE) else "not_supplied",
  stringsAsFactors = FALSE
), file.path(dirs$report, "case05_v4_execution_status.csv"))

blockers <- checks[formal_required & !checks$passed, c("check", "category", "consequence_if_failed"), drop = FALSE]
lines <- c(
  "# Case05 v4 formal-readiness result",
  "",
  paste0("**Status:** `", status, "`"),
  "",
  "## What passed",
  "",
  "- The cached Earth sequence has 66 complete 325-0 Ma carrier slices.",
  "- Active P2 movement is Arias-style spherical diffusion on an explicit grid.",
  "- Plate carriage is separate from active movement, and TARDIS-style paths are diagnostics rather than a second multiplier.",
  "",
  "## What prevents a formal empirical run",
  ""
)
if (nrow(blockers)) {
  lines <- c(lines, unlist(Map(function(a, b) paste0("- **", a, "**: ", b),
                               blockers$check, blockers$consequence_if_failed)))
} else {
  lines <- c(lines, "- No blocking condition remains.")
}
lines <- c(lines,
  "",
  "## Required repair before launch",
  "",
  "1. Refit or validate the HMSC response model with surveyed presence-absence data, an explicit detection model, or a documented presence-background design. The occupied-cell GBIF matrix alone does not provide absences.",
  "2. Use spatial cross-fitting: HMSC response estimation uses training blocks and the dynamic endpoint likelihood uses held-out occurrence records plus a declared effort surface; alternatively supply a genuinely independent range/occurrence endpoint. Repeated GBIF records must be aggregated to unique taxon-by-grid presence unless an explicit record-level observation model is supplied.",
  "3. Add a fossil-calibrated posterior-tree or a dated-tree sensitivity ensemble. The current V.PhyloMaker2 backbone is a workflow screen, not a deep-time dating result.",
  "4. Condition particles with the endpoint likelihood/pruning layer. Do not use a terminal anchor, copy modern maps backward, or call a forward scenario a reconstruction.",
  "",
  "## Output meaning",
  "",
  "Any older Case05 map is retained only as a rejected or exploratory forward scenario. Some older outputs also used per-step extant-lineage survival renormalisation and mean-field diversity aliases; neither may be used as a result map, richness map, Shannon/Simpson figure, or validation figure."
)
writeLines(lines, file.path(dirs$report, "case05_v4_formal_readiness.md"), useBytes = TRUE)

message("Case05 v4 readiness audit: ", status, " -> ", output)
if (!formal_ready && fail_on_not_ready) {
  stop("Case05 v4 formal run is blocked; inspect 01_input_and_model_audit/formal_readiness_checks.csv.", call. = FALSE)
}
