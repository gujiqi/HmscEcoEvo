#!/usr/bin/env Rscript

# Run one reference-calibrated HMSC response draw through Case05 v3.
#
# Case05 v3 retains the dynamic palaeogeographic carrier state, Arias-style
# spherical P2 kernel, and no-BioGeoBEARS global domain of v2. It changes P1:
# beta-only environmental response is centred on an observed-modern
# X-beta reference reconstructed along the dated tree. This prevents a
# zero-intercept probit transformation from being misreported as occurrence
# richness. The run remains a process scenario unless a separate endpoint
# likelihood is supplied; it never applies a terminal anchor.

`%||%` <- function(x, y) if (is.null(x) || !length(x) || identical(x, "") || is.na(x)) y else x
parse_args <- function(args) {
  out <- list()
  for (arg in args) if (grepl("^--", arg)) {
    z <- sub("^--", "", arg); at <- regexpr("=", z, fixed = TRUE)
    if (at < 0L) out[[z]] <- TRUE else out[[substr(z, 1L, at - 1L)]] <- substr(z, at + 1L, nchar(z))
  }
  out
}
as_num <- function(x, default) { z <- suppressWarnings(as.numeric(x)[1L]); if (is.finite(z)) z else default }
as_int <- function(x, default) { z <- as_num(x, default); if (is.finite(z)) as.integer(z) else default }
safe_dir <- function(x) { if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE); normalizePath(x, winslash = "/", mustWork = FALSE) }
write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(x, path)
  } else {
    utils::write.csv(x, path, row.names = FALSE)
  }
  invisible(path)
}
time_slug <- function(x) formatC(as.numeric(x), format = "f", digits = 6, drop0trailing = TRUE)

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
file_arg <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
if (!length(script)) script <- "scripts/case05_v3_run_reference_calibrated_draw.R"
pkg_root <- normalizePath(cfg$pkg_root %||% file.path(dirname(normalizePath(script[[1L]], winslash = "/", mustWork = FALSE)), ".."), winslash = "/", mustWork = TRUE)
workspace <- normalizePath(file.path(pkg_root, ".."), winslash = "/", mustWork = TRUE)
carrier_root <- cfg$carrier_root %||% file.path(pkg_root, "derived_inputs", "case05_v2_paleomap_h3_carrier_325_0Ma_1deg")
movement_root <- cfg$movement_root %||% file.path(
  pkg_root, "derived_inputs", "case05_v3_arias_target_landscape_cache"
)
beta_path <- cfg$beta_path %||% file.path(pkg_root, "derived_inputs", "case04_plant200_hmsc_beta_posterior_25draws_20260919", "case04_hmsc_beta_posterior_25draws.csv")
reference_path <- cfg$reference_path %||% file.path(
  pkg_root, "derived_inputs", "case05_v3_reference_environmental_support",
  "ancestral_response_reference_thresholds.csv"
)
plant_input_dir <- cfg$plant_input_dir %||% file.path(workspace, "HmscEcoEvo_API8_20260818", "data_external", "plant_200_multifamily_extant_20260820", "prepared_inputs", "case03_data_raw", "Plant200_multifamily_extant_1deg_traits")
tree_path <- cfg$tree %||% file.path(plant_input_dir, "tree.tre")
draw_id <- cfg$response_draw %||% ""
output <- safe_dir(cfg$output %||% stop("--output is required.", call. = FALSE))
internal_dt_myr <- as_num(cfg$internal_dt_myr, 1)
root_max_carriers <- as_int(cfg$root_max_carriers, 120L)
root_occupancy <- as_num(cfg$root_occupancy, 0.8)
root_prior_mode <- match.arg(
  tolower(cfg$root_prior_mode %||% "environment_weighted"),
  c("environment_weighted", "max_environment", "uniform_land")
)
root_support_power <- as_num(cfg$root_support_power, 1)
establishment_intercept <- as_num(cfg$establishment_intercept, -0.25)
establishment_slope <- as_num(cfg$establishment_slope, 0.85)
persistence_intercept <- as_num(cfg$persistence_intercept, 1.25)
persistence_slope <- as_num(cfg$persistence_slope, 0.75)
movement_rate_multiplier <- as_num(cfg$movement_rate_multiplier, 1)
diversity_min_expected_richness <- as_num(cfg$diversity_min_expected_richness, 0.05)
save_all_carrier_states <- tolower(cfg$save_all_carrier_states %||% "true") %in% c("true", "t", "1", "yes")
condition_extant_tree_survival <- tolower(cfg$condition_extant_tree_survival %||% "false") %in% c("true", "t", "1", "yes")
seed <- as_int(cfg$seed, 20260922L)

if (!nzchar(draw_id) || !file.exists(beta_path) || !file.exists(tree_path) ||
    !file.exists(reference_path)) {
  stop("response_draw, beta_path, reference_path and tree are required and must exist.", call. = FALSE)
}
if (!requireNamespace("pkgload", quietly = TRUE) || !requireNamespace("ape", quietly = TRUE)) {
  stop("pkgload and ape are required.", call. = FALSE)
}
if (!is.finite(internal_dt_myr) || internal_dt_myr <= 0 || root_max_carriers < 1L ||
    !is.finite(root_occupancy) || root_occupancy <= 0 || root_occupancy > 1 ||
    !is.finite(root_support_power) || root_support_power < 0 ||
    !is.finite(movement_rate_multiplier) || movement_rate_multiplier < 0 ||
    diversity_min_expected_richness < 0) {
  stop("Invalid numerical process settings.", call. = FALSE)
}
if (condition_extant_tree_survival) {
  stop(
    "A forward Case05 scenario must not condition every time step on extant-tree survival. ",
    "That operation reinflates low-probability states and is valid only after ",
    "complete particles are weighted by a declared endpoint/node event. Use ",
    "--condition_extant_tree_survival=false; formal smoothing requires the v4 ",
    "particle workflow and independent or held-out endpoint evidence.",
    call. = FALSE
  )
}
pkgload::load_all(pkg_root, quiet = TRUE)
set.seed(seed)

carrier_index_path <- cfg$carrier_index %||% file.path(carrier_root, "case05_v2_plate_carrier_index.csv")
movement_index_path <- cfg$movement_index %||% file.path(movement_root, "case05_v3_movement_cache_index.csv")
if (!file.exists(carrier_index_path) || !file.exists(movement_index_path)) {
  stop("Run carrier and movement-cache preparation before a draw.", call. = FALSE)
}
carrier_index <- utils::read.csv(carrier_index_path, stringsAsFactors = FALSE)
movement_index <- utils::read.csv(movement_index_path, stringsAsFactors = FALSE)
carrier_index$time_ma <- suppressWarnings(as.numeric(carrier_index$time_ma))
movement_index$time_ma <- suppressWarnings(as.numeric(movement_index$time_ma))
carrier_index <- carrier_index[order(carrier_index$time_ma, decreasing = TRUE), , drop = FALSE]
if (!all(carrier_index$time_ma %in% movement_index$time_ma) || any(!file.exists(carrier_index$file))) {
  stop("Carrier and movement indexes must cover the same valid time slices.", call. = FALSE)
}
# A PALEOMAP track has stable identity but can be absent from an individual
# reconstructed slice.  Keep the union as the latent state index; a missing
# coordinate is a hard geographic absence for that slice, not an instruction
# to recycle a different target cell or to reorder state by row position.
carrier_universe_ids <- unique(unlist(lapply(carrier_index$file, function(path) {
  as.character(readRDS(path)$carriers$track_id)
}), use.names = FALSE))
if (!length(carrier_universe_ids) || any(!nzchar(carrier_universe_ids)) ||
    anyDuplicated(carrier_universe_ids)) {
  stop("Carrier slices must define a non-empty, unique stable track universe.",
       call. = FALSE)
}
env_times <- carrier_index$time_ma
tree <- ape::read.tree(tree_path)
if (!ape::is.ultrametric(tree) || is.null(tree$edge.length)) {
  stop("Case05 v3 requires an ultrametric dated tree with branch lengths in Ma.", call. = FALSE)
}
beta <- utils::read.csv(beta_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("lineage", "response_draw", "intercept") %in% names(beta))) {
  stop("Beta posterior table must contain lineage, response_draw and intercept.", call. = FALSE)
}
beta <- beta[as.character(beta$response_draw) == draw_id, , drop = FALSE]
if (nrow(beta) != length(tree$tip.label) || !setequal(beta$lineage, tree$tip.label)) {
  stop("Selected HMSC response draw must contain exactly one Beta row per dated-tree tip.", call. = FALSE)
}
basis_cols <- setdiff(names(beta)[vapply(beta, is.numeric, logical(1))],
                       c("chain", "iteration", "response_weight", "intercept", "hmsc_intercept_original"))
if (!length(basis_cols)) stop("No environmental Beta axes found in beta_path.", call. = FALSE)

# Node events are inserted exactly; Earth state remains piecewise constant
# within the enclosing environmental interval, rather than inventing a new
# palaeogeographic raster at a fractional node age.
n_tip <- length(tree$tip.label)
depth <- ape::node.depth.edgelength(tree)
root_age <- max(depth[seq_len(n_tip)])
node_age <- root_age - depth
root_node <- setdiff(tree$edge[, 1L], tree$edge[, 2L])[[1L]]
node_labels <- tree$node.label
if (is.null(node_labels)) node_labels <- rep("", tree$Nnode)
lineage_name <- function(node) {
  vapply(as.integer(node), function(one_node) {
    if (one_node <= n_tip) return(tree$tip.label[[one_node]])
    label <- node_labels[[one_node - n_tip]]
    if (!is.na(label) && nzchar(label)) label else paste0("node_", one_node)
  }, character(1))
}
node_times <- node_age[(n_tip + 1L):(n_tip + tree$Nnode)]
node_times <- node_times[node_times <= max(env_times) + 1e-8 & node_times >= min(env_times) - 1e-8]
master <- HmscEcoEvo::hee_master_time_axis(env_times = env_times,
                                             tree = tree, tolerance = 1e-6)
# At a speciation age the response helper returns the incoming (parent)
# branch by convention.  Add a tiny, explicitly recorded younger-side query
# for each node, so a state that has just been split can retrieve its daughter
# response without changing the Earth slice or inventing a new time raster.
min_positive_branch_myr <- min(tree$edge.length[tree$edge.length > 0])
event_response_epsilon_myr <- min(1e-4, min_positive_branch_myr / 20)
post_node_response_times <- node_times - event_response_epsilon_myr
post_node_response_times <- post_node_response_times[
  post_node_response_times >= min(env_times) - 1e-8
]
response_times <- sort(unique(c(
  master$time_ma[master$time_ma <= max(env_times) + 1e-8],
  post_node_response_times
)), decreasing = TRUE)
tip_beta <- beta[, c("lineage", "response_draw", "intercept", basis_cols), drop = FALSE]
names(tip_beta)[names(tip_beta) == "lineage"] <- "species"
responses <- HmscEcoEvo::hee_evolution_ancestral_response_direct(
  tip_beta_draws = tip_beta, tree = tree, times = response_times,
  basis_cols = basis_cols, reconstruct_intercept = FALSE,
  branch_uncertainty = "none"
)
references <- utils::read.csv(reference_path, check.names = FALSE,
                              stringsAsFactors = FALSE)
if (!all(c("lineage", "response_draw", "time_ma", "eta_reference", "eta_scale") %in%
         names(references))) {
  stop("reference_path must contain lineage, response_draw, time_ma, eta_reference, and eta_scale.",
       call. = FALSE)
}
references <- references[as.character(references$response_draw) == draw_id,
                         c("lineage", "time_ma", "eta_reference", "eta_scale"), drop = FALSE]
if (!nrow(references) || any(!is.finite(references$eta_reference)) ||
    any(!is.finite(references$eta_scale) | references$eta_scale <= 1e-8)) {
  stop("reference_path has no finite thresholds for response draw ", draw_id,
       call. = FALSE)
}

config_dir <- safe_dir(file.path(output, "00_config"))
state_dir <- safe_dir(file.path(output, "16_state_space_inference", "carrier_states"))
summary_dir <- safe_dir(file.path(output, "16_state_space_inference", "carrier_metrics"))
process_dir <- safe_dir(file.path(output, "16_state_space_inference", "carrier_interval_process"))
gradient_dir <- safe_dir(file.path(output, "09_environmental_filtering", "reference_contrast_audit"))
diversity_status_dir <- safe_dir(file.path(output, "18_diversity_maps"))
write_csv(data.frame(
  status = "NOT_RUN_REQUIRES_BINARY_POSTERIOR_OCCUPANCY_DRAWS",
  requested_metric = c("posterior_incidence_shannon_entropy", "posterior_incidence_gini_simpson", "posterior_phylogenetic_diversity", "posterior_functional_diversity"),
  reason = "This v3 runner propagates deterministic mean-field occupancy matrices per HMSC response draw, not binary dynamic-occupancy particles. Use hee_posterior_occupancy_diversity() after a particle/smoothing workflow; do not relabel mean-field composition summaries as posterior diversity.",
  stringsAsFactors = FALSE
), file.path(diversity_status_dir, "posterior_diversity_status.csv"))
write_csv(data.frame(
  parameter = c("response_draw", "root_age_ma", "earth_time_representation", "plate_transport", "root_prior_mode", "root_support_power", "environmental_reference", "movement_rate_multiplier", "terminal_anchor", "process_mode", "extant_tree_survival_conditioning", "diversity_min_expected_richness", "scientific_boundary"),
  value = c(draw_id, root_age,
            "piecewise_constant_with_exact_tree_node_splitting",
            "stable_PALEOMAP_H3_carrier_identity_not_target_centred_raster_copy",
            root_prior_mode, root_support_power,
            "modern_X_beta_reference_and_robust_scale_reconstructed_on_dated_tree; see reference_support_metadata.csv",
            movement_rate_multiplier, "none",
            "reference_calibrated_forward_mean_field_process_scenario",
            "disabled_in_forward_scenario", diversity_min_expected_richness,
            "No independent endpoint likelihood is applied here. P1 uses a modern-reference X-beta contrast with a reconstructed robust response scale; it is a relative environmental-support covariate, not an occurrence probability. Extant-tree survival is not reimposed during the forward recursion: a lineage with zero global support is recorded as incompatible for later particle or endpoint weighting, not artificially reinflated. The draw remains a forward process scenario, not a posterior palaeodistribution reconstruction."),
  stringsAsFactors = FALSE
), file.path(config_dir, "draw_configuration.csv"))
write_csv(master, file.path(config_dir, "master_time_axis.csv"))
write_csv(responses, file.path(output, "04_evolution", "ancestral_beta_response_master_time.csv"))

load_slice <- function(time_ma) {
  cidx <- match(time_ma, carrier_index$time_ma)
  midx <- match(time_ma, movement_index$time_ma)
  carrier_blob <- readRDS(carrier_index$file[[cidx]])
  movement_blob <- readRDS(movement_index$file[[midx]])
  carrier <- as.data.frame(movement_blob$carrier)
  if (!identical(as.character(carrier$track_id), as.character(carrier_blob$carriers$track_id))) {
    stop("Carrier ordering differs between carrier and movement cache at ", time_ma, " Ma.")
  }
  if (anyDuplicated(carrier$track_id) || any(!carrier$track_id %in% carrier_universe_ids)) {
    stop("Carrier slice contains invalid stable track IDs at ", time_ma, " Ma.")
  }
  list(carrier = carrier, movement = movement_blob$movement,
       map_coverage = carrier_blob$output_coverage,
       output_coverage = carrier_blob$diagnostics$output_coverage_fraction)
}
first_slice <- load_slice(env_times[[1L]])
carrier_ids <- as.character(first_slice$carrier$track_id)
if (anyDuplicated(carrier_ids)) stop("Carrier IDs must be unique.")
# The dated root is 0.05 Myr older than the oldest available Earth slice.
# At the first Earth slice the root has already divided into its two child
# branches, so there is no physical ``root lineage'' column to simulate. Both
# daughter branches inherit one compact ancestral patch; this is the explicit
# root-boundary convention rather than a fabricated pre-root environment.
root_children <- tree$edge[tree$edge[, 1L] == root_node, 2L]
if (length(root_children) < 2L) {
  stop("The dated tree root must have at least two child branches.", call. = FALSE)
}
# Retain the first-slice ordering at the top so that root patch indices stay
# explicit; append tracks that only occur in younger Earth reconstructions.
carrier_universe_ids <- unique(c(carrier_ids, carrier_universe_ids))
q <- matrix(0, nrow = length(carrier_universe_ids), ncol = length(root_children),
            dimnames = list(carrier_universe_ids, lineage_name(root_children)))

response_for <- function(time_ma, lineages) {
  tol <- 1e-8
  out <- responses[abs(responses$time_ma - time_ma) <= tol &
                     responses$lineage %in% lineages, , drop = FALSE]
  missing <- setdiff(lineages, out$lineage)
  if (length(missing)) {
    post_node_time <- time_ma - event_response_epsilon_myr
    post <- responses[abs(responses$time_ma - post_node_time) <= tol &
                        responses$lineage %in% missing, , drop = FALSE]
    out <- rbind(out, post)
  }
  if (nrow(out) != length(lineages) || anyDuplicated(out$lineage)) {
    stop("Missing or duplicate ancestral Beta response at ", time_ma,
         " Ma for: ", paste(setdiff(lineages, out$lineage), collapse = ", "), call. = FALSE)
  }
  out[match(lineages, out$lineage), , drop = FALSE]
}

reference_for <- function(time_ma, lineages) {
  tol <- 1e-8
  out <- references[abs(references$time_ma - time_ma) <= tol &
                      references$lineage %in% lineages, , drop = FALSE]
  missing <- setdiff(lineages, out$lineage)
  if (length(missing)) {
    post_node_time <- time_ma - event_response_epsilon_myr
    post <- references[
      abs(references$time_ma - post_node_time) <= tol &
        references$lineage %in% missing, , drop = FALSE
    ]
    out <- rbind(out, post)
  }
  if (nrow(out) != length(lineages) || anyDuplicated(out$lineage) ||
      any(!is.finite(out$eta_reference)) ||
      any(!is.finite(out$eta_scale) | out$eta_scale <= 1e-8)) {
    stop("Missing or duplicate environmental reference at ", time_ma,
         " Ma for: ", paste(setdiff(lineages, out$lineage), collapse = ", "),
         call. = FALSE)
  }
  out[match(lineages, out$lineage), , drop = FALSE]
}

eta_for <- function(slice, time_ma, lineages) {
  carrier <- slice$carrier
  if (!all(basis_cols %in% names(carrier))) {
    stop("Carrier environmental table lacks Beta axis: ",
         paste(setdiff(basis_cols, names(carrier)), collapse = ", "), call. = FALSE)
  }
  resp <- response_for(time_ma, lineages)
  X <- as.matrix(carrier[, basis_cols, drop = FALSE])
  B <- as.matrix(resp[, basis_cols, drop = FALSE])
  raw_eta <- X %*% t(B)
  # Matrix multiplication may legally simplify in downstream arithmetic when
  # a time slice has one active lineage. Restore explicit cell-by-lineage
  # dimensions before every vectorised reference operation; otherwise a
  # one-branch interval can be mistaken for a cell vector and silently break
  # the grid state at a speciation boundary.
  raw_eta <- matrix(raw_eta, nrow = nrow(X), ncol = nrow(B),
                    dimnames = list(as.character(carrier$track_id), lineages))
  ref <- reference_for(time_ma, lineages)
  reference_eta <- sweep(raw_eta, 2L, as.numeric(ref$eta_reference), "-")
  # The reconstructed response deliberately has no HMSC prevalence intercept.
  # P1 is therefore a transparent modern-reference contrast with a positive
  # modern response scale reconstructed along the dated tree. Unlike a
  # per-time-slice centering operation, this retains genuine palaeoclimate
  # shifts in absolute environmental support. It remains a relative support
  # covariate, never a historical occurrence probability on its own.
  response_scale <- as.numeric(ref$eta_scale)
  process_eta <- sweep(reference_eta, 2L, response_scale, "/")
  process_eta <- pmax(-8, pmin(8, process_eta))
  reference_eta <- matrix(reference_eta, nrow = nrow(raw_eta), ncol = ncol(raw_eta),
                          dimnames = dimnames(raw_eta))
  process_eta <- matrix(process_eta, nrow = nrow(raw_eta), ncol = ncol(raw_eta),
                        dimnames = dimnames(raw_eta))
  process_eta[!is.finite(process_eta)] <- NA_real_
  reference_eta[!is.finite(reference_eta)] <- NA_real_
  rownames(process_eta) <- rownames(reference_eta) <- carrier$track_id
  colnames(process_eta) <- colnames(reference_eta) <- lineages
  list(
    process_eta = process_eta,
    reference_eta = reference_eta,
    audit = data.frame(
      lineage = lineages,
      time_ma = time_ma,
      eta_reference = as.numeric(ref$eta_reference),
      eta_scale = response_scale,
      process_eta_available_q025 = vapply(seq_along(lineages), function(j) {
        z <- process_eta[carrier$active_land, j]
        if (any(is.finite(z))) as.numeric(stats::quantile(z, .025, na.rm = TRUE)) else NA_real_
      }, numeric(1)),
      process_eta_available_q50 = vapply(seq_along(lineages), function(j) {
        z <- process_eta[carrier$active_land, j]
        if (any(is.finite(z))) as.numeric(stats::quantile(z, .5, na.rm = TRUE)) else NA_real_
      }, numeric(1)),
      process_eta_available_q975 = vapply(seq_along(lineages), function(j) {
        z <- process_eta[carrier$active_land, j]
        if (any(is.finite(z))) as.numeric(stats::quantile(z, .975, na.rm = TRUE)) else NA_real_
      }, numeric(1)),
      stringsAsFactors = FALSE
    )
  )
}

apply_habitat <- function(q, slice) {
  active <- !is.na(slice$carrier$active_land) & slice$carrier$active_land
  if (!identical(rownames(q), as.character(slice$carrier$track_id))) {
    stop("Occupancy and carrier state have incompatible row order.", call. = FALSE)
  }
  q[!active, ] <- 0
  q
}

compact_root <- function(q, slice, eta) {
  root_eta <- rowMeans(eta, na.rm = TRUE)
  active <- which(!is.na(slice$carrier$active_land) &
                    slice$carrier$active_land & is.finite(root_eta))
  if (!length(active)) stop("No active root carriers have complete environmental support.")
  support <- stats::plogis(root_eta[active])
  center <- switch(
    root_prior_mode,
    max_environment = active[[which.max(support)]],
    uniform_land = sample(active, size = 1L),
    environment_weighted = {
      weight <- pmax(support, .Machine$double.eps)^root_support_power
      active[[sample.int(length(active), size = 1L, prob = weight)]]
    }
  )
  d <- HmscEcoEvo::hee_great_circle_distance_km(
    slice$carrier$paleo_lon[center], slice$carrier$paleo_lat[center],
    slice$carrier$paleo_lon[active], slice$carrier$paleo_lat[active]
  )
  chosen <- active[order(d, -support)][seq_len(min(root_max_carriers, length(active)))]
  # The root prior selects a compact environmentally supported patch. Its
  # occupancy mass is a declared root-range prior, not a second transformation
  # of the reference-calibrated P1 covariate.
  q[chosen, ] <- root_occupancy
  list(occupancy = q, center = center, chosen = chosen, support = support)
}

speciate_at <- function(q, time_ma) {
  hits <- which(abs(node_age - time_ma) <= 1e-5 * max(1, abs(time_ma)))
  hits <- hits[hits > n_tip & hits != root_node]
  if (!length(hits)) return(q)
  for (parent in hits) {
    parent_name <- lineage_name(parent)
    if (!parent_name %in% colnames(q)) next
    children <- tree$edge[tree$edge[, 1L] == parent, 2L]
    if (length(children) < 2L) next
    inherited <- q[, parent_name, drop = FALSE]
    q <- q[, setdiff(colnames(q), parent_name), drop = FALSE]
    for (child in children) {
      child_name <- lineage_name(child)
      q <- cbind(q, inherited)
      colnames(q)[ncol(q)] <- child_name
    }
  }
  q
}

write_metrics <- function(time_ma, slice, q, environment, process = NULL) {
  process_eta <- environment$process_eta
  reference_eta <- environment$reference_eta
  env_support <- stats::plogis(process_eta)
  reference_support <- stats::plogis(reference_eta)
  env_support[!is.finite(env_support)] <- NA_real_
  reference_support[!is.finite(reference_support)] <- NA_real_
  diversity <- HmscEcoEvo::hee_occupancy_weighted_diversity(
    q, zero_tolerance = diversity_min_expected_richness
  )
  arrival <- colonisation <- persistence <- rep(NA_real_, nrow(q))
  if (!is.null(process)) {
    arrival <- rowMeans(process$arrival_probability_reporting_interval, na.rm = TRUE)
    colonisation <- rowMeans(process$colonisation_probability_reporting_interval, na.rm = TRUE)
    persistence <- rowMeans(process$persistence_probability, na.rm = TRUE)
  }
  tab <- data.frame(
    track_id = rownames(q), time_ma = time_ma,
    paleo_lon = slice$carrier$paleo_lon, paleo_lat = slice$carrier$paleo_lat,
    active_land = slice$carrier$active_land,
    expected_sampled_surviving_lineage_richness = rowSums(q),
    # This is a threshold count on a relative P1 contrast, not expected
    # occupancy or a posterior species-richness estimate.
    reference_threshold_environmental_support_count = rowSums(
      is.finite(reference_eta) & reference_eta >= 0
    ),
    relative_environmental_support_mass = rowSums(env_support, na.rm = TRUE),
    mean_relative_environmental_support = rowMeans(env_support, na.rm = TRUE),
    mean_reference_calibrated_environmental_support = rowMeans(reference_support, na.rm = TRUE),
    mean_field_compositional_shannon_entropy = diversity$mean_field_compositional_shannon_entropy,
    mean_field_compositional_gini_simpson = diversity$mean_field_compositional_gini_simpson,
    mean_field_effective_lineages_shannon = diversity$mean_field_effective_lineages_shannon,
    mean_field_effective_lineages_simpson = diversity$mean_field_effective_lineages_simpson,
    # Backwards-compatible aliases only.  They are not posterior Shannon or
    # Simpson diversity and must not be placed in a result figure as such.
    occupancy_weighted_shannon = diversity$occupancy_weighted_shannon_entropy,
    occupancy_weighted_simpson = diversity$occupancy_weighted_gini_simpson,
    diversity_defined = diversity$diversity_defined,
    diversity_scientific_boundary = diversity$scientific_boundary,
    mean_arrival_probability_1myr_over_active_lineages = arrival,
    mean_colonisation_probability_1myr_over_active_lineages = colonisation,
    mean_conditional_persistence_probability_1myr_over_active_lineages = persistence,
    # Legacy aliases retained for old map renderers.
    mean_arrival_probability_1myr = arrival,
    mean_colonisation_probability_1myr = colonisation,
    mean_persistence_probability = persistence,
    n_active_lineages = ncol(q),
    stringsAsFactors = FALSE
  )
  saveRDS(tab, file.path(summary_dir, paste0("carrier_metrics_", time_slug(time_ma), "Ma.rds")), compress = "gzip")
  write_csv(environment$audit, file.path(
    gradient_dir, paste0("reference_contrast_", time_slug(time_ma), "Ma.csv")
  ))
  if (save_all_carrier_states) {
    saveRDS(list(occupancy = q, lineages = colnames(q), time_ma = time_ma,
                 state_space = "stable_plate_carrier"),
            file.path(state_dir, paste0("occupancy_", time_slug(time_ma), "Ma.rds")), compress = "gzip")
  }
  invisible(tab)
}

# A dated extant tree tells us that a branch eventually reaches its next
# observed node, but conditioning each numerical state on that event changes
# the forward transition kernel. Keep the state unconditioned here. A future
# particle/smoothing workflow scores full trajectories at the appropriate
# node or endpoint rather than reinflating their cell probabilities in place.
survival_audit <- function(q) {
  p_any <- 1 - exp(colSums(log1p(-q)))
  data.frame(
    lineage = colnames(q),
    expected_occupancy_mass = colSums(q),
    survival_probability_mean_field = p_any,
    zero_global_support = p_any <= 1e-300,
    stringsAsFactors = FALSE
  )
}

# Arrival, colonisation, and persistence are interval quantities.  They are
# written on the older slice where their movement graph and environmental
# state were evaluated, never pasted onto a younger grid with a different
# carrier set.
write_interval_process <- function(time_from_ma, time_to_ma, slice, report) {
  tab <- data.frame(
    track_id = as.character(slice$carrier$track_id),
    time_from_ma = time_from_ma, time_to_ma = time_to_ma,
    paleo_lon = slice$carrier$paleo_lon, paleo_lat = slice$carrier$paleo_lat,
    active_land = slice$carrier$active_land,
    mean_arrival_probability_1myr_over_active_lineages = report$arrival / report$total_duration_myr,
    mean_colonisation_probability_1myr_over_active_lineages = report$colonisation / report$total_duration_myr,
    mean_conditional_persistence_probability_1myr_over_active_lineages = report$persistence / report$total_duration_myr,
    # Legacy aliases retained for old map renderers.
    mean_arrival_probability_1myr = report$arrival / report$total_duration_myr,
    mean_colonisation_probability_1myr = report$colonisation / report$total_duration_myr,
    mean_persistence_probability = report$persistence / report$total_duration_myr,
    n_internal_process_steps = report$n_steps,
    stringsAsFactors = FALSE
  )
  saveRDS(tab, file.path(
    process_dir,
    paste0("carrier_interval_process_", time_slug(time_from_ma),
           "_to_", time_slug(time_to_ma), "Ma.rds")
  ), compress = "gzip")
  invisible(tab)
}

# Initial state at the oldest available environmental time. The dated-tree root
# is 0.05 Myr older; that short unrepresented interval is recorded rather than
# filled with a fabricated palaeoenvironment.
current_time <- env_times[[1L]]
current_slice <- first_slice
q_current <- q[carrier_ids, , drop = FALSE]
environment <- eta_for(current_slice, current_time, colnames(q_current))
root <- compact_root(q_current, current_slice, environment$process_eta)
q_current <- apply_habitat(root$occupancy, current_slice)
q[carrier_ids, colnames(q_current)] <- q_current
write_csv(data.frame(track_id = rep(carrier_ids[root$chosen], ncol(q)),
                     lineage = rep(colnames(q), each = length(root$chosen)),
                     time_ma = current_time,
                     root_center_track_id = carrier_ids[[root$center]],
                     root_prior_mode = root_prior_mode,
                     root_initial_occupancy = as.numeric(q_current[root$chosen, ]),
                     stringsAsFactors = FALSE),
          file.path(output, "08_root_initialisation", "root_compact_patch.csv"))
write_metrics(current_time, current_slice, q_current, environment)

interval_audit <- list(); audit_i <- 0L
for (ii in seq_len(length(env_times) - 1L)) {
  older <- env_times[[ii]]; younger <- env_times[[ii + 1L]]
  current_slice <- load_slice(older)
  current_ids <- as.character(current_slice$carrier$track_id)
  q[setdiff(rownames(q), current_ids), ] <- 0
  q_current <- apply_habitat(q[current_ids, , drop = FALSE], current_slice)
  breaks <- sort(unique(c(older,
                           node_times[node_times < older - 1e-6 & node_times >= younger - 1e-6],
                           younger)), decreasing = TRUE)
  process_report <- list(
    arrival = rep(0, nrow(q_current)),
    colonisation = rep(0, nrow(q_current)),
    persistence = rep(0, nrow(q_current)),
    n_steps = 0L,
    total_duration_myr = 0,
    survival_probability_mean_field = numeric()
  )
  for (bb in seq_len(length(breaks) - 1L)) {
    seg_start <- breaks[[bb]]; seg_end <- breaks[[bb + 1L]]; remaining <- seg_start - seg_end
    while (remaining > 1e-10) {
      dt <- min(internal_dt_myr, remaining)
      environment <- eta_for(current_slice, seg_start, colnames(q_current))
      eta <- environment$process_eta
      eta[!is.finite(eta)] <- -Inf
      process <- HmscEcoEvo::hee_projection_global_grid_step(
        occupancy = q_current, movement_matrix = current_slice$movement, eta = eta,
        delta_t = dt, reporting_interval_myr = 1,
        persistence_reference_interval = 1,
        establishment_intercept = establishment_intercept,
        establishment_slope = establishment_slope,
        persistence_intercept = persistence_intercept,
        persistence_slope = persistence_slope,
        movement_multiplier = movement_rate_multiplier,
        habitat_state = as.numeric(current_slice$carrier$active_land), stochastic = FALSE
      )
      process_report$arrival <- process_report$arrival + dt * rowMeans(
        process$arrival_probability_reporting_interval, na.rm = TRUE
      )
      process_report$colonisation <- process_report$colonisation + dt * rowMeans(
        process$colonisation_probability_reporting_interval, na.rm = TRUE
      )
      process_report$persistence <- process_report$persistence + dt * rowMeans(
        process$persistence_probability, na.rm = TRUE
      )
      process_report$n_steps <- process_report$n_steps + 1L
      process_report$total_duration_myr <- process_report$total_duration_myr + dt
      q_current <- apply_habitat(process$occupancy_next, current_slice)
      survival <- survival_audit(q_current)
      process_report$survival_probability_mean_field <- c(
        process_report$survival_probability_mean_field,
        survival$survival_probability_mean_field
      )
      remaining <- remaining - dt
    }
    q[current_ids, colnames(q_current)] <- q_current
    q <- speciate_at(q, seg_end)
    q_current <- q[current_ids, , drop = FALSE]
  }
  write_interval_process(older, younger, current_slice, process_report)
  next_slice <- load_slice(younger)
  next_ids <- as.character(next_slice$carrier$track_id)
  q[setdiff(rownames(q), next_ids), ] <- 0
  q_next <- apply_habitat(q[next_ids, , drop = FALSE], next_slice)
  transport_survival <- survival_audit(q_next)
  q[next_ids, colnames(q_next)] <- q_next
  environment_next <- eta_for(next_slice, younger, colnames(q_next))
  write_metrics(younger, next_slice, q_next, environment_next)
  audit_i <- audit_i + 1L
  interval_audit[[audit_i]] <- data.frame(
    time_from_ma = older, time_to_ma = younger,
    interval_myr = older - younger, n_active_lineages = ncol(q_next),
    mean_expected_lineage_richness = mean(rowSums(q_next)),
    n_active_carriers = sum(next_slice$carrier$active_land),
    output_coverage_fraction = next_slice$output_coverage,
    min_lineage_survival_probability_mean_field = if (length(process_report$survival_probability_mean_field)) min(process_report$survival_probability_mean_field) else NA_real_,
    median_lineage_survival_probability_mean_field = if (length(process_report$survival_probability_mean_field)) stats::median(process_report$survival_probability_mean_field) else NA_real_,
    min_lineage_survival_probability_mean_field_after_transport = if (nrow(transport_survival)) min(transport_survival$survival_probability_mean_field) else NA_real_,
    scientific_boundary = "unconditioned_forward_mean_field_survival_diagnostic_not_posterior_conditioning",
    stringsAsFactors = FALSE
  )
  message("Completed ", draw_id, ": ", older, " -> ", younger,
          " Ma; active lineages ", ncol(q), ".")
}
write_csv(do.call(rbind, interval_audit), file.path(output, "16_state_space_inference", "interval_audit.csv"))
write_csv(data.frame(
  response_draw = draw_id, status = "complete", n_environment_slices = length(env_times),
  n_final_tip_lineages = ncol(q), terminal_anchor_mode = "none",
  endpoint_role = "not_applied_no_independent_endpoint_data",
  extant_tree_survival_conditioning = "disabled_in_forward_scenario",
  stringsAsFactors = FALSE
), file.path(output, "16_state_space_inference", "draw_completion.csv"))
message("Case05 v3 draw complete: ", draw_id)
