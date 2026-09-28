#!/usr/bin/env Rscript

# Case05 trait-informed process-history constructor
#
# This script creates an explicit lineage x time x HMSC-draw process scenario
# for Case04/Case05.  It deliberately does not estimate historical demographic
# rates from the modern HMSC.  Instead, it propagates observed/imputed tip
# traits through the dated tree and applies transparent, user-set coefficients
# as a prior/sensitivity scenario.  This keeps the five/eight-process roles
# separate: beta is P1 environmental filtering; trait-informed movement,
# establishment, and persistence are P2/Persistence inputs.

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(x, y) if (is.null(x)) y else x
parse_args <- function(x) {
  out <- list()
  for (arg in x) {
    if (!startsWith(arg, "--")) next
    key_value <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    key <- key_value[[1L]]
    value <- if (length(key_value) > 1L) paste(key_value[-1L], collapse = "=") else "TRUE"
    out[[key]] <- value
  }
  out
}
cfg <- parse_args(args)
need <- c("template", "tree", "traits", "output")
missing_args <- need[!vapply(need, function(x) nzchar(cfg[[x]] %||% ""), logical(1))]
if (length(missing_args)) {
  stop("Missing required command-line options: ",
       paste(paste0("--", missing_args, "="), collapse = ", "), call. = FALSE)
}

if (!requireNamespace("HmscEcoEvo", quietly = TRUE)) {
  stop("The current HmscEcoEvo package is not installed. Run devtools::install_local() first.", call. = FALSE)
}

arg_num <- function(name, default) {
  value <- suppressWarnings(as.numeric(cfg[[name]] %||% default))
  if (length(value) != 1L || !is.finite(value)) {
    stop("--", name, " must be one finite numeric value.", call. = FALSE)
  }
  value
}

template_path <- normalizePath(cfg$template, mustWork = TRUE)
tree_path <- normalizePath(cfg$tree, mustWork = TRUE)
traits_path <- normalizePath(cfg$traits, mustWork = TRUE)
output_path <- normalizePath(cfg$output, mustWork = FALSE)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

# Coefficients are intentionally explicit sensitivity assumptions.  They act
# on z-standardised reconstructed traits and never turn traits into beta.
base_movement <- arg_num("base_movement", 1)
movement_seed <- arg_num("movement_seed", -0.20)
movement_height <- arg_num("movement_height", 0.15)
movement_min <- arg_num("movement_min", 0.35)
movement_max <- arg_num("movement_max", 2.50)
base_establishment <- arg_num("base_establishment", -1.50)
establishment_seed <- arg_num("establishment_seed", 0.10)
establishment_ldmc <- arg_num("establishment_ldmc", 0.08)
establishment_slope <- arg_num("establishment_slope", 0.85)
base_persistence <- arg_num("base_persistence", 0)
persistence_wood <- arg_num("persistence_wood", 0.16)
persistence_lifespan <- arg_num("persistence_lifespan", 0.16)
persistence_ldmc <- arg_num("persistence_ldmc", 0.08)
persistence_slope <- arg_num("persistence_slope", 0.75)

if (movement_min <= 0 || movement_max < movement_min) {
  stop("movement_min must be > 0 and movement_max must be >= movement_min.", call. = FALSE)
}

template <- utils::read.csv(template_path, stringsAsFactors = FALSE, check.names = FALSE)
required_template <- c("lineage", "time_ma", "response_draw")
if (!all(required_template %in% names(template))) {
  stop("template must contain: ", paste(required_template, collapse = ", "), call. = FALSE)
}
template$lineage <- as.character(template$lineage)
template$response_draw <- as.character(template$response_draw)
template$time_ma <- suppressWarnings(as.numeric(template$time_ma))
if (any(!nzchar(template$lineage)) || any(!nzchar(template$response_draw)) ||
    any(!is.finite(template$time_ma)) || anyDuplicated(template[required_template])) {
  stop("template has empty, non-finite, or duplicated lineage/time/draw keys.", call. = FALSE)
}

traits <- utils::read.csv(traits_path, stringsAsFactors = FALSE, check.names = FALSE)
trait_axes <- c(
  "seed_mass", "maximum_whole_plant_height", "stem_wood_density",
  "leaf_life_span", "leaf_dry_mass_per_leaf_fresh_mass"
)
if (!all(c("species", trait_axes) %in% names(traits))) {
  stop("traits must contain species plus: ", paste(trait_axes, collapse = ", "), call. = FALSE)
}
traits$species <- as.character(traits$species)
if (any(!nzchar(traits$species)) || anyDuplicated(traits$species)) {
  stop("traits$species must be unique and non-empty.", call. = FALSE)
}
for (axis in trait_axes) {
  traits[[axis]] <- suppressWarnings(as.numeric(traits[[axis]]))
  if (any(!is.finite(traits[[axis]]))) {
    stop("traits contains non-finite values for ", axis,
         ". Supply a multiple-imputation draw with finite values.", call. = FALSE)
  }
  sd_axis <- stats::sd(traits[[axis]])
  if (!is.finite(sd_axis) || sd_axis <= 0) {
    stop("trait axis has zero/non-finite variance: ", axis, call. = FALSE)
  }
  traits[[axis]] <- as.numeric(scale(traits[[axis]]))
}

tree <- ape::read.tree(tree_path)
if (is.null(tree$edge.length) || any(!is.finite(tree$edge.length)) || any(tree$edge.length < 0)) {
  stop("tree must be dated and have non-negative finite branch lengths.", call. = FALSE)
}

# The dynamic engine first drops all tree tips not represented by the selected
# HMSC posterior.  Reproduce that exact pruning here: internal node numbers
# and their labels otherwise differ even when the source Newick is identical.
engine_tips <- intersect(unique(template$lineage), tree$tip.label)
if (length(engine_tips) < 2L) {
  stop("Could not recover at least two HMSC tip names from the engine template.",
       " Use the template emitted by the same Case04/Case05 run.", call. = FALSE)
}
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, engine_tips))
missing_trait_tips <- setdiff(tree$tip.label, traits$species)
if (length(missing_trait_tips)) {
  stop("traits is missing selected HMSC tree tips: ",
       paste(missing_trait_tips, collapse = ", "), call. = FALSE)
}
traits <- traits[match(tree$tip.label, traits$species), , drop = FALSE]

# The same dated tree and every requested template time are used, so the
# result has precisely the lineage identifiers that the dynamic engine needs.
tip_traits <- traits[, c("species", trait_axes), drop = FALSE]
tip_traits$response_draw <- "trait_imputation_001"
tip_traits <- tip_traits[, c("species", "response_draw", trait_axes), drop = FALSE]
trait_history <- HmscEcoEvo::hee_evolution_ancestral_response_direct(
  tip_beta_draws = tip_traits,
  tree = tree,
  times = sort(unique(template$time_ma), decreasing = TRUE),
  basis_cols = trait_axes,
  species_col = "species",
  response_draw_col = "response_draw",
  reconstruct_intercept = FALSE,
  branch_uncertainty = "none"
)
trait_history <- trait_history[, c("lineage", "time_ma", trait_axes), drop = FALSE]
if (anyDuplicated(trait_history[c("lineage", "time_ma")])) {
  stop("Trait reconstruction unexpectedly created duplicate lineage/time keys.", call. = FALSE)
}

key_template <- paste(template$lineage, format(template$time_ma, scientific = FALSE, trim = TRUE), sep = "|")
key_history <- paste(trait_history$lineage, format(trait_history$time_ma, scientific = FALSE, trim = TRUE), sep = "|")
idx <- match(key_template, key_history)
if (anyNA(idx)) {
  missing <- unique(template[is.na(idx), c("lineage", "time_ma"), drop = FALSE])
  missing_path <- file.path(dirname(output_path), "missing_trait_history_keys.csv")
  utils::write.csv(missing, missing_path, row.names = FALSE)
  stop("Trait reconstruction does not cover every engine lineage/time key. See: ",
       missing_path, call. = FALSE)
}
history <- trait_history[idx, trait_axes, drop = FALSE]
rownames(history) <- NULL

# Re-centre reconstructed traits within each active lineage-time slice.  The
# scalar baselines are calibrated separately to the 0 Ma endpoint; traits must
# express relative lineage differences rather than silently change the entire
# global demographic rate because ancestors have a different mean trait value
# from modern tips.  This is the appropriate first-order trait scenario until
# independent historical demographic observations identify absolute effects.
centered_history <- history
time_groups <- split(seq_len(nrow(template)), template$time_ma)
for (rows in time_groups) {
  for (axis in trait_axes) {
    centered_history[[axis]][rows] <- centered_history[[axis]][rows] -
      mean(centered_history[[axis]][rows])
  }
}

movement_deviation <- movement_seed * centered_history$seed_mass +
  movement_height * centered_history$maximum_whole_plant_height
movement_multiplier <- numeric(nrow(template))
for (rows in time_groups) {
  # The geometric tendency is centred to an arithmetic mean of base_movement,
  # avoiding an unintended log-normal inflation of total regional propagule
  # pressure while retaining trait-based relative movement differences.
  raw_multiplier <- exp(movement_deviation[rows])
  movement_multiplier[rows] <- base_movement * raw_multiplier / mean(raw_multiplier)
}
movement_multiplier <- pmin(pmax(movement_multiplier, movement_min), movement_max)
establishment_intercept <- base_establishment +
  establishment_seed * centered_history$seed_mass +
  establishment_ldmc * centered_history$leaf_dry_mass_per_leaf_fresh_mass
persistence_intercept <- base_persistence +
  persistence_wood * centered_history$stem_wood_density +
  persistence_lifespan * centered_history$leaf_life_span +
  persistence_ldmc * centered_history$leaf_dry_mass_per_leaf_fresh_mass

out <- template[, required_template, drop = FALSE]
out$movement_multiplier <- movement_multiplier
out$establishment_intercept <- establishment_intercept
out$establishment_slope <- establishment_slope
out$persistence_intercept <- persistence_intercept
out$persistence_slope <- persistence_slope
if (any(!is.finite(as.matrix(out[, -(1:3), drop = FALSE]))) || any(out$movement_multiplier <= 0)) {
  stop("Generated process parameters are invalid.", call. = FALSE)
}
utils::write.csv(out, output_path, row.names = FALSE)

trait_output <- cbind(template[, required_template, drop = FALSE], history)
names(centered_history) <- paste0("centred_", names(centered_history))
trait_output <- cbind(trait_output, centered_history)
utils::write.csv(trait_output,
                 file.path(dirname(output_path), "trait_history_used_for_process_scenario.csv"),
                 row.names = FALSE)
summary_rows <- data.frame(
  parameter = c("movement_multiplier", "establishment_intercept", "persistence_intercept"),
  minimum = c(min(out$movement_multiplier), min(out$establishment_intercept), min(out$persistence_intercept)),
  median = c(stats::median(out$movement_multiplier), stats::median(out$establishment_intercept), stats::median(out$persistence_intercept)),
  maximum = c(max(out$movement_multiplier), max(out$establishment_intercept), max(out$persistence_intercept)),
  stringsAsFactors = FALSE
)
utils::write.csv(summary_rows,
                 file.path(dirname(output_path), "trait_process_parameter_summary.csv"),
                 row.names = FALSE)
meta <- data.frame(
  field = c("scientific_status", "trait_input", "tree", "template", "trait_imputation", "movement_model", "establishment_model", "persistence_model"),
  value = c(
    "trait_informed_prior_scenario_not_empirically_estimated",
    traits_path, tree_path, template_path,
    "single imputed trait table; trait-imputation uncertainty is not yet propagated",
    "within-time centred exp(b_seed*z(seed_mass)+b_height*z(height)); rescaled to the supplied base_movement",
    "base_establishment+b_seed*centred_z(seed_mass)+b_ldmc*centred_z(LDMC)",
    "base_persistence+b_wood*centred_z(wood_density)+b_lifespan*centred_z(leaf_life_span)+b_ldmc*centred_z(LDMC)"
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(meta,
                 file.path(dirname(output_path), "trait_process_parameter_metadata.csv"),
                 row.names = FALSE)
cat("Wrote trait-informed lineage-time process scenario:\n", output_path, "\n", sep = "")
