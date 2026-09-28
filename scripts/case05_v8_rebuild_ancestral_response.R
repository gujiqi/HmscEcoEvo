#!/usr/bin/env Rscript

# Rebuild only the Case05 ancestral response. The geographic pruning model
# does not depend on beta, so its saved location distributions remain intact.

root <- "C:/Users/Google/Documents/HMSC-HIST"
base <- file.path(root, "outputs", "HmscEcoEvo")
run <- file.path(base, "case05_v8_full200_terrain_resistance_20260926")
comparison <- file.path(base, "case05_v8_four_scheme_comparison_20260926")
out <- file.path(comparison, "09_corrected_ancestral_response")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

input_tree <- file.path(root, "HmscEcoEvo_API8_20260818", "data_external",
  "plant_200_multifamily_extant_20260820", "prepared_inputs",
  "case03_data_raw", "Plant200_multifamily_extant_1deg_traits", "tree.tre")
beta_file <- file.path(base,
  "case05_plant200_4deg_effortaware_formal25_all66_20260923",
  "02_hmsc_4deg", "hmsc_environmental_beta_posterior_draws.rds")
selection <- utils::read.csv(file.path(run, "01_input_audit", "species_selection.csv"))
times <- utils::read.csv(file.path(run, "01_input_audit", "master_time_axis.csv"))$time_ma
draw_ids <- utils::read.csv(file.path(run, "04_phylogeographic_pruning",
  "pruning_log_likelihoods.csv"))$response_draw
tips <- selection$species[selection$selected]
tree <- ape::reorder.phylo(ape::keep.tip(ape::read.tree(input_tree), tips),
  "cladewise")

# Match the exact v8 dated-tree convention before changing any beta values.
depth <- ape::node.depth.edgelength(tree)
n_tip <- length(tree$tip.label)
root_age <- max(depth[seq_len(n_tip)])
tip_lag <- root_age - depth[seq_len(n_tip)]
if (max(abs(tip_lag)) > 1e-4) stop("Tree is not near-ultrametric.")
terminal_edge <- match(seq_len(n_tip), tree$edge[, 2L])
tree$edge.length[terminal_edge] <- tree$edge.length[terminal_edge] + tip_lag
root_age <- max(ape::node.depth.edgelength(tree)[seq_len(n_tip)])
if (root_age > max(times) + 1e-4) {
  tree$edge.length <- tree$edge.length * (max(times) / root_age)
}

beta <- readRDS(beta_file)
axes <- c("MAT_pohl_C", "MAP_pohl_mm_yr", "elevation_m",
  "T_seasonality_pohl_sd", "P_seasonality_pohl_sd",
  "moisture_availability_index_z", "wetland_potential_index_z",
  "relief_3x3_m")
if (!all(axes %in% names(beta))) stop("HMSC beta axes are incomplete.")
beta <- beta[beta$lineage %in% tree$tip.label &
  beta$response_draw %in% draw_ids,
  c("lineage", "response_draw", axes), drop = FALSE]
names(beta)[names(beta) == "lineage"] <- "species"
if (nrow(beta) != n_tip * length(draw_ids)) {
  stop("Tip beta posterior is incomplete for the stored v8 draws.")
}

response <- HmscEcoEvo::hee_evolution_ancestral_response_direct(
  tip_beta_draws = beta, tree = tree, times = times, basis_cols = axes,
  reconstruct_intercept = FALSE, branch_uncertainty = "brownian_bridge",
  seed = 20260923L)
old <- readRDS(file.path(run, "02_ancestral_response",
  "ancestral_beta_response_draws.rds"))
if (nrow(response) != nrow(old) || any(!is.finite(as.matrix(response[, axes])))) {
  stop("Corrected response rows or values do not match the v8 contract.")
}

root_new <- response[response$lineage_type == "root", , drop = FALSE]
root_old <- old[old$lineage_type == "root", , drop = FALSE]
root_comparison <- do.call(rbind, lapply(draw_ids, function(draw_id) {
  tip <- beta[beta$response_draw == draw_id, , drop = FALSE]
  new <- root_new[root_new$response_draw == draw_id, , drop = FALSE]
  previous <- root_old[root_old$response_draw == draw_id, , drop = FALSE]
  data.frame(response_draw = draw_id, axis = axes,
    old_root_beta = as.numeric(previous[1L, axes]),
    corrected_root_beta = as.numeric(new[1L, axes]),
    unweighted_tip_mean = colMeans(tip[, axes, drop = FALSE]),
    stringsAsFactors = FALSE)
}))
root_comparison$old_minus_tip_mean <- root_comparison$old_root_beta -
  root_comparison$unweighted_tip_mean
root_comparison$corrected_minus_tip_mean <-
  root_comparison$corrected_root_beta - root_comparison$unweighted_tip_mean

saveRDS(response, file.path(out, "ancestral_beta_response_draws.rds"),
  compress = "gzip")
utils::write.csv(root_comparison, file.path(out, "root_beta_correction_audit.csv"),
  row.names = FALSE)
utils::write.csv(data.frame(
  check = c("same_row_count_as_original", "all_corrected_beta_finite",
    "all_old_roots_equal_tip_mean", "corrected_roots_not_all_tip_mean",
    "five_draws", "all_tip_species"),
  pass = c(nrow(response) == nrow(old),
    all(is.finite(as.matrix(response[, axes]))),
    all(abs(root_comparison$old_minus_tip_mean) < 1e-8),
    any(abs(root_comparison$corrected_minus_tip_mean) > 1e-3),
    length(unique(response$response_draw)) == 5L,
    length(tree$tip.label) == 200L)
), file.path(out, "correction_quality_gates.csv"), row.names = FALSE)
message("Corrected ancestral response: ", nrow(response), " rows at ", out)
