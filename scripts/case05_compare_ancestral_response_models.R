#!/usr/bin/env Rscript

# Compare BM, OU and EB on the saved Plant200 HMSC Beta draws and dated tree.
# This reconstructs response coefficients only; no geographic occupancy is fit.

root <- "C:/Users/Google/Documents/HMSC-HIST"
base <- file.path(root, "outputs", "HmscEcoEvo")
run <- file.path(base, "case05_v8_full200_terrain_resistance_20260926")
comparison <- file.path(base, "case05_v8_four_scheme_comparison_20260926")
quick <- identical(Sys.getenv("HMSCEE_RESPONSE_QUICK"), "1")
out <- file.path(comparison,
  if (quick) "11_ancestral_response_model_quick_check" else
    "11_ancestral_response_BM_OU_EB")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

input_tree <- file.path(root, "HmscEcoEvo_API8_20260818", "data_external",
  "plant_200_multifamily_extant_20260820", "prepared_inputs",
  "case03_data_raw", "Plant200_multifamily_extant_1deg_traits", "tree.tre")
beta_file <- file.path(base,
  "case05_plant200_4deg_effortaware_formal25_all66_20260923",
  "02_hmsc_4deg", "hmsc_environmental_beta_posterior_draws.rds")
selection <- utils::read.csv(file.path(run, "01_input_audit",
  "species_selection.csv"))
times <- utils::read.csv(file.path(run, "01_input_audit",
  "master_time_axis.csv"))$time_ma
draw_ids <- utils::read.csv(file.path(run, "04_phylogeographic_pruning",
  "pruning_log_likelihoods.csv"))$response_draw
tips <- selection$species[selection$selected]
tree <- ape::reorder.phylo(ape::keep.tip(ape::read.tree(input_tree), tips),
  "cladewise")
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
  stop("Tip beta posterior is incomplete for the saved Case05 draws.")
}
if (quick) {
  beta <- beta[beta$response_draw == draw_ids[1L],
               c("species", "response_draw", axes[1L]), drop = FALSE]
  axes <- axes[1L]
  times <- sort(unique(c(max(times), 200, 100, 0)), decreasing = TRUE)
}

models <- if (quick) c("BM", "OU", "EB", "AICc") else
  c("BM", "OU", "EB", "AICc")
quality <- list()
comparison_rows <- list()
for (model in models) {
  message("Reconstructing ", model, ": ", n_tip, " tips, ",
          length(unique(beta$response_draw)), " posterior draw(s), ",
          length(axes), " axis/axes, ", length(times), " times")
  elapsed <- system.time({
    response <- HmscEcoEvo::hee_evolution_ancestral_response_direct(
      tip_beta_draws = beta, tree = tree, times = times, basis_cols = axes,
      evolution_model = model, reconstruct_intercept = FALSE,
      branch_uncertainty = "process_bridge", seed = 20260926L)
  })[["elapsed"]]
  fit <- attr(response, "model_comparison")
  fit$requested_model <- model
  comparison_rows[[model]] <- fit
  saveRDS(response, file.path(out, paste0("ancestral_beta_", model, ".rds")),
          compress = "gzip")
  utils::write.csv(fit, file.path(out, paste0("model_fit_", model, ".csv")),
                   row.names = FALSE)
  tip_rows <- response[response$time_ma == 0 &
    response$lineage %in% tree$tip.label, , drop = FALSE]
  original <- beta[match(paste(tip_rows$response_draw, tip_rows$lineage),
                         paste(beta$response_draw, beta$species)),
                   axes, drop = FALSE]
  tip_match <- isTRUE(all.equal(as.matrix(tip_rows[, axes, drop = FALSE]),
                                as.matrix(original), tolerance = 1e-8,
                                check.attributes = FALSE))
  quality[[model]] <- data.frame(
    model = model, n_tips = n_tip, n_draws = length(unique(response$response_draw)),
    n_times = length(unique(response$time_ma)), n_rows = nrow(response),
    n_axes = length(axes), all_beta_finite = all(is.finite(as.matrix(
      response[, axes, drop = FALSE]))),
    modern_tip_beta_recovered = tip_match,
    model_fit_rows = nrow(fit), elapsed_seconds = elapsed,
    stringsAsFactors = FALSE)
  message(model, " done in ", round(elapsed, 1), " s; rows=", nrow(response))
}
utils::write.csv(do.call(rbind, quality), file.path(out, "model_run_quality.csv"),
                 row.names = FALSE)
utils::write.csv(do.call(rbind, comparison_rows),
                 file.path(out, "all_model_fits.csv"), row.names = FALSE)
if ("AICc" %in% names(comparison_rows)) {
  selected <- comparison_rows$AICc[comparison_rows$AICc$selected, , drop = FALSE]
  counts <- table(selected$model, selected$boundary_status)
  lines <- c(
    "# Case05 ancestral environmental-response model comparison",
    "",
    "These outputs reconstruct HMSC response coefficients on the dated Plant200 tree.",
    "They are not palaeogeographic occupancy or historical dispersal maps.",
    "",
    paste0("Tips: ", n_tip, "; HMSC draws: ", length(unique(beta$response_draw)),
           "; axes: ", length(axes), "; requested times: ", length(times), "."),
    "",
    "BM is the backward-compatible baseline. OU is a single-optimum-like",
    "Gaussian process; EB allows the evolutionary rate to decline from the root.",
    "AICc is computed per environmental axis and HMSC posterior draw.",
    "It treats each tip Beta draw as exact and ignores tree uncertainty; do not",
    "interpret these scores as joint posterior model probabilities.",
    "OU support is not evidence of stabilizing selection.",
    "Fitted OU/EB parameters on a bound signal weak identifiability.",
    "",
    "## AICc selections by boundary status",
    "",
    paste(capture.output(print(counts)), collapse = "\n"),
    "",
    "See model_fit_AICc.csv for each axis/draw/parameter/AICc/boundary flag.",
    "The RDS files contain full lineage-time response trajectories."
  )
  writeLines(lines, file.path(out, "MODEL_INTERPRETATION.md"), useBytes = TRUE)
}
if (!all(vapply(quality, function(x) x$all_beta_finite &&
                x$modern_tip_beta_recovered, logical(1)))) {
  stop("At least one ancestral-response model failed quality gates.")
}
message("Model comparison saved to ", out)
