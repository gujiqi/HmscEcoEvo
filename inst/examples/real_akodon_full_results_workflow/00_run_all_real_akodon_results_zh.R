# Real Akodon full HMSC results workflow -------------------------------------
# This runs the S1-S7 result workflow with the real Akodon data bundled in
# inst/extdata/real_akodon.
#
# Output directory:
#   hmscHist_real_akodon_full_results/

localDir <- file.path(getwd(), "hmscHist_real_akodon_full_results")
if (!dir.exists(localDir)) dir.create(localDir, recursive = TRUE)

# Default MCMC settings for this real-data example. Increase these for formal
# analysis. You can override them before sourcing this file with:
# options(hmscHist.mcmc_settings = list(samples = 1000, transient = 1000, ...))
default_real_akodon_mcmc <- list(
  thin = 1,
  samples = 300,
  transient = 300,
  nChains = 2,
  nParallel = 1,
  verbose = 100
)
real_akodon_mcmc <- utils::modifyList(
  default_real_akodon_mcmc,
  getOption("hmscHist.mcmc_settings", list())
)
options(hmscHist.mcmc_settings = real_akodon_mcmc)

local_real_dir <- file.path("inst", "examples", "real_akodon_full_results_workflow")
if (dir.exists(local_real_dir)) {
  real_dir <- local_real_dir
} else {
  real_dir <- system.file("examples", "real_akodon_full_results_workflow", package = "HmscEcoEvo")
}
if (!nzchar(real_dir) || !dir.exists(real_dir)) {
  stop("Cannot find real_akodon_full_results_workflow scripts.", call. = FALSE)
}

local_common_dir <- file.path("inst", "examples", "full_results_workflow")
if (dir.exists(local_common_dir)) {
  common_dir <- local_common_dir
} else {
  common_dir <- system.file("examples", "full_results_workflow", package = "HmscEcoEvo")
}
if (!nzchar(common_dir) || !dir.exists(common_dir)) {
  stop("Cannot find full_results_workflow scripts.", call. = FALSE)
}

steps <- list(
  file.path(real_dir, "S1_define_real_akodon_model.R"),
  file.path(common_dir, "S2_fit_models_short.R"),
  file.path(common_dir, "S3_evaluate_convergence.R"),
  file.path(common_dir, "S4_compute_model_fit.R"),
  file.path(common_dir, "S5_show_model_fit_and_variance.R"),
  file.path(common_dir, "S6_show_parameter_estimates.R"),
  file.path(common_dir, "S7_make_predictions.R")
)

for (s in steps) {
  cat("\n============================================================\n")
  cat("Running", basename(s), "\n")
  cat("============================================================\n")
  source(s, local = TRUE)
}

cat("\nReal Akodon full workflow finished. Results directory:\n")
cat(localDir, "\n")
cat("Start with results/ALL_RESULTS_INDEX.csv.\n")
