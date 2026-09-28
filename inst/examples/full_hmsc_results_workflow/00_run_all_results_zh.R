# Toy full HMSC results workflow ---------------------------------------------
# This workflow uses the bundled toy data. For the real Akodon workflow, run:
#   inst/examples/run_real_akodon_full_results_zh.R

localDir <- file.path(getwd(), "hmscHist_full_results_example")
if (!dir.exists(localDir)) dir.create(localDir, recursive = TRUE)

scripts <- c(
  "S1_define_models_hmscHist.R",
  "S2_fit_models_short.R",
  "S3_evaluate_convergence.R",
  "S4_compute_model_fit.R",
  "S5_show_model_fit_and_variance.R",
  "S6_show_parameter_estimates.R",
  "S7_make_predictions.R"
)

local_script_dir <- file.path("inst", "examples", "full_hmsc_results_workflow")
if (dir.exists(local_script_dir)) {
  script_dir <- local_script_dir
} else {
  script_dir <- system.file("examples", "full_hmsc_results_workflow", package = "HmscEcoEvo")
}
if (!nzchar(script_dir) || !dir.exists(script_dir)) {
  stop("Cannot find full_hmsc_results_workflow scripts.", call. = FALSE)
}

for (s in scripts) {
  cat("\n============================================================\n")
  cat("Running", s, "\n")
  cat("============================================================\n")
  source(file.path(script_dir, s), local = TRUE)
}

cat("\nToy full workflow finished. Results directory:\n")
cat(localDir, "\n")
cat("Start with results/ALL_RESULTS_INDEX.csv.\n")
