# Run the real Akodon full HMSC results workflow ------------------------------
# From the package root:
#   source("inst/examples/run_real_akodon_full_results_zh.R")
#
# From an installed package:
#   source(system.file("examples/run_real_akodon_full_results_zh.R", package = "HmscEcoEvo"))

this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
if (!is.na(this_file) && file.exists(this_file)) {
  example_dir <- dirname(this_file)
} else if (dir.exists(file.path("inst", "examples"))) {
  example_dir <- file.path("inst", "examples")
} else {
  example_dir <- system.file("examples", package = "HmscEcoEvo")
}

if (!nzchar(example_dir) || !dir.exists(example_dir)) {
  stop("Cannot find hmscHist examples directory.", call. = FALSE)
}

source(file.path(example_dir, "real_akodon_full_results_workflow", "00_run_all_real_akodon_results_zh.R"))
