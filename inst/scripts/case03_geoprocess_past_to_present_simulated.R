#!/usr/bin/env Rscript

# Compatibility wrapper. The old Case03 product-model script was removed.
# Use the rebuilt nested HmscEE demo instead.

args_full <- commandArgs(FALSE)
cmd_file <- args_full[grep("^--file=", args_full)]
this_file <- if (length(cmd_file)) {
  normalizePath(sub("^--file=", "", cmd_file[[1]]), winslash = "/",
                mustWork = FALSE)
} else {
  normalizePath("scripts/case03_geoprocess_past_to_present_simulated.R",
                winslash = "/", mustWork = FALSE)
}
script <- file.path(dirname(this_file), "case03_geoprocess_simulated_demo.R")
if (!file.exists(script)) script <- "scripts/case03_geoprocess_simulated_demo.R"
status <- system2(file.path(R.home("bin"), "Rscript"),
                  c(normalizePath(script, winslash = "/", mustWork = TRUE),
                    "--analysis_mode=demo", commandArgs(trailingOnly = TRUE)))
quit(save = "no", status = status)
