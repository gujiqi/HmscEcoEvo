test_that("Case 03 nested demo quick workflow runs end to end", {
  script_candidates <- c(
    system.file("scripts", "case03_geoprocess_past_to_present_simulated.R",
                package = "HmscEcoEvo", mustWork = FALSE),
    file.path("scripts", "case03_geoprocess_past_to_present_simulated.R"),
    file.path("..", "..", "scripts", "case03_geoprocess_past_to_present_simulated.R")
  )
  script <- script_candidates[file.exists(script_candidates)][1]
  skip_if(is.na(script), "Case 03 compatibility wrapper is unavailable")

  out_dir <- file.path(tempdir(), paste0("case03_nested_wrapper_", Sys.getpid()))
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(normalizePath(script, winslash = "/", mustWork = TRUE),
      "--quick=true",
      paste0("--output=", normalizePath(out_dir, winslash = "/",
                                        mustWork = FALSE))),
    stdout = TRUE,
    stderr = TRUE
  )
  expect_null(attr(status, "status"), info = paste(status, collapse = "\n"))
  expect_true(file.exists(file.path(out_dir, "metadata", "analysis_mode.txt")))
  expect_equal(readLines(file.path(out_dir, "metadata", "analysis_mode.txt")),
               "demo")
  expect_true(file.exists(file.path(out_dir, "tables",
                                    "00.1_hmscee_nested_formula_catalog.csv")))
  expect_true(file.exists(file.path(out_dir, "tables",
                                    "01.1_paleo_earth_X_H_state.csv")))
  expect_true(file.exists(file.path(out_dir, "tables",
                                    "02.1_bsm_region_history_R.csv")))
  expect_true(file.exists(file.path(out_dir, "tables",
                                    "03.1_within_region_movement_kernel_K.csv")))
  expect_true(file.exists(file.path(out_dir, "tables",
                                    "04.2_cell_probability_summary_p.csv")))
  expect_true(file.exists(file.path(out_dir, "rds",
                                    "case03_hmscee_nested_demo_fit.rds")))

  p <- utils::read.csv(file.path(out_dir, "tables",
                                 "04.2_cell_probability_summary_p.csv"))
  expect_true(all(is.finite(p$probability_mean)))
  expect_true(all(p$probability_mean >= 0 & p$probability_mean <= 1))
  expect_true(all(p$probability_mean[p$H_state <= 0] == 0))

  R <- utils::read.csv(file.path(out_dir, "tables",
                                 "02.1_bsm_region_history_R.csv"))
  expect_true(all(R$R_region %in% c(0, 1)))
  K <- utils::read.csv(file.path(out_dir, "tables",
                                 "03.1_within_region_movement_kernel_K.csv"))
  expect_true(all(K$K_movement >= 0 & K$K_movement <= 1))
  expect_true(all(K$movement_scale == "within_BioGEOBEARS_region_only") ||
                all(K$movement_scale == "within_BioGeoBEARS_region_only"))

  catalog <- utils::read.csv(file.path(out_dir, "tables",
                                       "00.1_hmscee_nested_formula_catalog.csv"),
                             stringsAsFactors = FALSE)
  expect_true(all(c("regional_history", "movement_kernel",
                    "final_probability") %in% catalog$formula_id))
  expect_false(any(grepl("S_HMSC.*A_hist", catalog$formula)))
})
