case03_script_path <- function(name) {
  candidates <- c(
    system.file("scripts", name, package = "HmscEcoEvo", mustWork = FALSE),
    file.path("scripts", name),
    file.path("..", "..", "scripts", name)
  )
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) return(NA_character_)
  normalizePath(hit[[1]], winslash = "/", mustWork = TRUE)
}

test_that("Case 03 scripts use the nested framework and not the obsolete product core", {
  real_script <- case03_script_path("case03_geoprocess_realdata_past_to_present.R")
  demo_script <- case03_script_path("case03_geoprocess_simulated_demo.R")
  wrapper_script <- case03_script_path("case03_geoprocess_past_to_present_simulated.R")
  expect_false(is.na(real_script))
  expect_false(is.na(demo_script))
  expect_false(is.na(wrapper_script))

  real_code <- paste(readLines(real_script, warn = FALSE), collapse = "\n")
  demo_code <- paste(readLines(demo_script, warn = FALSE), collapse = "\n")
  expect_match(real_code, "hee_bsm_region_history", fixed = TRUE)
  expect_match(demo_code, "hee_nested_region_cell_occupancy", fixed = TRUE)
  expect_match(demo_code, "hee_within_region_movement_kernel", fixed = TRUE)
  expect_false(grepl("hee_dynamic_earth_biota_probability\\s*\\(",
                     real_code, perl = TRUE))
  expect_false(grepl("hee_dynamic_earth_biota_probability\\s*\\(",
                     demo_code, perl = TRUE))
  expect_false(grepl("S_HMSC\\s*\\*\\s*A_hist", real_code))
  expect_false(grepl("S_HMSC\\s*\\*\\s*A_hist", demo_code))
})

test_that("real Case 03 writes a missing-input manifest and stops", {
  real_script <- case03_script_path("case03_geoprocess_realdata_past_to_present.R")
  skip_if(is.na(real_script), "real Case 03 script is not available")
  out_dir <- file.path(tempdir(), paste0("case03_nested_real_missing_", Sys.getpid()))
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  result <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(real_script, "--analysis_mode=real",
      paste0("--output=", normalizePath(out_dir, winslash = "/",
                                        mustWork = FALSE))),
    stdout = TRUE, stderr = TRUE
  ))
  expect_false(is.null(attr(result, "status")))
  manifest <- file.path(out_dir, "missing_required_inputs.csv")
  expect_true(file.exists(manifest), info = paste(result, collapse = "\n"))
  missing <- utils::read.csv(manifest, stringsAsFactors = FALSE)
  expect_true(all(c("comm", "sites", "env_rds", "region_cube", "tree",
                    "bgb_or_bsm_region_history") %in% missing$input))
  expect_true(all(missing$status == "MISSING"))
})

test_that("simulated Case 03 demo runs the nested HmscEE functions", {
  demo_script <- case03_script_path("case03_geoprocess_simulated_demo.R")
  skip_if(is.na(demo_script), "simulated Case 03 script is unavailable")
  out_dir <- file.path(tempdir(), paste0("case03_nested_demo_", Sys.getpid()))
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  result <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(demo_script, "--analysis_mode=demo", "--quick=true",
      paste0("--output=", normalizePath(out_dir, winslash = "/",
                                        mustWork = FALSE))),
    stdout = TRUE, stderr = TRUE
  ))
  expect_null(attr(result, "status"), info = paste(result, collapse = "\n"))
  expect_true(file.exists(file.path(out_dir, "tables",
                                    "00.1_hmscee_nested_formula_catalog.csv")))
  expect_true(file.exists(file.path(out_dir, "tables",
                                    "04.2_cell_probability_summary_p.csv")))
  p <- utils::read.csv(file.path(out_dir, "tables",
                                 "04.2_cell_probability_summary_p.csv"))
  expect_true(all(p$probability_mean >= 0 & p$probability_mean <= 1))
  expect_true(all(p$H_state[p$probability_mean > 0] > 0))
  f <- utils::read.csv(file.path(out_dir, "tables",
                                 "00.1_hmscee_nested_formula_catalog.csv"),
                       stringsAsFactors = FALSE)
  expect_true("final_probability" %in% f$formula_id)
  expect_false(any(grepl("A_hist.*D_dynamic.*Q_extrap", f$formula)))
})
