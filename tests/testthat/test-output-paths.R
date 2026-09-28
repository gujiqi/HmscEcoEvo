test_that("output paths are deterministic and file-safe", {
  p <- hee_output_path("probability", time_ma = 2.58, species_id = "sp one/a",
                       stat = "mean", root = "outputs")
  expect_match(p, "species=sp_one_a")
  expect_match(p, "time=002.6Ma")
  expect_match(p, "stat=mean.csv")
})

test_that("Case 2 template stops clearly when real data are missing", {
  script <- file.path("scripts", "case02_realdata_full_540Ma_template.R")
  skip_if_not(file.exists(script), "case02 script is not installed in this test context")
  out <- system2(file.path(R.home("bin"), "Rscript"),
                 c(script, "--data_root=does_not_exist_for_test"),
                 stdout = TRUE, stderr = TRUE)
  expect_true(any(grepl("Missing required file: data_raw/comm.csv", out, fixed = TRUE)))
})

test_that("Case 1 projection task table records completed output when available", {
  task_file <- file.path("outputs", "case01_simulated_full_540Ma",
                         "tables", "projection_task_table.csv")
  skip_if_not(file.exists(task_file), "formal Case 1 output is not installed in this test context")
  tasks <- read.csv(task_file, stringsAsFactors = FALSE)
  expect_true(all(c("task_id", "tile_id", "row_range", "status",
                    "output_file", "input_hash") %in% names(tasks)))
  expect_false(any(tasks$status == "pending"))
  expect_true(all(file.exists(unique(tasks$output_file))))
})
