test_that("Word report renderer writes a docx when rmarkdown and pandoc are available", {
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc is not available")
  rmd <- tempfile(fileext = ".Rmd")
  writeLines(c(
    "---",
    "title: 'Tiny HmscEcoEvo report'",
    "output: word_document",
    "params:",
    "  value: 1",
    "---",
    "",
    "# Table",
    "",
    "```{r}",
    "data.frame(metric = 'ok', value = params$value)",
    "```"
  ), rmd)
  out <- tempfile(fileext = ".docx")
  rendered <- hee_render_word_report(rmd, out, params = list(value = 2))
  expect_true(file.exists(rendered))
  expect_true(file.info(rendered)$size > 0)
})

test_that("complete workflow report indexes all figures and tables", {
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc is not available")
  root <- tempfile("hee-report-")
  dir.create(file.path(root, "figures"), recursive = TRUE)
  dir.create(file.path(root, "tables"), recursive = TRUE)
  dir.create(file.path(root, "reports"), recursive = TRUE)

  grDevices::png(file.path(root, "figures", "01_Y_heatmap.png"),
                 width = 320, height = 240)
  graphics::plot(1, 1, main = "Y")
  grDevices::dev.off()
  grDevices::png(file.path(root, "figures", "geoprocess_network_fragmentation_through_time.png"),
                 width = 320, height = 240)
  graphics::plot(1, 2, main = "Network")
  grDevices::dev.off()
  utils::write.csv(data.frame(a = 1), file.path(root, "tables", "Y_matrix.csv"),
                   row.names = FALSE)
  utils::write.csv(data.frame(a = 2), file.path(root, "tables", "geoprocess_network_summary.csv"),
                   row.names = FALSE)

  out <- hee_write_full_workflow_report(root)
  expect_true(file.exists(out))
  idx <- utils::read.csv(file.path(root, "tables",
                                   "full_workflow_report_completeness_index.csv"))
  expect_equal(sum(idx$type == "figure"), 2)
  expect_equal(sum(idx$type == "table"), 2)
  expect_true(all(idx$included_in_report[idx$type == "figure"]))
  expect_false(anyNA(idx$section))
})

test_that("complete workflow report states heuristic and non-causal boundaries explicitly", {
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc is not available")
  root <- tempfile("hee-report-boundaries-")
  dir.create(file.path(root, "figures"), recursive = TRUE)
  dir.create(file.path(root, "tables"), recursive = TRUE)
  dir.create(file.path(root, "reports"), recursive = TRUE)
  utils::write.csv(data.frame(file = "example.png"),
                   file.path(root, "tables", "output_file_inventory.csv"),
                   row.names = FALSE)

  hee_write_full_workflow_report(root)
  rmd <- readLines(file.path(root, "reports",
                            "hmscecoevo_full_540Ma_workflow.Rmd"),
                   warn = FALSE)
  txt <- paste(rmd, collapse = " ")
  expect_match(txt, "heuristic", ignore.case = TRUE)
  expect_match(txt, "non-causal|not causal|非因果", ignore.case = TRUE)
  expect_match(txt, "mock|MOCK_ONLY", ignore.case = TRUE)
  expect_match(txt, "proxy", ignore.case = TRUE)
})
