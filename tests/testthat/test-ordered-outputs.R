test_that("ordered output catalog and copies are created", {
  root <- tempfile("hee-ordered-")
  dir.create(file.path(root, "figures"), recursive = TRUE)
  dir.create(file.path(root, "tables"), recursive = TRUE)
  dir.create(file.path(root, "reports"), recursive = TRUE)
  writeBin(charToRaw("png"), file.path(root, "figures", "01_Y_heatmap.png"))
  writeBin(charToRaw("pdf"), file.path(root, "figures", "01_Y_heatmap.pdf"))
  write.csv(data.frame(a = 1), file.path(root, "tables", "Y_matrix.csv"), row.names = FALSE)
  writeLines("report", file.path(root, "reports", "hmscecoevo_full_540Ma_workflow.Rmd"))

  catalog <- hee_write_ordered_outputs(root, include_rds = TRUE, overwrite = TRUE, copy_mode = "copy")

  expect_true(file.exists(file.path(root, "ordered_results", "00_ordered_output_index.csv")))
  expect_true(any(grepl("1.1.01_Y_community_matrix_heatmap[.]png$", catalog$ordered_relpath)))
  expect_true(any(grepl("1.1.01_Y_community_matrix_heatmap[.]pdf$", catalog$ordered_relpath)))
  expect_true(file.exists(file.path(root, "tables", "ordered_output_index.csv")))
  expect_equal(
    unique(catalog$order_number[basename(catalog$original_relpath) %in% c("01_Y_heatmap.png", "01_Y_heatmap.pdf")]),
    "1.1.1"
  )
})
