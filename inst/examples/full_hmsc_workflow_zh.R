# hmscHist 完整 HMSC 结果示例 -------------------------------------------------
# 这个脚本现在调用 full_results_workflow 中的 7 步完整流程。

this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
script_dir <- if (!is.na(this_file)) dirname(this_file) else "inst/examples"
source(file.path(script_dir, "full_results_workflow", "00_run_all_results_zh.R"), local = FALSE)
