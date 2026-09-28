# S5_show_model_fit.R ---------------------------------------------------------
# 目标：把模型拟合结果整理成直观表格和 PDF 图。

if (!exists("script_dir") || !file.exists(file.path(script_dir, "_common.R"))) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
  script_dir <- if (!is.na(this_file)) dirname(this_file) else "inst/examples/full_results_workflow"
}
source(file.path(script_dir, "_common.R"))

load(file.path(dataDir, "analysis_data.RData"))
fit_results <- readRDS(file.path(rdsDir, "04_fit_results.rds"))

for (nm in names(fit_results)) {
  MF_expl_df <- as_metric_table(fit_results[[nm]]$MF_expl)
  MF_cv_df <- as_metric_table(fit_results[[nm]]$MF_cv)

  pdf(file.path(figureDir, paste0("05_model_fit_summary_", nm, ".pdf")), width = 10, height = 8)
  par(mfrow = c(2, 2), mar = c(6, 4, 3, 1))

  if ("AUC" %in% names(MF_expl_df)) {
    barplot(MF_expl_df$AUC, names.arg = MF_expl_df$species, las = 2,
            main = "Explanatory AUC", ylab = "AUC")
  }
  if ("TjurR2" %in% names(MF_expl_df)) {
    barplot(MF_expl_df$TjurR2, names.arg = MF_expl_df$species, las = 2,
            main = "Explanatory Tjur R2", ylab = "Tjur R2")
  }
  if ("AUC" %in% names(MF_cv_df)) {
    barplot(MF_cv_df$AUC, names.arg = MF_cv_df$species, las = 2,
            main = "Cross-validation AUC", ylab = "AUC")
  }
  if ("TjurR2" %in% names(MF_cv_df)) {
    barplot(MF_cv_df$TjurR2, names.arg = MF_cv_df$species, las = 2,
            main = "Cross-validation Tjur R2", ylab = "Tjur R2")
  }
  dev.off()

  pred_mean <- prediction_mean_matrix(fit_results[[nm]]$pred_expl, Y)
  if (!is.null(pred_mean)) {
    obs_rich <- rowSums(Y, na.rm = TRUE)
    pred_rich <- rowSums(pred_mean, na.rm = TRUE)
    rich_df <- data.frame(site = rownames(Y), observed_richness = obs_rich, predicted_richness = pred_rich)
    write.csv(rich_df, file.path(tableDir, paste0("05_observed_vs_predicted_richness_", nm, ".csv")), row.names = FALSE)

    pdf(file.path(figureDir, paste0("05_observed_vs_predicted_richness_", nm, ".pdf")), width = 6, height = 6)
    plot(obs_rich, pred_rich, pch = 19, xlab = "Observed richness", ylab = "Predicted richness",
         main = paste(nm, "richness"))
    abline(0, 1, lty = 2)
    dev.off()
  }
}

cat("S5 done. Model fit figures saved.\n")
