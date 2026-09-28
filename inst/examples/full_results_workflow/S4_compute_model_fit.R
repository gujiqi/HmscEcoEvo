# S4_compute_model_fit.R ------------------------------------------------------
# Input: fitted HMSC models.
# Output: explanatory fit, cross-validation fit, WAIC, and predicted values.

library(Hmsc)

if (!exists("localDir")) localDir <- file.path(getwd(), "hmscHist_full_results_example")
modelDir <- file.path(localDir, "models")
resultDir <- file.path(localDir, "results")
tableDir <- file.path(resultDir, "tables")
for (d in c(resultDir, tableDir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

model_file <- list.files(modelDir, pattern = "models_thin_.*\\.RData$", full.names = TRUE)
if (length(model_file) == 0) {
  stop("No fitted model file found. Please run S2 first.", call. = FALSE)
}
model_file <- model_file[length(model_file)]
load(model_file)

fit_to_df <- function(fit) {
  if (is.null(fit) || length(fit) == 0) return(data.frame())
  maxlen <- max(vapply(fit, length, integer(1)))
  out <- data.frame(species_index = seq_len(maxlen))
  for (nm in names(fit)) {
    v <- fit[[nm]]
    out[[nm]] <- rep(NA_real_, maxlen)
    out[[nm]][seq_along(v)] <- as.numeric(v)
  }
  out
}

pred_mean_matrix <- function(pred, Y) {
  d <- dim(pred)
  if (length(d) == 2) return(pred)
  if (length(d) != 3) return(NULL)
  ny <- nrow(Y)
  ns <- ncol(Y)
  if (d[1] == ny && d[2] == ns) return(apply(pred, c(1, 2), mean, na.rm = TRUE))
  if (d[2] == ny && d[3] == ns) return(apply(pred, c(2, 3), mean, na.rm = TRUE))
  if (d[1] == ny && d[3] == ns) return(apply(pred, c(1, 3), mean, na.rm = TRUE))
  NULL
}

for (nm in names(models)) {
  m <- models[[nm]]
  predY <- computePredictedValues(m, expected = TRUE)
  MF <- evaluateModelFit(hM = m, predY = predY)
  saveRDS(predY, file.path(resultDir, paste0("predicted_values_array_", nm, ".rds")))
  saveRDS(MF, file.path(resultDir, paste0("model_fit_explanatory_", nm, ".rds")))
  write.csv(fit_to_df(MF), file.path(tableDir, paste0("model_fit_explanatory_", nm, ".csv")), row.names = FALSE)

  pm <- pred_mean_matrix(predY, m$Y)
  if (!is.null(pm)) {
    rownames(pm) <- rownames(m$Y)
    colnames(pm) <- colnames(m$Y)
    write.csv(pm, file.path(tableDir, paste0("predicted_probability_mean_", nm, ".csv")))
    richness <- data.frame(
      site = rownames(pm),
      observed_richness = rowSums(m$Y, na.rm = TRUE),
      predicted_richness = rowSums(pm, na.rm = TRUE)
    )
    write.csv(richness, file.path(tableDir, paste0("predicted_richness_", nm, ".csv")), row.names = FALSE)
  }

  cvMF <- try({
    if (!is.null(m$studyDesign) && "sample" %in% names(m$studyDesign)) {
      partition <- createPartition(m, nfolds = 2, column = "sample")
    } else {
      partition <- createPartition(m, nfolds = 2)
    }
    cvpred <- computePredictedValues(m, partition = partition, nParallel = 1, expected = TRUE)
    evaluateModelFit(hM = m, predY = cvpred)
  }, silent = TRUE)
  if (!inherits(cvMF, "try-error")) {
    saveRDS(cvMF, file.path(resultDir, paste0("model_fit_cross_validation_", nm, ".rds")))
    write.csv(fit_to_df(cvMF), file.path(tableDir, paste0("model_fit_cross_validation_", nm, ".csv")), row.names = FALSE)
  } else {
    writeLines(as.character(cvMF), con = file.path(resultDir, paste0("cross_validation_failed_", nm, ".txt")))
  }

  WAIC <- try(computeWAIC(m), silent = TRUE)
  if (!inherits(WAIC, "try-error")) {
    write.csv(data.frame(model = nm, WAIC = as.numeric(WAIC)), file.path(tableDir, paste0("WAIC_", nm, ".csv")), row.names = FALSE)
  }
}

cat("S4 done: model fit, cross-validation, WAIC, and predictions were written.\n")
