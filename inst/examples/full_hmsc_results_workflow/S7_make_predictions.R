# S7_make_predictions.R -------------------------------------------------------
# Input: fitted HMSC models.
# Output: gradient predictions, plots, and a compact result index.

library(Hmsc)

if (!exists("localDir")) localDir <- file.path(getwd(), "hmscHist_full_results_example")
modelDir <- file.path(localDir, "models")
resultDir <- file.path(localDir, "results")
plotDir <- file.path(resultDir, "plots")
tableDir <- file.path(resultDir, "tables")
for (d in c(resultDir, plotDir, tableDir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

model_file <- list.files(modelDir, pattern = "models_thin_.*\\.RData$", full.names = TRUE)
if (length(model_file) == 0) {
  stop("No fitted model file found. Please run S2 first.", call. = FALSE)
}
model_file <- model_file[length(model_file)]
load(model_file)
if (!exists("models")) stop("Model file does not contain an object named 'models'.", call. = FALSE)

choose_covariates <- function(m) {
  xn <- colnames(m$XData)
  if (is.null(xn)) return(character())
  priority <- c(
    "temp", "colonization_age", "historical_turnover_predictor",
    "prob_weighted_source_diversity", "precip", "soil_pH"
  )
  out <- intersect(priority, xn)
  if (length(out) < 3) out <- unique(c(out, xn))[seq_len(min(3, length(unique(c(out, xn)))))]
  out
}

safe_pdf <- function(file, expr) {
  grDevices::pdf(file)
  on.exit(grDevices::dev.off(), add = TRUE)
  try(capture.output(force(expr)), silent = TRUE)
}

for (nm in names(models)) {
  m <- models[[nm]]
  covariates <- choose_covariates(m)
  if (length(covariates) == 0) next

  ex.sp <- which.min(abs(colMeans(m$Y, na.rm = TRUE) - 0.5))
  if (length(ex.sp) == 0 || is.na(ex.sp)) ex.sp <- 1

  for (covariate in covariates) {
    Gradient <- try(constructGradient(m, focalVariable = covariate), silent = TRUE)
    if (inherits(Gradient, "try-error")) next
    predY <- try(predict(m, Gradient = Gradient, expected = TRUE), silent = TRUE)
    if (inherits(predY, "try-error")) next

    saveRDS(Gradient, file.path(resultDir, paste0("gradient_", nm, "_", covariate, ".rds")))
    saveRDS(predY, file.path(resultDir, paste0("gradient_prediction_", nm, "_", covariate, ".rds")))

    safe_pdf(
      file.path(plotDir, paste0("gradient_richness_", nm, "_", covariate, ".pdf")),
      plotGradient(m, Gradient, pred = predY, measure = "S", showData = TRUE,
                   main = paste0(nm, ": richness over ", covariate))
    )

    safe_pdf(
      file.path(plotDir, paste0("gradient_species_", nm, "_", covariate, ".pdf")),
      plotGradient(m, Gradient, pred = predY, measure = "Y", index = ex.sp,
                   showData = TRUE, yshow = c(-0.1, 1.1),
                   main = paste0(nm, ": species ", ex.sp, " over ", covariate))
    )

    if (!is.null(m$Tr) && length(m$trNames) > 1) {
      safe_pdf(
        file.path(plotDir, paste0("gradient_trait_", nm, "_", covariate, ".pdf")),
        plotGradient(m, Gradient, pred = predY, measure = "T", index = 1,
                     showData = TRUE,
                     main = paste0(nm, ": community-weighted trait over ", covariate))
      )
    }
  }
}

describe_result_file <- function(files) {
  desc <- rep("Result file", length(files))
  desc[grepl("XData", files)] <- "HMSC XData table"
  desc[grepl("TrData", files)] <- "HMSC TrData table"
  desc[grepl("history", files)] <- "hmscHist historical variables or diagnostics"
  desc[grepl("convergence", files)] <- "MCMC convergence diagnostics"
  desc[grepl("traceplots", files)] <- "MCMC trace plots"
  desc[grepl("model_fit", files)] <- "Model fit metrics"
  desc[grepl("predicted", files)] <- "Predicted values or predicted richness"
  desc[grepl("variance_partitioning", files)] <- "Variance partitioning output"
  desc[grepl("association", files)] <- "Species association output"
  desc[grepl("plotBeta", files)] <- "Beta parameter heatmap"
  desc[grepl("plotGamma", files)] <- "Gamma parameter heatmap"
  desc[grepl("gradient", files)] <- "Gradient prediction output"
  desc
}

files <- list.files(resultDir, recursive = TRUE, full.names = FALSE)
index <- data.frame(
  file = files,
  description = describe_result_file(files),
  stringsAsFactors = FALSE
)
write.csv(index, file.path(resultDir, "ALL_RESULTS_INDEX.csv"), row.names = FALSE)

cat("S7 done: prediction files and ALL_RESULTS_INDEX.csv were written.\n")
