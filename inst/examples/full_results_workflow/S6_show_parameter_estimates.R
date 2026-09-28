# S6_show_parameter_estimates.R ----------------------------------------------
# Input: fitted HMSC models.
# Output: Beta/Gamma/Rho posterior estimate tables and plots.

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

safe_pdf <- function(file, expr) {
  grDevices::pdf(file)
  on.exit(grDevices::dev.off(), add = TRUE)
  try(capture.output(force(expr)), silent = TRUE)
}

write_post_components <- function(post, prefix) {
  if (is.null(post)) return(invisible(NULL))
  saveRDS(post, paste0(prefix, ".rds"))
  for (nm in names(post)) {
    x <- post[[nm]]
    if (is.matrix(x) || is.data.frame(x)) {
      write.csv(x, paste0(prefix, "_", nm, ".csv"))
    } else if (is.array(x)) {
      saveRDS(x, paste0(prefix, "_", nm, ".rds"))
    } else if (is.atomic(x)) {
      write.csv(data.frame(value = as.vector(x)), paste0(prefix, "_", nm, ".csv"), row.names = FALSE)
    }
  }
}

for (nm in names(models)) {
  m <- models[[nm]]

  betaPost <- try(getPostEstimate(m, parName = "Beta"), silent = TRUE)
  if (!inherits(betaPost, "try-error")) {
    write_post_components(betaPost, file.path(tableDir, paste0("post_Beta_", nm)))
    safe_pdf(file.path(plotDir, paste0("plotBeta_support_", nm, ".pdf")),
             plotBeta(m, post = betaPost, param = "Support", supportLevel = 0.95))
    safe_pdf(file.path(plotDir, paste0("plotBeta_mean_", nm, ".pdf")),
             plotBeta(m, post = betaPost, param = "Mean", plotTree = !is.null(m$phyloTree)))
  }

  gammaPost <- try(getPostEstimate(m, parName = "Gamma"), silent = TRUE)
  if (!inherits(gammaPost, "try-error")) {
    write_post_components(gammaPost, file.path(tableDir, paste0("post_Gamma_", nm)))
    safe_pdf(file.path(plotDir, paste0("plotGamma_support_", nm, ".pdf")),
             plotGamma(m, post = gammaPost, param = "Support", supportLevel = 0.95))
    safe_pdf(file.path(plotDir, paste0("plotGamma_mean_", nm, ".pdf")),
             plotGamma(m, post = gammaPost, param = "Mean"))
  }

  rhoPost <- try(getPostEstimate(m, parName = "Rho"), silent = TRUE)
  if (!inherits(rhoPost, "try-error")) {
    write_post_components(rhoPost, file.path(tableDir, paste0("post_Rho_", nm)))
  }
}

cat("S6 done: posterior estimate tables and plots were written.\n")
