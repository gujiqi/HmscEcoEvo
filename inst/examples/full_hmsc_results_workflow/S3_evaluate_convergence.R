# S3_evaluate_convergence.R ---------------------------------------------------
# Input: fitted HMSC models.
# Output: MCMC convergence tables and trace/histogram plots.

library(Hmsc)
if (!requireNamespace("coda", quietly = TRUE)) {
  stop("Package 'coda' is required. Install it with install.packages('coda').", call. = FALSE)
}
library(coda)

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

safe_mcmc_diag <- function(x, label) {
  if (is.null(x)) return(NULL)
  out <- try({
    ess <- coda::effectiveSize(x)
    psrf <- coda::gelman.diag(x, multivariate = FALSE)$psrf
    psrf_df <- as.data.frame(psrf)
    data.frame(
      parameter = names(ess),
      ESS = as.numeric(ess),
      PSRF_point = psrf_df[, 1],
      PSRF_upper = psrf_df[, 2],
      group = label,
      row.names = NULL
    )
  }, silent = TRUE)
  if (inherits(out, "try-error")) return(NULL)
  out
}

safe_pdf <- function(file, expr) {
  grDevices::pdf(file)
  on.exit(grDevices::dev.off(), add = TRUE)
  try(capture.output(force(expr)), silent = TRUE)
}

plot_mcmc_summary <- function(diag_df, file, title) {
  if (is.null(diag_df) || nrow(diag_df) == 0) return(invisible(NULL))
  safe_pdf(file, {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op), add = TRUE)
    hist(diag_df$ESS, main = paste0(title, ": ESS"), xlab = "Effective sample size")
    hist(diag_df$PSRF_point, main = paste0(title, ": PSRF"), xlab = "PSRF point estimate")
  })
}

for (nm in names(models)) {
  m <- models[[nm]]
  mpost <- convertToCodaObject(m)
  saveRDS(mpost, file.path(resultDir, paste0("posterior_coda_", nm, ".rds")))

  diag_list <- list(
    Beta = safe_mcmc_diag(mpost$Beta, "Beta"),
    Gamma = safe_mcmc_diag(mpost$Gamma, "Gamma"),
    Rho = safe_mcmc_diag(mpost$Rho, "Rho"),
    V = safe_mcmc_diag(mpost$V, "V")
  )
  if (!is.null(mpost$Omega) && length(mpost$Omega) > 0) {
    diag_list$Omega1 <- safe_mcmc_diag(mpost$Omega[[1]], "Omega_random_level_1")
  }
  if (!is.null(mpost$Alpha) && length(mpost$Alpha) > 0) {
    diag_list$Alpha1 <- safe_mcmc_diag(mpost$Alpha[[1]], "Alpha_random_level_1")
  }
  keep <- !vapply(diag_list, is.null, logical(1))
  diag_df <- if (any(keep)) do.call(rbind, diag_list[keep]) else NULL
  if (!is.null(diag_df)) {
    write.csv(diag_df, file.path(tableDir, paste0("convergence_diagnostics_", nm, ".csv")), row.names = FALSE)
    plot_mcmc_summary(diag_df, file.path(plotDir, paste0("convergence_histograms_", nm, ".pdf")), nm)
  }

  safe_pdf(file.path(plotDir, paste0("traceplots_Beta_", nm, ".pdf")), plot(mpost$Beta))

  if (!is.null(mpost$Gamma)) {
    safe_pdf(file.path(plotDir, paste0("traceplots_Gamma_", nm, ".pdf")), plot(mpost$Gamma))
  }
  if (!is.null(mpost$Omega) && length(mpost$Omega) > 0) {
    safe_pdf(file.path(plotDir, paste0("traceplots_Omega1_", nm, ".pdf")), plot(mpost$Omega[[1]]))
  }
}

cat("S3 done: MCMC convergence diagnostics were written.\n")
