# S2_fit_models_short.R -------------------------------------------------------
# Input: models/unfitted_models.RData.
# Output: models/models_thin_...RData.

set.seed(1)
library(Hmsc)

if (!exists("localDir")) localDir <- file.path(getwd(), "hmscHist_full_results_example")
modelDir <- file.path(localDir, "models")
resultDir <- file.path(localDir, "results")
if (!dir.exists(modelDir)) dir.create(modelDir, recursive = TRUE)
if (!dir.exists(resultDir)) dir.create(resultDir, recursive = TRUE)

load(file.path(modelDir, "unfitted_models.RData"))

get_mcmc_value <- function(settings, name, default) {
  value <- settings[[name]]
  if (is.null(value)) default else value
}

settings <- getOption("hmscHist.mcmc_settings", list())
thin <- get_mcmc_value(settings, "thin", 1)
samples <- get_mcmc_value(settings, "samples", 30)
transient <- get_mcmc_value(settings, "transient", 15)
nChains <- get_mcmc_value(settings, "nChains", 2)
nParallel <- get_mcmc_value(settings, "nParallel", 1)
verbose <- get_mcmc_value(settings, "verbose", 0)

cat("MCMC settings:",
    "thin =", thin,
    "samples =", samples,
    "transient =", transient,
    "nChains =", nChains,
    "nParallel =", nParallel,
    "\n")

computational_time <- list()
fitted_models <- models
for (nm in names(models)) {
  cat("Fitting model:", nm, "\n")
  ptm <- proc.time()
  fitted_models[[nm]] <- sampleMcmc(
    models[[nm]],
    samples = samples,
    transient = transient,
    thin = thin,
    nChains = nChains,
    nParallel = nParallel,
    verbose = verbose
  )
  computational_time[[nm]] <- proc.time() - ptm
}
models <- fitted_models

filename <- file.path(modelDir, paste0("models_thin_", thin, "_samples_", samples, "_chains_", nChains, ".RData"))
save(models, computational_time, thin, samples, transient, nChains, nParallel, file = filename)
writeLines(c(
  paste("thin =", thin),
  paste("samples =", samples),
  paste("transient =", transient),
  paste("nChains =", nChains),
  paste("nParallel =", nParallel),
  "NOTE: Increase samples/transient/thin/nChains for formal inference."
), con = file.path(resultDir, "mcmc_settings.txt"))

cat("S2 done: fitted model file:", filename, "\n")
