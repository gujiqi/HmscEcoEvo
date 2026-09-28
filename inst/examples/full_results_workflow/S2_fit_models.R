# S2_fit_models.R -------------------------------------------------------------
# Input: models/unfitted_models.RData.
# Output: models/models_thin_...RData.

if (!exists("script_dir") || !file.exists(file.path(script_dir, "_common.R"))) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
  script_dir <- if (!is.na(this_file)) dirname(this_file) else "inst/examples/full_results_workflow"
}
source(file.path(script_dir, "_common.R"))

set.seed(1)
library(Hmsc)

load(file.path(modelDir, "unfitted_models.RData"))

get_mcmc_value <- function(settings, name, default) {
  value <- settings[[name]]
  if (is.null(value)) default else value
}

settings <- getOption("hmscHist.mcmc_settings", list())
thin <- get_mcmc_value(settings, "thin", 10)
samples <- get_mcmc_value(settings, "samples", 1000)
transient <- get_mcmc_value(settings, "transient", 1000)
nChains <- get_mcmc_value(settings, "nChains", 4)
nParallel <- get_mcmc_value(settings, "nParallel", 1)
verbose <- get_mcmc_value(settings, "verbose", 100)

mcmc_settings <- data.frame(
  samples = samples,
  transient = transient,
  thin = thin,
  nChains = nChains,
  nParallel = nParallel,
  verbose = verbose
)
write.csv(mcmc_settings, file.path(tableDir, "02_mcmc_settings.csv"), row.names = FALSE)

cat("MCMC settings:",
    "thin =", thin,
    "samples =", samples,
    "transient =", transient,
    "nChains =", nChains,
    "nParallel =", nParallel,
    "\n")

fitted_models <- vector("list", length(models))
names(fitted_models) <- names(models)
computational_time <- numeric(length(models))
names(computational_time) <- names(models)

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
  computational_time[nm] <- (proc.time() - ptm)[3]
}

models <- fitted_models

out_file <- file.path(
  modelDir,
  paste0("models_thin_", thin, "_samples_", samples, "_chains_", nChains, ".RData")
)
save(models, computational_time, mcmc_settings, file = out_file)
write.csv(data.frame(model = names(computational_time), seconds = computational_time),
          file.path(tableDir, "02_computational_time.csv"), row.names = FALSE)
writeLines(c(
  paste("thin =", thin),
  paste("samples =", samples),
  paste("transient =", transient),
  paste("nChains =", nChains),
  paste("nParallel =", nParallel),
  "NOTE: Increase samples/transient/thin/nChains for formal inference."
), con = file.path(resultDir, "mcmc_settings.txt"))

cat("S2 done. Fitted models saved to:\n")
cat(out_file, "\n")
