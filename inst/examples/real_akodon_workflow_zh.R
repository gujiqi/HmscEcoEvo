# Real Akodon hmscHist example
# This script uses real Akodon occurrence, phylogeny, geography and trait data.
# It does not require speciation_events or history_events.

library(HmscEcoEvo)

# If you want to fit HMSC too, install Hmsc and ape:
# install.packages(c("Hmsc", "ape", "coda", "corrplot"))
# 3. 检查是不是新版
history_dictionary()$variable
setwd("C:/Users/Google/Documents/HMSC-HIST/hmscHist_0.1.4_real_akodon_source/hmscHist")


unlink("real_akodon_hmscHist_outputs", recursive = TRUE)

source("inst/examples/real_akodon_workflow_zh.R")
data_dir <- "inst/extdata/real_akodon"
if (!dir.exists(data_dir)) {
  # If running from this data folder directly:
  data_dir <- "."
}

comm <- read.csv(file.path(data_dir, "comm.csv"), row.names = 1, check.names = FALSE)
site_region <- read.csv(file.path(data_dir, "site_region.csv"), stringsAsFactors = FALSE)
env <- read.csv(file.path(data_dir, "env.csv"), row.names = 1, check.names = FALSE)
coords <- read.csv(file.path(data_dir, "coords.csv"), row.names = 1, check.names = FALSE)
traits <- read.csv(file.path(data_dir, "traits.csv"), row.names = 1, check.names = FALSE)
species_region_history <- read.csv(file.path(data_dir, "species_region_history.csv"), stringsAsFactors = FALSE)
trait_history <- read.csv(file.path(data_dir, "trait_history.csv"), stringsAsFactors = FALSE)

phy <- NULL
if (requireNamespace("ape", quietly = TRUE)) {
  phy <- ape::read.tree(file.path(data_dir, "phylo.tre"))
}

proj <- hmscHist_data(
  comm = comm,
  site_region = site_region,
  env = env,
  coords = coords,
  traits = traits,
  phy = phy
)

proj <- build_history_tables(
  proj,
  species_region_history = species_region_history,
  trait_history = trait_history,
  method = "manual"
)

hist <- calc_history_indices(proj)
summary(hist)

XData <- as_hmsc_xdata(hist, env = env, strategy = "minimal")
TrData <- as_hmsc_trdata(hist, traits = traits, strategy = "minimal")
random_hist <- as_hmsc_random(hist, strategy = "minimal")

diag <- diagnose_history_indices(hist, env = env)
print(diag$warnings)

# Save hmscHist outputs
out_dir <- "real_akodon_hmscHist_outputs"
dir.create(out_dir, showWarnings = FALSE)
write.csv(hist$XData_history, file.path(out_dir, "XData_history.csv"))
if (!is.null(hist$TrData_history)) write.csv(hist$TrData_history, file.path(out_dir, "TrData_history.csv"))
write.csv(random_hist, file.path(out_dir, "random_history.csv"))
write.csv(XData, file.path(out_dir, "XData_for_HMSC.csv"))
if (!is.null(TrData)) write.csv(TrData, file.path(out_dir, "TrData_for_HMSC.csv"))

# Optional short HMSC run
if (requireNamespace("Hmsc", quietly = TRUE)) {
  hmsc_formals <- names(formals(Hmsc::Hmsc))
  y_arg <- if ("YData" %in% hmsc_formals) "YData" else if ("Y" %in% hmsc_formals) "Y" else NA_character_
  if (is.na(y_arg)) {
    message("Skipping optional HMSC fit: installed Hmsc::Hmsc() does not expose a Y/YData argument.")
  } else {
    library(Hmsc)
    Y <- as.matrix(comm)
    XFormula <- ~ .
    TrFormula <- if (!is.null(TrData) && ncol(TrData) > 0) ~ . else ~ 1
    studyDesign <- data.frame(sample = factor(rownames(Y)))
    rL <- HmscRandomLevel(units = levels(studyDesign$sample))
    hmsc_args <- list(
      XData = XData,
      XFormula = XFormula,
      TrData = TrData,
      TrFormula = TrFormula,
      phyloTree = phy,
      studyDesign = studyDesign,
      ranLevels = list(sample = rL),
      distr = "probit"
    )
    hmsc_args[[y_arg]] <- Y
    fit_result <- tryCatch({
      model0 <- do.call(Hmsc, hmsc_args)
      model <- sampleMcmc(model0, samples = 10, transient = 5, thin = 1,
                          nChains = 2, nParallel = 1, verbose = 0)
      pred <- computePredictedValues(model)
      fit <- evaluateModelFit(hM = model, predY = pred)
      saveRDS(model, file.path(out_dir, "toy_hmsc_model.rds"))
      saveRDS(fit, file.path(out_dir, "toy_hmsc_fit.rds"))
      fit
    }, error = function(e) {
      message("Skipping optional HMSC fit after Hmsc error: ", conditionMessage(e))
      NULL
    })
    if (!is.null(fit_result)) print(fit_result)
  }
}
