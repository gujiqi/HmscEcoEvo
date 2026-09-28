# S1_define_real_akodon_model.R ----------------------------------------------
# Input: real Akodon data in inst/extdata/real_akodon.
# Output: data/prepared_data.RData and models/unfitted_models.RData.

set.seed(1)

if (!requireNamespace("Hmsc", quietly = TRUE)) {
  stop("Please install Hmsc first.", call. = FALSE)
}
if (!requireNamespace("HmscEcoEvo", quietly = TRUE)) {
  stop("Please install or load hmscHist first.", call. = FALSE)
}

library(Hmsc)
library(HmscEcoEvo)

if (!exists("localDir")) localDir <- file.path(getwd(), "hmscHist_real_akodon_full_results")
dataOutDir <- file.path(localDir, "data")
modelDir <- file.path(localDir, "models")
resultDir <- file.path(localDir, "results")
plotDir <- file.path(resultDir, "plots")
tableDir <- file.path(resultDir, "tables")
for (d in c(dataOutDir, modelDir, resultDir, plotDir, tableDir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

dataDir <- system.file("extdata", "real_akodon", package = "HmscEcoEvo")
if (!nzchar(dataDir)) {
  dataDir <- file.path("inst", "extdata", "real_akodon")
}
if (!dir.exists(dataDir)) {
  stop("Cannot find real Akodon data directory: ", dataDir, call. = FALSE)
}

comm <- read.csv(file.path(dataDir, "comm.csv"), row.names = 1, check.names = FALSE)
Y <- as.matrix(comm)
storage.mode(Y) <- "numeric"

site_region <- read.csv(file.path(dataDir, "site_region.csv"), stringsAsFactors = FALSE)
env <- read.csv(file.path(dataDir, "env.csv"), row.names = 1, check.names = FALSE)
coords <- read.csv(file.path(dataDir, "coords.csv"), row.names = 1, check.names = FALSE)
traits <- read.csv(file.path(dataDir, "traits.csv"), row.names = 1, check.names = FALSE)
species_region_history <- read.csv(file.path(dataDir, "species_region_history.csv"), stringsAsFactors = FALSE)
trait_history <- read.csv(file.path(dataDir, "trait_history.csv"), stringsAsFactors = FALSE)

phy <- NULL
phy_file <- file.path(dataDir, "phylo.tre")
if (requireNamespace("ape", quietly = TRUE) && file.exists(phy_file)) {
  phy <- ape::read.tree(phy_file)
}

proj <- hmscHist_data(
  comm = Y,
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
hist_diag <- diagnose_history_indices(hist, env = env)

strategy <- "minimal"
XData <- as_hmsc_xdata(hist, env = env, strategy = strategy)
TrData <- as_hmsc_trdata(hist, traits = traits, strategy = strategy)
random_hist <- as_hmsc_random(hist, strategy = strategy)

XFormula <- stats::as.formula("~ .")
TrFormula <- if (!is.null(TrData) && ncol(TrData) > 0) stats::as.formula("~ .") else NULL

studyDesign <- data.frame(sample = factor(rownames(Y)), row.names = rownames(Y))
ranLevels <- list(sample = HmscRandomLevel(units = levels(studyDesign$sample)))

if (!is.null(random_hist) && ncol(random_hist) > 0) {
  studyDesign_with_history <- cbind(studyDesign, random_hist[rownames(studyDesign), , drop = FALSE])
} else {
  studyDesign_with_history <- studyDesign
}

make_hmsc <- function(args) {
  hmsc_formals <- names(formals(Hmsc::Hmsc))
  y_arg <- if ("YData" %in% hmsc_formals) "YData" else if ("Y" %in% hmsc_formals) "Y" else NA_character_
  if (is.na(y_arg)) stop("Installed Hmsc::Hmsc() does not expose a Y/YData argument.", call. = FALSE)
  args[[y_arg]] <- Y
  do.call(Hmsc, args)
}

model_notes <- character()
base_args <- list(
  XData = XData,
  XFormula = XFormula,
  TrData = TrData,
  TrFormula = TrFormula,
  phyloTree = phy,
  studyDesign = studyDesign,
  ranLevels = ranLevels,
  distr = "probit"
)

model0 <- try(make_hmsc(base_args), silent = TRUE)

if (inherits(model0, "try-error")) {
  model_notes <- c(model_notes, "Full model with phyloTree/TrData/sample random effect failed; trying without phyloTree.")
  args2 <- base_args
  args2$phyloTree <- NULL
  model0 <- try(make_hmsc(args2), silent = TRUE)
}

if (inherits(model0, "try-error")) {
  model_notes <- c(model_notes, "Model with TrData/sample random effect failed; trying XData + sample random effect only.")
  args3 <- list(
    XData = XData,
    XFormula = XFormula,
    studyDesign = studyDesign,
    ranLevels = ranLevels,
    distr = "probit"
  )
  model0 <- try(make_hmsc(args3), silent = TRUE)
}

if (inherits(model0, "try-error")) {
  model_notes <- c(model_notes, "Model with sample random effect failed; trying fixed-effects-only model.")
  args4 <- list(
    XData = XData,
    XFormula = XFormula,
    distr = "probit"
  )
  model0 <- make_hmsc(args4)
}

models <- list(real_akodon = model0)

write.csv(Y, file.path(tableDir, "01_Y_comm.csv"))
write.csv(XData, file.path(tableDir, "02_XData_for_HMSC.csv"))
if (!is.null(TrData)) write.csv(TrData, file.path(tableDir, "03_TrData_for_HMSC.csv"))
if (!is.null(random_hist)) write.csv(random_hist, file.path(tableDir, "04_random_history_candidates.csv"))
write.csv(hist$XData_history, file.path(tableDir, "05_all_XData_history.csv"))
if (!is.null(hist$TrData_history)) write.csv(hist$TrData_history, file.path(tableDir, "06_all_TrData_history.csv"))
write.csv(coords, file.path(tableDir, "07_coords.csv"))
write.csv(site_region, file.path(tableDir, "08_site_region.csv"), row.names = FALSE)

capture.output({
  print(hist)
  cat("\nGroups:\n")
  print(hist$groups)
  cat("\nDiagnostics:\n")
  print(hist_diag)
}, file = file.path(resultDir, "history_indices_summary.txt"))
capture.output(hist_diag, file = file.path(resultDir, "history_indices_diagnostics.txt"))
writeLines(model_notes, con = file.path(resultDir, "model_setup_notes.txt"))

save(Y, XData, TrData, random_hist, hist, hist_diag, env, traits, coords, phy,
     studyDesign, studyDesign_with_history, ranLevels, XFormula, TrFormula,
     file = file.path(dataOutDir, "prepared_data.RData"))
save(models, file = file.path(modelDir, "unfitted_models.RData"))

cat("S1 done: real Akodon model defined, inputs written, model not yet fitted.\n")
