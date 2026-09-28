# S1_define_models_hmscHist.R -------------------------------------------------
# 输入：包内 toy 数据
# 输出：models/unfitted_models.RData 和 data/prepared_data.RData

set.seed(1)

if (!requireNamespace("Hmsc", quietly = TRUE)) stop("请先安装 Hmsc：install.packages('Hmsc')")
if (!requireNamespace("HmscEcoEvo", quietly = TRUE)) stop("请先安装或 load_all hmscHist")

library(Hmsc)
library(HmscEcoEvo)

# localDir 由 00_run_all_results_zh.R 创建；单独运行时也可自动创建
if (!exists("localDir")) localDir <- file.path(getwd(), "hmscHist_full_results_example")
dataDir <- file.path(localDir, "data")
modelDir <- file.path(localDir, "models")
resultDir <- file.path(localDir, "results")
plotDir <- file.path(resultDir, "plots")
tableDir <- file.path(resultDir, "tables")
for (d in c(dataDir, modelDir, resultDir, plotDir, tableDir)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

pkg_file <- function(...) {
  p <- system.file(..., package = "HmscEcoEvo")
  if (nzchar(p)) return(p)
  file.path("inst", ...)
}

toyDir <- pkg_file("extdata", "full_hmsc")

# 1. 读取数据 -----------------------------------------------------------------
comm <- read.csv(file.path(toyDir, "comm.csv"), row.names = 1, check.names = FALSE)
Y <- as.matrix(comm)
storage.mode(Y) <- "numeric"

site_region <- read.csv(file.path(toyDir, "site_region.csv"), stringsAsFactors = FALSE)
env <- read.csv(file.path(toyDir, "env.csv"), row.names = 1, check.names = FALSE)
coords <- read.csv(file.path(toyDir, "coords.csv"), row.names = 1, check.names = FALSE)
traits <- read.csv(file.path(toyDir, "traits.csv"), row.names = 1, check.names = FALSE)

species_region_history <- read.csv(file.path(toyDir, "species_region_history.csv"), stringsAsFactors = FALSE)
speciation_events <- read.csv(file.path(toyDir, "speciation_events.csv"), stringsAsFactors = FALSE)
history_events <- read.csv(file.path(toyDir, "history_events.csv"), stringsAsFactors = FALSE)
trait_history <- read.csv(file.path(toyDir, "trait_history.csv"), stringsAsFactors = FALSE)

phy <- NULL
phy_file <- file.path(toyDir, "phylo.tre")
if (requireNamespace("ape", quietly = TRUE) && file.exists(phy_file)) {
  phy <- ape::read.tree(phy_file)
}

# 2. hmscHist 历史指数 ---------------------------------------------------------
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
  speciation_events = speciation_events,
  history_events = history_events,
  trait_history = trait_history,
  method = "manual"
)

hist <- calc_history_indices(proj)
hist_diag <- diagnose_history_indices(hist, env = env)

# 3. 选择 HMSC 输入 ------------------------------------------------------------
strategy <- "minimal"
XData <- as_hmsc_xdata(hist, env = env, strategy = strategy)
TrData <- as_hmsc_trdata(hist, traits = traits, strategy = strategy)
random_hist <- as_hmsc_random(hist, strategy = strategy)

# 因为 toy 数据很小，XData 中变量不能太多。minimal 一般最稳。
XFormula <- stats::as.formula("~ .")
TrFormula <- if (!is.null(TrData) && ncol(TrData) > 0) stats::as.formula("~ .") else NULL

# 4. studyDesign 和随机效应 ---------------------------------------------------
# 为了展示 computeAssociations，加入 sample-level random effect。
studyDesign <- data.frame(sample = factor(rownames(Y)), row.names = rownames(Y))
ranLevels <- list(sample = HmscRandomLevel(units = levels(studyDesign$sample)))

# 历史随机效应候选变量另外保存；正式分析时可把 dominant_source_region 加入 studyDesign。
if (!is.null(random_hist) && ncol(random_hist) > 0) {
  studyDesign_with_history <- cbind(studyDesign, random_hist[rownames(studyDesign), , drop = FALSE])
} else {
  studyDesign_with_history <- studyDesign
}

# 5. 定义 HMSC 模型 ------------------------------------------------------------
# 优先使用 phyloTree + TrData；如果失败，则退回简化模型。
model_notes <- character()
model0 <- try(
  Hmsc(
    Y = Y,
    XData = XData,
    XFormula = XFormula,
    TrData = TrData,
    TrFormula = TrFormula,
    phyloTree = phy,
    studyDesign = studyDesign,
    ranLevels = ranLevels,
    distr = "probit"
  ),
  silent = TRUE
)

if (inherits(model0, "try-error")) {
  model_notes <- c(model_notes, "Full model with phyloTree/TrData/sample random effect failed; trying without phyloTree.")
  model0 <- try(
    Hmsc(
      Y = Y,
      XData = XData,
      XFormula = XFormula,
      TrData = TrData,
      TrFormula = TrFormula,
      studyDesign = studyDesign,
      ranLevels = ranLevels,
      distr = "probit"
    ),
    silent = TRUE
  )
}

if (inherits(model0, "try-error")) {
  model_notes <- c(model_notes, "Model with TrData/sample random effect failed; trying XData + sample random effect only.")
  model0 <- try(
    Hmsc(
      Y = Y,
      XData = XData,
      XFormula = XFormula,
      studyDesign = studyDesign,
      ranLevels = ranLevels,
      distr = "probit"
    ),
    silent = TRUE
  )
}

if (inherits(model0, "try-error")) {
  model_notes <- c(model_notes, "Model with sample random effect failed; trying fixed-effects-only model.")
  model0 <- Hmsc(
    Y = Y,
    XData = XData,
    XFormula = XFormula,
    distr = "probit"
  )
}

models <- list(main = model0)

# 6. 保存 ---------------------------------------------------------------------
write.csv(Y, file.path(tableDir, "01_Y_comm.csv"))
write.csv(XData, file.path(tableDir, "02_XData_for_HMSC.csv"))
if (!is.null(TrData)) write.csv(TrData, file.path(tableDir, "03_TrData_for_HMSC.csv"))
if (!is.null(random_hist)) write.csv(random_hist, file.path(tableDir, "04_random_history_candidates.csv"))
write.csv(hist$XData_history, file.path(tableDir, "05_all_XData_history.csv"))
if (!is.null(hist$TrData_history)) write.csv(hist$TrData_history, file.path(tableDir, "06_all_TrData_history.csv"))

capture.output(summary(hist), file = file.path(resultDir, "history_indices_summary.txt"))
capture.output(hist_diag, file = file.path(resultDir, "history_indices_diagnostics.txt"))
writeLines(model_notes, con = file.path(resultDir, "model_setup_notes.txt"))

save(Y, XData, TrData, random_hist, hist, hist_diag, env, traits, coords, phy, studyDesign, studyDesign_with_history,
     ranLevels, XFormula, TrFormula, file = file.path(dataDir, "prepared_data.RData"))
save(models, file = file.path(modelDir, "unfitted_models.RData"))

cat("S1 完成：模型已定义但尚未拟合。\n")
