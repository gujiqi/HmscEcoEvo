# hmscHist 快速示例 ---------------------------------------------------------
# 这个示例只需要一个历史主表 species_region_history.csv。
# 不需要准备各种事件表。事件表是可选扩展。

library(HmscEcoEvo)

dat <- system.file("extdata/basic", package = "HmscEcoEvo")

comm <- read.csv(file.path(dat, "comm.csv"), row.names = 1, check.names = FALSE)
site_region <- read.csv(file.path(dat, "site_region.csv"), stringsAsFactors = FALSE)
env <- read.csv(file.path(dat, "env.csv"), row.names = 1, check.names = FALSE)
traits <- read.csv(file.path(dat, "traits.csv"), row.names = 1, check.names = FALSE)
species_region_history <- read.csv(file.path(dat, "species_region_history.csv"), stringsAsFactors = FALSE)
trait_history <- read.csv(file.path(dat, "trait_history.csv"), stringsAsFactors = FALSE)

# 1. 创建项目对象
proj <- hmscHist_data(
  comm = comm,
  site_region = site_region,
  env = env,
  traits = traits
)

# 2. 构建/标准化历史表
#    这里只给 species_region_history 和 trait_history。
#    没有 speciation_events 和 history_events 也可以跑。
proj <- build_history_tables(
  proj,
  species_region_history = species_region_history,
  trait_history = trait_history,
  method = "manual"
)

summary(proj)

# 3. 计算历史指数
hist <- calc_history_indices(proj)
summary(hist)

# 4. 选择最少、最代表性的历史变量
X_min <- select_history_variables(hist, strategy = "minimal", component = "XData")
Tr_min <- select_history_variables(hist, strategy = "minimal", component = "TrData")
Ran_min <- select_history_variables(hist, strategy = "minimal", component = "random")

X_min
Tr_min
Ran_min

# 5. 生成 HMSC 可用输入
XData <- as_hmsc_xdata(hist, env = env, strategy = "minimal")
TrData <- as_hmsc_trdata(hist, traits = traits, strategy = "minimal")
random_hist <- as_hmsc_random(hist, strategy = "minimal")

# 6. 诊断
#    查看历史变量之间是否重复、是否和环境变量高度相关。
diag <- diagnose_history_indices(hist, env = env)
diag$warnings
diag$missing_summary
diag$high_correlations
diag$environment_overlap

# 后续用户自己正常跑 HMSC，例如：
# model <- Hmsc::Hmsc(
#   YData = comm,
#   XData = XData,
#   XFormula = ~ temp + precip + soil_pH + colonization_age + prob_weighted_source_diversity,
#   TrData = TrData,
#   TrFormula = ~ body_size + dispersal_trait + trait_conservatism,
#   distr = "probit"
# )
