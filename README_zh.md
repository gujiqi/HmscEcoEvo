# HmscEcoEvo

HmscEcoEvo 是在原 `hmscHist` 历史生物地理 predictor 工作流基础上扩展出的 HMSC 生态-进化诊断包。

它的定位不是替代 `Hmsc`、BioGeoBEARS、bayou、SURFACE、`phytools`、`geiger` 或 `phylosignal`。它做三件事：

1. 保留原 `hmscHist` 的历史表标准化、历史指数计算和 HMSC 输入整理功能。
2. 从 HMSC 或用户提供的结果中统一提取 Beta、Gamma、Rho、Omega、posterior draws、traits、phylogeny、historical predictors、population/time-series 数据。
3. 构建统一 S3 对象 `hmsc_ecoevo`，计算并绘制七类生态-进化诊断。

## 当前六过程框架

`hee_core_process_catalog()` 的核心过程只有 Environmental filtering、Dispersal、Biotic filtering、Evolution、Speciation 和 Extinction。地学历史是外部驱动；化石、现代末端、BioGeoBEARS/BSM 是可选约束；避难所、通道和多样性是结果或诊断，不是新增过程。

旧版 `hee_colonisation_*()` 和 `hee_persistence_probability()` 仍为可复现历史情景保留兼容接口，但不属于当前六过程主模型，不能将其未经独立校准的值写成 Case05 的实测定殖率或存续率。当前 Case05 用定年树、逐期古陆地、板块携带和主动移动推断谱系位置；环境支持与位置分布分开报告，位置分布不是局地种群占据概率。

## 快速示例

```r
library(HmscEcoEvo)

evo <- hmsc_ecoevo(
  beta = beta_matrix,
  beta_draws = beta_draws,
  gamma = gamma_matrix,
  rho = rho_summary,
  traits = traits,
  phylo = phy,
  history = hist_indices
)

niche <- calc_niche_metrics(evo)
plot_niche_summary(niche)
```

## 七大模块

- `calc_niche_metrics()` / `plot_niche_summary()`：Beta 热图、posterior support、不确定性、响应曲线、最适值、生态位宽度、Beta 向量强度/方向、cosine similarity、Beta-PCA。
- `calc_phylo_signal_metrics()` / `plot_phylo_signal()`：HMSC rho、Pagel lambda、Blomberg K、Abouheif Cmean、Moran I、correlogram、尺度依赖信号、LIPA/local Moran、node-level signal、CCI、residual Beta map。
- `calc_evo_transition_metrics()` / `plot_evo_transition()`：branch shift、shift magnitude、rate、variance ratio、生态位扩张/收缩、DTT、MDI、EB、convergence、peak reuse、phylogenetic distance vs Beta distance、Beta cosine、integration/modularity。
- `calc_trait_mediation_metrics()` / `plot_trait_mediation()`：Gamma 热图、trait/history R2、residual rho、TMNS/HMNS/THMNS/HPNS、missing trait risk、trait phylogenetic redundancy、trait omission sensitivity、residual Beta heatmap。
- `calc_gamma_evolution_metrics()` / `plot_gamma_evolution()`：clade/regime-specific Gamma、Gamma-shift probability、sign flip、Gamma clade variance、turnover matrix、specialization、trait-function network。
- `calc_population_evolution_metrics()` / `plot_population_evolution()`：population Beta、种内生态位分化、local adaptation、plasticity、LA:PL ratio、genetic niche signal、reaction norms、GxE、eco-evolutionary feedback、trait evolution contribution、eco-evo turnover。
- `calc_validation_metrics()` / `plot_validation_dashboard()`：posterior uncertainty propagation、tree uncertainty、trait omission、spatial/environment confounding、prior sensitivity、simulation recovery、block CV、posterior predictive check、ESS/Rhat/trace/MCMC。

## 旧功能兼容

原来的函数仍然导出：

```r
proj <- hmscHist_data(comm, site_region, env = env, traits = traits)
proj <- build_history_tables(proj, species_region_history = species_region_history)
hist <- calc_history_indices(proj)
XData <- as_hmsc_xdata(hist, env = env)
TrData <- as_hmsc_trdata(hist, traits = traits)
```

新推荐写法是：

```r
proj <- HmscEcoEvo_data(comm, site_region, env = env, traits = traits)
```

## 重要边界

HmscEcoEvo 不会从 HMSC 的 Beta 里伪造 bayou/SURFACE/OU/BM/EB 结果。branch shift probability、rate shift probability、adaptive peak、convergence regime、tree uncertainty、prior sensitivity、simulation recovery 和严格 posterior predictive check 等指标，需要用户传入外部预计算表。

historical predictors 和 residual association 只能作为诊断线索，不能直接解释为因果证据。
