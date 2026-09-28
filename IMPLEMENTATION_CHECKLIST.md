# HmscEcoEvo Implementation Checklist

Status labels:

- Direct: computed by HmscEcoEvo from `hmsc_ecoevo` data.
- Optional package/precomputed: computed when an optional package is installed or a user table is supplied.
- Precomputed required: standardized and plotted by HmscEcoEvo, but must come from an external model/workflow.
- Optional data: computed only when population, time-series, coordinate, tree-set, or validation data are supplied.

## Package Migration

- Package name changed to `HmscEcoEvo`: Direct.
- Historical workflow retained: Direct.
- Old exports retained: Direct.
- `hmscHist_data()` backward-compatible alias: Direct.
- New `HmscEcoEvo_data()` constructor: Direct.
- `hmsc_history_indices` class retained with new `HmscEcoEvo_history_indices` class: Direct.

## Unified Object and Extractors

- `hmsc_ecoevo()`: Direct.
- `as_hmsc_ecoevo()`: Direct.
- `extract_hmsc_ecoevo()`: Direct.
- Species/trait/phylogeny/history alignment checks: Direct.
- Beta orientation standardized to species x axis: Direct.
- Gamma orientation standardized to trait x axis: Direct.
- Posterior draw summaries with mean, median, lwr, upr, support: Direct.

## 1. Niche Summary

- `calc_niche_metrics()`: Direct.
- `plot_niche_summary()`: Direct.
- Beta coefficient: Direct.
- Beta posterior support: Direct when `beta_draws` supplied.
- Beta uncertainty: Direct when `beta_draws` supplied.
- Niche response curves: Direct.
- Niche optimum: Direct for valid quadratic terms.
- Niche breadth: Direct for valid quadratic terms.
- Beta-vector magnitude: Direct.
- Beta-vector direction: Direct.
- Beta-vector angle / cosine similarity: Direct.
- Beta-PCA: Direct.

## 2. Phylogenetic Signal

- `calc_phylo_signal_metrics()`: Direct/optional.
- `plot_phylo_signal()`: Direct.
- HMSC rho: Direct when supplied/extractable.
- Pagel's lambda: Optional package/precomputed (`phytools` or `precomputed$pagel_lambda`).
- Blomberg's K: Optional package/precomputed (`phytools` or `precomputed$blomberg_k`).
- Abouheif's Cmean: Precomputed required (`precomputed$abouheif_cmean` or external `phylosignal` workflow).
- Moran's I: Direct with `ape` phylogeny.
- Phylogenetic correlogram: Direct with `ape` phylogeny.
- Scale-dependent signal: Direct descriptive bins from correlogram.
- Beta-LIPA / local Moran's I: Direct descriptive local Moran.
- Node-level signal: Direct descriptive clade-node variance/signal with `ape`.
- Clade Conservatism Index: Direct when clade labels supplied.
- Residual phylogenetic map: Direct residual Beta after traits/history.

## 3. Evolutionary Transition

- `calc_evo_transition_metrics()`: Direct + precomputed interfaces.
- `plot_evo_transition()`: Direct.
- Branch shift probability: Precomputed required (`precomputed$branch_shifts`).
- Branch contrast outlier score: Direct descriptive.
- Number of shifts: Direct descriptive threshold; model-based counts if precomputed.
- Delta theta shift magnitude: Direct descriptive branch delta; model-based if precomputed.
- Axis-specific shift: Direct descriptive branch delta.
- Beta-rate sigma2 / axis-specific rate: Direct descriptive PIC/variance.
- Rate ratio: Direct descriptive by clade when clades supplied.
- P_rate_shift: Precomputed required (`precomputed$rate_shifts`).
- Axis-specific rate: Direct descriptive.
- Clade niche variance: Direct when clades supplied.
- Variance ratio: Direct when clades supplied.
- Niche expansion/contraction: Direct from variance ratio.
- Beta-disparity: Direct descriptive variance.
- DTT: Precomputed required (`precomputed$dtt`).
- MDI: Precomputed required (`precomputed$mdi` or DTT workflow).
- Early-burst parameter: Precomputed required (`precomputed$early_burst`).
- Convergence score: Direct descriptive far-phylo/low-Beta-distance screen or precomputed.
- Peak reuse index: Precomputed required (`precomputed$peaks`).
- Phylogenetic distance vs Beta distance: Direct with phylogeny.
- Beta cosine similarity: Direct via niche metrics.
- Evolutionary integration network: Direct Beta-axis correlation network.
- Evolutionary modularity network: Direct descriptive clustering of Beta-axis correlations.

## 4. Trait / History Mediation

- `calc_trait_mediation_metrics()`: Direct.
- `plot_trait_mediation()`: Direct.
- Gamma heatmap: Direct when Gamma supplied.
- Trait-explained R2: Direct with traits.
- History-explained R2: Direct with species-level history.
- Residual rho: Direct descriptive residual Moran's I with phylogeny.
- TMNS: Direct descriptive trait-mediated niche signal.
- HMNS: Direct descriptive history-mediated niche signal.
- THMNS: Direct descriptive trait+history-mediated niche signal.
- HPNS: Direct descriptive hidden/residual phylogenetic niche signal.
- Missing trait risk index: Direct descriptive.
- Trait phylogenetic redundancy: Direct with traits; Moran component needs phylogeny.
- Trait omission sensitivity: Direct leave-one-trait-out R2 or precomputed.
- Residual Beta tree heatmap: Direct.

## 5. Gamma Evolution

- `calc_gamma_evolution_metrics()`: Direct + precomputed interfaces.
- `plot_gamma_evolution()`: Direct.
- Gamma matrix: Direct when supplied.
- Clade-specific Gamma: Direct descriptive within-clade Beta ~ traits or precomputed.
- Regime-specific Gamma: Direct descriptive within-regime Beta ~ traits or precomputed.
- Gamma-shift probability: Precomputed required (`precomputed$gamma_shift_probability`).
- Gamma sign-flip index: Direct when clade Gamma exists.
- Gamma variance among clades: Direct when clade Gamma exists.
- Trait-function turnover index: Direct when clade Gamma exists.
- Pairwise turnover matrix: Direct when clade Gamma exists.
- Trait-axis specialization: Direct when Gamma exists.
- Trait-function network: Direct when Gamma exists.

## 6. Population Evolution

- `calc_population_evolution_metrics()`: Optional data.
- `plot_population_evolution()`: Optional data.
- Population Beta: Direct when population data supplied.
- Within-species niche divergence: Direct when population Beta supplied.
- Local adaptation contribution: Optional data/precomputed.
- Plasticity contribution: Optional data/precomputed.
- LA:PL ratio: Optional data/precomputed.
- Genetic niche signal: Optional data/precomputed.
- Reaction norms: Optional data/precomputed.
- G x E: Optional data/precomputed.
- Eco-evolutionary feedback path: Precomputed required for causal interpretation; descriptive time-series summaries supported.
- Beta_t shift rate: Optional time-series data/precomputed.
- Trait evolution contribution: Optional data/precomputed.
- Community-driven selection index: Optional data/precomputed.
- Feedback strength: Optional descriptive correlation/precomputed; not causal by itself.
- Eco-evo turnover ratio: Optional time-series data/precomputed.
- Lagged niche shift: Optional data/precomputed.
- Interaction-mediated evolution: Optional data/precomputed.

## 7. Validation Dashboard

- `calc_validation_metrics()`: Direct + optional/precomputed.
- `plot_validation_dashboard()`: Direct.
- Posterior uncertainty propagation: Direct with posterior draws.
- Tree uncertainty sensitivity: Precomputed required.
- Trait omission sensitivity: Direct via trait mediation or precomputed.
- Spatial confounding check: Direct with `XData` and coordinates.
- Environment omission check: Precomputed required.
- Prior sensitivity: Precomputed required.
- Simulation recovery: Precomputed required.
- Block cross-validation: Precomputed required.
- Posterior predictive check: Precomputed required unless user supplies HMSC workflow output.
- ESS/Rhat/trace/MCMC diagnostics: ESS direct from posterior draws with `coda`; Rhat requires multi-chain precomputed diagnostics.
