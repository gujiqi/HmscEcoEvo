# HmscEcoEvo example-case and diagnostic checklist

This checklist was generated for the seven small runnable example cases in
`inst/extdata/ecoevo_cases/`. The cases are toy data for exercising code paths,
not empirical evidence.

## Example data

- `case01_niche_summary`
- `case02_phylo_signal`
- `case03_evo_transition`
- `case04_trait_history_mediation`
- `case05_gamma_evolution`
- `case06_population_evolution`
- `case07_validation_dashboard`

Each case contains `evo.rds`, `beta.csv`, `gamma.csv`, `beta_draws.rds`,
`gamma_draws.rds`, `traits.csv`, `history_species.csv`, `phylo.tre`,
`omega.csv`, `population.csv`, `timeseries.csv`, `XData.csv`, `coords.csv`,
`Y.csv`, `species_groups.csv`, and `precomputed.rds`.

The `hmsc_ecoevo()` constructor also extracts, when available from an HMSC
object or HMSC-like `postList`, Beta/Gamma means, Beta/Gamma posterior draw
arrays, rho posterior draws, Omega posterior draw arrays reconstructed from
Lambda, species/covariate/trait names, TrData, phylogeny, XData, and Y. If an
`hmsc_history_indices` object is supplied, `TrData_history` is aligned to
species and `XData_history` is merged into `XData`.

## Scripts

- `generate_ecoevo_cases.R`: regenerates the seven case datasets.
- `run_all_case_plots.R`: runs all seven calc functions and all seven plot
  functions for every case, saves metrics and plot outputs.

## Diagnostic coverage

### 1. Niche summary

Implemented directly from Beta/posterior draws:
`beta heatmap`, `beta posterior support/uncertainty`, `niche response curves`,
`niche optimum`, `niche breadth`, `beta-vector magnitude/direction`,
`beta-PCA`, `beta cosine similarity`.

### 2. Phylogenetic signal

Implemented directly when phylogeny is supplied:
`Hmsc rho`, `Moran's I`, binned Moran's-I `phylogenetic correlogram`,
`scale-dependent signal`, `beta-LIPA/local Moran's I`,
`node-level signal`, `clade conservatism index`,
`residual phylogenetic map`.

Precomputed or optional-package backed:
`Pagel's lambda`, `Blomberg's K`, `Abouheif's Cmean`, and full
`phylosignal::phyloCorrelogram` with confidence bands when `phylosignal` and
`phylobase` are installed.

### 3. Evolutionary transition

Implemented descriptively or from precomputed branch/rate tables:
`branch shift probability`, `number of shifts`, `delta theta shift magnitude`,
`axis-specific shift`, `rate ratio`, `P_rate_shift`, `axis-specific rate`,
`variance ratio`, `niche expansion/contraction`, `DTT`,
`early-burst parameter`, `MDI`, `convergence score`, `peak reuse index`,
`phylogenetic distance vs beta distance`, `beta cosine similarity`,
`evolutionary integration/modularity network`, `multivariate branch shift`,
`trait-niche covariance`.

Requires external precomputed model output for model-based inference:
BM/OU/EB model fits, bayou/SURFACE-style branch shifts, adaptive peaks,
rate-shift models, convergence-regime models.

### 4. Trait/history mediation

Implemented:
`Gamma heatmap`, `trait-explained R2`, `history-explained R2`,
`combined trait-history R2`, `residual rho`, `TMNS`, `HMNS`, `THMNS`, `HPNS`,
`missing trait risk index`, `trait phylogenetic redundancy`,
`trait omission sensitivity`, `residual beta tree heatmap`.

Interpretation guard:
historical predictors and residual association are descriptive diagnostics, not
causal evidence.

### 5. Gamma evolution

Implemented:
`clade-specific Gamma`, `regime-specific Gamma`, `Gamma-shift probability`,
`Gamma sign-flip index`, `Gamma variance among clades`,
`trait-function turnover index`, `pairwise turnover matrix`,
`trait-axis specialization`, `trait-function network`.

Requires external precomputed model output for model-based Gamma shifts.

### 6. Population evolution

Implemented when population/time-series data are supplied:
`population beta`, `within-species niche divergence`,
`local adaptation contribution`, `plasticity contribution`, `LA:PL ratio`,
`genetic niche signal`, `reaction norms`, `GxE`,
`eco-evolutionary feedback path`, `trait evolution contribution`,
`eco-evo turnover ratio`, `beta(t) shift rate`,
`community-driven selection index`, `feedback strength`,
`lagged niche shift`, `interaction-mediated evolution`.

Interpretation guard:
feedback paths are descriptive unless supplied from a causal/mechanistic model.

### 7. Validation dashboard

Implemented or standardized:
`posterior uncertainty propagation`, `tree uncertainty sensitivity`,
`trait omission sensitivity`, `spatial confounding check`,
`environment omission check`, `prior sensitivity`, `simulation recovery`,
`block cross-validation`, `posterior predictive check`,
`ESS/Rhat/trace/MCMC diagnostics`.

Requires external precomputed workflows for tree-set sensitivity, prior
sensitivity, environment omission model comparison, simulation recovery, block
CV, and full PPC when those are not already supplied.

## Fresh verification

The example run produced 49 diagnostic PDFs and 406 PNG panels under
`inst/extdata/ecoevo_cases/outputs/`. The outputs directory is ignored during
`R CMD build` to avoid bloating the installed package; regenerate it with
`run_all_case_plots.R`.
