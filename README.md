# HmscEcoEvo

HmscEcoEvo 1.0.4 implements HmscEE with six biological processes while
preserving the historical-predictor workflow from `hmscHist`.

The package is not a replacement for `Hmsc`, BioGeoBEARS, bayou, SURFACE,
`phytools`, `geiger`, `phylosignal`, plate-reconstruction software, or causal
eco-evolutionary modelling packages. Its core role is to keep scale and process
boundaries explicit: HMSC and phylogenetic response reconstruction supply
ancestral environmental responses; palaeogeography and plate histories supply
hard habitat existence and movement landscapes; palaeoclimate supplies
time-matched environmental predictors. Current global-grid Case05 infers
lineage-location distributions with time-ordered movement and extant endpoint
conditioning. Optional BioGeoBEARS/BSM analyses constrain regional histories;
older scenario-calibrated occupancy workflows remain reproducible but are not
the default empirical Case05. Extrapolation diagnostics record transferability
risk; they are not biological absence probabilities. The package deliberately
does not use `e` for Earth scenario, because BioGeoBEARS already uses `e` for
range contraction/local extinction.

In one sentence:

```text
HmscEE = ancestral environmental response
       + dynamic Earth-constrained lineage movement
       + dated-tree lineage history
       | modern endpoints and optional historical constraints
```

This replaces the older product-form scenario
`S_HMSC * A_hist * E_phylo * L_land * D_dynamic * Q_extrap` as the core model.
That product is retained only as an obsolete compatibility guard; it stops
with an explanation rather than producing a new core HmscEE result.

HmscEE uses six biological process names in the current public catalog:

| Process | Preferred prefix | Question | Main implementation |
|---|---|---|---|
| Environmental filtering | `hee_environmental_filtering_*()` | How does a lineage respond to the environment at this time? | HMSC plus ancestral environmental response produces environmental support `S_env`; this is not historical occupancy. |
| Dispersal | `hee_dispersal_*()` | How does a lineage move across a supplied land-only cell graph? | Spherical or resistance-aware movement is kept distinct from plate carriage; time-ordered pruning conditions locations on the dated tree and modern endpoint. |
| Biotic filtering | `hee_biotic_filtering_*()` | What independently supported effects do other organisms exert? | Optional external interaction term; not estimated from HMSC residual association alone. |
| Evolution | `hee_evolution_*()` | How do traits and environmental responses change through the tree? | Ancestral traits/responses, including `beta = Gamma T + u` when traits support it. |
| Speciation | `hee_speciation_*()` | When do lineages split? | The dated tree supplies lineage identity and node times; optional regional histories can add geography. |
| Extinction | `hee_extinction_*()` | What disappears at geographic, local-population, or lineage scale? | Geographic land loss is not automatically biological extinction; complete lineage extinction needs extinct-tip/fossil evidence. |

Colonisation and persistence are **not** seventh and eighth processes. Older
`hee_colonisation_*()` and `hee_persistence_probability()` functions are kept
only as legacy compatibility for explicitly calibrated or scenario-based
dynamic-occupancy examples. They are not called by the current formal Case05,
and their outputs must not be labelled empirical historical rates.

See [CASE_WORKFLOWS.md](CASE_WORKFLOWS.md) for the current versus legacy status
of Case01-Case07 entry points.

### Landscape-explicit dispersal

Version 1.0.1 makes the preferred dispersal API explicit and sparse. It
adopts a landscape-weighted spherical diffusion kernel for active movement and
a separate time-ordered graph for plate transport. For source cell `i` and
target cell `j`, the directional component is:

```text
K[i -> j, lineage, time] proportional to
  cell_area[j] * movement_landscape_weight[j] *
  exp(-effective_movement_distance_rad[i,j]^2 /
      (2 * diffusion_variance_rad2_per_myr * delta_t_myr))
```

Weights are normalised only across valid destinations of one source cell. A
separate `source_emigration_rate_per_myr` turns `K` into a movement-rate matrix;
arrival then feeds colonisation, and establishment remains conditional on
arrival. This is not a product of independent distance, connectivity, barrier,
accessibility and suitability scores. `effective_movement_distance` is the
great-circle distance unless the sparse candidate edges supply a documented
resistance-expanded `effective_cost_km` from the continuous topographic graph.
Use either that edge cost or the same topographic quantity as a cell-level
movement weight, not both.

```r
# Dynamic Earth: hard land state + one documented movement permeability.
movement_landscape <- hee_dispersal_landscape_weights(
  paleo_cells,
  hard_mask_col = "land_exists",
  movement_weight_col = "movement_habitat_weight",
  cell_area_col = "cell_area_km2"
)

# Sparse 4/8-neighbour candidates can come from the existing topographic graph.
resistance <- hee_dispersal_resistance_surface(paleo_cells)
local_graph <- hee_dispersal_connectivity_graph(resistance)

# Active biological movement on one historical time slice.
K <- hee_dispersal_spherical_kernel(
  movement_landscape,
  edges = local_graph$edges,
  time_ma = 40, delta_t_myr = 0.25,
  diffusion_variance_rad2_per_myr = 0.002,
  source_emigration_rate_per_myr = 0.1,
  spatial_domain = "global_explicit"
)
movement_rate <- hee_dispersal_transition_matrix(K, cell_order = cell_ids)
```

Use `spatial_domain = "within_region"` only when a time-stratified
BioGeoBEARS/BSM reconstruction already owns cross-region range evolution. In
that mode, active HmscEE movement is zero across regions. In a no-BioGeoBEARS
model, use `"global_explicit"` and make the global palaeogeographic graph,
root prior and endpoint/fossil conditioning explicit. `hee_plate_grid_transport()`
and `hee_dispersal_spacetime_graph()` represent movement of the Earth reference
frame before active dispersal; plate carriage is not organismal dispersal.

`hee_dispersal_path_diagnostic()` is for one interpretable least-cost route or
phylogeographic diagnostic, not for building an all-pairs occupancy matrix.
`hee_dispersal_pruning_likelihood()` implements a finite-state, explicit-grid
Felsenstein constraint for narrow-range phylogeographic data; it does not by
itself model wide ranges, full extinction, or community assembly. The former
`hee_dispersal_kernel()`, `hee_dispersal_edge_kernel()`, and
`hee_dispersal_particle_kernel()` remain available for reproducibility, but
new analyses should use `hee_dispersal_spherical_kernel()` and, when finite
propagules are needed, `hee_dispersal_propagule_particles()`.

Everything else is a layer around those engines: climate/geology/geography and
habitat are external drivers; BioGeoBEARS, dated trees and plate models are
historical constraints or inference tools; fossils, pollen, aDNA,
cross-validation, no-analog climates and phyloregion analyses are observation,
uncertainty or result layers. Use `hee_core_process_catalog()` and
`hee_function_catalog()` to inspect this map in machine-readable tables.

The active-kernel parameterisation follows the explicit-grid spherical
diffusion formulation of Arias (2024, [doi:10.1093/sysbio/syae051](https://doi.org/10.1093/sysbio/syae051)).
The separate time-ordered landscape graph is informed by the
landscape-explicit phylogeographic workflow of Flannery-Sutherland et al.
(2025, [doi:10.1038/s41559-025-02739-y](https://doi.org/10.1038/s41559-025-02739-y)).
Hmsc supplies the hierarchical modern community/environmental-response layer;
it does not, by itself, estimate a deep-time dispersal kernel from one modern
community snapshot.

Full workflow scripts use optional helper packages for plotting, report
generation, local source loading, and raster-ready objects (`gridExtra`,
`scales`, `pkgload`, `officer`, and `terra`). These are declared in `Suggests`;
missing helper packages should be reported as dependency limitations, not as
scientific evidence that a workflow result is valid.

## Core Workflow

```r
library(HmscEcoEvo)

evo <- hmsc_ecoevo(
  beta = beta_matrix,            # species x environmental/niche axis
  beta_draws = beta_draws,       # optional draw x species x axis
  gamma = gamma_matrix,          # optional trait x niche axis
  rho = rho_summary,             # optional HMSC rho
  omega = omega_matrix,          # optional residual associations
  traits = traits,               # species x trait
  phylo = phy,                   # ape::phylo
  history = hist_indices,        # hmsc_history_indices or species history table
  population = population_beta,  # optional
  timeseries = ecoevo_ts         # optional
)

niche <- calc_niche_metrics(evo)
plot_niche_summary(niche)

phylo_signal <- calc_phylo_signal_metrics(evo)
plot_phylo_signal(phylo_signal)
```

## Vignette Navigation

The vignettes have been consolidated into a numbered manual. They are written
as a combined user guide, teaching text, and reviewer checklist. Lightweight
chunks can run during package checks; full MCMC, BioGeoBEARS fitting, plate
reconstruction, and full 540-0 Ma workflows are shown with `eval = FALSE` and
should be run through the scripts.

| Vignette | Best for | Coverage | Runtime | Quick mode | External data needed | R CMD check friendly | Main output |
|---|---|---|---|---|---|---|---|
| `v00_index_roadmap` | Everyone | Reading order, interpretation classes, function map | Light | Yes | No | Yes | Roadmap |
| `v01_basic_hmsc_ecoevo_workflow` | HMSC users and reviewers | `hmsc_ecoevo` object, `calc_*`/`plot_*`, quick example | Light | Yes | No | Yes | First working object and metrics |
| `v02_seven_diagnostic_figures` | Figure readers | Figures 1-7, panel logic, formulas, misreads | Light | Yes | Optional phylogeny/traits/precomputed tables | Yes | Diagnostic plot guide |
| `v03_history_to_hmsc` | hmscHist users | Historical-process tables, trait-history indices, HMSC inputs, leakage guard and bundled trait-history mediation figures | Light | Yes | User history tables for real use | Yes | HMSC-ready historical predictors plus result-reading gallery |
| `v04_deep_time_hindcasting_540Ma` | Palaeoenvironment users | 0 Ma recipe, full time axis, masks, M1-M5 | Mostly documented | Yes; full via script | env_cube, landmask, origin times | Yes | Deep-time workflow recipe |
| `v05_geological_processes_dynamic_assembly` | Geological-process users | Legacy heuristic connectivity, colonisation, extinction, assembly and reliability scenarios; not the current six-process Case05 | Light | Yes | Palaeogeographic layers for real use | Yes | Process formula and scenario guide |
| `v06_biogeography_plate_correction` | Biogeography users | BioGeoBEARS, accessibility, BSM, plate tracks | Light | Yes | BGB/plate outputs for real use | Yes | External input schemas |
| `v07_case01_full_simulated_results` | Reproducing and reviewing Case 1 | Simulated inputs, HMSC/MCMC, seven diagnostics, 540-0 Ma full axis, M1-M5 formulas, output tables/figures, Word-report checks | Reads saved output plus runnable mini summary | Script supports quick/full | No real external data | Yes | Complete simulated result-reading manual |
| `v08_case02_realdata_template` | Real projects | `data_raw/` schemas, units/ranges, stop rules, real BGB/plate/fossil inputs, HMSC/CV/M1-M5 code, missing-data decisions | Template with runnable schema check | Script supports quick/full | Yes | Yes | Real-data template and limitation checklist |
| `v09_advanced_diagnostics_and_reporting` | Comparative-method and manuscript users | Posterior, niche, phylogenetic-signal, transition, geoprocess, validation formulas and paper-safe wording | Light formula guide | Yes | Precomputed external results for real use | Yes | Formula catalog and reporting templates |
| `v10_case03_geoprocess_past_to_present` | Geological-process workflow users | Separate real-data and simulated-demo entry points; real 0 Ma HMSC fitting, full palaeo grids, posterior uncertainty, time-specific regions, km connectivity, M1-M5 legality, geological diagnostics and reports | Quick smoke test plus full command | Yes; full via script | Real mode requires project community, sites and environment; M2-M5 need dated tree/origins and M3-M5 need external BioGeoBEARS accessibility | Yes | Real workflow plus explicitly labelled teaching demo |
| `v11_case04_plant200_ancestral_suitability` | Plant200 ancestor-response users | HMSC Beta posterior, dated tree, active lineages, ancestor-Beta trajectories and palaeoenvironment support | Reads the available development output; full run documented | Yes | Plant200 HMSC posterior, dated tree and historical environment cube | Yes | Ancestor-response result reader; not a final occupancy reconstruction |
| `v12_case05_global_dispersal_diversity` | Plant200 result readers | Formal 4-degree plate-carriage, spherical-diffusion and terminal-calibrated lineage-support maps | Reads formal PNG/GeoTIFF inventory; full run documented | Yes | Formal Case05 result directory | Yes | Complete map-category reader with strict interpretation boundaries |
| `v13_case05_patch_history_results` | Case05 patch-history readers | Connected supported patches, plate-carriage links and candidate refugia | Reads saved derived outputs | Yes | Case05 v8 patch-history results | Yes | Separates potential corridors from unobserved movement flux |
| `v14_case06_bm_ou_eb_plant200` | Ancestral-response model comparison | Matched BM, OU and EB fits, ancestor-Beta trajectories and time-matched palaeoenvironmental support | Reads eight real-data figures; runner documented | Yes | Plant200 HMSC posterior, dated tree and PALEOMAP carrier | Yes | Boundary-aware model comparison, not a historical occupancy reconstruction |

## Seven Diagnostic Modules

1. `calc_niche_metrics()` / `plot_niche_summary()`: Beta heatmap, posterior support and uncertainty, response curves, niche optimum, niche breadth, Beta-vector magnitude/direction, cosine similarity, and Beta-PCA.
2. `calc_phylo_signal_metrics()` / `plot_phylo_signal()`: HMSC rho, Pagel's lambda, Blomberg's K, Abouheif's Cmean, Moran's I, correlogram, scale-dependent signal, local Moran/LIPA, node-level signal, CCI, and residual Beta map.
3. `calc_evo_transition_metrics()` / `plot_evo_transition()`: branch shifts, shift magnitude, rates, variance shifts, expansion/contraction, DTT/MDI/EB placeholders or precomputed results, convergence, peak reuse, phylogenetic distance vs Beta distance, cosine similarity, and integration/modularity.
4. `calc_trait_mediation_metrics()` / `plot_trait_mediation()`: Gamma heatmap, trait/history R2, residual rho, TMNS/HMNS/THMNS/HPNS, missing-trait risk, trait phylogenetic redundancy, trait omission sensitivity, and residual Beta heatmap.
5. `calc_gamma_evolution_metrics()` / `plot_gamma_evolution()`: clade/regime-specific Gamma, Gamma-shift probability, sign flips, among-clade variance, turnover matrix, specialization, and trait-function network.
6. `calc_population_evolution_metrics()` / `plot_population_evolution()`: population Beta, within-species divergence, local adaptation, plasticity, LA:PL ratio, genetic niche signal, reaction norms, GxE, feedback paths, trait-evolution contribution, and eco-evo turnover.
7. `calc_validation_metrics()` / `plot_validation_dashboard()`: posterior uncertainty propagation, tree uncertainty sensitivity, trait omission, spatial/environment confounding checks, prior sensitivity, simulation recovery, block CV, posterior predictive checks, and MCMC diagnostics.

## Historical Predictor Workflow

All original `hmscHist` functions remain exported, including:

```r
proj <- hmscHist_data(comm, site_region, env = env, traits = traits)
proj <- build_history_tables(proj, species_region_history = species_region_history)
hist <- calc_history_indices(proj)
XData <- as_hmsc_xdata(hist, env = env)
TrData <- as_hmsc_trdata(hist, traits = traits)
```

The preferred new constructor name is:

```r
proj <- HmscEcoEvo_data(comm, site_region, env = env, traits = traits)
```

`hmscHist_data()` is kept as a backward-compatible alias.

## Deep-Time Hindcasting Workflow

HmscEcoEvo also includes a deep-time HMSC hindcasting layer for Phanerozoic
paleoenvironmental time cubes. The current HmscEE core is not a simple
time-slice multiplier. It is a conditional decomposition:

```text
P(cell occupied) =
  P(region occupied by BioGeoBEARS/BSM) *
  P(cell occupied inside that region | region occupied)
```

The main update is performed by `hee_nested_region_cell_occupancy()`:

```r
E <- hee_paleo_earth_state(earth, env_cols = c("MAT", "MAP"))
R <- hee_bsm_region_history(bsm_region_history)
ancestral_response_draws <- hee_evolution_ancestral_response_direct(
  tip_beta_draws = hmsc_beta_draws,
  tree = dated_tree,
  times = sort(unique(E$time_ma), decreasing = TRUE),
  basis_cols = c("MAT", "MAP")
)
S <- hee_lineage_suitability(E, ancestral_response_draws,
                             basis_cols = c("MAT", "MAP"))
K <- hee_within_region_movement_kernel(least_cost_paths,
                                       traits = dispersal_traits)

fit <- hee_nested_region_cell_occupancy(
  suitability = S,
  region_history = R,
  earth_state = E,
  movement_kernel = K
)
```

Do **not** reconstruct ancestors by averaging modern prediction maps. The
scientific path is modern HMSC tip `Beta` posterior -> dated-tree ancestral
environmental-response `Beta` -> palaeoenvironmental support `S_env`. The
shortcut wrapper is:

```r
S_env <- hee_environmental_filtering_ancestral_suitability(
  earth_state = E,
  tip_beta_draws = hmsc_beta_draws,
  tree = dated_tree,
  basis_cols = c("MAT", "MAP"),
  link = "probit"
)
```

`S_env$summary$suitability_mean` is ancestral environmental support, not final
historical occupancy. BioGeoBEARS region history and within-region dynamic
occupancy are applied after this step.

The core equations are:

```text
beta[l,k] = Gamma T[l,k] + u[l,k]
K[l,u'->u,k] = 1 - exp(-lambda_mov[l,u'->u,k] * delta_t[k])
A_arr[l,u,k] = 1 - prod(1 - q[l,u',k] K[l,u'->u,k])
gamma[l,u,k] = A_arr[l,u,k] * logit^-1(a_l + b_S eta[l,u,k])
phi[l,u,k] = logit^-1(a_l^P + b_P eta[l,u,k])
q[l,u,k+1] = q[l,u,k] phi[l,u,k] + (1 - q[l,u,k]) gamma[l,u,k]
p[l,u,k] = sum_g w[g] sum_h w[h|g] sum_c w[c|g] sum_s w[s]
           R[l,r(u,k),k]^(h,g) q[l,u,k]^(h,g,c,s)
```

`H_state = 0` or `R_region = 0` is a hard zero. Cross-region movement is not
estimated by HmscEE during ordinary time steps; BioGeoBEARS/BSM supplies the
range expansion, contraction, vicariance, subset and founder-event history.
No-analog and extrapolation scores are attached with
`hee_extrapolation_uncertainty_report()` and do not downweight `p`.
Palaeogeographic scenarios `g` are common upstream conditions: `g` supplies
the structural constraints used to obtain BioGeoBEARS histories `h | g`, and
the same `g` supplies HmscEE hard habitat, cell-region maps and movement
resistance. Palaeoclimate scenarios are written `c | g` because climate models
often depend on the palaeogeography. Do not average BioGeoBEARS histories,
plate/geography scenarios, and palaeoclimate scenarios as independent Cartesian
products when the BioGeoBEARS model used time-stratified geography. Ancestral
traits `T[l,k]`, ancestral environmental response `beta[l,k]`, and ancestral
geographic range `R[l,r,k]` are separate objects; BioGeoBEARS ranges must not
be used to define the ancestral niche.
In these equations `u` is a spatial grid cell; `c` is reserved for the
palaeoclimate scenario to avoid conflict with the BioGeoBEARS parameter `e`.

### Applying the Framework to the Phanerozoic Environment Factors

The palaeoenvironment adapter translates the 51-variable land-harmonized
Phanerozoic products into process-ready inputs. It prefers the PALEOMAP/DEM
geographic frame for the core land term (`land_mask_dem`, `land_area_km2`,
`elevation_m`) and treats Li/model-mean layers as sensitivity or uncertainty
alternatives unless the user explicitly chooses them.

```r
cube <- hee_load_timecube(
  "Phanerozoic_Environment_540_0Ma_5Myr_4deg_landharmonized_v2"
)

roles <- hee_paleoenv_variable_roles(cube)
var_map <- hee_paleoenv_default_variable_map(cube)

paleo_inputs <- hee_prepare_paleoenv_geoprocess_inputs(
  cube,
  times = c(540, 0),
  hmsc_variables = c("MAT_pohl_C", "MAP_pohl_mm_yr",
                     "elevation_m", "moisture_availability_index_z"),
  habitat_strategy = "land_only"
)

arena <- paleo_inputs$arena              # L_land / G_arena / geography
landscape <- paleo_inputs$landscape_state # elevation, climate, area, moisture
opportunity <- paleo_inputs$ecological_opportunity
```

Default role mapping:

| Environment factor | Default framework use | Main output/component |
|---|---|---|
| `land_mask_dem`, `land_area_km2` | Geographic arena and area change | `L_land`, `G_arena`, `cell_area_km2` |
| `MAT_pohl_C`, `MAP_pohl_mm_yr` | HMSC predictors, climatic connectivity and extrapolation checks | `S_HMSC`, climate mismatch, `Q_extrap` diagnostics |
| `elevation_m`, `slope_deg`, `relief_3x3_m`, `tpi_3x3_m` | Topographic resistance, heterogeneity and opportunity | connectivity barriers, `ecological_opportunity` |
| `distance_to_coast_km` | Coastal access/isolation proxy | structural connectivity context |
| `moisture_availability_index_z`, `wetland_potential_index_z`, `bryophyte_moisture_score_z` | Habitat/opportunity proxies; optional habitat gate only when stated | `habitat_availability` or opportunity proxy |
| `MAT_gridded_model_sd_C`, `MAP_gridded_model_sd_mm_yr`, `delta_*` | Climate-model disagreement and reliability diagnostics | uncertainty/reliability components |
| `koppen_*`, `hydro_class_simple` | Categorical summaries and stratified validation | reporting/context layers |

By default `habitat_strategy = "land_only"` so moisture and wetland indices are
not multiplied into `L_land` if they are already used by HMSC. Set
`habitat_strategy = "bryophyte"`, `"wetland"`, or `"moisture"` only for an
explicit habitat-gated sensitivity scenario.

```r
cube <- hee_load_timecube("Phanerozoic_Environment_540_0Ma_5Myr_4deg_landharmonized_v2")
recipe <- hee_lock_recipe(env_now, formula = ~ MAT_pohl_C + MAP_pohl_mm_yr)
mask <- hee_phylo_time_mask(times = c(100, 50, 0), species_origin = origin_ma)
suit <- hee_project_hmsc_timecube(beta = beta, env_cube = cube,
                                  recipe = recipe, times = c(100, 50, 0),
                                  variables = c("MAT_pohl_C", "MAP_pohl_mm_yr"))
proj <- hee_combine(suit, phylo_mask = mask, accessibility = A_BGB)
rich <- hee_richness(proj)
```

Two fast 20-species deep-time examples can be generated with the scripts under
`inst/extdata/deep_time_cases`: a paleo-map case and a tectonic-track case.
Large full-axis RDS outputs are treated as local analysis products and are not
included in source-package builds. Plate reconstruction and BioGeoBEARS fitting
are not performed internally; pass plate-corrected coordinates and
model-averaged accessibility as external inputs.

## Geological Process Diagnostics

The current HmscEE core separates geological arena change, BioGeoBEARS/BSM
regional history, ancestral environmental response, within-region movement,
arrival, establishment, persistence, and derived richness/range summaries.
Older geoprocess helpers remain for compatibility and sensitivity diagnostics,
but they are not the main HmscEE probability model.

| Process | Main function(s) | Required input | Main output | Directly computed by package | Needs external data | Proxy indicator | Main limitation |
|---|---|---|---|---|---|---|---|
| Six-process framework map | `hee_core_process_catalog()`, `hee_function_catalog()` | None | Six-process map plus driver, history, constraint, uncertainty and result layers | Yes | No | Conceptual contract | Use these tables as the first reference for "who is responsible for what"; older occupancy helpers are labelled legacy. |
| Core nested HmscEE probability | `hee_earth_state()`, `hee_bgb_region_history()`, `hee_environmental_filtering_suitability()`, `hee_dispersal_kernel()`, `hee_projection_nested_occupancy()` | Palaeo-Earth `X/H`, BioGeoBEARS/BSM region history `R`, ancestral response draws, within-region least-cost movement edges | Cell occupancy draws `q`, biological probability summary `p`, transition labels, and formula catalog | Yes | HMSC/ancestral response estimates, BSM region histories, palaeoenvironment and movement-resistance layers | Nested conditional model | BioGeoBEARS handles between-region history; HmscEE only refines occupancy inside occupied regions. No `A_hist * D_dynamic * Q_extrap` product is used. |
| Obsolete product-form core | `hee_dynamic_earth_biota_probability()`, `hee_dynamic_state_transition()`, `hee_regional_species_pool_transition()`, `hee_static_snapshot_approximation()` | Former `S_HMSC`, `A_hist`, `E_phylo`, `L_land`, `D_dynamic`, `Q_extrap` components | Stops with an explanatory error | No | Historical only | Obsolete | Kept as compatibility guards so old scripts fail loudly instead of double-counting dispersal. |
| Geographic stage | `hee_geographic_stage()` | Land/geographic existence plus optional habitat, not-ice, availability masks | `G_arena`, `L_land`, `geographic_existence` | Yes | Palaeogeographic layers | External geography-derived mask | Supplied `NA` is treated as absent geography; omitted optional columns are neutral. |
| Landscape gain/loss | `hee_landscape_events()`, `hee_land_age()` | Cell/region land, habitat, area, elevation, or component through `time_ma` | Event labels and continuous land age | Yes | Paleo land/area/elevation layers | Yes | Detects state changes in supplied maps; it is not a tectonic model. |
| Ecological opportunity | `hee_ecological_opportunity()` | Land age plus paleoenvironmental variables | `ecological_opportunity` in `[0, 1]` | Yes | Paleoenvironmental cube | Yes | A configurable opportunity index, not a speciation or diversification rate. |
| Structural connectivity | `hee_structural_connectivity()` | Region/cell pair distances, barriers, shared boundaries | `structural_connectivity` in `[0, 1]` | Yes | Paleodistance/barrier reconstruction | Yes | Distance and barrier kernels must be calibrated for the study system. |
| Functional connectivity | `hee_functional_connectivity()` | Connectivity edges plus species dispersal traits | `functional_connectivity` in `[0, 1]` | Yes | Trait/dispersal data | Yes | Species expansion is explicit only when requested; no uncontrolled Cartesian joins. If no trait table is supplied, the functional layer is neutral (`1`) so structural distance/barrier effects are not double counted. |
| Climatic connectivity | `hee_climatic_connectivity()` | Connectivity edges plus climatic mismatch | `climatic_connectivity` in `[0, 1]` | Yes | Paleoenvironmental contrasts | Yes | Measures environmental continuity, not realised dispersal. |
| Dynamic network metrics | `hee_network_metrics()` | Region-pair connectivity graph through `time_ma` | Node centrality, component count, edge density, fragmentation, corridor persistence | Yes | Dynamic region/cell connectivity | Yes | Graph proxy only; it does not prove realised dispersal routes. |
| Colonisation pressure/probability | `hee_colonisation_pressure()`, `hee_colonisation_probability()` | Previous occupancy, connectivity, suitability, accessibility | Source pressure and `colonisation_probability` | Yes | Occupancy/accessibility inputs | Heuristic state-transition term | Missing process inputs are neutral; explicit zero source evidence gives zero colonisation. |
| Rescue and extinction | `hee_rescue_effect()`, `hee_extinction_probability()` | Previous occupancy/connectivity, suitability, land loss, isolation | `rescue_effect`, `extinction_probability` | Yes | Connectivity and habitat-state inputs | Heuristic state-transition term | Rescue lowers extinction by construction; empirical rates require calibration. Explicit zero-area rows are treated as absent arena and force local extinction. |
| Layered extinction | `hee_extinction_layers()` | Dynamic occupancy and optional extinction probability | Local, regional, and lineage-level extinction proxies | Yes | Independent extinction events for observed interpretation | Yes | Regional and lineage terms are proxies unless calibrated by fossil or lineage-extinction evidence. |
| Speciation-opportunity diagnostics | `hee_speciation_opportunity()` and related helpers | Regional isolation, ecological opportunity, area, lineage presence | `speciation_opportunity` in `[0, 1]` | Yes | Region histories and lineage information | Yes | Post hoc diagnostic only. Core HmscEE does not create new species; dated nodes and BSM inheritance define lineage splitting. |
| Cradle/museum/grave/source/sink roles | `hee_cradle_museum_grave()` | Species-pool and optional source-sink summaries | Regional role proxy scores and dominant status | Yes | Species-pool/source-sink histories | Yes | Descriptive regional role labels, not causal classification. |
| Dynamic assembly | `hee_dynamic_assembly()` | Previous state, suitability, colonisation, extinction, phylo/geographic masks | Occupancy probability through time | Yes | Optional masks and process tables | Semi-mechanistic recursion | Time is Ma and sorted old-to-young internally; deep-time species-level interpretation is limited. |
| Dynamic model family | `hee_dynamic_model_family()` | None | M0-M8 legacy scenario-diagnostic definitions and required components | Yes | External data for fitting higher models | Scenario/proxy map | Backward-compatible sensitivity guide; not the HmscEE core process list. |
| Mechanism hypotheses | `hee_geoprocess_hypotheses()` | None | H1-H8 geological-driver hypotheses with proxy outputs and limits | Yes | External validation for formal hypothesis tests | Driver hypothesis map | Geological hypotheses condition environmental filtering, dispersal, colonisation, biotic filtering, persistence, evolution, speciation and extinction; they are not extra biological processes. |
| Process model comparison | `hee_geoprocess_model_comparison()` | Named list of model outputs plus optional validation summaries | Mean probability, impossible-geography fraction, richness and validation summaries | Yes | Model outputs and validation data | Diagnostic comparison | Descriptive comparison unless paired with formal validation or likelihood. |
| Legacy static projection combination | `hee_combine_static()`, `hee_combine()` | Suitability plus `geographic_existence`, `D_static`, `weight`, accessibility, connectivity terms | Final probability and component columns | Yes | Suitability, masks, optional BGB/dispersal | Legacy multiplicative scenario model | Useful for older Case 1/Case 2 sensitivity contrasts, not the main nested HmscEE core. |
| Mechanism attribution and limitation | `hee_mechanism_attribution()`, `hee_limitation_map()` | Dynamic or static projection components | Dominant limiting process labels | Yes | Projection components | Diagnostic summary | Describes bottlenecks in the model, not causal proof. |
| Prediction reliability | `hee_prediction_reliability()` | Extrapolation, model agreement, sample support, MCMC quality, external confidence | `prediction_reliability` and component scores | Yes | Diagnostics from workflow and external confidence if available | Epistemic diagnostic | It is not a posterior credible interval or frequentist confidence interval. |
| Plate correction | `hee_import_plate_points()`, `hee_reconstruct_points()`, `hee_extract_paleoenv_track()` | Present and paleo coordinates or external reconstruction backend | Paleo coordinates and corrected-vs-uncorrected contrasts | Partly | Plate model, user table, `rgplates`, or `pygplates` | External reconstruction input | The package does not infer plate motion from modern points alone. |
| BioGeoBEARS accessibility | `hee_bgb_prepare()`, `hee_bgb_run_models()`, `hee_bgb_import()`, `hee_bgb_accessibility()` | Ranges, tree, and fitted/imported BioGeoBEARS outputs | Model weights and accessibility through time | Import/runner interface | BioGeoBEARS or precomputed tables | Model-derived external input | Case 1 may use labelled mock data; real-data templates must not treat mocks as evidence. AICc weights are reported only when `n > k + 1` for every candidate model. |
| Scenario stress tests | `hee_geoprocess_scenario_suite()`, `plot_geoprocess_scenario_suite()` | Built-in deterministic normal, extreme, missing, duplicate-key, single-time, uneven-time, shuffled-time, and land appearance/loss scenarios | Input tables, intermediate result tables, validation table, formula catalog, interpretation table, and plot list | Yes | No real data; optional only for method QA | Simulated software-audit cases | These results validate process logic and edge cases; they are not palaeogeographic estimates. |
| Case 03 real-data workflow | `scripts/case03_geoprocess_realdata_past_to_present.R` | Strict project input audit for observed `comm.csv`/`sites.csv`, supplied environment cube, real region history, dated tree and optional external BSM/plate/fossil tables | Missing-input manifest, metadata, and validated `R` region-history tables when inputs exist | Yes | No simulated fallback. Real use must provide project data and BioGeoBEARS/BSM or equivalent region history | Real-data scaffold | `analysis_mode=real` is the default. Missing required files write `missing_required_inputs.csv` and stop. |
| Case 03 simulated demo | `scripts/case03_geoprocess_simulated_demo.R` | Small nested HmscEE teaching workflow with virtual palaeo-Earth `X/H`, BSM-like `R`, within-region `K`, `q`, `p`, and region-pool aggregation | `outputs/case03_hmscee_nested_demo_*` tables, RDS and one probability map figure | Yes | Uses generated teaching inputs only | Simulation | Demonstrates the nested equations; do not cite demo outputs as empirical palaeodistributions, dispersal routes, speciation rates or causal geological effects. The old script name remains a compatibility wrapper. |
| Case 04 Plant200 ancestral suitability | `scripts/case04_plant200_ancestral_suitability_no_bgb.R` | HMSC tip-response draws, dated tree and complete palaeoenvironment grid | Active-lineage tables and ancestor-response trajectories through time | Yes | Plant200 HMSC posterior, dated tree and historical environment cube | Development result | Reconstructs environment-response functions before they are projected to historical climate. It does not average modern maps backward or claim final palaeo-occupancy. |
| Case 05 Plant200 three-dispersal comparison (legacy) | `scripts/case05_plant200_eight_process_three_dispersal.R` | Earlier three-scenario forward comparison | Scenario tables, maps and differences | Yes | Plant200 input bundle and 1° land/habitat/elevation time cube | Legacy sensitivity workflow | Retained for reproducibility; its colonisation, persistence and occupancy-derived maps are not the current formal Case05 result. |
| Case 05 Plant200 geographic support diversity | `scripts/case05_v7_arias_global_sphere_4deg.R` | HMSC response draws, dated tree, PALEOMAP carrier tracks and a global 4° palaeo-Earth grid | 66 time slices of PNG and GeoTIFF maps for geographic location density, environment predictor, destination weight and terminal-calibrated lineage-support diversity | Yes | Plant200 posterior/tree/environment inputs and plate-carriage tracks | Formal descriptive geographic result | Uses plate carriage plus sparse spherical diffusion and modern terminal conditioning. It deliberately does not claim colonisation, persistence, local extinction, final occupancy, species richness, PD or FD without separately calibrated demographic evidence. |
| Case 06 Plant200 BM/OU/EB sensitivity | `scripts/case06_plant200_response_model_study.R` | The same five HMSC Beta draws, 200-tip dated tree and 66 native palaeo-Earth slices for each evolutionary model | Per-axis AICc/boundary audit, focal ancestral paths, global support curves, matched 4-degree maps and difference maps | Yes | Case05 Plant200 fitted posterior and carrier cache | Empirical-input response-model sensitivity | OU fits at its upper parameter bound for all 40 axis-draw combinations; EB nearly collapses to BM. Environmental support is not occupancy. |

Case 03 now follows the nested framework. The real-data script is deliberately
strict: it audits required project files, refuses simulated fallback, validates
BioGeoBEARS/BSM region histories with `hee_bsm_region_history()`, and records
what is missing. The simulated demo is separate and uses small virtual data only
to show how `X/H`, `R`, `K`, `q`, and `p` connect.

The old M3-M5 product interpretation is not used in Case 03 core outputs.
Between-region history belongs to BioGeoBEARS/BSM. Within-region movement uses
`hee_within_region_movement_kernel()` and is set to zero across BSM regions
unless a user explicitly constructs a node/entry event. Extrapolation is joined
with `hee_extrapolation_uncertainty_report()` and is not a biological
downweight.

If `species_origin.csv` is omitted, origin ages can be derived from the
terminal-parent ages of `tree.tre` only when `--tree_time_unit=Ma` is declared
explicitly. A tree measured in substitutions must not be read as geological
time. The derived tip-parent age is an operational lineage-origin constraint,
not direct fossil evidence for a species' first appearance.

BSM events are diagnostic by default. The optional
`--apply_bsm_accessibility_adjustment=true` switch enables a labelled heuristic
sensitivity analysis and records its coefficients; it is not the default
real-data estimate. Fossil validation requires externally reconstructed
palaeogeographic coordinates (`paleo_lon`, `paleo_lat`) and records temporal and
spatial matching distances.

```r
system2(file.path(R.home("bin"), "Rscript"), c(
  "scripts/case03_geoprocess_realdata_past_to_present.R",
  "--analysis_mode=real",
  "--comm=data_raw/comm.csv", "--sites=data_raw/sites.csv",
  "--traits=data_raw/traits.csv", "--tree=data_raw/tree.tre",
  "--env_rds=data_raw/environment.rds",
  "--region_cube=data_raw/region_cube.csv",
  "--bgb_accessibility=data_raw/bgb_accessibility.csv",
  "--plate_points=data_raw/plate_points.csv",
  "--bsm_events=data_raw/bsm_events.csv",
  "--fossils=data_raw/fossils.csv",
  "--output=outputs/case03_hmscee_realdata_audit"
))
```

```r
system2(file.path(R.home("bin"), "Rscript"), c(
  "scripts/case03_geoprocess_simulated_demo.R",
  "--analysis_mode=demo", "--quick=true",
  "--output=outputs/case03_hmscee_nested_demo_quick"
))
```

```r
events <- hee_landscape_events(landscape, area_col = "area_km2")
land_age <- hee_land_age(landscape)
opp <- hee_ecological_opportunity(land_age, env_cols = c("bio1", "bio12"))
conn <- hee_build_connectivity_cube(distances, traits = dispersal_traits)
conn2 <- hee_connectivity_cube(distances = distances, traits = dispersal_traits)
net <- hee_network_metrics(conn2)
iso <- hee_isolation_history(conn)
src <- hee_colonisation_pressure(previous_state, conn)
gamma <- hee_colonisation_probability(suitability, source_pressure = src)
rescue <- hee_rescue_effect(previous_state, conn)
epsilon <- hee_extinction_probability(state, delta_t = 5)
layers <- hee_extinction_layers(dynamic_state, epsilon)
spec <- hee_speciation_opportunity(region_state)
pool <- hee_update_species_pool(initial_pool, lineage_events, times)
roles <- hee_cradle_museum_grave(pool)
dyn <- hee_dynamic_assembly(previous_state, suitability,
                            colonisation = gamma, extinction = epsilon)
why <- hee_mechanism_attribution(dyn)
lim <- hee_limitation_map(why)
pri <- hee_prediction_reliability(extrapolation, model_agreement)
hee_dynamic_model_family()
hee_geoprocess_hypotheses()
suite <- hee_geoprocess_scenario_suite()
suite$validation_table
plot_geoprocess_scenario_suite(suite)
```

These helpers can compute arena creation/loss, structural-functional-climatic
connectivity, isolation duration, source pressure, colonisation probability,
rescue effects, local extinction probability, speciation opportunity, regional
pool updates, mechanism-limitation maps, and the dynamic update
`p[t+1] = G * E * (p[t] * (1 - epsilon) + (1 - p[t]) * gamma)`.
The scenario suite writes the same logic into compact audit cases and explicitly
labels each output as a real estimate, proxy indicator, heuristic scenario, or
external-data-dependent process.
Extrapolation risk and prediction reliability are reported as epistemic
diagnostics by default; they are not silently multiplied into biological
process probability unless the user explicitly chooses a downweighting
sensitivity analysis. True plate reconstruction, fossil-calibrated extinction,
and speciation-rate estimation still require external data or specialist
models.

### Safety defaults for deep-time state transitions

The geological-process helpers separate omitted components from incomplete
component tables. If accessibility, geographic existence, phylogenetic validity,
static dispersal, or dynamic dispersal is not supplied at all, the corresponding
component is neutral and the simpler model remains interpretable. Once a
component table is supplied, however, unmatched species/cell/region/time rows
are treated as zero evidence for that component rather than full accessibility.
This prevents sparse process tables from silently widening deep-time ranges.
Missing colonisation and extinction tables in
`hee_dynamic_assembly()` are neutral (`gamma = 0`, `epsilon = 0`); derive them
explicitly with `hee_colonisation_probability()` and
`hee_extinction_probability()` when you want those processes in the dynamic
update. Species-time phylogenetic masks can be supplied either as long tables
or as species-by-time matrices; unmatched entries mean lineage absence, not
unknown presence. In `hee_colonisation_probability()`, an omitted
`accessibility` argument is neutral, while a supplied accessibility table with
unmatched or missing rows means not historically accessible. Explicit or
unmatched source pressure of zero means no propagule evidence and gives zero
colonisation probability. Negative or non-finite time intervals now error
instead of being silently changed to zero.

For connectivity, omitted barrier/passability columns are neutral, but supplied
barrier, route, or passability columns with `NA` are conservative. Provide
explicit `1` values for routes known to be open. Dispersal scales must be finite
and non-negative; a scale of zero means no movement across positive distance.

Joins between process tables are checked by explicit keys. Duplicate
species/cell/region/time keys now produce a readable error, because silent row
expansion can change richness, source pressure, accessibility, and reliability
maps. `hee_combine()` only uses named process columns such as
`geographic_existence`, `land_mask_dem`, `D_static`, `structural_connectivity`,
`functional_connectivity`, `climatic_connectivity`, `D_dynamic`, and
`Q_extrap`; unrelated numeric columns such as temperature are not guessed as
weights. Dynamic dispersal uses the user-supplied first-slice default only when
there is no previous time slice; later time slices with no source evidence have
`D_dynamic = 0`.
Species-level accessibility and phylogenetic masks must include a `species`
key when more than one species is projected. This prevents a time-only or
region-only table from being silently broadcast to every species. Use
`allow_species_broadcast = TRUE` only for an explicitly shared neutral scenario.
Species-level deep-time masks also require explicit `species_origin` by
default. Letting modern species exist from the oldest requested time slice can
inflate palaeo-distribution and richness estimates; use
`missing_origin = "oldest"` only for a stated lineage/clade summary or
backwards-compatible sensitivity run. If an extrapolation-weight table is
supplied, unmatched cell/time rows are treated as `Q_extrap = 0`; omitting the
whole extrapolation component is the only neutral case.

Deep-time projection helpers now fail early when key scientific inputs are
ambiguous. Requested paleoenvironment variables must exist in the cube,
`land_only = TRUE` requires the named landmask layer, every non-intercept Beta
axis must be present after recipe application, model-averaged BioGeoBEARS
weights are checked and normalised, and BioGeoBEARS LRT is restricted to the
standard nested +J comparisons unless custom nested models are explicitly
supplied with a warning. Matrix accessibility inputs are supported only when
their dimnames identify species/cells/regions and time or another projection
key.
`hee_combine_projection_models()` is strict by default: it refuses to label
M1-M5 as complete unless `L_land`, `E_phylo`, `A_BGB`, `D_static`, and
`D_dynamic` inputs are supplied. Set `strict = FALSE` only for transparent
formula demonstrations or sensitivity runs with intentionally neutral
components.

## Boundary Conditions

HmscEcoEvo does not infer ancestral ranges, does not fit BioGeoBEARS/BSM, and does not infer branch shifts, rate shifts, adaptive peaks, convergence regimes, or causal eco-evolutionary feedback paths from HMSC outputs alone. For those analyses, pass precomputed tables from the appropriate external model into `hmsc_ecoevo(precomputed = list(...))`.

Historical predictors and residual species associations are diagnostic summaries. They should not be interpreted as causal evidence without additional design, experiments, temporal data, or mechanistic modelling.

The full workflow report uses the same boundary language: colonisation,
extinction, rescue, connectivity, speciation-opportunity, reliability, and
M1-M5 dispersal contrasts are heuristic/proxy scenario diagnostics unless
calibrated with independent process data, and they are non-causal summaries of
the supplied model components.
