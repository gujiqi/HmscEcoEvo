# HmscEcoEvo 1.0.1

* Added the landscape-explicit active-dispersal API:
  `hee_dispersal_landscape_weights()`,
  `hee_dispersal_spherical_kernel()`,
  `hee_dispersal_transition_matrix()`,
  `hee_dispersal_spacetime_graph()`,
  `hee_dispersal_path_diagnostic()`, and
  `hee_dispersal_pruning_likelihood()`. The preferred P2 kernel is now
  area-aware spherical diffusion over sparse historical-grid edges. It can
  use an explicit resistance-expanded `effective_cost_km` once, while plate
  carriage remains a separate Dynamic-Earth transport operation.

* Reclassified `hee_dispersal_kernel()`, `hee_dispersal_edge_kernel()`, and
  `hee_dispersal_particle_kernel()` as legacy compatibility interfaces.
  Existing scripts remain valid; new work should use the spherical kernel and
  `hee_dispersal_propagule_particles()`. The eight-process catalog and README
  now distinguish active P2 movement, conditional P3 colonisation, Dynamic
  Earth, historical constraints, uncertainty, and derived biodiversity maps.

# HmscEcoEvo 1.0.0

* Added `hee_dispersal_particle_kernel()` and
  `hee_projection_global_grid_step()`. The first converts a continuous local
  movement kernel into a reproducible finite-propagule Monte Carlo
  approximation while preserving every source-cell kernel total; it is not a
  posterior-particle reconstruction or a gene-flow estimate. The second
  exposes the no-BioGeoBEARS global-grid CTMC update, keeping environmental
  filtering, arrival, conditional establishment, optional biotic terms, and
  local persistence in their respective process roles.

* Added `scripts/case05_plant200_eight_process_three_dispersal.R`. It runs a
  controlled Plant200 comparison with identical inputs and eight-process
  settings under pure-distance forward movement, continuous topographic
  movement, and finite-propagule topographic movement. It writes complete
  1-degree GIS GeoTIFF/PNG layers and fixed-scale difference maps only after
  calculating a shared legend contract across all scenarios and time slices.
  The case explicitly excludes BioGeoBEARS, BSM, regional accessibility and
  `A_hist`; it remains a process-constrained forward scenario unless modern
  endpoints and fossils are independently used to condition trajectories.

* Added the public topographic dispersal API: `hee_dispersal_resistance_surface()`,
  `hee_dispersal_connectivity_graph()`, `hee_dispersal_edge_kernel()`, and
  `hee_dispersal_arrival_pressure()`. These functions turn raw continuous
  palaeo-elevation into a land-only 4/8-neighbour structural-resistance layer,
  edge-specific effective movement costs, a local movement kernel, and source
  arrival pressure. The default uses local relief and elevation steps rather
  than absolute elevation; it excludes ocean edges and keeps climate suitability
  outside the resistance layer to avoid double-counting environmental filtering.
  The result is a generic topographic P2 scenario unless calibrated for the
  focal taxon, not an empirical gene-flow or realised-dispersal estimate.

* Added the direct ancestral environmental-response path:
  `hee_evolution_ancestral_response_direct()` reconstructs ancestral Beta
  response trajectories from HMSC tip Beta posterior draws on a dated tree, and
  `hee_environmental_filtering_ancestral_suitability()` projects those
  ancestral responses onto palaeoenvironmental cells. This prevents the
  scientifically unsafe shortcut of averaging modern prediction maps backward
  in time. HMSC intercepts are not reconstructed by default because they can
  encode prevalence, sampling, detection and accessibility rather than niche
  shape.

* Public API navigation now uses eight explicit biological process names:
  environmental filtering, dispersal, colonisation, biotic filtering,
  persistence, evolution, speciation, and extinction. Added
  `hee_core_process_catalog()` and `hee_function_catalog()` as the preferred
  machine-readable maps of canonical function names, legacy aliases, driver
  layers, historical-predictor functions, BioGeoBEARS/tree constraints,
  uncertainty layers, and result diagnostics. Older `hee_hmscee_*`,
  `hee_geoprocess_*`, `hee_dynamic_*`, `hee_lineage_*`, and hmscHist names are
  retained for compatibility, but new documentation should prefer the explicit
  eight-process names.

* Rebuilt the HmscEE core as a two-scale nested Dynamic Earth-Biota Assembly
  Framework. The former multiplicative equation
  `S_HMSC * A_hist * E_phylo * L_land * D_dynamic * Q_extrap` is no longer the
  main model because BioGeoBEARS/BSM regional range history and an HmscEE
  dynamic dispersal filter can otherwise explain the same dispersal chain
  twice.
* Fixed the theoretical contract around process responsibility. Climate,
  geology, geography, habitat history, BioGeoBEARS/BSM, dated trees, fossils,
  extrapolation and phyloregion summaries are now consistently documented as
  drivers, inference constraints, observations, uncertainty layers or derived
  diagnostics rather than additional biological processes.
* Added the current 1.0 core functions:
  `hee_paleo_earth_state()`, `hee_bsm_region_history()`,
  `hee_trait_mediated_ancestral_response()`, `hee_lineage_suitability()`,
  `hee_within_region_movement_kernel()`, `hee_arrival_pressure()`,
  `hee_establishment_probability()`, `hee_colonisation_from_arrival()`,
  `hee_persistence_probability()`, `hee_entry_initialization()`,
  `hee_nested_region_cell_occupancy()`, `hee_region_pool_from_cell_occupancy()`,
  `hee_extrapolation_uncertainty_report()`, and
  `hee_hmscee_process_catalog()` / `hee_hmscee_formula_catalog()`.
* Palaeogeographic/plate scenarios are now written `g`, palaeoclimate
  scenarios are written `c|g`, BioGeoBEARS/BSM histories are `h|g`, and
  HMSC/ancestral response posterior draws are `s`. The package no longer uses
  `e` for Earth scenario because BioGeoBEARS already uses parameter `e` for
  range contraction/local extinction. `hee_bsm_region_history()` accepts
  `geography_scenario`, `geography_weight`, and conditional
  `history_weight_given_geography`; `hee_paleo_earth_state()` accepts
  `climate_scenario` and `climate_weight_given_geography`;
  `hee_nested_region_cell_occupancy()` summarizes as
  `sum_g w[g] sum_h w[h|g] sum_c w[c|g] sum_s w[s] R^(h,g) q^(h,g,c,s)`
  instead of using an independent Cartesian product over geography, climate,
  BioGeoBEARS histories, and response draws.
* Added `hee_trait_mediated_ancestral_response()` to keep ancestral traits,
  ancestral environmental response, and ancestral geographic range separate.
  The implemented response formula is `beta[lineage,time] = Gamma *
  T[lineage,time] + u[lineage,time]`; BioGeoBEARS range columns are rejected as
  niche traits.
* `hee_dynamic_earth_biota_probability()`,
  `hee_dynamic_state_transition()`,
  `hee_regional_species_pool_transition()`, and
  `hee_static_snapshot_approximation()` are now explicit obsolete-core guards.
  They stop with a message directing users to
  `hee_nested_region_cell_occupancy()` instead of silently producing the old
  product-form probability.
* In the new core, BioGeoBEARS/BSM answers "which region", HMSC plus ancestral
  response reconstruction answers "which environment", and HmscEE answers
  "which cells inside an occupied region". Geography scenarios create hard
  habitat state `H`, cell-region maps, and movement resistance/least-cost paths
  `W`; climate scenarios create environmental predictors `X`.
  Extrapolation/no-analog scores are reported beside probability as
  uncertainty diagnostics, not multiplied into biological probability.
* Rebuilt Case 03 into separate entry points:
  `scripts/case03_geoprocess_realdata_past_to_present.R` is a strict real-data
  input-audit workflow with no simulated fallback, and
  `scripts/case03_geoprocess_simulated_demo.R` is an explicitly labelled nested
  teaching demo. The old `case03_geoprocess_past_to_present_simulated.R` file is
  a compatibility wrapper that runs the demo only.
* `hee_combine()` now recognises the 1.0 component aliases `S_HMSC`, `A_hist`,
  `E_phylo`, `L_land`, `G_arena`, `D_dynamic`, and `Q_extrap`, and reports
  explicit component columns in the combined projection table. Missing HMSC
  suitability in a supplied suitability table is now conservative
  (`probability = 0`) and flagged with `missing_suitability`.
* Static and dynamic dispersal/connectivity helpers now use kilometre
  great-circle distances for longitude-latitude coordinates instead of treating
  decimal-degree differences as distances. This can materially change
  `D_static`, `D_dynamic`, structural connectivity, regional filters, and
  downstream M4/M5 outputs compared with older selected-cell demonstrations.
* README and package-level documentation now describe HmscEcoEvo 1.0 as a
  Dynamic Earth-Biota Assembly Framework: Earth dynamics create/remove the
  arena, historical biogeography and connectivity constrain the regional pool,
  HMSC evaluates local establishment, and extrapolation diagnostics report
  transferability risk.
* Added palaeoenvironment adapter helpers for applying the framework directly
  to land-harmonized Phanerozoic environment factors:
  `hee_paleoenv_variable_roles()`,
  `hee_paleoenv_default_variable_map()`, and
  `hee_prepare_paleoenv_geoprocess_inputs()`. The defaults prioritise
  `land_mask_dem`, `land_area_km2`, `elevation_m`, `MAT_pohl_C`, and
  `MAP_pohl_mm_yr`; moisture, wetland, and bryophyte indices are treated as
  opportunity or explicit habitat-gating sensitivity layers rather than silent
  land multipliers.
* Case 03 is now aligned with the nested HmscEE 1.0 core. The simulated demo
  uses geography-specific `H/W`, climate-specific `X`, conditional BSM-like
  region histories `h|g`, trait-mediated ancestral responses `s`, and
  within-region movement in `hee_nested_region_cell_occupancy()`. The old M5
  product-form validation is no longer described as the core Case 03 result.

# HmscEcoEvo 0.2.0

* The Case 03 simulated demonstration now defaults to a study-scale design:
  36 species, 240 complete longitude-latitude grid cells, 100 spatially
  stratified modern HMSC sites, six regions, and 25 Myr intervals from 540 to
  0 Ma. It no longer probes a user-specific Windows palaeoenvironment path or
  silently changes between external and simulated environments. Fixed map-tile
  widths were replaced by grid-derived geometry. These changes alter demo
  figures, sampling coverage, and all derived demo-only summaries; they do not
  alter real Case 03 calculations.

* Case 03 now has separate real-data and simulated-demo entry points. The real
  entry point defaults to `analysis_mode=real`, has no generated biological or
  geographical fallback, writes `missing_required_inputs.csv` for absent core
  inputs, and records module-specific `NOT_RUN_*` states for missing dated-tree,
  region, BioGeoBEARS, plate, BSM, or fossil inputs. The former mixed workflow
  is retained as `scripts/case03_geoprocess_simulated_demo.R`; the old filename
  is a backwards-compatible wrapper.
* The real Case 03 uses full time-specific palaeo grids, locked 0 Ma predictor
  scaling, kilometre great-circle distances, DEM-neighbourhood topographic
  metrics, HMSC Beta posterior uncertainty summaries, actual task states and
  content hashes. These changes can materially alter results relative to the
  former selected-cell/posterior-mean demonstration.
* Case 03 real mode now writes GIS-ready single-species and diversity map
  products on the complete palaeogeographic grid: posterior mean/sd/quantile
  suitability, binary and historically accessible suitability, environment-
  suitable-but-inaccessible areas, expected and binary richness, Shannon,
  Simpson, weighted endemism, phylogenetic and functional diversity summaries,
  temporal transition maps, per-raster metadata sidecars, and all-time
  NetCDF/SpatRaster stack indexes. The quick workflow accepts
  `--quick_times=...` for explicit small-scale smoke tests, while full mode
  still processes every time slice supplied by the environment cube.
* Case 03 now applies externally supplied least-cost paths directly when they
  are available. Structural connectivity combines kilometre path distance,
  explicit route status, and barrier strength; otherwise the output is labelled
  as a region-centroid approximation. BSM events are diagnostic by default and
  no longer modify historical accessibility unless the user explicitly enables
  the documented heuristic sensitivity scenario.
* Independent fossil validation now matches the nearest available age first and
  then the nearest palaeogeographic grid cell by great-circle distance. Modern
  coordinates are rejected unless externally reconstructed. This can change
  fossil-validation probabilities relative to the former degree-distance
  matcher.
* Case 03 M5 no longer multiplies continuous BioGeoBEARS accessibility twice
  when initialising the oldest time slice or crossing a species-origin time.
  `D_dynamic` is now a binary external-support state controlled by
  `source_accessibility_threshold`, while continuous `A_BGB` remains the
  separate historical-accessibility term in M5. This can increase M5 relative
  to the former unintended `A_BGB^2` result where accessibility was between
  zero and one. Zero accessibility never seeds a source, even when the
  threshold is set to zero.
* Real Case 03 now refuses to derive `origin_ma` from a tree unless
  `--tree_time_unit=Ma` is explicitly declared. An external
  `species_origin.csv` remains preferred when species first-appearance or
  lineage-origin evidence is available. HMSC variance partitioning and each
  planned cross-validation analysis now receive explicit tabular/figure or
  failure-status outputs; full mode stops if a planned CV analysis fails.

## Final full-chain audit fixes

* Optional phylogenetic signal diagnostics now correctly extract numeric
  Pagel's lambda and Blomberg's K estimates from `phytools::phylosig()` list
  returns. This fixes a dependency-complete failure that was hidden on systems
  where `phytools` was not installed.
* DESCRIPTION now declares script/report helper packages used by the full
  workflow examples (`gridExtra`, `scales`, `pkgload`, `officer`, and `terra`)
  in `Suggests` for reproducible local validation.

## Geological-process safety fixes

* The Case 03 simulated demo includes a gen3sis-style observer layer. The demo writes
  tables and figures for dynamic-landscape mapping, ancestral initialization,
  dispersal/ecological-filter mapping, carrying-capacity and crowding proxies,
  environmental adaptation scores, species range-size observers,
  abundance-proxy weighted trait observers, and
  colonisation-extinction-speciation process balance. These are explicitly
  labelled as proxy/heuristic observers, not gen3sis engine outputs and not
  causal or true abundance/speciation estimates.
* Case 03 output files are now written with hierarchical numeric prefixes
  (`01.1`, `02.1`, `05.1.01`, etc.) and `--clean_output=TRUE` by default
  removes stale figures/tables from known output subdirectories before a run.
  Gen3sis-style origin seeding is restricted to each species' ancestral region,
  adaptation observers use standardized palaeoenvironmental optima, and the
  process-balance table now preserves a signed `net_assembly_balance` in
  `[-1, 1]` plus a separate positive-only opportunity proxy. This can change
  Case 03 observer trajectories and file paths compared with earlier outputs.
* Case 03 now fits a real `Hmsc::Hmsc()` model and samples it with
  `Hmsc::sampleMcmc()` on simulated modern communities extracted from the 0 Ma
  environment layer, then projects the fitted HMSC environmental suitability
  through the 540-0 Ma time axis. When the local 4-degree Phanerozoic
  environment RDS is available it is used directly; otherwise the script falls
  back to a clearly labelled simulated environment for portable tests.
* Case 03 now writes explicit input-inventory, process-source, internal-helper
  coverage, palaeoenvironment-role, region-environment, and final M5
  probability/richness map outputs. Region connectivity is derived from
  palaeoenvironmental summaries (regional land area, land existence,
  elevation, temperature, precipitation, habitat and disturbance) rather than
  a disconnected artificial route table. External BioGeoBEARS/accessibility,
  plate-point, and BSM/event tables can be supplied with command-line
  arguments; supplied but missing files now stop with clear errors.
* Case 03 M4/M5 dispersal filters now multiply the static/dynamic distance
  filter by structural-functional-climatic geological connectivity, and M5
  includes an explicit origin-time seed controlled by `--origin_seed`. This can
  change previous Case 03 M4/M5 richness and area because connectivity is now
  part of the projection formula rather than only a diagnostic side table.
* `hee_ecological_opportunity()` now masks opportunity to zero when supplied
  geography or habitat availability is absent or unknown. `hee_colonisation_
  probability()` now forces colonisation probability to zero when suitability,
  source pressure, or supplied accessibility is explicitly zero. Missing
  palaeodistance in `hee_colonisation_pressure()` is treated as no route rather
  than a zero-distance route.
* `hee_geoprocess_diagnostics()` now recomputes total
  `speciation_opportunity` after the allopatric, founder, in-situ and
  radiation opportunity subcomponents are updated, so the reported total and
  component columns are internally consistent.
* `hee_extrapolation_risk()` now treats palaeo rows with missing environmental
  predictors as high-risk/no-analog rows rather than reliable rows. This avoids
  overconfidence when palaeoenvironmental layers are incomplete and can reduce
  projections when `Q_extrap` is used as a weighting or exclusion scenario.
* `hee_phylo_time_mask()` now requires explicit `species_origin` by default.
  The previous fallback assumed every supplied species existed from the oldest
  requested time slice, which can overestimate species-level deep-time
  distributions. Set `missing_origin = "oldest"` only for backwards-compatible
  neutral demonstrations, lineage/clade summaries, or sensitivity scenarios
  where this assumption is stated.
* `hee_combine()` now treats unmatched rows in a supplied `extrapolation_weight`
  table as conservative (`Q_extrap = 0`) rather than fully reliable
  (`Q_extrap = 1`). Omitting the whole extrapolation component remains neutral.
  This can make old downweighted projections smaller in cells or times where
  extrapolation risk was not evaluated.
* Full workflow Word reports now state explicitly that colonisation,
  extinction, rescue, connectivity, speciation-opportunity, reliability, and
  M1-M5 dispersal contrasts are heuristic/proxy scenario diagnostics and
  non-causal unless calibrated with independent process data. This changes only
  report wording, not numerical results.
* `hee_combine()` now applies hard-zero constraints after component lookup:
  if phylogenetic existence, land/geographic existence, accessibility, static
  dispersal, dynamic dispersal, or extrapolation weight is explicitly zero, the
  final probability is zero even when environmental suitability is `NA`.
* `hee_combine()` now rejects species-level accessibility or phylogenetic-mask
  tables that omit `species` while projecting multiple species. Set
  `allow_species_broadcast = TRUE` only for an explicit shared neutral
  scenario.
* `hee_accessibility_to_cells()` now validates already cell-level accessibility
  tables, rejects duplicate `species + cell_id + time_ma` keys, clips values to
  `[0, 1]`, and treats explicit `NA` accessibility as zero evidence.
* `hee_combine_projection_models()` is strict by default and requires the
  formula components needed to label M1-M5 as complete (`L_land`, `E_phylo`,
  `A_BGB`, `D_static`, and `D_dynamic`). Use `strict = FALSE` only for
  transparent neutral-component demonstrations.
* Formula/default audit update: `hee_colonisation_pressure()` and
  `hee_rescue_effect()` now exclude `from_region == to_region` self-links by
  default (`include_self = FALSE`). Local persistence is handled by the
  occupancy survival term, not counted as external colonisation or demographic
  rescue. Set `include_self = TRUE` only for backwards-compatible scenarios
  where within-region source pressure is intentional.
* `hee_extinction_probability()` and `hee_extinction_risk()` now treat explicit
  `NA` values in `geographic_existence` conservatively as absent habitat and
  force local extinction. Omitting the entire geography column remains neutral.
* `hee_build_connectivity_cube()` now derives passability as
  `1 - barrier_strength` when a barrier table supplies only barrier strength.
  Explicitly supplied but missing `barrier_passability` remains conservative
  (`0`).
* Static dispersal helpers now return zero reachability when explicitly
  supplied accessibility is zero for all candidate source regions, instead of
  selecting a "best" zero-accessibility source.
* Speciation-opportunity helpers now use conservative defaults: missing
  fragmentation, persistence, ecological opportunity, land age, or connectivity
  evidence no longer creates positive opportunity indices. The new
  `speciation_opportunity_lambda` column is the preferred name; the old
  `speciation_rate_lambda` column is retained as a backwards-compatible alias
  and should not be interpreted as a calibrated speciation rate.
* `hee_bgb_accessibility()` now validates direct data-frame inputs, rejects
  duplicate accessibility keys, checks non-negative `time_ma`, and clips
  accessibility to `[0, 1]`. This can change old results that contained
  duplicated or out-of-range BioGeoBEARS accessibility rows.
* `hee_bgb_compare()` now rejects duplicated BioGeoBEARS model names before
  computing AICc weights, preventing duplicated candidate models from changing
  model-averaged accessibility.
* The Case 1 and Case 2 workflow scripts now accept both `--key=value` and
  `key=value` command-line arguments. Case 2 region joins now reject duplicated
  `cell_id + time_ma` keys instead of silently expanding projection rows.
* Case 1 M1-M5 projection probabilities now use the package-level
  `hee_combine()` formula inside the tiled projection loop so example outputs
  follow the same masking and probability-bound rules as package users.
* `hee_dynamic_assembly()` now treats missing colonisation and extinction
  tables as neutral (`gamma = 0`, `epsilon = 0`). Previous versions used
  suitability-derived defaults, which could increase or decrease occupancy even
  when a process table was absent.
* `hee_colonisation_probability()` and `hee_extinction_probability()` now reject
  negative or non-finite `delta_t` values instead of silently replacing them
  with zero.
* `hee_colonisation_probability()` now treats an omitted, zero, or unmatched
  `source_pressure` term as no propagule evidence and yields zero colonisation
  probability. This is a conservative behaviour change: source pressure must be
  supplied explicitly when colonisation should occur.
* `hee_extinction_probability()` and `hee_extinction_risk()` now force local
  extinction when `geographic_existence <= 0`, matching the deep-time
  interpretation that a taxon cannot persist in a non-existent arena.
* `hee_combine()` now uses only explicit process columns for land,
  accessibility, static dispersal, dynamic dispersal, and extrapolation
  weights. It no longer guesses a weight from the first numeric column.
* Process joins now check duplicate keys and row-count changes. Ambiguous
  species/cell/region/time joins error with the offending key values.
* `hee_functional_connectivity()` now errors on duplicated trait species and
  requires explicit `expand_species = TRUE` when expanding non-species
  connectivity to all species.
* `hee_prediction_reliability()` now records raw extrapolation risk, applies a
  saturating transform to unbounded risk scores, supports component weights, and
  treats missing optional components as neutral.
* Single-row and constant-richness inputs now produce finite refugia scores.
* `hee_combine()` now honours named matrix accessibility or weight inputs
  instead of treating matrices as neutral weights. Matrices must have
  interpretable row/column names, for example species by time (`0Ma`) or
  cell/region by time.
* `hmsc_ecoevo()` now reorders `Gamma` by trait and Beta-axis names, matching
  the stricter handling already used for `gamma_draws`.
* `hee_phylo_time_mask()` now correctly names unnamed `species_origin` vectors
  when a matching `species` vector is supplied, and rejects invalid ages.
* `hee_make_paleo_grid()`, `hee_extract_paleoenv_track()`, and HMSC projection
  helpers now reject missing requested paleo variables or missing Beta axes
  instead of silently dropping them.
* `hee_bgb_accessibility()` now checks model-average weights, normalises them
  to sum to one, rejects duplicate accessibility keys, and clips final
  accessibility to `[0, 1]`.
* `hee_bgb_lrt()` now checks standard nested BioGeoBEARS pairs
  (`DEC`/`DEC+J`, `DIVALIKE`/`DIVALIKE+J`, `BAYAREALIKE`/`BAYAREALIKE+J`) and
  rejects invalid likelihood-ratio tests.
* Static dispersal helpers now reject duplicated accessibility/barrier keys
  that could otherwise change row counts during joins.
* `hee_geoprocess_model_comparison()` now requires a named list of model
  outputs so every summary has a stable `model_id`.
* Formula/model-definition audit fixes: `hee_functional_connectivity()` is now
  neutral when no trait table is supplied; species-specific connectivity now
  errors when a species is missing from the trait table; explicit zero-area rows
  now force local geographic extinction; extreme speciation-rate proxy values are
  clamped to finite rates; undefined small-sample AICc weights now error;
  extrapolation thresholds are applied once; and dynamic dispersal defaults are
  used only for the first time slice, not for later missing source evidence.
* Supplied projection-component tables are now conservative when incomplete:
  unmatched accessibility, landmask, static-filter, dynamic-filter, and
  phylogenetic-mask rows are treated as zero evidence rather than full
  accessibility. Explicit `NA` land or composite static-dispersal weights are
  also conservative. Omitting the whole optional component remains neutral.
* Tip-level local phylogenetic signal now keeps Moran/LIPA values matched to
  species names even when the phylogeny tip order differs from the Beta matrix
  row order.
* `hee_update_species_pool()` now honours its documented time rule: when
  `initial_pool` contains `time_ma`, the oldest row per region is used as the
  initial pool, independent of input row order.
* `hee_compare_plate_corrected()` now requires one-to-one `point_id` matches
  and rejects duplicate plate-comparison keys before computing displacement.
* `calc_niche_metrics()` now reports quadratic niche breadth as the
  Gaussian-equivalent `sqrt(-1/(2 * beta_quadratic))`, consistent with
  `eta = const - (x - optimum)^2 / (2 * sigma^2)`.
* `hee_ecological_opportunity()` now ignores and warns on unknown weight names
  instead of adding them to the denominator and depressing opportunity scores.
* `hee_geoprocess_diagnostics()` now applies the same conservative unmatched-row
  logic for supplied accessibility and phylogenetic masks as `hee_combine()`;
  its internal metric lookup no longer assumes a pre-existing `.row_id` column.
* Connectivity functions now reject negative palaeodistance or environmental
  distance values instead of clipping them to zero and inflating connectivity.
* Supplied metric tables now need explicit recognised value columns; HmscEcoEvo
  no longer guesses from the first numeric column. Multi-row metric tables with
  no shared join keys now error instead of being silently ignored.
* `hee_extrapolation_weight(mode = "downweight")` now treats `Inf`
  extrapolation score as extreme risk (`Q_extrap = 0`) while missing, negative,
  or `NaN` scores remain neutral zero risk.
* `hee_dynamic_assembly()` now accepts species-by-time matrix `phylo_mask`
  inputs as well as long tables, and recognises `phylo_existence` as a lineage
  mask alias. Supplied but unmatched matrix/table rows are treated as lineage
  absence.
* `hee_dynamic_assembly()` now converts source-pressure-only colonisation
  inputs through `hee_colonisation_probability()` using the internal
  `interval_myr`, instead of using `suitability * source_pressure` directly.
* `hee_colonisation_probability()` now treats unmatched or explicit `NA`
  accessibility rows as not historically accessible when an accessibility table
  is supplied. Omitting accessibility remains neutral.
* `hee_update_occupancy()`, `hee_combine_static()`, and `hee_combine()` now use
  conservative values for explicit or unmatched geographic/phylogenetic masks.
  Omitted optional masks remain neutral.
* Structural, functional, static, and dynamic connectivity helpers now treat
  explicit missing barriers/routes/passability as conservative rather than
  fully passable. Users should provide explicit `1` values for routes known to
  be open.
* Dispersal-scale inputs are now checked for finite non-negative values in
  static/dynamic dispersal helpers and dispersal matrices. A scale of zero is
  allowed and means no movement across positive distances.
# HmscEcoEvo 1.0.4

* The current HmscEE process and formula catalogs now name six biological
  processes: environmental filtering, dispersal, biotic filtering, evolution,
  speciation, and extinction. Colonisation and persistence are no longer
  listed as independent core processes. Their older exported helpers remain
  compatibility functions for labelled, calibrated or scenario-based
  dynamic-occupancy examples, not the empirical Case05 workflow.
* Process-attribution state budgets now label retained occupancy, new
  occupancy, recolonisation, and local loss as state-transition components
  rather than extra biological processes.
