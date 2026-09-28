# HmscEcoEvo geoprocess final audit

Audit date: 2026-07-27  
Package version: 0.2.0  
Workspace: `C:/Users/Google/Documents/HMSC-HIST/HmscEcoEvo`

This document records the final acceptance audit for the geological-process,
deep-time projection, dynamic assembly, connectivity, colonisation, extinction,
speciation-opportunity, and prediction-reliability layer of HmscEcoEvo.

## 0A. 2026-08-11 formula and function audit addendum

This addendum records the follow-up audit of whether process formulas are
actually embedded in package calculations and whether defaults can alter
scientific conclusions.

### 0A.1 Critical fixes

| Item | Old behaviour | New behaviour | Files |
|---|---|---|---|
| Colonisation without source pressure | Missing `source_pressure` was treated as `M = 1`, so suitability alone could produce colonisation. | Missing, explicit zero, or unmatched `source_pressure` is no propagule evidence and gives `colonisation_probability = 0`. | `R/geological_processes.R`, `man/hee_colonisation_probability.Rd`, `NEWS.md` |
| Negative palaeodistance or climate distance | Negative distances were clipped or transformed in ways that could inflate connectivity. | Negative finite palaeodistance/environment distance now errors before connectivity is calculated. | `R/geological_processes.R` |
| Supplied process tables with wrong value column | Some tables could be silently ignored or replaced by defaults if they lacked a recognised value column. | Supplied component/metric tables now require recognised columns or error with the table name. | `R/deep_time_workflow.R`, `R/geological_processes.R` |
| Infinite extrapolation risk | `Inf` extrapolation score was treated as neutral zero risk in `hee_extrapolation_weight()`. | `Inf` is extreme extrapolation risk; in `mode = "downweight"` it yields `Q_extrap = 0`. | `R/combine_projection.R`, `man/hee_extrapolation_weight.Rd` |

### 0A.2 Formula status after fix

| Process | Formula embedded in calculation? | Interpretation |
|---|---|---|
| Source pressure | Yes: `M[j,c,t] = 1 - product(1 - p[j,c',t] * K[d] * C[j,c',c,t])`. | Proxy for propagule pressure from previous occupancy and connectivity. |
| Colonisation | Yes: if `M <= 0` or `A <= 0`, rate is zero; otherwise `gamma = 1 - exp(-exp(eta_col) * delta_t)`. | Heuristic/calibratable transition probability, not direct observed colonisation. |
| Extinction | Yes: `epsilon = 1 - exp(-exp(eta_ext) * delta_t)`, with geographic absence forcing `epsilon = 1`. | Heuristic/calibratable transition probability. |
| Dynamic occupancy | Yes: `p_next = G * E * (p_prev * (1 - epsilon) + (1 - p_prev) * gamma)`. | Ensures no occupancy before lineage origin or outside geographic existence. |
| Structural connectivity | Yes, distance/passability decay with non-negative distances only. | Process proxy requiring palaeogeographic input. |
| Functional connectivity | Yes, species traits modify distance decay; species expansion must be explicit. | Trait-mediated proxy, not measured dispersal unless calibrated. |
| Climatic connectivity | Yes, environmental mismatch reduces connectivity. | Corridor proxy unless real corridor reconstructions are supplied. |
| Prediction reliability | Yes, combines extrapolation, model agreement, phylo validity, and calibration with explicit weights. | Reliability index, not a confidence interval. |

### 0A.3 Verification

- Function/export scan: 191 public-prefix function definitions, 191 exports,
  0 duplicate definitions, 0 `.hee_*` internal exports, 0 exported functions
  without definitions.
- Targeted tests passed: `test-geological-processes.R` 96 pass,
  `test-geological-processes-systematic.R` 128 pass, `test-combine.R` 30 pass,
  `test-projection-components.R` 31 pass.
- Full `testthat::test_local()` passed with 2 expected skips for case-output
  tests that are not installed in the package test context.
- `roxygen2::roxygenise()` completed and updated the affected Rd files.
- `R CMD build .` produced `HmscEcoEvo_0.2.0.tar.gz`.
- `_R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual
  HmscEcoEvo_0.2.0.tar.gz` completed with `Status: OK`; missing Suggests were
  `phytools`, `geiger`, `phylosignal`, and `phylobase`.
- Quick Case 1 rerun after formula fixes succeeded at
  `outputs/_audit_case01_quick_after_formula_fix`, producing 173 PNG figures,
  173 PDF figures, 141 CSV tables, 508 ordered-result files, and the Word
  report `reports/hmscecoevo_full_540Ma_workflow.docx`.

### 0A.4 Remaining interpretation limits

- The geoprocess transition formulas are configurable heuristic/process-proxy
  models unless coefficients are calibrated with independent fossil,
  experimental, or palaeobiogeographic evidence.
- Source pressure, opportunity, reliability, and limitation maps should not be
  reported as causal estimates without explicit model calibration.
- Quick-mode case outputs are software verification outputs. Scientific
  interpretation should use `quick = FALSE`; the 2026-08-11 formula changes
  mean older case outputs should be regenerated before interpretation.

## 0. 2026-07-27 scenario-suite extension audit

This addendum records the current implementation round based on
`地学过程-20260727.docx`. It supersedes older output counts elsewhere in this
file where they differ.

### 0.1 Newly added public functions

| Function | Purpose | Inputs | Required columns | Outputs | Interpretation boundary |
|---|---|---|---|---|---|
| `hee_geoprocess_scenario_suite()` | Runs deterministic normal and stress scenarios for the geological-process layer. | Scenario names, species names, seed. | Built-in generated tables contain `species`, `cell_id`, `region`, `time_ma`, `geographic_existence`, `suitability`, `probability`, connectivity, accessibility, phylo mask, traits, and pseudo BSM events. | `input_tables`, `result_tables`, `validation_table`, `scenario_summary`, `formula_catalog`, `interpretation_table`. | Software-audit and scientific-logic validation; not a real palaeogeographic estimate. |
| `plot_geoprocess_scenario_suite()` | Plots scenario means, validation status, landscape events, and dynamic occupancy. | Output from `hee_geoprocess_scenario_suite()`. | Suite result tables. | Named printable ggplot list. | Visual QA of process logic and edge cases. |

Internal helpers added in `R/geoprocess_scenarios.R`: `.hee_run_one_geoprocess_scenario()`,
`.hee_geoprocess_scenario_inputs()`, `.hee_geoprocess_scenario_derived()`,
`.hee_validate_geoprocess_scenario()`, `.hee_geoprocess_formula_catalog()`,
`.hee_geoprocess_interpretation_table()`, and `.hee_bind_fill()`. These are not
exported.

### 0.2 Scenario coverage

The suite covers: normal, extreme, missing optional inputs, duplicate key,
single time slice, uneven time intervals, shuffled time order, and sudden land
appearance/loss. The duplicate-key case is expected to stop with a clear
duplicate-key error and is recorded as `EXPECTED_ERROR`; it is not silently
repaired.

Current Case 1 scenario summary:

| Scenario | Status | Validation failures |
|---|---|---|
| normal | PASS | 0 |
| extreme | PASS | 0 |
| missing_optional | PASS | 0 |
| duplicate_key | EXPECTED_ERROR | 0 |
| single_time | PASS | 0 |
| uneven_time | PASS | 0 |
| shuffled_time | PASS | 0 |
| land_appearance_loss | PASS | 0 |

### 0.3 Formula catalog added to outputs

The new `geoprocess_formula_catalog.csv` records the formulas used by the
scenario suite:

- Static projection M1: `P = S_HMSC * L_land`.
- Static projection M2: `P = S_HMSC * E_phylo * L_land`.
- Static projection M3: `P = S_HMSC * E_phylo * A_BGB * L_land`.
- Static projection M4: `P = S_HMSC * E_phylo * A_BGB * D_static * L_land`.
- Static projection M5: `P = S_HMSC * E_phylo * A_BGB * D_dynamic * L_land`.
- Source pressure: `M[j,c,t] = 1 - product(1 - p[j,c',t] * K[d] * C[j,c',c,t])`.
- Colonisation: `gamma = 1 - exp(-exp(eta_col) * delta_t)`.
- Extinction: `epsilon = 1 - exp(-exp(eta_ext) * delta_t)`.
- Dynamic occupancy: `p[t+1] = G * E * (p[t] * (1 - epsilon) + (1 - p[t]) * gamma)`.
- Structural connectivity: `C_struct = route_open * exp(-alpha * distance) * exp(-barrier_weight * barrier)`.
- Functional connectivity: `C_func = plogis(k0 + k_trait * trait - k_distance * distance - k_barrier * barrier)`.
- Climatic connectivity: `C_clim = exp(-phi * climate_mismatch)` or supplied corridor suitability.
- Prediction reliability: configurable multiplicative/weighted diagnostic; not
  a credible interval.

### 0.4 Fixes made in this round

| Issue | Old behavior | New behavior | Scientific effect |
|---|---|---|---|
| Ecological-opportunity name collision | If `env_cols` included `habitat_heterogeneity`, an internal merge created `.x/.y` columns and then attempted to assign a zero-length vector. | Internal computed heterogeneity now uses `.hee_habitat_heterogeneity` and removes it before return. | Prevents normal geoprocess diagnostics from failing when user supplies habitat heterogeneity as an environmental variable. |
| Zero-row scenario tables | Missing optional inputs could produce 0-row tables; adding scalar `scenario` failed. | Scenario collectors now add 0-length scenario columns and `.hee_bind_fill()` pads missing columns with `rep(NA, nrow(x))`. | Missing-data stress tests are explicit and auditable. |
| Refugia score validation | `refugia_score` was treated like a probability and negative standardized scores failed validation. | `refugia_score` is finite-only; probability-like columns remain constrained to `[0, 1]`. | Avoids false scientific failures for standardized refugia scores. |
| Landscape appearance/loss evidence | Scenario suite did not expose the package's event detector in result tables. | `hee_landscape_events()` is now run inside the suite and written as `landscape_events`. | Sudden emergence and submergence are tested and reported. |

### 0.5 Current verification results

| Check | Result |
|---|---|
| Code consistency scan | PASS: 163 public `hee_*`/scenario plot functions, 0 duplicate public definitions, 0 missing exports, 0 exports without definitions, 0 exported functions missing Rd aliases. |
| `roxygen2::roxygenise('.')` | PASS. |
| Targeted new scenario tests | PASS: 39 assertions. |
| Geological-process tests | PASS: `test-geological-processes.R` and `test-geological-processes-systematic.R`. |
| Full `testthat::test_local('.', reporter='summary')` | PASS with 2 expected skips for non-installed external case-output context. |
| Case 1 full-axis verification run | PASS: `quick = FALSE`, 109 time slices, 30 species, 80 sites, 157 PNG figures, 121 CSV tables, Word report generated. MCMC was intentionally shortened to `samples = 10`, `transient = 10` for code verification. |
| Case 2 missing-file behavior | PASS as expected: stops with `Missing required file: data_raw/comm.csv`. |
| `R CMD build .` | PASS: `HmscEcoEvo_0.2.0.tar.gz` created. |
| `_R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual HmscEcoEvo_0.2.0.tar.gz` | PASS, Status OK. INFO: suggested packages unavailable for checking: `phytools`, `geiger`, `phylosignal`, `phylobase`. |

### 0.6 Current output locations

- Case 1 output root:
  `C:/Users/Google/Documents/HMSC-HIST/HmscEcoEvo/outputs/case01_simulated_full_540Ma`.
- Word report:
  `outputs/case01_simulated_full_540Ma/reports/hmscecoevo_full_540Ma_workflow.docx`.
- Ordered output folders:
  `outputs/case01_simulated_full_540Ma/ordered_results`.
- Scenario RDS:
  `outputs/case01_simulated_full_540Ma/rds/geoprocess_scenario_suite.rds`.
- Scenario figures:
  `geoprocess_scenario_process_heatmap.png`,
  `geoprocess_scenario_validation_status.png`,
  `geoprocess_scenario_landscape_events.png`, and
  `geoprocess_scenario_dynamic_occupancy.png`.
- Scenario tables:
  `geoprocess_scenario_summary.csv`,
  `geoprocess_scenario_validation.csv`,
  `geoprocess_module_evidence_index.csv`,
  `geoprocess_formula_catalog.csv`,
  `geoprocess_interpretation_table.csv`,
  12 `geoprocess_scenario_input_*.csv` files, and
  22 `geoprocess_scenario_result_*.csv` files.

### 0.7 Remaining scientific limits

- The scenario suite is simulated. It validates process logic, boundary
  handling, and data-join safety, but it is not evidence for a real geological
  history.
- Plate positions, palaeodistance, palaeobarriers, palaeoclimate corridors,
  BioGeoBEARS accessibility, and fossil validations remain external-data
  dependent for real studies.
- Colonisation, extinction, rescue, and speciation-opportunity formulas are
  transparent configurable process models or proxy indices unless calibrated
  against independent data.
- Prediction reliability is an epistemic diagnostic and must not be reported as
  posterior probability, credible interval, or frequentist confidence.
- Deep-time species-level maps become increasingly lineage/clade-level
  hypotheses with age; the report explicitly flags this limitation.

## 1. Verification commands

The following commands were executed from the package root.

| Check | Command | Result |
|---|---|---|
| Function/export/Rd consistency | Custom R scan of `R/*.R`, `NAMESPACE`, and `man/*.Rd` | PASS: 209 `hee_*` definition rows, 161 public `hee_*` functions, 199 exports, 0 duplicate definitions, 0 accidental `.hee_*` exports, 0 public definitions missing from `NAMESPACE`, 0 exports without definitions, 0 exports missing Rd aliases. |
| Documentation generation | `roxygen2::roxygenise(".")` | PASS. |
| Case 1 quick workflow | `Rscript scripts/case01_simulated_full_540Ma.R --quick=TRUE --output=outputs/final_acceptance_case01_quick` | PASS: HMSC MCMC executed, Word report produced, ordered outputs written. |
| Time-cube full-axis check | `hee_time_axis_from_env_cube()` plus `hee_assert_full_time_axis(..., quick = FALSE, min_time_slices = 50, required_range = c(540, 0))` | PASS: 109 cube time slices, range 540 to 0 Ma. |
| Geoprocess targeted tests | `testthat::test_file()` for geoprocess, systematic geoprocess, combine, richness | PASS. |
| Full tests | `testthat::test_local(reporter = "summary")` | PASS with 2 expected skips for non-installed external case-output context. |
| Build | `R CMD build .` | PASS: `HmscEcoEvo_0.2.0.tar.gz` created. |
| Check | `_R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual HmscEcoEvo_0.2.0.tar.gz` | PASS: `Status: OK`. Suggested packages unavailable for checking: `phytools`, `geiger`, `phylosignal`, `phylobase`. |
| Case 2 missing-file behavior | `Rscript scripts/case02_realdata_full_540Ma_template.R --quick=TRUE ...` | PASS as expected: stopped with `Missing required file: data_raw/comm.csv`. |

## 2. Output acceptance

Clean quick-case output directory:
`outputs/final_acceptance_case01_quick`

| Output check | Result |
|---|---|
| PNG figures | 139 |
| CSV tables | 82 |
| PNG files smaller than 5 KB | 0 |
| Empty or unreadable CSV files | 0 |
| Word report | `outputs/final_acceptance_case01_quick/reports/hmscecoevo_full_540Ma_workflow.docx` |
| Word report size | 2,496,478 bytes |
| Embedded media in Word report | 22 |
| Projection task table | 5 rows, all `done` |
| Quick time range | 540 to 0 Ma |
| Quick time slices | 10 |

Existing full-case output directory:
`outputs/case01_simulated_full_540Ma`

| Full output check | Result |
|---|---|
| PNG figures | 138 |
| CSV tables | 58 |
| Word report | `outputs/case01_simulated_full_540Ma/reports/hmscecoevo_full_540Ma_workflow.docx` |
| Word report size | 6,348,176 bytes |
| Full time slices | 109 |
| Species | 30 |
| Sites | 80 |
| quick | FALSE |

## 3. Function inventory

### 3.1 Geological arena and opportunity

| Function | Purpose | Required structure | Key return columns | Time assumption | Status |
|---|---|---|---|---|---|
| `hee_define_geoprocess_models()` | Defines named M1-M5 projection scenarios and process formulas. | None. | `model_id`, `formula`, `required_components`. | Not time-indexed. | PASS. |
| `hee_landscape_events()` | Detects land/habitat/area/elevation/component changes through time. | Cell or region table with `time_ma` and at least one state column. | `landscape_event`, `area_change`, `habitat_change`. | Ma; sorted internally from old to young for interval comparisons. | PASS. |
| `hee_land_age()` | Computes persistence age of land/habitat state. | Landscape state through `time_ma`. | `land_age_myr`. | Ma; adjacent absolute intervals. | PASS. |
| `hee_ecological_opportunity()` | Builds a bounded opportunity index from land age, novelty, environmental change, and optional area. | Cell/region-time table plus selected environmental columns. | `ecological_opportunity`. | Ma only for inherited time columns; no propagation. | PASS. |

Formula notes:

- `hee_ecological_opportunity()` is a configurable index. Inputs already in
  `[0, 1]` are kept; non-negative intensities use stable positive-index
  transforms so a single positive value does not collapse to zero.
- This is a proxy for ecological opportunity, not a direct estimate of
  diversification or speciation rate.

### 3.2 Connectivity

| Function | Purpose | Required structure | Key return columns | Join behavior | Status |
|---|---|---|---|---|---|
| `hee_structural_connectivity()` | Physical route permeability from palaeodistance, barriers, and route openness. | Pair table with distance column; optional `from_region`, `to_region`, `time_ma`. | `structural_connectivity`. | Checks duplicate `from_region`-`to_region`-`time_ma` keys when present. | PASS. |
| `hee_functional_connectivity()` | Species-specific permeability from dispersal traits plus route properties. | Connectivity table plus optional trait table with unique `species`. | `functional_connectivity`, `connectivity`. | If connectivity already has `species`, only matching species are joined; no uncontrolled expansion. | PASS. |
| `hee_climatic_connectivity()` | Climatic traversability from environmental distance, corridor suitability, or precomputed climate connectivity. | Connectivity table. | `climate_connectivity`, `connectivity`. | Preserves rows; missing optional climate component is neutral (`1`). | PASS. |
| `hee_build_connectivity_cube()` / `hee_connectivity_cube()` | Builds structural, functional, climatic, and total connectivity across times. | Distances or landscape-derived region coordinates; optional traits. | Pair/species/time connectivity cube. | Species expansion must be explicit when needed. | PASS. |
| `hee_isolation_history()` | Summarises isolation through time from connectivity. | Connectivity cube. | `isolation`, `isolation_duration`. | Ma; sorted internally. | PASS. |

Formula notes:

- Structural connectivity:
  `C_struct = route_open * exp(-alpha * max(paleodistance, 0)) * exp(-barrier_weight * barrier_strength)`, clipped to `[0, 1]`.
- Functional connectivity:
  `C_func = plogis(k0 + k_trait * trait_index - k_distance * distance_index - k_barrier * barrier)`.
- Climatic connectivity:
  `C_clim = exp(-phi * max(environment_distance, 0))`, or supplied corridor suitability/precomputed connectivity, clipped to `[0, 1]`.
- Total connectivity multiplies available structural, functional, and climatic terms and is clipped to `[0, 1]`.

### 3.3 Colonisation, rescue, extinction, and dynamic assembly

| Function | Purpose | Required structure | Key return columns | Time assumption | Status |
|---|---|---|---|---|---|
| `hee_colonisation_pressure()` | Aggregates source occupancy through connectivity. | Previous occupancy plus pair/species/time connectivity. | `source_pressure`. | Uses explicit keys; duplicate keys error. | PASS. |
| `hee_rescue_effect()` | Computes rescue from occupied sources through connectivity. | Previous occupancy plus connectivity. | `rescue_effect`. | Uses explicit keys; duplicate keys error. | PASS. |
| `hee_colonisation_probability()` | Converts suitability, source pressure, accessibility, and opportunity to interval colonisation probability. | Suitability table or vector; optional source/accessibility/opportunity tables. | `colonisation_rate_lambda`, `colonisation_probability`. | `delta_t` must be finite and non-negative. | PASS. |
| `hee_extinction_risk()` | Bounded instantaneous extinction risk from unsuitability, land loss, isolation, area, rescue, disturbance, and vulnerability. | Species/cell/region-time table. | `extinction_probability`, `forced_geographic_extinction`. | Not interval-scaled. | PASS. |
| `hee_extinction_probability()` | Converts extinction drivers to interval extinction probability. | Species/cell/region-time table. | `extinction_rate_lambda`, `extinction_probability`, `forced_geographic_extinction`. | `delta_t` must be finite and non-negative. | PASS. |
| `hee_update_occupancy()` | One-step occupancy update. | Previous probability plus colonisation/extinction/geographic/phylo masks. | Occupancy probability. | One interval. | PASS. |
| `hee_dynamic_assembly()` | Propagates occupancy from old to young time slices. | Suitability plus optional previous state, colonisation, extinction, masks. | `probability`, process components. | Sorts `time_ma` old to young; repeated times within species-cell error. | PASS. |

Formula notes:

- Colonisation linear predictor:
  `eta = b0 + bS * logit(S) + bM * log(M) + bA * log(A) + bO * O`.
- Colonisation interval probability:
  `gamma = 1 - exp(-exp(eta) * delta_t)`.
- Omitted, explicit zero, or unmatched `source_pressure` means no propagule
  evidence and sets colonisation rate to zero. Source pressure must be supplied
  explicitly when colonisation should occur.
- Extinction linear predictor:
  `eta_ext = b0 + bU * (1 - S) + bL * land_loss + bI * isolation + bA * area_inverse + bR * rescue + bD * disturbance + bV * vulnerability`.
- Extinction interval probability:
  `epsilon = 1 - exp(-exp(eta_ext) * delta_t)`.
- `land_loss >= 1` or `geographic_existence <= 0` forces `epsilon = 1`.
- Dynamic occupancy update:
  `p[t+1] = G * E * (p[t] * (1 - epsilon) + (1 - p[t]) * gamma)`,
  where `G` is geographic existence and `E` is the phylogenetic time mask.

Interpretation:

- `gamma` and `epsilon` are scenario probabilities unless coefficients are
  calibrated externally.
- They are not automatically inferred colonisation or extinction rates.

### 3.4 Speciation and species pools

| Function | Purpose | Required structure | Key return columns | Scientific interpretation | Status |
|---|---|---|---|---|---|
| `hee_speciation_opportunity()` | General speciation-opportunity index. | Region/species/time drivers such as isolation, area, ecological opportunity. | `speciation_opportunity`. | Proxy only. | PASS. |
| `hee_allopatric_speciation_opportunity()` | Opportunity emphasising isolation/barriers. | Regional isolation/barrier drivers. | Opportunity score. | Proxy only. | PASS. |
| `hee_founder_speciation_opportunity()` | Opportunity emphasising colonisation/founder events. | Colonisation/accessibility/source context. | Opportunity score. | Proxy only. | PASS. |
| `hee_insitu_speciation_opportunity()` | Opportunity emphasising persistent in-situ arena. | Persistence/opportunity/context table. | Opportunity score. | Proxy only. | PASS. |
| `hee_radiation_opportunity()` | Opportunity for clade radiation. | Area/opportunity/lineage context. | Opportunity score. | Proxy only. | PASS. |
| `hee_update_species_pool()` | Updates available regional species pools from lineage events and times. | Initial pool, lineage events, `time_ma`. | Pool membership/status. | Uses explicit events; does not infer true speciation. | PASS. |

### 3.5 Static projection, limitation, and reliability

| Function | Purpose | Required structure | Key return columns | Boundary behavior | Status |
|---|---|---|---|---|---|
| `hee_combine_static()` | Backward-compatible static projection combiner. | Suitability plus optional accessibility, phylo, geo, dispersal, extrapolation weights. | `probability`, component columns. | Only named process columns are used. | PASS. |
| `hee_combine()` | Main M1-M5 style multiplicative combiner. | Suitability and optional process tables with explicit join keys. | `S_HMSC`, `E_phylo`, `A_BGB`, `L_land`, `D_static`, `D_dynamic`, `Q_extrap`, `probability`. | Join row-count changes error. | PASS. |
| `hee_mechanism_attribution()` | Identifies the lowest or most limiting component per row. | Projection components. | `dominant_mechanism`. | Diagnostic, not causal proof. | PASS. |
| `hee_limitation_map()` | Aggregates limitation labels to map/community level. | Mechanism attribution table. | `dominant_community_mechanism`, `mechanism_fraction`. | Handles weights; does not infer causation. | PASS. |
| `hee_prediction_reliability()` | Combines epistemic reliability components. | Scalar, vector, or data-frame diagnostics. | `prediction_reliability` and component scores. | Missing optional components are neutral or weights re-normalise by exponentiation. | PASS. |
| `hee_geoprocess_diagnostics()` | Convenience wrapper for full geoprocess diagnostics. | Landscape, connectivity, suitability, masks, optional diagnostics. | List of process tables. | Uses the same safe helpers. | PASS. |

## 2026-07-27 Addendum: Remaining Geological-Process Gaps Closed

The follow-up audit against `地学过程-20260727.docx` found that several process
concepts were implemented as lower-level tables but not yet exposed as named
framework functions, scenario outputs, figures, or report-ready tables. The
following additions close those workflow gaps.

| Function | Formula or rule | Output | Interpretation boundary | Status |
|---|---|---|---|---|
| `hee_network_metrics()` | `fragmentation = 1 - largest_component / n_nodes`; corridor persistence is the fraction of time slices with an open edge. | Node centrality, component count, edge density, largest-component fraction, fragmentation, corridor persistence. | Dynamic graph proxy, not observed dispersal. | PASS |
| `hee_extinction_layers()` | Local risk from supplied `epsilon`; regional proxy `1 - mean occupancy`; lineage proxy `I(sum occupancy <= 0)`. | Local, regional, and lineage extinction proxy tables. | Observed extinction requires independent fossil or lineage-extinction evidence. | PASS |
| `hee_cradle_museum_grave()` | Dominant role is the largest of cradle, museum, grave, source, and sink proxy scores. | Regional role scores and label. | Descriptive proxy; not causal classification. | PASS |
| `hee_dynamic_model_family()` | Defines M0-M8 from climate-only to feedback-enabled Earth-Biota assembly. | Model id, formula, required components, implementation status. | Design map; does not automatically fit all M0-M8 models. | PASS |
| `hee_geoprocess_hypotheses()` | Defines H1-H8 mechanism hypotheses and proxy outputs. | Hypothesis, prediction, proxy outputs, interpretation limits. | Hypothesis map; not evidence by itself. | PASS |
| `hee_geoprocess_model_comparison()` | Descriptive summaries across named model outputs; flags `geographic_existence == 0 & probability > 0`. | Mean probability, impossible fraction, richness/validation summaries when supplied. | Diagnostic comparison unless paired with formal validation/likelihood. | PASS |

These functions are now integrated into:

- `hee_geoprocess_scenario_suite()` result tables.
- `plot_geoprocess_scenario_suite()` plot list.
- `scripts/case01_simulated_full_540Ma.R` full Case 1 output tables, figures,
  evidence index, and Word report.
- `scripts/case02_realdata_full_540Ma_template.R` real-data template outputs and
  report placeholders.
- `README.md` geological-process overview table.
- `tests/testthat/test-geoprocess-framework.R` and
  `tests/testthat/test-geoprocess-scenarios.R`.

Remaining scientific limitations are explicit: true plate reconstruction,
BioGeoBEARS fitting, fossil-calibrated extinction, and actual speciation-rate
estimation still require external data or specialist models. The package
computes transparent proxies and state transitions, not causal proof.

Formula notes:

- Static/five-model projection:
  - M1: `P = S_HMSC * L_land`.
  - M2: `P = S_HMSC * E_phylo * L_land`.
  - M3: `P = S_HMSC * E_phylo * A_BGB * L_land`.
  - M4: `P = S_HMSC * E_phylo * A_BGB * D_static * L_land`.
  - M5: `P = S_HMSC * E_phylo * A_BGB * D_dynamic * L_land`.
- Prediction reliability:
  `PRI = (1 - extrapolation_risk)^wE * model_agreement^wM * phylo_validity^wP * calibration^wC`.
- Reliability is an epistemic diagnostic, not a posterior credible interval.

### 3.6 Deep-time HMSC, BioGeoBEARS, plate correction, and validation interfaces

| Function group | Main functions | Directly computed by package | External requirement | Status |
|---|---|---|---|---|
| Time cube and recipe | `hee_load_timecube()`, `hee_index_timecube()`, `hee_time_axis_from_env_cube()`, `hee_lock_recipe()`, `hee_apply_recipe()` | Yes | Paleoenvironment cube. | PASS. |
| Leakage/time mask | `hee_detect_leakage()`, `hee_phylo_time_mask()` | Yes | Species origins or dated tree. | PASS. |
| HMSC projection | `hee_project_hmsc_timecube()`, `hee_project_hmsc_table()` | Yes for matrix/posterior prediction adapter; real HMSC fitting is through `Hmsc`. | Fitted HMSC/evo object or Beta-like parameters. | PASS. |
| BioGeoBEARS | `hee_bgb_prepare()`, `hee_bgb_run_models()`, `hee_bgb_import()`, `hee_bgb_compare()`, `hee_bgb_lrt()`, `hee_bgb_bsm()`, `hee_bgb_accessibility()` | Import/runner/model-averaging helpers; not an internal BioGeoBEARS engine. | BioGeoBEARS or precomputed accessibility/events. | PASS with external dependency boundary. |
| Plate correction | `hee_import_plate_points()`, `hee_reconstruct_points()`, `hee_extract_paleoenv_track()`, `hee_compare_plate_corrected()` | User-table and comparison helpers. | External plate reconstruction table/backend for real paleo coordinates. | PASS with external dependency boundary. |
| Dynamic filters | `hee_static_dispersal()`, `hee_static_region_dispersal()`, `hee_dynamic_region_filter()`, `hee_dynamic_filter()` | Yes, from supplied distances/connectivity/projections. | Calibrated dispersal scale if used for inference. | PASS. |
| Community summaries | `hee_richness()`, `hee_turnover()`, `hee_refugia()`, `hee_cwm()`, `hee_phylo_diversity()` | Yes. | Traits/tree for CWM/PD. | PASS. |
| Validation | `hee_validate_fossils()`, `hee_cv_random()`, `hee_cv_spatial_block()`, `hee_cv_environment_block()`, `hee_simulate_recovery()` | Yes for summaries and simulation recovery; real fossil validation needs fossil records. | Independent fossils/pseudo-hindcast design for real inference. | PASS. |

## 4. High-risk checks and final status

| Risk item | Status | Evidence |
|---|---|---|
| Duplicate function definitions or old functions overwriting new ones | PASS | Automated scan found 0 duplicate `hee_*` definitions. |
| Internal helper functions accidentally exported | PASS | 0 `.hee_*` exports. Three unexported public-prefix helpers were renamed to `.hee_*`. |
| Exported functions missing from `NAMESPACE` or `.Rd` aliases | PASS | 0 missing exports, 0 export-without-definition, 0 missing Rd aliases. |
| Time direction mismatch | PASS | Dynamic functions validate `time_ma`, reject negative/repeated invalid times, and sort old-to-young internally where propagation occurs. |
| Negative or non-finite `delta_t` | PASS | `.hee_validate_delta_t()` errors on negative or non-finite intervals. |
| Probability outside `[0, 1]` | PASS | Core functions use `.hee_clip01()` and tests cover extreme, zero, one, and non-finite inputs. |
| Neutral defaults changing biological process strength | PASS | Missing colonisation/extinction tables in dynamic assembly are neutral; omitted, explicit zero, or unmatched source pressure is zero propagule evidence. |
| Uncontrolled Cartesian joins | PASS | Safe join helpers check keys and row counts; species expansion in functional connectivity must be explicit. |
| Duplicate join keys silently retained | PASS | Duplicate species, cell-time, region-time, and metric keys error with offending key values. |
| Single positive value collapsed to zero | PASS | Opportunity/reliability transforms avoid unsafe min-max behavior for one-row and constant positive inputs. |
| Geographic non-existence ignored | PASS | Geographic existence zero forces extinction/projection probability to zero. |
| Static combination silently ignoring `D_static` or `weight` | PASS | Static combiner recognises `geographic_existence`, `D_static`, `weight`, `accessibility`, `structural_connectivity`, `functional_connectivity`, and `climatic_connectivity`. |
| Reliability unable to handle scalar/vector/data frame | PASS | Tests cover numeric scalar, vector, and table inputs. |
| Case 2 fake real data | PASS | Missing files stop with explicit `Missing required file` error; no mock real-data fallback. |
| Full-axis support | PASS | Cube has 109 time slices from 540 to 0 Ma; full-mode assertion passes. |

## 5. Scientific interpretation boundaries

Implemented as direct calculations:

- Matrix and table validation, recipe locking/application, time-axis extraction,
  phylogenetic time mask, HMSC projection adapters, static/dynamic projection
  combination, richness, turnover, refugia, CWM, simple PD summaries, output
  cataloguing, and report generation.

Implemented as transparent proxy or heuristic indices:

- Ecological opportunity, isolation history, colonisation probability,
  extinction probability, speciation opportunity, mechanism attribution,
  limitation maps, and prediction reliability.

Implemented as external-data interfaces:

- BioGeoBEARS model fitting and BSM, model-averaged historical accessibility,
  plate reconstructions, palaeodistance/barrier reconstructions, fossil
  validation, and real dispersal-rate calibration.

Not causal without additional design:

- Historical accessibility, residual association, colonisation pressure,
  rescue effect, opportunity indices, limitation maps, and prediction
  reliability should not be interpreted as causal evidence.

Parameters needing calibration for real inference:

- Distance-decay coefficients (`alpha`, `phi`), barrier weights, functional
  connectivity coefficients (`kappa`), colonisation/extinction coefficients,
  dispersal scales, fossil-detection thresholds, extrapolation thresholds, and
  reliability weights.

Functions most suitable for real data:

- Data checks and alignment, recipe functions, phylogenetic time mask,
  projection combination, richness/turnover/refugia/CWM/PD summaries, and
  reporting functions, provided the user supplies real HMSC fits, real
  paleoclimate layers, real plate reconstructions, and real historical
  accessibility.

Functions mainly for simulation and sensitivity analysis unless calibrated:

- Colonisation/extinction transition probabilities, speciation-opportunity
  indices, dynamic assembly, prediction reliability, and mock BioGeoBEARS
  objects in Case 1.

## 6. Documentation audit

| Requirement | Status |
|---|---|
| Exported functions have Rd aliases | PASS. |
| Roxygen usage matches current implementation | PASS: `R CMD check` code/documentation mismatch check is OK. |
| Time direction documented for time-axis and dynamic functions | PASS. |
| Probability bounds documented for process helpers | PASS. |
| Neutral defaults documented | PASS in README and roxygen. |
| Real-data template requirements documented | PASS in scripts and README boundaries. |
| Process proxies vs external reconstructions vs simulation outputs distinguished | PASS in README and this audit. |

## 7. Test audit

Total `test_that()` blocks: 88.  
Geological-process focused `test_that()` blocks: 36.

Covered cases include:

- Empty, zero-row, one-row, mixed-type, NA/non-finite, repeated-key, unordered,
  and extra-column inputs.
- Time axes `540, 500, 0`, reversed input order, unequal intervals, single time
  slices, repeated times, illegal negative ages, and `delta_t`.
- Probability bounds for colonisation, extinction, rescue, structural,
  functional, climatic connectivity, speciation opportunity, ecological
  opportunity, occupancy, and prediction reliability.
- Directional monotonicity: suitability/source/accessibility do not lower
  colonisation, barriers and distance do not increase connectivity, dispersal
  ability does not lower functional connectivity, climate mismatch does not
  increase climatic connectivity, rescue does not increase extinction, and
  geographic absence forces final probability to zero.
- Join behavior: extra species, species order changes, duplicate species,
  duplicate edge keys, species-expanded vs species-specific connectivity, and
  row-count preservation.
- Cross-function consistency from structural connectivity to functional and
  climatic connectivity, colonisation, extinction, dynamic assembly, static
  combination, and reliability.

## 8. Fixes made in this final acceptance pass

| File | Change | Reason | API impact |
|---|---|---|---|
| `R/write_outputs.R` | Renamed unexported helpers `hee_sanitize_ordered_name()`, `hee_order_key()`, and `hee_filename_order_number()` to `.hee_*`. | Avoid public-prefix internal helpers being mistaken for missing exports. | No exported API change. |
| `README.md` | Added a geological-process module overview table. | Documents direct calculations, external-data needs, proxy status, and limitations. | Documentation only. |
| `NEWS.md` | Changed development heading to `HmscEcoEvo 0.2.0`. | Removes R CMD check NOTE about unparseable news entries. | Documentation only. |
| `.Rbuildignore` | Added `GEOPROCESS_AUDIT.md`. | Keeps the local audit deliverable out of source-package checks. | No API change. |

## 9. Remaining limitations and next development priorities

Remaining limitations:

- The package does not infer plate positions from present coordinates alone.
  Real plate correction requires user-supplied paleo coordinates or an external
  backend such as `rgplates`/`pygplates`.
- BioGeoBEARS fitting/BSM is not implemented as a full internal engine. Real
  analyses should import fitted BioGeoBEARS outputs or provide a runner.
- Colonisation, extinction, opportunity, and reliability formulas are
  configurable scenario models unless calibrated to independent data.
- Quick workflow output is for development validation; formal inference should
  use full mode and adequate HMSC MCMC settings.
- Full Case 1 output currently uses the available 5 Myr environmental cube
  resolution: 109 time slices from 540 to 0 Ma, not every 1 Myr.

Next priorities:

1. Add optional wrappers around real `rgplates`/`pygplates` outputs with strict
   CRS and plate-model metadata.
2. Add a BioGeoBEARS import validator that checks model names, log-likelihoods,
   parameters, model weights, and accessibility array dimensions against the
   tree and regions.
3. Add calibration workflows for colonisation/extinction coefficients using
   independent occurrence or fossil data.
4. Add richer Word-report diagnostics for full HMSC cross-validation that refits
   models by random, spatial, environmental, and clade folds.
5. Add visual snapshot tests or image hash checks for representative figures.

## 10. Addendum: Remaining formula, mask, and connectivity audit (2026-08-11)

This follow-up audit rechecked the most common way deep-time projections can be
over-optimistic: treating incomplete supplied process evidence as full
accessibility, full passability, or full lineage/geographic existence. The
package now uses one consistent rule: omitted optional components are neutral,
but supplied incomplete components are conservative.

Critical and major fixes made:

- `hee_dynamic_assembly()` now accepts species-by-time matrix `phylo_mask`
  inputs, recognises `phylo_existence`, and treats unmatched entries as lineage
  absence.
- Source-pressure-only dynamic assembly now calls
  `hee_colonisation_probability()` with `interval_myr` instead of using
  `suitability * source_pressure` directly.
- `hee_colonisation_probability()` now treats unmatched or missing
  accessibility rows as zero accessibility when an accessibility table is
  supplied.
- `hee_update_occupancy()`, `hee_combine_static()`, and `hee_combine()` now use
  conservative values for supplied but missing geographic/phylogenetic masks.
- Structural, functional, static, dynamic, and cube connectivity helpers now
  treat explicit missing barrier, route, or passability information as
  conservative instead of fully passable.
- Dispersal scales are now checked for finite non-negative values. Scale zero
  is allowed and means no movement across positive distance.
- `hee_make_dispersal_matrices()` now handles `scale = 0` without `NaN` and
  rejects finite negative distances.

Regression tests added or extended:

- `tests/testthat/test-geological-processes.R`
- `tests/testthat/test-geological-processes-systematic.R`
- `tests/testthat/test-geoprocess-scenarios.R`
- `tests/testthat/test-projection-components.R`
- `tests/testthat/test-combine.R`

The detailed machine-readable table is saved as
`outputs/audit/geoprocess_remaining_fix_20260811.csv`; the narrative audit note
is saved as `outputs/audit/geoprocess_remaining_fix_20260811.md`.
