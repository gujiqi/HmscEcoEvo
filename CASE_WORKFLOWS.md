# Case workflow status (HmscEcoEvo 1.0.4)

The current conceptual model has six biological processes: environmental
filtering, dispersal, biotic filtering, evolution, speciation and extinction.
The cases below do not all estimate every process. Dynamic Earth is an external
driver; maps of richness, patches, refugia and corridors are derived results.

| Case | Current entry point | Role and interpretation |
| --- | --- | --- |
| 01 | `scripts/case01_simulated_full_540Ma.R` | Simulated integration and legacy M1-M5 sensitivity only; generated history is not empirical evidence. |
| 02 | `scripts/case02_realdata_full_540Ma_template.R` | Strict real-input audit/template; it stops when required inputs are absent. |
| 03 | `scripts/case03_geoprocess_realdata_past_to_present.R` | Optional nested regional-history and historical-driver workflow; uncalibrated transition summaries are scenario diagnostics. |
| 04 | `scripts/case04_plant200_ancestral_suitability_no_bgb.R` | HMSC tip Beta posterior, dated-tree ancestral response and time-matched palaeoenvironmental support; no historical occupancy is inferred. |
| 05 | `scripts/case05_v8_terrestrial_three_scheme_4deg.R` and `scripts/case05_v8_compare_four_schemes.R` | Global terrestrial lineage-location inference and dispersal-scheme comparison; plate carriage is separate from active movement. No demographic establishment or persistence rates are fitted. |
| 06 | `scripts/case06_plant200_response_model_study.R` | Matched BM, OU and EB ancestral-response comparison; not a palaeo-occupancy reconstruction. |
| 07 | `scripts/case07_patch_evidence_ladder.R`, `scripts/case07_real_time_map_atlas.R` and `scripts/case07_plant200_source_corridor_maps.R` | Real environmental patch/candidate-refugium summaries and explicitly labelled forward-predicted movement-opportunity maps. Controlled joint-edge examples are demonstrations, not empirical historical paths. |

All current entry points require the installed six-process package version
1.0.4 or load that source tree. The older `case04_plant200_global_dynamic_*`
and `case05_plant200_eight_process_three_dispersal.R` scripts remain as
reproducible, scenario-calibrated occupancy experiments. Their colonisation
and persistence terms are **legacy state-transition components**, not seventh
and eighth HmscEE processes or outputs of the current Case05.

The package test suite parses every `case01` through `case07` R script and
checks the current Case04-Case07 entry points for obsolete core helper calls.
Full real-data runs require their original external data and may take hours;
parsing and unit tests do not substitute for rerunning those analyses.

## Verification on 2026-09-27

The installed package was rebuilt as version 1.0.4 for R 4.5.3 and R 4.4.3.
The source test suite passed without failures. A clean-source `R CMD check`
passed with one expected NOTE because the check used a source directory rather
than a package tarball; tests, examples and vignettes were excluded from that
check and validated separately where noted.

| Case | Run performed | Result |
| --- | --- | --- |
| 01 | Quick simulation, 10 times | Completed; PNG figures and a Word report were generated. Duplicate tutorial/report renderings are not counted as independent results. |
| 02 | Script parsing and strict input audit | Real run not possible: the required complete input bundle, including `plate_points.csv`, is absent. |
| 03 | Quick nested demo | Completed; a real-data fit still requires its external inputs. |
| 04 | Reduced ancestral-suitability smoke run | Completed for 135 and 0 Ma, five tips and one response draw. Not a full atlas. |
| 05 | Reduced current-v8 terrestrial smoke run | Completed for 66 slices and four tips; all 13 quality gates passed. Not the full 200-tip four-scheme rerun. |
| 06 | Full response-model comparison | Completed for 200 tips, 66 times, five draws and BM/OU/EB; all eight gates passed. |
| 07 | Patch analysis and real-time map atlases | Completed; 28 patch gates passed, 594 time-map PNG/GeoTIFF pairs and 132 forward-predictive source/corridor PNG/GeoTIFF pairs. The latter are not posterior realized flows. |

Fresh verification outputs are under `HmscEcoEvo/outputs/six_process_case*`.
These runs demonstrate software execution, not that uncalibrated deep-time
rates or missing historical observations have become empirically identified.
