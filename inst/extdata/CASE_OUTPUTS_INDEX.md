# HmscEcoEvo Case Output Index

This directory contains regenerated, runnable example cases and saved workflow
outputs.

## Seven eco-evolutionary diagnostic cases

Output root:

`inst/extdata/ecoevo_cases/outputs`

Each of the seven case folders contains:

- `metrics.rds`: all calculated metric objects for the seven diagnostic modules.
- `niche_summary.pdf/png`
- `phylo_signal.pdf/png`
- `evo_transition.pdf/png`
- `trait_mediation.pdf/png`
- `gamma_evolution.pdf/png`
- `population_evolution.pdf/png`
- `validation_dashboard.pdf/png`

Summary table:

`inst/extdata/ecoevo_cases/outputs/case_run_summary.csv`

Regeneration scripts:

- `inst/extdata/ecoevo_cases/generate_ecoevo_cases.R`
- `inst/extdata/ecoevo_cases/run_all_case_plots.R`

## Deep-time hindcasting cases

Output root:

`inst/extdata/deep_time_cases/outputs`

The old small-slice cases were removed and replaced by two full-axis cases.
Both use the 4-degree Phanerozoic land-harmonized compact RDS product and the
full time axis `540, 535, ..., 0 Ma` (109 slices).

Runtime coverage check:

- Word `hee_*` functions detected: 108
- Exported `hee_*` functions: 111
- Missing detected Word functions: 0

Coverage table:

`inst/extdata/deep_time_cases/outputs/word_function_coverage_runtime_check.csv`

### case01_global540_20

Global paleo-map mode using the full 540-0 Ma time axis, 20 species, a dated
540 Ma phylogeny, unequal species-origin times, HMSC Beta projection,
phylogenetic time masking, extrapolation-risk weighting, richness, turnover,
refugia, modern validation, simulated fossil validation, and BioGeoBEARS-style
precomputed accessibility/event summaries.

Output folder:

`inst/extdata/deep_time_cases/outputs/c01`

Saved outputs include:

- `case_input_and_all_results.rds`
- `suitability_full_540_0Ma.rds`
- `projection_full_540_0Ma.rds`
- `richness_full_540_0Ma.rds`
- `turnover_full_540_0Ma.rds`
- `refugia_full_540_0Ma.rds`
- `source_and_size_summary.csv`
- `beta.csv`
- `traits.csv`
- `species_origin.csv`
- `timecube_index_full.csv`
- `summary_predictors_by_time.csv`
- `summary_extrapolation_by_time.csv`
- `summary_land_area_by_time.csv`
- `summary_richness_by_time.csv`
- `summary_turnover_by_time.csv`
- `refugia_summary.csv`
- `simulated_fossils_full_axis.csv`
- `fossil_validation.csv`
- `summary_fossil_validation.csv`
- `summary_fossil_validation_by_time.csv`
- `modern_validation_metrics.csv`
- `bgb_model_compare.csv`
- `bgb_events.csv`
- `workflow_steps.csv`
- `plot_dated_phylogeny_540Ma.pdf/png`
- `plot_beta_heatmap.pdf/png`
- `plot_source_land_area_full_axis.pdf/png`
- `plot_predictor_trajectories_full_axis.pdf/png`
- `plot_extrapolation_risk_full_axis.pdf/png`
- `plot_richness_full_axis.pdf/png`
- `plot_representative_richness_maps.pdf/png`
- `plot_turnover_full_axis.pdf/png`
- `plot_refugia_map.pdf/png`
- `plot_fossil_validation_full_axis.pdf/png`
- `plot_modern_validation.pdf/png`
- `plot_bgb_accessibility_events.pdf/png`
- `plot_workflow_steps.pdf/png`

### case02_tracks540_20

Plate-corrected track mode using the full 540-0 Ma time axis, 20 species, a
dated 540 Ma phylogeny, six moving paleo-regions/tracks, nearest-land
paleoenvironment extraction, model-averaged accessibility, dynamic dispersal
filtering, extrapolation-risk weighting, richness, turnover, refugia, modern
validation, and BioGeoBEARS-style precomputed events.

Output folder:

`inst/extdata/deep_time_cases/outputs/c02`

Saved outputs include:

- `case_input_and_all_results.rds`
- `projection_tracks_full_540_0Ma.rds`
- `track_environment_full_540_0Ma.rds`
- `richness_tracks_full_540_0Ma.rds`
- `turnover_tracks_full_540_0Ma.rds`
- `refugia_tracks_full_540_0Ma.rds`
- `source_and_size_summary.csv`
- `beta.csv`
- `traits.csv`
- `species_origin.csv`
- `plate_corrected_tracks_full_axis.csv`
- `track_environment_full_axis.csv`
- `track_environment_long.csv`
- `track_projection_full_axis.csv`
- `track_richness_full_axis.csv`
- `track_turnover_full_axis.csv`
- `track_refugia_summary.csv`
- `summary_richness_by_time.csv`
- `summary_component_weights_by_time.csv`
- `summary_accessibility_by_region_time.csv`
- `modern_validation_metrics.csv`
- `bgb_model_compare.csv`
- `bgb_events.csv`
- `workflow_steps.csv`
- `plot_dated_phylogeny_540Ma.pdf/png`
- `plot_beta_heatmap.pdf/png`
- `plot_plate_corrected_tracks.pdf/png`
- `plot_track_environment_full_axis.pdf/png`
- `plot_track_richness_full_axis.pdf/png`
- `plot_track_turnover_full_axis.pdf/png`
- `plot_projection_component_weights.pdf/png`
- `plot_accessibility_heatmap.pdf/png`
- `plot_track_refugia.pdf/png`
- `plot_modern_validation.pdf/png`
- `plot_bgb_accessibility_events.pdf/png`
- `plot_workflow_steps.pdf/png`

Summary table:

`inst/extdata/deep_time_cases/outputs/deep_time_case_run_summary.csv`

Regeneration scripts:

- `inst/extdata/deep_time_cases/generate_deep_time_cases.R`
- `inst/extdata/deep_time_cases/run_deep_time_workflows.R`
