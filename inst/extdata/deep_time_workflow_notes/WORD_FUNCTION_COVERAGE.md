# Word Function Coverage Audit

This audit was generated from the three deep-time workflow manuscripts stored as
Word documents and extracted into this directory. HmscEcoEvo now exports every
valid `hee_*` function name detected in those manuscripts.

Coverage summary:

- Valid `hee_*` function names detected from Word text: 108
- Exported `hee_*` functions in `NAMESPACE`: 111
- Missing detected Word functions: 0

Implementation classes:

- Direct computation: data checks, recipes, time-cube indexing/extraction,
  HMSC Beta projection, phylogenetic time masks, suitability combination,
  richness, turnover, refugia, fossil validation, thresholds, summaries,
  cross-validation folds, taxon tables, tiling, coordinate normalization,
  collinearity filtering, trait imputation, community-weighted means, and
  output writing.
- Wrapper/import/precomputed interface: HMSC fitting/prediction,
  BioGeoBEARS model fitting and BSM import, model-averaged accessibility,
  plate reconstruction, long MCMC diagnostics, scenario orchestration, and
  uncertainty budgets.
- Explicitly not fabricated: true BioGeoBEARS inference, GPlates-style plate
  rotation, causal eco-evolutionary feedback, and model-based adaptive
  shift/rate/peak inference.

The extra exported `hee_*` functions not literally named in the manuscripts are
small glue functions needed by the implemented workflow:

- `hee_load_timecube`
- `hee_project_hmsc_table`
- `hee_run_deep_time_pipeline`
