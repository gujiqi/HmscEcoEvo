# Case 03 real-data input contract

These files are schemas, not observations. The real workflow never fills
missing inputs with simulated values.

Core inputs:

- `comm.csv`: `site_id` followed by one 0/1 occurrence column per species.
- `sites.csv`: `site_id`, decimal-degree `lon`, and `lat`.
- `env_rds`: compact HmscEcoEvo palaeoenvironment RDS supplied through the
  command line. Its 0 Ma layer must contain every HMSC predictor and its full
  time axis supplies the palaeo prediction grids.

Module inputs:

- `traits.csv`: one row per species. `dispersal_distance_km` is required for
  M4/M5 connectivity; other numeric traits are optional HMSC trait predictors.
- `tree.tre`: phylogenetic Newick tree. If its branch lengths are used to
  derive `origin_ma`, declare `--tree_time_unit=Ma`; substitution-scale branch
  lengths are not geological ages.
- `species_origin.csv`: optional explicit species origin ages. If absent, the
  workflow derives terminal-parent ages only from a tree explicitly declared
  to be in Ma and records that operational choice.
- `region_cube.csv` or a compatible RDS: unique `cell_id x time_ma` region
  assignments. Region identifiers are time-specific through `region_time_id`.
- `bgb_accessibility.csv`: model-derived species-region-time accessibility.
  M3-M5 are not run without this table.
- `connectivity_paths.csv`: optional externally calculated least-cost paths,
  unique by `from_region`, `to_region`, and `time_ma`. Distances are kilometres;
  `barrier_strength`, `route_open`, and `corridor_suitability` are in `[0,1]`.
  If absent, Case 03 labels its region-centroid great-circle calculation as an
  approximation rather than a realised route.
- `plate_points.csv`: externally reconstructed modern-site tracks. The package
  does not invent palaeocoordinates.
- `bsm_events.csv`: BioGeoBEARS stochastic-map events.
- `fossils.csv`: independent palaeontological validation records. Supply
  `age_min`, `age_max`, `paleo_lon`, and `paleo_lat`, or explicitly label
  `lon`/`lat` as `coordinate_frame=palaeogeographic`. Modern coordinates are
  rejected because they cannot be compared directly with ancient grids.

Example invocation:

```r
system2(file.path(R.home("bin"), "Rscript"), c(
  "scripts/case03_geoprocess_realdata_past_to_present.R",
  "--analysis_mode=real",
  "--comm=data_raw/comm.csv",
  "--sites=data_raw/sites.csv",
  "--traits=data_raw/traits.csv",
  "--tree=data_raw/tree.tre",
  "--tree_time_unit=Ma",
  "--env_rds=data_raw/environment.rds",
  "--region_cube=data_raw/region_cube.csv",
  "--bgb_accessibility=data_raw/bgb_accessibility.csv",
  "--connectivity_paths=data_raw/connectivity_paths.csv",
  "--output=outputs/case03_real"
))
```

`time_ma` is age before present in Ma; larger values are older. Coordinates in
the full environmental cube are palaeogeographic coordinates. `track_id`
identifies a reconstructed modern-point trajectory and must not be used as the
identity of the full palaeo prediction grid.
