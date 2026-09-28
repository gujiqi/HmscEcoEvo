# HmscEcoEvo Example Data And Real-Data Schemas

This directory is for lightweight example data and schema documentation. Large
full-axis outputs are local analysis products and are not bundled into source
package builds.

## Case 1 Simulated Full Output

The full simulated 540-0 Ma workflow writes to an `outputs/` directory, for
example:

`outputs/case01_full_final_paper_audit_20260812_run2/`

That output should contain ordered figures, tables, RDS files, metadata, task
tables, and a Word report. It is a simulated teaching and validation case.
BioGeoBEARS/accessibility, plate points, geological processes, colonisation,
extinction, speciation opportunity, and prediction reliability are mock,
proxy, or heuristic unless replaced by real external data.

## Case 3 Real And Demo Workflows

The scientific main entry point is:

`scripts/case03_geoprocess_realdata_past_to_present.R`

It defaults to `analysis_mode=real`, uses observed community/site inputs, the
complete grids in a user-supplied palaeoenvironment cube, and no simulated
fallback. Core missing files create `missing_required_inputs.csv`. Missing tree,
region, BioGeoBEARS, plate, BSM, or fossil inputs create module-specific
`NOT_RUN_*` states. Schemas are in `case03_realdata_template/`.

The teaching simulation is retained separately as:

`scripts/case03_geoprocess_simulated_demo.R`

It writes to:

`outputs/case03_geoprocess_demo_quick/`

or, for the denser run:

`outputs/case03_geoprocess_demo_full/`

Run it with:

```r
system2(file.path(R.home("bin"), "Rscript"), c(
  "scripts/case03_geoprocess_simulated_demo.R",
  "--analysis_mode=demo",
  "--quick=true",
  "--output=outputs/case03_geoprocess_demo_quick"
))
```

The demo first looks for the local 4-degree land-harmonized Phanerozoic
environment RDS and uses it when present. If that file is not available, it
falls back to a clearly labelled portable simulated environment so package
checks and teaching examples can still run. Modern `XData` are extracted from
the 0 Ma layer, a real `Hmsc::Hmsc()` model is fitted to simulated modern
community data with `Hmsc::sampleMcmc()`, and the fitted HMSC environmental
suitability is projected from 540 Ma to 0 Ma.

Species traits, species origin times, the dated teaching phylogeny and modern
community observations are simulated for the example. Historical accessibility,
plate tracks and process events are mock-only unless the user supplies external
tables via:

```r
system2(file.path(R.home("bin"), "Rscript"), c(
  "scripts/case03_geoprocess_simulated_demo.R",
  "--analysis_mode=demo",
  "--quick=false",
  "--output=outputs/case03_geoprocess_full_external_inputs",
  "--bgb_accessibility=data_raw/bgb_accessibility.csv",
  "--plate_points=data_raw/plate_points.csv",
  "--bsm_events=data_raw/bsm_events.csv"
))
```

If one of those external file paths is supplied but missing, the script stops
with a clear error instead of generating fake real-data evidence. The case is
intended to demonstrate the entire geological-process chain from 540 Ma to
0 Ma:

- ecological stage generation and disappearance;
- land age, area, heterogeneity, and ecological opportunity;
- structural, functional, and climatic connectivity;
- isolation duration, source pressure, rescue, colonisation, and extinction;
- species-pool updates and speciation opportunity proxies;
- dynamic occupancy, M1-M5 projection contrasts, richness, turnover, refugia,
  limitation maps, reliability, validation checks, ordered outputs, and a Word
  report.

The outputs are software and methods demonstrations. They are not real
palaeodistribution estimates, real dispersal routes, causal geological
mechanisms, or true speciation/extinction rates.

## Case 2 Required `data_raw/` Files

Real projects should provide these files:

| File | Required columns or contents |
|---|---|
| `data_raw/comm.csv` | `site` plus one column per species |
| `data_raw/sites.csv` | `site`, `lon`, `lat`, optional blocks or regions |
| `data_raw/traits.csv` | `species` plus trait columns |
| `data_raw/tree.tre` | dated Newick tree with branch lengths |
| `data_raw/tip_ranges.csv` | `species`, `region`, range or presence columns |
| `data_raw/env_cube_index.csv` | `time_ma`, file path, variable metadata |
| `data_raw/region_cube_index.csv` | `time_ma`, file path or layer id |
| `data_raw/landmask_cube_index.csv` | `time_ma`, file path or layer id |
| `data_raw/plate_points.csv` | `track_id`, `time_ma`, `lon`, `lat`, `plate_id`, `paleo_lon`, `paleo_lat` as applicable |
| `data_raw/connectivity_paths.csv` | `from_region`, `to_region`, `time_ma`, `least_cost_distance_km`, `barrier_strength`, `route_open`, `corridor_suitability` |
| `data_raw/fossils.csv` | `species`, `age_min`, `age_max`, `paleo_lon`, `paleo_lat`; palaeogeographic coordinates are required |

Missing required files should stop with a clear error such as:

```r
stop("Missing required file: data_raw/comm.csv")
```

The real-data template must not generate fake real results. If independent
fossil validation, BioGeoBEARS output, or plate reconstruction is missing, the
corresponding section should be skipped or explicitly labelled as absent.

## Time And Interpretation Rules

`time_ma` means Ma before present; larger values are older. `origin_ma` is used
by `hee_phylo_time_mask()` so species-level probability is zero before a
species origin (`time_ma > origin_ma`). M1-M5 projection models are scenario
products combining HMSC suitability, land, phylogenetic time, historical
accessibility, and dispersal filters. They are not observed historical ranges.

Use the numbered vignette sequence beginning with
`vignettes/v00_index_roadmap.Rmd` for the complete workflow map.
