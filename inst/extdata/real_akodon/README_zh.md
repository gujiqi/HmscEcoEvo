# Real Akodon dataset for hmscHist

This folder contains a real-data example derived from the public Herodotools Akodon example data.

## Source

The source files were downloaded from the Herodotools GitHub repository:

- `Table_Akodon_coords_pa.txt`: Akodon occurrence/presence-absence grid data with longitude and latitude.
- `akodon.new`: Akodon phylogenetic tree in Newick format.
- `geo_area_akodon.data`: BioGeoBEARS-style species-area geography file for Akodon species.
- `size_akodon.txt`: body size data for Akodon species.

See Herodotools documentation and vignette for the original example workflow.

## Files for hmscHist

- `comm.csv`: site x species presence/absence matrix, subset of real Akodon grid cells.
- `coords.csv`: site coordinates.
- `env.csv`: WorldClim 2.1 bioclimatic variables extracted from the site coordinates. Columns are `bio1`, `bio2`, `bio3`, `bio9`, `bio10`, `bio11`, `bio15`, and `bio19`.
- `site_region.csv`: site-level region assigned from dominant Akodon species-area membership at the site.
- `site_region_scores.csv`: scores used to assign site regions.
- `traits.csv`: real Akodon body-size trait data.
- `phylo.tre`: Akodon phylogeny.
- `tip_ranges.csv`: species x biogeographic area ranges from `geo_area_akodon.data`.
- `species_region_history.csv`: hmscHist-format species-region history table derived from current area ranges and terminal branch lengths.
- `trait_history.csv`: hmscHist-format trait-history placeholder using real body size as `dispersal_trait`.

## Important caveat

The occurrence matrix, tree, body-size trait table, coordinates, and species-area geography file are real public data.

The environmental table is derived from local WorldClim 2.1 2.5 arc-min GeoTIFF files:
`wc2.1_2.5m_bio_1.tif`, `wc2.1_2.5m_bio_2.tif`, `wc2.1_2.5m_bio_3.tif`,
`wc2.1_2.5m_bio_9.tif`, `wc2.1_2.5m_bio_10.tif`, `wc2.1_2.5m_bio_11.tif`,
`wc2.1_2.5m_bio_15.tif`, and `wc2.1_2.5m_bio_19.tif`. Thirteen coastal
sites fell on NA raster cells and were filled from the nearest valid land cell
within 250 km so the example can run without missing environmental covariates.

The file `species_region_history.csv` is a standardized hmscHist demonstration table derived from real species-area ranges and the real tree. Its `entry_time` column is a terminal-branch-length proxy, not a published BioGeoBEARS/BSM colonization-time estimate. For formal historical inference, replace `species_region_history.csv` with a table extracted from your own BioGeoBEARS/BSM analysis.

No `speciation_events.csv` or `history_events.csv` are included in this no-event real-data example, because you requested a real example without event tables.
