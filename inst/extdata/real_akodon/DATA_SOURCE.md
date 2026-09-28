# Data sources and citations

This real example is derived from the public Herodotools Akodon example dataset hosted on GitHub.

Use and cite:

- Herodotools GitHub repository: https://github.com/GabrielNakamura/Herodotools
- Herodotools vignette: https://gabrielnakamura.github.io/Herodotools/articles/Intro_Herodotools_vignette.html

The Herodotools repository describes the package as a toolset for historical biogeography analysis and provides the Akodon occurrence data, phylogeny, species-area file, and trait data used here.

`env.csv` was regenerated from the site coordinates using local WorldClim 2.1
2.5 arc-min BIO GeoTIFF files for BIO1, BIO2, BIO3, BIO9, BIO10, BIO11, BIO15,
and BIO19. Direct point extraction was used where possible. Thirteen coastal
sites fell on NA raster cells and were filled from the nearest valid land cell
within 250 km.

Important: `species_region_history.csv` is a standardized demonstration table derived from real Akodon species-area ranges and phylogenetic branch lengths. Replace it with a BioGeoBEARS/BSM-derived table for formal analyses.
