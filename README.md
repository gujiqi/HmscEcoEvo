# HmscEcoEvo

**HmscEcoEvo** is an R package for exploring how present-day community patterns relate to species traits, phylogeny, historical environments and changing geography. It builds on fitted HMSC results and retains the historical-predictor tools from `hmscHist`.

## What can you do with it?

- **Prepare historical predictors for HMSC.** Standardize species-by-region histories, calculate historical indices, and create HMSC-ready site and trait tables.
- **Explore eco-evolutionary patterns in an HMSC fit.** Plot environmental response shapes, phylogenetic signal, evolutionary shifts, trait and history contributions, clade differences, population-level variation and model checks.
- **Estimate ancestral environmental responses.** Reconstruct response coefficients along a dated tree (BM, OU or EB) and evaluate them against the environment at each historical time slice. This produces ancestral *environmental support*, not a claim that a lineage occupied every suitable place.
- **Represent a changing Earth.** Prepare historical land and habitat masks, movement-resistance landscapes and plate-carriage links. Active biological dispersal is kept separate from movement of the land itself.
- **Compare historical movement assumptions.** Run and compare grid-based dispersal models, including landscape-aware and rare long-distance movement, with lineage locations conditioned on the dated tree and available endpoint evidence.
- **Summarize spatial results.** Make time-slice maps and derived summaries of lineage richness, diversity, community change, connected patches and candidate refugia. Source and corridor maps are labelled by the evidence used to produce them.

The functions are organized around six biological processes: environmental filtering, dispersal, biotic filtering, evolution, speciation and extinction. **A case study need not estimate all six.** In particular, the current Plant200 Case05 does not estimate historical establishment, persistence or complete lineage-extinction rates from the available data.

## Install and explore

```r
install.packages("remotes") # if needed
remotes::install_github("gujiqi/HmscEcoEvo", build_vignettes = FALSE)
library(HmscEcoEvo)

hee_function_catalog()     # find functions by task
hee_core_process_catalog() # see the six-process organization
```

Start with the [basic workflow](vignettes/v01_basic_hmsc_ecoevo_workflow.Rmd), the [historical-predictor guide](vignettes/v03_history_to_hmsc.Rmd), or the [case-study guide](CASE_WORKFLOWS.md). The [Chinese README](README_zh.md) and numbered [vignettes](vignettes/) provide more detail.

Small example inputs are included. Full empirical analyses require the fitted HMSC objects, dated trees and time-resolved Earth data described in their case scripts; the large input datasets and generated maps are not stored in this repository.
