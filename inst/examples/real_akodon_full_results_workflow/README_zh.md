# Real Akodon full results workflow

Run this from the package root or after installing `hmscHist`:

```r
source(system.file("examples/real_akodon_full_results_workflow/00_run_all_real_akodon_results_zh.R",
                   package = "HmscEcoEvo"))
```

It uses `inst/extdata/real_akodon` and writes all outputs to:

```text
hmscHist_real_akodon_full_results/
```

Start with:

```text
hmscHist_real_akodon_full_results/results/ALL_RESULTS_INDEX.csv
```

The MCMC settings are intentionally short for a runnable example. Increase
`samples`, `transient`, `thin`, and `nChains` in `S2_fit_models_short.R` for
formal inference.
