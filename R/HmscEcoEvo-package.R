#' HmscEcoEvo: Eco-evolutionary diagnostics for HMSC outputs
#'
#' HmscEcoEvo preserves the historical-predictor workflow originally developed
#' in hmscHist and adds a post-processing framework for HMSC eco-evolutionary
#' diagnostics. It extracts Beta, Gamma, Rho, Omega, posterior draws, traits,
#' phylogeny, historical predictors, and optional population/time-series data
#' into a unified [hmsc_ecoevo()] object, then computes seven groups of metrics
#' and diagnostic plots.
#'
#' Version 1.0 also implements HmscEE as a Dynamic Earth-Biota Assembly
#' Framework for deep-time projections. The public API names eight explicit
#' biological processes: environmental filtering, dispersal, colonisation,
#' biotic filtering, persistence, evolution, speciation, and extinction.
#' Climate, geology, geography and habitat history are external drivers, not
#' additional biological processes. BioGeoBEARS/BSM and dated trees are
#' historical inference tools or constraints, not an extra process; fossils,
#' pollen, aDNA, extrapolation diagnostics and phyloregion summaries are
#' observation, uncertainty or result layers.
#'
#' In the current nested implementation, BioGeoBEARS/BSM supplies macro-scale
#' lineage-by-region history, HMSC and phylogenetic response reconstruction
#' supply environmental response, palaeo-Earth layers supply environment,
#' hard habitat existence and movement resistance, and HmscEE updates
#' conditional cell occupancy inside occupied regions. The core decomposition
#' is:
#' \deqn{P(Z_{\ell,u,k}=1) =
#' \sum_g w_g \sum_h w_{h|g}\sum_c w_{c|g}\sum_s w_s
#' R_{\ell,r(u,k),k}^{(h,g)}q_{\ell,u,k}^{(h,g,c,s)}.}
#' Here `g` is a palaeogeographic/plate scenario, `c` is a palaeoclimate
#' scenario conditional on `g`, `h` is a BioGeoBEARS/BSM history conditional on
#' `g`, `s` is an HMSC/ancestral trait/ancestral response posterior draw, and
#' `u` is a spatial grid cell.
#' HmscEE does not use `e` for Earth scenario because BioGeoBEARS already uses
#' parameter `e` for range contraction/local extinction.
#' Geography scenarios are common upstream conditions for both BioGeoBEARS/BSM
#' regional histories and HmscEE hard habitat/resistance state; BSM histories
#' are not averaged independently of the geography scenario that generated
#' them.
#' Ancestral traits, ancestral environmental response, and ancestral geographic
#' range are separate objects: BioGeoBEARS supplies range `R`, not the
#' trait-mediated response \eqn{\beta_{\ell,k} = \Gamma T_{\ell,k}+u_{\ell,k}}.
#' HMSC residual associations are not interpreted as causal biotic interactions
#' unless independent interaction evidence is supplied.
#' Extrapolation and no-analog diagnostics are reported as epistemic
#' uncertainty rather than multiplied into biological probability. The package
#' does not infer true plate motion, causal dispersal routes, realised
#' speciation rates, or true historical species distributions by itself.
#'
#' The package is not a replacement for Hmsc, BioGeoBEARS, bayou, SURFACE,
#' phylogenetic comparative packages, or causal eco-evolutionary models.
#' Model-based branch shifts, rate shifts, adaptive peaks, convergence regimes,
#' tree uncertainty, prior sensitivity, simulation recovery, and full
#' posterior-predictive checks should be supplied as precomputed results from
#' appropriate external workflows.
#'
#' @importFrom stats na.omit setNames
#' @importFrom utils modifyList
#' @keywords internal
"_PACKAGE"

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "LA_PL_ratio", "PC1", "PC2", "P_shift", "axis", "axis1", "axis2",
    "axis_specific_rate", "beta_distance", "beta_magnitude", "block",
    "branch", "branch_contrast_outlier_score", "breadth_index", "clade",
    "component", "contribution", "correlation", "delta_r2", "delta_theta",
    "disparity", "draw",
    "estimate", "estimated_value", "feedback_strength",
    "GxE", "gamma_shift_probability", "gamma_sign_flip_index",
    "gamma_variance_among_clades", "genetic_niche_signal", "genotype",
    "group1", "group2", "local_moran_i", "lwr", "metric",
    "max_P_shift", "missing_trait_risk_index", "moran_i",
    "multivariate_shift_magnitude", "n_species", "node",
    "node_signal", "number_of_niche_shifts", "omitted_trait", "optimum",
    "optimum_type", "parameter", "path", "phylo_distance", "population",
    "predictor_set", "r2", "r2_with_coords", "redundancy",
    "residual_beta", "response", "response_axis", "response_lwr",
    "response_upr", "signal", "similarity", "source_group",
    "species", "support", "time", "trait", "trait_axis_specialization",
    "true_value", "turnover", "upr", "value", "variable",
    "variance_ratio", "weight",
    "within_species_niche_divergence"
  ))
}
