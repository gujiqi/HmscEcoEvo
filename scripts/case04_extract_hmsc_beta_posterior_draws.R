#!/usr/bin/env Rscript

# Extract genuine retained HMSC Beta samples for Case04/Case05.  This script
# intentionally works from the fitted modern HMSC object, never from modern
# prediction maps or posterior means.

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x)) y else x

parse_args <- function(x) {
  out <- list()
  for (arg in x) {
    if (!grepl("^--", arg)) next
    kv <- sub("^--", "", arg)
    pos <- regexpr("=", kv, fixed = TRUE)
    if (pos < 0L) out[[kv]] <- TRUE else {
      out[[substr(kv, 1L, pos - 1L)]] <- substr(kv, pos + 1L, nchar(kv))
    }
  }
  out
}

as_int <- function(x, default) {
  z <- suppressWarnings(as.integer(x)[1L])
  if (is.na(z)) default else z
}

safe_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

write_csv <- function(x, path) {
  safe_dir(dirname(path))
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(x, path)
  } else {
    utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
  }
  invisible(path)
}

cfg <- parse_args(args)
pkg_root <- normalizePath(cfg$pkg_root %||%
                            "C:/Users/Google/Documents/HMSC-HIST/HmscEcoEvo",
                          winslash = "/", mustWork = FALSE)
previous_case04 <- normalizePath(
  cfg$previous_case04 %||%
    "C:/Users/Google/Documents/HMSC-HIST/HmscEcoEvo_API8_20260818/outputs/case04_standard200_true_bgb_all109_stage1_20260821",
  winslash = "/", mustWork = FALSE
)
model_path <- normalizePath(
  cfg$model %||% file.path(previous_case04, "03_hmsc_environment", "case04_hmsc_model.rds"),
  winslash = "/", mustWork = FALSE
)
recipe_path <- normalizePath(
  cfg$recipe %||% file.path(previous_case04, "03_hmsc_environment", "case04_modern_0Ma_environment_recipe.rds"),
  winslash = "/", mustWork = FALSE
)
output <- safe_dir(cfg$output %||%
                     file.path(pkg_root, "derived_inputs", "case04_plant200_hmsc_beta_posterior"))
n_draws <- as_int(cfg$n_draws, 25L)
seed <- as_int(cfg$seed, 20260919L)

if (!requireNamespace("HmscEcoEvo", quietly = TRUE)) {
  stop("HmscEcoEvo must be installed before extracting HMSC posterior draws.", call. = FALSE)
}
if (!file.exists(model_path)) stop("Missing HMSC model: ", model_path, call. = FALSE)
if (!file.exists(recipe_path)) stop("Missing locked modern recipe: ", recipe_path, call. = FALSE)

model <- readRDS(model_path)
recipe <- readRDS(recipe_path)
basis_cols <- setdiff(as.character(recipe$variables), "log_sampling_effort")
if (!all(basis_cols %in% model[["covNames"]])) {
  stop("The locked recipe contains variables missing from the HMSC model: ",
       paste(setdiff(basis_cols, model[["covNames"]]), collapse = ", "), call. = FALSE)
}

draws <- HmscEcoEvo::hee_hmsc_beta_posterior_draws(
  model = model,
  n_draws = n_draws,
  axes = basis_cols,
  seed = seed,
  include_intercept = FALSE
)
draw_file <- file.path(output,
                       paste0("case04_hmsc_beta_posterior_", n_draws, "draws.csv"))
write_csv(draws, draw_file)

manifest <- unique(draws[, c("response_draw", "chain", "iteration", "response_weight")])
manifest <- manifest[order(manifest$chain, manifest$iteration), , drop = FALSE]
write_csv(manifest, file.path(output, "case04_hmsc_beta_posterior_draw_manifest.csv"))

summary <- do.call(rbind, lapply(basis_cols, function(axis) {
  x <- draws[[axis]]
  data.frame(
    axis = axis,
    posterior_mean = mean(x, na.rm = TRUE),
    posterior_sd = stats::sd(x, na.rm = TRUE),
    q025 = stats::quantile(x, 0.025, na.rm = TRUE, names = FALSE),
    q50 = stats::quantile(x, 0.5, na.rm = TRUE, names = FALSE),
    q975 = stats::quantile(x, 0.975, na.rm = TRUE, names = FALSE),
    stringsAsFactors = FALSE
  )
}))
write_csv(summary, file.path(output, "case04_hmsc_beta_posterior_axis_summary.csv"))

metadata <- data.frame(
  field = c("model", "recipe", "model_md5", "recipe_md5", "n_available_posterior_samples",
            "n_selected_posterior_samples", "sampling_design", "n_species",
            "n_environment_axes", "include_hmsc_intercept_in_projection", "seed",
            "HmscEcoEvo_version", "Hmsc_version", "scientific_boundary"),
  value = c(
    model_path, recipe_path,
    unname(tools::md5sum(model_path)), unname(tools::md5sum(recipe_path)),
    attr(draws, "n_available_posterior_samples"),
    attr(draws, "n_selected_posterior_samples"),
    attr(draws, "sampling_design"), length(unique(draws$lineage)),
    length(basis_cols), FALSE, seed,
    as.character(utils::packageVersion("HmscEcoEvo")),
    if (requireNamespace("Hmsc", quietly = TRUE)) as.character(utils::packageVersion("Hmsc")) else NA_character_,
    "actual_retained_HMSC_Beta_samples;_intercept_retained_for_audit_but_not_projected_as_ancestral_niche"
  ),
  stringsAsFactors = FALSE
)
write_csv(metadata, file.path(output, "case04_hmsc_beta_posterior_extraction_metadata.csv"))

writeLines(c(
  "# Case04 actual HMSC Beta posterior extraction",
  "",
  "Each response_draw is a retained MCMC Beta sample from the fitted Plant200 HMSC object, selected with stratification across retained chains.",
  "The HMSC intercept is retained as hmsc_intercept_original for audit, but the projection intercept is set to zero by design: modern prevalence and sampling are not silently interpreted as ancestral niche responses.",
  "Use this CSV as --response_history in the global dynamic engine. The engine reconstructs each draw's tip Beta along the dated tree before projecting it through each palaeoenvironmental time slice."
), file.path(output, "README.md"), useBytes = TRUE)

message("Case04 HMSC posterior Beta extraction complete: ", draw_file)
