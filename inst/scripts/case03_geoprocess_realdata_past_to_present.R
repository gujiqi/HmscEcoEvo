#!/usr/bin/env Rscript

# Case 03 real-data workflow: HmscEE nested macro-region / micro-cell model.
# Real mode never simulates missing biological, geographical, BSM, or plate data.

args_full <- commandArgs(FALSE)
cmd_file <- args_full[grep("^--file=", args_full)]
script_path <- if (length(cmd_file)) {
  normalizePath(sub("^--file=", "", cmd_file[[1]]), winslash = "/",
                mustWork = FALSE)
} else {
  normalizePath("scripts/case03_geoprocess_realdata_past_to_present.R",
                winslash = "/", mustWork = FALSE)
}
pkg_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/",
                          mustWork = FALSE)
load_hmscecoevo <- function(script_path) {
  if (requireNamespace("HmscEcoEvo", quietly = TRUE)) {
    library(HmscEcoEvo)
    return(invisible(TRUE))
  }
  roots <- unique(normalizePath(c(
    file.path(dirname(script_path), ".."),
    file.path(dirname(script_path), "..", "..")
  ), winslash = "/", mustWork = FALSE))
  roots <- roots[file.exists(file.path(roots, "DESCRIPTION")) &
                   dir.exists(file.path(roots, "R"))]
  if (length(roots) && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(roots[[1]], quiet = TRUE)
    return(invisible(TRUE))
  }
  stop("HmscEcoEvo is not installed and a source package root could not be found.",
       call. = FALSE)
}
load_hmscecoevo(script_path)

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = "") {
  hit <- grep(paste0("^(--)?", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^(--)?", name, "="), "", hit[[1]])
}

analysis_mode <- tolower(arg_value("analysis_mode", "real"))
if (!analysis_mode %in% c("real", "demo")) {
  stop("--analysis_mode must be `real` or `demo`.", call. = FALSE)
}
if (analysis_mode == "demo") {
  demo <- file.path(dirname(script_path), "case03_geoprocess_simulated_demo.R")
  if (!file.exists(demo)) stop("Missing demo script: ", demo, call. = FALSE)
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c(normalizePath(demo, winslash = "/", mustWork = TRUE),
                      "--analysis_mode=demo",
                      args[!grepl("^(--)?analysis_mode=", args)]))
  quit(save = "no", status = status)
}

output <- normalizePath(arg_value("output",
                                  "outputs/case03_hmscee_nested_realdata"),
                        winslash = "/", mustWork = FALSE)
dirs <- file.path(output, c("tables", "metadata", "rds", "figures",
                            "reports", "tasks"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
tab_dir <- file.path(output, "tables")
meta_dir <- file.path(output, "metadata")
writeLines("real", file.path(meta_dir, "analysis_mode.txt"))
utils::write.csv(hee_hmscee_process_catalog(),
                 file.path(tab_dir, "00.0_hmscee_five_process_catalog.csv"),
                 row.names = FALSE)

paths <- list(
  comm = arg_value("comm"),
  sites = arg_value("sites"),
  traits = arg_value("traits"),
  tree = arg_value("tree"),
  tip_ranges = arg_value("tip_ranges"),
  env_rds = arg_value("env_rds"),
  region_cube = arg_value("region_cube"),
  landmask_cube = arg_value("landmask_cube"),
  plate_points = arg_value("plate_points"),
  bgb_accessibility = arg_value("bgb_accessibility"),
  bsm_events = arg_value("bsm_events"),
  fossils = arg_value("fossils")
)

region_history_path <- if (nzchar(paths$bsm_events)) paths$bsm_events else
  paths$bgb_accessibility
required <- data.frame(
  input = c("comm", "sites", "env_rds", "region_cube", "tree",
            "bgb_or_bsm_region_history"),
  argument = c("--comm", "--sites", "--env_rds", "--region_cube", "--tree",
               "--bgb_accessibility or --bsm_events"),
  path = c(paths$comm, paths$sites, paths$env_rds, paths$region_cube,
           paths$tree, region_history_path),
  reason = c(
    "modern site x species matrix for P1 HMSC abiotic filtering",
    "modern survey coordinates for 0 Ma XData extraction and validation",
    "external driver layer for P1 X, hard habitat H and P2 movement resistance W",
    "time-specific region map for P2/P5 BioGeoBEARS region states r(c,k)",
    "P4/P5 dated phylogeny for active lineages and ancestral response reconstruction",
    "P2/P5 BioGeoBEARS/BSM regional history R[lineage,region,time]"
  ),
  stringsAsFactors = FALSE
)
required$status <- ifelse(nzchar(required$path) & file.exists(required$path),
                          "AVAILABLE", "MISSING")
utils::write.csv(required, file.path(tab_dir, "00.1_required_real_inputs.csv"),
                 row.names = FALSE)
if (any(required$status == "MISSING")) {
  missing <- required[required$status == "MISSING", , drop = FALSE]
  utils::write.csv(missing, file.path(output, "missing_required_inputs.csv"),
                   row.names = FALSE)
  stop("Missing required Case 03 real inputs. See ",
       file.path(output, "missing_required_inputs.csv"), call. = FALSE)
}

comm <- utils::read.csv(paths$comm, check.names = FALSE,
                        stringsAsFactors = FALSE)
sites <- utils::read.csv(paths$sites, check.names = FALSE,
                         stringsAsFactors = FALSE)
region_cube <- utils::read.csv(paths$region_cube, stringsAsFactors = FALSE)
region_history_raw <- utils::read.csv(region_history_path,
                                      stringsAsFactors = FALSE)

if (!"site_id" %in% names(comm)) stop("comm must contain site_id.", call. = FALSE)
if (!all(c("site_id", "lon", "lat") %in% names(sites))) {
  stop("sites must contain site_id, lon and lat.", call. = FALSE)
}
if (!setequal(comm$site_id, sites$site_id)) {
  stop("comm and sites must contain exactly the same site_id values.",
       call. = FALSE)
}
species <- setdiff(names(comm), "site_id")
Y <- as.matrix(comm[, species, drop = FALSE])
storage.mode(Y) <- "numeric"

if (!all(c("cell_id", "time_ma", "region") %in% names(region_cube))) {
  stop("region_cube must contain cell_id, time_ma and region.", call. = FALSE)
}
if (!all(c("lineage", "region", "time_ma") %in% names(region_history_raw))) {
  stop("BioGeoBEARS/BSM history must contain lineage, region and time_ma.",
       call. = FALSE)
}

region_history <- hee_bsm_region_history(region_history_raw)
utils::write.csv(region_history, file.path(tab_dir,
                                           "01.1_bsm_region_history_R.csv"),
                 row.names = FALSE)

metadata <- data.frame(
  item = c("analysis_mode", "n_sites", "n_species", "n_region_history_rows",
           "core_framework", "biological_process_engines",
           "external_layers_not_processes"),
  value = c("real", nrow(sites), length(species), nrow(region_history),
            "Five-process HmscEE with BioGeoBEARS/BSM regional history + within-region HmscEE occupancy",
            "P1 abiotic filtering; P2 dispersal-colonisation; P3 biotic interactions; P4 niche/trait evolution; P5 speciation-lineage extinction",
            "Climate, geology, geography, habitat, dated trees, BioGeoBEARS, fossils, extrapolation and phyloregion diagnostics"),
  stringsAsFactors = FALSE
)
utils::write.csv(metadata, file.path(meta_dir, "case03_hmscee_real_metadata.csv"),
                 row.names = FALSE)
saveRDS(list(Y = Y, sites = sites, region_history = region_history),
        file.path(output, "rds", "case03_real_verified_inputs.rds"))

message("Case 03 real input audit completed. No simulated fallback was used.")
