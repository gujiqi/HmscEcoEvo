# Common helper functions for the full hmscHist + HMSC example.

set.seed(1)

pkg_file <- function(...) {
  p <- system.file(..., package = "HmscEcoEvo")
  if (nzchar(p)) return(p)
  file.path("inst", ...)
}

localDir <- file.path(getwd(), "hmscHist_HMSC_results")
dataDir <- file.path(localDir, "data")
modelDir <- file.path(localDir, "models")
resultDir <- file.path(localDir, "results")
tableDir <- file.path(resultDir, "tables")
figureDir <- file.path(resultDir, "figures")
rdsDir <- file.path(resultDir, "rds")
textDir <- file.path(resultDir, "text")

dir_create <- function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE)
}
invisible(lapply(c(localDir, dataDir, modelDir, resultDir, tableDir, figureDir, rdsDir, textDir), dir_create))

read_csv_rownames <- function(file) {
  read.csv(file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
}

write_list_as_tables <- function(x, prefix, dir = tableDir) {
  if (is.null(x)) return(invisible(NULL))
  if (is.data.frame(x) || is.matrix(x)) {
    write.csv(x, file.path(dir, paste0(prefix, ".csv")))
    return(invisible(NULL))
  }
  if (is.list(x)) {
    for (nm in names(x)) {
      obj <- x[[nm]]
      if (is.data.frame(obj) || is.matrix(obj)) {
        write.csv(obj, file.path(dir, paste0(prefix, "_", nm, ".csv")))
      } else {
        saveRDS(obj, file.path(rdsDir, paste0(prefix, "_", nm, ".rds")))
      }
    }
  }
  invisible(NULL)
}

as_metric_table <- function(MF) {
  if (is.null(MF) || !is.list(MF)) return(data.frame())
  all_species <- unique(unlist(lapply(MF, names)))
  if (length(all_species) == 0) {
    max_len <- max(vapply(MF, length, integer(1)))
    all_species <- paste0("sp", seq_len(max_len))
  }
  out <- data.frame(species = all_species, stringsAsFactors = FALSE)
  for (nm in names(MF)) {
    val <- MF[[nm]]
    if (is.null(names(val))) names(val) <- all_species[seq_along(val)]
    out[[nm]] <- as.numeric(val[match(out$species, names(val))])
  }
  out
}

prediction_mean_matrix <- function(predY, Y) {
  if (is.null(predY)) return(NULL)
  d <- dim(predY)
  if (is.null(d)) return(NULL)
  if (length(d) == 2) {
    M <- predY
  } else if (length(d) == 3) {
    if (d[1] == nrow(Y) && d[2] == ncol(Y)) {
      M <- apply(predY, c(1, 2), mean, na.rm = TRUE)
    } else if (d[2] == nrow(Y) && d[3] == ncol(Y)) {
      M <- apply(predY, c(2, 3), mean, na.rm = TRUE)
    } else {
      # fallback: average over the last dimension
      M <- apply(predY, c(1, 2), mean, na.rm = TRUE)
    }
  } else {
    return(NULL)
  }
  rownames(M) <- rownames(Y)
  colnames(M) <- colnames(Y)
  M
}

find_latest_model_file <- function(modelDir = modelDir) {
  files <- list.files(modelDir, pattern = "^models_thin_.*_samples_.*_chains_.*\\.RData$", full.names = TRUE)
  if (length(files) == 0) stop("No fitted model file found in: ", modelDir, call. = FALSE)
  files[which.max(file.info(files)$mtime)]
}

safe_pdf <- function(file, expr, width = 9, height = 7) {
  grDevices::pdf(file, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  try(force(expr), silent = TRUE)
}

safe_write_csv <- function(x, file) {
  if (is.null(x)) return(invisible(FALSE))
  if (is.data.frame(x) || is.matrix(x)) {
    write.csv(x, file)
    return(invisible(TRUE))
  }
  invisible(FALSE)
}

get_mcmc_diag_one <- function(x, label) {
  if (!inherits(x, "mcmc.list")) return(NULL)
  ess <- try(coda::effectiveSize(x), silent = TRUE)
  psrf <- try(coda::gelman.diag(x, multivariate = FALSE)$psrf, silent = TRUE)
  par_names <- unique(c(names(ess), rownames(psrf)))
  if (length(par_names) == 0) par_names <- paste0(label, "_", seq_along(ess))
  out <- data.frame(parameter = par_names, block = label, stringsAsFactors = FALSE)
  out$effective_size <- NA_real_
  out$psrf_point <- NA_real_
  out$psrf_upper <- NA_real_
  if (!inherits(ess, "try-error")) out$effective_size <- as.numeric(ess[match(out$parameter, names(ess))])
  if (!inherits(psrf, "try-error")) {
    out$psrf_point <- as.numeric(psrf[match(out$parameter, rownames(psrf)), 1])
    out$psrf_upper <- as.numeric(psrf[match(out$parameter, rownames(psrf)), 2])
  }
  out
}

collect_mcmc_diagnostics <- function(obj, prefix = "") {
  out <- list()
  if (inherits(obj, "mcmc.list")) {
    out[[prefix]] <- get_mcmc_diag_one(obj, prefix)
  } else if (is.list(obj)) {
    nms <- names(obj)
    if (is.null(nms)) nms <- as.character(seq_along(obj))
    for (i in seq_along(obj)) {
      lab <- if (nzchar(prefix)) paste(prefix, nms[i], sep = "_") else nms[i]
      out <- c(out, collect_mcmc_diagnostics(obj[[i]], lab))
    }
  }
  out
}

plot_mcmc_blocks <- function(obj, prefix = "") {
  if (inherits(obj, "mcmc.list")) {
    try(plot(obj, main = prefix), silent = TRUE)
  } else if (is.list(obj)) {
    nms <- names(obj)
    if (is.null(nms)) nms <- as.character(seq_along(obj))
    for (i in seq_along(obj)) {
      lab <- if (nzchar(prefix)) paste(prefix, nms[i], sep = "_") else nms[i]
      plot_mcmc_blocks(obj[[i]], lab)
    }
  }
}

classify_covariate_group <- function(cov_names, env_cols, hist_groups) {
  labels <- rep("Other", length(cov_names))
  labels[cov_names %in% c("(Intercept)", "Intercept")] <- "Intercept"
  labels[cov_names %in% env_cols] <- "Environment"
  for (g in names(hist_groups)) {
    if (g %in% c("trait_history", "random_history")) next
    labels[cov_names %in% hist_groups[[g]]] <- paste0("History: ", g)
  }
  labels
}

save_post_estimate <- function(post, prefix) {
  if (is.null(post)) return(invisible(NULL))
  if (!is.null(post$mean)) write.csv(post$mean, file.path(tableDir, paste0(prefix, "_mean.csv")))
  if (!is.null(post$support)) write.csv(post$support, file.path(tableDir, paste0(prefix, "_support.csv")))
  saveRDS(post, file.path(rdsDir, paste0(prefix, ".rds")))
}
