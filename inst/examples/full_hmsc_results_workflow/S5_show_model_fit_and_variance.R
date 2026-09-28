# S5_show_model_fit_and_variance.R -------------------------------------------
# Input: fitted HMSC models.
# Output: variance partitioning, associations, and fit plots.

library(Hmsc)

if (!exists("localDir")) localDir <- file.path(getwd(), "hmscHist_full_results_example")
modelDir <- file.path(localDir, "models")
dataDir <- file.path(localDir, "data")
resultDir <- file.path(localDir, "results")
plotDir <- file.path(resultDir, "plots")
tableDir <- file.path(resultDir, "tables")
for (d in c(resultDir, plotDir, tableDir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

model_file <- list.files(modelDir, pattern = "models_thin_.*\\.RData$", full.names = TRUE)
if (length(model_file) == 0) {
  stop("No fitted model file found. Please run S2 first.", call. = FALSE)
}
model_file <- model_file[length(model_file)]
load(model_file)
load(file.path(dataDir, "prepared_data.RData"))

make_vp_group <- function(m, env_names, hist_names) {
  Xn <- colnames(m$X)
  group <- rep(NA_integer_, length(Xn))
  groupnames <- c("Intercept", "Environment", "History")
  group[grepl("Intercept", Xn, fixed = TRUE) | Xn == "(Intercept)"] <- 1
  for (i in seq_along(Xn)) {
    if (is.na(group[i])) {
      if (any(vapply(env_names, function(z) grepl(z, Xn[i], fixed = TRUE), logical(1)))) group[i] <- 2
      if (any(vapply(hist_names, function(z) grepl(z, Xn[i], fixed = TRUE), logical(1)))) group[i] <- 3
    }
  }
  group[is.na(group)] <- 2
  list(group = group, groupnames = groupnames)
}

safe_pdf <- function(file, expr) {
  grDevices::pdf(file)
  on.exit(grDevices::dev.off(), add = TRUE)
  try(capture.output(force(expr)), silent = TRUE)
}

for (nm in names(models)) {
  m <- models[[nm]]

  VP_default <- try(computeVariancePartitioning(m), silent = TRUE)
  if (!inherits(VP_default, "try-error")) {
    saveRDS(VP_default, file.path(resultDir, paste0("variance_partitioning_default_", nm, ".rds")))
    write.csv(VP_default$vals, file.path(tableDir, paste0("variance_partitioning_default_vals_", nm, ".csv")))
    if (!is.null(VP_default$R2T$Beta)) {
      write.csv(VP_default$R2T$Beta, file.path(tableDir, paste0("variance_partitioning_R2T_Beta_", nm, ".csv")))
    }
    if (!is.null(VP_default$R2T$Y)) {
      write.csv(VP_default$R2T$Y, file.path(tableDir, paste0("variance_partitioning_R2T_Y_", nm, ".csv")))
    }
    safe_pdf(
      file.path(plotDir, paste0("variance_partitioning_default_", nm, ".pdf")),
      plotVariancePartitioning(m, VP = VP_default, main = paste("Variance partitioning:", nm))
    )
  } else {
    writeLines(as.character(VP_default), file.path(resultDir, paste0("variance_partitioning_default_failed_", nm, ".txt")))
  }

  env_names <- colnames(env)
  hist_names <- setdiff(colnames(XData), env_names)
  gh <- make_vp_group(m, env_names, hist_names)
  VP_grouped <- try(computeVariancePartitioning(m, group = gh$group, groupnames = gh$groupnames), silent = TRUE)
  if (!inherits(VP_grouped, "try-error")) {
    saveRDS(VP_grouped, file.path(resultDir, paste0("variance_partitioning_env_history_", nm, ".rds")))
    write.csv(VP_grouped$vals, file.path(tableDir, paste0("variance_partitioning_env_history_vals_", nm, ".csv")))
    safe_pdf(
      file.path(plotDir, paste0("variance_partitioning_env_history_", nm, ".pdf")),
      plotVariancePartitioning(m, VP = VP_grouped, main = paste("Environment vs history:", nm))
    )
  }

  assoc <- try(computeAssociations(m), silent = TRUE)
  if (!inherits(assoc, "try-error")) {
    saveRDS(assoc, file.path(resultDir, paste0("species_associations_", nm, ".rds")))
    for (k in seq_along(assoc)) {
      if (!is.null(assoc[[k]]$mean)) {
        write.csv(assoc[[k]]$mean, file.path(tableDir, paste0("association_mean_", nm, "_level", k, ".csv")))
      }
      if (!is.null(assoc[[k]]$support)) {
        write.csv(assoc[[k]]$support, file.path(tableDir, paste0("association_support_", nm, "_level", k, ".csv")))
      }
      mat <- assoc[[k]]$mean
      if (!is.null(mat)) {
        if (!is.null(assoc[[k]]$support)) {
          supportLevel <- 0.95
          mat <- ((assoc[[k]]$support > supportLevel) + (assoc[[k]]$support < (1 - supportLevel)) > 0) * mat
        }
        safe_pdf(file.path(plotDir, paste0("association_heatmap_", nm, "_level", k, ".pdf")), {
          if (requireNamespace("corrplot", quietly = TRUE)) {
            corrplot::corrplot(mat, method = "color",
                               col = grDevices::colorRampPalette(c("blue", "white", "red"))(200),
                               title = paste("Association:", nm, "level", k), mar = c(0, 0, 2, 0))
          } else {
            graphics::image(t(mat[nrow(mat):1, ]), main = paste("Association:", nm, "level", k), axes = FALSE)
          }
        })
      }
    }
  }

  pred_file <- file.path(tableDir, paste0("predicted_richness_", nm, ".csv"))
  if (file.exists(pred_file)) {
    rich <- read.csv(pred_file)
    safe_pdf(file.path(plotDir, paste0("observed_vs_predicted_richness_", nm, ".pdf")), {
      plot(rich$observed_richness, rich$predicted_richness, pch = 19,
           xlab = "Observed richness", ylab = "Predicted richness",
           main = paste("Observed vs predicted richness:", nm))
      abline(0, 1, lty = 2)
    })
  }
}

cat("S5 done: variance partitioning, associations, and fit plots were written.\n")
