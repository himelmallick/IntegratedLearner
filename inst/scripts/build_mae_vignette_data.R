# Usage: Rscript inst/scripts/build_mae_vignette_data.R SOURCE_DIR [OUTPUT_DIR]
# SOURCE_DIR contains the seven original MAE .RData files (or existing .rds
# examples). Historical originals can be recovered from Git before data/ removal.
# Selection uses training assays only: first 25 finite, nonconstant rows in their
# original order. All samples and assay values of the retained rows are kept.
args <- commandArgs(trailingOnly = TRUE)
source_dir <- if (length(args)) args[[1]] else "inst/extdata"
output_dir <- if (length(args) > 1L) args[[2]] else "inst/extdata"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_source <- function(name) {
  rds <- file.path(source_dir, paste0(name, ".rds"))
  if (file.exists(rds)) return(readRDS(rds))
  env <- new.env(parent = emptyenv())
  load(file.path(source_dir, paste0(name, ".RData")), envir = env)
  env[[name]]
}

compact_pair <- function(train_name, valid_name = NULL, max_features = 25L) {
  train <- read_source(train_name)
  valid <- if (!is.null(valid_name)) read_source(valid_name) else NULL
  assays <- MultiAssayExperiment::experiments(train)
  valid_assays <- if (!is.null(valid)) MultiAssayExperiment::experiments(valid) else NULL
  for (layer in names(assays)) {
    matrix <- SummarizedExperiment::assay(assays[[layer]])
    eligible <- apply(matrix, 1L, function(x) all(is.finite(x)) && stats::var(x) > 0)
    ids <- head(rownames(matrix)[eligible], max_features)
    stopifnot(length(ids) > 0L)
    assays[[layer]] <- assays[[layer]][ids, , drop = FALSE]
    if (!is.null(valid)) {
      stopifnot(all(ids %in% rownames(valid_assays[[layer]])))
      valid_assays[[layer]] <- valid_assays[[layer]][ids, , drop = FALSE]
    }
  }
  MultiAssayExperiment::experiments(train) <- assays
  stopifnot(methods::validObject(train))
  saveRDS(train, file.path(output_dir, paste0(train_name, ".rds")), compress = "xz")
  if (!is.null(valid)) {
    MultiAssayExperiment::experiments(valid) <- valid_assays
    stopifnot(methods::validObject(valid))
    saveRDS(valid, file.path(output_dir, paste0(valid_name, ".rds")), compress = "xz")
  }
}

compact_pair("PRISM_MAE", "NLIBD_MAE")
compact_pair("pregnancy_MAE")
compact_pair("Franzosa_MAE_train", "Franzosa_MAE_valid")
compact_pair("TCGA_survival_MAE_train", "TCGA_survival_MAE_valid")
