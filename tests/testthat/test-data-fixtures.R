test_that("tutorial fixtures are small valid MAEs with usable outcomes", {
  files <- list.files(system.file("extdata", package = "IntegratedLearner"),
                      pattern = "\\.rds$", full.names = TRUE)
  expect_length(files, 7L)
  expect_lt(sum(file.info(files)$size), 500000)
  for (path in files) {
    mae <- readRDS(path)
    expect_s4_class(mae, "MultiAssayExperiment")
    expect_true(methods::validObject(mae))
    metadata <- MultiAssayExperiment::colData(mae)
    expect_gt(nrow(metadata), 10L)
    for (se in as.list(MultiAssayExperiment::experiments(mae))) {
      mat <- SummarizedExperiment::assay(se)
      expect_true(nrow(mat) > 0L && nrow(mat) <= 25L)
      expect_true(all(is.finite(mat)))
      if (!grepl("NLIBD|_valid", basename(path))) {
        expect_true(all(apply(mat, 1L, stats::var) > 0))
      }
      expect_identical(colnames(mat), rownames(metadata))
    }
    outcome <- if ("Y" %in% names(metadata)) metadata$Y else if ("diseaseCat" %in% names(metadata)) {
      metadata$diseaseCat
    } else metadata$event
    expect_false(anyNA(outcome))
    expect_gt(length(unique(outcome)), 1L)
  }
})

test_that("paired MAEs have identical assay feature IDs in training order", {
  pairs <- list(c("PRISM_MAE", "NLIBD_MAE"),
                c("Franzosa_MAE_train", "Franzosa_MAE_valid"),
                c("TCGA_survival_MAE_train", "TCGA_survival_MAE_valid"))
  for (pair in pairs) {
    train <- MultiAssayExperiment::experiments(readRDS(system.file(
      "extdata", paste0(pair[1], ".rds"), package = "IntegratedLearner", mustWork = TRUE
    )))
    valid <- MultiAssayExperiment::experiments(readRDS(system.file(
      "extdata", paste0(pair[2], ".rds"), package = "IntegratedLearner", mustWork = TRUE
    )))
    expect_identical(names(train), names(valid))
    for (layer in names(train)) {
      expect_identical(rownames(train[[layer]]), rownames(valid[[layer]]))
    }
  }
})
