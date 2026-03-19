test_that("IntegratedLearner runs on toy PCL data (binomial)", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  
  pcl <- make_toy_pcl(n_samples = 24, n_features = 12, seed = 10, binary = TRUE)
  
  fit <- IntegratedLearner::IntegratedLearner(
    PCL_train = pcl,
    folds = 2,
    seed = 2026,
    base_learner = "SL.glm",
    run_stacked = FALSE,
    run_concat = FALSE,
    print_learner = FALSE,
    verbose = FALSE,
    family = stats::binomial()
  )
  
  expect_identical(fit$input_mode, "PCL")
  expect_identical(fit$family, "binomial")
  expect_true(is.matrix(fit$yhat.train) || is.data.frame(fit$yhat.train))
  expect_equal(NROW(fit$yhat.train), nrow(pcl$sample_metadata))
  expect_true(all(is.finite(fit$AUC.train)))
  expect_true(all(fit$AUC.train >= 0 & fit$AUC.train <= 1))
})

test_that("IntegratedLearner runs on an iHMP subset with validation set", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  
  path <- testthat::test_path("..", "..", "data", "iHMP.RData")
  skip_if_not(file.exists(path), "Missing fixture: data/iHMP.RData")
  
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  pcl <- get("pcl", envir = env)
  
  y_raw <- as.character(pcl$sample_metadata$Y)
  class_map <- split(rownames(pcl$sample_metadata), y_raw)
  if (length(class_map) < 2) {
    skip("Need at least two classes for binomial smoke test.")
  }
  
  class_sizes <- vapply(class_map, length, integer(1))
  if (any(class_sizes < 6)) {
    skip("Not enough samples per class for train/validation split.")
  }
  
  keep_samples <- sort(unique(unlist(
    lapply(class_map, function(ids) ids[seq_len(min(10, length(ids)))]),
    use.names = FALSE
  )))
  
  layers <- unique(pcl$feature_metadata$featureType)
  keep_features <- sort(unique(unlist(
    lapply(layers, function(layer) head(rownames(pcl$feature_metadata)[pcl$feature_metadata$featureType == layer], 12)),
    use.names = FALSE
  )))
  
  pcl_small <- subset_pcl(pcl, feature_ids = keep_features, sample_ids = keep_samples)
  pcl_small$sample_metadata <- pcl_small$sample_metadata[, c("Y", "subjectID"), drop = FALSE]
  pcl_small$feature_metadata <- pcl_small$feature_metadata[, c("featureID", "featureType"), drop = FALSE]
  
  class_map_small <- split(
    rownames(pcl_small$sample_metadata),
    as.character(pcl_small$sample_metadata$Y)
  )
  train_ids <- sort(unique(unlist(
    lapply(class_map_small, function(ids) {
      ids[seq_len(max(1L, floor(0.75 * length(ids))))]
    }),
    use.names = FALSE
  )))
  valid_ids <- setdiff(rownames(pcl_small$sample_metadata), train_ids)
  
  pcl_train <- subset_pcl(pcl_small, feature_ids = keep_features, sample_ids = train_ids)
  pcl_valid <- subset_pcl(pcl_small, feature_ids = keep_features, sample_ids = valid_ids)
  
  fit <- IntegratedLearner::IntegratedLearner(
    PCL_train = pcl_train,
    PCL_valid = pcl_valid,
    folds = 2,
    seed = 2026,
    base_learner = "SL.glm",
    run_stacked = FALSE,
    run_concat = FALSE,
    print_learner = FALSE,
    verbose = FALSE,
    family = stats::binomial()
  )
  
  expect_true(isTRUE(fit$test))
  expect_true(is.matrix(fit$yhat.test) || is.data.frame(fit$yhat.test))
  expect_equal(NROW(fit$yhat.test), nrow(pcl_valid$sample_metadata))
  expect_true(all(is.finite(fit$AUC.test)))
  expect_true(all(fit$AUC.test >= 0 & fit$AUC.test <= 1))
})
