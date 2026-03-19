test_that("IntegratedLearner binomial stacked + concat returns expected structure", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  
  pcl <- make_toy_pcl(n_samples = 28, n_features = 12, seed = 41, binary = TRUE)
  
  fit <- suppressWarnings(IntegratedLearner::IntegratedLearner(
    PCL_train = pcl,
    folds = 2,
    seed = 2026,
    base_learner = "SL.glm",
    run_stacked = TRUE,
    run_concat = TRUE,
    print_learner = FALSE,
    verbose = FALSE,
    family = stats::binomial()
  ))
  
  expected_cols <- c(unique(pcl$feature_metadata$featureType), "stacked", "concatenated")
  expect_identical(colnames(fit$yhat.train), expected_cols)
  expect_identical(names(fit$AUC.train), expected_cols)
  expect_true(all(is.finite(fit$AUC.train)))
  
  expect_true(!is.null(fit$weights))
  expect_equal(sum(fit$weights), 1, tolerance = 1e-6)
  expect_identical(sort(names(fit$weights)), sort(unique(pcl$feature_metadata$featureType)))
  
  expect_true(is.numeric(fit$feature_importance_signed))
  expect_true(is.list(fit$feature_importance_signed_by_layer))
})

test_that("IntegratedLearner gaussian path computes R2 train and validation", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  
  pcl <- make_toy_pcl(n_samples = 30, n_features = 12, seed = 42, binary = FALSE)
  train_ids <- rownames(pcl$sample_metadata)[1:22]
  valid_ids <- rownames(pcl$sample_metadata)[23:30]
  
  pcl_train <- subset_pcl(pcl, feature_ids = rownames(pcl$feature_table), sample_ids = train_ids)
  pcl_valid <- subset_pcl(pcl, feature_ids = rownames(pcl$feature_table), sample_ids = valid_ids)
  
  fit <- suppressWarnings(IntegratedLearner::IntegratedLearner(
    PCL_train = pcl_train,
    PCL_valid = pcl_valid,
    folds = 2,
    seed = 2026,
    base_learner = "SL.glm",
    run_stacked = FALSE,
    run_concat = FALSE,
    print_learner = FALSE,
    verbose = FALSE,
    family = stats::gaussian()
  ))
  
  expect_identical(fit$family, "gaussian")
  expect_true(all(names(fit$R2.train) %in% colnames(fit$yhat.train)))
  expect_true(all(names(fit$R2.test) %in% colnames(fit$yhat.test)))
  expect_true(all(is.finite(fit$R2.train)))
  expect_true(all(is.finite(fit$R2.test)))
})

test_that("IntegratedLearner runs in MAE mode for binomial outcome", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  
  mae <- make_toy_mae(
    n_samples = 18,
    n_features_layer1 = 7,
    n_features_layer2 = 5,
    seed = 43,
    outcome = "binomial"
  )
  
  fit <- suppressWarnings(IntegratedLearner::IntegratedLearner(
    MAE_train = mae,
    experiment = c("taxonomy", "pathway"),
    assay.type = c("relative_abundance", "pathway_abundance"),
    folds = 2,
    seed = 2026,
    base_learner = "SL.glm",
    run_stacked = FALSE,
    run_concat = FALSE,
    print_learner = FALSE,
    verbose = FALSE,
    family = stats::binomial()
  ))
  
  expect_identical(fit$input_mode, "MAE")
  expect_identical(fit$family, "binomial")
  expect_equal(NROW(fit$yhat.train), 18)
  expect_equal(NCOL(fit$yhat.train), 2)
})

test_that("predict.learner validates inputs and returns aligned predictions", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  
  pcl <- make_toy_pcl(n_samples = 26, n_features = 12, seed = 44, binary = TRUE)
  
  fit <- suppressWarnings(IntegratedLearner::IntegratedLearner(
    PCL_train = pcl,
    folds = 2,
    seed = 2026,
    base_learner = "SL.glm",
    run_stacked = TRUE,
    run_concat = TRUE,
    print_learner = FALSE,
    verbose = FALSE,
    family = stats::binomial()
  ))
  
  valid_ids <- rownames(pcl$sample_metadata)[1:8]
  feature_table_valid <- pcl$feature_table[, valid_ids, drop = FALSE]
  sample_metadata_valid <- pcl$sample_metadata[valid_ids, , drop = FALSE]
  
  pred <- suppressWarnings(IntegratedLearner:::predict.learner(
    object = fit,
    feature_table_valid = feature_table_valid,
    sample_metadata_valid = sample_metadata_valid,
    feature_metadata = pcl$feature_metadata
  ))
  
  expect_true(is.data.frame(pred$yhat.test) || is.matrix(pred$yhat.test))
  expect_equal(NROW(pred$yhat.test), length(valid_ids))
  expect_true(all(c("stacked", "concatenated") %in% colnames(pred$yhat.test)))
  expect_true(all(is.finite(pred$AUC.test)))
  
  bad_ft <- feature_table_valid
  rownames(bad_ft)[1] <- "bad_feature"
  expect_error(
    IntegratedLearner:::predict.learner(
      object = fit,
      feature_table_valid = bad_ft,
      sample_metadata_valid = sample_metadata_valid,
      feature_metadata = pcl$feature_metadata
    ),
    "Both feature_table and feature_table_valid should have the same rownames.",
    fixed = TRUE
  )
})

test_that("update.learner handles full and partial layer overlap", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  
  pcl <- make_toy_pcl(n_samples = 24, n_features = 12, seed = 45, binary = TRUE)
  
  fit <- suppressWarnings(IntegratedLearner::IntegratedLearner(
    PCL_train = pcl,
    folds = 2,
    seed = 2026,
    base_learner = "SL.glm",
    run_stacked = FALSE,
    run_concat = FALSE,
    print_learner = FALSE,
    verbose = FALSE,
    family = stats::binomial()
  ))
  
  valid_ids <- rownames(pcl$sample_metadata)[1:6]
  feature_table_valid <- pcl$feature_table[, valid_ids, drop = FALSE]
  sample_metadata_valid <- pcl$sample_metadata[valid_ids, , drop = FALSE]
  
  pred_full <- suppressWarnings(IntegratedLearner:::predict.learner(
    object = fit,
    feature_table_valid = feature_table_valid,
    sample_metadata_valid = sample_metadata_valid,
    feature_metadata = pcl$feature_metadata
  ))
  
  capture.output({
    fit_updated_full <- suppressWarnings(IntegratedLearner:::update.learner(
      object = fit,
      feature_table_valid = feature_table_valid,
      sample_metadata_valid = sample_metadata_valid,
      feature_metadata_valid = pcl$feature_metadata,
      verbose = FALSE
    ))
  })
  
  expect_identical(colnames(fit_updated_full$yhat.test), colnames(pred_full$yhat.test))
  expect_equal(NROW(fit_updated_full$yhat.test), NROW(pred_full$yhat.test))
  
  layer_a_features <- rownames(pcl$feature_metadata)[pcl$feature_metadata$featureType == "omicsA"]
  ft_partial <- feature_table_valid[layer_a_features, , drop = FALSE]
  fm_partial <- pcl$feature_metadata[layer_a_features, , drop = FALSE]
  
  capture.output({
    fit_updated_partial <- suppressWarnings(IntegratedLearner:::update.learner(
      object = fit,
      feature_table_valid = ft_partial,
      sample_metadata_valid = sample_metadata_valid,
      feature_metadata_valid = fm_partial,
      verbose = FALSE
    ))
  })
  
  expect_identical(colnames(fit_updated_partial$yhat.test), "omicsA")
  expect_equal(NCOL(fit_updated_partial$yhat.test), 1)
  
  fm_no_overlap <- fm_partial
  fm_no_overlap$featureType <- "no_overlap_layer"
  expect_error(
    IntegratedLearner:::update.learner(
      object = fit,
      feature_table_valid = ft_partial,
      sample_metadata_valid = sample_metadata_valid,
      feature_metadata_valid = fm_no_overlap,
      verbose = FALSE
    ),
    "Validation set has no layers in common with model fit.",
    fixed = TRUE
  )
})

test_that("plot.learner returns plotting payloads for binomial and gaussian fits", {
  skip_if_not_installed("SuperLearner")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("cowplot")
  skip_if_not_installed("stringr")
  suppressPackageStartupMessages(library(SuperLearner))
  
  tmp_plot <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp_plot)
  on.exit({
    if (grDevices::dev.cur() > 1L) grDevices::dev.off()
    unlink(tmp_plot)
  }, add = TRUE)
  
  pcl_bin <- make_toy_pcl(n_samples = 22, n_features = 10, seed = 46, binary = TRUE)
  fit_bin <- suppressWarnings(IntegratedLearner::IntegratedLearner(
    PCL_train = pcl_bin,
    folds = 2,
    seed = 2026,
    base_learner = "SL.glm",
    run_stacked = FALSE,
    run_concat = FALSE,
    print_learner = FALSE,
    verbose = FALSE,
    family = stats::binomial()
  ))
  
  out_bin <- suppressWarnings(IntegratedLearner:::plot.learner(fit_bin))
  expect_true(is.list(out_bin))
  expect_true(all(c("plot", "ROC_table") %in% names(out_bin)))
  expect_true(is.data.frame(out_bin$ROC_table))
  
  pcl_gauss <- make_toy_pcl(n_samples = 22, n_features = 10, seed = 47, binary = FALSE)
  fit_gauss <- suppressWarnings(IntegratedLearner::IntegratedLearner(
    PCL_train = pcl_gauss,
    folds = 2,
    seed = 2026,
    base_learner = "SL.glm",
    run_stacked = FALSE,
    run_concat = FALSE,
    print_learner = FALSE,
    verbose = FALSE,
    family = stats::gaussian()
  ))
  
  out_gauss <- suppressWarnings(IntegratedLearner:::plot.learner(fit_gauss))
  expect_true(is.list(out_gauss))
  expect_true(all(c("plot", "R2_table") %in% names(out_gauss)))
  expect_true(is.data.frame(out_gauss$R2_table))
})
