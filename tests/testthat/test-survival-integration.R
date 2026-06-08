test_that("IntegratedLearner survival mode works with TCGA fixture and validation", {
  skip_if_not_installed("survival")
  skip_if_not_installed("timeROC")

  suppressPackageStartupMessages(library(survival))

  tcga <- make_tcga_survival_pcl(
    n_samples = 220, n_train = 160, n_gene_features = 10,
    n_mirna_features = 8, seed = 2026
  )

  fit <- suppressWarnings(IntegratedLearner::IntegratedLearner(
    PCL_train = tcga$train,
    PCL_valid = tcga$valid, folds = 2, seed = 2026, base_learner = "surv.coxph",
    weight_method = "COX", verbose = FALSE
  ))

  expect_identical(fit$family, "survival")
  expect_identical(fit$input_mode, "PCL")

  expect_true(is.list(fit$train_out))
  expect_true(all(c("single", "early", "late", "intermediate") %in% names(fit$train_out)))
  expect_true(is.list(fit$valid_out))
  expect_true(all(c("single", "early", "late", "intermediate") %in% names(fit$valid_out)))

  expected_layers <- unique(as.character(tcga$train$feature_metadata$featureType))
  expect_equal(length(expected_layers), 2)
  expect_true(all(names(fit$train_out$late$weights) %in% expected_layers))
  expect_equal(sum(fit$train_out$late$weights), 1, tolerance = 1e-06)
  expect_true(is.finite(fit$valid_out$late$valid_cindex))
  expect_true(fit$valid_out$late$valid_cindex >= 0 && fit$valid_out$late$valid_cindex <=
    1)

  expect_true(is.list(fit$train_out$single$metrics))
  expect_equal(length(fit$train_out$single$metrics), 2)
  expect_true(is.list(fit$valid_out$single$valid_auc))
  expect_equal(length(fit$valid_out$single$valid_auc), 2)
})

test_that("IntegratedLearner survival IBS branch returns train metrics without validation", {
  skip_if_not_installed("survival")
  skip_if_not_installed("timeROC")

  suppressPackageStartupMessages(library(survival))

  tcga <- make_tcga_survival_pcl(
    n_samples = 180, n_train = 140, n_gene_features = 8,
    n_mirna_features = 8, seed = 2027
  )

  fit <- suppressWarnings(IntegratedLearner::IntegratedLearner(
    PCL_train = tcga$train,
    folds = 2, seed = 2027, base_learner = "surv.coxph", weight_method = "IBS",
    verbose = FALSE
  ))

  expect_identical(fit$family, "survival")
  expect_null(fit$valid_out)
  expect_true(is.finite(fit$train_out$late$train_cindex))
  expect_true(is.data.frame(fit$train_out$late$train_auc))
  expect_true(all(c("time", "AUC") %in% colnames(fit$train_out$late$train_auc)))
})

test_that("plot.learner returns survival AUC and KM payloads", {
  skip_if_not_installed("survival")
  skip_if_not_installed("timeROC")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("cowplot")
  skip_if_not_installed("stringr")

  suppressPackageStartupMessages(library(survival))

  tmp_plot <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp_plot)
  on.exit(
    {
      if (grDevices::dev.cur() > 1L) grDevices::dev.off()
      unlink(tmp_plot)
    },
    add = TRUE
  )

  tcga <- make_tcga_survival_pcl(
    n_samples = 200, n_train = 150, n_gene_features = 8,
    n_mirna_features = 6, seed = 2030
  )

  fit <- suppressWarnings(IntegratedLearner::IntegratedLearner(
    PCL_train = tcga$train,
    PCL_valid = tcga$valid, folds = 2, seed = 2030,
    base_learner = "surv.coxph", weight_method = "COX", verbose = FALSE
  ))

  out <- suppressWarnings(IntegratedLearner:::plot.learner(fit))
  expect_true(is.list(out))
  expect_true(all(c("plot", "AUC_table_train", "KM_table_train") %in% names(out)))
  expect_true(is.data.frame(out$AUC_table_train))
  expect_true(is.data.frame(out$KM_table_train))
  expect_true(all(c("time", "AUC", "model") %in% colnames(out$AUC_table_train)))
  expect_true(all(c("time", "surv", "strata") %in% colnames(out$KM_table_train)))
})
