make_cooperative_toy_multiclass_pcl <- function(
  n_samples = 36, n_features = 12, n_classes = 3,
  seed = 91L
) {
  set.seed(seed)

  sample_ids <- sprintf("S%03d", seq_len(n_samples))
  feature_ids <- sprintf("F%03d", seq_len(n_features))

  mat <- matrix(stats::rnorm(n_features * n_samples),
    nrow = n_features, ncol = n_samples,
    dimnames = list(feature_ids, sample_ids)
  )

  signal_a <- colMeans(mat[1:4, , drop = FALSE])
  signal_b <- colMeans(mat[5:8, , drop = FALSE])
  score <- cbind(signal_a, signal_b, -signal_a - signal_b)
  y <- apply(score, 1, function(z) paste0("C", which.max(z)))

  sample_metadata <- data.frame(
    Y = y, subjectID = sample_ids, row.names = sample_ids,
    stringsAsFactors = FALSE
  )

  feature_metadata <- data.frame(
    featureID = feature_ids,
    featureType = rep(c("omicsA", "omicsB"), length.out = n_features),
    row.names = feature_ids,
    stringsAsFactors = FALSE
  )

  list(
    feature_table = as.data.frame(mat),
    sample_metadata = sample_metadata,
    feature_metadata = feature_metadata
  )
}

test_that("IntegratedLearner adds cooperative predictions for gaussian outcomes", {
  skip_if_not_installed("multiview")

  pcl <- make_toy_pcl(n_samples = 36, n_features = 12, seed = 301, binary = FALSE)

  fit <- suppressWarnings(IntegratedLearner::IntegratedLearner(
    PCL_train = pcl,
    folds = 3,
    seed = 2026,
    base_learner = "not_a_real_base_learner",
    meta_learner = "not_a_real_meta_learner",
    run_intermediate = TRUE,
    cooperative_rho = c(0, 0.25),
    print_learner = FALSE,
    family = stats::gaussian()
  ))

  expect_true(isTRUE(fit$run_intermediate))
  expect_identical(colnames(fit$yhat.train), "cooperative")
  expect_false(isTRUE(fit$run_stacked))
  expect_false(isTRUE(fit$run_concat))
  expect_null(fit$base_learner)
  expect_null(fit$meta_learner)
  expect_true(is.finite(fit$R2.train[["cooperative"]]))
  expect_true(fit$cooperative_rho %in% fit$cooperative_rho_grid)
  expect_true(is.numeric(fit$cooperative_feature_importance))
  expect_true(is.numeric(fit$cooperative_feature_importance_signed))
  expect_equal(length(fit$cooperative_feature_importance), nrow(pcl$feature_table))
  expect_equal(length(fit$cooperative_feature_importance_signed), nrow(pcl$feature_table))
  expect_setequal(names(fit$cooperative_feature_importance_by_layer), unique(pcl$feature_metadata$featureType))
  expect_setequal(names(fit$cooperative_feature_importance_signed_by_layer), unique(pcl$feature_metadata$featureType))
  expect_identical(fit$feature_importance_signed, fit$cooperative_feature_importance_signed)
  expect_identical(fit$feature_importance_signed_by_layer, fit$cooperative_feature_importance_signed_by_layer)
})

test_that("IntegratedLearner adds cooperative predictions for binary outcomes", {
  skip_if_not_installed("multiview")

  pcl <- make_toy_pcl(n_samples = 36, n_features = 12, seed = 302, binary = TRUE)

  fit <- suppressWarnings(IntegratedLearner::IntegratedLearner(
    PCL_train = pcl,
    PCL_valid = pcl,
    folds = 3,
    seed = 2026,
    base_learner = "not_a_real_base_learner",
    meta_learner = "not_a_real_meta_learner",
    run_intermediate = TRUE,
    cooperative_rho = c(0, 0.25),
    print_learner = FALSE,
    family = stats::binomial()
  ))

  expect_identical(colnames(fit$yhat.train), "cooperative")
  expect_identical(colnames(fit$yhat.test), "cooperative")
  expect_null(fit$model_fits$model_layers)
  expect_null(fit$model_fits$model_stacked)
  expect_null(fit$model_fits$model_concat)
  expect_true(is.finite(fit$AUC.train[["cooperative"]]))
  expect_true(fit$AUC.train[["cooperative"]] >= 0)
  expect_true(fit$AUC.train[["cooperative"]] <= 1)
  expect_true(is.numeric(fit$cooperative_feature_importance))
  expect_true(is.numeric(fit$cooperative_feature_importance_signed))
  expect_equal(length(fit$cooperative_feature_importance), nrow(pcl$feature_table))
  expect_equal(length(fit$cooperative_feature_importance_signed), nrow(pcl$feature_table))
  expect_setequal(names(fit$cooperative_feature_importance_by_layer), unique(pcl$feature_metadata$featureType))
  expect_setequal(names(fit$cooperative_feature_importance_signed_by_layer), unique(pcl$feature_metadata$featureType))
})

test_that("IntegratedLearner ignores cooperative learning for multiclass outcomes", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("multiview")

  pcl <- make_cooperative_toy_multiclass_pcl(n_samples = 36, n_features = 12, seed = 303)

  expect_message(
    fit <- suppressWarnings(IntegratedLearner::IntegratedLearner(
      PCL_train = pcl,
      PCL_valid = pcl,
      folds = 3,
      seed = 2026,
      base_learner = "SL.glmnet",
      run_stacked = FALSE,
      run_concat = FALSE,
      run_intermediate = TRUE,
      cooperative_rho = c(0, 0.25),
      print_learner = FALSE,
      family = stats::binomial()
    )),
    "not supported for multiclass outcomes"
  )

  expect_false(isTRUE(fit$run_intermediate))
  expect_false("cooperative" %in% names(fit$prob.train))
  expect_false("cooperative" %in% names(fit$prob.test))
  expect_false("cooperative" %in% fit$metrics.train$model)
})

test_that("ILsurv adds cooperative survival risk output", {
  skip_if_not_installed("survival")
  skip_if_not_installed("timeROC")
  skip_if_not_installed("multiview")

  tcga <- make_tcga_survival_pcl(
    n_samples = 120, n_train = 90, n_gene_features = 6,
    n_mirna_features = 6, seed = 304
  )

  fit <- suppressWarnings(IntegratedLearner::IntegratedLearner(
    PCL_train = tcga$train,
    PCL_valid = tcga$valid,
    folds = 3,
    seed = 2026,
    run_intermediate = TRUE,
    cooperative_rho = c(0, 0.25),
    verbose = FALSE
  ))

  expect_true(isTRUE(fit$run_intermediate))
  expect_null(fit$base_learner)
  expect_null(fit$train_out$single)
  expect_null(fit$train_out$early)
  expect_null(fit$train_out$late)
  expect_true(is.list(fit$train_out$cooperative))
  expect_true(is.finite(fit$train_out$cooperative$train_cindex))
  expect_equal(length(fit$train_out$cooperative$train_risk), nrow(tcga$train$sample_metadata))
  expect_null(fit$valid_out$single)
  expect_null(fit$valid_out$early)
  expect_null(fit$valid_out$late)
  expect_true(is.list(fit$valid_out$cooperative))
  expect_equal(length(fit$valid_out$cooperative$valid_risk), nrow(tcga$valid$sample_metadata))
  expect_true(is.numeric(fit$train_out$cooperative$feature_importance))
  expect_true(is.numeric(fit$train_out$cooperative$feature_importance_signed))
  expect_equal(length(fit$train_out$cooperative$feature_importance), nrow(tcga$train$feature_table))
  expect_equal(length(fit$train_out$cooperative$feature_importance_signed), nrow(tcga$train$feature_table))
  expect_setequal(
    names(fit$train_out$cooperative$feature_importance_by_layer),
    unique(tcga$train$feature_metadata$featureType)
  )
  expect_identical(fit$cooperative_feature_importance, fit$train_out$cooperative$feature_importance)
  expect_identical(
    fit$cooperative_feature_importance_signed_by_layer,
    fit$train_out$cooperative$feature_importance_signed_by_layer
  )
})
