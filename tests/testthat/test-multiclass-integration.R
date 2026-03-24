make_toy_multiclass_pcl <- function(
    n_samples = 36,
    n_features = 12,
    n_classes = 3,
    seed = 91L
) {
  set.seed(seed)

  sample_ids <- sprintf("S%03d", seq_len(n_samples))
  feature_ids <- sprintf("F%03d", seq_len(n_features))

  mat <- matrix(
    stats::rnorm(n_features * n_samples),
    nrow = n_features,
    ncol = n_samples,
    dimnames = list(feature_ids, sample_ids)
  )

  signal_a <- colMeans(mat[1:4, , drop = FALSE])
  signal_b <- colMeans(mat[5:8, , drop = FALSE])
  score <- cbind(signal_a, signal_b, -signal_a - signal_b)
  y <- apply(score, 1, function(z) paste0("C", which.max(z)))

  sample_metadata <- data.frame(
    Y = y,
    subjectID = sample_ids,
    row.names = sample_ids,
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

test_that("IntegratedLearner multiclass path runs with native multiclass backend", {
  skip_if_not_installed("glmnet")

  pcl <- make_toy_multiclass_pcl(n_samples = 30, n_features = 12, seed = 2026)

  fit <- suppressWarnings(IntegratedLearner::IntegratedLearner(
    PCL_train = pcl,
    folds = 2,
    seed = 2026,
    base_learner = "SL.glmnet",
    run_stacked = TRUE,
    run_concat = TRUE,
    print_learner = FALSE,
    family = stats::binomial()
  ))

  expect_identical(fit$family, "multinomial")
  expect_true(is.list(fit$prob.train))
  expect_true(all(c("omicsA", "omicsB", "stacked", "concatenated") %in% names(fit$prob.train)))
  expect_true(is.data.frame(fit$metrics.train))
  expect_true(all(c("model", "accuracy", "balanced_accuracy", "logloss") %in% colnames(fit$metrics.train)))

  p <- fit$prob.train[["stacked"]]
  expect_true(is.matrix(p) || is.data.frame(p))
  expect_equal(nrow(p), nrow(pcl$sample_metadata))
  expect_true(all(abs(rowSums(as.matrix(p)) - 1) < 1e-6))
})

test_that("predict.learner returns multiclass probabilities and metrics", {
  skip_if_not_installed("glmnet")

  pcl <- make_toy_multiclass_pcl(n_samples = 33, n_features = 12, seed = 2027)

  fit <- suppressWarnings(IntegratedLearner::IntegratedLearner(
    PCL_train = pcl,
    folds = 2,
    seed = 2027,
    base_learner = "SL.glmnet",
    run_stacked = TRUE,
    run_concat = TRUE,
    print_learner = FALSE,
    family = stats::binomial()
  ))

  pred <- IntegratedLearner:::predict.learner(
    object = fit,
    feature_table_valid = pcl$feature_table,
    sample_metadata_valid = pcl$sample_metadata,
    feature_metadata = pcl$feature_metadata
  )

  expect_true(is.list(pred$prob.test))
  expect_true(is.data.frame(pred$class.test))
  expect_true(is.data.frame(pred$metrics.test))
  expect_true(all(c("accuracy", "balanced_accuracy", "logloss") %in% colnames(pred$metrics.test)))
})

test_that("update.learner is blocked with explicit message for multiclass fits", {
  skip_if_not_installed("glmnet")

  pcl <- make_toy_multiclass_pcl(n_samples = 27, n_features = 10, seed = 2028)

  fit <- suppressWarnings(IntegratedLearner::IntegratedLearner(
    PCL_train = pcl,
    folds = 2,
    seed = 2028,
    base_learner = "SL.glmnet",
    run_stacked = TRUE,
    run_concat = FALSE,
    print_learner = FALSE,
    family = stats::binomial()
  ))

  expect_error(
    IntegratedLearner:::update.learner(
      object = fit,
      feature_table_valid = pcl$feature_table,
      sample_metadata_valid = pcl$sample_metadata,
      feature_metadata_valid = pcl$feature_metadata
    ),
    "update.learner\\(\\) for multiclass fits is not implemented yet",
    perl = TRUE
  )
})

test_that("multiclass learner aliases map to native backends", {
  expect_identical(IntegratedLearner:::.map_multiclass_learner("SL.glmnet"), "glmnet")
  expect_identical(IntegratedLearner:::.map_multiclass_learner("SL.ranger"), "ranger")
  expect_identical(IntegratedLearner:::.map_multiclass_learner("SL.randomForest"), "randomforest")
  expect_identical(IntegratedLearner:::.map_multiclass_learner("randomforest"), "randomforest")
  expect_identical(IntegratedLearner:::.map_multiclass_learner("xgboost"), "xgboost")
  expect_identical(IntegratedLearner:::.map_multiclass_learner("mbart"), "mbart")
  expect_identical(IntegratedLearner:::.map_multiclass_learner("multinom"), "multinom")

  expect_identical(IntegratedLearner:::.map_multiclass_meta_learner("SL.nnls.auc"), "glmnet")
  expect_identical(IntegratedLearner:::.map_multiclass_meta_learner("randomforest"), "randomforest")
  expect_identical(IntegratedLearner:::.map_multiclass_meta_learner("xgboost"), "xgboost")
  expect_identical(IntegratedLearner:::.map_multiclass_meta_learner("mbart"), "mbart")
  expect_identical(IntegratedLearner:::.map_multiclass_meta_learner("multinom"), "multinom")
})

test_that("IntegratedLearner multiclass supports xgboost backend", {
  skip_if_not_installed("xgboost")

  pcl <- make_toy_multiclass_pcl(n_samples = 30, n_features = 12, seed = 2029)

  fit <- suppressWarnings(IntegratedLearner::IntegratedLearner(
    PCL_train = pcl,
    folds = 2,
    seed = 2029,
    base_learner = "xgboost",
    meta_learner = "xgboost",
    run_stacked = TRUE,
    run_concat = FALSE,
    print_learner = FALSE,
    family = stats::binomial()
  ))

  expect_identical(fit$base_learner_used, "xgboost")
  expect_identical(fit$meta_learner_used, "xgboost")
  expect_true(is.data.frame(fit$metrics.train))
  expect_true(all(c("accuracy", "balanced_accuracy", "logloss") %in% colnames(fit$metrics.train)))
})

test_that("IntegratedLearner multiclass supports randomForest backend", {
  skip_if_not_installed("randomForest")

  pcl <- make_toy_multiclass_pcl(n_samples = 30, n_features = 12, seed = 2031)

  fit <- suppressWarnings(IntegratedLearner::IntegratedLearner(
    PCL_train = pcl,
    folds = 2,
    seed = 2031,
    base_learner = "randomforest",
    meta_learner = "randomforest",
    run_stacked = TRUE,
    run_concat = FALSE,
    print_learner = FALSE,
    family = stats::binomial()
  ))

  expect_identical(fit$base_learner_used, "randomforest")
  expect_identical(fit$meta_learner_used, "randomforest")
  expect_true(is.data.frame(fit$metrics.train))
  expect_true(all(c("accuracy", "balanced_accuracy", "logloss") %in% colnames(fit$metrics.train)))
})

test_that("IntegratedLearner multiclass supports mbart backend", {
  skip_if_not_installed("BART")

  pcl <- make_toy_multiclass_pcl(n_samples = 24, n_features = 10, seed = 2033)

  fit <- suppressWarnings(IntegratedLearner::IntegratedLearner(
    PCL_train = pcl,
    folds = 2,
    seed = 2033,
    base_learner = "mbart",
    meta_learner = "glmnet",
    run_stacked = FALSE,
    run_concat = FALSE,
    print_learner = FALSE,
    family = stats::binomial(),
    ntree = 10,
    ndpost = 20,
    nskip = 10,
    keepevery = 1
  ))

  expect_identical(fit$base_learner_used, "mbart")
  expect_true(is.data.frame(fit$metrics.train))
  expect_true(all(c("accuracy", "balanced_accuracy", "logloss") %in% colnames(fit$metrics.train)))
})

test_that("mbart prediction aligns columns to training features", {
  skip_if_not_installed("BART")

  set.seed(2040)
  X <- as.data.frame(matrix(rnorm(30 * 6), nrow = 30, ncol = 6))
  colnames(X) <- paste0("f", seq_len(ncol(X)))
  y <- factor(rep(c("A", "B", "C"), length.out = nrow(X)), levels = c("A", "B", "C"))

  fit_obj <- suppressWarnings(IntegratedLearner:::.fit_multiclass_model_impl(
    X = X,
    y = y,
    learner_id = "mbart",
    seed = 2040,
    model_args = list(
      ntree = 10L,
      ndpost = 20L,
      nskip = 10L,
      keepevery = 1L,
      printevery = 1000L
    )
  ))

  X_new <- X[, c("f3", "f1", "f2", "f4", "f5", "f6"), drop = FALSE]
  X_new$extra_col <- rnorm(nrow(X_new))

  prob <- suppressWarnings(IntegratedLearner:::.predict_multiclass_model_impl(
    fit_obj = fit_obj,
    newX = X_new,
    class_levels = levels(y)
  ))

  expect_equal(nrow(prob), nrow(X_new))
  expect_equal(ncol(prob), length(levels(y)))
  expect_identical(colnames(prob), levels(y))
})

test_that("mbart handles constant columns without predict-time dimension mismatch", {
  skip_if_not_installed("BART")

  set.seed(2041)
  X <- as.data.frame(matrix(rnorm(45 * 8), nrow = 45, ncol = 8))
  colnames(X) <- paste0("f", seq_len(ncol(X)))
  X$f7 <- 1
  X$f8 <- 0
  y <- factor(rep(c("A", "B", "C"), length.out = nrow(X)), levels = c("A", "B", "C"))

  fit_obj <- suppressWarnings(IntegratedLearner:::.fit_multiclass_model_impl(
    X = X,
    y = y,
    learner_id = "mbart",
    seed = 2041,
    model_args = list(
      ntree = 10L,
      ndpost = 20L,
      nskip = 10L,
      keepevery = 1L,
      printevery = 1000L
    )
  ))

  prob <- suppressWarnings(IntegratedLearner:::.predict_multiclass_model_impl(
    fit_obj = fit_obj,
    newX = X,
    class_levels = levels(y)
  ))

  expect_equal(nrow(prob), nrow(X))
  expect_equal(ncol(prob), length(levels(y)))
  expect_identical(colnames(prob), levels(y))
})
