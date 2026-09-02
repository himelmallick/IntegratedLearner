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
  skip_if_not_installed("SuperLearner")
  skip_if_not_installed("multiview")

  pcl <- make_toy_pcl(n_samples = 36, n_features = 12, seed = 301, binary = FALSE)

  fit <- suppressWarnings(IntegratedLearner::IntegratedLearner(
    PCL_train = pcl,
    folds = 3,
    seed = 2026,
    base_learner = "SL.mean",
    run_stacked = FALSE,
    run_concat = FALSE,
    run_intermediate = TRUE,
    cooperative_rho = c(0, 0.25),
    print_learner = FALSE,
    family = stats::gaussian()
  ))

  expect_true(isTRUE(fit$run_intermediate))
  expect_true("cooperative" %in% colnames(fit$yhat.train))
  expect_true(is.finite(fit$R2.train[["cooperative"]]))
  expect_true(fit$cooperative_rho %in% fit$cooperative_rho_grid)
})

test_that("IntegratedLearner adds cooperative predictions for binary outcomes", {
  skip_if_not_installed("SuperLearner")
  skip_if_not_installed("multiview")

  pcl <- make_toy_pcl(n_samples = 36, n_features = 12, seed = 302, binary = TRUE)

  fit <- suppressWarnings(IntegratedLearner::IntegratedLearner(
    PCL_train = pcl,
    PCL_valid = pcl,
    folds = 3,
    seed = 2026,
    base_learner = "SL.mean",
    run_stacked = FALSE,
    run_concat = FALSE,
    run_intermediate = TRUE,
    cooperative_rho = c(0, 0.25),
    print_learner = FALSE,
    family = stats::binomial()
  ))

  expect_true("cooperative" %in% colnames(fit$yhat.train))
  expect_true("cooperative" %in% colnames(fit$yhat.test))
  expect_true(is.finite(fit$AUC.train[["cooperative"]]))
  expect_true(fit$AUC.train[["cooperative"]] >= 0)
  expect_true(fit$AUC.train[["cooperative"]] <= 1)
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
    base_learner = "surv.coxph",
    run_intermediate = TRUE,
    cooperative_rho = c(0, 0.25),
    verbose = FALSE
  ))

  expect_true(isTRUE(fit$run_intermediate))
  expect_true(is.list(fit$train_out$cooperative))
  expect_true(is.finite(fit$train_out$cooperative$train_cindex))
  expect_equal(length(fit$train_out$cooperative$train_risk), nrow(tcga$train$sample_metadata))
  expect_true(is.list(fit$valid_out$cooperative))
  expect_equal(length(fit$valid_out$cooperative$valid_risk), nrow(tcga$valid$sample_metadata))
})
