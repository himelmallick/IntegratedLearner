test_that("utility helpers for outcomes and metrics behave correctly", {
  expect_identical(IntegratedLearner:::coerce_binary_truth(c(0, 1, 1, 0)), c(0L, 1L, 1L, 0L))
  expect_identical(IntegratedLearner:::coerce_binary_truth(c(FALSE, TRUE)), c(0L, 1L))
  expect_identical(IntegratedLearner:::coerce_binary_truth(factor(c("no", "yes", "yes"))), c(0L, 1L, 1L))
  expect_error(
    IntegratedLearner:::coerce_binary_truth(c("a", "b", "c")),
    "exactly two outcome classes"
  )

  preds <- cbind(layer1 = c(0.1, 0.2, 0.8, 0.9), layer2 = c(0.3, 0.4, 0.7, 0.6))
  metrics <- IntegratedLearner:::binary_model_metrics(preds, c(0, 0, 1, 1))
  expect_true(all(c("auc", "accuracy", "balanced_accuracy", "metrics") %in% names(metrics)))
  expect_true(all(is.finite(metrics$accuracy)))

  auc_vals <- IntegratedLearner:::binary_auc_values_raw(preds, c(0, 0, 1, 1))
  expect_true(all(is.finite(auc_vals)))

  scores_bin <- IntegratedLearner:::single_layer_fusion_scores(preds, c(0, 0, 1, 1), "binomial")
  expect_identical(scores_bin$metric, "AUC")
  expect_identical(names(scores_bin$scores), colnames(preds))

  preds_g <- cbind(layer1 = c(1, 2, 3, 4), layer2 = c(4, 3, 2, 1))
  scores_g <- IntegratedLearner:::single_layer_fusion_scores(preds_g, c(1, 2, 3, 4), "gaussian")
  expect_identical(scores_g$metric, "R2")
  expect_true(all(is.finite(scores_g$scores)))

  expect_error(
    IntegratedLearner:::single_layer_fusion_scores(preds, c(0, 1, 0, 1), "poisson"),
    "Unsupported family"
  )

  fm <- data.frame(
    featureID = c("f1", "f2", "f3"),
    featureType = c("A", "B", "A"),
    row.names = c("f1", "f2", "f3"),
    stringsAsFactors = FALSE
  )
  expect_identical(IntegratedLearner:::feature_ids_for_layers(fm, "A"), c("f1", "f3"))
  expect_identical(IntegratedLearner:::feature_ids_for_layers(fm, character(0)), character(0))
  expect_equal(IntegratedLearner:::expit(0), 0.5)
})

test_that("metadata and screening utility helpers validate inputs", {
  sm <- data.frame(
    outcome = c("case", "control", "case"),
    idcol = c("S1", "S2", "S3"),
    stringsAsFactors = FALSE
  )
  norm <- IntegratedLearner:::normalize_sample_metadata_columns(
    sample_metadata = sm,
    outcome_col = "outcome",
    subject_id_col = "idcol"
  )
  expect_true(all(c("Y", "subjectID") %in% colnames(norm)))

  coerced_g <- IntegratedLearner:::coerce_outcome_by_family(
    sample_metadata = data.frame(Y = c("1.5", "2.5"), subjectID = c("A", "B")),
    family_name = "gaussian"
  )
  expect_identical(coerced_g$sample_metadata$Y, c(1.5, 2.5))

  coerced_b <- IntegratedLearner:::coerce_outcome_by_family(
    sample_metadata = data.frame(Y = c("yes", "no", "yes"), subjectID = c("A", "B", "C")),
    family_name = "binomial"
  )
  expect_identical(sort(unique(coerced_b$sample_metadata$Y)), c(0, 1))
  expect_identical(coerced_b$binary_levels, c("no", "yes"))

  expect_error(
    IntegratedLearner:::coerce_outcome_by_family(
      sample_metadata = data.frame(Y = c("a", "", "b"), subjectID = c("A", "B", "C")),
      family_name = "binomial"
    ),
    "missing/empty"
  )

  expect_identical(IntegratedLearner:::safe_family_name(stats::gaussian()), "gaussian")
  expect_identical(IntegratedLearner:::safe_family_name("BINOMIAL"), "binomial")
  expect_true(is.na(IntegratedLearner:::safe_family_name(1)))

  rs <- suppressWarnings(IntegratedLearner:::resolve_screening_args(screen_pct = 0.25))
  expect_true(rs$enabled)
  expect_equal(rs$screen_pct, 25)
  expect_error(
    IntegratedLearner:::resolve_screening_args(run_screening = TRUE, screen_pct = NULL),
    "screening is enabled"
  )

  expect_warning(
    IntegratedLearner:::resolve_feature_filter_args(filter_method = "prevalence", filter_pct = 0.2),
    "interpreting as a proportion"
  )
  rf <- suppressWarnings(
    IntegratedLearner:::resolve_feature_filter_args(filter_method = "prevalence", filter_pct = 0.2)
  )
  expect_true(rf$enabled)
  expect_identical(rf$filter_method, "prevalence")
  expect_equal(rf$filter_pct, 20)
  expect_error(
    IntegratedLearner:::resolve_feature_filter_args(filter_method = "bad", filter_pct = 10),
    "Unknown 'filter_method'"
  )

  expect_true(IntegratedLearner:::is_survival_outcome("cox", data.frame(Y = 1)))
  expect_true(IntegratedLearner:::is_survival_outcome("gaussian", data.frame(time = 1:2, event = c(1, 0), Y = 1:2)))
  expect_false(IntegratedLearner:::is_survival_outcome("gaussian", data.frame(Y = 1:2)))

  surv_df <- data.frame(subjectID = c("A", "B"))
  surv_df$Y <- I(survival::Surv(c(1, 2), c(1, 0)))
  surv_out <- IntegratedLearner:::ensure_survival_metadata(surv_df)
  expect_identical(surv_out$time, c(1, 2))
  expect_identical(surv_out$event, c(1, 0))
  expect_error(
    IntegratedLearner:::ensure_survival_metadata(data.frame(Y = c(1, 2))),
    "Survival mode requires"
  )

  expect_true(is.null(IntegratedLearner:::require_package("stats")))
  expect_error(IntegratedLearner:::require_package("definitelyNotAPackage"), "package not found")
  jv <- IntegratedLearner:::java_major_version()
  expect_true(is.na(jv) || (is.numeric(jv) && length(jv) == 1L))

  options(IntegratedLearner.java_warned = FALSE)
  expect_true(is.null(IntegratedLearner:::warn_if_java_bart_machine_mismatch()))
})

test_that("feature filtering and univariate importance helpers work", {
  ft <- as.data.frame(rbind(
    f1 = c(1, 0, 1, 0, 1, 0),
    f2 = c(0, 0, 0, 0, 1, 1),
    f3 = c(1, 2, 3, 4, 5, 6),
    f4 = c(1, 1, 1, 1, 1, 1)
  ))
  colnames(ft) <- paste0("S", 1:6)
  fm <- data.frame(
    featureID = rownames(ft),
    featureType = c("A", "A", "B", "B"),
    row.names = rownames(ft),
    stringsAsFactors = FALSE
  )
  ft_valid <- ft

  filt_prev <- IntegratedLearner:::filter_features_by_method(
    feature_table = ft,
    feature_metadata = fm,
    feature_table_valid = ft_valid,
    filter_method = "prevalence",
    filter_pct = 50
  )
  expect_equal(nrow(filt_prev$feature_table), 2)
  expect_true(all(rownames(filt_prev$feature_table) %in% rownames(ft)))

  filt_var <- IntegratedLearner:::filter_features_by_method(
    feature_table = ft,
    feature_metadata = fm,
    feature_table_valid = ft_valid,
    filter_method = "variance",
    filter_pct = 50
  )
  expect_equal(nrow(filt_var$feature_table), 2)

  filt_legacy <- expect_warning(
    IntegratedLearner:::filter_features_by_prevalence(
      feature_table = ft,
      feature_metadata = fm,
      prevalence_pct = 50
    ),
    NA
  )
  expect_equal(nrow(filt_legacy$feature_table), 2)

  ft_missing <- ft_valid[setdiff(rownames(ft_valid), rownames(filt_prev$feature_table)[1]), , drop = FALSE]
  expect_error(
    IntegratedLearner:::filter_features_by_method(
      feature_table = ft,
      feature_metadata = fm,
      feature_table_valid = ft_missing,
      filter_method = "prevalence",
      filter_pct = 50
    ),
    "Validation feature table is missing"
  )

  sm <- data.frame(Y = c(0, 0, 1, 1, 1, 1), subjectID = colnames(ft), row.names = colnames(ft))
  imp_bin <- suppressWarnings(
    IntegratedLearner:::compute_signed_univariate_importance(ft, sm, fm, stats::binomial())
  )
  expect_true(is.numeric(imp_bin$all))
  expect_true(is.list(imp_bin$by_layer))

  sm$Y <- c(0.2, 0.4, 0.7, 1.2, 1.5, 1.8)
  imp_g <- IntegratedLearner:::compute_signed_univariate_importance(ft, sm, fm, stats::gaussian())
  expect_true(is.numeric(imp_g$all))

  expect_identical(IntegratedLearner:::infer_input_mode(MAE_train = NULL, PCL_train = list(a = 1)), "PCL")
  expect_identical(IntegratedLearner:::infer_input_mode(MAE_train = list(a = 1), PCL_train = NULL), "MAE")
  expect_error(
    IntegratedLearner:::infer_input_mode(MAE_train = list(a = 1), PCL_train = list(a = 1)),
    "either MAE_\\* inputs or PCL_\\* inputs"
  )
})

test_that("learner wrappers and print method return expected structures", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("glmnetUtils")
  skip_if_not_installed("quadprog")
  skip_if_not_installed("nloptr")
  skip_if_not_installed("SuperLearner")

  set.seed(100)
  X <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
  Yg <- X$x1 - 0.3 * X$x2 + rnorm(20, sd = 0.2)
  Yb <- c(rep(0, 10), rep(1, 10))
  wts <- rep(1, nrow(X))

  fit_glmnet2 <- IntegratedLearner::sl_glmnet2(
    Y = Yg, X = X, newX = X[1:4, ], family = stats::gaussian(),
    obsWeights = wts, id = seq_len(nrow(X)), nfolds = 3, nlambda = 5
  )
  expect_equal(nrow(fit_glmnet2$pred), 4)

  fit_lasso <- suppressWarnings(IntegratedLearner::sl_lasso(
    Y = Yb, X = X, newX = X[1:4, ], family = stats::binomial(),
    obsWeights = wts, id = seq_len(nrow(X)), nfolds = 3, nlambda = 5
  ))
  expect_equal(nrow(fit_lasso$pred), 4)

  fit_enet <- IntegratedLearner::sl_enet(
    Y = Yg, X = X, newX = X[1:4, ], family = stats::gaussian(),
    obsWeights = wts, id = seq_len(nrow(X)), alpha = c(0, 0.5), nfolds = 3, nlambda = 5
  )
  expect_equal(nrow(fit_enet$pred), 4)
  expect_true(is.finite(IntegratedLearner:::get_cvm(fit_enet$fit$object)))

  expect_true(is.finite(IntegratedLearner::auc_obj(
    b = c(0.5, 0.5),
    X = data.frame(m1 = c(0.1, 0.4, 0.8, 0.9), m2 = c(0.2, 0.3, 0.7, 0.95)),
    Y = c(0, 0, 1, 1)
  )))

  nnls_fit <- IntegratedLearner::NNLS(
    x = cbind(c(0.1, 0.4, 0.6, 0.9), c(0.2, 0.5, 0.7, 0.8)),
    y = c(0.15, 0.45, 0.65, 0.85),
    wt = rep(1, 4)
  )
  expect_equal(sum(nnls_fit$solution), 1, tolerance = 1e-6)

  sl_nnls_g <- IntegratedLearner::sl_nnls_auc(
    Y = Yg,
    X = data.frame(m1 = runif(20), m2 = runif(20)),
    newX = data.frame(m1 = runif(4), m2 = runif(4)),
    family = stats::gaussian(),
    obsWeights = wts
  )
  expect_equal(length(sl_nnls_g$pred), 4)
  pred_nnls_g <- IntegratedLearner:::predict.sl_nnls_auc(sl_nnls_g$fit, newdata = data.frame(m1 = runif(3), m2 = runif(3)))
  expect_equal(length(pred_nnls_g), 3)

  sl_nnls_b <- IntegratedLearner::sl_nnls_auc(
    Y = Yb,
    X = data.frame(m1 = runif(20), m2 = runif(20)),
    newX = data.frame(m1 = runif(4), m2 = runif(4)),
    family = stats::binomial(),
    obsWeights = wts
  )
  expect_equal(length(sl_nnls_b$pred), 4)

  pcl <- make_toy_pcl(n_samples = 20, n_features = 8, seed = 55, binary = TRUE)
  fit <- suppressWarnings(IntegratedLearner::IntegratedLearner(
    PCL_train = pcl,
    folds = 2, seed = 2026, base_learner = "SL.glm", run_stacked = TRUE,
    run_concat = TRUE, print_learner = FALSE, verbose = FALSE, family = stats::binomial()
  ))
  expect_true(length(capture.output(print(fit))) > 0)
})
