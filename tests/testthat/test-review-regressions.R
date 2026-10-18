test_that("prediction rejects incomplete, repeated, and reordered feature IDs", {
  feature_ids <- c("f1", "f2", "f3", "f4")
  fit <- structure(list(family = "gaussian", feature.names = feature_ids), class = "learner")
  metadata <- data.frame(featureID = feature_ids, row.names = feature_ids)
  for (ids in list(feature_ids[1:2], rep(feature_ids[1:2], 2), rev(feature_ids), character())) {
    ft <- matrix(0, length(ids), 2, dimnames = list(ids, c("s1", "s2")))
    expect_error(stats::predict(fit, feature_table_valid = ft, feature_metadata = metadata),
                 "identical row names in the same order", fixed = TRUE)
  }
})

test_that("glmnet wrappers preserve predictions and the elastic net search", {
  withr::local_seed(15)
  x <- data.frame(a = rnorm(36), b = rnorm(36), c = rnorm(36))
  y <- 2 * x$a + rnorm(36)
  args <- list(Y = y, X = x, newX = x[1:4, ], family = stats::gaussian(),
               obsWeights = rep(1, 36), id = 1:36, nfolds = 3, nlambda = 10)
  lasso <- withr::with_seed(42, do.call(sl_lasso, args))
  reference <- withr::with_seed(42, do.call(sl_glmnet2, c(args, list(alpha = 1))))
  expect_equal(lasso$pred, reference$pred)
  expect_s3_class(lasso$fit, "SL.glmnet")
  enet <- withr::with_seed(42, do.call(sl_enet, c(args, list(alpha = c(0, 0.5, 1)))))
  expect_true(enet$fit$object$alpha %in% c(0, 0.5, 1))
  expect_true(all(is.finite(enet$pred)))
})

test_that("local CV seeding restores the caller RNG state", {
  withr::local_seed(321)
  before <- .Random.seed
  metadata <- data.frame(subjectID = paste0("s", 1:12))
  first <- IntegratedLearner:::build_subject_cv_partition(metadata, folds = 3, seed = 14)
  expect_identical(.Random.seed, before)
  second <- IntegratedLearner:::build_subject_cv_partition(metadata, folds = 3, seed = 14)
  expect_identical(first$fold_id, second$fold_id)
  expect_identical(.Random.seed, before)
})

test_that("vectorized Brier calculations equal the scalar formula", {
  withr::local_seed(53)
  times <- c(5, 1, 4, 2, 3)
  events <- c(1, 0, 1, 1, 0)
  grid <- c(1, 2, 4)
  csurv <- c(1, 0.8, Inf, 0.6, 0.5)
  csurv_grid <- c(1, 0.7, Inf)
  matrices <- lapply(1:3, function(i) matrix(runif(15), 5, 3))
  params <- c(0.4, -0.3)
  layer_weights <- exp(c(params, 0)) / sum(exp(c(params, 0)))
  combined <- Reduce(`+`, Map(`*`, matrices, layer_weights))
  reference <- vapply(seq_along(grid), function(j) {
    mean(combined[, j]^2 * (times <= grid[j] & events == 1) / csurv +
         (1 - combined[, j])^2 * (times > grid[j]) / csurv_grid[j])
  }, numeric(1))
  weights <- IntegratedLearner:::ipcw_brier_weights(times, events, grid, csurv, csurv_grid)
  expect_equal(IntegratedLearner:::weighted_brier_values(combined, weights), reference)
  objective <- IntegratedLearner:::make_ibs_objective(matrices, grid, times, events, csurv, csurv_grid)
  expect_equal(objective(params), IntegratedLearner:::integrated_brier(grid, reference))
  expect_true(is.finite(objective(c(1000, -1000))))
  # Preserve matrix dimensions for a single subject or evaluation time.
  expect_length(IntegratedLearner:::weighted_brier_values(
    matrix(0.5, 1, 1), IntegratedLearner:::ipcw_brier_weights(1, 1, 2, 1, 1)), 1L)
})
