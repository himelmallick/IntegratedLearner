make_survival_calibration_data <- function() {
  set.seed(306)
  n <- 96L
  x <- matrix(rnorm(n * 4), nrow = 4,
    dimnames = list(paste0("f", 1:4), paste0("s", seq_len(n))))
  latent <- 0.6 * x[1, ] - 0.4 * x[3, ]
  event_time <- rexp(n, rate = exp(latent) / 5)
  censor_time <- rexp(n, rate = 0.08)
  metadata <- data.frame(time = pmin(event_time, censor_time),
    event = as.integer(event_time <= censor_time), subjectID = colnames(x),
    row.names = colnames(x))
  list(x = as.data.frame(x), metadata = metadata,
    features = data.frame(featureID = rownames(x), featureType = rep(c("a", "b"), each = 2),
      row.names = rownames(x)), train = 1:72, valid = 73:96)
}

fit_calibration_fixture <- function(data, validation, cooperative = FALSE) {
  IntegratedLearner::ILsurv(
    feature_table = data$x[, data$train, drop = FALSE],
    sample_metadata = data$metadata[data$train, , drop = FALSE], feature_metadata = data$features,
    valid_feature_table = data$x[, data$valid, drop = FALSE], valid_sample_metadata = validation,
    folds = 3, seed = 306, run_intermediate = cooperative, cooperative_rho = 0
  )
}

test_that("calibration predicts from a frozen training fit, including after serialization", {
  data <- make_survival_calibration_data()
  marker <- as.numeric(data$x[1, data$train])
  calibration <- IntegratedLearner:::fit_survival_calibration(
    marker, data$metadata$time[data$train], data$metadata$event[data$train], c(1, 3, 5))
  risk_new <- c(-1, 0, 1, NA_real_)
  predicted <- IntegratedLearner:::predict_survival_calibration(calibration, risk_new)
  expect_equal(dim(predicted), c(4, 3))
  expect_true(all(is.finite(predicted)))
  expect_true(all(predicted >= 0 & predicted <= 1))
  expect_equal(IntegratedLearner:::predict_survival_calibration(
    unserialize(serialize(calibration, NULL)), risk_new), predicted)
  expect_equal(IntegratedLearner:::predict_survival_calibration(calibration, 0),
    predicted[2, , drop = FALSE])
  testthat::local_mocked_bindings(
    fit_survival_calibration = function(...) stop("Prediction must not refit calibration"),
    .package = "IntegratedLearner"
  )
  expect_equal(IntegratedLearner:::predict_survival_calibration(calibration, risk_new), predicted)
})

test_that("survival scoring consumes probabilities without fitting calibration", {
  times <- c(1, 2, 3, 4, 5, 6)
  events <- c(1, 1, 0, 1, 0, 1)
  marker <- c(.9, .7, .4, .6, .2, .1)
  grid <- c(1, 3, 5)
  probs <- outer(seq(.2, .8, length.out = 6), c(.9, .7, .4))
  expected <- IntegratedLearner:::compute_ipcw_brier(times, events, probs, grid)
  testthat::local_mocked_bindings(
    fit_survival_calibration = function(...) stop("Scoring must not fit calibration"),
    risk_to_surv_matrix = function(...) stop("Scoring must not fit calibration"),
    .package = "IntegratedLearner"
  )
  scored <- IntegratedLearner:::compute_auc_cindex(times, events, marker,
    surv_mat = probs, brier_time_grid = grid)
  expect_equal(scored$brier, expected)
  expect_equal(scored$ibs, IntegratedLearner:::integrated_brier(expected$time, expected$Brier))
  without_probs <- IntegratedLearner:::compute_auc_cindex(times, events, marker)
  expect_true(all(is.na(without_probs$brier$Brier)))
  expect_true(is.na(without_probs$ibs))
  expect_equal(without_probs$cindex, scored$cindex)
  expect_equal(without_probs$auc, scored$auc)
  expect_error(IntegratedLearner:::compute_auc_cindex(times, events, marker,
    surv_mat = probs[-1, ], brier_time_grid = grid), "must match")
})

test_that("Kaplan-Meier calibration fallback also uses only its training fit", {
  calibration <- IntegratedLearner:::fit_survival_calibration(
    rep(1, 6), 1:6, c(1, 0, 1, 0, 1, 0), c(1, 3, 5))
  expect_null(calibration$model)
  predicted <- IntegratedLearner:::predict_survival_calibration(calibration, c(-100, 0, 100))
  expect_equal(predicted[1, ], predicted[2, ])
  expect_equal(predicted[2, ], predicted[3, ])
})

test_that("standard survival probabilities and calibration ignore validation labels", {
  data <- make_survival_calibration_data()
  validation <- data$metadata[data$valid, , drop = FALSE]
  fit <- fit_calibration_fixture(data, validation)
  changed <- validation
  changed$time <- rev(validation$time) * 1.7
  changed$event <- 1 - validation$event
  altered <- fit_calibration_fixture(data, changed)
  expect_identical(fit$brier_time_grid, altered$brier_time_grid)
  expect_equal(fit$valid_out$single$valid_risk, altered$valid_out$single$valid_risk)
  expect_equal(fit$valid_out$single$valid_surv, altered$valid_out$single$valid_surv)
  expect_equal(fit$valid_out$early$valid_surv, altered$valid_out$early$valid_surv)
  expect_false(isTRUE(all.equal(fit$valid_out$single$valid_brier,
    altered$valid_out$single$valid_brier)))
  for (layer in names(fit$survival_calibration$layers)) {
    expect_equal(fit$valid_out$single$valid_surv[[layer]],
      IntegratedLearner:::predict_survival_calibration(fit$survival_calibration$layers[[layer]],
        fit$valid_out$single$valid_risk[[layer]]))
  }
  expect_equal(fit$valid_out$early$valid_surv,
    IntegratedLearner:::predict_survival_calibration(fit$survival_calibration$early,
      fit$valid_out$early$valid_risk))
  for (method in c("IBS", "COX")) {
    expect_equal(fit$valid_out$late[[method]]$valid_surv,
      altered$valid_out$late[[method]]$valid_surv)
    expected <- IntegratedLearner:::compute_ipcw_brier(validation$time, validation$event,
      fit$valid_out$late[[method]]$valid_surv, fit$brier_time_grid)
    expect_equal(fit$valid_out$late[[method]]$valid_brier, expected)
    expect_equal(fit$valid_out$late[[method]]$valid_ibs,
      IntegratedLearner:::integrated_brier(expected$time, expected$Brier))
  }
  expect_equal(fit$valid_out$late$COX$valid_surv,
    IntegratedLearner:::predict_survival_calibration(fit$survival_calibration$late$COX,
      fit$valid_out$late$COX$valid_risk))
  weights <- fit$train_out$late$IBS$weights
  expected_surv <- Reduce(`+`, lapply(names(weights), function(layer) {
    fit$valid_out$single$valid_surv[[layer]] * weights[[layer]]
  }))
  expect_equal(fit$valid_out$late$IBS$valid_surv, expected_surv)
  train <- data$metadata[data$train, , drop = FALSE]
  expected_train <- IntegratedLearner:::compute_ipcw_brier(train$time, train$event,
    fit$train_out$late$IBS$train_surv, fit$brier_time_grid)
  expect_equal(fit$train_out$late$IBS$train_brier, expected_train)
  expect_equal(fit$train_out$late$IBS$train_ibs,
    IntegratedLearner:::integrated_brier(expected_train$time, expected_train$Brier))
})

test_that("cooperative survival validation uses frozen training calibration", {
  skip_if_not_installed("multiview")
  data <- make_survival_calibration_data()
  validation <- data$metadata[data$valid, , drop = FALSE]
  fit <- fit_calibration_fixture(data, validation, cooperative = TRUE)
  changed <- validation
  changed$time <- rev(validation$time) * 1.7
  changed$event <- 1 - validation$event
  altered <- fit_calibration_fixture(data, changed, cooperative = TRUE)
  expect_identical(fit$brier_time_grid, altered$brier_time_grid)
  expect_equal(fit$valid_out$cooperative$valid_risk, altered$valid_out$cooperative$valid_risk)
  expect_equal(fit$valid_out$cooperative$valid_surv, altered$valid_out$cooperative$valid_surv)
  expected <- IntegratedLearner:::predict_survival_calibration(
    fit$survival_calibration$cooperative, fit$valid_out$cooperative$valid_risk)
  expect_equal(fit$valid_out$cooperative$valid_surv, expected)
  brier <- IntegratedLearner:::compute_ipcw_brier(validation$time, validation$event,
    expected, fit$brier_time_grid)
  expect_equal(fit$valid_out$cooperative$valid_brier, brier)
})
