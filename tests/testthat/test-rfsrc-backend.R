test_that("required randomForestSRC backend fits and predicts survival risk", {
  withr::local_seed(29)
  x <- matrix(rnorm(160), ncol = 4,
              dimnames = list(paste0("s", 1:40), paste0("f", 1:4)))
  times <- rexp(40, rate = exp(x[, 1] / 3)) + 0.1
  events <- rep(c(1, 1, 0, 1), 10)
  fit <- IntegratedLearner:::fit_surv_model(
    "surv.rfsrc", x, times, events,
    method_args = list(ntree = 25, nodesize = 3, importance = FALSE, seed = -29)
  )
  expect_s3_class(fit$model, "rfsrc")
  expect_identical(fit$feature_names, colnames(x))
  risk <- IntegratedLearner:::predict_surv_risk("surv.rfsrc", fit, x[1:5, , drop = FALSE])
  expect_length(risk, 5)
  expect_true(all(is.finite(risk)))
})
