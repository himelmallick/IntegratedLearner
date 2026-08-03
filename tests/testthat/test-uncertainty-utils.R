test_that("VarImp.learner rejects unsupported base learners", {
  skip_if_not_installed("bartMachine")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tibble")
  skip_if_not_installed("stringr")

  fit <- list(
    base_learner = "SL.glm",
    meta_learner = "sl_nnls_auc",
    weights = c(omicsA = 0.6, omicsB = 0.4),
    model_fits = list(model_layers = list())
  )

  expect_error(
    IntegratedLearner:::VarImp.learner(fit),
    "currently available only for BART base learner",
    fixed = TRUE
  )
})

test_that("VarImp.learner validates requested layer names before BART calls", {
  skip_if_not_installed("bartMachine")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tibble")
  skip_if_not_installed("stringr")

  fit <- list(
    base_learner = "sl_bart",
    meta_learner = "sl_nnls_auc",
    weights = c(omicsA = 1),
    model_fits = list(model_layers = list(omicsA = list(dummy = TRUE)))
  )

  expect_error(
    IntegratedLearner:::VarImp.learner(fit, layer.names = "missing_layer"),
    "Invalid layer name(s): missing_layer.",
    fixed = TRUE
  )
})

test_that("VarImp.learner returns a ggplot for mocked BART fits", {
  skip_if_not_installed("bartMachine")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tibble")
  skip_if_not_installed("stringr")

  fit <- list(
    base_learner = "sl_bart",
    meta_learner = "sl_nnls_auc",
    weights = c(layer1 = 0.7, layer2 = 0.3),
    model_fits = list(
      model_layers = list(
        layer1 = structure(list(), class = "bartMachine"),
        layer2 = structure(list(), class = "bartMachine")
      )
    )
  )

  testthat::local_mocked_bindings(
    investigate_var_importance = function(object, plot = FALSE) {
      list(
        avg_var_props = c(f_1 = 0.6, f_2 = 0.4),
        sd_var_props = c(f_1 = 0.05, f_2 = 0.03)
      )
    },
    .package = "bartMachine"
  )

  p <- IntegratedLearner:::VarImp.learner(fit, num.var = 2)
  expect_s3_class(p, "ggplot")
})

test_that("omicsEye_theme returns a ggplot theme object", {
  skip_if_not_installed("ggplot2")

  th <- IntegratedLearner:::omicsEye_theme()
  expect_s3_class(th, "theme")
})

test_that("credint.learner rejects unsupported non-BART configurations", {
  skip_if_not_installed("bayesplot")
  skip_if_not_installed("ggplot2")

  fit <- list(base_learner = "SL.glm", meta_learner = "sl_nnls_auc")

  expect_error(
    IntegratedLearner:::credint.learner(fit = fit),
    "currently only available for BART as base learner and NNLS/AUC as the meta learner",
    fixed = TRUE
  )
})

test_that("credint.learner multiclass path requires at least one mbart model", {
  skip_if_not_installed("BART")
  skip_if_not_installed("bayesplot")
  skip_if_not_installed("ggplot2")

  fit <- list(
    family = "multinomial",
    test = FALSE,
    X_train_layers = list(omicsA = data.frame(x = 1:3)),
    Y_train = factor(c("A", "B", "C")),
    prob.train = list(omicsA = matrix(
      c(0.8, 0.1, 0.1,
        0.2, 0.7, 0.1,
        0.1, 0.2, 0.7),
      nrow = 3, byrow = TRUE,
      dimnames = list(c("S1", "S2", "S3"), c("A", "B", "C"))
    )),
    model_fits = list(model_layers = list(omicsA = list(learner_id = "glmnet")))
  )

  expect_error(
    IntegratedLearner:::credint.learner(fit = fit, test = FALSE),
    "at least one fitted model must use 'mbart'",
    fixed = TRUE
  )
})

test_that("credint.learner multiclass path validates requested class early", {
  skip_if_not_installed("BART")
  skip_if_not_installed("bayesplot")
  skip_if_not_installed("ggplot2")

  fit <- list(
    family = "multinomial",
    test = FALSE,
    X_train_layers = list(omicsA = data.frame(x = c(1, 2, 3), row.names = c("S1", "S2", "S3"))),
    Y_train = factor(c("A", "B", "C"), levels = c("A", "B", "C")),
    prob.train = list(omicsA = matrix(
      c(0.8, 0.1, 0.1,
        0.2, 0.7, 0.1,
        0.1, 0.2, 0.7),
      nrow = 3, byrow = TRUE,
      dimnames = list(c("S1", "S2", "S3"), c("A", "B", "C"))
    )),
    class_levels = c("A", "B", "C"),
    model_fits = list(model_layers = list(omicsA = list(learner_id = "mbart")))
  )

  expect_error(
    IntegratedLearner:::credint.learner(
      fit = fit,
      test = FALSE,
      model = "omicsA",
      class = "Z"
    ),
    "Requested class 'Z' is not among class levels: A, B, C",
    fixed = TRUE
  )
})
