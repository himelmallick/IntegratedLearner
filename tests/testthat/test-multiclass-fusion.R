make_fusion_probabilities <- function() {
  good <- matrix(c(.8, .1, .1, .1, .8, .1, .1, .1, .8),
    nrow = 3, byrow = TRUE, dimnames = list(c("s1", "s2", "s3"), c("A", "B", "C")))
  list(good = good, bad = good[, c("B", "C", "A")])
}

test_that("majority voting counts hard votes and resolves ties in class order", {
  probs <- make_fusion_probabilities()
  # Keep labels fixed while moving high probabilities to the wrong classes.
  colnames(probs$bad) <- colnames(probs$good)
  vote <- IntegratedLearner:::multiclass_majority_probs(probs)
  expect_equal(unname(rowSums(vote)), rep(1, 3))
  expect_equal(vote[1, ], c(A = .5, B = 0, C = .5))
  expect_identical(IntegratedLearner:::class_from_prob(vote), c("A", "A", "B"))
  probs$third <- probs$good
  expect_identical(IntegratedLearner:::class_from_prob(
    IntegratedLearner:::multiclass_majority_probs(probs)), c("A", "B", "C"))
  expect_equal(unname(IntegratedLearner:::multiclass_majority_probs(probs[1])), diag(3))
})

test_that("learned averaging favors the useful layer and aligns layer and class names", {
  probs <- make_fusion_probabilities()
  colnames(probs$bad) <- colnames(probs$good)
  y <- c("A", "B", "C")
  weights <- IntegratedLearner:::learn_multiclass_weights(probs, y)
  expect_named(weights, names(probs))
  expect_true(all(weights >= 0))
  expect_equal(sum(weights), 1)
  expect_gt(weights[["good"]], .99)
  combined <- IntegratedLearner:::multiclass_weighted_probs(probs, weights)
  expect_gt(mean(log(diag(combined))), mean(log(diag((probs$good + probs$bad) / 2))))
  expect_equal(IntegratedLearner:::multiclass_weighted_probs(
    list(bad = probs$bad[, c("C", "A", "B")], good = probs$good),
    weights[c("bad", "good")])[, c("A", "B", "C")], combined)
  expect_equal(IntegratedLearner:::learn_multiclass_weights(probs[1], y), c(good = 1))
  expect_equal(IntegratedLearner:::multiclass_weighted_probs(probs[1], c(good = 1)), probs$good)
})

test_that("learned averaging handles complementary, identical and zero probabilities", {
  p <- matrix(c(1, 0, 0, 0, 1, 0), 2, byrow = TRUE,
    dimnames = list(NULL, c("A", "B", "C")))
  q <- p[2:1, , drop = FALSE]
  weights <- IntegratedLearner:::learn_multiclass_weights(list(p = p, q = q), c("A", "A"))
  expect_equal(unname(weights), c(.5, .5), tolerance = 1e-4)
  identical_weights <- IntegratedLearner:::learn_multiclass_weights(
    list(p = p, q = p, r = p), c("A", "B"))
  expect_equal(unname(identical_weights), rep(1 / 3, 3), tolerance = 1e-4)
})

test_that("fusion rejects mismatched inputs and invalid weights", {
  probs <- make_fusion_probabilities()
  expect_error(IntegratedLearner:::multiclass_majority_probs(list()), "nonempty")
  broken <- probs
  broken$bad <- broken$bad[1:2, ]
  expect_error(IntegratedLearner:::multiclass_majority_probs(broken), "matching samples")
  broken <- probs
  broken$bad[1, 1] <- NA_real_
  expect_error(IntegratedLearner:::multiclass_majority_probs(broken), "finite")
  expect_error(IntegratedLearner:::multiclass_weighted_probs(probs, c(good = 1)), "weights")
  expect_error(IntegratedLearner:::multiclass_weighted_probs(probs, c(good = 2, bad = -1)), "weights")
  expect_error(IntegratedLearner:::learn_multiclass_weights(probs, c("A", "B", "D")), "outcomes")
})

make_multiclass_fusion_pcl <- function() {
  set.seed(2042)
  y <- rep(c("A", "B", "C"), each = 20)
  x <- matrix(rnorm(12 * length(y)), 12,
    dimnames = list(paste0("f", 1:12), paste0("s", seq_along(y))))
  x[1, ] <- x[1, ] + as.integer(factor(y))
  list(feature_table = as.data.frame(x),
    sample_metadata = data.frame(Y = y, subjectID = colnames(x), row.names = colnames(x)),
    feature_metadata = data.frame(featureID = rownames(x),
      featureType = rep(c("rna", "protein"), each = 6), row.names = rownames(x)))
}

test_that("fusion integrates with fitting, cross-fitting, prediction and serialization", {
  skip_if_not_installed("nnet")
  pcl <- make_multiclass_fusion_pcl()
  fit <- IntegratedLearner::IntegratedLearner(
    PCL_train = pcl, PCL_valid = pcl, folds = 3, seed = 2042,
    base_learner = "multinom", run_stacked = FALSE, run_concat = FALSE,
    family = stats::binomial(), print_learner = FALSE)
  methods <- c("majority", "weighted_average")
  expect_true(all(methods %in% names(fit$prob.train)))
  expect_true(all(methods %in% fit$metrics.train$model))
  expect_true(all(methods %in% fit$metrics.test$model))
  expect_equal(sum(fit$model_fits$model_weighted$weights), 1)
  for (f in unique(fit$fold_id)) {
    held_out <- fit$fold_id == f
    layers <- fit$prob.train[names(fit$model_fits$model_layers)]
    weights <- IntegratedLearner:::learn_multiclass_weights(
      lapply(layers, function(p) p[!held_out, , drop = FALSE]), fit$Y_train[!held_out])
    expect_equal(fit$prob.train$weighted_average[held_out, , drop = FALSE],
      IntegratedLearner:::multiclass_weighted_probs(
        lapply(layers, function(p) p[held_out, , drop = FALSE]), weights))
  }
  restored <- unserialize(serialize(fit, NULL))
  pred <- IntegratedLearner:::predict_multiclass.learner(restored,
    feature_table_valid = pcl$feature_table, sample_metadata_valid = pcl$sample_metadata,
    feature_metadata = pcl$feature_metadata)
  expect_equal(pred$prob.test$weighted_average,
    IntegratedLearner:::multiclass_weighted_probs(
      pred$prob.test[names(fit$model_fits$model_layers)],
      fit$model_fits$model_weighted$weights))
  for (method in methods) {
    expect_equal(pred$prob.test[[method]], fit$prob.test[[method]])
    expect_equal(unname(rowSums(pred$prob.test[[method]])), rep(1, 60))
    expect_identical(pred$class.test[[method]], fit$class.test[[method]])
  }
  # Older serialized fits without these flags retain their original outputs.
  restored$run_majority <- restored$run_weighted <- NULL
  old_pred <- IntegratedLearner:::predict_multiclass.learner(restored,
    feature_table_valid = pcl$feature_table, feature_metadata = pcl$feature_metadata)
  expect_false(any(methods %in% names(old_pred$prob.test)))
  # Adding the methods does not change layer or stacked outputs in older fits.
  expect_equal(old_pred$prob.test, pred$prob.test[names(old_pred$prob.test)])
})

test_that("fusion flags can be disabled independently", {
  skip_if_not_installed("nnet")
  pcl <- make_multiclass_fusion_pcl()
  for (majority in c(FALSE, TRUE)) {
    fit <- IntegratedLearner::IntegratedLearner(
      PCL_train = pcl, folds = 3, seed = 2042, base_learner = "multinom",
      run_stacked = FALSE, run_concat = FALSE, run_majority = majority,
      run_weighted = !majority, family = stats::binomial(), print_learner = FALSE)
    expect_identical("majority" %in% names(fit$prob.train), majority)
    expect_identical("weighted_average" %in% names(fit$prob.train), !majority)
  }
  fit <- IntegratedLearner::IntegratedLearner(
    PCL_train = pcl, folds = 3, seed = 2042, base_learner = "multinom",
    run_stacked = FALSE, run_concat = FALSE, run_majority = FALSE,
    run_weighted = FALSE, family = stats::binomial(), print_learner = FALSE)
  expect_setequal(names(fit$prob.train), c("rna", "protein"))
  expect_null(fit$model_fits$model_weighted)
})

test_that("a single layer supports both fusion methods", {
  skip_if_not_installed("nnet")
  pcl <- make_multiclass_fusion_pcl()
  pcl$feature_metadata$featureType <- "rna"
  fit <- IntegratedLearner::IntegratedLearner(
    PCL_train = pcl, PCL_valid = pcl, folds = 3, seed = 2042,
    base_learner = "multinom", run_stacked = FALSE, run_concat = FALSE,
    family = stats::binomial(), print_learner = FALSE)
  expect_equal(fit$model_fits$model_weighted$weights, c(rna = 1))
  expect_equal(fit$prob.train$weighted_average, fit$prob.train$rna)
  expect_equal(fit$prob.test$weighted_average, fit$prob.test$rna)
  expect_identical(fit$class.train$majority, fit$class.train$rna)
  expect_identical(fit$class.test$majority, fit$class.test$rna)
})

test_that("fusion flags require single logical values", {
  expect_error(IntegratedLearner::IL_multiclass(run_majority = NA), "run_majority")
  expect_error(IntegratedLearner::IL_multiclass(run_weighted = c(TRUE, FALSE)), "run_weighted")
  expect_error(IntegratedLearner::IL_multiclass(run_weighted = 1), "run_weighted")
})

test_that("all multiclass methods have survival-style grouped outputs", {
  skip_if_not_installed("nnet")
  pcl <- make_multiclass_fusion_pcl()
  fit <- IntegratedLearner::IntegratedLearner(
    PCL_train = pcl, PCL_valid = pcl, folds = 3, seed = 2042,
    base_learner = "multinom", meta_learner = "multinom",
    family = stats::binomial(), print_learner = FALSE)
  methods <- c("stacked", "majority", "weighted_average")
  expect_identical(fit$late_methods, methods)
  for (split in c("train", "valid")) {
    grouped <- fit[[paste0(split, "_out")]]
    suffix <- if (split == "train") "train" else "test"
    expect_named(grouped, c("single", "early", "late"))
    expect_named(grouped$late, methods)
    expect_equal(grouped$single[[paste0(split, "_prob")]],
      fit[[paste0("prob.", suffix)]][names(fit$model_fits$model_layers)])
    expect_equal(grouped$early[[paste0(split, "_prob")]],
      fit[[paste0("prob.", suffix)]]$concatenated)
    for (method in methods) {
      out <- grouped$late[[method]]
      expect_equal(out[[paste0(split, "_prob")]], fit[[paste0("prob.", suffix)]][[method]])
      expect_identical(out[[paste0(split, "_class")]], fit[[paste0("class.", suffix)]][[method]])
      metrics <- fit[[paste0("metrics.", suffix)]]
      expect_equal(out$metrics, metrics[metrics$model == method, , drop = FALSE])
      expect_equal(out[[paste0(split, "_logloss")]], out$metrics$logloss)
    }
    expect_equal(grouped$late$weighted_average$weights, fit$model_fits$model_weighted$weights)
  }
  expect_named(fit$train_out$late$weighted_average$weights_by_fold,
    as.character(sort(unique(fit$fold_id))))
  for (f in names(fit$train_out$late$weighted_average$weights_by_fold)) {
    expect_equal(sum(fit$train_out$late$weighted_average$weights_by_fold[[f]]), 1)
  }
  pred <- IntegratedLearner:::predict_multiclass.learner(fit,
    feature_table_valid = pcl$feature_table, sample_metadata_valid = pcl$sample_metadata,
    feature_metadata = pcl$feature_metadata)
  expect_equal(pred$valid_out, fit$valid_out)
  expect_identical(pred$late_methods, methods)
  unlabeled <- IntegratedLearner:::predict_multiclass.learner(fit,
    feature_table_valid = pcl$feature_table, feature_metadata = pcl$feature_metadata)
  expect_equal(unlabeled$valid_out$late$weighted_average$valid_prob,
    fit$prob.test$weighted_average)
  expect_null(unlabeled$valid_out$late$weighted_average$metrics)
  expect_null(unlabeled$valid_out$late$weighted_average$valid_accuracy)
})

test_that("grouped outputs respect disabled methods and missing validation", {
  skip_if_not_installed("nnet")
  pcl <- make_multiclass_fusion_pcl()
  fit <- IntegratedLearner::IntegratedLearner(
    PCL_train = pcl, folds = 3, seed = 2042, base_learner = "multinom",
    run_stacked = FALSE, run_concat = FALSE, run_majority = FALSE, run_weighted = FALSE,
    family = stats::binomial(), print_learner = FALSE)
  expect_null(fit$valid_out)
  expect_true("valid_out" %in% names(fit))
  expect_null(fit$train_out$early)
  expect_length(fit$train_out$late, 0)
  expect_length(fit$late_methods, 0)
  expect_setequal(fit$train_out$single$metrics$model, c("rna", "protein"))
})
