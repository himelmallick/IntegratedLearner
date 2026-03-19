test_that("infer_mae_defaults selects experiments and assays", {
  mae <- make_toy_mae(n_samples = 12, n_features_layer1 = 6, n_features_layer2 = 4, seed = 31)
  
  defaults <- IntegratedLearner:::.infer_mae_defaults(mae_train = mae)
  expect_identical(defaults$experiment, c("taxonomy", "pathway"))
  expect_identical(unname(defaults$assay.type), c("relative_abundance", "pathway_abundance"))
  
  defaults_idx <- IntegratedLearner:::.infer_mae_defaults(
    mae_train = mae,
    experiment = c(2L, 1L),
    assay.type = c("pathway_abundance", "relative_abundance")
  )
  expect_identical(defaults_idx$experiment, c("pathway", "taxonomy"))
  expect_identical(defaults_idx$assay.type, c("pathway_abundance", "relative_abundance"))
})

test_that("infer_mae_defaults validates experiment and assay arguments", {
  mae <- make_toy_mae(n_samples = 10, n_features_layer1 = 5, n_features_layer2 = 5, seed = 32)
  
  expect_error(
    IntegratedLearner:::.infer_mae_defaults(mae_train = mae, experiment = c(99L)),
    "'experiment' indices are out of range for mae_train.",
    fixed = TRUE
  )
  
  expect_error(
    IntegratedLearner:::.infer_mae_defaults(mae_train = mae, experiment = c("missing_layer")),
    "Some 'experiment' entries are not found in mae_train.",
    fixed = TRUE
  )
  
  mae_multi_assay <- make_toy_mae(
    n_samples = 10,
    n_features_layer1 = 5,
    n_features_layer2 = 5,
    seed = 33,
    multi_assay_layer1 = TRUE
  )
  
  expect_error(
    IntegratedLearner:::.infer_mae_defaults(mae_train = mae_multi_assay),
    "Multiple assays found for an experiment",
    fixed = TRUE
  )
})

test_that("get_data_from_MAE handles survival fallback and required colData columns", {
  mae <- make_toy_mae(
    n_samples = 12,
    n_features_layer1 = 5,
    n_features_layer2 = 4,
    seed = 34,
    outcome = "survival",
    include_time_event = TRUE
  )
  
  exps <- MultiAssayExperiment::experiments(mae)
  for (nm in names(exps)) {
    cd <- as.data.frame(SummarizedExperiment::colData(exps[[nm]]))
    cd$Y <- NULL
    SummarizedExperiment::colData(exps[[nm]]) <- S4Vectors::DataFrame(cd)
  }
  MultiAssayExperiment::experiments(mae) <- exps
  
  long_df <- IntegratedLearner:::.get_data_from_MAE(
    mae = mae,
    experiment = c("taxonomy", "pathway"),
    assay.type = c("relative_abundance", "pathway_abundance"),
    outcome.col = "Y"
  )
  
  expect_true(is.data.frame(long_df))
  expect_true(all(c("sampleID", "subjectID", "Y", "time", "event", "featureID", "experiment") %in% colnames(long_df)))
  expect_true(nrow(long_df) > 0)
  
  mae_missing_subject <- make_toy_mae(n_samples = 10, n_features_layer1 = 5, n_features_layer2 = 4, seed = 35)
  exps2 <- MultiAssayExperiment::experiments(mae_missing_subject)
  for (nm in names(exps2)) {
    cd <- as.data.frame(SummarizedExperiment::colData(exps2[[nm]]))
    cd$subjectID <- NULL
    SummarizedExperiment::colData(exps2[[nm]]) <- S4Vectors::DataFrame(cd)
  }
  MultiAssayExperiment::experiments(mae_missing_subject) <- exps2
  
  expect_error(
    IntegratedLearner:::.get_data_from_MAE(
      mae = mae_missing_subject,
      experiment = c("taxonomy", "pathway"),
      assay.type = c("relative_abundance", "pathway_abundance"),
      outcome.col = "Y"
    ),
    "colData(se) must contain a 'subjectID' column.",
    fixed = TRUE
  )
})

test_that("prepare_from_MAE returns canonical objects and validates feature identity", {
  mae <- make_toy_mae(n_samples = 14, n_features_layer1 = 6, n_features_layer2 = 5, seed = 36)
  
  out <- IntegratedLearner:::.prepare_from_MAE(
    mae_train = mae,
    experiment = c("taxonomy", "pathway"),
    assay.type = c("relative_abundance", "pathway_abundance"),
    na.rm = FALSE
  )
  
  expect_true(is.data.frame(out$feature_table))
  expect_true(is.data.frame(out$sample_metadata))
  expect_true(is.data.frame(out$feature_metadata))
  expect_identical(rownames(out$feature_table), rownames(out$feature_metadata))
  expect_identical(colnames(out$feature_table), rownames(out$sample_metadata))
  
  mae_valid <- mae
  exps_valid <- MultiAssayExperiment::experiments(mae_valid)
  se_tax <- exps_valid$taxonomy
  new_names <- rownames(se_tax)
  new_names[1] <- "different_feature"
  rownames(se_tax) <- new_names
  exps_valid$taxonomy <- se_tax
  MultiAssayExperiment::experiments(mae_valid) <- exps_valid
  
  expect_error(
    IntegratedLearner:::.prepare_from_MAE(
      mae_train = mae,
      mae_valid = mae_valid,
      experiment = c("taxonomy", "pathway"),
      assay.type = c("relative_abundance", "pathway_abundance")
    ),
    "Feature sets differ between mae_train and mae_valid in experiment 'taxonomy' - they must match exactly.",
    fixed = TRUE
  )
})
