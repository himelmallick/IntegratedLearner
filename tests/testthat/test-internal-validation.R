test_that("infer_input_mode enforces mutually exclusive inputs", {
  pcl <- make_toy_pcl(n_samples = 10, n_features = 8, seed = 1)
  
  expect_identical(
    IntegratedLearner:::.infer_input_mode(MAE_train = list(dummy = TRUE), PCL_train = NULL),
    "MAE"
  )
  
  expect_identical(
    IntegratedLearner:::.infer_input_mode(MAE_train = NULL, PCL_train = pcl),
    "PCL"
  )
  
  expect_error(
    IntegratedLearner:::.infer_input_mode(MAE_train = list(dummy = TRUE), PCL_train = pcl),
    "Provide either MAE_* inputs or PCL_* inputs, not both.",
    fixed = TRUE
  )
  
  expect_error(
    IntegratedLearner:::.infer_input_mode(MAE_train = NULL, PCL_train = NULL),
    "You must supply either MAE_train or PCL_train.",
    fixed = TRUE
  )
})

test_that("prepare_from_PCL validates structure and alignment", {
  pcl <- make_toy_pcl(n_samples = 12, n_features = 10, seed = 2)
  
  out <- IntegratedLearner:::.prepare_from_PCL(PCL_train = pcl)
  expect_identical(out$feature_table, pcl$feature_table)
  expect_identical(out$sample_metadata, pcl$sample_metadata)
  expect_identical(out$feature_metadata, pcl$feature_metadata)
  
  missing_component <- pcl[c("feature_table", "feature_metadata")]
  expect_error(
    IntegratedLearner:::.prepare_from_PCL(PCL_train = missing_component),
    "PCL_train is missing required components: sample_metadata",
    fixed = TRUE
  )
  
  bad_rows <- pcl
  rownames(bad_rows$feature_metadata)[1] <- "not_a_feature"
  bad_rows$feature_metadata$featureID[1] <- "not_a_feature"
  expect_error(
    IntegratedLearner:::.prepare_from_PCL(PCL_train = bad_rows),
    "Row names of PCL_train$feature_table must equal row names of PCL_train$feature_metadata.",
    fixed = TRUE
  )
})

test_that("prepare_from_PCL drops incomplete features with na.rm", {
  pcl <- make_toy_pcl(n_samples = 10, n_features = 8, seed = 3)
  pcl$feature_table[1, 1] <- NA_real_
  
  out <- IntegratedLearner:::.prepare_from_PCL(PCL_train = pcl, na.rm = TRUE)
  
  expect_equal(nrow(out$feature_table), nrow(pcl$feature_table) - 1L)
  expect_identical(rownames(out$feature_table), rownames(out$feature_metadata))
  expect_false(anyNA(out$feature_table))
})

test_that("prepare_from_PCL checks validation feature identity", {
  pcl_train <- make_toy_pcl(n_samples = 12, n_features = 8, seed = 4)
  pcl_valid <- make_toy_pcl(n_samples = 6, n_features = 8, seed = 5)
  
  rownames(pcl_valid$feature_table)[1] <- "different_feature"
  rownames(pcl_valid$feature_metadata)[1] <- "different_feature"
  pcl_valid$feature_metadata$featureID[1] <- "different_feature"
  
  expect_error(
    IntegratedLearner:::.prepare_from_PCL(PCL_train = pcl_train, PCL_valid = pcl_valid),
    "Validation feature_table must have identical rownames as training feature_table (same features, same order).",
    fixed = TRUE
  )
})

test_that("validate_IL_inputs enforces non-survival constraints", {
  pcl <- make_toy_pcl(n_samples = 14, n_features = 10, seed = 6)
  
  expect_invisible(
    IntegratedLearner:::.validate_IL_inputs(
      feature_table = pcl$feature_table,
      sample_metadata = pcl$sample_metadata,
      feature_metadata = pcl$feature_metadata,
      family_name = "binomial",
      is_survival = FALSE
    )
  )
  
  non_finite <- pcl
  non_finite$sample_metadata$Y[1] <- Inf
  expect_error(
    IntegratedLearner:::.validate_IL_inputs(
      feature_table = non_finite$feature_table,
      sample_metadata = non_finite$sample_metadata,
      feature_metadata = non_finite$feature_metadata,
      family_name = "binomial",
      is_survival = FALSE
    ),
    "sample_metadata$Y contains non-finite values; please clean or impute the outcome.",
    fixed = TRUE
  )
  
  multiclass <- pcl
  multiclass$sample_metadata$Y <- rep(c(0, 1, 2), length.out = nrow(multiclass$sample_metadata))
  expect_error(
    IntegratedLearner:::.validate_IL_inputs(
      feature_table = multiclass$feature_table,
      sample_metadata = multiclass$sample_metadata,
      feature_metadata = multiclass$feature_metadata,
      family_name = "binomial",
      is_survival = FALSE
    ),
    "Multiclass outcomes are not supported",
    fixed = TRUE
  )
  
  non01 <- pcl
  non01$sample_metadata$Y <- rep(c(1, 2), length.out = nrow(non01$sample_metadata))
  expect_warning(
    IntegratedLearner:::.validate_IL_inputs(
      feature_table = non01$feature_table,
      sample_metadata = non01$sample_metadata,
      feature_metadata = non01$feature_metadata,
      family_name = "binomial",
      is_survival = FALSE
    ),
    "Binomial outcome is not coded as {0,1}",
    fixed = TRUE
  )
  
  low_unique <- pcl
  low_unique$sample_metadata$Y <- rep_len(1:5, nrow(low_unique$sample_metadata))
  expect_warning(
    IntegratedLearner:::.validate_IL_inputs(
      feature_table = low_unique$feature_table,
      sample_metadata = low_unique$sample_metadata,
      feature_metadata = low_unique$feature_metadata,
      family_name = "gaussian",
      is_survival = FALSE
    ),
    "Gaussian family with <=5 unique Y values - are you sure?",
    fixed = TRUE
  )
})

test_that("validate_IL_inputs enforces survival constraints", {
  pcl <- make_toy_pcl(n_samples = 16, n_features = 10, seed = 7)
  sample_ids <- rownames(pcl$sample_metadata)
  surv_meta <- make_toy_survival_metadata(sample_ids = sample_ids, seed = 77)
  
  expect_invisible(
    IntegratedLearner:::.validate_IL_inputs(
      feature_table = pcl$feature_table,
      sample_metadata = surv_meta,
      feature_metadata = pcl$feature_metadata,
      family_name = "survival",
      is_survival = TRUE
    )
  )
  
  bad_time <- surv_meta
  bad_time$time[1] <- -1
  expect_error(
    IntegratedLearner:::.validate_IL_inputs(
      feature_table = pcl$feature_table,
      sample_metadata = bad_time,
      feature_metadata = pcl$feature_metadata,
      family_name = "survival",
      is_survival = TRUE
    ),
    "Survival time must be non-negative.",
    fixed = TRUE
  )
  
  bad_event <- surv_meta
  bad_event$event[1] <- 2
  expect_error(
    IntegratedLearner:::.validate_IL_inputs(
      feature_table = pcl$feature_table,
      sample_metadata = bad_event,
      feature_metadata = pcl$feature_metadata,
      family_name = "survival",
      is_survival = TRUE
    ),
    "Survival event indicator must be in {0,1}.",
    fixed = TRUE
  )
})

test_that("ensure_survival_metadata derives time/event from Y", {
  sample_ids <- sprintf("S%03d", 1:8)
  times <- seq_along(sample_ids)
  events <- rep(c(0, 1), length.out = length(sample_ids))
  
  from_surv <- data.frame(
    subjectID = sample_ids,
    row.names = sample_ids,
    stringsAsFactors = FALSE
  )
  from_surv$Y <- survival::Surv(times, events)
  
  out_surv <- IntegratedLearner:::.ensure_survival_metadata(from_surv, context = "training")
  expect_equal(out_surv$time, as.numeric(times))
  expect_equal(out_surv$event, as.numeric(events))
  
  from_matrix <- data.frame(
    subjectID = sample_ids,
    row.names = sample_ids,
    stringsAsFactors = FALSE
  )
  from_matrix$Y <- I(cbind(time = times, event = events))
  
  out_matrix <- IntegratedLearner:::.ensure_survival_metadata(from_matrix, context = "training")
  expect_equal(out_matrix$time, as.numeric(times))
  expect_equal(out_matrix$event, as.numeric(events))
  
  no_surv_cols <- data.frame(
    subjectID = sample_ids,
    row.names = sample_ids,
    stringsAsFactors = FALSE
  )
  
  expect_error(
    IntegratedLearner:::.ensure_survival_metadata(no_surv_cols, context = "training"),
    "Survival mode requires 'time' and 'event' columns in sample_metadata (training).",
    fixed = TRUE
  )
})
