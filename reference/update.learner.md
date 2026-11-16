# Update IntegratedLearner fit object based on layers available in the test set

Allow update of IntegratedLearner if only a subset of omics layers are
available in test set. If all layers and features match, it calls
predict.learner()

## Usage

``` r
# S3 method for class 'learner'
update(
  object,
  feature_table_valid,
  sample_metadata_valid = NULL,
  feature_metadata_valid,
  outcome_col = NULL,
  subject_id_col = NULL,
  ...
)
```

## Arguments

- object:

  fitted 'IntegratedLearner' object

- feature_table_valid:

  Feature table from validation set. It should be a data frame with
  features in rows and samples in columns. Feature names should be a
  subset of training data feature names.

- sample_metadata_valid:

  OPTIONAL (can provide feature_table_valid and not this):
  Sample-specific metadata table from independent validation set. If
  provided, it must have the exact same structure as sample_metadata.
  Default is NULL.

- feature_metadata_valid:

  Matrix containing feature names and their corresponding layers. Must
  be subset of feature_metadata provided in IntegratedLearner object.

- outcome_col:

  Optional outcome column name in `sample_metadata_valid`. If `NULL`,
  uses `object$column_map$outcome_col` (or `"Y"`).

- subject_id_col:

  Optional subject ID column name in `sample_metadata_valid`. If `NULL`,
  uses `object$column_map$subject_id_col` (or `"subjectID"`).

- ...:

  Additional arguments (unused)

## Value

An updated IntegratedLearner prediction object with `yhat.test` and,
when validation outcomes are supplied, the observed outcomes and
corresponding performance summaries.

## Progress, summaries, and reproducibility

Call [`set.seed()`](https://rdrr.io/r/base/Random.html) before fitting
to reproduce the data splits and stochastic fitting steps. To display
progress, wrap the call in
[`progressr::with_progress()`](https://progressr.futureverse.org/reference/with_progress.html).
Use `print(fit)` to show a fitted model summary.

## Examples

``` r
set.seed(1)
sample_ids <- paste0("S", seq_len(18))
layer1 <- matrix(rnorm(3 * length(sample_ids)),
  nrow = 3,
  dimnames = list(paste0("L1_F", 1:3), sample_ids)
)
layer2 <- matrix(rnorm(2 * length(sample_ids)),
  nrow = 2,
  dimnames = list(paste0("L2_F", 1:2), sample_ids)
)
signal <- colMeans(layer1)
sample_metadata_valid <- data.frame(
  Y = as.numeric(signal > stats::median(signal)),
  subjectID = sample_ids,
  row.names = sample_ids
)
feature_table_valid <- as.data.frame(rbind(layer1, layer2))
feature_metadata_valid <- data.frame(
  featureID = rownames(feature_table_valid),
  featureType = c(rep("Layer1", 3), rep("Layer2", 2)),
  row.names = rownames(feature_table_valid)
)
fit <- IntegratedLearner(
  PCL_train = list(
    feature_table = feature_table_valid,
    sample_metadata = sample_metadata_valid,
    feature_metadata = feature_metadata_valid
  ),
  folds = 2, base_learner = "SL.mean",
  run_stacked = FALSE, run_concat = FALSE,
  family = stats::binomial()
)
upd <- update(
  object = fit,
  feature_table_valid = feature_table_valid,
  sample_metadata_valid = sample_metadata_valid,
  feature_metadata_valid = feature_metadata_valid
)
names(upd)
#> [1] "Y_test"                 "yhat.test"              "AUC.test"              
#> [4] "accuracy.test"          "balanced_accuracy.test" "metrics.test"          
```
