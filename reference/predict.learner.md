# Make predictions using a trained 'IntegratedLearner' model

This function makes predictions using a trained 'IntegratedLearner'
model for new samples for which predictions are to be made

## Usage

``` r
# S3 method for class 'learner'
predict(
  object,
  feature_table_valid = NULL,
  sample_metadata_valid = NULL,
  feature_metadata = NULL,
  outcome_col = NULL,
  subject_id_col = NULL,
  ...
)
```

## Arguments

- object:

  Fitted 'IntegratedLearner' object

- feature_table_valid:

  Feature table from validation set. Must have the exact same structure
  as feature_table.

- sample_metadata_valid:

  OPTIONAL (can provide feature_table_valid and not this):
  Sample-specific metadata table from independent validation set. If
  provided, it must have the exact same structure as sample_metadata.

- feature_metadata:

  Matrix containing feature names and their corresponding layers. Must
  be same as that provided in IntegratedLearner object.

- outcome_col:

  Optional outcome column name in `sample_metadata_valid`. If `NULL`,
  uses the mapping stored in `object$column_map` (or falls back to
  `"Y"`).

- subject_id_col:

  Optional subject identifier column name in `sample_metadata_valid`. If
  `NULL`, uses the mapping stored in `object$column_map` (or falls back
  to `"subjectID"`).

- ...:

  Additional arguments (currently unused)

## Value

For gaussian and binomial fits, a list with at least

- `yhat.test`:

  A data frame or matrix of predictions for each single layer and, when
  requested during training, the `stacked` and `concatenated` fusion
  predictions.

- `Y_test`:

  Observed validation outcomes, returned only when
  `sample_metadata_valid` is supplied.

- `AUC.test`, `accuracy.test`, `balanced_accuracy.test`, `metrics.test`:

  Binary-outcome performance summaries, returned only when validation
  outcomes are supplied for binomial fits.

- `R2.test`:

  Named vector of validation-set R\\(^2\\) values, returned only when
  validation outcomes are supplied for gaussian fits.

For multiclass fits, returns the output of
`predict_multiclass.learner()`, which contains class probabilities,
predicted classes, and optional validation metrics. The flat outputs are
`prob.test`, `class.test`, and `yhat.test`, with `metrics.test` when
outcomes are available. Grouped results are returned in
`valid_out$single`, `valid_out$early`, and `valid_out$late`. The latter
includes all methods enabled during fitting: `stacked`, `majority`, and
`weighted_average`, as listed in `late_methods`. Each fusion method
includes `valid_prob`, `valid_class`, and a `metrics` table, plus
`valid_accuracy`, `valid_balanced_accuracy`, `valid_auc`, and
`valid_logloss` when labels are available. Without labels, `metrics` is
`NULL` and scalar metrics are absent. Learned averaging also includes
the training-derived layer `weights`; prediction does not refit them.
Majority outputs are vote fractions; ties use training class order, and
the fractions are not calibrated probabilities. Disabled early fusion
returns `early = NULL`, and disabled late methods are omitted from
`late`.

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
sample_metadata <- data.frame(
  Y = as.numeric(signal > stats::median(signal)),
  subjectID = sample_ids,
  row.names = sample_ids
)
feature_table <- as.data.frame(rbind(layer1, layer2))
feature_metadata <- data.frame(
  featureID = rownames(feature_table),
  featureType = c(rep("Layer1", 3), rep("Layer2", 2)),
  row.names = rownames(feature_table)
)
fit <- IntegratedLearner(
  PCL_train = list(
    feature_table = feature_table,
    sample_metadata = sample_metadata,
    feature_metadata = feature_metadata
  ),
  folds = 2, base_learner = "SL.mean",
  run_stacked = FALSE, run_concat = FALSE,
  family = stats::binomial()
)
pred <- predict(
  fit,
  feature_table_valid = feature_table,
  sample_metadata_valid = sample_metadata,
  feature_metadata = feature_metadata
)
names(pred)
#> [1] "Y_test"                 "yhat.test"              "AUC.test"              
#> [4] "accuracy.test"          "balanced_accuracy.test" "metrics.test"          
```
