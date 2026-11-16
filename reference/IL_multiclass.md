# Integrated machine learning for multi-omics multiclass classification

Native multiclass backend used by
[`IntegratedLearner()`](https://himelmallick.github.io/IntegratedLearner/reference/IntegratedLearner.md)
when `family = binomial()` and the outcome has more than two classes.
This backend preserves the existing API while using multiclass-native
probability modeling and multiclass stacking.

## Usage

``` r
IL_multiclass(
  feature_table,
  sample_metadata,
  feature_metadata,
  feature_table_valid = NULL,
  sample_metadata_valid = NULL,
  folds = 5,
  base_learner = "glmnet",
  base_screener = "All",
  run_screening = FALSE,
  screen_pct = NULL,
  filter_method = NULL,
  filter_pct = NULL,
  prevalence_pct = NULL,
  meta_learner = "glmnet",
  run_concat = TRUE,
  run_stacked = TRUE,
  run_intermediate = FALSE,
  cooperative_rho = c(0, 0.1, 0.25, 0.5, 1),
  cooperative_s = "lambda.min",
  cooperative_type_measure = NULL,
  family = stats::binomial(),
  eps = 1e-15,
  run_majority = TRUE,
  run_weighted = TRUE,
  ...
)
```

## Arguments

- feature_table:

  An R data frame containing multi-layer features (in rows) and samples
  (in columns). Column names of `feature_metadata` must match the row
  names of `sample_metadata`.

- sample_metadata:

  An R data frame of metadata variables (in columns). Must have a column
  named `subjectID` describing per-subject unique identifiers. For
  longitudinal designs, this variable may be non-unique. Additionally, a
  column named `Y` must be present which is the continuous or binary
  outcome of interest. Row names of `sample_metadata` must match the
  column names of `feature_table`.

- feature_metadata:

  An R data frame of feature-specific metadata across views (in columns)
  and features (in rows). Must have a column named `featureID` as a
  unique per-feature identifier and a column named `featureType`
  describing the source layers. Row names of `feature_metadata` must
  match the row names of `feature_table`.

- feature_table_valid:

  Feature table from validation set for which prediction is desired.
  Must have the exact same structure as `feature_table`. If missing,
  uses `feature_table`.

- sample_metadata_valid:

  Sample-specific metadata table from independent validation set when
  available. Must have the exact same structure as `sample_metadata`.

- folds:

  Number of subject-level folds for out-of-fold predictions. Default 5.
  Fully nested evaluation of the fusion pipeline is separate.

- base_learner:

  Native multiclass base learner. Default `"glmnet"`.

- base_screener:

  Deprecated for this backend; currently ignored and kept only for
  backward compatibility.

- run_screening:

  Logical; if `TRUE`, run supervised screening within each training fold
  (and again on full training data) before fitting base models.

- screen_pct:

  Percentage of features to retain during screening (`(0,100]`). Applied
  after any optional filtering.

- filter_method:

  Optional feature-filter method before model fitting. Supported values
  are `'prevalence'` (top detected features by prevalence) and
  `'variance'` (top features by caret-based variance ranking via
  `nearZeroVar` + empirical variance). If `NULL`, defaults to
  `'prevalence'` when filtering is requested.

- filter_pct:

  Optional retention percentage (in `(0,100]`) for the selected
  `filter_method`. Keeps the top `filter_pct` percent features.

- prevalence_pct:

  Optional retention percentage (in `(0,100]`) for prevalence-based
  filtering before model fitting. Deprecated alias of `filter_pct` with
  `filter_method = 'prevalence'`.

- meta_learner:

  Native multiclass stacking learner. Default `"glmnet"`. Legacy
  `sl_nnls_auc` and `SL.nnls` aliases map to multinomial glmnet, not
  scalar layer-weight optimization.

- run_concat:

  Should early fusion be run? Default is TRUE.

- run_stacked:

  Should stacked model (late fusion) be run? Default is TRUE.

- run_intermediate:

  Ignored for multiclass outcomes with a message; cooperative multiview
  learning is not currently supported for multiclass fits.

- cooperative_rho:

  Non-negative agreement-penalty values to evaluate for cooperative
  learning.

- cooperative_s:

  Penalty value used when predicting from the selected cooperative
  model. Defaults to `"lambda.min"`.

- cooperative_type_measure:

  Optional multiview CV loss/score. If `NULL`, an outcome-appropriate
  default is used.

- family:

  Allows [`gaussian()`](https://rdrr.io/r/stats/family.html) for
  continuous and [`binomial()`](https://rdrr.io/r/stats/family.html) for
  binary outcomes. Survival outcomes must be handled via IL_survival().

- eps:

  Small positive constant used to stabilize probabilities.

- run_majority:

  Logical; return multiclass majority voting. Default TRUE; ignored for
  other outcomes. Ties use the first class in training class order.

- run_weighted:

  Logical; return multiclass probability averaging with nonnegative
  layer weights summing to one, learned by minimizing out-of-fold
  multiclass log loss. Default TRUE; ignored for other outcomes.

- ...:

  Additional arguments passed to native base-model fitting; these are
  not passed to the stacked meta-model or fusion-weight optimizer.

## Value

A fitted multiclass IntegratedLearner object containing class
probabilities, predicted classes, fitted single-layer and fusion models,
and training/validation multiclass performance summaries when available.

## Details

Majority voting returns vote fractions under `majority`; these are not
calibrated class probabilities. Ties use the first class in training
class order. Learned averaging returns probabilities under
`weighted_average` and stores final layer weights under
`model_fits$model_weighted$weights`. Weights minimize multiclass log
loss on out-of-fold layer probabilities. Training outputs cross-fit the
weights across the existing folds. As with stacking, fully nested
evaluation is needed for an unbiased assessment of the entire fitted
pipeline. Results are also grouped under `train_out` and `valid_out`,
each with `single`, `early`, and `late` entries. The `late` list
contains all enabled methods: `stacked`, `majority`, and
`weighted_average`. Each method includes split-specific probabilities,
classes, and metrics; learned averaging also includes final weights and,
for training, the weights fitted outside each fold. Existing flat
outputs remain available. `late_methods` lists the returned late-fusion
methods. Fusion entries contain `train_prob`, `train_class`, a `metrics`
table, and scalar `train_accuracy`, `train_balanced_accuracy`,
`train_auc`, and `train_logloss`; validation entries use the
corresponding `valid_` names. AUC is the macro average of available
one-vs-rest class AUCs. Single-layer entries contain probability lists
and class data frames. Without validation data, `valid_out` is `NULL`.
Disabled early fusion gives `early = NULL`, and disabled late methods
are omitted from `late`.

Feature importance retains a univariate class-association calculation:
for each feature and class, the mean in the remaining classes is
subtracted from the mean in that class and divided by the feature's
overall training standard deviation plus `1e-8`. The signed scores form
a feature-by-class matrix in `feature_importance_signed_by_class`.
`feature_importance_signed_by_layer_by_class` groups these matrices by
omics layer and ranks features by mean absolute score across classes.
`feature_importance_global` is the unsigned mean absolute score across
classes, ranked highest first. These associations use the training
features after optional filtering, including features excluded by
screening; they are shared across fusion methods and are not
model-specific attributions.

## Progress, summaries, and reproducibility

Call [`set.seed()`](https://rdrr.io/r/base/Random.html) before fitting
to reproduce the data splits and stochastic fitting steps. To display
progress, wrap the call in
[`progressr::with_progress()`](https://progressr.futureverse.org/reference/with_progress.html).
Use `print(fit)` to show a fitted model summary.

## Examples

``` r
set.seed(1)
n <- 24
feature_table <- as.data.frame(rbind(
  matrix(rnorm(3 * n), nrow = 3, dimnames = list(paste0("L1_F", 1:3), paste0("S", 1:n))),
  matrix(rnorm(2 * n), nrow = 2, dimnames = list(paste0("L2_F", 1:2), paste0("S", 1:n)))
))
y <- rep(c("A", "B", "C"), length.out = n)
sample_metadata <- data.frame(
  subjectID = paste0("ID", 1:n), Y = y,
  row.names = colnames(feature_table)
)
feature_metadata <- data.frame(
  featureID = rownames(feature_table),
  featureType = c(rep("Layer1", 3), rep("Layer2", 2)),
  row.names = rownames(feature_table)
)
fit <- IL_multiclass(
  feature_table = feature_table,
  sample_metadata = sample_metadata,
  feature_metadata = feature_metadata,
  folds = 3, base_learner = "multinom", meta_learner = "multinom",
  run_stacked = FALSE, run_concat = FALSE
)
fit$family
#> [1] "multinomial"
fit$late_methods
#> [1] "majority"         "weighted_average"
fit$train_out$late$weighted_average$weights
#>    Layer1    Layer2 
#> 0.2901878 0.7098122 
head(fit$feature_importance_signed_by_layer_by_class[[1]])
#>                A          B          C
#> L1_F1 -0.5194429 -0.1502983  0.6697412
#> L1_F3  0.5698738 -0.3140530 -0.2558207
#> L1_F2  0.3711250 -0.1626989 -0.2084261
```
