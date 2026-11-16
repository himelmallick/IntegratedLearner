# Integrated machine learning for multi-omics prediction (continuous/binary outcomes)

Performs integrated machine learning to predict a **binary or
continuous** outcome based on two or more omics layers (views). This
function implements the core IntegratedLearner engine for **non-survival
outcomes**.

## Usage

``` r
IL_conbin(
  feature_table,
  sample_metadata,
  feature_metadata,
  feature_table_valid = NULL,
  sample_metadata_valid = NULL,
  folds = 5,
  base_learner = "sl_bart",
  base_screener = "All",
  run_screening = FALSE,
  screen_pct = NULL,
  filter_method = NULL,
  filter_pct = NULL,
  prevalence_pct = NULL,
  meta_learner = "sl_nnls_auc",
  run_concat = TRUE,
  run_stacked = TRUE,
  run_intermediate = FALSE,
  cooperative_rho = c(0, 0.1, 0.25, 0.5, 1),
  cooperative_s = "lambda.min",
  cooperative_type_measure = NULL,
  drop_poor_performing_layers = FALSE,
  refit.stack = FALSE,
  family = stats::gaussian(),
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

  How many folds in the V-fold nested cross-validation? Default is 5.

- base_learner:

  Base learner for late fusion and early fusion. Default: `sl_bart`.

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

  Meta-learner for late fusion (stacked generalization). Defaults to
  `sl_nnls_auc`.

- run_concat:

  Should early fusion be run? Default is TRUE.

- run_stacked:

  Should stacked model (late fusion) be run? Default is TRUE.

- run_intermediate:

  Should standalone cooperative learning via multiview be run? Default
  is FALSE. When TRUE, only the cooperative multiview model is fit; base
  learners, early fusion, and late fusion are skipped.

- cooperative_rho:

  Non-negative agreement-penalty values to evaluate for cooperative
  learning.

- cooperative_s:

  Penalty value used when predicting from the selected cooperative
  model. Defaults to `"lambda.min"`.

- cooperative_type_measure:

  Optional multiview CV loss/score. If `NULL`, an outcome-appropriate
  default is used.

- drop_poor_performing_layers:

  Logical; if `TRUE`, layers with single-layer performance below the
  screening threshold are removed from early and late fusion only. The
  threshold is AUC \< 0.5 for binary and R\\(^2\\) \< 0.5 for continuous
  outcomes. Single-layer outputs are still retained and reported.

- refit.stack:

  logical; For late fusion, refit predictions on the entire data are
  returned if specified. Default FALSE.

- family:

  Allows [`gaussian()`](https://rdrr.io/r/stats/family.html) for
  continuous and [`binomial()`](https://rdrr.io/r/stats/family.html) for
  binary outcomes. Survival outcomes must be handled via IL_survival().

- ...:

  Additional arguments (currently unused).

## Value

A list-like IntegratedLearner object containing fitted layer-specific,
stacked, and concatenated models, cross-validated performance (AUC or
R\302\262), and predictions for training and validation sets. When
`run_intermediate = TRUE`, coefficient-based cooperative
feature-importance scores from the selected multiview model are returned
in `cooperative_feature_importance`,
`cooperative_feature_importance_by_layer`,
`cooperative_feature_importance_signed`, and
`cooperative_feature_importance_signed_by_layer`.

## Details

`IL_conbin` takes a training set
`(feature_table, sample_metadata, feature_metadata)` and, optionally, a
corresponding validation set, and returns predicted values based on the
validation set. It also performs V-fold nested cross-validation to
estimate the prediction accuracy of various fusion algorithms.

Two integration paradigms are supported: early and late. The software
includes multiple ML models based on the
[`SuperLearner`](https://rdrr.io/pkg/SuperLearner/man/SuperLearner.html)
R package as well as several data exploration capabilities and
visualization modules in a unified estimation framework.

Although `IL_conbin()` is typically called internally by
[`IntegratedLearner()`](https://himelmallick.github.io/IntegratedLearner/reference/IntegratedLearner.md)
after extracting multi-view feature tables from `MultiAssayExperiment`
objects, advanced users may call `IL_conbin()` directly when they
already have multi-layer feature tables and metadata in
matrix/data.frame form.

## Progress, summaries, and reproducibility

Call [`set.seed()`](https://rdrr.io/r/base/Random.html) before fitting
to reproduce the data splits and stochastic fitting steps. To display
progress, wrap the call in
[`progressr::with_progress()`](https://progressr.futureverse.org/reference/with_progress.html).
Use `print(fit)` to show a fitted model summary.

## See also

[`IntegratedLearner`](https://himelmallick.github.io/IntegratedLearner/reference/IntegratedLearner.md),
IL_survival()

## Author

Himel Mallick, <him4004@med.cornell.edu>

## Examples

``` r
set.seed(1)
n <- 20
feature_table <- as.data.frame(rbind(
  matrix(rnorm(3 * n), nrow = 3, dimnames = list(paste0("L1_F", 1:3), paste0("S", 1:n))),
  matrix(rnorm(2 * n), nrow = 2, dimnames = list(paste0("L2_F", 1:2), paste0("S", 1:n)))
))
sample_metadata <- data.frame(
  subjectID = paste0("ID", 1:n), Y = rnorm(n),
  row.names = colnames(feature_table)
)
feature_metadata <- data.frame(
  featureID = rownames(feature_table),
  featureType = c(rep("Layer1", 3), rep("Layer2", 2)),
  row.names = rownames(feature_table)
)
fit <- IL_conbin(
  feature_table = feature_table,
  sample_metadata = sample_metadata,
  feature_metadata = feature_metadata,
  folds = 3, base_learner = "SL.mean", run_stacked = FALSE,
  run_concat = FALSE, family = stats::gaussian()
)
#> Loading required package: nnls
#> Warning: All algorithms have zero weight
#> Warning: All metalearner coefficients are zero, predictions will all be equal to 0
#> Warning: All algorithms have zero weight
#> Warning: All metalearner coefficients are zero, predictions will all be equal to 0
names(fit)
#>  [1] "model_fits"                         "SL_fits"                           
#>  [3] "X_train_layers"                     "Y_train"                           
#>  [5] "yhat.train"                         "base_learner"                      
#>  [7] "meta_learner"                       "base_screener"                     
#>  [9] "base_screener_used"                 "screening_used"                    
#> [11] "run_concat"                         "run_stacked"                       
#> [13] "run_intermediate"                   "cooperative_rho_grid"              
#> [15] "cooperative_s"                      "drop_poor_performing_layers"       
#> [17] "fusion_layers_retained"             "fusion_layers_removed"             
#> [19] "fusion_layer_scores"                "fusion_layer_metric"               
#> [21] "fusion_layer_threshold"             "family"                            
#> [23] "feature.names"                      "test"                              
#> [25] "R2.train"                           "feature_importance_signed"         
#> [27] "feature_importance_signed_by_layer" "folds"                             
#> [29] "cvControl"                          "time"                              
```
