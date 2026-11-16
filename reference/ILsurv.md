# IntegratedLearner Survival Engine

BioC-friendly survival backend used by
[`IntegratedLearner()`](https://himelmallick.github.io/IntegratedLearner/reference/IntegratedLearner.md)
for time-to-event outcomes. This preserves the historical `ILsurv()`
interface while using model fitting and fusion routines that do not
depend on mlr3proba/mlr3extralearners.

## Usage

``` r
ILsurv(
  feature_table,
  sample_metadata,
  feature_metadata,
  valid_feature_table = NULL,
  valid_sample_metadata = NULL,
  base_learner = "surv.coxph",
  folds = 5,
  run_screening = FALSE,
  screen_pct = NULL,
  drop_poor_performing_layers = FALSE,
  do_early_fusion = TRUE,
  t_vec = NULL,
  t_vec_probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
  layer_score = c("sum", "mean", "l2"),
  eps = 1e-12,
  weight_lambda = 0.02,
  weight_penalty = c("l2_to_uniform", "entropy"),
  weight_cap = 1,
  optim_maxit_cox = 4000,
  optim_maxit_ibs = 300,
  ibs_shrink_to_uniform = 0,
  run_intermediate = FALSE,
  cooperative_rho = c(0, 0.1, 0.25, 0.5, 1),
  cooperative_s = "lambda.min",
  cooperative_type_measure = NULL,
  ...
)

IL_survival(
  feature_table,
  sample_metadata,
  feature_metadata,
  valid_feature_table = NULL,
  valid_sample_metadata = NULL,
  base_learner = "surv.coxph",
  folds = 5,
  run_screening = FALSE,
  screen_pct = NULL,
  drop_poor_performing_layers = FALSE,
  do_early_fusion = TRUE,
  t_vec = NULL,
  t_vec_probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
  layer_score = c("sum", "mean", "l2"),
  eps = 1e-12,
  weight_lambda = 0.02,
  weight_penalty = c("l2_to_uniform", "entropy"),
  weight_cap = 1,
  optim_maxit_cox = 4000,
  optim_maxit_ibs = 300,
  ibs_shrink_to_uniform = 0,
  run_intermediate = FALSE,
  cooperative_rho = c(0, 0.1, 0.25, 0.5, 1),
  cooperative_s = "lambda.min",
  cooperative_type_measure = NULL,
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

- valid_feature_table:

  Validation feature table (features x samples).

- valid_sample_metadata:

  Validation sample metadata containing `time` and `event`.

- base_learner:

  Survival base learner.

- folds:

  How many folds in the V-fold nested cross-validation? Default is 5.

- run_screening:

  Logical; if `TRUE`, run supervised screening within each training fold
  (and again on full training data) before fitting base models.

- screen_pct:

  Percentage of features to retain during screening (`(0,100]`). Applied
  after any optional filtering.

- drop_poor_performing_layers:

  Logical; if `TRUE`, layers with single-layer performance below the
  screening threshold are removed from early and late fusion only. The
  threshold is AUC \< 0.5 for binary and R\\(^2\\) \< 0.5 for continuous
  outcomes. Single-layer outputs are still retained and reported.

- do_early_fusion:

  Logical; run early fusion model.

- t_vec:

  Optional explicit time points for COX-based late-fusion feature
  construction.

- t_vec_probs:

  Quantiles used to derive `t_vec` when `t_vec` is `NULL`.

- layer_score:

  Layer summary on cumhaz increments: `'sum'`, `'mean'`, or `'l2'`.

- eps:

  Numerical lower bound for survival probabilities.

- weight_lambda:

  COX-weight optimization regularization strength.

- weight_penalty:

  COX-weight optimization regularizer: `'l2_to_uniform'` or `'entropy'`.

- weight_cap:

  Optional cap on any one layer weight (COX method).

- optim_maxit_cox:

  Maximum iterations for COX-weight optimization.

- optim_maxit_ibs:

  Maximum iterations for IBS-weight optimization.

- ibs_shrink_to_uniform:

  Optional convex shrinkage of IBS weights toward uniform.

- run_intermediate:

  Should standalone cooperative learning via multiview be run? Default
  is FALSE. When TRUE, only the cooperative multiview Cox model is fit;
  survival base learners, early fusion, and late fusion are skipped.

- cooperative_rho:

  Non-negative agreement-penalty values to evaluate for cooperative
  learning.

- cooperative_s:

  Penalty value used when predicting from the selected cooperative
  model. Defaults to `"lambda.min"`.

- cooperative_type_measure:

  Optional multiview CV loss/score. If `NULL`, Harrell's C is used.

- ...:

  Additional base-learner hyperparameters passed to `base_learner`. You
  may also pass `model_args = list(...)` as a named list where each
  entry is keyed by learner ID.

## Value

List with `train_out` and `valid_out`. Single-layer survival outputs
include `train_out$single$metrics`, `train_out$single$train_risk`, and
`train_out$single$feature_importance_signed_by_layer`, a named list of
unweighted signed feature-importance scores by layer. Late-fusion
results are returned under `train_out$late$IBS`, `train_out$late$COX`,
and the analogous validation entries. When `run_intermediate = TRUE`,
`train_out$cooperative` includes coefficient-based cooperative
feature-importance scores from the selected multiview Cox model.
`brier_time_grid` gives the shared training-derived grid used for
training and validation Brier/IBS evaluation. `survival_calibration`
stores the training-fitted risk-to-survival mappings. Single, early,
late, and cooperative entries include `train_surv`/`valid_surv`
probabilities when the corresponding model is fitted (lists for single
layers, matrices for fusion models).

## Details

Risk-to-survival Cox calibration is fitted using training outcomes and
OOF risk scores, then held fixed for validation, including coefficients,
baseline survival, missing-risk replacement, and the Brier time grid. A
training Kaplan-Meier curve is used if Cox calibration is not estimable.
Validation outcomes are used only for performance evaluation, including
estimating censoring weights, and never to fit a prediction calibration
model. IBS fusion is scored directly on its weighted survival curves;
single-layer, early-fusion, Cox-fusion, and cooperative models use their
frozen training calibration. Reported training probability-error metrics
remain internal summaries because calibration and fusion fitting use
training outcomes; they are not fully nested evaluation estimates.

## Progress, summaries, and reproducibility

Call [`set.seed()`](https://rdrr.io/r/base/Random.html) before fitting
to reproduce the data splits and stochastic fitting steps. To display
progress, wrap the call in
[`progressr::with_progress()`](https://progressr.futureverse.org/reference/with_progress.html).
Use `print(fit)` to show a fitted model summary.

## Examples

``` r
set.seed(1)
n <- 20
feature_table <- as.data.frame(rbind(
  matrix(rnorm(3 * n), nrow = 3, dimnames = list(paste0("L1_F", 1:3), paste0("S", 1:n))),
  matrix(rnorm(2 * n), nrow = 2, dimnames = list(paste0("L2_F", 1:2), paste0("S", 1:n)))
))
sample_metadata <- data.frame(
  subjectID = paste0("ID", 1:n),
  time = rexp(n, rate = 0.1),
  event = rbinom(n, 1, 0.6),
  row.names = colnames(feature_table)
)
sample_metadata$Y <- survival::Surv(sample_metadata$time, sample_metadata$event)
feature_metadata <- data.frame(
  featureID = rownames(feature_table),
  featureType = c(rep("Layer1", 3), rep("Layer2", 2)),
  row.names = rownames(feature_table)
)
fit <- ILsurv(
  feature_table = feature_table,
  sample_metadata = sample_metadata,
  feature_metadata = feature_metadata,
  folds = 3
)
names(fit)
#>  [1] "family"                      "train_out"                  
#>  [3] "valid_out"                   "screening_used"             
#>  [5] "brier_time_grid"             "survival_calibration"       
#>  [7] "screen_method"               "screen_pct"                 
#>  [9] "surv_plot_data"              "backend"                    
#> [11] "base_learner"                "supported_learners"         
#> [13] "fold_id"                     "folds"                      
#> [15] "test"                        "drop_poor_performing_layers"
#> [17] "late_methods"                "run_intermediate"           
#> [19] "cooperative_rho"             "cooperative_rho_grid"       
#> [21] "cooperative_rho_scores"      "cooperative_s"              
#> [23] "cooperative_type_measure"    "fusion_layers_retained"     
#> [25] "fusion_layers_removed"       "fusion_layer_scores"        
#> [27] "fusion_layer_metric"         "fusion_layer_threshold"     
fit$brier_time_grid
#>  [1]  0.2025041  0.3715227  0.6941235  0.9110887  1.2556587  1.7856945
#>  [7]  2.1369196  2.3565563  2.9915892  4.3587055  7.1861760  7.4304362
#> [13]  7.6318001  9.8181941 13.6457221 14.1954920 15.7271671 17.0992583
#> [19] 18.3564098 32.5933710
```
