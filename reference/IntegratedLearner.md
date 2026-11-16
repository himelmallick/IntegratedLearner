# Integrated machine learning for multi-omics prediction and classification

Performs integrated machine learning to predict a binary, continuous, or
time-to-event outcome based on two or more omics layers (views). The
`IntegratedLearner` function takes a training `MultiAssayExperiment`
and, optionally, a validation `MultiAssayExperiment`, inspects standard
Bioconductor containers and accessors, extracts multi-layer feature
tables, and returns predicted values based on the validation set. It
also performs V-fold nested cross-validation to estimate the prediction
accuracy of various fusion algorithms. Two integration paradigms are
supported: early and late. The software includes multiple ML models
based on the
[`SuperLearner`](https://rdrr.io/pkg/SuperLearner/man/SuperLearner.html)
R package as well as several data exploration capabilities and
visualization modules in a unified estimation framework.

## Usage

``` r
IntegratedLearner(
  MAE_train = NULL,
  MAE_valid = NULL,
  PCL_train = NULL,
  PCL_valid = NULL,
  experiment = NULL,
  assay.type = NULL,
  outcome_col = "Y",
  subject_id_col = "subjectID",
  na.rm = FALSE,
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
  run_majority = TRUE,
  run_weighted = TRUE,
  ...
)
```

## Arguments

- MAE_train:

  A `MultiAssayExperiment` containing the training data. Each experiment
  corresponds to one view (omics layer), usually stored as a
  `SummarizedExperiment` or `TreeSummarizedExperiment`. The `colData` of
  every selected experiment must contain a column named by `outcome_col`
  describing the outcome of interest (binary, continuous, or
  survival-encoded), and a subject identifier column named by
  `subject_id_col`.

- MAE_valid:

  Optional `MultiAssayExperiment` containing an independent
  validation/test set for which prediction is desired. It must contain
  the same set of experiments used for training, with identical assay
  names and identical feature (row) names for each experiment.
  Additional experiments in `mae_test` are ignored.

- PCL_train:

  Optional list of per-layer feature matrices (legacy PCL mode for
  backward compatibility).

- PCL_valid:

  Optional validation list of per-layer feature matrices.

- experiment:

  Optional character or integer vector specifying which experiments
  (layers) to extract from `MAE_train` (and `MAE_valid`, if supplied).
  If `NULL`, all experiments in `MAE_train` are used.

- assay.type:

  Optional character vector of assay names, one per experiment. If
  `NULL`, each selected experiment must contain exactly one assay;
  otherwise, the user must supply `assay.type`.

- outcome_col:

  Outcome column name in MAE/PCL sample metadata. Defaults to `"Y"` for
  backward compatibility.

- subject_id_col:

  Subject identifier column name in MAE/PCL sample metadata. Defaults to
  `"subjectID"` for backward compatibility.

- na.rm:

  Logical; if `TRUE`, features with missing values are dropped after
  extraction.

- folds:

  How many folds in the V-fold nested cross-validation? Default is 5.

- base_learner:

  Base learner for late fusion and early fusion. Check out the
  [SuperLearner package
  page](https://CRAN.R-project.org/package=SuperLearner) for all
  available options. Default is `` `sl_bart` ``.

- base_screener:

  Deprecated for `IL_conbin` and `IL_multiclass`; kept for backward
  compatibility and currently ignored in those backends.

- run_screening:

  Logical; if `TRUE`, run supervised screening on training data within
  each CV fold (and again on full training data before final fitting),
  then apply the same selected features to fold-validation or external
  validation data.

- screen_pct:

  Percentage of features to retain during screening (`(0,100]`).

- filter_method:

  Optional feature-filter method for non-survival runs. Supported
  values: `'prevalence'` and `'variance'`.

- filter_pct:

  Optional retention percentage (in `(0,100]`) for the selected
  `filter_method`. Keeps the top `filter_pct` percent features.

- prevalence_pct:

  Optional retention percentage (in `(0,100]`) for prevalence-based
  feature filtering before model fitting. Deprecated alias of
  `filter_pct` with `filter_method = 'prevalence'`.

- meta_learner:

  Meta-learner for late fusion (stacked generalization). Defaults to
  `` `sl_nnls_auc` ``. For multiclass outcomes this legacy default maps
  to multinomial `glmnet`; use `run_weighted` for learned scalar layer
  weights. Check out the [SuperLearner package
  page](https://CRAN.R-project.org/package=SuperLearner) for all
  available options.

- run_concat:

  Should early fusion be run? Default is TRUE. Uses the specified
  `base_learner` as the learning algorithm.

- run_stacked:

  Should stacked model (late fusion) be run? Default is TRUE.

- run_intermediate:

  Should standalone cooperative learning via multiview be run? Default
  is FALSE. When TRUE, IntegratedLearner runs only the cooperative
  multiview model and skips base learners, early fusion, and late
  fusion. Supported for continuous, binary, and survival outcomes;
  ignored with a message for multiclass outcomes.

- cooperative_rho:

  Non-negative agreement-penalty values to evaluate for cooperative
  learning. A separate `cv.multiview()` fit is trained for each value,
  and the best value is selected by cross-validation.

- cooperative_s:

  Penalty value used when predicting from the selected cooperative
  model. Defaults to `"lambda.min"`.

- cooperative_type_measure:

  Optional multiview CV loss/score. If `NULL`, an outcome-appropriate
  default is used.

- drop_poor_performing_layers:

  Logical; if `TRUE`, layers with poor single-layer performance are
  removed from early and late fusion only. The thresholds are AUC \< 0.5
  for binary, R\\(^2\\) \< 0.5 for continuous, and C-index \< 0.5 for
  survival outcomes. Single-layer outputs are still retained and
  reported.

- refit.stack:

  logical; For late fusion, post-refit predictions on the entire data
  are returned if specified. Default is FALSE.

- family:

  Currently allows `` `gaussian()` `` for continuous, `` `binomial()` ``
  for binary, or an appropriate Cox/survival family for time-to-event
  outcomes. Determines whether `IL_conbin` or `IL_survival` is used
  internally.

- run_majority:

  Logical; return multiclass majority voting. Default TRUE; ignored for
  other outcomes. Ties use the first class in training class order.

- run_weighted:

  Logical; return multiclass probability averaging with nonnegative
  layer weights summing to one, learned by minimizing out-of-fold
  multiclass log loss. Default TRUE; ignored for other outcomes.

- ...:

  Additional arguments passed to the underlying engines.

## Value

A list-like IntegratedLearner object containing the trained model fits
(layer-specific, stacked, and/or concatenated models), cross-validated
performance estimates, and predicted values for training and, if
supplied, validation data. When `run_intermediate = TRUE`, the object
also includes coefficient-based cooperative feature-importance scores
from the selected multiview model. Multiclass fits also provide
`train_out` and `valid_out` grouped into `single`, `early`, and `late`
results. All enabled late methods (`stacked`, `majority`,
`weighted_average`) include probabilities, predicted classes, and
metrics; learned averaging also reports layer weights. Existing flat
multiclass outputs remain available.

## Details

In a typical Bioconductor workflow, users construct assay-level data in
`SummarizedExperiment`, combine those assays in `MultiAssayExperiment`,
inspect the object with `experiments()`, `assay()`, `colData()`, and
`sampleMap()`, and then pass the MAE into `IntegratedLearner()` for
supervised multi-omics prediction. Internally, `IntegratedLearner()`
converts the `MultiAssayExperiment` into tabular multi-view matrices and
then calls `IL_conbin` for continuous/binary outcomes or `IL_survival`
for time-to-event outcomes, depending on the specified `family`. Legacy
PCL-style inputs remain available for backward compatibility, but the
primary user-facing workflow is MAE-based.

## Progress, summaries, and reproducibility

Call [`set.seed()`](https://rdrr.io/r/base/Random.html) before fitting
to reproduce the data splits and stochastic fitting steps. To display
progress, wrap the call in
[`progressr::with_progress()`](https://progressr.futureverse.org/reference/with_progress.html).
Use `print(fit)` to show a fitted model summary.

## Author

Himel Mallick, <him4004@med.cornell.edu>

## Examples

``` r
PRISM_MAE <- readRDS(system.file(
  "extdata", "PRISM_MAE.rds", package = "IntegratedLearner", mustWork = TRUE
))
mae <- PRISM_MAE
names(MultiAssayExperiment::experiments(mae))
#> [1] "species"     "metabolites"
head(as.data.frame(S4Vectors::DataFrame(MultiAssayExperiment::colData(mae)))[, 1:4])
#>        Diagnosis dysbiosis_score Y subjectID
#> G35127        CD       0.9341207 1    G35127
#> G35128        CD       0.5962602 1    G35128
#> G35152        CD       0.9505732 1    G35152
#> G36347        CD       0.9966957 1    G36347
#> G36348        CD       0.8475403 1    G36348
#> G36349        CD       0.9569372 1    G36349
head(MultiAssayExperiment::sampleMap(mae))
#> DataFrame with 6 rows and 3 columns
#>      assay     primary     colname
#>   <factor> <character> <character>
#> 1  species      G35127      G35127
#> 2  species      G35128      G35128
#> 3  species      G35152      G35152
#> 4  species      G36347      G36347
#> 5  species      G36348      G36348
#> 6  species      G36349      G36349
se <- MultiAssayExperiment::experiments(mae)[[1]]
SummarizedExperiment::assayNames(se)
#> [1] "abundance"
fit <- IntegratedLearner(
  MAE_train = mae,
  experiment = names(MultiAssayExperiment::experiments(mae))[1:2],
  assay.type = rep("abundance", 2),
  folds = 2, base_learner = "SL.mean",
  run_stacked = FALSE, run_concat = FALSE,
  family = stats::binomial()
)
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
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
#> [25] "train_out"                          "AUC.train"                         
#> [27] "accuracy.train"                     "balanced_accuracy.train"           
#> [29] "metrics.train"                      "feature_importance_signed"         
#> [31] "feature_importance_signed_by_layer" "folds"                             
#> [33] "cvControl"                          "time"                              
#> [35] "input_mode"                         "column_map"                        
```
