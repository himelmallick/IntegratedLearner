# IntegratedLearner 0.99.5

## Multiclass stacking

- Grouped multiclass results under survival-style `train_out` and `valid_out`
  entries for `single`, `early`, and all enabled `late` methods. Each method
  includes probabilities, classes, and metrics, with final and fold-specific
  weights for learned averaging. Prediction returns the same validation layout;
  existing flat outputs are preserved. Added `late_methods` to enumerate outputs.
- Added majority voting (`majority`) and an out-of-fold log-loss-weighted
  probability average (`weighted_average`) to multiclass fitting, prediction,
  and metrics. Both default to enabled and can be disabled independently with
  `run_majority` and `run_weighted`. Majority voting returns vote fractions and
  resolves ties by training class order. Final averaging weights are stored in
  `model_fits$model_weighted$weights`; training predictions cross-fit these
  weights across the existing folds.
- Fit the final multiclass meta-model on out-of-fold layer probabilities,
  matching the inputs used to train the cross-validated stack and reducing
  overfitting from in-sample base predictions.
- Updated the README, vignette, and fitting, prediction, and plotting help for
  all multiclass fusion methods, grouped results, and evaluation caveats.
- Clarified the existing multiclass class-association importance calculation
  and per-layer grouping. Corrected documentation to describe feature-by-class
  matrices and an unsigned global ranking; the calculation is unchanged.

# IntegratedLearner 0.99.4

## Tutorial data and documentation

- Removed the `data/` directory, raw tables, duplicate PCL datasets, and their
  obsolete dataset help pages. Seven MAE-only RDS fixtures in `inst/extdata`
  total 257,684 bytes, down from 13,720,424 bytes of packaged datasets.
- Each assay retains at most 25 finite, variable features selected using only
  training assays in original feature order. Paired validation assays use the
  same feature IDs and order. All original samples, subject IDs, outcomes,
  sample maps, and retained assay values are preserved.
- Updated the README, vignette, examples, and tests to access the tutorial
  RDS files directly with `readRDS(system.file("extdata", ...))`, replacing the
  former `data()` calls without adding a package-specific loader function.
  PCL inputs remain supported without bundling duplicate PCL datasets.
- Added a deterministic fixture-generation script and an inventory with source
  provenance and regeneration instructions in `inst/extdata/README.md`.
- Removed the outdated README project-status paragraph and documented the
  public entry points and naming conventions. Existing low-level exported
  names remain available for compatibility; package learner wrappers use `sl_*`.
- The default tutorial uses non-Java learners, including ranger with multinomial
  stacking for the small multiclass example. BART and its uncertainty examples
  remain optional extensions enabled with `IntegratedLearner.tutorial_bart`.

## Dependencies and S3 methods

- Registered `update.learner` as an S3 method and imported the `stats::update`
  generic. Updated its example and tests to call `update(fit, ...)`.
- Declared `randomForestSRC` in Imports as a required survival backend and
  replaced its runtime availability guard and dynamic namespace lookup with
  a direct `randomForestSRC::rfsrc()` call.
- Removed redundant availability checks for packages declared in Imports and
  the unused dependency-check helpers.
- Consolidated global-variable declarations in `R/globals.R` and removed the
  obsolete R-version check.

## Readability and shared implementation

- Documented why `make_sl_env()` provides package-local learner lookup for
  SuperLearner without attaching the package or modifying upstream namespaces.
- Removed `set_seed_internal()` and called `withr::local_seed()` directly.
  Centralized and documented the historical cooperative seed offset.
- Consolidated final glmnet fitting across `sl_glmnet2()`, `sl_lasso()`, and
  `sl_enet()`, preserving the elastic-net alpha search. Shared roxygen sections
  are inherited; each learner retains its own example.
- Shared cooperative matrix construction, model tuning, and result metadata
  across outcome families while retaining their distinct responses, folds,
  filtering, metrics, and output structures. Reused computed feature importance.
- Reused `slice_feature_table()` for the six cited multiclass transpositions
  and shared numeric screener-option parsing.
- Consolidated prediction feature checks with exact checks of feature names,
  lengths, and order to avoid accepting mismatches through vector recycling.
- Removed redundant numeric conversions, multiplication by 1, and no-op null
  fallbacks. Numeric percentage results may use integer or double storage.
- Removed unused `extract_multiclass_screen_args()`, `map_multiclass_screener()`,
  and `ibs_measure()` helpers and tests that only exercised deleted helpers.
- Flattened the nesting identified by the reviewer, simplified fragmented
  comment banners, and clarified section comments without splitting source files.
- Reused MAE preparation in tests and removed the legacy survival-fixture
  conversion/fallback code and redundant fixture-data loader.

## Performance and validation

- Replaced NNLS transpose/matrix products with `crossprod()`.
- Shared and vectorized IPCW Brier calculations. IBS optimization now caches
  its layer design matrix and censoring weights rather than constructing a
  three-dimensional array at every evaluation, and uses stable softmax weights.
- Added regression coverage for RNG restoration and reproducible folds,
  glmnet-wrapper predictions, exact feature alignment, fixture size and
  structure, vectorized Brier/IBS equivalence, and randomForestSRC fitting and
  prediction using the required backend.
- Verified retained MAE contents against their source objects. Executed the
  complete default tutorial across binary, continuous, multiclass, and survival
  workflows without warnings or errors; optional BART extensions were excluded.
- The revised implementation passes the installed-package test suite (627
  assertions, no failures or skips) and the requested unnecessary-nesting
  linter. 

# IntegratedLearner 0.99.3

## Cooperative multiview learning

- Added `run_intermediate` as a standalone cooperative multiview learning mode
  powered by `multiview`. In this mode, continuous/binary and survival fits use
  direct feature-level cooperative learning and skip base learners, early
  fusion, late fusion, and supervised screening.
- Multiclass outcomes do not support cooperative multiview learning; if
  `run_intermediate = TRUE`, a message is printed and the option is ignored.
- Added cooperative tuning controls `cooperative_rho`, `cooperative_s`, and
  `cooperative_type_measure`. `rho` is selected across repeated
  `cv.multiview()` fits with shared folds, while lambda selection is handled by
  `cv.multiview()`.
- Added coefficient-based cooperative feature-importance outputs, including
  absolute and signed global scores plus per-layer score lists for
  continuous/binary and survival cooperative fits.
- Updated the README, vignette, and Rd files to document standalone
  intermediate/cooperative mode and its output structure.
- Added tests for gaussian, binomial, survival, and multiclass
  `run_intermediate` behavior, including checks that cooperative mode does not
  depend on base or meta learners.

# IntegratedLearner 0.99.2

## Review-driven updates

- Revised package metadata for the current Bioconductor review cycle.
- Added Himel Mallick's ORCID to the `Authors@R` field in `DESCRIPTION`.
- Replaced the previous unformatted news file with a versioned `NEWS.md`.

## Bioconductor workflow and documentation

- Reworked the package documentation and vignette to make
  `MultiAssayExperiment` the primary user-facing workflow.
- Expanded the vignette introduction to describe how
  `IntegratedLearner` fits into a typical Bioconductor pipeline.
- Updated examples and vignette code to use standard Bioconductor containers
  and accessors such as `MultiAssayExperiment`, `SummarizedExperiment`,
  `experiments()`, `assay()`, `colData()`, and `sampleMap()`.
- Added package-level documentation and clarified return values for exported
  methods, including `predict.learner()`.
- Revised man page examples to be runnable under package checking and removed
  blank `if (FALSE)` example blocks.
- Removed disabled vignette install chunks and switched installation snippets
  to plain code blocks to reduce the risk of documentation drift.

## Dependencies and installation

- Reorganized dependencies so that packages needed at runtime are listed in
  `Imports`, while vignette/example/test-only packages are kept in `Suggests`
  where appropriate.
- Removed the runtime package-code dependence on `stringr` and replaced those
  uses with base R equivalents.
- Updated README installation instructions to include Bioconductor usage.

## API and naming cleanup

- Continued aligning package-specific exported helper/learner names with
  `snake_case` conventions and reserving `.` for S3 methods.
- Preserved backward compatibility for older learner identifiers while
  supporting the newer naming style.

## Examples, testing, and package checks

- Updated exported examples for `IntegratedLearner`, `predict.learner()`,
  `update.learner()`, and `plot.learner()` to improve reproducibility under
  `R CMD check`.
- Removed rendered vignette artifacts from `vignettes/` and added ignore rules
  to avoid recommitting generated HTML outputs.
- Expanded and refreshed tests around model methods, MAE workflows,
  prediction/update behavior, and utility helpers.

## Internal code quality

- Replaced the previous internal seed handling approach with a cleaner
  `withr`-based implementation.
- Updated `print.learner()` to return invisibly, consistent with base `print()`
  conventions.
- Reduced code repetition in the continuous/binary training, prediction, and
  update paths by introducing shared helper utilities for validation-layer
  prediction assembly, fusion prediction assembly, and family-specific metric
  calculation.
