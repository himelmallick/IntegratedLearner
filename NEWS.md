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
