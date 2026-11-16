# Contributing to IntegratedLearner

Thank you for your interest in contributing. `IntegratedLearner` is an
open source R package for multi-omics prediction and classification, and
we welcome contributions from the community — bug reports, documentation
improvements, new learners, and code.

By participating in this project, you agree to abide by our [Code of
Conduct](https://himelmallick.github.io/IntegratedLearner/CODE_OF_CONDUCT.md).

## Ways to contribute

- **Report a bug** — open an
  [issue](https://github.com/himelmallick/IntegratedLearner/issues)
- **Request a feature** — open an issue describing the use case
- **Improve documentation** — vignettes, roxygen docs, and examples
- **Contribute code** — bug fixes, new learners, performance
  improvements

## Reporting bugs

Before opening an issue, please search [existing
issues](https://github.com/himelmallick/IntegratedLearner/issues) to
avoid duplicates.

A good bug report includes:

1.  A **minimal reproducible example** (a
    [reprex](https://reprex.tidyverse.org/) is ideal) — the smallest
    code that reproduces the problem.
2.  The **expected** behavior and the **actual** behavior, including the
    full error message.
3.  Your session information:
    [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html) or
    [`sessioninfo::session_info()`](https://sessioninfo.r-lib.org/reference/session_info.html).
4.  Whether the issue involves the `PCL` or `MAE` input mode, and which
    outcome type (continuous, binary, multiclass, or survival).

Note that BART-based workflows (`SL.BART`, `bartMachine`) require a
working Java installation. If you hit a Java-related error, please
confirm Java is configured before filing.

## Development setup

\
[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``"devtools"``)`\
`devtools``::`[`install_github`](https://devtools.r-lib.org/reference/install-deprecated.html)`(``"himelmallick/IntegratedLearner"``, dependencies ``=`` ``TRUE``)`

To work on the package locally:

\
`# from the repository root`\
`devtools``::`[`load_all`](https://devtools.r-lib.org/reference/load_all.html)`(``)``     ``# load the package`\
`devtools``::`[`document`](https://devtools.r-lib.org/reference/document.html)`(``)``     ``# regenerate roxygen docs and NAMESPACE`\
`devtools``::`[`test`](https://devtools.r-lib.org/reference/test.html)`(``)``         ``# run the testthat suite`\
`devtools``::`[`check`](https://devtools.r-lib.org/reference/check.html)`(``)``        ``# full R CMD check`

`IntegratedLearner` targets Bioconductor, so please also run:

\
`BiocCheck``::``BiocCheck``(``)`

## Pull request process

1.  **Open an issue first** for anything beyond a trivial fix, so we can
    agree on the approach before you invest time.
2.  **Fork** the repository and create a branch from `master`:
    `git checkout -b feature/short-description`
3.  **Make your changes**, following the coding standards below.
4.  **Add tests** covering the new behavior (`tests/testthat/`).
    Contributions that change behavior should not reduce test coverage.
5.  **Update documentation** — roxygen comments, and the vignette if
    user-facing behavior changes. Run
    [`devtools::document()`](https://devtools.r-lib.org/reference/document.html).
6.  **Ensure
    [`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
    and `BiocCheck::BiocCheck()` pass** with no new errors, warnings, or
    notes.
7.  **Open a pull request** against `master`, describing what changed
    and why, and referencing the related issue.

A maintainer will review your PR. We aim to respond to pull requests and
issues promptly, and we may ask for changes before merging.

## Coding standards

- Follow the [Bioconductor coding
  style](https://contributions.bioconductor.org/r-code.html): 4-space
  indentation, lines under 80 characters, and `camelCase` or
  `snake_case` used consistently with the surrounding code.
- Document all exported functions with roxygen2, including `@param`,
  `@return`, and a runnable `@examples` block.
- Use `MultiAssayExperiment` / `SummarizedExperiment` classes for
  multi-omics containers rather than introducing new bespoke structures
  — interoperability with the Bioconductor ecosystem is a core design
  goal.
- Avoid adding heavy dependencies. New packages should go in `Suggests`
  unless they are essential to core functionality.
- Do not commit large data files. Tutorial data belongs in
  `inst/extdata/` as small MAE RDS files; do not duplicate it as raw
  tables or PCL objects.

## Finding your way around `R/`

`IntegratedLearner.R` prepares the inputs and dispatches by outcome.
`IL_conbin.R` fits continuous and binary outcomes, `IL_multiclass.R`
fits multiclass outcomes, `IL_survival.R` fits time-to-event outcomes,
and `IL_cooperative.R` implements cooperative learning.

| If you are changing… | Start in… |
|----|----|
| MAE or PCL input preparation, outcome conversion, or input checks | `IL_PrepareInputs.R`, `IL_Outcomes.R`, or `IL_validate.R` |
| Feature selection or filtering | `IL_Features.R` |
| SuperLearner adapters | `IL_SuperLearner.R` |
| Fold setup or fit controls | `IL_FitSupport.R` |
| Shared fusion results, metrics, or layer prediction | `IL_Fusion.R` or `IL_PredictLayers.R` |
| Numeric parsing shared by outcomes and feature options | `IL_Numeric.R` |
| Fit summaries | `IL_Print.R` |

## Adding a new learner

New learners are a welcome contribution. Please:

- For binary/continuous outcomes, follow the `SuperLearner` `SL.*`
  convention.
- For multiclass outcomes, register the learner with the native
  multiclass backend (see existing learners such as `glmnet`, `ranger`,
  `xgboost`).
- For survival outcomes, follow the `surv.*` convention used by the
  `ILsurv` engine.
- Include tests demonstrating the learner runs end-to-end in at least
  one fusion mode, and document it in the README’s supported-model list.

## Attribution

Contributors are credited in the package `DESCRIPTION` (`Authors@R`) and
in the GitHub contributor list. Substantial contributions may warrant
co-authorship on resulting publications — please raise this with the
maintainers.

## Questions

Open an
[issue](https://github.com/himelmallick/IntegratedLearner/issues) or
contact the maintainers:

- Nalin Arora — <naa4050@med.cornell.edu>
- Himel Mallick — <him4004@med.cornell.edu>
