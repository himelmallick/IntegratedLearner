# Wrapper for bartMachine learner

Support bayesian additive regression trees via the bartMachine package.

## Usage

``` r
sl_bart(
  Y,
  X,
  newX,
  family,
  obsWeights,
  id,
  num_trees = 50,
  num_burn_in = 250,
  verbose = FALSE,
  alpha = 0.95,
  beta = 2,
  k = 2,
  q = 0.9,
  nu = 3,
  num_iterations_after_burn_in = 1000,
  serialize = TRUE,
  seed = 5678,
  ...
)
```

## Arguments

- Y:

  Outcome variable

- X:

  Covariate dataframe

- newX:

  Optional dataframe to predict the outcome

- family:

  'gaussian' for regression, 'binomial' for binary classification

- obsWeights:

  Optional observation-level weights (supported but not tested)

- id:

  Optional id to group observations from the same unit (not used
  currently).

- num_trees:

  The number of trees to be grown in the sum-of-trees model.

- num_burn_in:

  Number of MCMC samples to be discarded as 'burn-in'.

- verbose:

  Prints information about progress of the algorithm to the screen.

- alpha:

  Base hyperparameter in tree prior for whether a node is nonterminal or
  not.

- beta:

  Power hyperparameter in tree prior for whether a node is nonterminal
  or not.

- k:

  For regression, k determines the prior probability that E(Y\|X) is
  contained in the interval \\((y\_{min}, y\_{max})\\), based on a
  normal distribution. For example, when k=2, the prior probability is
  95 percent. For classification, k determines the prior probability
  that E(Y\|X) is between (-3,3). Note that a larger value of k results
  in more shrinkage and a more conservative fit.

- q:

  Quantile of the prior on the error variance at which the data-based
  estimate is placed. Note that the larger the value of q, the more
  aggressive the fit as you are placing more prior weight on values
  lower than the data-based estimate. Not used for classification.

- nu:

  Degrees of freedom for the inverse chi^2 prior. Not used for
  classification.

- num_iterations_after_burn_in:

  Number of MCMC samples to draw from the posterior distribution of
  f(x).

- serialize:

  If TRUE, bartMachine results can be saved to a file, but will require
  additional RAM.

- seed:

  Seed for reproducibility passed to bartMachine.

- ...:

  Additional arguments (not used)

## Value

A list with elements `pred` (predictions for `newX`) and `fit` (the
fitted model object).

## Examples

``` r
if (
  requireNamespace("bartMachine", quietly = TRUE) &&
  nzchar(Sys.which("java")) &&
  any(grepl("jdk.incubator.vector", getOption("java.parameters", character()), fixed = TRUE))
) {
  set.seed(1)
  X <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
  Y <- rnorm(20)
  fit <- sl_bart(
    Y = Y, X = X, newX = X[1:3, ],
    family = stats::gaussian(),
    obsWeights = rep(1, nrow(X)), id = seq_len(nrow(X)),
    num_trees = 5, num_burn_in = 5, num_iterations_after_burn_in = 20
  )
  fit$pred
}
```
