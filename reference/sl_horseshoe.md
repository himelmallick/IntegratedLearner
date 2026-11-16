# Horseshoe regression

Horseshoe regression

## Usage

``` r
sl_horseshoe(
  Y,
  X,
  newX,
  family,
  prior = "horseshoe",
  N = 20000L,
  burnin = 1000L,
  thinning = 1L,
  ...
)
```

## Arguments

- Y:

  Outcome variable

- X:

  Covariate data frame

- newX:

  Dataframe to predict the outcome

- family:

  'gaussian' for regression, 'binomial' for binary classification.
  Untested options: 'poisson' for for integer or count data

- prior:

  prior for regression coefficients to use. 'Horseshoe' by default.
  Untested options: ridge regression (prior='rr' or prior='ridge'),
  lasso regression (prior='lasso') and horseshoe+ regression
  (prior='hs+' or prior='horseshoe+')

- N:

  Number of posterior samples to generate.

- burnin:

  Number of burn-in samples.

- thinning:

  Desired level of thinning.

- ...:

  other parameters passed to bayesreg function

## Value

A list with elements `pred` (predictions for `newX`) and `fit` (the
fitted Bayesian regression object).

## Examples

``` r
if (
  requireNamespace("bayesreg", quietly = TRUE) &&
  identical(Sys.getenv("INTEGRATEDLEARNER_RUN_BAYESREG_EXAMPLES"), "true")
) {
  set.seed(1)
  X <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
  Y <- rnorm(20)
  fit <- sl_horseshoe(
    Y = Y, X = X, newX = X[1:3, ],
    family = stats::gaussian()
  )
  fit$pred
}
```
