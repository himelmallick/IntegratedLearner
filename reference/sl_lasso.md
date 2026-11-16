# Lasso regression with a fixed alpha

Uses the shared glmnet fitting implementation with `alpha = 1` by
default.

## Usage

``` r
sl_lasso(
  Y,
  X,
  newX,
  family,
  obsWeights,
  id,
  alpha = 1,
  nfolds = 10,
  nlambda = 100,
  useMin = TRUE,
  loss = "deviance",
  ...
)
```

## Arguments

- Y:

  Outcome variable

- X:

  Covariate dataframe

- newX:

  Dataframe to predict the outcome

- family:

  'gaussian' for regression, 'binomial' for binary classification.
  Untested options: 'multinomial' for multiple classification or
  'mgaussian' for multiple response, 'poisson' for non-negative outcome
  with proportional mean and variance, 'cox'.

- obsWeights:

  Optional observation-level weights

- id:

  Optional id to group observations from the same unit (not used
  currently).

- alpha:

  Elastic net mixing parameter, range 0 to 1. 0 = ridge regression and 1
  = lasso.

- nfolds:

  Number of folds for internal cross-validation to optimize lambda.

- nlambda:

  Number of lambda values to check, recommended to be 100 or more.

- useMin:

  If TRUE use lambda that minimizes risk, otherwise use 1 standard-error
  rule which chooses a higher penalty with performance within one
  standard error of the minimum (see Breiman et al. 1984 on CART for
  background).

- loss:

  Loss function, can be 'deviance', 'mse', or 'mae'. If family =
  binomial can also be 'auc' or 'class' (misclassification error).

- ...:

  Any additional arguments are passed through to cv.glmnet.

## Value

A list with elements `pred` (predictions for `newX`) and `fit`
(cross-validated glmnet fit metadata).

## References

Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization paths
for generalized linear models via coordinate descent. Journal of
statistical software, 33(1), 1.

Hoerl, A. E., & Kennard, R. W. (1970). Ridge regression: Biased
estimation for nonorthogonal problems. Technometrics, 12(1), 55-67.

Tibshirani, R. (1996). Regression shrinkage and selection via the lasso.
Journal of the Royal Statistical Society. Series B (Methodological),
267-288.

Zou, H., & Hastie, T. (2005). Regularization and variable selection via
the elastic net. Journal of the Royal Statistical Society: Series B
(Statistical Methodology), 67(2), 301-320.

## See also

`predict.SL.glmnet`
[`cv.glmnet`](https://glmnet.stanford.edu/reference/cv.glmnet.html)
[`glmnet`](https://glmnet.stanford.edu/reference/glmnet.html)

## Examples

``` r
set.seed(1)
X <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
linpred <- 0.5 * X$x1 - 0.25 * X$x2
Y <- stats::rbinom(nrow(X), 1, stats::plogis(linpred))
fit <- sl_lasso(
  Y = Y, X = X, newX = X[1:3, ], family = stats::binomial(),
  obsWeights = rep(1, nrow(X)), id = seq_len(nrow(X)),
  nfolds = 3, nlambda = 10
)
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
head(fit$pred)
#>   lambda.min
#> 1        0.5
#> 2        0.5
#> 3        0.5
```
