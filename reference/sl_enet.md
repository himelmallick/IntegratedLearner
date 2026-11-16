# Elastic net regression with optimized alpha and lambda

Selects alpha by cross-validation, then uses the shared glmnet fitting
implementation.

## Usage

``` r
sl_enet(
  Y,
  X,
  newX,
  family,
  obsWeights,
  id,
  alpha = seq(0, 1, 0.1),
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

  Numeric grid of elastic net mixing parameters between 0 and 1.

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
Y <- rnorm(20)
fit <- sl_enet(
  Y = Y, X = X, newX = X[1:3, ], family = stats::gaussian(),
  obsWeights = rep(1, nrow(X)), id = seq_len(nrow(X)),
  alpha = c(0, 0.5, 1), nfolds = 3, nlambda = 10
)
head(fit$pred)
#>   lambda.min
#> 1 -0.3440358
#> 2 -0.1192026
#> 3 -0.2179840
```
