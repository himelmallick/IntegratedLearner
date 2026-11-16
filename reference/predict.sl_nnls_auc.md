# Predict function for sl_nnls_auc

Predict function for sl_nnls_auc

## Usage

``` r
# S3 method for class 'sl_nnls_auc'
predict(object, newdata, ...)
```

## Arguments

- object:

  Fitted sl_nnls_auc model

- newdata:

  Validation layer-level predictions

- ...:

  Additional arguments (unused)

## Value

Prediction from the meta-learner

## Examples

``` r
set.seed(1)
X <- data.frame(m1 = runif(20), m2 = runif(20))
Y <- rnorm(20)
fit <- sl_nnls_auc(
  Y = Y, X = X, newX = X[1:4, ],
  family = stats::gaussian(), obsWeights = rep(1, nrow(X))
)
head(predict(fit$fit, newdata = X[1:4, ]))
#>        [,1]
#> 1 0.9347052
#> 2 0.2121425
#> 3 0.6516738
#> 4 0.1255551
```
