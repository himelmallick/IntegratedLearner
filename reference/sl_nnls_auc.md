# Combined SuperLearner function for both NNLS/AUC maximization

Combined SuperLearner function for both NNLS/AUC maximization

## Usage

``` r
sl_nnls_auc(Y, X, newX, family, obsWeights, bounds = c(0, Inf), ...)
```

## Arguments

- Y:

  Outcome matrix from metalearner

- X:

  Layer-level predictions used to train the metalearner

- newX:

  Layer-level predictions for validation data

- family:

  Family object

- obsWeights:

  Observation weights

- bounds:

  Lower/upper bounds for weights (binomial case)

- ...:

  Additional arguments passed through

## Value

Estimated meta-learner coefficients and predictions

## Examples

``` r
set.seed(1)
X <- data.frame(m1 = runif(20), m2 = runif(20))
Y <- rnorm(20)
fit <- sl_nnls_auc(
  Y = Y, X = X, newX = X[1:4, ],
  family = stats::gaussian(), obsWeights = rep(1, nrow(X))
)
head(fit$pred)
#>        [,1]
#> 1 0.9347052
#> 2 0.2121425
#> 3 0.6516738
#> 4 0.1255551
```
