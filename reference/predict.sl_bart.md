# Predict function for sl_bart

Predict function for sl_bart

## Usage

``` r
# S3 method for class 'sl_bart'
predict(object, newdata, family, X = NULL, Y = NULL, ...)
```

## Arguments

- object:

  Fitted sl_bart model object

- newdata:

  Data frame for prediction

- family:

  Family object passed through (unused)

- X:

  Training design matrix (unused)

- Y:

  Training outcome (unused)

- ...:

  Additional arguments (unused)

## Value

Prediction from the sl_bart model

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
  head(predict(fit$fit, newdata = X[1:3, ], family = stats::gaussian()))
}
```
