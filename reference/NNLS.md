# NNLS function to optimize weights of several base learners

NNLS function to optimize weights of several base learners

## Usage

``` r
NNLS(x, y, wt)
```

## Arguments

- x:

  x

- y:

  y

- wt:

  wt

## Value

Solution of the quadratic programming problem

## Examples

``` r
x <- cbind(c(0.1, 0.4, 0.6, 0.9), c(0.2, 0.5, 0.7, 0.8))
y <- c(0.15, 0.45, 0.65, 0.85)
fit <- NNLS(x = x, y = y, wt = rep(1, nrow(x)))
fit$solution
#> [1] 0.5 0.5
```
