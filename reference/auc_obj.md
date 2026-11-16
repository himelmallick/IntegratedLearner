# Title Meta level Objective function: NNLS for gaussian; Rank loss for binary observations

Title Meta level Objective function: NNLS for gaussian; Rank loss for
binary observations

## Usage

``` r
auc_obj(b, X, Y)
```

## Arguments

- b:

  Weights vector

- X:

  Design matrix (data frame)

- Y:

  Outcome variable

## Value

1 - AUC

## Examples

``` r
X <- data.frame(m1 = c(0.1, 0.4, 0.8, 0.9), m2 = c(0.2, 0.3, 0.7, 0.95))
Y <- c(0, 0, 1, 1)
auc_obj(b = c(0.5, 0.5), X = X, Y = Y)
#> [1] 0
```
