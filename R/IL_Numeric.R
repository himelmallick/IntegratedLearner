# Parse decimal numbers supplied as values, strings, or factor labels.
parse_numeric_values <- function(x) {
  x <- trimws(as.character(x))
  decimal <- grepl("^[+-]?(?:[0-9]+[.]?[0-9]*|[.][0-9]+)(?:[eE][+-]?[0-9]+)?$", x)
  out <- rep(NA_real_, length(x))
  out[decimal] <- as.numeric(x[decimal])
  out
}

parse_numeric_scalar <- function(x) {
  if (length(x) == 0L) NA_real_ else parse_numeric_values(x[[1L]])
}

expit <- function(x) {
  1 / (1 + exp(-x))
}
