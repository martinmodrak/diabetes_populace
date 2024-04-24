phano <- function(x) {
  mean(x) / var(x)
}

prop_zero <- function(x) {
  mean(x == 0)
}
