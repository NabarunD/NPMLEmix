#' Evaluating Exponents
#'
#' Raising one number to the power of another number.
#' @param x A numeric value which is to be raised to an exponent.
#' @param y The numeric exponent to which x needs to be raised.
#' @return x^y from the input.
#' @export

powerxy <- function(x,y) {
  return(x^y)
}

#' Evaluating Logistic Function
#'
#' Evaluating the value of the logistic function at an input real number x.
#' @param x A numeric value at which the logistic function is to be evaluated.
#' @return The value of the logistic function at x.
#' @export
#'
logis <- function(x) {
  return(exp(x)/(1+exp(x)))
}
