#' Squared exponential kernel for GP
#'
#' @param x1 a numeric vector.
#' @param x2 a numeric vector.
#' @param theata a numeric vetor.
#'
#' @export
kernel_gp_squared_exponential <- function(x1, x2, theta) {
  kernel_gp_exponential(x1, x2, theta, power = 2)
}

#' Exponential kernel for GP
#'
#' @param x1 a numeric vector.
#' @param x2 a numeric vector.
#' @param theata a numeric vetor.
#' @param power a scalar.
#'
#' @export
kernel_gp_exponential <- function(x1, x2, theta, power = 2) {
  exp(- sum(((x1 - x2)/theta) ^ power) ^ (1 / power))
}
