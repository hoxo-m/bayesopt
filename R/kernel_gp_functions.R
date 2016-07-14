#' Squared exponential kernel for GP
#'
#' @param x1 a numeric vector.
#' @param x2 a numeric vector.
#' @param theta a numeric vetor.
#'
#' @export
kernel_gp_squared_exponential <- function(x1, x2, theta) {
  exp(- (1/2) * sum(((x1 - x2)/theta)^2) )
}
