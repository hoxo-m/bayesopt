#' Squared exponential kernel for GP
#'
#' @param x1 a numeric vector.
#' @param x2 a numeric vector.
#' @param theta a numeric vetor.
#'
#' @importFrom assertthat assert_that
#'
#' @export
kernel_gp_squared_exponential <- function(x1, x2, theta) {
  assert_that(length(x1) == length(x2))
  assert_that(length(theta) == 1 || length(x1) == length(theta))
  exp(- (1/2) * sum(((x1 - x2)/theta)^2) )
}
