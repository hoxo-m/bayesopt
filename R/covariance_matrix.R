#' Compute covaiance matrix
#'
#' @param x_mat a matrix.
#' @param nu a scalar.
#' @param theta a numeric vector.
#' @param theta0 a numeric vector.
#' @param kernel_func a function.
#'
#' @importFrom assertthat assert_that
#'
#' @return covariance matrix
#'
compute_covariance_matrix <- function(x_mat, nu, theta, theta0,
                              kernel_func = kernel_gp_squared_exponential) {
  if(!is.matrix(x_mat)) x_mat <- matrix(x_mat)
  n <- nrow(x_mat)
  assert_that(length(theta) == 1 || ncol(x_mat) == length(theta))
  result_mat <- matrix(0, nrow = n, ncol = n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      value <- kernel_func(x_mat[i, ], x_mat[j, ], theta = theta)
      result_mat[i, j] <- value
      result_mat[j, i] <- value
    }
  }
  result_mat <- theta0 * result_mat
  diag(result_mat) <- nu
  result_mat
}
