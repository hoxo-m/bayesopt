#' Compute covaiance matrix
#'
#' @param x_mat a matrix.
#' @param nu a scalar.
#' @param theta a numeric vector.
#' @param kernel_func a function.
#'
#' @importFrom assertthat assert_that
#'
#' @return covariance matrix
#'
compute_covariance_matrix <- function(x_mat, nu, theta, theta0,
                              kernel_func = kernel_gp_squared_exponential) {
  n <- nrow(x_mat)
  assert_that(length(theta) == 1 || ncol(x_mat) == length(theta))
  result_mat <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in i:n) {
      value <- kernel_func(x_mat[i, ], x_mat[j, ], theta)
      if (i == j) {
        result_mat[i, i] <- value + nu
      } else {
        result_mat[i, j] <- value
        result_mat[j, i] <- value
      }
    }
  }
  theta0 * result_mat
}
