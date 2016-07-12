#' Compute covaiance matrix
#'
#' @param x_mat a matrix.
#' @param theta a numeric vector.
#' @param kernel_func a function.
#'
#' @return covariance matrix
#'
compute_covariance_matrix <- function(x_mat, theta,
                              kernel_func = kernel_gp_squared_exponential) {
  n <- nrow(x_mat)
  result_mat <- matrix(0, nrow = n, ncol = n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      value <- kernel_func(x_mat[i, ], x_mat[j, ], theta)
      result_mat[i, j] <- value
      result_mat[j, i] <- value
    }
  }
  result_mat
}
