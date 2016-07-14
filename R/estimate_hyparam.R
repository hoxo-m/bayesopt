#' Estimate Hyper Parameters of GP by Empirical Bayes Method
#'
#' @param x_mat a matrix.
#' @param y a numeric vector.
#' @param kernel_func a function.
#'
#' @importFrom assertthat assert_that
#' @importFrom stats optim
#'
estimate_hyperparam <- function(x_mat, y, kernel_func = kernel_gp_squared_exponential) {
  n_dim <- ncol(x_mat)

  trans_func <- function(x) x

  obj_func <- function(params) {
    nu <- trans_func(params[1])
    theta <- trans_func(params[-1])
    if (nu < 0 || any(theta < 0)) return(-Inf)
    cov_mat <- compute_covariance_matrix(x_mat, nu = nu, theta = theta,
                                         kernel_func = kernel_func)
    y_colvec <- matrix(y)
    exponent <- - (1/2) * ( t(y_colvec) %*% solve(cov_mat) %*% y_colvec )
    denominator <- (1/2) * n_dim * log(2 * pi) + (1/2) * log(det(cov_mat))
    log_likelihood <- exponent - denominator
    log_likelihood
  }

  initial_vars <- c(nu = 1, rep(1, n_dim))

  # pre <- optim(initial_vars, obj_func, control = list(trace = 1, fnscale = -1, maxit = 500), method = "BFGS")
  # assert_that(pre$convergence == 0)
  result <- optim(initial_vars, obj_func, control = list(trace = 0, fnscale = -1, maxit = 500))
  assert_that(result$convergence == 0)
  opt_nu <- unname(trans_func(result$par[1]))
  opt_theta <- unname(trans_func(result$par[-1]))
  list(nu = opt_nu, theta = opt_theta)
}

