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

  obj_func <- function(params) {
    nu <- params[1]
    theta0 <- params[2]
    theta <- params[-(1:2)]
    if (nu < 0 || theta0 < 0 || any(theta < 0)) return(-Inf)
    cov_mat <- compute_covariance_matrix(x_mat, nu = nu, theta = theta, theta0 = theta0,
                                         kernel_func = kernel_func)
    if (cov_mat[1,1] > var(y)) return(-Inf)
    y_colvec <- matrix(y)
    exponent <- - (1/2) * ( t(y_colvec) %*% solve(cov_mat) %*% y_colvec )
    denominator <- (1/2) * n_dim * log(2 * pi) + (1/2) * log(det(cov_mat))
    log_likelihood <- exponent - denominator
    log_likelihood
  }

  initial_vars <- c(nu = 1, theta0 = 1, rep(1, n_dim))
  # initial_vars <- c(theta0 = 1, rep(1, n_dim))

  # pre <- optim(initial_vars, obj_func, control = list(trace = 1, fnscale = -1, maxit = 500), method = "BFGS")
  # assert_that(pre$convergence == 0)
  result <- optim(initial_vars, obj_func, control = list(trace = 1, fnscale = -1, maxit = 5000))
  assert_that(result$convergence == 0)
  opt_nu <- unname(result$par[1])
  opt_theta0 <- unname(result$par[2])
  opt_theta <- unname(result$par[-(1:2)])
  list(nu = opt_nu, theta0 = opt_theta0, theta = opt_theta)
}


estimate_hyperparam2 <- function(x_mat, y, kernel_func = kernel_gp_squared_exponential) {
  n_dim <- ncol(x_mat)

  obj_func <- function(params) {
    theta0 <- params[1]
    theta <- params[-(1)]
    if (theta0 < 0 || any(theta < 0)) return(-Inf)
    cov_mat <- compute_covariance_matrix(x_mat, nu = 0, theta = theta, theta0 = theta0,
                                         kernel_func = kernel_func)
    y_colvec <- matrix(y)
    exponent <- - (1/2) * ( t(y_colvec) %*% solve(cov_mat) %*% y_colvec )
    denominator <- (1/2) * n_dim * log(2 * pi) + (1/2) * log(det(cov_mat))
    log_likelihood <- exponent - denominator
    log_likelihood
  }

  # initial_vars <- c(nu = 1, theta0 = 1, rep(1, n_dim))
  initial_vars <- c(theta0 = 1, rep(1, n_dim))

  # pre <- optim(initial_vars, obj_func, control = list(trace = 1, fnscale = -1, maxit = 500), method = "BFGS")
  # assert_that(pre$convergence == 0)
  result <- optim(initial_vars, obj_func, control = list(trace = 1, fnscale = -1, maxit = 500))
  assert_that(result$convergence == 0)
  opt_theta0 <- unname(result$par[1])
  opt_theta <- unname(result$par[-(1)])
  list(theta0 = opt_theta0, theta = opt_theta)
}

