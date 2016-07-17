#' Estimate Hyper Parameters of GP by Empirical Bayes Method
#'
#' @param x_mat a matrix.
#' @param y a numeric vector.
#' @param noise logical.
#' @param kernel_func a function.
#' @param debug logical.
#'
#' @importFrom assertthat assert_that
#' @importFrom Matrix rcond
#' @importFrom stats optim
#'
estimate_hyperparam <- function(x_mat, y, noise = TRUE,
                                kernel_func = kernel_gp_squared_exponential,
                                debug = FALSE) {
  if(!is.matrix(x_mat)) x_mat <- matrix(x_mat)
  n_dim <- ncol(x_mat)
  n_sample <- nrow(x_mat)
  opt_m <- mean(y)
  y <- y - opt_m

  trans_func <- function(x) exp(x) + 1
  # trans_func <- function(x) x

  obj_func <- function(params) {
    if (noise) {
      nu <- trans_func(params[1])
      theta <- params[-1]
      theta0 <- 1
      if (nu < 1 || any(theta <= 0)) return(min_log_likelihood)
    } else {
      nu <- 1
      theta <- params
      theta0 <- 1
      if (any(theta <= 0)) return(min_log_likelihood)
    }
    cov_mat <- compute_covariance_matrix(x_mat, nu = nu, theta = theta, theta0 = theta0,
                                         kernel_func = kernel_func)

    if (rcond(cov_mat) < 1e-10) return(min_log_likelihood)

    inv_cov_mat <- solve(cov_mat)
    y_colvec <- matrix(y)
    exponent <- - ( t(y_colvec) %*% inv_cov_mat %*% y_colvec )
    denominator <- log(det(cov_mat))
    log_likelihood <- exponent - denominator
    log_likelihood
  }

  theta_uppers <- apply(x_mat, 2, function(x) {
    (diff(range(x)) / 2) ^ 2
  })

  if (noise) {
    initial_vars <- c(nu = 0, theta = theta_uppers / 2)
  } else {
    initial_vars <- c(theta = theta_uppers / 2)
  }

  initial_value <- obj_func(initial_vars)
  min_log_likelihood <- initial_value

  if (length(initial_vars) == 1) {
    result <- optim(initial_vars, obj_func, control = list(trace = debug, fnscale = -1, maxit = 5000), method = "Brent", lower = 0, upper = theta_uppers)
  } else {
    result <- optim(initial_vars, obj_func, control = list(trace = debug, fnscale = -1, maxit = 5000))
    # result <- optim(result$par, obj_func, control = list(trace = debug, fnscale = -1, maxit = 1000), method="BFGS")
  }
  assert_that(result$convergence == 0)

  if (noise) {
    opt_nu <- unname(trans_func(result$par[1]))
    opt_theta <- unname(result$par[-1])
    opt_theta0 <- 1
  } else {
    opt_nu <- 1
    opt_theta <- unname(result$par)
    opt_theta0 <- 1
  }
  list(m = opt_m, nu = opt_nu, theta0 = opt_theta0, theta = opt_theta)
}
