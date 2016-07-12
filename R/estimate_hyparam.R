#' Estimate Hyper Parameters of GP by Empirical Bayes Method
#'
#' @importFrom assertthat assert_that
#' @importFrom stats optim
#'
estimate_hyperparam <- function(y, x_mat, kernel_func = kernel_gp_squared_exponential) {
  n_dim <- ncol(x_mat)
  vars <- paste0("theta", seq_len(n_dim), collapse = ",")

  obj_func <- function(params) {
    m <- params[1]
    nu <- exp(params[2])
    theta0 <- exp(params[3])
    for (i in 1:n_dim) {
      var_name <- paste0("theta", i)
      assign(var_name, exp(params[3+i]))
    }
    theta <- eval(parse(text=sprintf("c(%s)", vars)))
    cov_mat <- theta0 * compute_covariance_matrix(x_mat, theta, kernel_func)
    diag(cov_mat) <- nu
    y_mat <- matrix(y - m)
    exponent <- - (1/2) * ( t(y_mat) %*% solve(cov_mat) %*% y_mat )
    denominator <- - (1/2) * log(det(cov_mat))
    log_likelihood <- exponent + denominator
    log_likelihood
  }

  vars <- paste0("theta", seq_len(n_dim), "=0", collapse = ",")
  initial_theta <- eval(parse(text = sprintf("c(%s)", vars)))
  initial_vars <- c(m = 1, nu = 0, theta0 = 0, initial_theta)

  result <- optim(initial_vars, obj_func, control = list(trace = 0, fnscale = -1, maxit = 500))
  assert_that(result$convergence == 0)
  opt_m <- unname(result$par["m"])
  opt_nu <- unname(exp(result$par["nu"]))
  opt_theta0 <- unname(exp(result$par["theta0"]))
  opt_theta <- unname(exp(result$par[4:length(result$par)]))
  list(m = opt_m, nu = opt_nu, theta0 = opt_theta0, theta = opt_theta)
}

