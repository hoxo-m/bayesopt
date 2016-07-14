gp_fit <- function(x_mat, y, kernel_func = kernel_gp_squared_exponential) {
  y <- scale(y, center = TRUE, scale = FALSE)
  hyparam <- estimate_hyperparam2(x_mat, y, kernel_func = kernel_func)
  n_sample <- nrow(x_mat)
  n_dim <- ncol(x_mat)
  cov_mat <- compute_covariance_matrix(x_mat, nu = 0, theta0 = hyparam$theta0,
                                       theta = hyparam$theta,
                                       kernel_func = kernel_func)
  inv_cov_mat <- solve(cov_mat)

  gp_pred <- function(new_x_mat) {
    k_t <- apply(new_x_mat, 1, function(new_x) {
      k_vec <- apply(x_mat, 1, function(x) {
        hyparam$theta0 * kernel_func(x, new_x, hyparam$theta)
      })
      matrix(k_vec)
    })
    k <- t(k_t)
    mu <- k %*% inv_cov_mat %*% matrix(y)
    sigma2 <- apply(k, 1, function(k_vec) {
      k_colvec <- matrix(k_vec)
      hyparam$theta0 - t(k_colvec) %*% inv_cov_mat %*% k_colvec
    })
    res <- data.frame(mu = mu, sigma2 = sigma2)
    res
  }
  gp_pred
}
