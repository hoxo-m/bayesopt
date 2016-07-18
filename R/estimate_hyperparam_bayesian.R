estimate_hyperparam_bayesian <- function(x_mat, y, noise = TRUE,
                                         kernel_func = kernel_gp_squared_exponential,
                                         debug = FALSE) {
  stan_model <- readRDS("R/bayesian_inference_hyparam.rds")
  stan_data <- list(x_mat = x_mat, y = y, n_sample = nrow(x_mat), n_dim = ncol(x_mat))
  stan_fit <- rstan::sampling(stan_model, data = stan_data, chains = 1)
  rstan::extract(stan_fit)

  hyparam <- rstan::extract(stanfit, pars = c("mu", "nu", "theta0", "theta"))
  hyparam
}
