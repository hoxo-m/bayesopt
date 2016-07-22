#' Acquisition Function GP-UCB (Upper Confidence Bound)
#'
#' @param beta a scalar.
#' @param kappa a scalar.
#'
#' @return acquisition function GP-UCB
#'
#' @export
acq_GP_UCB <- function(beta = log(2/1e-6), kappa = sqrt(beta)) {
  if (length(kappa) > 1) {
    i <- 1
    acq_func <- function(mu, var) {
      acqs <- mu + kappa[i] * sqrt(var)
      i <- i + 1
      max_acq <- max(acqs)
      max_inds <- which(acqs == max_acq)
      sample(max_inds, 1)
    }
  } else {
    acq_func <- function(mu, var) {
      acqs <- mu + kappa * sqrt(var)
      max_acq <- max(acqs)
      max_inds <- which(acqs == max_acq)
      sample(max_inds, 1)
    }
  }
  acq_func
}


#' Acquisition Function GP-MI (Mutual Information)
#'
#' @param delta a scalar.
#' @param alpha a scalar.
#' @param kappa a scalar.
#'
#' @return acquisition function GP-UCB
#'
#' @export
acq_GP_MI <- function(delta = 1e-6, alpha = log(2/delta), kappa = sqrt(alpha)) {
  sum_prev_vars <- 0
  acq_func <- function(mu, var) {
    acqs <- mu + kappa * (sqrt(var + sum_prev_vars) - sqrt(sum_prev_vars))
    max_acq <- max(acqs)
    max_inds <- which(acqs == max_acq)
    ind <- sample(max_inds, 1)
    sum_prev_vars <- sum_prev_vars + var[ind]
    ind
  }
  acq_func
}

#' Acquisition Function GP-UCB (Upper Confidence Bound)
#'
#' @param beta a scalar.
#' @param kappa a scalar.
#'
#' @return acquisition function GP-UCB
#'
#' @export
acq_GP_UCB_full_bayes <- function(beta = log(2/1e-6), kappa = sqrt(beta)) {
  acq_func <- function(x_mat, y, new_x_mat) {
    if (!require(rstan)) stop("rstan is not found.")
    standata <- list(
      n_sample = nrow(x_mat),
      n_dim = ncol(x_mat),
      n_sample_new = nrow(new_x_mat),
      x_mat = x_mat,
      y = y,
      new_x_mat = new_x_mat,
      kappa = kappa
    )
    # print(file.exists("R/acq_ucb_se.stan"))
    # stanmodel <- rstan::stan_model("R/acq_ucb_se.stan")
    stanmodel <- readRDS("R/acq_ucb_se.rds")
    stanfit <- rstan::sampling(stanmodel, data = standata, chains = 1, verbose = FALSE, warmup = 200, iter = 700)
    acqs <- colMeans(rstan::extract(stanfit, pars = "acq")$acq)
    max_acq <- max(acqs)
    max_inds <- which(acqs == max_acq)
    sample(max_inds, 1)
  }
  acq_func
}
