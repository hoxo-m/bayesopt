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
