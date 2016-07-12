#' Acquisition Function GP-UCB (Upper Confidence Bound)
#'
#' @export
acq_GP_UCB <- function(beta = log(2/1e-6), kappa = sqrt(beta)) {
  if (length(kappa) > 1) {
    i <- 1
    acq_func <- function(mu, var) {
      acqs <- mu + kappa[i] * sqrt(var)
      i <- i + 1
      which.max(acqs)
    }
  } else {
    acq_func <- function(mu, var) {
      acqs <- mu + kappa * sqrt(var)
      which.max(acqs)
    }
  }
  acq_func
}


#' Acquisition Function GP-MI (Mutual Information)
#'
#' @export
acq_GP_MI <- function(delta = 1e-6, alpha = log(2/delta), kappa = sqrt(alpha)) {
  sum_prev_vars <- 0
  acq_func <- function(mu, var) {
    acqs <- mu + kappa * (sqrt(var + sum_prev_vars) - sqrt(sum_prev_vars))
    ind <- which.max(acqs)
    sum_prev_vars <- sum_prev_vars + var[ind]
    ind
  }
  acq_func
}
