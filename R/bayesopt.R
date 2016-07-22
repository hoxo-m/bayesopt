#' Bayesian Optimization
#'
#' @param objective_func a function.
#' @param ... list of parameters.
#' @param iter an integer.
#' @param kernel a character string.
#' @param acq_func a function.
#' @param plot logical.
#'
#' @return a list
#'
#' @importFrom GPfit GP_fit predict.GP plot.GP
#' @importFrom graphics legend
#' @export
bayesopt <- function(objective_func, ..., iter = 10, init_size = 3,
                     kernel = c("matern5_2", "matern3_2", "square_exp"),
                     acq_func = acq_GP_UCB(), plot = FALSE) {
  # Prepare -----------------------------------------------------------------
  objective_func <- match.fun(objective_func)
  kernel <- match.arg(kernel)
  parameter_list <- list(...)
  dimension <- length(parameter_list)
  grid <- unname(expand.grid(parameter_list, KEEP.OUT.ATTRS = FALSE))
  evaluate <- function(i) do.call(objective_func, args = as.list(grid[i,]))
  if(kernel == "matern5_2") {
    gpfit <- function(inds, ys) {
      GP_fit(grid[inds, ], ys, corr = list(type = "matern", nu = 5/2))
    }
  } else if(kernel == "matern3_2") {
    gpfit <- function(inds, ys) {
      GP_fit(grid[inds, ], ys, corr = list(type = "matern", nu = 3/2))
    }
  } else {
    gpfit <- function(inds, ys) {
      GP_fit(grid[inds, ], ys, corr = list(type = "exponential", power = 2))
    }
  }
  # Initialize --------------------------------------------------------------
  inds <- sample(nrow(grid), size = init_size)
  ys <- c()
  for(i in inds) {
    y <- evaluate(i)
    ys <- c(ys, y)
    message(sprintf("input: %s, output: %f", paste(as.character(grid[i,]), collapse = ", "), y))
  }
  # Search ------------------------------------------------------------------
  gps <- vector("list", length = iter - init_size)
  for(i in seq_len(iter - init_size)) {
    gp <- gpfit(inds, ys)
    gps[[i]] <- gp
    if(dimension == 1 && plot) {
      plot.GP(gp, range=range(grid))
      plot(objective_func, add=TRUE, col=3)
      legend("topleft", legend = c("true", "pred"), col = c(3, 4), lty = 1)
    }
    # pred <- predict.GP(gp, grid[-inds, ])
    pred <- predict.GP(gp, grid)
    next_ind <- acq_func(pred$Y_hat, pred$MSE)
    # next_ind <- next_ind + sum(inds < seq_len(nrow(grid))[-inds][next_ind])
    y <- evaluate(next_ind)
    message(sprintf("input: %s, output: %f", paste(as.character(grid[next_ind,]), collapse = ", "), y))
    inds <- c(inds, next_ind)
    ys <- c(ys, y)
  }
  # Result ------------------------------------------------------------------
  max_y_ind <- which.max(ys)
  max_y <- ys[max_y_ind]
  max_x <- grid[inds[max_y_ind], ]
  list(opt_x = max_x, opt_y = max_y, gp = gps, ys = ys)
}
