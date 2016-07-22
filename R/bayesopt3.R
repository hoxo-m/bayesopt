#' Bayesian Optimization
#'
#' @param objective_func a function.
#' @param ... list of parameters.
#' @param iter an integer.
#' @param noise logical.
#' @param kernel a character string.
#' @param acq_func a function.
#'
#' @return a list
#'
#' @importFrom methods formalArgs
#'
#' @export
bayesopt3 <- function(objective_func, ..., iter = 10) {
  # Prepare -----------------------------------------------------------------
  objective_func <- match.fun(objective_func)
  parameter_list <- list(...)
  dimension <- length(parameter_list)
  grid <- as.matrix(expand.grid(parameter_list, KEEP.OUT.ATTRS = FALSE))
  colnames(grid) <- formalArgs(objective_func)
  evaluate <- function(i) do.call(objective_func, args = as.list(grid[i,]))
  acq_func <- acq_GP_UCB_full_bayes()
  # Initialize --------------------------------------------------------------
  initial_size <- dimension + 1
  inds <- sample(nrow(grid), size = initial_size)
  ys <- c()
  for(i in inds) {
    y <- evaluate(i)
    ys <- c(ys, y)
    message(sprintf("input: %s, output: %f", paste(as.character(grid[i,]), collapse = ", "), y))
  }
  # Search ------------------------------------------------------------------
  for(i in seq_len(iter - initial_size)) {
    next_ind <- acq_func(grid[inds, , drop=FALSE], ys, grid)
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
  list(opt_x = max_x, opt_y = max_y, ys = ys)
}
