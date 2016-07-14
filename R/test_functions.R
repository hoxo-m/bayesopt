#' Himmelblau's function
#'
#' @param x a numeric vector.
#' @param y a numeric vector.
#'
#' @return value
#'
#' @export
Himmelblau <- function(x, y) {
  value <- (x^2 + y - 11)^2 + (x + y^2 - 7)^2
  -value
}
