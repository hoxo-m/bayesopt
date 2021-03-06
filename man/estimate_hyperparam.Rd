% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_hyparam.R
\name{estimate_hyperparam}
\alias{estimate_hyperparam}
\title{Estimate Hyper Parameters of GP by Empirical Bayes Method}
\usage{
estimate_hyperparam(x_mat, y, noise = TRUE,
  kernel_func = kernel_gp_squared_exponential, debug = FALSE)
}
\arguments{
\item{x_mat}{a matrix.}

\item{y}{a numeric vector.}

\item{noise}{logical.}

\item{kernel_func}{a function.}

\item{debug}{logical.}
}
\description{
Estimate Hyper Parameters of GP by Empirical Bayes Method
}

