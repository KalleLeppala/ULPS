% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ULPS.R
\name{ULPS}
\alias{ULPS}
\title{Uniform longitudinal phenotype simulator}
\usage{
ULPS(X, beta, residual, C = NA)
}
\arguments{
\item{X}{A centered and scaled N x P genotype matrix.}

\item{beta}{A list of three matrices, the result created by the function \code{create_effects()}.}

\item{residual}{A T x T covariance matrix, can be created by the function \code{residual_kernel()}.}

\item{C}{An optional N x T matrix of other effects affecting the phenotype, by column t at time point t.}
}
\description{
The main function that simulates the phenotypes over time when given the effect sizes over time
and the covariance matrix describing the autocorrelation of residuals over time.
}
\seealso{
\code{\link[create_effects()]{create_effects}}, \code{\link[residual_kernel()]{residual_kernel}}
}
