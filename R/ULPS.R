#' Uniform longitudinal phenotype simulator
#'
#' The main function that simulates the phenotypes over time when given the effect sizes over time
#' and the covariance matrix describing the autocorrelation of residuals over time.
#'
#' @param X A centered and scaled N x P genotype matrix.
#' @param beta A list of three matrices, the result created by the function `create_effects()`.
#' @param residual A T x T covariance matrix, can be created by the function `residual_kernel()`.
#' @param C An optional N x T matrix of other effects affecting the phenotype, by column t at time point t.
#' @result A N x T matrix of simulated phenotypes, each column corresponding to a time point.
#'
#' @seealso \code{\link[create_effects()]{create_effects}}, \code{\link[residual_kernel()]{residual_kernel}}
#'
#' @export
ULPS <- function(X, beta, residual, C = NA) {
  Y <- X %*% beta$beta
  suppressWarnings(if (is.na(C) == FALSE) {Y <- Y + C})
  epsilon <- matrix(rnorm(NROW(X)*NCOL(beta$beta)), nrow = NROW(X))
  epsilon <- epsilon %*% chol(residual)
  Y <- Y + epsilon
  return(Y)
}
