#' @export
ULPS <- function(X, beta, residual, C = NA) {
  Y <- X %*% beta$beta
  suppressWarnings(if (is.na(C) == FALSE) {Y <- Y + C})
  epsilon <- matrix(rnorm(NROW(X)*NCOL(beta$beta)), nrow = NROW(X))
  epsilon <- epsilon %*% chol(residual)
  Y <- Y + epsilon
  return(Y)
}
