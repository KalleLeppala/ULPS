#' Stable covariance
#'
#' The sample covariance matrices needed by various functions of the package are sometimes not positive definite, for actual or numeric reasons.
#' We use the functions `is.positive.definite()` and `cov.shrink()` from the package `corpcor` to identify and fix the issue.
#'
#' @param M An A x B data matrix.
#' @return A B x B positive definite matrix that estimates the covariance matrix of the distribution that produced `M`.
#'
#' @export
stable_covariance <- function(M) {
	MMT <- cov(M)
	if (is.positive.definite(MMT) == FALSE) {
		warning("An LD matrix was not positive definite. Shrinking.")
		capture.output({MMT <- cov.shrink(M)})
	}
	return(MMT)
}
