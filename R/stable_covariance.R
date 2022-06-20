#' Description
#' The functions cov.shrink, cor.shrink, and pcor.shrink implement a shrinkage approach to estimate covariance and (partial) correlation matrices, as suggested in Schaefer and Strimmer (2005). The advantages of using this approach in comparison with the standard empirical estimates (cov and cor) are that the shrinkage estimates
#' are always positive definite,
#' well conditioned (so that the inverse always exists), and
#' exhibit (sometimes dramatically) better mean squared error.
#' Furthermore, they are inexpensive to compute and do not require any tuning parameters (the shrinkage intensity is analytically estimated from the data).
#' export
stable_covariance <- function(M) {
	# Is it better to only use shrinkage if the covariance matrix is not positive definite, or to always use it?
	MMT <- cov(M)
	if (is.positive.definite(MMT) == FALSE) {
		warning("An LD matrix was not positive definite. Shrinking.")
		capture.output({MMT <- cov.shrink(M)})
	}
	return(MMT)
}
