#' Covariance matrices for residual autocorrelation between time points
#'
#' Returns the covariance matrix of the normally distributed residuals of a single sample over time points.
#'
#' @param variances A T-vector of residual variance components over time points.
#'                  This forms the diagonal of the finished covariance matrix;
#'                  first a matrix K with unit diagonal elements is created and then all rows and columns are multiplied with the corresponding standard deviation.
#' @param type Determines the kernel function of the Gaussian process. Available options are:
#' * "noise" White Gaussian noise. The matrix K is the T x T identity matrix. The default option.
#' * "squared" Squared exponential kernel, also known as radial basis function kernel or Gaussian kernel
#'             The (i, j)-element of K is defined as
#'             \ifelse{html}{\out{exp(-0.5 d(i, j)<sup>2</sup> / l<sup>2</sup> ,}}{\eqn{\exp\left(-0.5 d(i, j)^2 / l^2 \right) \,,}}
#'             where \ifelse{html}{\out{d(i, j)}}{d(i, j)} = |`points[`i`]` - `points[`j`]`| is the distance between i:th and j:th time points
#'             and \ifelse{html}{\out{l}}{l} = `params[1]` is the range parameter describing how many units away points significantly influence others.
#' * "o-u" The Ornstein–Uhlenbeck -kernel.
#' @param params The parameters
#' * parameter 1 stuff
#' * parameter 2 stuff
#'
#' @export
residual_kernel <- function(variances, type = "noise", params = NA, points = NA, sample = FALSE) {
	if (is.na(points) == TRUE) {points <- seq(1, length(variances))}
	if (type == "noise") {
	  if (is.na(params) == FALSE) {stop("Type 'noise' doesn't use parameters.")}
    covariance <- diag(length(points))
	}
	if (type == "squared") {
    if (length(params) != 1) {stop("Type 'squared' uses 1 parameter; range.")}
    l <- params[1]
    D <- abs(outer(points, points, "-")) # Distance matrix.
    covariance <- exp(-0.5*(D/l)**2)
    # OK this is totally cheating but because of numerical issues the covariance matrix
    # doesn't look positive definite and the cov.shrink() -function introduces *massive*
    # changes for some reason, so I'll just add a little something to the diagonal after all.
    covariance <- 0.999*covariance + 0.001*diag(length(points))
	}
	if (type == "o-u") {
    if (length(params) != 1) {stop("Type 'o-u' uses 1 parameter; range.")}
    l <- params[1]
    D <- abs(outer(points, points, "-")) # Distance matrix.
    covariance <- exp(-D/l)
	}
	if (type == "matern") {
	  if (length(params) != 2) {stop("Type 'matern' uses 2 parameters; smoothness and range.")}
		nu <- params[1]
		l <- params[2]
    D <- abs(outer(points, points, "-")) # Distance matrix.
    covariance <- Matern(D, range = l, smoothness = nu)
	}
  if (is.positive.definite(covariance) == FALSE) {
		warning("A residual autocorrelation matrix was not positive definite. Shrinking.")
		capture.output({covariance <- cov.shrink(covariance)})
  }
	covariance <- diag(sqrt(variances)) %*% covariance %*% diag(sqrt(variances))
	if (sample == TRUE) {
		plot(points, t(chol(covariance)) %*% rnorm(length(points)), "l")
	}
	return(covariance)
}
