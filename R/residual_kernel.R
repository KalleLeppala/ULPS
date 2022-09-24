#' Residual kernel
#'
#' Creates a covariance matrix of the independent normally distributed residuals autocorrelated over time points.
#' Four kernel functions available for the generation of this matrix.
#' Also able to visualize independent samples using the matrix created.
#'
#' @param variances A T-vector of residual variance components over time points.
#'                  This forms the diagonal of the finished covariance matrix;
#'                  first a matrix K with unit diagonal elements is created and then all rows and columns are multiplied with the corresponding standard deviation.
#' @param type Determines the kernel function of the Gaussian process. Available options are:
#' * "noise" White Gaussian noise. The matrix K is the T x T identity matrix. No `params` needed. The default option.
#' * "squared" Squared exponential kernel, also known as radial basis function kernel or Gaussian kernel.
#'             The (i, j)-element of K is defined as \out{exp(-0.5 d(i, j)<sup>2</sup> / l<sup>2</sup>)},
#'             where d(i, j) is `points[j] - points[i]` and l is the range parameter `params[1]`, the only parameter needed for `type` "squared".
#' * "o-u" The Ornstein-Uhlenbeck -kernel, also known as the AR(1)-kernel.
#'         The (i, j)-element of K is defined as \out{exp(-0.5 |d(i, j)| / l)},
#'         where d(i, j) is `points[j] - points[i]` and l is the range parameter `params[1]`, the only parameter needed for `type` "o-u".
#' * "matern" The Mat\out{&#233;}rn-kernel.
#'            The (i, j)-element of K is defines as \out{	2<sup>1- &#957;</sup> (sqrt(2&#957;)|d(i, j)|/l) K<sub>&#957;</sub>((sqrt(2&#957;)|d(i, j)|/l))<sup>&#957;</sup>  / &#915;(&#957;)},
#'            where \out{K<sub>&#957;</sub>()} is the modified Bessel function, \out{&#915;()} is the Gamma-function, d(i, j) is `points[j] - points[i]`,
#'            \out{&#957;} is the smoothness parameter `params[1]` and the l is the range parameter `params[2]`.
#' @param params The parameters required, the amount and meaning of which depends on `type`.
#' @param points An optional T-vector if the time points are not evenly spaced.
#'               By default, distance between consecutive time points is one unit from the range parameters point of view.
#' @param sample The function can visualize samples of the Gaussian process.
#'               The positive integer given is the number of independent samples, the default value FALSE produces no plot.
#' @param palette The palette used for plotting. Default is the inclusive Okabe-Ito palette without black.
#' @return A T x T covariance matrix if `sample` is FALSE, a ggplot if it's a positive integer.
#'
#' @seealso \code{\link[ULPS()]{ULPS}}, \code{\link[create_effects()]{create_effects}}
#'
#' @export
residual_kernel <- function(variances, type = "noise", params = NA, points = NA, sample = FALSE, palette = palette.colors()[2:9]) {
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
	if (sample != FALSE) {
		okornot <- FALSE
		if (is.numeric(sample) == TRUE) {
			if (sample > 0 && sample%%1==0) {okornot <- TRUE}
		}
		if (okornot == FALSE) {stop("Parameter 'sample' is the number of residuals visualized.")}
		time = points
		figure <- ggplot() +
              theme_classic() +
			        ylab("residuals") +
	            theme(legend.position = "none")
    for (j in seq(0, sample - 1)) {
	    figure <- figure + geom_line(aes(x = time, y = t(chol(covariance)) %*% rnorm(length(time))), colour = palette[j%%length(palette) + 1])
    }
    suppressMessages(print(figure))
	}
	if (sample == FALSE) {return(covariance)}
	else {return(figure)}
}
