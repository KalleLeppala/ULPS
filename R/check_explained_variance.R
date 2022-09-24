#' Check explained variance
#'
#' Prints the amount \out{&#946;<sup>T</sup> R &#946;} of variance a set of QTL:s explain.
#' Here \out{&#946;} is the effect size vector of those variants, and R is their linkage disequilibrium matrix.
#' Also prints the sum of squared effect sizes (the amount of variance the QTL:s would explain if we disregarded LD) for comparison.
#'
#' @param LD Either a centered and scaled N x P genotype matrix or a P x P matrix of correlation coefficients between the P genetic variants from an external source.
#' @param effects A list of three matrices, the result created by the function `create_effects()`.
#' @param j The time point.
#' @param i The effectiveness class; if not provided will consider all classes lumped together.
#'
#' @seealso \code{\link[create_effects()]{create_effects}}
#'
#' @export
check_explained_variance <- function(LD, effects, j, i) {
	loci <- effects$loci; size <- effects$size
	if (mode(LD) == "numeric") {
		if (is.matrix(LD) == TRUE) {
      p <- NCOL(LD)
      if (NROW(LD) != p || sum(LD[, 1]) == 0) {LD_style <- "genotypes"}
      else {LD_style <- "correlation"}
		} else {stop("LD not a matrix")}
	} else {stop("LD needs to contain numerical values")}
	if (missing(i)) {
		all_loci <- numeric(0); all_size <- numeric(0)
		for (i in seq(1, length(loci))) {all_loci <- c(all_loci, loci[[i]][[j]]); all_size <- c(all_size, size[[i]][[j]])}
		if (length(all_loci) > 0) {
			if (LD_style == "correlation") {temp <- LD[all_loci, all_loci]}
      else {temp <- suppressWarnings(stable_covariance(as.matrix(LD[, all_loci])))}
      variance <- t(all_size) %*% temp %*% all_size
      print(paste("Time point ", j, ", all classes together: variance = ", variance, ", sum of squared effects = ", sum(all_size**2), sep = ""))
		} else {
      print(paste("Time point ", j, ", all classes together: variance = 0, sum of squared effects = 0 (empty set)", sep = ""))
		}
	} else if (length(loci[[i]][[j]]) > 0) {
    if (LD_style == "correlation") {temp <- LD[loci[[i]][[j]], loci[[i]][[j]]]}
    else {temp <- suppressWarnings(stable_covariance(as.matrix(LD[, loci[[i]][[j]]])))}
    variance <- t(size[[i]][[j]]) %*% temp %*% size[[i]][[j]]
    print(paste("Time point ", j, ", class ", i, ": variance = ", variance, ", sum of squared effects = ", sum(size[[i]][[j]]**2), sep = ""))
	} else {
		print(paste("Time point ", j, ", class ", i, ": variance = 0, sum of squared effects = 0 (empty set)", sep = ""))
	}
}
