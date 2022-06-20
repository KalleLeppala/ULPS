# Connections between variance components and other things.
convert <- function(varY = NA, varX = NA, varG = NA, varR = NA, covXG = NA, b2 = NA) {
  if (is.vector(varY) == FALSE) {stop("Argument 'varY' not a vector.")}
  if (is.vector(varX) == FALSE) {stop("Argument 'varX' not a vector.")}
  if (is.vector(varG) == FALSE) {stop("Argument 'varG' not a vector.")}
  if (is.vector(varR) == FALSE) {stop("Argument 'varR' not a vector.")}
  if (is.vector(covXG) == FALSE) {stop("Argument 'covXG' not a vector.")}
  if (is.vector(b2) == FALSE) {stop("Argument 'b2' not a vector.")}
	result <- list(varY = varY, varX = varX, varG = varG, varR = varR, covXG = covXG, b2 = b2)
	# If varY and varG are missing but the others are not:
	if (is.na(result$varY) && is.na(result$varG) && sum(is.na(result)) == 4) {
		result$varY <- (result$varX + 2*result$covXG + result$varR)/(1 - result$b2)
		result$varG <- result$b2*result$varY
	}
	if (is.na(result$b2) && !is.na(result$varG) && !is.na(result$varY)) {result$b2 <- result$varG/result$varY}
	if (!is.na(result$b2) && is.na(result$varG) && !is.na(result$varY)) {result$varG <- result$b2*result$varY}
	if (!is.na(result$b2) && !is.na(result$varG) && is.na(result$varY)) {result$varY <- result$varG/result$b2}
	if (is.na(result$varY) && !is.na(result$varX) && !is.na(result$varG) && !is.na(result$covXG) && !is.na(result$varR)) {
		result$varY <- result$varX + result$varG + 2*result$covXG + result$varR
	}
	if (!is.na(result$varY) && is.na(result$varX) && !is.na(result$varG) && !is.na(result$covXG) && !is.na(result$varR)) {
		result$varX <- result$varY - result$varG - 2*result$covXG - result$varR
	}
  if (!is.na(result$varY) && !is.na(result$varX) && is.na(result$varG) && !is.na(result$covXG) && !is.na(result$varR)) {
	  result$varG <- result$varY - result$varX - 2*result$covXG - result$varR
	}
  if (!is.na(result$varY) && !is.na(result$varX) && !is.na(result$varG) && is.na(result$covXG) && !is.na(result$varR)) {
		result$covXG <- 0.5*(result$varY - result$varX - result$varG - result$varR)
	}
  if (!is.na(result$varY) && !is.na(result$varX) && !is.na(result$varG) && !is.na(result$covXG) && is.na(result$varR)) {
		result$varR <- result$varY - result$varX - result$varG - 2*result$covXG
	}
	if (is.na(result$b2) && !is.na(result$varG) && !is.na(result$varY)) {result$b2 <- result$varG/result$varY}
	if (!is.na(result$b2) && is.na(result$varG) && !is.na(result$varY)) {result$varG <- result$b2*result$varY}
	if (!is.na(result$b2) && !is.na(result$varG) && is.na(result$varY)) {result$varY <- result$varG/result$b2}
}
