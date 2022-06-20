#' Auxiliary function for simulating effect sizes using c distinct ''effectiveness classes''
#'
#' @param LD Either a centered and scaled n x p genotype matrix or a p x p matrix of correlation coefficients between the p genetic variants from an external source.
#' @param amounts A c x t matrix of the number of QTL:s in each effectiveness class at each time point.
#' @param variances A c x t matrix of variances each effectiveness class explains ignoring the LD with other classes.
#' @param changes A c x (t - 1) matrix containing the number of QTL:s that change between time points in each effectiveness class, in addition to changes dictated by the matrix amounts.
#'                Must have changes[i, j] at most min(amounts[i, j], amounts[i, j + 1]). Not required when t = 1.
#' @param shuffles A c x (t - 1) logical matrix coding whether the QTL:s that are kept between time points have their effects re-randomized or not.
#'                 Not required when t = 1.
#' @param subsets A list of c subsets of the set {1, ..., p} to which indices of the QTL:s of each effectiveness class are restricted. By default there are no restrictions.
#' @param max_reject Maximum number of times the rejection sampling loop is executed. If reached, accepts the sample but gives a warning. Default value 100.
#' @param prioritize_total_variance If TRUE, will scale the effect sizes so that the variance explained by all the effectiveness classes jointly is the column sum of variances.
#'                                  At the cost of accuracy on the variance explained by classes individually. Default value FALSE.
#' @return A p x t matrix of effect sizes.
#'
#' @export
create_effects <- function(LD, amounts, variances, changes = NA, shuffles = NA, subsets = NA, max_reject = 100, prioritize_total_variance = FALSE) {
	# Check if everything is OK with the arguments.
	if (mode(LD) == "numeric") {
		if (is.matrix(LD) == TRUE) {
      p <- NCOL(LD)
      if (NROW(LD) != p || sum(LD[, 1]) == 0) {LD_style <- "genotypes"}
      else {LD_style <- "correlation"}
		} else {stop("Argument 'LD' not a matrix.")}
	} else {stop("Argument 'LD' needs to contain numerical values.")}
	if (mode(amounts) == "numeric") {
		if (is.matrix(amounts) == TRUE) {
		  c <- NROW(amounts)
		  t <- NCOL(amounts)
	  } else {stop("Argument 'amounts' not a matrix.")}
	} else {stop("Argument 'amounts' needs to contain numerical values.")}
	if (mode(variances) == "numeric") {
		if (is.matrix(variances) == TRUE) {
		  if (NROW(variances) != c) {stop("Mismatch in numbers of effectiveness classes.")}
			if (NCOL(variances) != t) {stop("Mismatch in numbers of time points.")}
		} else {stop("Argument 'variances' not a matrix.")}
	} else {stop("Argument 'variances' needs to contain numerical values.")}
	if (t > 1) {
		if (mode(changes) == "numeric") {
		  if (is.matrix(changes) == TRUE) {
			  if (NROW(changes) != c) {stop("Mismatch in numbers of effectiveness classes.")}
			  if (NCOL(changes) != t - 1) {stop("Mismatch in numbers of time points.")}
		  } else {stop("Argument 'changes' not a matrix.")}
	  } else {stop("Argument 'changes' needs to contain numerical values.")}
	  if (mode(shuffles) == "logical") {
		  if (is.matrix(shuffles) == TRUE) {
			  if (NROW(shuffles) != c) {stop("Mismatch in numbers of effectiveness classes.")}
			  if (NCOL(shuffles) != t - 1) {stop("Mismatch in numbers of time points.")}
		  } else {stop("Argument 'shuffles' not a matrix.")}
	  } else {stop("Argument 'shuffles' needs to contain logical values.")}
	}
  # Check if the special requirement of the argument changes is met.
  if (t > 1) { for (j in seq(1, t - 1)) {
  	for (i in seq(1, c)) {
  	  if (changes[i, j] > amounts[i, j] || changes[i, j] > amounts[i, j + 1]) {stop("changes contains too big values compared to amounts.")}
  	}
  } }
	# Creating the argument subset if it's not given.
	if (is.na(subsets) == TRUE) {
		subsets <- list()
	  for (i in seq(1, c)) {subsets[[length(subsets) + 1]] <- seq(1, p)}
	}
	# OK it seems all is well. Starting the program.
	free <- seq(1, p) # The variants currently free, that can be picked as a QTL at some effectiveness class.
	loci <- list(); size <- list() # Keeping track of the QTL locations and sizes, initializing as lists of lists containing empty vectors of correct length.
	for (i in seq(1, c)) {
		loci[[i]] <- list(); size[[i]] <- list()
		for (j in seq(1, t)) {loci[[i]][[j]] <- numeric(amounts[i, j]); size[[i]][[j]] <- numeric(amounts[i, j])}
	}
	# The first time point:
	for (i in seq(1, c)) { if (amounts[i, 1] > 0) {
 	  loci[[i]][[1]] <- sample(intersect(free, subsets[[i]]), length(loci[[i]][[1]])) # Randomizing QTL loci.
		free <- setdiff(free, loci[[i]][[1]]) # Removing the chosen indices from the pool of currently free indices.
    # temp is the positive definite LD matrix of the relevant variants.
		if (LD_style == "correlation") {temp <- LD[loci[[i]][[1]], loci[[i]][[1]]]}
		else {temp <- stable_covariance(as.matrix(LD[, loci[[i]][[1]]]))}
		eigen <- eigen(temp)
		U <- eigen$vectors
	  Sinv <- diag(1/sqrt(eigen$values))
		counter <- 1
		while (TRUE) { # Rejection sampling because ellipses are not spheres.
      sizes <- rnorm(amounts[i, 1]) # Randomizing QTL effect sizes, first they're just some multinormal.
      sizes <- sqrt(variances[i, 1])*sizes/sqrt(sum(sizes**2)) # Scaling.
      if (amounts[i, 1] == 1) {factor <- 1; maxfactor <- 1}
      else {
      	V <- cbind(qr.Q(qr(cbind(sizes, matrix(rnorm(amounts[i, 1]*(amounts[i, 1] - 1)), nrow = amounts[i, 1]))))[, -1]) # Some basis vectors orthogonal to alpha.
        factor <- sqrt(det(t(V) %*% Sinv %*% Sinv %*% V)) # The amount of distortion.
        maxfactor <- prod(diag(Sinv)[2:amounts[i, 1]]) # Maximal amount of distortion.
      }
      if (runif(1) < factor/maxfactor) {break()}
      counter <- counter + 1
      if (counter > max_reject) {warning(paste("Effects at class ", i, ", time 1 kept due to reaching the maximum number of rejections allowed.\n Final acceptance rate was ", factor/maxfactor, ".", sep = "")); break()}
		}
    sizes <- U %*% Sinv %*% sizes # Transformed so that the LD within the effectiveness class is accounted for.
    size[[i]][[1]] <- c(sizes) # Saved.
	}
	print(paste("Class ", i, ", time point 1 done.", sep = ""))
	}
	# The rest of the time points:
	if (t > 1) {
  	for (j in seq(2, t)) {
	  	for (i in seq(1, c)) { if (amounts[i, j] > 0) {
		  	# We keep min(amounts[i, j - 1], amounts[i, j]) - changes[i, j - 1] old loci.
  			keep <- min(amounts[i, j - 1], amounts[i, j]) - changes[i, j - 1]
	  		# And introduce amounts[i, j] - min(amounts[i, j - 1], amounts[i, j]) + changes[i, j - 1] new loci.
        introduce <- amounts[i, j] - keep
	   	  if (keep == 0) {
          loci[[i]][[j]] <- c(sample(intersect(free, subsets[[i]]), introduce)) # Randomizing QTL loci.
          free <- setdiff(c(free, na.omit(loci[[i]][[j - 1]][keep + 1:amounts[i, j - 1]])), loci[[i]][[j]][keep + 1:amounts[i, j]]) # Updating the pool of currently free indices.
        } else {
          loci[[i]][[j]] <- c(loci[[i]][[j - 1]][1:keep], sample(intersect(free, subsets[[i]]), introduce)) # Randomizing QTL loci.
          free <- setdiff(c(free, na.omit(loci[[i]][[j - 1]][keep + 1:amounts[i, j - 1]])), loci[[i]][[j]][keep + 1:amounts[i, j]]) # Updating the pool of currently free indices.
        }
        # temp is the positive definite LD matrix of all relevant variants, kept or introduced.
        if (LD_style == "correlation") {temp <- LD[loci[[i]][[j]], loci[[i]][[j]]]}
        else {temp <- stable_covariance(as.matrix(LD[, loci[[i]][[j]]]))}
        eigen <- eigen(temp)
		    U <- eigen$vectors
		    Sinv <- diag(1/sqrt(eigen$values))
		    counter <- 1
        while (TRUE) { # Rejection sampling because ellipses are not spheres.
  		    # Are we randomizing the effects of the kept variants again or keeping them?
	  	    if (shuffles[i, j - 1] == TRUE || keep == 0) { # Randomizing, easy.
	    	    sizes <- rnorm(amounts[i, j]) # Randomizing QTL effect sizes, first they're just some multinormal.
            sizes <- sqrt(variances[i, j])*sizes/sqrt(sum(sizes**2)) # Scaling.
     	    } else { # Keeping, a bit harder.
    	  	  if (introduce == 0) {
    		    	size[[i]][[j]] <- size[[i]][[j - 1]][1:keep]
  	  	    } else {
  		        choltemp <- diag(sqrt(eigen$values)) %*% t(U)
  		        Omega <- as.matrix(choltemp[, 1:keep]) %*% size[[i]][[j - 1]][1:keep] # A fixed point (or vector).
 	  		      QQ <- qr.Q(qr(as.matrix(choltemp[, (keep + 1):NCOL(choltemp)]))) # An orthogonal basis of the subspace V.
    		      center <- (diag(length(loci[[i]][[j]])) - QQ %*% t(QQ)) %*% Omega # Center of the little sphere.
    		      sizes <- rnorm(introduce) # Randomizing QTL effect sizes, first they're just some multinormal.
  	  	   	  sizes <- as.matrix(choltemp[, (keep + 1):NCOL(choltemp)]) %*% sizes # Rotating.
  	  	   	  if (variances[i, j] < sum(center**2)) {warning(paste("Because of the effects that are kept, class ", i, " at time ", j, " explains more variance than is allocated for it.\n  Consider re-running the simulation.", sep = ""))}
    		      sizes <- sqrt(max(variances[i, j] - sum(center**2), 0))*sizes/sqrt(sum(sizes**2)) # Scaling.
              sizes <- center + sizes # Translating.
  	  	    }
     	    }
        	if (amounts[i, j] == 1) {factor <- 1; maxfactor <- 1}
          else {
       	    V <- cbind(qr.Q(qr(cbind(sizes, matrix(rnorm(amounts[i, j]*(amounts[i, j] - 1)), nrow = amounts[i, j]))))[, -1]) # Some basis vectors orthogonal to alpha.
            factor <- sqrt(det(t(V) %*% Sinv %*% Sinv %*% V)) # The amount of distortion.
            maxfactor <- prod(diag(Sinv)[2:amounts[i, j]]) # Maximal amount of distortion.
          }
          if (runif(1) < factor/maxfactor) {break()}
        	counter <- counter + 1
        	if (counter > max_reject) {warning(paste("Effects at class ", i, ", time ", j, " kept due to reaching the maximum number of rejections allowed.\n Final acceptance rate was ", factor/maxfactor, ".", sep = "")); break()}
        	}
        sizes <- U %*% Sinv %*% sizes # Transformed so that the LD within the effectiveness class is accounted for.
        size[[i]][[j]] <- c(sizes) # Saved.
	  	  # The order of indices and effect sizes are randomized so that it's easy to sample a random subset if needed.
	  	  random_order <- sample(length(loci[[i]][[j]]), length(loci[[i]][[j]]))
        loci[[i]][[j]] <- loci[[i]][[j]][random_order]
        size[[i]][[j]] <- size[[i]][[j]][random_order]
	  	}
	  	print(paste("Class ", i, ", time point ", j, " done.", sep = ""))
	  	}
		}
	}
	if (prioritize_total_variance == TRUE) {
    for (j in seq(1, t)) {
    	all_loci <- numeric(0); all_size <- numeric(0)
		  for (i in seq(1, c)) {all_loci <- c(all_loci, loci[[i]][[j]]); all_size <- c(all_size, size[[i]][[j]])}
  		if (length(all_loci) > 0) {
	  		if (LD_style == "correlation") {temp <- LD[all_loci, all_loci]}
        else {temp <- suppressWarnings(stable_covariance(as.matrix(LD[, all_loci])))}
      	target <- sum(variances[, j])
        variance <- t(all_size) %*% temp %*% all_size
        for (i in seq(1, c)) {size[[i]][[j]] <- sqrt(target/variance[1, 1])*size[[i]][[j]]}
  		}
    }
	}
 	beta <- matrix(0, nrow = t, ncol = p) # Creating an effect size vector.
 	for (i in seq(1, c)) {
 		for (j in seq(1, t)) {
 			for (k in seq(1, length(loci[[i]][[j]]))) {
   			beta[j, loci[[i]][[j]][k]] <- size[[i]][[j]][k]
 			}
 		}
 	}
 	beta <- t(beta) # Whoops.
 	return(list(loci = loci, size = size, beta = beta))
}
