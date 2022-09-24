#' Plot effects
#'
#' Visualizes the effects of QTL:s from various effectiveness classes over time.
#'
#' @param effects A list of three matrices, the result created by the function `create_effects()`.
#' @param palette The palette used for plotting. Default is the inclusive Okabe-Ito palette without black.
#' @return Returns the ggplot.
#'
#' @seealso \code{\link[create_effects()]{create_effects}}
#'
#' @export
plot_effects <- function(effects, palette = palette.colors()[2:9]) {
	loci <- effects$loci
	size <- effects$size
	time <- numeric(0)
	class <- character(0)
	effect <- numeric(0)
  locus <- numeric(0)
	for (i in seq(1, length(loci))) {
		if (length(loci[[i]]) > 0) { for (j in seq(1, length(loci[[i]]))) {
			if (length(loci[[i]][[j]]) > 0) { for (k in seq(1, length(loci[[i]][[j]]))) {
				time[length(time) + 1] <- j
			  class[length(class) + 1] <- palette[i]
			  effect[length(effect) + 1] <- size[[i]][[j]][k]
			  locus[length(locus) + 1] <- loci[[i]][[j]][k]
			} }
		} }
	}
	df <- data.frame(time, class, effect, locus)
  df <- df[order(df$locus, df$time), ]
  odd <- c(1)
  even <- c(1)
  for (k in seq(2, NROW(df))) {
  	if (k%%2 == 0) {
  		odd[k] <- odd[k - 1] + 1
  		if (df$class[k] == df$class[k - 1] && df$locus[k] == df$locus[k - 1] && df$time[k] == df$time[k - 1] + 1) {
  			even[k] <- even[k - 1]
  		} else {
  			even[k] <- even[k - 1] + 1
  		}
  	} else {
  		even[k] <- even[k - 1] + 1
  		if (df$class[k] == df$class[k - 1] && df$locus[k] == df$locus[k - 1] && df$time[k] == df$time[k - 1] + 1) {
  			odd[k] <- odd[k - 1]
  		} else {
  			odd[k] <- odd[k - 1] + 1
  		}
  	}
  } # ...I'm not proud of this!
  df$odd <- odd
  df$even <- even + NROW(df)
  df <- df[order(df$class), ]
  come_on_now <- c(df$class, df$class)
  come_on_now <- setNames(come_on_now, c(as.character(df$odd), as.character(df$even)))
  df$odd <- factor(df$odd)
  df$even <- factor(df$even)
  df$effects <- df$effect
	figure <- ggplot(df, aes(x = time, y = effects)) +
            theme_classic() +
	          theme(legend.position = "none") +
	          geom_line(aes(group = odd, colour = odd)) +
            geom_line(aes(group = even, colour = even)) +
        	  scale_color_manual(values = come_on_now) +
            geom_point(colour = df$class)
	suppressMessages(print(figure))
	return(figure)
}
