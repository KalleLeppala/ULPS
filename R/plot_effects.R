#' @export
# Oletuspaletti tunnetaan nimellä Okabe-Ito (Okabe & Ito 2008) ja sen pitäisi olla inklusiivinen.
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
  df$even <- even
  df <- df[order(df$class), ]
	figure <- ggplot(df, aes(x = time, y = effect)) +
            theme_classic() +
	          theme(legend.position = "none") +
	          geom_line(aes(group = odd, color = class)) +
            geom_line(aes(group = even, color = class)) +
            geom_point(aes(color = class))
	suppressMessages(print(figure))
}
