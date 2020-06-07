.onLoad <- function (libname, pkgname)
{
	# Helpers.R
	utils::globalVariables(c('prev', 'like', 'post', 'val'))

	# Estimate.R
	utils::globalVariables(c('term'))

	plots.theme = ggplot2::theme_minimal() +
		ggplot2::theme(
			axis.text = ggplot2::element_text(size = 12, colour = 'dimgray'),
			#axis.text.x = element_text(face = "bold", size = 15, colour = 'dimgray'),
			axis.title = ggplot2::element_text(face = "bold", size = 15, colour = 'dimgray'),
			axis.title.x = ggplot2::element_text(margin = ggplot2::margin(20, 0, 0, 0)),
			axis.title.y = ggplot2::element_text(margin = ggplot2::margin(0, 20, 0, 0)),
			#axis.title.y = element_blank(),
			axis.ticks.y = ggplot2::element_blank(),
			panel.grid.major.y = ggplot2::element_line(colour = "lightgray", size = .5),
			panel.grid.minor.y = ggplot2::element_line(colour = "lightgray", size = .25),
			panel.grid.major.x = ggplot2::element_blank(),
			panel.grid.minor.x = ggplot2::element_blank(),
			panel.border = ggplot2::element_blank(),
			plot.margin = ggplot2::unit(c(1,1.5,1,1), 'cm'),
			plot.title = ggplot2::element_text(face = "bold"),
			#legend.position = 'none',
			legend.key.width = ggplot2::unit(3,"line"),
			legend.title = ggplot2::element_text(size = 14, face = "bold")
		)

	options("PooledPrevalence.ggtheme" = plots.theme)

	#rstan::rstan_options(auto_write = TRUE)

	invisible()
}
