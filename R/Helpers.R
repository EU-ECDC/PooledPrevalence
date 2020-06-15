#' Estimate Beta distibution parameters given a mean and uncertainty boundaries
#'
#' This tool may be helpful in building priors based on the Beta distribution.
#' The user need to pass the expected mean/median of the distibution and the
#' uncertainty boundaries. Then an optimizion algorithm is used to the find the
#' alpha and beta parameters of the beta that reproduce a shape as similar as
#' possible to the expected one. The method doesn't assure perfect coverage,
#' either for optimization limits or because the expected shape is not allowed
#' by the Beta distribution, but is a starting point for further exploration of
#' useful alfa and beta parameters.
#'
#' @param min.val Expected lower uncertainty bound of the ditribution.
#' @param mid.val Expected central value. Can represent the mean, the median or
#'   the mode of the distribution according to the \code{mid.val.type}
#'   parameter.
#' @param max.val Expected upper uncertainty bound of the ditribution.
#' @param levels Uncertainty level corresponding to the \code{min.val} and
#'   \code{max.val} parameters.
#' @param mid.val.type Define which parameter of the distribution \code{mid.val}
#'   refers to, either the mean, the median or the mode.
#' @param force.proper Discourage values of alpha and beta < 1 during the
#'   optimization process. Tests should be done turning this option on and off
#'   and evaluating the results
#' @param ... Pass parameters to the \code{\link{evaluate_beta_params}}
#'   function, like \code{plot = TRUE} or \code{cut.plot.at} to modify y axis
#'   limits.
#'
#' @seealso \code{\link{evaluate_beta_params}}
#'
#' @return A list with the proposed alpha and beta parameters, a description of
#'   the resulting distribution through \code{\link{evaluate_beta_params}}, and
#'   the actual output of the optimization algorithm.
#'
#' @export
#' @importFrom stats qbeta optim
#'
#' @examples
#'
#' # Estimate the alpha and beta parameter of a Beta distribution with 95%
#' # boundaries near 0.0001 and 0.1 and mean value 0.015
#'
#' get_beta_params(.0001, .015, .1, force.proper = TRUE, plot = TRUE)
#'
#' # Same inputs but now force.proper is FALSE
#'
#' get_beta_params(.0001, .015, .1, force.proper = FALSE, plot = TRUE)
#'

get_beta_params <- function(min.val, mid.val, max.val, levels = c(.025, .975), mid.val.type = c('mean', 'median', 'mode'), force.proper = TRUE, ...) {
	obj.fun <- function(params, q = levels, intended = c(min.val, mid.val, max.val)) {
		a = exp(params[1])
		b = exp(params[2])

		med.pred <- if (mid.val.type[1] == 'mean') {
			a / (a + b)
		} else if (mid.val.type[1] == 'median') {
			qbeta(.5, a, b)
		} else if (mid.val.type[1] == 'mode') {
			(a - 1) / (a + b + 2)
		}

		predicted <- c(qbeta(q, a, b), med.pred) %>% log

		res <- sum((log(intended) - predicted)^2)

		if (is.finite(res) & !is.nan(res) & if (force.proper) (a >= 1 & b >= 1) else TRUE) res else Inf
	}

	res <- optim(par = c(0, 0), fn = obj.fun, method = 'Nelder-Mead')

	list(
		alpha = res$par[1] %>% exp,
		beta = res$par[2] %>% exp,
		distr.parameters = evaluate_beta_params(res$par[1] %>% exp, res$par[2] %>% exp, ...),
		opt.results = res
	)
}


#' Inspect the shape of a beta distribution described by the \eqn{\alpha} and
#' \eqn{\beta} arguments
#'
#' This function is useful in defining priors based on the beta distribution.
#' Given the \eqn{\alpha} and \eqn{\beta} parameters, it returns relevant
#' statistics of the distribution and boundary points to help understand its
#' characteristics. Optionally the full PDF of the distribution is plotted.
#'
#' @param alpha \eqn{\alpha} Shape parameter of the Beta distribution.
#' @param beta \eqn{\beta} Shape parameter of the Beta distribution.
#' @param lower.bound Lower bound of the distribution to report.
#' @param upper.bound Upper bound of the distribution to report.
#' @param plot Whether to plot the PDF of the distribution.
#' @param cut.plot.at Since distribution defined by \code{alpha} or \code{beta}
#'   < 1 are strongly concentrated at 0 or 1, squeezing down the rest of the PDF
#'   in the plot, it may be useful to cut the y axis in order to show the
#'   interesting parts of the distribution.
#'
#'
#' @return A vector with the lower/upper bound chosen by the user, the median,
#'   the mean, the standard deviation and the mode of the distribution. If
#'   \code{alpha} or \code{beta} are below 1, the mode cannot be calculated
#'   through \code{(alpha - 1)/(alpha + beta - 2)} so it is approximated through
#'   grid search
#'
#' @export
#'
#' @importFrom magrittr set_names
#'
#' @examples
#'
#' # Alpha and beta chosen in order to have the mean
#' # alpha/(alpha + beta) equal to .05 and precision alpha + beta
#' # equal to 30
#'
#' evaluate_beta_params(.05 * 30, (1 - .05) * 30)
#'
#' # This time the alpha and beta are chosen to have a mode
#' # (alpha - 1)/(alpha + beta - 2) equal to 0.5, with the same precision
#'
#' evaluate_beta_params(.05 * (30 - 2) + 1, (1 - .05) * (30 - 2) + 1)
#'
#' # If we again center the mean at 0.05 but this time the precision is lowered to 15,
#' # alpha goes below 1 and this produces a distribution with most of the mass
#' # at zero. The mode value calculated analytically in this case would be wrong (negative)
#' # so the empirical mode is shown as an approximation
#'
#' evaluate_beta_params(.05 * 15, (1 - .05) * 15)
#'
#'
#'

evaluate_beta_params <- function(alpha, beta, lower.bound = .025, upper.bound = .975, plot = T, cut.plot.at = NA) {

	mid <- .5

	a <- alpha
	b <- beta

	precis <- .000001

	minq <- qbeta(precis, a, b)
	maxq <- qbeta(1 - precis, a, b)

	grid <- dplyr::tibble(
		prev = seq(minq, maxq, 0.0001),
		like = dbeta(prev, a, b),
		post = like/sum(like)
	)

	params <- c(
		qbeta(c(lower.bound, mid, upper.bound), a, b) %>% set_names(c(lower.bound, mid, upper.bound) %>% percent),
		mean = a/(a+b),
		sd = a * b / ((a + b)^2 * (a + b + 1)),
		mode = if (a < 1 | b < 1) grid$prev[which.max(grid$like)] else (a - 1)/(a + b - 2)
	)

	if (plot) {
		p <- ggplot(grid, aes(prev, post)) +
			geom_line() +
			geom_vline(xintercept = c(params[c(1, 3)]), linetype = 'dashed', alpha = .8) +
			geom_point(
				data = data.frame(
					val = params[-5],
					like = dbeta(params[-5], a, b),
					col = c('Interval', 'Median', 'Interval', 'Mean', 'Mode')
				),
				aes(x = val, y = like/sum(grid$like), color = col), size = 3, alpha = .5) +
			# geom_point(aes(x = params['50%'], y = dbeta(params['50%'], a, b), color = 'Median'), size = 3, alpha = .5) +
			# geom_point(aes(x = params['mode'], y = dbeta(params['mode'], a, b), color = 'Mode'), size = 3, alpha = .5) +
			# geom_point(aes(x = params['2.5%'], y = dbeta(params['2.5%'], a, b), color = 'Interval'), size = 3, alpha = .5) +
			# geom_point(aes(x = params['97.5%'], y = dbeta(params['97.5%'], a, b), color = 'Interval'), size = 3, alpha = .5) +
			scale_x_continuous(labels = percent) +
			scale_y_continuous(labels = percent) +
			coord_cartesian(xlim = c(NA, maxq), ylim = c(NA, cut.plot.at)) +
			theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
			labs(x = 'Value', y = 'PDF', color = NULL) +
			options("PooledPrevalence.ggtheme")

		print(p)
	}

	params
}
