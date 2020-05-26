#' Estimate Beta distibution parameters given a mean and uncertainty boundaries
#'
#' This tool may be helpful in building priors based on the Beta distribution.
#' The user need to pass the expected mean/median of the distibution and the
#' uncertainty boundaries. Then an optimizion algorithm is used to the find the
#' alpha and beta parameters of the beta that reproduce a shape as similar as
#' possible to the expected one. The method don't assure perfect coverage,
#' either for optimization limits or because the expected shape is not allowed
#' by the Beta distribution, but is a starting point for further esploration of
#' useful alfa and beta parameters.
#'
#' @param min.val Expected lower uncertainty bound of the ditribution.
#' @param mid.val Expected central value. Can represent the median or the mean
#'   of the distribution according to the \code{mid.val.type} parameter.
#' @param max.val Expected upper uncertainty bound of the ditribution.
#' @param levels Uncertainty level corresponding to the \code{min.val} and
#'   \code{max.val} parameters.
#' @param mid.val.type Define which parameter of the distribution \code{mid.val}
#'   refers to, either the mean or the median.
#' @param force.proper Discourage values of alpha and beta < 1 during the
#'   optimization process. Tests should be done turning this option on and off
#'   and evaluating the results
#' @param plot Whether to plot the resulting distribution using
#'   \code{\link{evaluate_beta_params}} fuction.
#'
#' @seealso \code{\link{evaluate_beta_params}}
#'
#' @return A list with the proposed alpha and beta parameters, a description of
#'   the resulting distribution through \code{\link{evaluate_beta_params}}, and
#'   the actual output of the optimization algorithm.
#' @export
#'
#' @examples
#'
#' # Estimate the alpha and beta parameter of a Beta distribution with 95% boundaries near 0.0001 and 0.1 and mean value 0.015
#'
#' get.beta.params(.0001, .015, .1, force.proper = T, plot = T)
#'
#' # Same inputs but now force.proper is FALSE
#'
#' get.beta.params(.0001, .015, .1, force.proper = F, plot = T)
#'

get_beta_params <- function(min.val, mid.val, max.val, levels = c(.025, .975), mid.val.type = c('mean', 'median'), force.proper = T, plot = F) {
	obj.fun <- function(params, q = levels, intended = c(min.val, mid.val, max.val)) {
		a = exp(params[1])
		b = exp(params[2])

		med.pred <- if (mid.val.type[1] == 'mean') a / (a + b) else if (mid.val.type[1] == 'median') qbeta(.5, a, b)

		predicted <- c(qbeta(q, a, b), med.pred) %>% log

		res <- sum((log(intended) - predicted)^2)

		if (is.finite(res) & !is.nan(res) & if (force.proper) (a >= 1 & b >= 1) else T) res else 10000
	}

	res <- optim(par = c(0, 0), fn = obj.fun, method = 'Nelder-Mead')

	list(
		alpha = res$par[1] %>% exp,
		beta = res$par[2] %>% exp,
		distr.parameters = evaluate.beta.params(res$par[1] %>% exp, res$par[2] %>% exp, plot = plot),
		opt.results = res
	)
}
