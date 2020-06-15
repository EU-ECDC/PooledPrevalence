#' Simulation of pooled prevalence studies
#'
#' Generate a sample of individuals with a given disease status informed by a
#' chosen underlying diseases prevalence and risk variability. Then simulate a
#' pooled test results from such sample. The process can be repeated for a
#' number of iteration to have a distribution of results. The prevalence in each
#' simulation is estimated through the \code{get_estimates(method = 'CB')}
#' function.
#'
#' @param w Number of pools to test (default \code{200})
#' @param s Number of sampled individuals per pool (default \code{12})
#' @param p Expected population prevalence (default \code{0.01})
#' @param iters Number of simulations (default \code{3000})
#' @param consider.sensitivity Consider sensitivity in sample acquisition in the
#'   simulation (default \code{FALSE})
#' @param asens \eqn{\alpha} Parameter of the Beta distribution for the sampling
#'   sensitivity (default \code{8.88})
#' @param bsens \eqn{\beta} Parameter of the Beta distribution for the sampling
#'   sensitivity (default \code{0.74})
#' @param simulate.results Add estimation results to the simulation (default
#'   \code{TRUE})
#' @param estimation.method Whether to use Bayesian Conjugate \\(default\\) or
#'   MLE for the prevalence estimation (default \code{'CB'})
#'
#' @return A dataframe with: the initial w, s, and p parameters; p.sample: the
#'   individual risk ditribution quantiles (0.025, .5, .975); cases: the number
#'   of positive individuals in the sample, after accounting for sampling
#'   sensitivity; pos & neg: the number of positive and negative pools out of w
#'   tested pools.
#' @export
#' @importFrom magrittr set_colnames
#'
#' @examples
#' # Simulate a study in a population with a disease risk of 1%
#'
#' set.seed(1234)
#' simulate_pool_test(s = 8, w = 250, p = .01, iter = 1, consider.sensitivity = FALSE)
#'
#' # Same simulation this time considering false negatives due to sampling sensitivity,
#' # with \code{asens} and \code{bsens} chosen to have a mean sensitivity of 95% and
#' # a lower bound around 70%
#'
#' set.seed(1234)
#' simulate_pool_test(s = 8, w = 250, p = .01, iter = 1, consider.sensitivity = TRUE)
#'
#' # use the \code{iters} argument to create more simulations
#'
#' simulate_pool_test(s = 8, w = 250, p = .01, iter = 3000, consider.sensitivity = FALSE)
#'

simulate_pool_test <- function(s = 12, w = 200, p = .01, iters = 3000, consider.sensitivity = FALSE, simulate.results = TRUE, asens = 8.88, bsens = 0.74, estimation.method = 'CB') {
	if (estimation.method %nin% c('CB', 'ML')) stop('Only Bayesian Conjugate and ML estimation methods are allowed')

	u <- 80 # Strangely, it does not seem to impact the simulations. to investigate

	a <- p * u
	b <- u - a
	sim.p <- rbeta(w * s * iters, a, b)
	cases <- rbinom(w * s * iters, 1, sim.p)
	positives <- if (consider.sensitivity) cases * rbinom(w * s * iters, 1, rbeta(w * s * iters, asens, bsens)) else cases

	pool.res <- (matrix(positives, nrow = s) %>% Rfast::colsums()) > 0
	test.pos.pools <- matrix(pool.res, nrow = w) %>% Rfast::colsums()

	test.cases <- matrix(positives, ncol = iters) %>% Rfast::colsums()
	p.distr <- matrix(sim.p, ncol = iters) %>% apply(2, quantile, c(.025, .5, .975)) %>% t %>% set_colnames(c('p.sample.l', 'p.sample', 'p.sample.u'))

	res <- data.frame(w, s, p, p.distr, cases = test.cases, pos = test.pos.pools, neg = w - test.pos.pools)

	if (simulate.results) {
		p.est <- get_estimates(s = s, w = w, k = test.pos.pools, method = estimation.method)$estimates

		p.est <- p.est[,colnames(p.est) != 'k']

		res <- data.frame(res, p.est)
	}

	res
}

#' Formats the results of a simulation
#'
#' Adds more readable names, transforms estimates to percentages and adds uncertainty
#' intervals around the simulated parameters. If the simulation are not already enriched
#' via \code{enrich_simulation} it does it automatically.
#'
#' @param simulations The results of a call to \code{simulate_pool_test}
#' @param level The uncertainty level around the estimates.
#' @param reverse Whether to show the results in vertical format instead of wide.
#'
#' @return A data frame with reformatted titles
#' @export
#'
#' @examples
#'
#' sim <- simulate_pool_test(s = 8, w = 250, p = .01, iter = 3000, consider.sensitivity = FALSE)
#'
#' formatted.sim <- format_simulation(sim)
#'

format_simulation <- function(simulations, level = .95, reverse = T) {
	lvl <- level

	if ('unc' %nin% colnames(simulations)) simulations <- enrich_simulation(simulations)

	out <- simulations %>%
		dplyr::select(
			'Estimated prevalence' = Est,
			'Estimation uncertainty' = unc,
			'Estimation error' = err,
			'Unpooled study prevalence' = base.prev,
			'Unpooled study uncertainty' = base.unc,
			'Unpooled study error' = base.err
		) %>%
		dplyr::summarise(dplyr::across(dplyr::everything(), sumfun, level = lvl))

	if (reverse) t(out) else out
}

#' Enrich simulated estimations with accuracy information
#'
#' Given some simulated estimates, compute the uncertainty interval size and
#' estimation error given a expected real prevalence. It is also computed the
#' expected prevalence in an unpooled study, plus the relative estimation error
#' and uncertainty size. The unpooled prevalence is computed using the Conjugate
#' Bayes method, with Beta prior parameters both equal to 0.3
#'
#' @param simulations A data frame with as columns the pool size \code{s}, the
#'   number of pools \code{w}, the number of positive cases \code{cases}, the
#'   prevalence estimate \code{Est}, the upper \code{Up} and lower \code{Lo}
#'   uncertainty the expected real prevalence \code{p}.
#' @param a \eqn{\alpha} parameter of the prior Beta distribution for the
#'   unpooled estimate.
#' @param b \eqn{\beta} parameter of the prior Beta distribution for the
#'   unpooled estimate.
#'
#' @return A data frame with added the uncertainty size \code{unc} and the
#'   estimation error \code{err}. Also the estimated prevalence in an unpooled
#'   study \code{base.prev} is reported, together with the unpooled uncertainty
#'   \code{base.unc} and error \code{base.err}.
#'
#' @export

enrich_simulation <- function(simulations, a = 0.3, b = 0.3) {

	unpooled.a <- a + simulations$cases
	unpooled.b <- with(simulations, b + w * s - cases)

	dplyr::mutate(simulations,
								unc = Up - Lo,
								err = abs(Est - p),
								base.prev = qbeta(.5, unpooled.a, unpooled.b),
								base.unc = qbeta(.975, unpooled.a, unpooled.b) - qbeta(.025, unpooled.a, unpooled.b),
								base.err = abs(base.prev - p)
	)
}
