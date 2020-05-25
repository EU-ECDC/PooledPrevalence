#' Simulation of pooled prevalence studies
#'
#' Generate a sample of individuals with a given disease status informed by a
#' chosen underlying diseases prevalence and risk variability. Then simulate a
#' pooled test results from such sample. The process can be repeated for a
#' number of iteration to have a distribution of results.
#'
#' @param w number of pools to test (default \code{200})
#' @param s number of sampled individuals per pool (default \code{12})
#' @param p expected population prevalence (default \code{0.01})
#' @param iters number of simulations (default \code{3000})
#' @param consider.sensitivity consider sensitivity in sample acquisition in the
#'   simulation (default \code{FALSE})
#' @param asens \alpha parameter of the Beta distribution for the sampling
#'   sensitivity (default \code{8.88})
#' @param bsens \beta parameter of the Beta distribution for the sampling
#'   sensitivity (default \code{0.74})
#' @param simulate.results add estimation results to the simulation (default
#'   \code{TRUE})
#' @param estimation.method whether to use Bayesian Conjugate \\(default\\) or
#'   MLE for the prevalence estimation (default \code{'B'})
#'
#' @return A dataframe with: the initial w, s, and p parameters; p.sample: the
#'   individual risk ditribution quantiles (0.025, .5, .975); cases: the number
#'   of positive individuals in the sample, after accounting for sampling
#'   sensitivity; pos & neg: the number of positive and negative pools out of w
#'   tested pools.
#'
#' @examples
#' # Simulate a study in a population with a disease risk of 1%
#'
#' set.seed(1234)
#' simulate.pool.test(s = 8, w = 250, p = .01, iter = 1, consider.sensitivity = F)
#'
#' # Same simulation this time considering false negatives due to sampling sensitivity,
#' # with \emph{asens} and \emph{bsens} chosen to have a mean sensitivity of 95% and
#' # a lower bound around 70%
#'
#' set.seed(1234)
#' simulate.pool.test(s = 8, w = 250, p = .01, iter = 1, consider.sensitivity = T)
#'
#' # use the \emph{iters} argument to create more simulations
#'
#' simulate.pool.test(s = 8, w = 250, p = .01, iter = 3000, consider.sensitivity = F)
#'
#' @export

simulate.pool.test <- function(w = 200, s = 12, p = .01, iters = 3000, consider.sensitivity = F, simulate.results = T, asens = 8.88, bsens = 0.74, estimation.method = 'B') {
	if (estimation.method %nin% c('B', 'MLE')) stop('Only Bayesian Conjugate and ML estimation methods are allowed')

	u <- 80 # Strangely, it does not seem to impact the simulations. to investigate

	a <- p * u
	b <- u - a
	sim.p <- rbeta(w * s * iters, a, b)
	cases <- rbernoulli(w * s * iters, sim.p)
	positives <- if (consider.sensitivity) cases * rbernoulli(w * s * iters, rbeta(w * s * iters, asens, bsens)) else cases

	pool.res <- (matrix(positives, nrow = s) %>% colsums()) > 0
	test.pos.pools <- matrix(pool.res, nrow = w) %>% colsums()

	test.cases <- matrix(positives, ncol = iters) %>% colsums()
	p.distr <- matrix(sim.p, ncol = iters) %>% apply(2, quantile, c(.025, .5, .975)) %>% t %>% set_colnames(c('p.sample.l', 'p.sample', 'p.sample.u'))

	res <- data.frame(w, s, p, p.distr, cases = test.cases, pos = test.pos.pools, neg = w - test.pos.pools)

	if (simulate.results) {
		p.est <- get.estimates(s = s, w = w, k = test.pos.pools)
		p.est['k'] <- NULL

		res <- data.frame(res, p.est)
	}

	res
}
