#' Prevalence estimation from a pooled study
#'
#' The function allow to estimate prevalence from the results of prevalence
#' study. It use three statistical frameworks: Maximum Likelihood estimation
#' (ML), Conjugate Bayesian (CB) estimation, Hierarchical MCMC Bayesian
#' estimation (HB).
#'
#' The CB and HB estimation allows the setting of priors, the first on the study
#' results, the second on the prevalence estimates themselves; the default prior
#' is a Beta(0.3, 0.3) distribution to account for the positive bias in the
#' prevalence estimation in case of low sample size. The HB method allows also
#' to take in account the sensitivity of the sample acquisition procedure,
#' distibuted as ~ Beta(asens, bsens). The default parameters describe a mean
#' sensitivity of 95% with a lower bound near 70%.
#'
#' The function takes in input the pool size, the number of tested pools and the
#' positive pools. In alternative, to the actual test results, users can provide
#' only the probability of positive pools in the study or the theoretical mean
#' prevalence; the expected number of positive pools is then generated throgh
#' Maximum Likelihood. All parameters study parameters \code{s, w, p_tests, p}
#' can be vectors; this permit to analyze multiple (real or simulated) studies
#' at once.
#'
#' The HB method requires the compilation of a STAN model at the first use. The
#' model will be compiled only the first time.
#'
#' @param s A vector of pool sizs (number of individual samples in a pool)
#' @param w A vector of numbers of tested pools in a study
#' @param k A vector of numbers of positive pools in a study
#' @param p_test A vector of probabilities of positive pools in a study. if
#'   \code{k} is not provided, it used to generate the expected number of
#'   positive pools
#' @param p A vector of prevalences of the disease. if \code{k} and
#'   \code{p_test} is not provided, it used to generate the expected number of
#'   positive pools
#' @param level The uncertainty level to use for the confidence/credibility
#'   interval. Must be strictly greater than 0 and less than 1. Defaults to
#'   0.95, which corresponds to a 95 percent confidence interval
#' @param method Either Maximum Likelihood Estimation (ML), Conjugate Bayesian
#'   (CB) estimation, Hierarchical MCMC Bayesian estimation (HB)
#' @param a \eqn{\alpha} parameter of the Beta prior. The prior is on the test
#'   results in case of \code{method = 'CB'} and on the prevalence estimation in
#'   case of \code{method = 'HB'}
#' @param b \eqn{\beta} parameter of the Beta prior. The prior is on the test
#'   results in case of \code{method = 'CB'} and on the prevalence estimation in
#'   case of \code{method = 'HB'}
#' @param iters Number of samples to draw when using methods CB e HB.
# @param hc.chains number of chains for the Hierarchical Bayesian method. By
#   default it's the number of cores determined by
#'   \code{parallel::detectCores()} or as stored in \code{getOption("mc.cores")}
#' @param asens \eqn{\alpha} parameter of the Beta distribution describing the
#'   sample acquisition sensitivity
#' @param bsens \eqn{\beta} parameter of the Beta distribution describing the
#'   sample acquisition sensitivity
#' @param consider.sensitivity whether to consider sample acquisition
#'   sensitivity in the HB estimation
#'
#' @return A list with two components: A named numeric matrix of height equal to
#'   the length of the parameter vectors. The matrix contains:
#'   \describe{\item{p_est}{The positivity rate for the pools in the study}
#'   \item{k}{The number of positive pools} \item{Est}{The estimated disease
#'   prevalence} \item{Lo}{The lower uncertainty bound} \item{Up}{The upper
#'   uncertainty bound} } For methods CB and HB, if just one set of parameters
#'   is passed, a data frame with the posterior samples is returned.
#'
#' @export
#'
#' @examples
#'
#' # Estimate prevalence from a study with 30 positive pools out of 200, with 10 samples per pool
#'
#' get_estimates(s = 10, w = 200, k = 30)$estimates
#'
#' # Same study now analyzed with a full bayesian estimation considering sampling sensitivity.
#' # Notice the slightly higher estimated prevalence, due to the incorporation of the false
#' # negatives due to sampling
#'
#' get_estimates(s = 10, w = 200, k = 30, method = 'HB')$estimates
#'
#' ## Through the Bayesian method is possible to perform Bayesian hypothesis testing
#'
#' # Which is the probability that the estimated prevalence from the previous example is
#' # lower than 2%?
#'
#' sampled.prev <- get_estimates(s = 10, w = 200, k = 30, iters = 20000)$samples
#'
#' mean(sampled.prev <= 0.02)
#'
#' # Which is the probability that the prevalence in region A is higher than that in region B,
#' # given just 5 positive pools more?
#'
#' sampled.prev.A <- get_estimates(s = 10, w = 200, k = 35)$samples
#'
#' sampled.prev.B <- get_estimates(s = 10, w = 200, k = 30)$samples
#'
#' mean(sampled.prev.A >= sampled.prev.B)
#'
#'


get_estimates <- function(s, w, k = NULL, p_test = NULL, p = NULL, level = .95,
													method = c('CB', 'HB', 'ML'), a = .3, b = .3, iters = 20000,
													#hb.chains = getOption("mc.cores", parallel::detectCores()),
													consider.sensitivity = T, asens = 8.88, bsens = 0.74) {

	if (is.null(p) & is.null(p_test) & is.null(k)) stop('Provide k, p or p_test')

	if (is.null(p_test) & is.null(k) & !is.null(p)) p_test <- 1 - (1 - p)^s

	if (is.null(p_test) & !is.null(k)) p_test <- round(k) / round(w)

	if (!is.null(p_test) & is.null(k)) k <- round(w) * p_test

	k <- round(k)
	w <- round(w)

	p_test <- k / w

	samples <- NULL

	if (method[1] == 'ML') {
		Z <- qnorm((1 - level)/2) %>% abs

		Up <- 1 - (1 - invlogit(logit(p_test) + Z * sqrt(1/(w * p_test * (1 - p_test)))))^(1/s)
		Lo <- 1 - (1 - invlogit(logit(p_test) - Z * sqrt(1/(w * p_test * (1 - p_test)))))^(1/s)

		Est <-  1 - (1 - p_test)^(1/s)

	}
	else if (method[1] == 'CB') {

		ql <- (1 - level)/2
		qu <- level + ql

		# a <- (1 - (1 - a / (a + b))^s) * (a + b)
		# b <- (1 - (1 - b / (a + b))^s) * (a + b)

		#browser()

		Up <- 1 - (1 - qbeta(qu, a + k, b + w - k))^(1/s)
		Lo <- 1 - (1 - qbeta(ql, a + k, b + w - k))^(1/s)
		Est <- 1 - (1 - qbeta(.5, a + k, b + w - k))^(1/s)

		if (length(Est) == 1) {
			samples <- data.frame(p_test = rbeta(iters, a + k, b + w - k))
			samples <- 1 - (1 - samples$p_test)^(1/s)
		}
	}
	else if (method[1] == 'HB') {

		#stop('Hierarchical Bayes method will be implemented soon!')
		if (nrow(data.frame(w, s, k)) > 1) stop('Only one set of parameters can be used with method HB')
		# process.args <- list(
		# 	stan_data = list(
		# 		w = w, # number of pools
		# 		k = k, # number of positive pools
		# 		s = s, # number of samples in each pool
		# 		a = a, b = b, # prior parameters for p_sample
		# 		asens = asens, bsens = bsens # prior for swabbing sensitivity
		# 	),
		# 	chains =  hc.chains
		# )

		fit.args <- list(
				w = w, # number of pools
				k = k, # number of positive pools
				s = s, # number of samples in each pool
				a = a, b = b, # prior parameters for p_sample
				asens = asens, bsens = bsens # prior for swabbing sensitivity
			)

		# if (!file.exists('Estimation_model.rds')) compile.estimation.model()
		#
		# fit <- callr::r(function(x) {
		# 	rstan::sampling(object = readRDS('Estimation_model.rds'), data = x$stan_data, iter = iters, chains = x$chains, cores = x$chains)
		# }, args = process.args)
		#
		# Estimates <- broom::tidy(conf.int = T, conf.level = level, estimate.method = "median") %>% dplyr::filter(term == 'p_sample')

		lk.fun <- function(p_sample, args) {

			if (p_sample < 0 | p_sample > 1) return(-Inf)

			if (consider.sensitivity) p_sample <- p_sample * rbeta(1,  args$asens,  args$bsens)

			p_test <- 1 - dbinom(0, args$s, p_sample)

			dbinom(args$k, args$w, p_test, log = T) + dbeta(p_sample, args$a, args$b, log = T)
		}

		ql <- (1 - level)/2
		qu <- level + ql


		fit <- MCMCmetrop1R(lk.fun,
								 theta.init = 1 - (1 - 30/200)^(1/10),
								 args = fit.args,
								 mcmc = iters
								 )

		# Est = Estimates$estimate
		# Lo = Estimates$conf.low
		# Up = Estimates$conf.high

		Est <- quantile(fit, .5)
		Lo <- quantile(fit, ql)
		Up <- quantile(fit, qu)

		# if (length(Est) == 1) {
		# 	samples <- rstan::extract(fit, pars = c('p_sample', 'p_test', 'sens')) %>% as.data.frame()
		# }

		if (length(Est) == 1) {
			samples <- fit[,1]
		}
	}
	else {
		stop('Unknown "method" value')
	}

	list(
		samples = samples,
		estimates = data.frame(p_test, k, Est, Lo, Up)
	)
}
