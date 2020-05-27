#' Evaluate the performance of alternative pooled testing designs
#'
#' This tool conduct a grid search among possible pairs of pool sizes (number of
#' samples per pool) and total number of tested pools, given some operational
#' limits and an hypothetical underlying disease prevalence. Then for each pair
#' the algorithm simulate a result using ML and back-estimate the prevalence.
#' From the estimates the function computes a log score based on the estimation
#' error and the estimation uncertainty both for the pooled and unpooled (pool
#' size = 1) estimation. The score is also increased by a term given by the
#' ratio of the number of pools to the pool size. This last term act as balance
#' to precise estimations and lab resource usage. Finally a decision tree model
#' partitions the simulated solution in optimization windows defined by a range
#' of pool sizes and total number of tests. These windows can help the users to
#' choose the best study design for their specific setting.
#'
#' The suggested optimal designs depend strongly on the expected disease
#' prevalence. In general a lower prevalence allows to have larger pool sizes
#' without loss of estimation accuracy and therefore larger saving of testing
#' resources. It is suggested to choose a prevalence higher than the expected
#' one (e.g. slightly higher than the expected higher possible value) in order
#' not to end up with an underpowered study in case the prevalence is higher
#' than expected,
#'
#'
#' @param w.max Maximum number of tests which are consiedered affordable for the study.
#' @param s.max Maximum number of individual samples in a pool. To be chosen
#'   after internal tests to evaluate possible loss in sensitivity due to
#'   pooling.
#' @param p Hypothisized prevalence in the target population. It is better to
#'   chose a value higher than what it is really expected.
#' @param n.max Total number of individual samples that can be collected
#' @param w.min Minimum number of tests which is possible to perform.
#' @param s.min Minimum pool size for pooling.
#'
#' @return
#' @export
#' @importFrom stats qbeta
#'
#' @examples
design_optimization <- function(w.max = 100, s.max = 32, p, n.max = w.max * s.max, w.min = 1, s.min = 1) {
	grid <- expand.grid(s = s.min:s.max, w = w.min:w.max, p = p) %>%
		dplyr::filter(w * s <= n.max) %>% {data.frame(., get_estimates(w = .$w, s = .$s, p = .$p))} %>%
		dplyr::mutate(
			unc = Up - Lo,
			err = abs(Est - p),
			base.prev = qbeta(.5, .05 + p * w * s, .05 + (1 - p) * w * s),
			base.unc = qbeta(.975, .05 + p * w * s, .05 + (1 - p) * w * s) - qbeta(.025, .05 + p * w * s, .05 + (1 - p) * w * s),
			d = log(err/p) + log(unc/p) + log(base.unc/p) + log(abs(base.prev - p)/p) + log(w/s), acc = T)
}
