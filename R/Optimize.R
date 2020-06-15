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
#' @param w.max Maximum number of tests which are considered affordable for the
#'   study.
#' @param s.max Maximum number of individual samples in a pool. To be chosen
#'   after internal tests to evaluate possible loss in sensitivity due to
#'   pooling.
#' @param p Maximal expected prevalence in the target population. It is better
#'   to chose a value higher than what it is really expected in order not to
#'   risk to design an underpowered study.
#' @param n.max Total number of individual samples that can be collected
#' @param w.min Minimum number of tests which is possible to perform.
#' @param s.min Minimum pool size for pooling.
#' @param max.groups The maximum number of rules that should be identified (the
#'   actual number may be lower).
#' @param score.fun The function that computes the score (lower is better). It
#'   takes as argument grid, wich is the simulated results plus enrichment via
#'   \code{enrich_simulation()}.
#' @param a The \eqn{alpha} parameter of the prior distribution.
#' @param b The \eqn{beta} parameter of the prior distribution.
#' @param ... Arguments to passed to \code{ctree()}, to control the creation of
#'   the rules
#'
#' @return
#' @export
#'
#' @examples
#'
#' # Explore possible test design given a maximum pool size of 12, 2000 maximum
#' # number of tests and 2000 sampled individuals and a 5% expected prevalence
#'
#' grid <- design_optimization(s.max = 12, w.max = 2000, p = .05, n.max = 2000)
#'
#' plot_optimization_grid(grid)

design_optimization <- function(s.max = 32, w.max, p, n.max = w.max * s.max,
																w.min = 1, s.min = 1, max.groups = 50, score.fun = function(grid) {
																	with(
																		grid,
																		log(err/p) + log(unc/p) + log(base.unc/p) + log(abs(base.prev - p)/p) + log(w/s)
																	)
																}, a = 0.3, b = 0.3, ...) {

	grid <- expand.grid(s = s.min:s.max, w = w.min:w.max, p = p) %>%
		dplyr::mutate(cases = round(p * w * s)) %>%
		dplyr::filter(w * s <= n.max) %>% {data.frame(., get_estimates(w = .$w, s = .$s, p = .$p, a = a, b = b)$estimates)} %>%
		enrich_simulation(a = a, b = b)

	grid$score <- score.fun(grid)

	if (any(grid$score == -Inf)) warning('-Inf in the optimization score, due to zero estimation error. Put a score 10 times greater than the minimum finite score. Probably you passed an extreme design to the algotrithm. Results may be biased.')

	grid$score[grid$score == -Inf] <- 10 * min(grid$score[is.finite(grid$score)])

	mod <- ctree(score ~ w + s, grid, minbucket = round(nrow(grid)/max.groups), ...)

	if (length(mod) == 1) {
		#grid$mean.score <- mean(grid$score)
		stop('No optimization was found')
	} else {
		grid$mean.score <- predict(mod, newdata = grid)
	}

	dplyr::group_by(grid, mean.score) %>%
		dplyr::mutate(
			dplyr::across(c(s,w), list(group.max = max, group.min = min)),
			dplyr::across(matches('s_'), ~ replace(.x, .x %in% c(s.min, s.max), NA)),
			dplyr::across(dplyr::matches('w_'), ~ replace(.x, .x %in% c(w.min, w.max), NA)),
			s.rule = paste(s_group.min, '<= s <=', s_group.max),
			w.rule = paste(w_group.min, '<= w <=', w_group.max),
			dplyr::across(s.rule:w.rule, ~ stringr::str_remove_all(.x, 'NA <= | <= NA') %>% stringr::str_remove('^\\w$')),
			rule = paste(s.rule, '&', w.rule) %>% stringr::str_remove_all('^ & | & $'),
		) %>%
		dplyr::arrange(mean.score) %>%
		dplyr::mutate(rule = factor(rule, levels = unique(rule))) %>%
		dplyr::select(!(s_group.max:w.rule)) %>%
		dplyr::ungroup()
}

#' Plot the optimization grid
#'
#' Take as input the grid of scores and the optimization windows produced by the
#' optimization algorithm and produce a plot which represent the possible
#' designs defined by pool size \code{s} and number of tested pools \code{w}.
#' The plot is color coded to show the optimization window rules with the
#' relative mean score.
#'
#' @param grid The grid produced by \code{design_optimization()}.
#'
#' @return A ggplot2 object
#' @export
#'
#' @importFrom scales pretty_breaks
#'
#' @examples
#'
#' grid <- design_optimization(w.max = 2000, p = .05, n.max = 2000)
#'
#' plot_optimization_grid(grid)

plot_optimization_grid <- function(grid) {
	dplyr::mutate(grid,
				 rule = sprintf('%s (%2g)', rule, round(mean.score, 2)),
				 rule = factor(rule, levels = unique(rule))
	) %>%
		ggplot(aes(s, w, color = rule)) +
		geom_point(size = 3, stroke = 0, shape = 16) +
		scale_color_discrete('Rules (score)') +
		scale_y_continuous(breaks = scales::pretty_breaks(10)) +
		scale_x_continuous(breaks = function(x) round(x[1]):round(x[2])) +
		labs(x = 'Pool size', y = 'Number of pools') +
		options("PooledPrevalence.ggtheme") +
		theme(legend.position = 'left') +
		guides(color = guide_legend(ncol = 1))
}

#
# summarise_optimization_grid <- function(grid) {
#
# }


