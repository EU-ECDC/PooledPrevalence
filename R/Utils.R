percent <- function(x) {
	#sapply(x, function(x) if (!is.na(x)) {if (abs(x * 100) < 1) sprintf('%.2g%%', x * 100) else sprintf('%.3g%%', x * 100)} else NA)
	sapply(x, function(x) if (!is.na(x)) {if (abs(x * 100) < 1) paste0(signif(x * 100, 2), '%') else paste0(signif(x * 100, 3), '%')} else NA)
}

'%nin%' <- function(x, y) {
	!(x %in% y)
}

invlogit <- function(x)
{
	1/(1 + exp(-x))
}

logit <- function(x)
{
	log(x/(1 - x))
}

my_round <- function(z) {
	dplyr::if_else(z < 1, signif(z, 2), round(z))
}


sumfun <- function(x, level = .95) {

	#, mode: {names(table(x))[which.max(table(x))]}
	sprintf('%s, %s%%CI: [%s, %s]', my_round(median(x)), 100 * level, my_round(quantile(x, .5 - level/2)), my_round(quantile(x, .5 + level/2))) %>%
		stringr::str_replace_all('0\\.\\d+', function(x) percent(as.numeric(x)))
}

# compile.estimation.model <- function() {
# 	invisible(rstan::stan_model('Estimation_model.stan', save_dso = T, auto_write = T))
# }
