percent <- function(x) {
	#sapply(x, function(x) if (!is.na(x)) {if (abs(x * 100) < 1) sprintf('%.2g%%', x * 100) else sprintf('%.3g%%', x * 100)} else NA)
	sapply(x, function(x) if (!is.na(x)) {if (abs(x * 100) < 1) glue('{signif(x * 100, 2)}%') else glue('{signif(x * 100, 3)}%')} else NA)
}

'%nin%' <- function(x, y) {
	!(x %in% y)
}
