---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# PooledPrevalence

<!-- badges: start -->
<!-- badges: end -->

A set of tools to help design a laboratory based prevalence study harnessing the pooling of laboratory samples to increase resource efficiency.
The tools provide guidance in choosing the number of samples to pool and the total number of tests to perform, focusing on achieving the maximum estimation performance (low uncertainty and estimation error) while keeping the total amount of tests low.
In the package there are functions to simulate study outcomes, optimize the study design, and finally analyze the results of a pooled study to retrieve the relative prevalence estimates.
The package is based on a Bayesian framework, with the possibility to perform a full hierarchical analysis including uncertainty due to the sensitivity in acquiring the testing material.
The package was developed at the occasion of the COVID-19 pandemic, but can be utilized in every laboratory based prevalence study.

!!DISCLAIMER: the Hierarchical Bayes method of estimation is still under development

## Installation

You can install the released version of PooledPrevalence with:

``` r
# install.packages("devtools")
devtools::install_github("EU-ECDC/PooledPrevalence")
```
## Basic usage

The main goal of the package is to estimate prevalence in a population/risk group given the results of a pooled laboratory test. Let's assume a test with the following design: 2000 individual biological samples divided in 10 samples per testing pool (pool size), for a total of 200 pools. In this experiment we observe 30 positive pools (15%). With get_estimates is possible to estimate the underlying prevalence:


```r
library(PooledPrevalence)

results <- get_estimates(10, 200, 30)

print(results$estimates)
#>   p_test  k        Est         Lo         Up
#> 1   0.15 30 0.01610737 0.01103635 0.02251263
```

The object `results` contains both the estimates and, if a Bayesian method is used, the posterior samples of the distribution. The samples allows the inspection of the full posterior distribution of the prevalence. We suggest to use a high number of iteration for visualization purpose when `method = 'BC'` (the default) is used.


```r

results <- get_estimates(10, 200, 30, iters = 200000)

plot(density(results$samples$p_sample))
```

<img src="man/figures/README-plot posterior-1.png" title="plot of chunk plot posterior" alt="plot of chunk plot posterior" width="100%" />

## Hypothesis testing

Having access to the posterior samples allows Bayesian hypothesis testing for decision making. For example a certain region may want to remove a certain control measure only if the average prevalence is below 2% with a > 90% probability:


```r

results <- get_estimates(10, 200, 30, iters = 200000)

mean(results$samples$p_sample < .02)
#> [1] 0.891305
```

In this case the probability of a mean prevalence below 2% is 89.1%, therefore there's not enough evidence in favor of removing the control measure. We may want to increase the sample size of the study or wait more. The Bayesian framework allows to simply add new samples the the already collected ones and simply update the analysis.


```r

old.results <- get_estimates(10, 200, 30, iters = 200000)

new.tested.pools <- 60

new.positve.pools <- 9

new.results <- get_estimates(10, 200 + new.tested.pools, 30 + new.positve.pools, iters = 200000)

mean(new.results$samples$p_sample < .02)
#> [1] 0.921165
```

The new probability is 92.1%, so we have enough evidence to remove the control measure.

Finally we want to ascertain the probability that region A, with 30 positive pools out of 200, has a lower prevalence or region B (38 positive pools out of 180).


```r

region.A <- get_estimates(10, 200, 30, iters = 200000)

region.B <- get_estimates(10, 180, 38, iters = 200000)

mean(region.A$samples$p_sample < region.B$samples$p_sample)
#> [1] 0.939005
```
