#' PooledPrevalence: Tools for Implementation and Analysis of Prevalence Studies
#' based Pooled Testing.
#'
#' A set of tools to help design a laboratory based prevalence study harnessing
#' the pooling of laboratory samples to increase resource efficiency. The tools
#' provide guidance in choosing the number of samples to pool and the total
#' number of tests to perform, focusing on achieve the maximum estimation
#' performance (low uncertainty and estimation error) while keeping the total
#' amount of tests low. In the package there functions to simulate study
#' outcomes, optimize the study design and finally analyze the results of a
#' pooled study to retrieve the relative prevalence estimates. The package is
#' based on a Bayesian framework, with the possibility to perform a full Baysian
#' analysis including also uncertainty due to the sensitivity in acquiring the
#' testing material. The package was developed in the occasion of the COVID-19
#' pandemic, but can be utilized in every laboratory based prevalence study.
#'
#'
#'
#' @importFrom dplyr %>%
#' @importFrom stats median qnorm quantile rbeta rbinom
#' @import ggplot2
#' @importFrom partykit ctree
#'
#' @docType package
#' @name PooledPrevalence
NULL
