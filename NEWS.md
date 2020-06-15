# PooledTest 0.1

* First commit!

# PooledTest 1.0

* Added optimization algorithm.
* Added hierarchical Bayesian estimation.
* Added proper Readme.
* Numberous bug fixing and implementation of helpers.

# PooledTest 1.0.1

* Fixed a bug into evaluate_beta_params(). cut.plot.at argument was not implemented.

# PooledTest 1.0.2

* Added ... to get_beta_params() to pass arguments to evaluate_beta_params()
* Fixed a bug into evaluate_beta_params() that was preventing the visualization or vertical lines on the distribution boundaries
* Fixed a bug into simulate_pool_test(), the estimation method was not really passed to get_estimates() and the previous method labeling was used.
* In extreme cases (e.g. prevalence = 50%) design_optimization() may give -Inf scores which crash the decision tree algorithm. Set those values to 10 times the minum value. Triggered also a warning.
* If no optimization windows are found by design_optimization() a meaningfull error is raised. 
