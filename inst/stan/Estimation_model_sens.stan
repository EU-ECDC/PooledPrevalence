data {
  //declare observed quantities
  int<lower=0> w; //number of pools
  int<lower=0, upper=w> k; //number of positive pools
  int<lower=0> s ; // number of samples in a pool
  real<lower=0> a;
  real<lower=0> b;
  real<lower=0> asens;
  real<lower=0> bsens;
}
transformed data {

}
parameters {
  //declare unobserved quantities
  real<lower=0, upper=1> p_sample; //
  real<lower=0, upper=1> sens;
}
transformed parameters {
  real<lower=0, upper=1> p_test;
  p_test = 1 - exp(binomial_lpmf( 0 | s, p_sample * sens));
}
model {
  p_sample ~ beta(a, b);
  sens ~ beta(asens, bsens);
  //p_test ~ beta(1, 1); // should this be added?
  target += binomial_lpmf( k | w , p_test );
}
generated quantities {

}
