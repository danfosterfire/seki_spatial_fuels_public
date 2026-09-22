
data {
  int<lower=1> N; // number of observations, 1 observation per plot
  int<lower=1> J; // number of fixed effect predictors; here 1 per mort class
  matrix[N,J] X; // fixed effect predictors
  real<lower=0> plot_area_ha; // plot area, used for standardizing the betas to per hectare
  array[N] int Y; // response vector, the count of trees (incl snags) or saplings (incl dead)
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[J] beta; // fixed effect coeffs
  real<lower=0> kappa; // NB dispersion
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // variable declarations
  vector[N] logmu; // linear predictor
  vector[N] mu; // NB mean
  // priors
  beta ~ normal(0, 5);
  kappa ~ normal(0, 5);
  
  // linear predictors
  logmu = X * beta;
  mu = exp(logmu);
  
  // likelihoods
  Y ~ neg_binomial_2(mu*plot_area_ha, kappa);
}











