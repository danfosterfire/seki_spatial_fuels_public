data{
  int<lower=1> N; //total number of observations
  int<lower=1> L; // sampling locations per plot
  int<lower=1> P; //number of plots
  int<lower=1> G; // number of groups
  int<lower=1, upper=P> plotID[N];
  int<lower=1, upper=L> locID[N];
  int<lower=1, upper=G> groupID[N];
  int<lower=0> y[N]; //response variable 
  vector[2] coords[L]; // coordinates
  //int Nsim; // number of distances to simulate cov.
  //vector[1] simDist[Nsim]; // distance vector to simulate covariance
}
parameters{
  vector[G] a0;
  real<lower=0> kappa;
  //hyperparameters for ap
  real<lower=0> sigma_ap;
  vector[P] z_ap;
  //GP parameters
  matrix[L,P] zGP; 
  real<lower=0> alpha;
  real<lower=0> rho;
}
transformed parameters{ 
  vector[N] logmu; 
  vector[P] ap;
  matrix[L,L] K;
  matrix[L,P] GP;
  
  // Define the GP matrix
  {
    matrix[L,L] KL;
    K = gp_exp_quad_cov(coords, alpha, rho); //cov_exp_quad
    KL = cholesky_decompose(add_diag(K, 1e-9));
    GP = KL * zGP;
  }
  

  //plot-level varying intercept, non-centered
  ap = sigma_ap*z_ap;

  //main model
  for(i in 1:N){
    logmu[i] = a0[groupID[i]] + ap[plotID[i]] + GP[locID[i], plotID[i]]; 
    }
}
model{
  //priors
  a0 ~ normal(1.5, 1);
  kappa ~ normal(0, 2);
  sigma_ap~ normal(0, 2);
  z_ap~ normal(0, 1);
  to_vector(zGP) ~ normal(0,1);
  alpha ~ normal(0, 0.5); //alpha ~ inv_gamma(5, 2);
  rho ~ inv_gamma(5,40); //rho ~ inv_gamma(3, 8);
  
  //likelihood
  y ~ neg_binomial_2_log(logmu, kappa);
}
generated quantities{
  vector[N] y_rep; 
  for(i in 1:N){
    y_rep[i] = neg_binomial_2_log_rng(logmu[i], kappa);
  }
}
