
data {
  int N; // number of observations
  int P; // number of plots
  int L; // number of distinct locations within each plot
  int G; // number of groups (fuel models, drought mortality classes, etc.)
  int J; // number of fixed effects covarites, including intercept and dummies 
         // for the G groups, plus any additional variables
  
  int<lower=0,upper=1> Y[N]; // observed responses; 1 if present 0 if absent
  int plot_id[N]; // plot IDs
  int location_id[N]; // location IDs
  int group_id[P]; // group IDs
  
  matrix[N,J] X; // fixed effects design matrix
  
  vector[2] coords[L]; // plot-relative coordinates of each subsample
}

parameters {
  
  vector[J] beta; // fixed effect coefficients 
  matrix[L,P] zGP; // standard-normal deviates for GP effect for each loction,plot,group
  
  real<lower=0> alpha[G]; // standard deviations of GPs
  real<lower=0> rho[G]; // length-scales of GPs
  
  real<lower=0> sigmaPlot[G]; // standard deviation of plot effect
  vector[P] z_plot; // standard normal deviates for plot effect
  
}

transformed parameters {

  vector[P] plot_effect; // atual plot effects
  
  for (i in 1:P){
    plot_effect[i] = sigmaPlot[group_id[i]] * z_plot[i];
  }
}

model {
  // variable declarations
  vector[N] logit_eta; // linear predictor for presence/absence
  matrix[L,P] GP; // actual GP effects for each plot:location
  matrix[L,L] K[G]; // covariance matrix for within-plot GP (prefire)
  matrix[L,L] KL[G]; // cholesky decomposition of K_pre
  
  
  
  // priors
  beta ~ normal(0, 2);
  sigmaPlot ~ normal(0, 0.5);
  alpha ~ normal(0, 0.5);
  rho ~ inv_gamma(5,40);
  
  // random effect realizations
  z_plot ~ std_normal();
  to_vector(zGP) ~ std_normal();
  
  // gaussian processses
  for (group in 1:G){
    K[group] = gp_exp_quad_cov(coords, alpha[group], rho[group]);
    K[group] = add_diag(K[group], 1e-9);
    KL[group] = cholesky_decompose(K[group]);
    
  }
  for (plot in 1:P){
    GP[,plot] = KL[group_id[plot]] * zGP[,plot];
  }
  
  // linear predictors
  for (i in 1:N){
    logit_eta[i] = 
      X[i,] * beta +
      GP[location_id[i], plot_id[i]] +
      plot_effect[plot_id[i]];
      
  }
  
  // likelihood
  Y ~ bernoulli_logit(logit_eta);
  
}



