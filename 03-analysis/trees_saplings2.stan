
data {
  int N; // number of observations
  int P; // number of plots
  int L; // number of distinct locations within each plot
  int G; // number of groups (fuel models, drought mortality classes, etc.)
  
  array[N] int Y; // observed responses
  array[N] int plot_id; // plot IDs
  array[N] int location_id; // location IDs
  array[P] int group_id; // group IDs
  
  vector[2] coords[L]; // plot-relative coordinates of each subsample
}

parameters {
  vector[G] beta; // fixed effect coefficients 
  matrix[L,P] zGP; // standard-normal deviates for GP effect for each loction,plot,group
  
  array[G] real<lower=0> alpha; // standard deviations of GPs
  array[G] real<lower=0> rho; // length-scales of GPs
  
  array[G] real<lower=0> sigmaPlot; // standard deviation of plot effect
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
  vector[N] logmu; // linear predictor
  matrix[L,P] GP; // actual GP effects for each plot:location
  matrix[L,L] K[G]; // covariance matrix for within-plot GP (prefire)
  matrix[L,L] KL[G]; // cholesky decomposition of K_pre
  
  // priors
  beta ~ normal(-1, 0.5);
  sigmaPlot ~ normal(0, 0.5);
  alpha ~ normal(0,0.5);
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
    logmu[i] = 
      beta[group_id[plot_id[i]]] +
      GP[location_id[i], plot_id[i]] +
      plot_effect[plot_id[i]];
  }
  
  // likelihood
  Y ~ poisson_log(logmu);
}



