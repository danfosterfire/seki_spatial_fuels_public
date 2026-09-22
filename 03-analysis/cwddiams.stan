//

data {
  int<lower=1> N; // number of unique particles whose diameter we have
  int<lower=1> P; // number of unique plots
  int<lower=1> G; // number of unique groups
  
  int<lower=1> group_id[P]; // group ids for each plot
  int<lower=1> plot_id[N]; // plot ids
  
  real<lower=0> lb; // lower bound for truncation of size
  real<lower=lb> ub; // upper bound for truncation
  
  real Y[N]; // observed diameters in cm
}


parameters {
  array[G] real<lower = lb, upper = ub> mu; // this needs a constraint; or should it be on logmu?
  array[G] real<lower=0> phi;
  
  //array[G] real<lower=0> SD_plot;
  //vector[P] z_plot;
  
}

//transformed parameters {
//  vector[P] effect_plot;
//  for(p in 1:P){
//    effect_plot[p] = SD_plot[group_id[p]] * z_plot[p];
//  }
//}

model {
  // variable declarations
  vector[N] shape;
  vector[N] scale;
  
  // priors (uniform omitted)
  phi ~ normal(0, 10);
  
  // linear predictor and reprarameterize
  for (i in 1:N){
    shape[i] = phi[group_id[plot_id[i]]]+2;
    
    scale[i] = mu[group_id[plot_id[i]]] * (1+phi[group_id[plot_id[i]]]);
  }

  // likelihood
  for (i in 1:N){
    Y[i] ~ inv_gamma(shape[i], scale[i])T[lb,ub];
  }
}

