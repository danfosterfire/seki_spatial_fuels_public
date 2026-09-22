//

data {
  int<lower=1> N; // number of unique particles whose diameter we have
  int<lower=1> P; // number of unique plots
  int<lower=1> G; // number of unique groups
  int<lower=1> TS; // number of unique transects
  
  int<lower=1> group_id[P]; // group ids for each plot
  int<lower=1> plot_id[TS]; // plot ids
  int<lower=1> transect_id[N]; // transect ids
  
  real<lower=0> lb; // lower bound for truncation of size
  real<lower=lb> ub; // upper bound for truncation
  
  real Y[N]; // observed diameters in cm
}

parameters {
  array[G] real intercept;
  array[G] real<lower=0> phi;
  
  array[G] real<lower=0> SD_plot;
  vector[P] z_plot;
  
  array[G] real<lower=0> SD_trans;
  vector[TS] z_trans;
}

transformed parameters {
  vector[P] effect_plot;
  vector[TS] effect_trans;
  for(p in 1:P){
    effect_plot[p] = SD_plot[group_id[p]] * z_plot[p];
  }
  for (ts in 1:TS){
    effect_trans[ts] = SD_trans[group_id[plot_id[ts]]]*z_trans[ts];
  }
}

model {
  // variable declarations
  vector[N] shape;
  vector[N] scale;
  vector[N] logmu;
  
  // priors (uniform omitted)
  intercept ~ normal(0, 1);
  SD_plot ~ normal(0, 0.25);
  SD_trans ~ normal(0, 0.25);
  phi ~ normal(0, 5);
  // random effect realizations
  z_plot ~ std_normal();
  z_trans ~ std_normal();
  
  // linear predictor and reprarameterize
  for (i in 1:N){
    logmu[i] = 
      intercept[group_id[plot_id[transect_id[i]]]] + 
      effect_plot[plot_id[transect_id[i]]] +
      effect_trans[transect_id[i]];
    
    shape[i] = phi[group_id[plot_id[transect_id[i]]]]+2;
    
    scale[i] = exp(logmu[i]) * (1+phi[group_id[plot_id[transect_id[i]]]]);
  }

  // likelihood
  for (i in 1:N){
    Y[i] ~ inv_gamma(shape[i], scale[i])T[lb,ub];
  }
}

