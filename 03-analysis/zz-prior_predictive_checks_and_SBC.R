
#### setup #####################################################################

library(here)
library(tidyverse)



#### function definitions ######################################################


simulate_data = 
  function(id){
    
    # load in the observed X data
    obs_data = readRDS(here::here('02-data',
                                  '05-for_analysis',
                                  'litterduff_data.rds'))
    
    # pull parameter values from the prior distribution
    # NOTE: because we're interested in the GP parameters for this study,
    # not as concerned with doing SBC for the fixed effects and plot SD; 
    # just setting these rather than drawing them from the prior facilitates 
    # creating realistic datasets to do SBC with
    #beta = rnorm(n = 6, mean = 0, sd = 2)
    #sigmaPlot = truncnorm::rtruncnorm(n = 1, a = 0, mean = 0, sd =0.5)
    beta = c(2, 0.5, -0.5)
    sigmaPlot = 0.25
    kappa = rcauchy(n = 1, location = 0, scale = 5)
    while (kappa <= 0){
      kappa = rcauchy(n = 1, location = 0, scale = 5)
    }
    alpha = truncnorm::rtruncnorm(n = obs_data$G, a = 0, mean = 0, sd = 0.5)
    rho = invgamma::rinvgamma(n = obs_data$G, shape = 5, rate = 40)
    
    # fixed effects
    XB = obs_data$X %*% beta
    
    # GP realizations
    D = dist(obs_data$coords) %>% as.matrix()
    Sigma =
      lapply(X = 1:obs_data$G,
             FUN = function(g){
               alpha[g]**2 * exp(-(1/(2*(rho[g]**2)))*(D**2))+
                 diag(x = 1e-9, nrow = obs_data$L, ncol = obs_data$L)
             })

    
    KL = 
      lapply(X = 1:obs_data$G,
             FUN = function(g){
               t(chol(Sigma[[g]]))
             })
    zGP = 
      matrix(nrow = obs_data$L, ncol = obs_data$P, byrow = FALSE,
                      data = sapply(X = 1:obs_data$P,
                                    FUN = function(plot){
                                      rnorm(n = obs_data$L, mean = 0, sd = 1)
                                      }))
      
    GP = 
      matrix(nrow = obs_data$L, ncol = obs_data$P, byrow = FALSE,
             data = sapply(X = 1:obs_data$P,
                           FUN = function(plot){
                             KL[[obs_data$group_id[plot]]] %*% zGP[,plot]
                           }))
    
    # plot random effect realization
    plotEffects = rnorm(n = obs_data$P, mean = 0, sd = sigmaPlot)
    
    # linear predictor for each observation
    logmu = 
      sapply(X = 1:obs_data$N,
             FUN = function(i){
               XB[i] +
                 GP[obs_data$location_id[i],obs_data$plot_id[i]] +
                 plotEffects[obs_data$plot_id[i]]
             })
    
    # draw responses
    Y = rnbinom(n = obs_data$N, mu = exp(logmu), size = kappa)
    
    sim_data = obs_data
    sim_data$Y = Y
    
    # save resulting data
    results = 
      list('id' = id,
           'true_params' = list('beta' = beta,
                                'alpha' = alpha,
                                'rho' = rho,
                                'sigmaPlot' = sigmaPlot,
                                'kappa' = kappa),
           'sim_data' = sim_data)
    
    saveRDS(results,
            here::here('02-data',
                       '06-results',
                       'sim_fits',
                       paste0('simdata_', id, '.rds')))
  }


estimate_parameters = 
  function(id){
    library(here)
    library(cmdstanr)
    
    # load the data and the stan model
    stan_model = 
      cmdstanr::cmdstan_model(here::here('03-analysis', 'litterduff_fwd.stan'))
    
    # fit the model
    sim_fit = 
      stan_model$sample(
        data = readRDS(here::here('02-data', '03-results', 'sim_fits', 
                                  paste0('simdata_',id,'.rds')))$sim_data,
        parallel_chains = 1,
        output_dir = here::here('02-data', '03-results', 'sim_fits'),
        output_basename = paste0('simfit_', id),
        thin = 4,
        iter_warmup = 1000,
        iter_sampling = 1000)
    
    sim_fit$save_object(here::here('02-data', '03-results', 'sim_fits',
                                   paste0('simfit_',id,'.rds')))
    
    print(paste0('successful run ', id))
    return(paste0('Successful run ', id))
  }

rank_parameters = 
  function(id){
    
    simulated_data = readRDS(here::here('02-data',
                                        '03-results',
                                        'sim_fits',
                                        paste0('simdata_',id,'.rds')))
    
    sim_fit = readRDS(here::here('02-data',
                                 '03-results',
                                 'sim_fits',
                                 paste0('simfit_',id,'.rds')))
    
    # compare posterior samples against true parameter values
    parameter_ranks = 
      
      # start with the (thinned) posterior samples
      as_draws_df(sim_fit$draws(variables = 
                                  c('beta', 'rho_pre', 'rho_post', 
                                    'alpha_pre', 'alpha_post', 'sigmaPlot', 
                                    'kappa'))) %>%
      
      # create indicator columns comparing posterior samples to true values
      mutate(
        Rbeta01 = as.integer(`beta[1]` < simulated_data$true_params$beta[1]),
        Rbeta02 = as.integer(`beta[2]` < simulated_data$true_params$beta[2]),
        Rbeta03 = as.integer(`beta[3]` < simulated_data$true_params$beta[3]),
        Rbeta04 = as.integer(`beta[4]` < simulated_data$true_params$beta[4]),
        Rbeta05 = as.integer(`beta[5]` < simulated_data$true_params$beta[5]),
        Rbeta06 = as.integer(`beta[6]` < simulated_data$true_params$beta[6]),
        
        Ralphapre = as.integer(alpha_pre < simulated_data$true_params$alpha_pre),
        Ralphapost = as.integer(alpha_post < simulated_data$true_params$alpha_post),
        Rrhopre = as.integer(rho_pre < simulated_data$true_params$rho_pre),
        Rrhopost = as.integer(rho_post < simulated_data$true_params$rho_post),
        
        Rsigmaplot = as.integer(sigmaPlot < simulated_data$true_params$sigmaPlot),
        
        Rkappa = as.integer(kappa < simulated_data$true_params$kappa)) %>%
      
      select(contains(c('Rbeta', 'Ralphapre', 'Ralphapost', 'Rrhopre', 
                        'Rrhopost', 'Rsigmaplot', 'Rkappa'))) %>%
      
      # sum the indicators to get the rank for each variable
      summarise_all(sum)
    
    return(parameter_ranks)
    
  }