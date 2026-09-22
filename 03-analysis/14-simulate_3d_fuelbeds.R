
#### setup #####################################################################


library(here)
library(tidyverse) 
library(posterior)
library(VoxR)

library(data.table)

set.seed(112188)

#### duff and litter ###########################################################




simulate_litterduff = 
  function(params, data){
    
    # extract parameters
    beta = 
      params %>% select(contains('beta')) %>% as.numeric()
    alpha = 
      params %>% select(contains('alpha')) %>% as.numeric() 
    rho = 
      params %>% select(contains('rho')) %>% as.numeric()
    sigmaPlot = 
      params %>% select(contains('sigmaPlot')) %>% as.numeric()
    kappa = 
      params %>% select(contains('kappa')) %>% as.numeric()
    
    # random effects realizations
    zPlot = rnorm(n = data$P, mean = 0, sd = 1)
    effectPlot = sigmaPlot[data$group_id] * zPlot
    
    zGP = matrix(nrow = data$L, ncol = data$P, 
                 data = rnorm(n = data$L * data$P, mean = 0, sd = 1))
    
    dists = dist(data$coords) %>% as.matrix()
      
    # construct GP effects
     K = 
       array(dim = c(data$L, data$L, data$G),
             dimnames = 
               list('l' = 1:data$L, 'lprime' = 1:data$L,'g' = 1:data$G),
             data = 
               sapply(X = 1:data$G,
                      FUN = function(g){
                        sapply(X = 1:data$L,
                               FUN = function(l){
                                 sapply(X = 1:data$L,
                                        FUN = function(lprime){
                                          ((alpha[g]**2)*
                                            exp(-(dists[lprime,l]**2)/
                                                  (2*rho[g]**2)))
                                        })
                               })
                      }))
     
     for (g in 1:data$G){
       K[,,g] = K[,,g]+diag(nrow = nrow(K[,,g]),
                            ncol = ncol(K[,,g]),
                            x = 1e-9)
     }
     
     KL = 
       array(dim = 
               c(data$L,data$L,data$G),
             dimnames = 
               list('l' = 1:data$L,'lprime' = 1:data$L,'g' = 1:data$G),
             data = 
               sapply(X = 1:data$G,
                      FUN = function(g){
                        t(chol(K[,,g]))
                      }))
   
     # realized random effects
     GP = 
       matrix(nrow = data$L,
              ncol = data$P,
              data = 
                sapply(X = 1:data$P,
                       FUN = function(p){
                         KL[,,data$group_id[p]] %*%
                           zGP[,p]
                       }))
  
     logmu = 
       beta[data$group_id[data$plot_id]]+
       effectPlot[data$plot_id]+
       GP[data$location_id, data$plot_id]
     
     Y = sapply(X = 1:data$N,
                FUN = function(i){
                  rnbinom(n = 1, 
                          mu = exp(logmu[i]), 
                          size = kappa[data$group_id[data$plot_id]])
                })
     
     return(Y)
     
     }


duff_sim_data = 
  list(
    G = 1,
    P = 1,
    coords = 
      expand.grid(x = seq(0.25, 29.75, 0.5), y = seq(0.25, 29.75, 0.5)) %>% 
      as.matrix()
  )

duff_sim_data$L = nrow(duff_sim_data$coords)
duff_sim_data$N = duff_sim_data$P * duff_sim_data$L
duff_sim_data$group_id = rep(1:duff_sim_data$G, times = duff_sim_data$P)
duff_sim_data$plot_id = rep(1:duff_sim_data$P, each = duff_sim_data$L)
duff_sim_data$location_id = rep(1:duff_sim_data$L, times = duff_sim_data$P)

duff_params = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'duff_fit.rds'))$draws() %>%
  as_draws_df() %>%
  summarise_all(median)

duff_sim_data$Y = simulate_litterduff(params = duff_params, data = duff_sim_data)

duff_sim_df = 
  data.frame(Y = duff_sim_data$Y,
             group_id = duff_sim_data$group_id[duff_sim_data$plot_id],
             plot_id = duff_sim_data$plot_id,
             location_id = duff_sim_data$location_id,
             x_coord = duff_sim_data$coords[duff_sim_data$location_id,1],
             y_coord = duff_sim_data$coords[duff_sim_data$location_id,2])




ggplot(data = duff_sim_df,
       aes(x = x_coord, y = y_coord, fill = Y))+
  geom_tile()+
  coord_fixed()+
  scale_fill_viridis_c(option = 'B', direction = -1)+
  theme_minimal()


litter_sim_data = 
  list(
    G = 1,
    P = 1,
    coords = 
      expand.grid(x = seq(0.25, 29.75, 0.5), y = seq(0.25, 29.75, 0.5)) %>% 
      as.matrix()
  )

litter_sim_data$L = nrow(litter_sim_data$coords)
litter_sim_data$N = litter_sim_data$P * litter_sim_data$L
litter_sim_data$group_id = rep(1:litter_sim_data$G, times = litter_sim_data$P)
litter_sim_data$plot_id = rep(1:litter_sim_data$P, each = litter_sim_data$L)
litter_sim_data$location_id = rep(1:litter_sim_data$L, times = litter_sim_data$P)

litter_params = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'litter_fit.rds'))$draws() %>%
  as_draws_df() %>%
  summarise_all(median)

litter_sim_data$Y = simulate_litterduff(params = litter_params, data = litter_sim_data)

litter_sim_df = 
  data.frame(Y = litter_sim_data$Y,
             group_id = litter_sim_data$group_id[litter_sim_data$plot_id],
             plot_id = litter_sim_data$plot_id,
             location_id = litter_sim_data$location_id,
             x_coord = litter_sim_data$coords[litter_sim_data$location_id,1],
             y_coord = litter_sim_data$coords[litter_sim_data$location_id,2])

ggplot(data = litter_sim_df,
       aes(x = x_coord, y = y_coord, fill = Y))+
  geom_tile()+
  coord_fixed()+
  scale_fill_viridis_c(option = 'B', direction = -1)+
  theme_minimal()

#### fine woody debris #########################################################


fwdtallies_sim_data = 
  list(
    G = 1,
    P = 1,
    coords = 
      expand.grid(x = seq(0.5, 29.5, 1), y = seq(0.5, 29.5, 1)) %>% 
      as.matrix()
  )

fwdtallies_sim_data$L = nrow(fwdtallies_sim_data$coords)
fwdtallies_sim_data$N = fwdtallies_sim_data$P * fwdtallies_sim_data$L
fwdtallies_sim_data$group_id = rep(1:fwdtallies_sim_data$G, times = fwdtallies_sim_data$P)
fwdtallies_sim_data$plot_id = rep(1:fwdtallies_sim_data$P, each = fwdtallies_sim_data$L)
fwdtallies_sim_data$location_id = rep(1:fwdtallies_sim_data$L, times = fwdtallies_sim_data$P)

fwdtallies_params = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'fwdtallies_fit.rds'))$draws() %>%
  as_draws_df() %>%
  summarise_all(median)

fwddiams_params = 
  readRDS(here::here('02-data',
                    '06-results',
                    'real_fits',
                    'fwddiams_fit.rds'))$draws() %>%
  as_draws_df() %>%
  summarise_all(median)

simulate_fwdtallies = 
  function(params, data){
    
    # extract parameters
    beta = 
      params %>% select(contains('beta')) %>% as.numeric()
    alpha = 
      params %>% select(contains('alpha')) %>% as.numeric() 
    rho = 
      params %>% select(contains('rho')) %>% as.numeric()
    sigmaPlot = 
      params %>% select(contains('sigmaPlot')) %>% as.numeric()
    tau = 
      params %>% select(contains('tau')) %>% as.numeric()
    kappa = 
      params %>% select(contains('kappa')) %>% as.numeric()
    
    # random effects realizations
    zPlot = rnorm(n = data$P, mean = 0, sd = 1)
    effectPlot = sigmaPlot[data$group_id] * zPlot
    
    zGP = matrix(nrow = data$L, ncol = data$P, 
                 data = rnorm(n = data$L * data$P, mean = 0, sd = 1))
    
    dists = dist(data$coords) %>% as.matrix()
      
    # construct GP effects
     K = 
       array(dim = c(data$L, data$L, data$G),
             dimnames = 
               list('l' = 1:data$L, 'lprime' = 1:data$L,'g' = 1:data$G),
             data = 
               sapply(X = 1:data$G,
                      FUN = function(g){
                        sapply(X = 1:data$L,
                               FUN = function(l){
                                 sapply(X = 1:data$L,
                                        FUN = function(lprime){
                                          ((alpha[g]**2)*
                                            exp(-(dists[lprime,l]**2)/
                                                  (2*rho[g]**2)))
                                        })
                               })
                      }))
     
     for (g in 1:data$G){
       K[,,g] = K[,,g]+diag(nrow = nrow(K[,,g]),
                            ncol = ncol(K[,,g]),
                            x = tau[g])
     }
     
     KL = 
       array(dim = 
               c(data$L,data$L,data$G),
             dimnames = 
               list('l' = 1:data$L,'lprime' = 1:data$L,'g' = 1:data$G),
             data = 
               sapply(X = 1:data$G,
                      FUN = function(g){
                        t(chol(K[,,g]))
                      }))
   
     # realized random effects
     GP = 
       matrix(nrow = data$L,
              ncol = data$P,
              data = 
                sapply(X = 1:data$P,
                       FUN = function(p){
                         KL[,,data$group_id[p]] %*%
                           zGP[,p]
                       }))
  
     logmu = 
       beta[data$group_id[data$plot_id]]+
       effectPlot[data$plot_id]+
       GP[data$location_id, data$plot_id]
     
     Y = sapply(X = 1:data$N,
                FUN = function(i){
                  rnbinom(n = 1, 
                          mu = exp(logmu[i]), 
                          size = kappa[data$group_id[data$plot_id]])
                })
     
     return(Y)
     
  }

fwdtallies_sim_data$Y = simulate_fwdtallies(params = fwdtallies_params, data = fwdtallies_sim_data)

fwdtallies_sim_df = 
  data.frame(Y = fwdtallies_sim_data$Y,
             group_id = fwdtallies_sim_data$group_id[fwdtallies_sim_data$plot_id],
             plot_id = fwdtallies_sim_data$plot_id,
             location_id = fwdtallies_sim_data$location_id,
             x_coord = fwdtallies_sim_data$coords[fwdtallies_sim_data$location_id,1],
             y_coord = fwdtallies_sim_data$coords[fwdtallies_sim_data$location_id,2])

ggplot(data = fwdtallies_sim_df,
       aes(x = x_coord, y = y_coord, fill = Y))+
  geom_tile()+
  coord_fixed()+
  scale_fill_viridis_c(option = 'B', direction = -1)+
  theme_minimal()

simulate_fwddiams = 
  function(params, data){
    data_long = 
      data %>%
      rowid_to_column('transect_id') %>%
      filter(Y > 0) %>%
      group_by(group_id, plot_id, location_id, x_coord, y_coord, transect_id) %>%
      mutate(particle_id = list(seq(1, Y, 1))) %>%
      unnest(cols = c(particle_id)) %>%
      ungroup() %>%
      select(-Y)
    
    intercept = params %>% select(contains('intercept')) %>% as.numeric()
    phi = params %>% select(contains('phi')) %>% as.numeric()
    SD_plot = params %>% select(contains('SD_plot')) %>% as.numeric()
    SD_trans = params %>% select(contains('SD_trans')) %>% as.numeric()

    z_plot = rnorm(n = length(unique(data$plot_id)), mean = 0, sd = 1)
    z_trans = rnorm(n = nrow(data), mean = 0, sd = 1)
    
    data_long$effect_plot = SD_plot[data_long$group_id] * z_plot[data_long$plot_id]
    
    data_long$effect_trans = SD_trans[data_long$group_id]*z_trans[data_long$transect_id]
    
    data_long$intercept = intercept[data_long$group_id]
    data_long$phi = phi[data_long$group_id]
    
    data_long$logmu = data_long$intercept + data_long$effect_plot + data_long$effect_trans
    
    data_long$shape = data_long$phi + 2
    
    data_long$scale = exp(data_long$logmu)*(1+data_long$phi)
    
    # again, comparing the pdf functions listed it looks like what stan and 
    # the invgamma paper i read call the "scale" the invgamma package is calling 
    # the "rate", so we want to use "rate" here even thought its "scale" above
    data_long$diam_cm = 
      sapply(X = 1:nrow(data_long),
             FUN = function(i){
               y = 0
               while(y <= 0 | y >= 7.55){
                 y = invgamma::rinvgamma(n = 1, shape = data_long$shape[i], rate = data_long$scale[i])
               }
               return(y)
             })
  
    data_long = 
      data_long %>%
      mutate(timelag_class = ifelse(diam_cm < 0.64,
                                    '1h',
                                    ifelse(diam_cm >= 0.64 & diam_cm < 2.54,
                                           '10h',
                                           '100h')))
  
    return(data_long)
  }

fwddiams_sim_data_long = 
  simulate_fwddiams(params = fwddiams_params, data = fwdtallies_sim_df)

fwddiams_sim_data = 
  fwdtallies_sim_df%>%
    select(group_id, plot_id, location_id, x_coord, y_coord) %>%
    expand(nesting(group_id, plot_id, location_id, x_coord, y_coord),
           timelag_class = c('1h', '10h', '100h')) %>%
    left_join(
      fwddiams_sim_data_long %>%
        group_by(group_id, plot_id, location_id, x_coord, y_coord, timelag_class) %>%
        summarise(count = n()) %>%
        ungroup()
    ) %>%
    mutate(count = ifelse(is.na(count), 0, count))

head(fwddiams_sim_data)

ggplot(fwddiams_sim_data,
       aes(x = x_coord, y = y_coord, fill = count))+
  geom_tile()+
  coord_fixed()+
  theme_minimal()+
  scale_fill_viridis_c(option = 'C', direction = -1)+
  facet_wrap(~timelag_class)


#### cwd #######################################################################


cwdtallies_sim_data = 
  list(
    G = 1,
    P = 1,
    coords = 
      expand.grid(x = seq(0.5, 29.5, 1), y = seq(0.5, 29.5, 1)) %>% 
      as.matrix()
  )

cwdtallies_sim_data$L = nrow(cwdtallies_sim_data$coords)
cwdtallies_sim_data$N = cwdtallies_sim_data$P * cwdtallies_sim_data$L
cwdtallies_sim_data$group_id = rep(1:cwdtallies_sim_data$G, times = cwdtallies_sim_data$P)
cwdtallies_sim_data$plot_id = rep(1:cwdtallies_sim_data$P, each = cwdtallies_sim_data$L)
cwdtallies_sim_data$location_id = rep(1:cwdtallies_sim_data$L, times = cwdtallies_sim_data$P)

cwdtallies_params = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'cwdtallies_fit.rds'))$draws() %>%
  as_draws_df() %>%
  summarise_all(median)

cwddiams_params = 
  readRDS(here::here('02-data',
                    '06-results',
                    'real_fits',
                    'cwddiams_fit.rds'))$draws() %>%
  as_draws_df() %>%
  summarise_all(median)

simulate_cwdtallies = 
  function(params, data){
    
    # extract parameters
    beta = 
      params %>% select(contains('beta')) %>% as.numeric()
    alpha = 
      params %>% select(contains('alpha')) %>% as.numeric() 
    rho = 
      params %>% select(contains('rho')) %>% as.numeric()
    sigmaPlot = 
      params %>% select(contains('sigmaPlot')) %>% as.numeric()
    
    # random effects realizations
    zPlot = rnorm(n = data$P, mean = 0, sd = 1)
    effectPlot = sigmaPlot[data$group_id] * zPlot
    
    zGP = matrix(nrow = data$L, ncol = data$P, 
                 data = rnorm(n = data$L * data$P, mean = 0, sd = 1))
    
    dists = dist(data$coords) %>% as.matrix()
      
    # construct GP effects
     K = 
       array(dim = c(data$L, data$L, data$G),
             dimnames = 
               list('l' = 1:data$L, 'lprime' = 1:data$L,'g' = 1:data$G),
             data = 
               sapply(X = 1:data$G,
                      FUN = function(g){
                        sapply(X = 1:data$L,
                               FUN = function(l){
                                 sapply(X = 1:data$L,
                                        FUN = function(lprime){
                                          ((alpha[g]**2)*
                                            exp(-(dists[lprime,l]**2)/
                                                  (2*rho[g]**2)))
                                        })
                               })
                      }))
     
     for (g in 1:data$G){
       K[,,g] = K[,,g]+diag(nrow = nrow(K[,,g]),
                            ncol = ncol(K[,,g]),
                            x = 1e-9)
     }
     
     KL = 
       array(dim = 
               c(data$L,data$L,data$G),
             dimnames = 
               list('l' = 1:data$L,'lprime' = 1:data$L,'g' = 1:data$G),
             data = 
               sapply(X = 1:data$G,
                      FUN = function(g){
                        t(chol(K[,,g]))
                      }))
   
     # realized random effects
     GP = 
       matrix(nrow = data$L,
              ncol = data$P,
              data = 
                sapply(X = 1:data$P,
                       FUN = function(p){
                         KL[,,data$group_id[p]] %*%
                           zGP[,p]
                       }))
  
     logmu = 
       beta[data$group_id[data$plot_id]]+
       effectPlot[data$plot_id]+
       GP[data$location_id, data$plot_id]
     
     Y = sapply(X = 1:data$N,
                FUN = function(i){
                  rpois(n = 1,lambda= exp(logmu[i]))
                })
     
     return(Y)
     
  }

set.seed(110819)
cwdtallies_sim_data$Y = simulate_cwdtallies(params = cwdtallies_params, data = cwdtallies_sim_data)

cwdtallies_sim_df = 
  data.frame(Y = cwdtallies_sim_data$Y,
             group_id = cwdtallies_sim_data$group_id[cwdtallies_sim_data$plot_id],
             plot_id = cwdtallies_sim_data$plot_id,
             location_id = cwdtallies_sim_data$location_id,
             x_coord = cwdtallies_sim_data$coords[cwdtallies_sim_data$location_id,1],
             y_coord = cwdtallies_sim_data$coords[cwdtallies_sim_data$location_id,2])

ggplot(data = cwdtallies_sim_df,
       aes(x = x_coord, y = y_coord, fill = Y))+
  geom_tile()+
  coord_fixed()+
  scale_fill_viridis_c(option = 'B', direction = -1)+
  theme_minimal()

simulate_cwddiams = 
  function(params, data){
    data_long = 
      data %>%
      rowid_to_column('transect_id') %>%
      filter(Y > 0) %>%
      group_by(group_id, plot_id, location_id, x_coord, y_coord, transect_id) %>%
      mutate(particle_id = list(seq(1, Y, 1))) %>%
      unnest(cols = c(particle_id)) %>%
      ungroup() %>%
      select(-Y)
    
    mu = params %>% select(contains('mu')) %>% as.numeric()
    phi = params %>% select(contains('phi')) %>% as.numeric()
    
    
    data_long$mu = mu[data_long$group_id]
    data_long$phi = phi[data_long$group_id]
    
    data_long$shape = data_long$phi + 2
    
    data_long$scale = data_long$mu*(1+data_long$phi)
    
    # again, comparing the pdf functions listed it looks like what stan and 
    # the invgamma paper i read call the "scale" the invgamma package is calling 
    # the "rate", so we want to use "rate" here even thought its "scale" above
    data_long$diam_cm = 
      sapply(X = 1:nrow(data_long),
             FUN = function(i){
               y = 0
               while(y <= 7.55 | y >= 110){
                 y = invgamma::rinvgamma(n = 1, shape = data_long$shape[i], rate = data_long$scale[i])
               }
               return(y)
             })
  
    return(data_long)
  }



cwddiams_sim_data_long = 
  simulate_cwddiams(params = cwddiams_params, data = cwdtallies_sim_df)

cwddiams_sim_data = 
  cwdtallies_sim_df %>%
    select(group_id, plot_id, location_id, x_coord, y_coord) %>%
    left_join(
      cwddiams_sim_data_long %>%
        group_by(group_id, plot_id, location_id, x_coord, y_coord) %>%
        summarise(ssd_cm2 = sum(diam_cm**2),
                  count = n()) %>%
        ungroup()
    ) %>%
    mutate(ssd_cm2 = ifelse(is.na(ssd_cm2), 0, ssd_cm2),
           count = ifelse(is.na(count), 0, count))

head(cwddiams_sim_data)

ggplot(cwddiams_sim_data,
       aes(x = x_coord, y = y_coord, fill = ssd_cm2))+
  geom_tile()+
  coord_fixed()+
  theme_minimal()+
  scale_fill_viridis_c(option = 'C', direction = -1)

#### understory veg ############################################################



simulate_vegpa = 
  function(params, data){
    
    # extract parameters
    beta = 
      params %>% select(contains('beta')) %>% as.numeric()
    alpha = 
      params %>% select(contains('alpha')) %>% as.numeric() 
    rho = 
      params %>% select(contains('rho')) %>% as.numeric()
    sigmaPlot = 
      params %>% select(contains('sigmaPlot')) %>% as.numeric()
    
    # random effects realizations
    zPlot = rnorm(n = data$P, mean = 0, sd = 1)
    effectPlot = sigmaPlot[data$group_id] * zPlot
    
    zGP = matrix(nrow = data$L, ncol = data$P, 
                 data = rnorm(n = data$L * data$P, mean = 0, sd = 1))
    
    dists = dist(data$coords) %>% as.matrix()
    
    # construct GP effects
    K = 
      array(dim = c(data$L, data$L, data$G),
            dimnames = 
              list('l' = 1:data$L, 'lprime' = 1:data$L,'g' = 1:data$G),
            data = 
              sapply(X = 1:data$G,
                     FUN = function(g){
                       sapply(X = 1:data$L,
                              FUN = function(l){
                                sapply(X = 1:data$L,
                                       FUN = function(lprime){
                                         ((alpha[g]**2)*
                                            exp(-(dists[lprime,l]**2)/
                                                  (2*rho[g]**2)))
                                       })
                              })
                     }))
    
    for (g in 1:data$G){
      K[,,g] = K[,,g]+diag(nrow = nrow(K[,,g]),
                           ncol = ncol(K[,,g]),
                           x = 1e-9)
    }
    
    KL = 
      array(dim = 
              c(data$L,data$L,data$G),
            dimnames = 
              list('l' = 1:data$L,'lprime' = 1:data$L,'g' = 1:data$G),
            data = 
              sapply(X = 1:data$G,
                     FUN = function(g){
                       t(chol(K[,,g]))
                     }))
    
    # realized random effects
    GP = 
      matrix(nrow = data$L,
             ncol = data$P,
             data = 
               sapply(X = 1:data$P,
                      FUN = function(p){
                        KL[,,data$group_id[p]] %*%
                          zGP[,p]
                      }))
    
    logmu = 
      beta[data$group_id[data$plot_id]]+
      effectPlot[data$plot_id]+
      GP[data$location_id, data$plot_id]
    
    Y = sapply(X = 1:data$N,
               FUN = function(i){
                 rbinom(n = 1,
                        size = 1,
                        prob = boot::inv.logit(logmu[i]))
               })
    
    return(Y)
    
  }


vegpa_sim_data = 
  list(
    G = 1,
    P = 1,
    coords = 
      expand.grid(x = seq(0.25, 29.75, 0.5), y = seq(0.25, 29.75, 0.5)) %>% 
      as.matrix()
  )

vegpa_sim_data$L = nrow(vegpa_sim_data$coords)
vegpa_sim_data$N = vegpa_sim_data$P * vegpa_sim_data$L
vegpa_sim_data$group_id = rep(1:vegpa_sim_data$G, times = vegpa_sim_data$P)
vegpa_sim_data$plot_id = rep(1:vegpa_sim_data$P, each = vegpa_sim_data$L)
vegpa_sim_data$location_id = rep(1:vegpa_sim_data$L, times = vegpa_sim_data$P)

vegpa_params = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'vegpa_fit.rds'))$draws() %>%
  as_draws_df() %>%
  summarise_all(median)

set.seed(110819)
vegpa_sim_data$Y = simulate_vegpa(params = vegpa_params, data = vegpa_sim_data)

vegpa_sim_df = 
  data.frame(Y = vegpa_sim_data$Y,
             group_id = vegpa_sim_data$group_id[vegpa_sim_data$plot_id],
             plot_id = vegpa_sim_data$plot_id,
             location_id = vegpa_sim_data$location_id,
             x_coord = vegpa_sim_data$coords[vegpa_sim_data$location_id,1],
             y_coord = vegpa_sim_data$coords[vegpa_sim_data$location_id,2])




ggplot(data = vegpa_sim_df,
       aes(x = x_coord, y = y_coord, fill = Y))+
  geom_tile()+
  coord_fixed()+
  scale_fill_viridis_c(option = 'B', direction = -1)+
  theme_minimal()

# assign a species and height to each pixel by drawing from the observed data
veg_observed = 
  readRDS(here::here('02-data',
                     '03-clean',
                     'seki_veg.rds')) %>%
  left_join(
    readRDS(here::here('02-data', '02-intermediate', 'plot_locations_true.rds')) %>%
      select(mort, plot_id)
  ) %>%
  mutate(group_id = ifelse(mort=='low',1,
                           ifelse(mort=='mid',2, 3))) %>%
  mutate(transect_length_m = abs(end_m - start_m)) %>%
  select(group_id, transect_length_m, spp, status, height_m) %>%
  group_by(group_id, spp, status, height_m) %>%
  summarise(transect_length_m = sum(transect_length_m)) %>%
  ungroup() %>%
  rowid_to_column('veg_id')

head(vegpa_sim_df)
head(veg_observed)  

vegpa_sim_df$veg_id = 
  sapply(X = 1:nrow(vegpa_sim_df),
       FUN = function(i){
         sample(x = veg_observed$veg_id[veg_observed$group_id==vegpa_sim_df$group_id[i]],
                size = 1,
                replace = TRUE)})

vegpa_sim_df = 
  vegpa_sim_df %>%
  mutate(veg_id = ifelse(Y==1, veg_id, NA)) %>%
  left_join(
    veg_observed %>%
      select(veg_id, spp, status, height_m))



ggplot(data = vegpa_sim_df,
       aes(x = x_coord, y = y_coord, fill = Y))+
  geom_tile()+
  coord_fixed()+
  scale_fill_viridis_c(option = 'B', direction = -1)+
  theme_minimal()

ggplot(data = vegpa_sim_df,
       aes(x = x_coord, y = y_coord, fill = height_m))+
  geom_tile()+
  coord_fixed()+
  scale_fill_viridis_c(option = 'B', direction = -1)+
  theme_minimal()

ggplot(data = vegpa_sim_df,
       aes(x = x_coord, y = y_coord, fill = spp))+
  geom_tile()+
  coord_fixed()+
  #scale_fill_viridis_d(option = 'B', direction = -1)+
  theme_minimal()



#### trees #####################################################################


trees_sim_data = 
  list(
    G = 1,
    P = 1,
    coords = 
      expand.grid(x = seq(0.5, 29.5, 1), y = seq(0.5, 29.5, 1)) %>% 
      as.matrix()
  )

trees_sim_data$L = nrow(trees_sim_data$coords)
trees_sim_data$N = trees_sim_data$P * trees_sim_data$L
trees_sim_data$group_id = rep(1:trees_sim_data$G, times = trees_sim_data$P)
trees_sim_data$plot_id = rep(1:trees_sim_data$P, each = trees_sim_data$L)
trees_sim_data$location_id = rep(1:trees_sim_data$L, times = trees_sim_data$P)

trees_params = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'trees_fit.rds'))$draws() %>%
  as_draws_df() %>%
  summarise_all(median)


simulate_trees = 
  function(params, data){
    
    # extract parameters
    beta = 
      params %>% select(contains('beta')) %>% as.numeric()
    alpha = 
      params %>% select(contains('alpha')) %>% as.numeric() 
    rho = 
      params %>% select(contains('rho')) %>% as.numeric()
    sigmaPlot = 
      params %>% select(contains('sigmaPlot')) %>% as.numeric()
    
    # random effects realizations
    zPlot = rnorm(n = data$P, mean = 0, sd = 1)
    effectPlot = sigmaPlot[data$group_id] * zPlot
    
    zGP = matrix(nrow = data$L, ncol = data$P, 
                 data = rnorm(n = data$L * data$P, mean = 0, sd = 1))
    
    dists = dist(data$coords) %>% as.matrix()
      
    # construct GP effects
     K = 
       array(dim = c(data$L, data$L, data$G),
             dimnames = 
               list('l' = 1:data$L, 'lprime' = 1:data$L,'g' = 1:data$G),
             data = 
               sapply(X = 1:data$G,
                      FUN = function(g){
                        sapply(X = 1:data$L,
                               FUN = function(l){
                                 sapply(X = 1:data$L,
                                        FUN = function(lprime){
                                          ((alpha[g]**2)*
                                            exp(-(dists[lprime,l]**2)/
                                                  (2*rho[g]**2)))
                                        })
                               })
                      }))
     
     for (g in 1:data$G){
       K[,,g] = K[,,g]+diag(nrow = nrow(K[,,g]),
                            ncol = ncol(K[,,g]),
                            x = 1e-9)
     }
     
     KL = 
       array(dim = 
               c(data$L,data$L,data$G),
             dimnames = 
               list('l' = 1:data$L,'lprime' = 1:data$L,'g' = 1:data$G),
             data = 
               sapply(X = 1:data$G,
                      FUN = function(g){
                        t(chol(K[,,g]))
                      }))
   
     # realized random effects
     GP = 
       matrix(nrow = data$L,
              ncol = data$P,
              data = 
                sapply(X = 1:data$P,
                       FUN = function(p){
                         KL[,,data$group_id[p]] %*%
                           zGP[,p]
                       }))
  
     logmu = 
       beta[data$group_id[data$plot_id]]+
       effectPlot[data$plot_id]+
       GP[data$location_id, data$plot_id]
     
     Y = sapply(X = 1:data$N,
                FUN = function(i){
                  rpois(n = 1,lambda= exp(logmu[i]))
                })
     
     return(Y)
     
  }

set.seed(101722)
trees_sim_data$Y = simulate_trees(params = trees_params, data = trees_sim_data)

trees_sim_df = 
  data.frame(Y = trees_sim_data$Y,
             group_id = trees_sim_data$group_id[trees_sim_data$plot_id],
             plot_id = trees_sim_data$plot_id,
             location_id = trees_sim_data$location_id,
             x_coord = trees_sim_data$coords[trees_sim_data$location_id,1],
             y_coord = trees_sim_data$coords[trees_sim_data$location_id,2])

ggplot(data = trees_sim_df,
       aes(x = x_coord, y = y_coord, fill = Y))+
  geom_tile()+
  coord_fixed()+
  scale_fill_viridis_c(option = 'B', direction = -1)+
  theme_minimal()

trees_observed = 
  readRDS(here::here('02-data', '03-clean', 'seki_trees.rds')) %>%
      filter(dbh_cm >= 11.4) %>%
  rowid_to_column('tree_id.i') %>%
  left_join(
    readRDS(here::here('02-data', '02-intermediate', 'plot_locations_true.rds')) %>%
    mutate(group_id = ifelse(mort == 'low', 
                             1,
                             ifelse(mort == 'mid',
                                    2,
                                    3))) %>%
    select(plot_id, group_id)
  )

head(trees_observed)
trees_sim_treelist = 
  
  # convert the continuous grid into a treelist
  trees_sim_df %>%
  filter(Y>0) %>%
  group_by(group_id, plot_id, location_id, x_coord, y_coord) %>%
  mutate(tree_i = list(seq(1, Y, 1))) %>%
  unnest(cols = c(tree_i)) %>%
  ungroup() %>%
  
  # fuzz the coordinates
  mutate(x_coord = x_coord + runif(n = n(), min = -0.5, max = 0.5),
         y_coord = y_coord + runif(n = n(), min = -0.5, max = 0.5)) %>%
  select(-Y, -tree_i)

trees_sim_treelist$tree_id.i = 
  sapply(X = 1:nrow(trees_sim_treelist),
         FUN = function(i){
           sample(x = 
                    trees_observed$tree_id.i[
                      trees_observed$group_id==trees_sim_treelist$group_id[i]],
                  size = 1,
                  replace = TRUE)
         })


trees_sim_treelist = 
  trees_sim_treelist %>%
  left_join(
    trees_observed %>%
      select(tree_id.i, status, spp, dbh_cm, height_m, htcb_m, decay_class),
    by = c('tree_id.i' = 'tree_id.i') 
  )

trees_sim_treelist

ggplot(trees_sim_treelist,
       aes(x = x_coord, y = y_coord, size = dbh_cm, color = spp))+
  geom_point()+
  theme_minimal()+
  coord_fixed()



#### saplings ##################################################################

saplings_sim_data = 
  list(
    G = 1,
    P = 1,
    coords = 
      expand.grid(x = seq(0.5, 29.5, 1), y = seq(0.5, 29.5, 1)) %>% 
      as.matrix()
  )

saplings_sim_data$L = nrow(saplings_sim_data$coords)
saplings_sim_data$N = saplings_sim_data$P * saplings_sim_data$L
saplings_sim_data$group_id = rep(1:saplings_sim_data$G, times = saplings_sim_data$P)
saplings_sim_data$plot_id = rep(1:saplings_sim_data$P, each = saplings_sim_data$L)
saplings_sim_data$location_id = rep(1:saplings_sim_data$L, times = saplings_sim_data$P)

saplings_params = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'saplings_fit.rds'))$draws() %>%
  as_draws_df() %>%
  summarise_all(median)


simulate_saplings = 
  function(params, data){
    
    # extract parameters
    beta = 
      params %>% select(contains('beta')) %>% as.numeric()
    alpha = 
      params %>% select(contains('alpha')) %>% as.numeric() 
    rho = 
      params %>% select(contains('rho')) %>% as.numeric()
    sigmaPlot = 
      params %>% select(contains('sigmaPlot')) %>% as.numeric()
    
    # random effects realizations
    zPlot = rnorm(n = data$P, mean = 0, sd = 1)
    effectPlot = sigmaPlot[data$group_id] * zPlot
    
    zGP = matrix(nrow = data$L, ncol = data$P, 
                 data = rnorm(n = data$L * data$P, mean = 0, sd = 1))
    
    dists = dist(data$coords) %>% as.matrix()
    
    # construct GP effects
    K = 
      array(dim = c(data$L, data$L, data$G),
            dimnames = 
              list('l' = 1:data$L, 'lprime' = 1:data$L,'g' = 1:data$G),
            data = 
              sapply(X = 1:data$G,
                     FUN = function(g){
                       sapply(X = 1:data$L,
                              FUN = function(l){
                                sapply(X = 1:data$L,
                                       FUN = function(lprime){
                                         ((alpha[g]**2)*
                                            exp(-(dists[lprime,l]**2)/
                                                  (2*rho[g]**2)))
                                       })
                              })
                     }))
    
    for (g in 1:data$G){
      K[,,g] = K[,,g]+diag(nrow = nrow(K[,,g]),
                           ncol = ncol(K[,,g]),
                           x = 1e-9)
    }
    
    KL = 
      array(dim = 
              c(data$L,data$L,data$G),
            dimnames = 
              list('l' = 1:data$L,'lprime' = 1:data$L,'g' = 1:data$G),
            data = 
              sapply(X = 1:data$G,
                     FUN = function(g){
                       t(chol(K[,,g]))
                     }))
    
    # realized random effects
    GP = 
      matrix(nrow = data$L,
             ncol = data$P,
             data = 
               sapply(X = 1:data$P,
                      FUN = function(p){
                        KL[,,data$group_id[p]] %*%
                          zGP[,p]
                      }))
    
    logmu = 
      beta[data$group_id[data$plot_id]]+
      effectPlot[data$plot_id]+
      GP[data$location_id, data$plot_id]
    
    Y = sapply(X = 1:data$N,
               FUN = function(i){
                 rpois(n = 1,lambda= exp(logmu[i]))
               })
    
    return(Y)
    
  }

set.seed(110819)
saplings_sim_data$Y = simulate_saplings(params = saplings_params, data = saplings_sim_data)

saplings_sim_df = 
  data.frame(Y = saplings_sim_data$Y,
             group_id = saplings_sim_data$group_id[saplings_sim_data$plot_id],
             plot_id = saplings_sim_data$plot_id,
             location_id = saplings_sim_data$location_id,
             x_coord = saplings_sim_data$coords[saplings_sim_data$location_id,1],
             y_coord = saplings_sim_data$coords[saplings_sim_data$location_id,2])

ggplot(data = saplings_sim_df,
       aes(x = x_coord, y = y_coord, fill = Y))+
  geom_tile()+
  coord_fixed()+
  scale_fill_viridis_c(option = 'B', direction = -1)+
  theme_minimal()

saplings_observed = 
  readRDS(here::here('02-data', '03-clean', 'seki_trees.rds')) %>%
      filter(dbh_cm < 11.4) %>%
  rowid_to_column('tree_id.i') %>%
  left_join(
    readRDS(here::here('02-data', '02-intermediate', 'plot_locations_true.rds')) %>%
    mutate(group_id = ifelse(mort == 'low', 
                             1,
                             ifelse(mort == 'mid',
                                    2,
                                    3))) %>%
    select(plot_id, group_id)
  )

head(saplings_observed)
saplings_sim_treelist = 
  
  # convert the continuous grid into a treelist
  saplings_sim_df %>%
  filter(Y>0) %>%
  group_by(group_id, plot_id, location_id, x_coord, y_coord) %>%
  mutate(tree_i = list(seq(1, Y, 1))) %>%
  unnest(cols = c(tree_i)) %>%
  ungroup() %>%
  
  # fuzz the coordinates
  mutate(x_coord = x_coord + runif(n = n(), min = -0.5, max = 0.5),
         y_coord = y_coord + runif(n = n(), min = -0.5, max = 0.5)) %>%
  select(-Y, -tree_i)

saplings_sim_treelist$tree_id.i = 
  sapply(X = 1:nrow(saplings_sim_treelist),
         FUN = function(i){
           sample(x = 
                    saplings_observed$tree_id.i[
                      saplings_observed$group_id==saplings_sim_treelist$group_id[i]],
                  size = 1,
                  replace = TRUE)
         })


saplings_sim_treelist = 
  saplings_sim_treelist %>%
  left_join(
    saplings_observed %>%
      select(tree_id.i, status, spp, dbh_cm, height_m, htcb_m, decay_class),
    by = c('tree_id.i' = 'tree_id.i') 
  )

saplings_sim_treelist

ggplot(saplings_sim_treelist,
       aes(x = x_coord, y = y_coord, size = dbh_cm, color = spp))+
  geom_point()+
  theme_minimal()+
  coord_fixed()




#### Convert surface fuels to biomass ##########################################

# load constants tables
species_codes = readRDS(here::here('02-data',
                                   '02-intermediate',
                                   'van_wagtendonk',
                                   'species_codes.rds'))

kvals = readRDS(here::here('02-data',
                           '02-intermediate',
                           'van_wagtendonk',
                           'kvals.rds'))


litterduff_coeffs = readRDS(here::here('02-data',
                                       '02-intermediate',
                                       'van_wagtendonk',
                                       'litterduff_coeffs.rds'))


QMDcm = readRDS(here::here('02-data',
                           '02-intermediate',
                           'van_wagtendonk',
                           'QMDcm.rds'))

SEC = readRDS(here::here('02-data',
                         '02-intermediate',
                         'van_wagtendonk',
                         'SEC.rds'))

SG = readRDS(here::here('02-data',
                        '02-intermediate',
                        'van_wagtendonk',
                        'SG.rds'))

sg_1000r = readRDS(here::here('02-data',
                              '02-intermediate',
                              'van_wagtendonk',
                              'vw96_sg_1000r.rds'))

overstory_composition = 
  trees_observed %>%
  
  # keep only big trees (>11.4cm DBH)
  filter(dbh_cm >= 11.4) %>%
  
  # map species wihout coefficients from van wagtendonk's work to 'OTHER' 
  # category, will get the 'All species' coefficient
  mutate(spp = ifelse(!is.element(spp, as.character(species_codes$species_code)),
                      'OTHER',
                      spp)) %>%
  
  mutate(ba_m2 = pi*((dbh_cm/100)/2), # basal area of each tree in m2
         # scaled by 0.05 ha plot size (10x30m + 10x10m + 10x10m)
         ba_m2ha = ba_m2*(1/0.05)) %>%
  
  group_by(group_id, plot_id, spp) %>%
  summarise(ba_m2ha = sum(ba_m2ha)) %>%
  ungroup() %>%
  
  # fill in 0s for species which werent present on a plot
  complete(nesting(group_id, plot_id), spp) %>%
  mutate(ba_m2ha = ifelse(is.na(ba_m2ha), 0, ba_m2ha)) %>%
  
  # get the total BA/ha on each plot
  left_join(x = .,
            y = 
              group_by(., group_id, plot_id) %>% 
              summarise(total_ba_m2ha = sum(ba_m2ha)) %>%
              ungroup()) %>%
  
  # get the proportion of plot total basal area occoupied by each species
  mutate(pBA = ba_m2ha / total_ba_m2ha) %>%
  
  # aggregate to mortality-class average composition
  group_by(group_id, spp) %>%
  summarise(pBA = mean(pBA)) %>%
  ungroup()

overstory_composition

# get BA-weighted average coefficient for fuel load (kg/m2) as a function of 
# litterduff_depth (cm)
litterduff_coeffs = 
  overstory_composition %>%
  left_join(litterduff_coeffs) %>%
  mutate(weighted_litterduff = pBA * litterduff_coeff,
         weighted_litter = pBA * litter_coeff,
         weighted_duff = pBA * duff_coeff) %>%
  group_by(group_id) %>%
  summarise(weighted_coeff_litterduff = sum(weighted_litterduff),
            weighted_coeff_litter = sum(weighted_litter),
            weighted_coeff_duff = sum(weighted_duff)) %>%
  ungroup()

litterduff_coeffs

# get BA-weighted SEC (secant of acute angle) for each timelag class
SEC = 
  overstory_composition %>%
  left_join(SEC) %>%
  select(group_id, spp, pBA, x1h, x10h, x100h, x1000s = x1000h, x1000r = x1000h) %>%
  pivot_longer(c(x1h, x10h, x100h, x1000s, x1000r),
               names_to = 'timelag_class',
               values_to = 'sec',
               names_prefix = 'x') %>%
  mutate(weighted = sec*pBA) %>%
  group_by(group_id, timelag_class) %>%
  summarise(weighted_sec = sum(weighted)) %>%
  ungroup()


# get the BA-weighted average specific gravity for each timelag class
SG = 
  overstory_composition %>%
  left_join(SG) %>%
  select(group_id, spp, pBA, x1h, x10h, x100h, x1000s) %>%
  pivot_longer(c(x1h, x10h, x100h, x1000s),
               names_to = 'timelag_class',
               values_to = 'sg',
               names_prefix = 'x') %>%
  mutate(weighted = sg*pBA) %>%
  group_by(group_id, timelag_class) %>%
  summarise(weighted_sg = sum(weighted)) %>%
  ungroup() %>%
  
  # and add an entry for 1000h rotten, which doesn't vary by species
  rbind(.,
        data.frame(group_id = c(1, 2, 3),
                   timelag_class = c('1000r', '1000r', '1000r'),
                   weighted_sg = rep(sg_1000r, times = 3)) %>%
          mutate(timelag_class = as.character(timelag_class)))

# get BA-weighted average QMD (cm2) for each FWD timelag class
QMDcm = 
  overstory_composition %>%
  left_join(QMDcm) %>%
  select(group_id, spp, pBA, x1h, x10h, x100h) %>%
  pivot_longer(c(x1h, x10h, x100h),
               names_to = 'timelag_class',
               values_to = 'qmd_cm',
               names_prefix = 'x') %>%
  mutate(weighted = qmd_cm*pBA) %>%
  group_by(group_id, timelag_class) %>%
  summarise(weighted_qmd = sum(weighted)) %>%
  ungroup()

head(litter_sim_df)

litter_sim_df = 
  litter_sim_df %>%
  rename(litter_cm = Y) %>%
  left_join(litterduff_coeffs %>%
              select(group_id, weighted_coeff_litter)) %>%
  mutate(litter_kgm2 = litter_cm*weighted_coeff_litter,
         litter_mgha = litter_kgm2 * 10) %>%
  select(group_id, plot_id, x_coord, y_coord, litter_cm, litter_kgm2)

duff_sim_df = 
  duff_sim_df %>%
  rename(duff_cm = Y) %>%
  left_join(litterduff_coeffs %>%
              select(group_id, weighted_coeff_duff)) %>%
  mutate(duff_kgm2 = duff_cm*weighted_coeff_duff,
         duff_mgha = duff_kgm2 * 10) %>%
  select(group_id, plot_id, x_coord, y_coord, duff_cm, duff_kgm2)


ggplot(litter_sim_df,
       aes(x = x_coord, y = y_coord, fill = litter_kgm2))+
  geom_tile()+
  coord_fixed()+
  theme_minimal()+
  scale_fill_viridis_c(option = 'D', direction = -1)

ggplot(duff_sim_df,
       aes(x = x_coord, y = y_coord, fill = duff_kgm2))+
  geom_tile()+
  coord_fixed()+
  theme_minimal()+
  scale_fill_viridis_c(option = 'D', direction = -1)

head(fwddiams_sim_data)

fwd_sim_df = 
  fwddiams_sim_data %>%
  left_join(QMDcm) %>%
  left_join(SEC) %>%
  left_join(SG) %>%
  mutate(transect_length_m = 1,
         slp_c = 1,
         k = 1.234) %>%
  mutate(fwd_mgha = (k*weighted_qmd*weighted_sg*weighted_sec*slp_c*count)/transect_length_m,
         fwd_kgm2 = fwd_mgha*0.1) %>%
  select(group_id, plot_id, x_coord, y_coord, timelag_class, count, fwd_kgm2)

head(fwd_sim_df)

ggplot(fwd_sim_df,
       aes(x = x_coord, y = y_coord, fill = fwd_kgm2))+
  geom_tile()+
  theme_minimal()+
  coord_fixed()+
  scale_fill_viridis_c(direction = -1, option = 'D')+
  facet_wrap(~timelag_class)

cwd_sim_df = 
  cwddiams_sim_data %>%
  mutate(timelag_class = '1000s') %>%
  left_join(SEC) %>%
  left_join(SG) %>%
  mutate(transect_length_m = 1, slp_c = 1, k = 1.234) %>%
  mutate(cwd_mgha = (k*ssd_cm2*weighted_sg*weighted_sec*slp_c)/transect_length_m,
         cwd_kgm2 = cwd_mgha*0.1) %>%
  select(group_id, plot_id, x_coord, y_coord,count, ssd_cm2, cwd_kgm2)

head(cwd_sim_df)

ggplot(cwd_sim_df,
       aes(x = x_coord, y = y_coord, fill = cwd_kgm2))+
  geom_tile()+
  theme_minimal()+
  coord_fixed()+
  scale_fill_viridis_c(option = 'D', direction = -1)

cwddiams_sim_data %>%
  arrange(desc(ssd_cm2)) %>%
  head()

# volume of a 110cm diameter by 100cmm length log
testvol_cm3 = ((110/2)**2)*pi*100

# douglas fir is .512 g/cm3, converted from grams to kg
testvol_cm3*0.512*(1/1000) # gives ~500 kg / m2


# case where completely overlapping
# treat given slices as a circle (or two, one for top and one for bottom)
# need area of a chord on a circle
# draw different ways for arc of circle to traverse a square, work out area of top and of bottom
# average them
# given the center of the cone and the radius of the cross section, work out a formulat 
# for overlap of each voxel with the circular slice of the crown

head(vegpa_sim_df)

#### biomass of trees and saplings #############################################

head(trees_sim_treelist)

table9 = read.csv(here::here('02-data',
                             '00-source',
                             'gill',
                             'table9.csv'))

lockhart_table_3 = 
  read.csv(here::here('02-data',
                      '00-source',
                      'lockhart',
                      'lockhart_table_3.csv'))



ref_species = 
  read.csv(here::here('02-data',
                      '00-source',
                      'fia',
                      'REF_SPECIES.csv')) %>%
  as_tibble() %>%
  
  select(GENUS, SPECIES,
         SPECIES_SYMBOL, 
         JENKINS_TOTAL_B1, JENKINS_TOTAL_B2,
         JENKINS_STEM_WOOD_RATIO_B1, JENKINS_STEM_WOOD_RATIO_B2,
         JENKINS_STEM_BARK_RATIO_B1, JENKINS_STEM_BARK_RATIO_B2,
         JENKINS_FOLIAGE_RATIO_B1, JENKINS_FOLIAGE_RATIO_B2,
         JENKINS_SAPLING_ADJUSTMENT, 
         MC_PCT_GREEN_WOOD,
         MC_PCT_GREEN_BARK,
         STANDING_DEAD_DECAY_RATIO1,
         STANDING_DEAD_DECAY_RATIO2,
         STANDING_DEAD_DECAY_RATIO3,
         STANDING_DEAD_DECAY_RATIO4,
         STANDING_DEAD_DECAY_RATIO5)

combo_sim_treelist = 
  trees_sim_treelist %>%
  bind_rows(saplings_sim_treelist) %>%
  # estimate the crown radius of each tree using the results of gill et al. 2000
  # first join in the regression coefficients
  left_join(table9, by = c('spp' = 'Species')) %>%
  
  # fill in unmatched values with the average of the hardwood species from the 
  # lockhart paper; which are different hardwood species but probably better than 
  # a generic conifer coefficient
  mutate(b0 = 
           ifelse(is.na(b0), 
                  lockhart_table_3[lockhart_table_3$species=='average','a'],
                  b0),
         b1 = 
           ifelse(is.na(b1), 
                  lockhart_table_3[lockhart_table_3$species=='average','b'], 
                  b1)) %>%
  
  # now, use the DBH and regression coefficients to estimate the crown radius
  mutate(crownrad_m = b0+(b1*dbh_cm)) %>%
  select(-b0, -b1) %>%
  
  # map the species codes we used to the official species codes used by the FIA
  # table; they don't appear to have generic coefficients to use for unknown 
  # species, so I'm just going to use the ABCO coefficients when species 
  # is unknown, based on the overstory composition data from seki
  left_join(data.frame(spp = c('ABCO', 'ACMA', 'ALRH', 'CADE',
                               'CONU', 'PILA', 'PIPO', 'QUCH', 'QUKE', 'TOCA',
                               'UNK'),
                       spp4 = c('ABCO', 'ACMA3', 'ALRH2', 'CADE27', 'CONU4',
                                'PILA', 'PIPO', 'QUCH2', 'QUKE', 'TOCA',
                                'ABCO')) %>%
              mutate(spp = as.character(spp),
                     spp4 = as.character(spp4))) %>%
  
  # join in the biomass coefficients
  left_join(ref_species %>%
              select(-GENUS, -SPECIES),
            by = c('spp4' = 'SPECIES_SYMBOL')) %>%
  
  # estimate the live biomass using equations from
  # PNW. (2015). A Data Dictionary and User Guide for the PNW-FIADB database. 
  # Pacific Northwest Research Station. 
  # https://www.fs.fed.us/pnw/rma/fia-topics/documentation/documents/PNW_FIADB_P2_Manual_2014.pdf
  # appendix J
  mutate(
    DIA = dbh_cm / 2.54,
    agbiomass_lb = 
      exp(JENKINS_TOTAL_B1+JENKINS_TOTAL_B2*log(DIA*2.54))*2.206,
    stem_ratio = 
      exp(JENKINS_STEM_WOOD_RATIO_B1+JENKINS_STEM_WOOD_RATIO_B2/(DIA*2.54)),
    bark_ratio = 
      exp(JENKINS_STEM_BARK_RATIO_B1+JENKINS_STEM_BARK_RATIO_B2/(DIA*2.54)),
    foliage_ratio = 
      exp(JENKINS_FOLIAGE_RATIO_B1+JENKINS_FOLIAGE_RATIO_B2/(DIA*2.54))) %>%
  
  mutate(
    # the foliage ratio equation breaks down for very small trees, giving 
    # ratios well over 1; cap it such that bark+stem+foliage=1
    foliage_ratio = 
      ifelse(bark_ratio+foliage_ratio+stem_ratio>1,
             1-bark_ratio-stem_ratio,
             foliage_ratio),
    stembiomass_lb = 
      agbiomass_lb*stem_ratio,
    barkbiomass_lb = 
      agbiomass_lb*bark_ratio,
    foliagebiomass_lb = 
      agbiomass_lb*foliage_ratio,
    branchbiomass_lb = 
      agbiomass_lb-stembiomass_lb-barkbiomass_lb-foliagebiomass_lb,
    
  ) %>%

  # adjust for snags
  mutate(
    
    # stem is woody, gets woody decay
    stembiomass_lb = 
      ifelse(status=='D',
             ifelse(decay_class=='1',
                    stembiomass_lb*STANDING_DEAD_DECAY_RATIO1,
                    ifelse(decay_class=='2',
                           stembiomass_lb*STANDING_DEAD_DECAY_RATIO2,
                           ifelse(decay_class=='3',
                                  stembiomass_lb*STANDING_DEAD_DECAY_RATIO3,
                                  ifelse(decay_class=='4',
                                         stembiomass_lb*STANDING_DEAD_DECAY_RATIO4,
                                         ifelse(decay_class=='5',
                                                stembiomass_lb*STANDING_DEAD_DECAY_RATIO5,
                                                stembiomass_lb))))),
             stembiomass_lb),
    
    # assume bark decays at the same rate as the wood (not great but also
    # not that important, and snag decay class provides little info about
    # how much bark is remaining)
    barkbiomass_lb = 
      ifelse(status=='D',
             ifelse(decay_class=='1',
                    barkbiomass_lb*STANDING_DEAD_DECAY_RATIO1,
                    ifelse(decay_class=='2',
                           barkbiomass_lb*STANDING_DEAD_DECAY_RATIO2,
                           ifelse(decay_class=='3',
                                  barkbiomass_lb*STANDING_DEAD_DECAY_RATIO3,
                                  ifelse(decay_class=='4',
                                         barkbiomass_lb*STANDING_DEAD_DECAY_RATIO4,
                                         ifelse(decay_class=='5',
                                                barkbiomass_lb*STANDING_DEAD_DECAY_RATIO5,
                                                barkbiomass_lb))))),
             barkbiomass_lb),
    
    # branches are woody, get woody decay
    branchbiomass_lb = 
      ifelse(status=='D',
             ifelse(decay_class=='1',
                    branchbiomass_lb*STANDING_DEAD_DECAY_RATIO1,
                    ifelse(decay_class=='2',
                           branchbiomass_lb*STANDING_DEAD_DECAY_RATIO2,
                           ifelse(decay_class=='3',
                                  branchbiomass_lb*STANDING_DEAD_DECAY_RATIO3,
                                  ifelse(decay_class=='4',
                                         branchbiomass_lb*STANDING_DEAD_DECAY_RATIO4,
                                         ifelse(decay_class=='5',
                                                branchbiomass_lb*STANDING_DEAD_DECAY_RATIO5,
                                                branchbiomass_lb))))),
             branchbiomass_lb),
    
    # assume all snags have no foliage
    foliagebiomass_lb = 
      ifelse(status=='D',
             0,
             foliagebiomass_lb)
  ) %>%
  
  # recalculate biomass aggregatsions after doing the decay adjustments
  mutate(
    agbiomass_lb = 
      stembiomass_lb + barkbiomass_lb + branchbiomass_lb + foliagebiomass_lb,
    
    # lump the bark in with the stem
    bolebiomass_lb = stembiomass_lb+barkbiomass_lb,
    
    crownbiomass_lb = branchbiomass_lb+foliagebiomass_lb) %>%
  
  # convert to metric for the values we intend to keep
  mutate(
    agbiomass_kg = agbiomass_lb *0.453592,
    bolebiomass_kg = bolebiomass_lb* 0.453592,
    crownbiomass_kg = crownbiomass_lb * 0.453592,
    foliagebiomass_kg = foliagebiomass_lb * 0.453592
  ) %>%
  
  # ditch all the intermediate columns
  select(
    group_id, plot_id,  x_coord, y_coord,
    status, spp, dbh_cm, height_m, htcb_m, decay_class,
    crownrad_m, agbiomass_kg, bolebiomass_kg, foliagebiomass_kg, crownbiomass_kg) %>%

  # calculate some standard summary statistics values
  mutate(
    tph = ifelse(dbh_cm>=11.4,
                 1/0.09,
                 1/((2*30+(2*28))/10000)),
    ba_m2 = pi*((dbh_cm/2/100)**2),
    ba_m2ha = ba_m2 * tph,
    agbiomass_mgha = 
      (agbiomass_kg/1000)*tph
  ) %>%
  
  # calculate intercept a and slope b for radius_m ~ a + b * height_m for 
  # both the bole and crown
  mutate(
    a_bole = (-dbh_cm*height_m)/(200*(1.37-height_m)),
    b_bole = dbh_cm / (200*(1.37-height_m)),
    a_crown = 
      ifelse(status=='L',
             -(crownrad_m*height_m)/(htcb_m-height_m),
             0),
    b_crown = 
      ifelse(status=='L',
             crownrad_m / (htcb_m-height_m),
             0)
  ) %>%
  mutate(bolevol_m3 = 
           pi*((a_bole+b_bole*0)**2)*(height_m/3),
         crownvol_m3 = 
           
           ifelse(is.element(spp, c('ABCO', 'CADE', 'PILA', 'PIPO')),
                  # cone for conifers
                  pi*(crownrad_m**2)*((height_m-htcb_m)/3),
                  # cylinder for hardwoods
                  pi*(crownrad_m**2)*((height_m-htcb_m))),
         
         boledensity_kgm3 = 
           bolebiomass_kg / bolevol_m3,
         
         crowndensity_kgm3 = 
           crownbiomass_kg / crownvol_m3)


ggplot(combo_sim_treelist,
       aes(x = x_coord, y = y_coord, color = spp, size = dbh_cm))+
  geom_point()+
  coord_fixed()+
  theme_minimal()

#### biomass of understory pixels ##############################################

head(vegpa_sim_df)


# pixel size on the veg height simulation is set such that 1 pixel =~ 1 plant
# use the generic all species formula from mcginnis et al. 2010 to estimate 
# the total biomass, biomass of foliage, and biomass of 1-hour woody branches 
# for each plant using an assumed total crown volume (assuming crown dimensions
# are 0.5x0.5xveg_height meters^3)
# I'm getting really high biomass numbers using the crown volume approach, 
# which seem too high based on the (admittedly scarce) literature. So instead 
# i'm going to estimate volume from the equation for crown diameter and height, 
# assuming a crown radius r such that pi*r^2 = 0.5^2 so that the crown area 
# of one individual is one pixel (where we picked the pixel size to have 
# approximately the area of the midpoint individual from the mcginnis data)
# so r = sqrt(0.5^2 / pi) = 0.282 so diameter = 0.56m

mcginnis_t1 = 
  read.csv(here::here('02-data',
                      '00-source',
                      'mcginnis',
                      'mcginnis_t1.csv')) %>%
  as_tibble()


# GO BACK AND ADD THE CORRECTION FACTORS
vegpa_sim_df = 
  vegpa_sim_df %>%
  mutate(crown_vol_m3 = 0.5*0.5*height_m) %>%
  mutate(spp_match = 
           ifelse(spp=='RIRO',
                  'RIB-',
                  ifelse(is.element(spp, mcginnis_t1$species),
                         spp,
                         'ALL'))) %>%
  left_join(mcginnis_t1 %>%
              filter(stock == 'total' & equation == 'cvol_m3'),
            by = c('spp_match' = 'species')) %>%
  mutate(biomass_g = exp(a+(b1*log(crown_vol_m3)))*cf,
         biomass_kgm3 = (biomass_g/1000)/crown_vol_m3) %>%
  select(group_id, plot_id, x_coord, y_coord, spp, status, height_m, 
         crown_vol_m3, biomass_g, biomass_kgm3, present = Y)
  
head(vegpa_sim_df)

ggplot(vegpa_sim_df,
       aes(x = x_coord, y = y_coord, fill = height_m))+
  geom_tile()+
  coord_fixed()+
  theme_minimal()+
  scale_fill_viridis_c()

vegpa_sim_df %>% 
  summarise(biomass_Mgha = (sum(biomass_g, na.rm = TRUE)/1000000)*(0.09))

#### voxelize trees  ###########################################################


# resolution cannot be finer than the resolution of the understory veg pixels
resolution = 1
fine_scale = 10
coarse_x_coords = seq(from = resolution/2, to = 30-(resolution/2), by = resolution)
coarse_y_coords = seq(from = resolution/2, to = 30-(resolution/2), by = resolution)
coarse_z_coords = 
  seq(from = resolution/2, 
      to = ceiling(max(combo_sim_treelist$height_m)/resolution)*resolution-(resolution/2), 
      by = resolution)



build_bulkdensity_array = 
  function(trees,
           x_coords,
           y_coords,
           z_coords){
    
    array(dim = c(length(x_coords),
                  length(y_coords),
                  length(z_coords)),
          dimnames = list('x' = x_coords, 'y' = y_coords, 'z' = z_coords),
          
          data = 
            sapply(X = z_coords,
                   FUN = function(z){
                   # calculate the bole radius of each tree at this z,
                    # eqn goes negative above height of each tree, 
                    # so multiply by 0 if z is greater than each trees height
                    bole_radii = 
                      (trees$a_bole+(trees$b_bole*z))*(z <= trees$height_m)
                      
                    # calculate the crown radius of each tree at this z
                    # multiply by 0 if z is greater than the tree's height,
                    # or below the trees htcb
                    crown_radii = 
                      (trees$a_crown+(trees$b_crown*z))*(z>=trees$htcb_m & z<= trees$height_m)
                    # set NAs to 0 (for snags)
                    #crown_radii[is.na(crown_radii)] = 0

                    # for every y
                    sapply(X = y_coords,
                          FUN = function(y){
                            
                            # for every x
                            sapply(X = x_coords,
                                   FUN = function(x){
                                     
                                     # get the distance between this 
                                     # point and each tree
                                     horiz_distance = 
                                       sqrt((x-trees$x_coord)**2+
                                              (y-trees$y_coord)**2)
                                     
                                     # get the (within plot) indices
                                     # of the trees which this point 
                                     # is within the crown
                                     within_crown = 
                                       (horiz_distance < crown_radii)
                                     
                                     within_bole = 
                                       (horiz_distance < bole_radii)
                                     
                                     # get the contribution of each 
                                     # tree to the total bulk density
                                     bulkdensity_kgm3 = 
                                       (within_bole*trees$boledensity_kgm3)+
                                       (within_crown*trees$crowndensity_kgm3)
                                     
                                     # sum all the tree's contributions
                                     # to get the overall bulk density 
                                     # for this voxel; nas are from snags which 
                                     # have NA crown density
                                     return(sum(bulkdensity_kgm3, na.rm = TRUE))
                                     
                                   })
                          })
                  })
    )}
           



bulkdensity_shrubs = 
  array(dim = c(length(coarse_x_coords),
                length(coarse_y_coords),
                length(coarse_z_coords),
                length(unique(vegpa_sim_df$plot_id))),
        dimnames = list('x' = coarse_x_coords,
                        'y' = coarse_y_coords,
                        'z' = coarse_z_coords),
        data = 
          sapply(X = unique(vegpa_sim_df$plot_id),
                 FUN = function(i){
                   
                   understory = vegpa_sim_df %>% filter(plot_id==i)
                   
                   understory_res = 
                     abs(unique(understory$x_coord)[order(unique(understory$x_coord))][2] -
                           unique(understory$x_coord)[order(unique(understory$x_coord))][1])
                     
                   max_veg_height = 
                     floor(max(understory$height_m,na.rm = TRUE)/resolution)*resolution+
                     (resolution/2)
                   
                   x_fine = seq(0+(resolution/fine_scale/2), 
                                30-resolution/fine_scale+(resolution/fine_scale/2), 
                                by = resolution/fine_scale)
                   y_fine = seq(0+(resolution/fine_scale/2), 
                                30-resolution/fine_scale+(resolution/fine_scale/2), 
                                by = resolution/fine_scale)
                    x_understory = 
                     floor(x_fine/understory_res)*understory_res+
                     understory_res/2
                 
                   y_understory = 
                    floor(y_fine/understory_res)*understory_res+
                     understory_res/2
                   
                   x_coarse = 
                     floor(x_fine/resolution)*resolution+
                     resolution/2
                   y_coarse = 
                     floor(x_fine/resolution)*resolution+
                     resolution/2
                   
                   fine_veg_mapping = 
                     data.frame(x_fine, x_understory,
                                x_coarse,
                                y_fine, y_understory,
                                y_coarse) %>%
                     expand(nesting(x_fine, x_understory, x_coarse),
                            nesting(y_fine, y_understory, y_coarse)) %>%
                     left_join(understory,
                               by = c('x_understory' = 'x_coord',
                                      'y_understory' = 'y_coord'))
                                       
                   
                   sapply(X = coarse_z_coords,
                          FUN = function(z_coarse){
                            
                            if (z_coarse > max_veg_height){
                              zero_array = 
                                array(dim = c(length(coarse_x_coords),
                                              length(coarse_y_coords),
                                              1),
                                      dimnames = list('x' = coarse_x_coords,
                                                      'y' = coarse_y_coords,
                                                      'z' = z_coarse),
                                      data = rep(0, 
                                                 times = 
                                                   length(coarse_x_coords)*
                                                   length(coarse_y_coords)))
                              
                              return(zero_array)
                            }
                            
                            # define the fine-resolution zs to iterate over
                            z_fine = seq(from = z_coarse-(resolution/2)+
                                           (resolution/fine_scale/2),
                                         to = z_coarse + (resolution/2)-
                                           resolution/fine_scale/2,
                                         by = resolution/fine_scale)
                            
                            fine_bd_slices =
                              
                              array(dim = c(length(coarse_x_coords),
                                            length(coarse_y_coords),
                                            length(z_fine)),
                                    dimnames = list('x' = coarse_x_coords,
                                                    'y' = coarse_y_coords,
                                                    'z' = z_fine),
                                    
                                    data = 
                                      sapply(X = z_fine,
                                             FUN = function(z){
                                               
                                               coarse_veg_mapping = 
                                                 fine_veg_mapping %>%
                                                 mutate(bulkdensity_kgm3 = 
                                                          ifelse(!is.na(biomass_kgm3)&
                                                                   !is.nan(biomass_kgm3),
                                                                 ifelse(height_m >= z,
                                                                        biomass_kgm3,
                                                                        0),
                                                                 0)) %>%
                                                 group_by(x_coarse, y_coarse) %>%
                                                 summarise(bulkdensity_kgm3 = 
                                                             mean(bulkdensity_kgm3)) %>%
                                                 ungroup() %>%
                                                 arrange(y_coarse, x_coarse)
                                         
                                               sapply(X = 1:length(coarse_y_coords),
                                                      FUN = function(y_i){
                                                        sapply(X = 1:length(coarse_x_coords),
                                                               FUN = function(x_i){
                                                                 
                                                                 coarse_veg_mapping$bulkdensity_kgm3[
                                                                   x_i+(length(coarse_x_coords)*(y_i-1))
                                                                        ]
                                                               })
                                                      })
                                               
                                             })
                                    
                                    
                                    )
                            
                            return(apply(X = fine_bd_slices,
                                         MARGIN = c(1,2),
                                         FUN = mean))

                            
                            })
                   
                   
                 }))

# more memory efficient version?
bulkdensity_trees = 
  
  array(dim = c(length(coarse_x_coords),
                length(coarse_y_coords),
                length(coarse_z_coords),
                length(unique(trees_sim_df$plot_id))),
        
        dimnames = 
          list('x' = coarse_x_coords,
               'y' = coarse_y_coords,
               'z' = coarse_z_coords,
               'plot_id' = unique(trees_sim_df$plot_id)),
        
        data = 
          
          # for every plot
          sapply(X = unique(trees_sim_df$plot_id),
                 FUN = function(i){
                   
                   # get the trees on this plot
                   treelist = combo_sim_treelist %>% filter(plot_id==i)
                   
                   
                   # for every z
                   sapply(X = coarse_z_coords,
                          FUN = function(z_coarse){
                            
                            
                            
                            # define the fine-resolution zs to iterate over
                            z_fine = seq(from = 
                                           z_coarse-(resolution/2)+(resolution/fine_scale/2),
                                         to = z_coarse + (resolution/2)-
                                           (resolution/fine_scale/2),
                                         by = resolution/fine_scale)
                            
                            # for every y
                            sapply(X = coarse_y_coords,
                                   FUN = function(y_coarse){
                                     
                                     
                                     # define the fineres ys to iterate over
                                     y_fine = seq(from = 
                                                    y_coarse-(resolution/2)+
                                                    (resolution/fine_scale/2),
                                                  to = y_coarse+(resolution/2)-
                                                    (resolution/fine_scale/2),
                                                  by = resolution/fine_scale)
                                     
                                     # for every x
                                     sapply(X = coarse_x_coords,
                                            FUN = function(x_coarse){
                                              
                                              
                                              x_fine = seq(x_coarse-(resolution/2)+
                                                             (resolution/fine_scale/2),
                                                           x_coarse+(resolution/2)-
                                                             (resolution/fine_scale/2),
                                                           by = resolution/fine_scale)
                                              
                                              print(paste0('Working on (',
                                                           x_coarse,
                                                           ', ',
                                                           y_coarse,
                                                           ', ',
                                                           z_coarse,
                                                           ')'))
                                              
                                              # define a fine-resolution array of subvoxels
                                              array_fine = 
                                                build_bulkdensity_array(trees = treelist,
                                                                        x_coords = x_fine,
                                                                        y_coords = y_fine,
                                                                        z_coords = z_fine)
                                              
                                              # the bulk density for this voxel is the 
                                              # mean of the fine resolution array
                                              result = mean(array_fine)
                                              
                                              return(result)
                                              
                                              
                                            })
                                     
                                   })
                          })
                 })
  )


bulkdensity_shrubs.df = 
  as.data.table(bulkdensity_shrubs) %>%
  as_tibble() %>%
  mutate(x = as.numeric(x),
         y = as.numeric(y),
         z = as.numeric(z)) %>%
  rename(plot_id = V1,
         bulkdensity_kgm3 = value)

head(bulkdensity_shrubs.df)

bulkdensity_trees.df = 
  as.data.table(bulkdensity_trees) %>%
  as_tibble() %>%
  mutate(x = as.numeric(x),
         y = as.numeric(y),
         z = as.numeric(z)) %>%
  rename(bulkdensity_kgm3 = value)

head(bulkdensity_trees.df)

#### make 2d figures ###########################################################

duff_depth_surface = 
  ggplot(duff_sim_df,
         aes(x = x_coord, y = y_coord, fill = duff_cm))+
  geom_tile()+
  coord_fixed()+
  theme_minimal()+
  scale_fill_viridis_c(direction = -1, option = 'D')+
  labs(fill = 'Duff depth (cm)',
       x = 'Easting (m)',
       y = 'Northing (m)')+
  theme(legend.position = 'bottom')

duff_depth_surface

ggsave(duff_depth_surface+
         theme(text = element_text(size = 16)),
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'duff_depth_surface.png'),
       height = 6.5, width = 6.5, units = 'in')

duff_load_surface = 
  ggplot(duff_sim_df,
         aes(x = x_coord, y = y_coord, fill = duff_kgm2))+
  geom_tile()+
  coord_fixed()+
  theme_minimal()+
  scale_fill_viridis_c(direction = -1, option = 'D')+
  labs(fill = 'Duff loading (kg/m^2)',
       x = 'Easting (m)',
       y = 'Northing (m)')+
  theme(legend.position = 'bottom')

duff_load_surface

ggsave(duff_load_surface+
         theme(text = element_text(size = 16)),
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'duff_load_surface.png'),
       height = 6.5, width = 6.5, units = 'in')


litter_depth_surface = 
  ggplot(litter_sim_df,
         aes(x = x_coord, y = y_coord, fill = litter_cm))+
  geom_tile()+
  coord_fixed()+
  theme_minimal()+
  scale_fill_viridis_c(direction = -1, option = 'D')+
  labs(fill = 'Litter depth (cm)',
       x = 'Easting (m)',
       y = 'Northing (m)')+
  theme(legend.position = 'bottom')

litter_depth_surface

ggsave(litter_depth_surface+
         theme(text = element_text(size = 16)),
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'litter_depth_surface.png'),
       height = 6.5, width = 6.5, units = 'in')


litter_load_surface = 
  ggplot(litter_sim_df,
         aes(x = x_coord, y = y_coord, fill = litter_kgm2))+
  geom_tile()+
  coord_fixed()+
  theme_minimal()+
  scale_fill_viridis_c(direction = -1, option = 'D')+
  labs(fill = 'Litter loading (kg/m^2)',
       x = 'Easting (m)',
       y = 'Northing (m)')+
  theme(legend.position = 'bottom')

litter_load_surface

ggsave(litter_load_surface+
         theme(text = element_text(size = 16)),
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'litter_load_surface.png'),
       height = 6.5, width = 6.5, units = 'in')

fwd_count_surface = 
  ggplot(fwd_sim_df %>%
           group_by(group_id, plot_id, x_coord, y_coord) %>%
           summarise(count = sum(count)) %>%
           ungroup(),
         aes(x = x_coord, y = y_coord, fill = count))+
  geom_tile()+
  coord_fixed()+
  theme_minimal()+
  scale_fill_viridis_c(direction = -1, option = 'D')+
  labs(fill = 'FWD count',
       x = 'Easting (m)',
       y = 'Northing (m)')+
  theme(legend.position = 'bottom')

fwd_count_surface

ggsave(fwd_count_surface+
         theme(text = element_text(size = 16)),
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'fwd_count_surface.png'),
       height = 6.5, width = 6.5, units = 'in')


fwd_diams_histogram = 
  ggplot(fwddiams_sim_data_long,
         aes(x = diam_cm))+
  geom_histogram()+
  geom_vline(lty = 2, xintercept = 0.635, lwd = 1)+
  geom_vline(lty = 2, xintercept = 2.54, lwd = 1)+
  theme_minimal()+
  labs(x = 'Particle diameter (cm)', y = 'Count')

fwd_diams_histogram

ggsave(fwd_diams_histogram+
         theme(text = element_text(size = 16)),
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'fwd_diams_histogram.png'),
       height = 6.5, width = 6.5, units = 'in')


fwd_load_surfaces = 
  ggplot(fwd_sim_df,
         aes(x = x_coord, y = y_coord, fill = fwd_kgm2))+
  geom_tile()+
  coord_fixed()+
  theme_minimal()+
  scale_fill_viridis_c(direction = -1, option = 'D')+
  labs(fill = 'FWD Loading (kg/m^2)',
       x = 'Easting (m)',
       y = 'Northing (m)')+
  facet_wrap(timelag_class~.,
             ncol = 2)

fwd_load_surfaces

ggsave(fwd_load_surfaces+
         theme(text = element_text(size = 16)),
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'fwd_load_surface.png'),
       height = 6.5, width = 6.5, units = 'in')

cwd_count_surface = 
  ggplot(cwd_sim_df %>%
           group_by(group_id, plot_id, x_coord, y_coord) %>%
           summarise(count = sum(count)) %>%
           ungroup(),
         aes(x = x_coord, y = y_coord, fill = count))+
  geom_tile()+
  coord_fixed()+
  theme_minimal()+
  scale_fill_viridis_c(direction = -1, option = 'D')+
  labs(fill = 'CWD count',
       x = 'Easting (m)',
       y = 'Northing (m)')+
  theme(legend.position = 'bottom')

cwd_count_surface

ggsave(cwd_count_surface+
         theme(text = element_text(size = 16)),
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'cwd_count_surface.png'),
       height = 6.5, width = 6.5, units = 'in')


cwd_diams_histogram = 
  ggplot(cwddiams_sim_data_long,
         aes(x = diam_cm))+
  geom_histogram()+
  theme_minimal()+
  labs(x = 'Particle diameter (cm)', y = 'Count')

cwd_diams_histogram

ggsave(cwd_diams_histogram+
         theme(text = element_text(size = 16)),
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'cwd_diams_histogram.png'),
       height = 6.5, width = 6.5, units = 'in')


cwd_load_surfaces = 
  ggplot(cwd_sim_df,
         aes(x = x_coord, y = y_coord, fill = cwd_kgm2))+
  geom_tile()+
  coord_fixed()+
  theme_minimal()+
  scale_fill_viridis_c(direction = -1, option = 'D')+
  labs(fill = 'cwd Loading (kg/m^2)',
       x = 'Easting (m)',
       y = 'Northing (m)')+
  theme(legend.position = 'bottom')

cwd_load_surfaces

ggsave(cwd_load_surfaces+
         theme(text = element_text(size = 16)),
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'cwd_load_surface.png'),
       height = 6.5, width = 6.5, units = 'in')


veg_presence_surface = 
  ggplot(vegpa_sim_df,
         aes(x = x_coord, y = y_coord, fill = as.logical(present)))+
  geom_tile()+
  coord_fixed()+
  theme_minimal()+
  scale_fill_viridis_d(direction = -1, option = 'D')+
  labs(fill = 'Veg Present',
       x = 'Easting (m)',
       y = 'Northing (m)')+
  theme(legend.position = 'bottom')

veg_presence_surface

ggsave(veg_presence_surface+
         theme(text = element_text(size = 16)),
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'veg_presence_surface.png'),
       height = 6.5, width = 6.5, units = 'in')


trees_count_surface = 
  ggplot(trees_sim_df,
         aes(x = x_coord, y = y_coord, fill = as.integer(Y)))+
  geom_tile()+
  coord_fixed()+
  theme_minimal()+
  scale_fill_viridis_c(direction = -1, option = 'D')+
  labs(fill = 'Tree count',
       x = 'Easting (m)',
       y = 'Northing (m)')+
  theme(legend.position = 'bottom')

trees_count_surface

ggsave(trees_count_surface+
         theme(text = element_text(size = 16)),
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'trees_count_surface.png'),
       height = 6.5, width = 6.5, units = 'in')


saplings_count_surface = 
  ggplot(saplings_sim_df,
         aes(x = x_coord, y = y_coord, fill = Y))+
  geom_tile()+
  coord_fixed()+
  theme_minimal()+
  scale_fill_viridis_c(direction = -1, option = 'D')+
  labs(fill = 'Saplings count',
       x = 'Easting (m)',
       y = 'Northing (m)')+
  theme(legend.position = 'bottom')

saplings_count_surface

ggsave(saplings_count_surface+
         theme(text = element_text(size = 16)),
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'saplings_count_surface.png'),
       height = 6.5, width = 6.5, units = 'in')


combo_stem_map = 
  ggplot(combo_sim_treelist,
         aes(x = x_coord, y = y_coord, color = spp, size = dbh_cm))+
  geom_point()+
  coord_fixed(xlim = c(0, 30), ylim = c(0, 30))+
  theme_minimal()+
  labs(color = 'Species',
       size = 'DBH (cm)',
       x = 'Easting (m)',
       y = 'Northing (m)')

combo_stem_map

ggsave(combo_stem_map+
         theme(text = element_text(size = 16)),
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'combo_stem_map.png'),
       height = 6.5, width = 6.5, units = 'in')


#### make 3d figures ###########################################################

# make nice viridis color pallette
map_viridis = 
  function(x, res, begin = 0, end = 1, opt = 'D'){
    
    # define the sequence
    minval = min(x)
    maxval = max(x)
    
    scaled_values = 
      as.integer(round((x-minval) / (maxval-minval), res) * 10 * res)
    
    viridis_colors = viridis(n = 10*res+1, begin = begin, end = end, option = opt)
    
    return(viridis_colors[scaled_values+1])
    
  }

#rgl.open()
rgl::open3d()
#rgl.bg(color = 'white')

# axes
segments3d(x = c(60, 60), z = c(0, 0), y = c(60, 0), color = 'black')
segments3d(x = c(60, 60), z = c(0, 130), y = c(60, 60), color = 'black')
segments3d(x = c(60, 0), z = c(0, 0), y = c(60, 60), color = 'black')


# axes labels
rgl.texts(x = 60, y = -5, z = 0, text = '30m', color = 'black')
rgl.texts(z = 0, y = 65, x = 0, text = '30m', color = 'black')
rgl.texts(x = 60, z= 130, y= 60, text = '65m', color = 'black')

veg_voxels_mesh = 
  VoxR::plot_voxels(bulkdensity_shrubs.df %>%
                      filter(bulkdensity_kgm3>0), 
                    plot = FALSE, 
                    alpha = 0.5,
                    res = 1)

veg_colors = map_viridis(veg_voxels_mesh$additionnal$bulkdensity_kgm3,
                     res = 2, begin = 1, end = 0)

#rgl::open3d()
rgl::shade3d(veg_voxels_mesh$mesh,
             col = veg_colors,
             #alpha = 0.5,
             lit = FALSE)
rgl::wire3d(veg_voxels_mesh$mesh,
             #col = veg_colors,
             alpha = 0.5,
             lit = FALSE)

snapshot3d(here::here('04-communication',
                      'figures',
                      'powerpoint',
                      'understory_voxels.png'),
           width = 900, height = 900)

rgl.clear()

rgl.close()

rgl::open3d()
#rgl.bg(color = 'white')

# axes
segments3d(x = c(60, 60), z = c(0, 0), y = c(60, 0), color = 'black')
segments3d(x = c(60, 60), z = c(0, 130), y = c(60, 60), color = 'black')
segments3d(x = c(60, 0), z = c(0, 0), y = c(60, 60), color = 'black')


# axes labels
rgl.texts(x = 60, y = -5, z = 0, text = '30m', color = 'black')
rgl.texts(z = 0, y = 65, x = 0, text = '30m', color = 'black')
rgl.texts(x = 60, z= 130, y= 60, text = '65m', color = 'black')

trees_voxels_mesh = 
  VoxR::plot_voxels(bulkdensity_trees.df %>%
                      filter(bulkdensity_kgm3>0), 
                    plot = FALSE, 
                    alpha = 0.5,
                    res = 1)

trees_colors = map_viridis(trees_voxels_mesh$additionnal$bulkdensity_kgm3,
                     res = 2, begin = 1, end = 0)

#rgl::open3d()
rgl::shade3d(trees_voxels_mesh$mesh,
             col = trees_colors,
             #alpha = 0.1,
             lit = TRUE)
rgl::wire3d(trees_voxels_mesh$mesh,
             #col = trees_colors,
             alpha = 0.5,
             lit = TRUE)


snapshot3d(here::here('04-communication',
                      'figures',
                      'powerpoint',
                      'overstory_voxels.png'),
           width = 900, height = 900)

rgl.clear()

rgl.close()

#### write data products #######################################################

write.csv(
  duff_sim_df,
  here::here('02-data',
             '06-results',
             'simulated_fuelbeds',
             'duff.csv'),
  row.names = FALSE
)


write.csv(
  litter_sim_df,
  here::here('02-data',
             '06-results',
             'simulated_fuelbeds',
             'litter.csv'),
  row.names = FALSE
)


write.csv(
  fwd_sim_df,
  here::here('02-data',
             '06-results',
             'simulated_fuelbeds',
             'fwd.csv'),
  row.names = FALSE
)


write.csv(
  cwd_sim_df,
  here::here('02-data',
             '06-results',
             'simulated_fuelbeds',
             'cwd.csv'),
  row.names = FALSE
)

write.csv(
  bulkdensity_shrubs.df,
  here::here('02-data',
             '06-results',
             'simulated_fuelbeds',
             'bulkdensity_understory.csv'),
  row.names = FALSE
)


write.csv(
  bulkdensity_trees.df,
  here::here('02-data',
             '06-results',
             'simulated_fuelbeds',
             'bulkdensity_treessaplings.csv'),
  row.names = FALSE
)

