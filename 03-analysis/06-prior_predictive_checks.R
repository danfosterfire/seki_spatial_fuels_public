

#### setup #####################################################################

library(here)
library(tidyverse)


#### litter and duff ###########################################################

# just to get the data structure, not actually looking at the real data yet
litterduff_data = 
  list(
    G = 1,
    P = 10,
    L = readRDS(here::here('02-data','05-for_analysis','litter_training_data.rds'))$L,
    coords = readRDS(here::here('02-data','05-for_analysis','litter_training_data.rds'))$coords

  )

litterduff_data$N = litterduff_data$P*litterduff_data$L
litterduff_data$location_id = rep(1:litterduff_data$L, times = litterduff_data$P)
litterduff_data$plot_id = rep(1:litterduff_data$P, each = litterduff_data$L)
litterduff_data$group_id = rep(1, times = litterduff_data$P)

litterduff_prior = 
  data.frame(
    'beta' = rnorm(n = 4000, mean = 0, sd = 2),
    'alpha' = truncnorm::rtruncnorm(n = 4000, mean = 0, sd = 2, a = 0),
    'rho' = invgamma::rinvgamma(n = 4000, shape = 5, rate = 40),
    'sigmaPlot' = truncnorm::rtruncnorm(n = 4000, mean = 0, sd = 2, a = 0),
    'kappa' = truncnorm::rtruncnorm(n = 4000, mean = 0, sd = 2, a = 0),
    draw = 1:4000
  )

set.seed(110819)
litterduff_sim = 
  do.call('bind_rows',
          lapply(X = 1:100,
                 FUN = function(draw){
                   
                   # extract parameters for this draw
                   beta = 
                     litterduff_prior %>%
                     select(contains('beta')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   alpha = 
                     litterduff_prior %>%
                     select(contains('alpha')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   rho = 
                     litterduff_prior %>%
                     select(contains('rho')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   sigmaPlot = 
                     litterduff_prior %>%
                     select(contains('sigmaPlot')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   plotEffects = 
                     rnorm(n = litterduff_data$P,
                           mean = 0,
                           sd = sigmaPlot)
                     
                   zGP = 
                     matrix(nrow = litterduff_data$L,
                            ncol = litterduff_data$P,
                            byrow = FALSE,
                            data = 
                              rnorm(n = 
                                      litterduff_data$L*
                                      litterduff_data$P,
                                    mean = 0,
                                    sd = 1))
                   
                   kappa = 
                     litterduff_prior %>%
                     select(contains('kappa')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   dists = 
                     litterduff_data$coords %>%
                     dist() %>%
                     as.matrix()
                   
                   # construct GP effects
                   K = 
                     array(dim = 
                             c(litterduff_data$L,
                               litterduff_data$L,
                               litterduff_data$G),
                           dimnames = 
                             list('l' = 1:litterduff_data$L,
                                  'lprime' = 1:litterduff_data$L,
                                  'g' = 1:litterduff_data$G),
                           data = 
                             sapply(X = 1:litterduff_data$G,
                                    FUN = function(g){
                                      sapply(X = 1:litterduff_data$L,
                                             FUN = function(l){
                                               sapply(X = 1:litterduff_data$L,
                                                      FUN = function(lprime){
                                                        ((alpha[g]**2)*
                                                          exp(-(dists[lprime,l]**2)/
                                                                (2*rho[g]**2)))
                                                      })
                                             })
                                    }))
                   
                   for (g in 1:litterduff_data$G){
                     K[,,g] = K[,,g]+diag(nrow = nrow(K[,,g]),
                                          ncol = ncol(K[,,g]),
                                          x = 1e-9)
                   }
                   
                   KL = 
                     array(dim = 
                             c(litterduff_data$L,
                               litterduff_data$L,
                               litterduff_data$G),
                           dimnames = 
                             list('l' = 1:litterduff_data$L,
                                  'lprime' = 1:litterduff_data$L,
                                  'g' = 1:litterduff_data$G),
                           data = 
                             sapply(X = 1:litterduff_data$G,
                                    FUN = function(g){
                                      t(chol(K[,,g]))
                                    }))
                 
                   # realized random effects
                   GP = 
                     matrix(nrow = litterduff_data$L,
                            ncol = litterduff_data$P,
                            data = 
                              sapply(X = 1:litterduff_data$P,
                                     FUN = function(p){
                                       KL[,,litterduff_data$group_id[p]] %*%
                                         zGP[,p]
                                     }))
                   
                   # linear predictors
                   logmu = 
                     sapply(X = 1:litterduff_data$N,
                            FUN = function(i){
                              
                               beta[litterduff_data$group_id[litterduff_data$plot_id[i]]] +
                                plotEffects[litterduff_data$plot_id[i]]+
                                as.numeric(GP[litterduff_data$location_id[i],
                                              litterduff_data$plot_id][i])
                              
                            })
                    
                   # realizations
                   Y = 
                     sapply(X = 1:litterduff_data$N,
                            FUN = function(i){
                              rnbinom(n = 1,
                                      mu = exp(logmu[i]),
                                      size = 
                                        kappa[litterduff_data$group_id[litterduff_data$plot_id[i]]])
                            })
                   
                   
                   results = 
                     data.frame(Y_sim = Y,
                                group_id = 
                                  litterduff_data$group_id[litterduff_data$plot_id],
                                plot_id = litterduff_data$plot_id,
                                location_id = litterduff_data$location_id,
                                logmu = logmu,
                                draw = draw) %>%
                     as_tibble() %>%
                     rowid_to_column('obs_id')
                   
                   return(results)
                 }))

# plot mean as a function of beta
litterduff_sim %>%
  left_join(litterduff_prior) %>%
  ggplot(aes(x = Y_sim))+
  geom_histogram()+
  scale_x_log10()+
  theme_minimal()

#### fwd tallies ###############################################################


fwdtallies_data = 
  list(
    G = 1,
    P = 10,
    L =   readRDS(here::here('02-data', 
                     '05-for_analysis',
                     'fwdtallies_training_data.rds'))$L,
    coords =   readRDS(here::here('02-data', 
                     '05-for_analysis',
                     'fwdtallies_training_data.rds'))$coords
  )

fwdtallies_data$N = fwdtallies_data$P*fwdtallies_data$L
fwdtallies_data$location_id = rep(1:fwdtallies_data$L, times = fwdtallies_data$P)
fwdtallies_data$plot_id = rep(1:fwdtallies_data$P, each = fwdtallies_data$L)
fwdtallies_data$group_id = rep(1, times = fwdtallies_data$P)

fwdtallies_prior = 
  data.frame(
    'beta' = rnorm(n = 4000, mean = 0, sd = 2),
    'alpha' = truncnorm::rtruncnorm(n = 4000, mean = 0, sd = 2, a = 0),
    'rho' = invgamma::rinvgamma(n = 4000, shape = 5, rate = 40),
    'sigmaPlot' = truncnorm::rtruncnorm(n = 4000, mean = 0, sd = 2, a = 0),
    'tau' = truncnorm::rtruncnorm(n = 4000, mean = 0, sd = 2, a = 0),
    'kappa' = truncnorm::rtruncnorm(n = 4000, mean = 0, sd = 2, a = 0),
    draw = 1:4000
  )
  

set.seed(110819)
fwdtallies_sim = 
  do.call('bind_rows',
          lapply(X = 1:100,
                 FUN = function(draw){
                   
                   # extract parameters for this draw
                   beta = 
                     fwdtallies_prior %>%
                     select(contains('beta')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   alpha = 
                     fwdtallies_prior %>%
                     select(contains('alpha')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   rho = 
                     fwdtallies_prior %>%
                     select(contains('rho')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   sigmaPlot = 
                     fwdtallies_prior %>%
                     select(contains('sigmaPlot')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   plotEffects = 
                     rnorm(n = fwdtallies_data$P,
                           mean = 0,
                           sd = sigmaPlot)
                   
                   zGP = 
                     matrix(nrow = fwdtallies_data$L,
                            ncol = fwdtallies_data$P,
                            byrow = FALSE,
                            data = 
                              rnorm(n = fwdtallies_data$L*fwdtallies_data$P,
                                    mean = 0,
                                    sd = 1))
                   
                   tau = 
                     fwdtallies_prior %>%
                     select(contains('tau')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   kappa = 
                     fwdtallies_prior %>%
                     select(contains('kappa')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   dists = 
                     fwdtallies_data$coords %>%
                     dist() %>%
                     as.matrix()
                   
                   # construct GP effects
                   K = 
                     array(dim = 
                             c(fwdtallies_data$L,
                               fwdtallies_data$L,
                               fwdtallies_data$G),
                           dimnames = 
                             list('l' = 1:fwdtallies_data$L,
                                  'lprime' = 1:fwdtallies_data$L,
                                  'g' = 1:fwdtallies_data$G),
                           data = 
                             sapply(X = 1:fwdtallies_data$G,
                                    FUN = function(g){
                                      sapply(X = 1:fwdtallies_data$L,
                                             FUN = function(l){
                                               sapply(X = 1:fwdtallies_data$L,
                                                      FUN = function(lprime){
                                                        ((alpha[g]**2)*
                                                          exp(-(dists[lprime,l]**2)/
                                                                (2*rho[g]**2)))
                                                      })
                                             })
                                    }))
                   
                   for (g in 1:fwdtallies_data$G){
                     K[,,g] = K[,,g]+diag(nrow = nrow(K[,,g]),
                                          ncol = ncol(K[,,g]),
                                          x = tau[g])
                   }
                   
                   KL = 
                     array(dim = 
                             c(fwdtallies_data$L,
                               fwdtallies_data$L,
                               fwdtallies_data$G),
                           dimnames = 
                             list('l' = 1:fwdtallies_data$L,
                                  'lprime' = 1:fwdtallies_data$L,
                                  'g' = 1:fwdtallies_data$G),
                           data = 
                             sapply(X = 1:fwdtallies_data$G,
                                    FUN = function(g){
                                      t(chol(K[,,g]))
                                    }))
                 
                   # realized random effects
                   GP = 
                     matrix(nrow = fwdtallies_data$L,
                            ncol = fwdtallies_data$P,
                            data = 
                              sapply(X = 1:fwdtallies_data$P,
                                     FUN = function(p){
                                       KL[,,fwdtallies_data$group_id[p]] %*%
                                         zGP[,p]
                                     }))
                   
                   # linear predictors
                   logmu = 
                     sapply(X = 1:fwdtallies_data$N,
                            FUN = function(i){
                              
                               beta[fwdtallies_data$group_id[fwdtallies_data$plot_id[i]]] +
                                plotEffects[fwdtallies_data$plot_id[i]]+
                                as.numeric(GP[fwdtallies_data$location_id[i],
                                              fwdtallies_data$plot_id][i])
                              
                            })
                    
                   # realizations
                   Y = 
                     sapply(X = 1:fwdtallies_data$N,
                            FUN = function(i){
                              rnbinom(n = 1,
                                      mu = exp(logmu[i]),
                                      size = 
                                        kappa[fwdtallies_data$group_id[fwdtallies_data$plot_id[i]]])
                            })
                   
                   
                   results = 
                     data.frame(Y_sim = Y,
                                group_id = 
                                  fwdtallies_data$group_id[fwdtallies_data$plot_id],
                                plot_id = fwdtallies_data$plot_id,
                                location_id = fwdtallies_data$location_id,
                                logmu = logmu,
                                draw = draw) %>%
                     as_tibble() %>%
                     rowid_to_column('obs_id')
                   
                   return(results)
                 }))


fwdtallies_sim %>%
  left_join(fwdtallies_prior) %>%
  ggplot(aes(x = Y_sim))+
  geom_histogram()+
  scale_x_log10()+
  theme_minimal()


#### fwd diams #################################################################


fwddiams_data = 
  list(
    N = readRDS(here::here('02-data',
                           '05-for_analysis',
                           'fwddiams_training_data.rds'))$N,
    G = 1,
    P = readRDS(here::here('02-data',
                           '05-for_analysis',
                           'fwddiams_training_data.rds'))$P,
    TS = readRDS(here::here('02-data',
                           '05-for_analysis',
                           'fwddiams_training_data.rds'))$TS,
    group_id = rep(1, 
                   times = readRDS(here::here('02-data',
                           '05-for_analysis',
                           'fwddiams_training_data.rds'))$P),
    plot_id = readRDS(here::here('02-data',
                           '05-for_analysis',
                           'fwddiams_training_data.rds'))$plot_id,
    transect_id = readRDS(here::here('02-data',
                           '05-for_analysis',
                           'fwddiams_training_data.rds'))$transect_id,
    lb = readRDS(here::here('02-data',
                           '05-for_analysis',
                           'fwddiams_training_data.rds'))$lb,
    ub = readRDS(here::here('02-data',
                           '05-for_analysis',
                           'fwddiams_training_data.rds'))$ub
  )

set.seed(110819)
fwddiams_prior = 
  data.frame(
    intercept = rnorm(n = 4000, mean = 0, sd = 1),
    SD_plot = truncnorm::rtruncnorm(n = 4000, mean = 0, sd = 0.25, a = 0),
    SD_trans = truncnorm::rtruncnorm(n = 4000, mean = 0, sd = 0.25, a = 0),
    phi = truncnorm::rtruncnorm(n = 4000, mean = 0, sd = 5, a = 0),
    draw = 1:4000
  )

fwddiams_sim = 
  do.call('rbind',
          lapply(X = 1:100,
                 FUN = 
                   function(draw){
                     print(paste0('Draw ',draw))
                     intercept = 
                       fwddiams_prior %>%
                       select(contains('intercept')) %>%
                       slice(draw) %>%
                       as.data.frame() %>%
                       as.numeric()
                     
                     phi = 
                       fwddiams_prior %>%
                       select(contains('phi')) %>%
                       slice(draw) %>%
                       as.data.frame() %>%
                       as.numeric()
                   
                     SD_plot = 
                       fwddiams_prior %>%
                       select(contains('SD_plot')) %>%
                       slice(draw) %>%
                       as.data.frame() %>%
                       as.numeric()
                     
                     effect_plot = 
                       rnorm(n = fwddiams_data$P,
                             mean = 0,
                             sd = SD_plot)
                   
                     SD_trans = 
                       fwddiams_prior %>%
                       select(contains('SD_trans')) %>%
                       slice(draw) %>%
                       as.data.frame() %>%
                       as.numeric()
                     effect_trans =
                       rnorm(n = fwddiams_data$TS,
                             mean = 0,
                             sd = SD_trans)
                     N = fwddiams_data$N
                     group_id = fwddiams_data$group_id
                     plot_id = fwddiams_data$plot_id
                     transect_id = fwddiams_data$transect_id
                   
                     logmu = 
                       sapply(X = 1:fwddiams_data$N,
                              FUN = function(i){
                                intercept[group_id[plot_id[transect_id[i]]]]+
                                  effect_plot[plot_id[transect_id[i]]]+
                                  effect_trans[transect_id[i]]
                              })
                    
                     shape = 
                       sapply(X = 1:N,
                              FUN = function(i){
                                phi[group_id[plot_id[transect_id[i]]]]+2
                              })
                   
                     rate = 
                       sapply(X = 1:N,
                              FUN = function(i){
                                exp(logmu[i])*(1+phi[group_id[plot_id[transect_id[i]]]])
                              })
                   
                     Ysim = 
                       sapply(X = 1:N,
                              FUN = function(i){
                                y = 0
                                while (y<=fwddiams_data$lb |
                                       y>=fwddiams_data$ub){
                                  y = 
                                    invgamma::rinvgamma(n = 1,
                                                        shape = shape[i],
                                                        rate = rate[i])
                                }
                                return(y)
                              })
                   
                     results = 
                       data.frame(transect_id = transect_id,
                                  plot_id = plot_id[transect_id],
                                  group_id = group_id[plot_id[transect_id]],
                                  Ysim = Ysim,
                                  draw = draw) %>%
                       rowid_to_column('obs_id')
                         
                     
                     return(results)
                   }))

fwddiams_sim %>%
  left_join(fwddiams_prior) %>%
  ggplot(aes(x = Ysim))+
  geom_histogram()+
  theme_minimal()


#### cwd tallies ###############################################################


cwdtallies_data = 
  list(
    N = readRDS(here::here('02-data',
                           '05-for_analysis',
                           'cwdtallies_training_data.rds'))$N,
    G = 1,
    P =   readRDS(here::here('02-data', 
                     '05-for_analysis',
                     'cwdtallies_training_data.rds'))$P,
    TR =   readRDS(here::here('02-data', 
                     '05-for_analysis',
                     'cwdtallies_training_data.rds'))$TR,
    L =   readRDS(here::here('02-data', 
                     '05-for_analysis',
                     'cwdtallies_training_data.rds'))$L,
    coords =   readRDS(here::here('02-data', 
                     '05-for_analysis',
                     'cwdtallies_training_data.rds'))$coords,
    location_id =   readRDS(here::here('02-data', 
                     '05-for_analysis',
                     'cwdtallies_training_data.rds'))$location_id,
    transect_id =   readRDS(here::here('02-data', 
                     '05-for_analysis',
                     'cwdtallies_training_data.rds'))$transect_id,
    plot_id =   readRDS(here::here('02-data', 
                     '05-for_analysis',
                     'cwdtallies_training_data.rds'))$plot_id,
    group_id = rep(1, times =   readRDS(here::here('02-data', 
                     '05-for_analysis',
                     'cwdtallies_training_data.rds'))$P)
  )

set.seed(110819)
cwdtallies_prior = 
  data.frame(
    draw = 1:4000,
    beta = rnorm(n = 4000, mean = 0, sd = 0.5),
    sigmaTrans = truncnorm::rtruncnorm(n = 4000, mean = 0, sd = 0.5, a = 0),
    sigmaPlot = truncnorm::rtruncnorm(n = 4000, mean = 0, sd = 0.5, a = 0),
    alpha = truncnorm::rtruncnorm(n = 4000, mean = 0, sd = 0.5, a = 0),
    rho = invgamma::rinvgamma(n = 4000, shape = 5, rate = 40)
  )

cwdtallies_sim = 
  do.call('bind_rows',
          lapply(X = 1:100,
                 FUN = function(draw){
                   
                   # extract parameters for this draw
                   beta = 
                     cwdtallies_prior %>%
                     select(contains('beta')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   alpha = 
                     cwdtallies_prior %>%
                     select(contains('alpha')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   rho = 
                     cwdtallies_prior %>%
                     select(contains('rho')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   sigmaPlot = 
                     cwdtallies_prior %>%
                     select(contains('sigmaPlot')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   plotEffects = 
                     rnorm(n = cwdtallies_data$P,
                           mean = 0,
                           sd = sigmaPlot)
                   
                   sigmaTrans = 
                     cwdtallies_prior %>%
                     select(contains('sigmaTrans')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   transectEffects = 
                     rnorm(n = cwdtallies_data$TR,
                           mean = 0, 
                           sd = sigmaTrans)
                   
                   zGP = 
                     matrix(nrow = cwdtallies_data$L,
                            ncol = cwdtallies_data$P,
                            byrow = FALSE,
                            data = 
                              rnorm(n = cwdtallies_data$L*cwdtallies_data$P,
                                    mean = 0,
                                    sd = 1))
                   
                   
                   dists = 
                     cwdtallies_data$coords %>%
                     dist() %>%
                     as.matrix()
                   
                   # construct GP effects
                   K = 
                     array(dim = 
                             c(cwdtallies_data$L,
                               cwdtallies_data$L,
                               cwdtallies_data$G),
                           dimnames = 
                             list('l' = 1:cwdtallies_data$L,
                                  'lprime' = 1:cwdtallies_data$L,
                                  'g' = 1:cwdtallies_data$G),
                           data = 
                             sapply(X = 1:cwdtallies_data$G,
                                    FUN = function(g){
                                      sapply(X = 1:cwdtallies_data$L,
                                             FUN = function(l){
                                               sapply(X = 1:cwdtallies_data$L,
                                                      FUN = function(lprime){
                                                        ((alpha[g]**2)*
                                                          exp(-(dists[lprime,l]**2)/
                                                                (2*rho[g]**2)))
                                                      })
                                             })
                                    }))
                   
                   for (g in 1:cwdtallies_data$G){
                     K[,,g] = K[,,g]+diag(nrow = nrow(K[,,g]),
                                          ncol = ncol(K[,,g]),
                                          x = 1e-9)
                   }
                   
                   KL = 
                     array(dim = 
                             c(cwdtallies_data$L,
                               cwdtallies_data$L,
                               cwdtallies_data$G),
                           dimnames = 
                             list('l' = 1:cwdtallies_data$L,
                                  'lprime' = 1:cwdtallies_data$L,
                                  'g' = 1:cwdtallies_data$G),
                           data = 
                             sapply(X = 1:cwdtallies_data$G,
                                    FUN = function(g){
                                      t(chol(K[,,g]))
                                    }))
                 
                   # realized random effects
                   GP = 
                     matrix(nrow = cwdtallies_data$L,
                            ncol = cwdtallies_data$P,
                            data = 
                              sapply(X = 1:cwdtallies_data$P,
                                     FUN = function(p){
                                       KL[,,cwdtallies_data$group_id[p]] %*%
                                         zGP[,p]
                                     }))
                   
                   # linear predictors
                   logmu = 
                     sapply(X = 1:cwdtallies_data$N,
                            FUN = function(i){
                              
                               beta[
                                 cwdtallies_data$group_id[
                                   cwdtallies_data$plot_id[
                                     cwdtallies_data$transect_id[i]]]] +
                                plotEffects[
                                  cwdtallies_data$plot_id[
                                    cwdtallies_data$transect_id[i]]]+
                                transectEffects[
                                  cwdtallies_data$transect_id[i]]+
                                as.numeric(GP[cwdtallies_data$location_id[i],
                                              cwdtallies_data$plot_id][
                                                cwdtallies_data$transect_id[i]
                                              ])
                              
                            })
                    
                   # realizations
                   Y = 
                     sapply(X = 1:cwdtallies_data$N,
                            FUN = function(i){
                              rpois(n = 1,
                                    lambda = exp(logmu[i]))})
                   
                   
                   results = 
                     data.frame(Y_sim = Y,
                                group_id = 
                                  cwdtallies_data$group_id[
                                    cwdtallies_data$plot_id[
                                      cwdtallies_data$transect_id]],
                                plot_id = 
                                  cwdtallies_data$plot_id[
                                    cwdtallies_data$transect_id],
                                transect_id = 
                                  cwdtallies_data$transect_id,
                                location_id = cwdtallies_data$location_id,
                                logmu = logmu,
                                draw = draw) %>%
                     as_tibble() %>%
                     rowid_to_column('obs_id')
                   
                   return(results)
                 }))

cwdtallies_sim %>%
  ggplot(aes(x = Y_sim))+
  geom_histogram()+
  theme_minimal()+
  scale_x_log10()

#### cwddiams ##################################################################

cwddiams_data = 
  list(
    N = readRDS(here::here('02-data',
                     '05-for_analysis',
                     'cwddiams_training_data.rds'))$N,
    P =   readRDS(here::here('02-data',
                     '05-for_analysis',
                     'cwddiams_training_data.rds'))$P,
    G = 1,
    group_id = rep(1, 
                   times = readRDS(here::here('02-data',
                                              '05-for_analysis',
                                              'cwddiams_training_data.rds'))$P),
    plot_id = readRDS(here::here('02-data',
                                 '05-for_analysis',
                                 'cwddiams_training_data.rds'))$plot_id,
    lb =   readRDS(here::here('02-data',
                     '05-for_analysis',
                     'cwddiams_training_data.rds'))$lb,
    ub =   readRDS(here::here('02-data',
                     '05-for_analysis',
                     'cwddiams_training_data.rds'))$ub)

set.seed(110819)
cwddiams_prior = 
  data.frame(
    mu = runif(n = 4000, min = cwddiams_data$lb, max = cwddiams_data$ub),
    phi = truncnorm::rtruncnorm(n= 4000, mean = 0, sd = 10, a = 0)
  )

cwddiams_sim = 
  do.call('rbind',
          lapply(X = 1:100,
                 FUN = 
                   function(draw){
                     print(paste0('Working on draw ', draw))
                     mu = 
                       cwddiams_prior %>%
                       select(contains('mu')) %>%
                       slice(draw) %>%
                       as.data.frame() %>%
                       as.numeric()
                     
                     phi = 
                       cwddiams_prior %>%
                       select(contains('phi')) %>%
                       slice(draw) %>%
                       as.data.frame() %>%
                       as.numeric()
                   
                     N = cwddiams_data$N
                     group_id = cwddiams_data$group_id
                     plot_id = cwddiams_data$plot_id
                   
                     shape = 
                       sapply(X = 1:N,
                              FUN = function(i){
                                phi[group_id[plot_id[i]]]+2
                              })
                   
                     rate = 
                       sapply(X = 1:N,
                              FUN = function(i){
                                mu[group_id[plot_id[i]]]*(1+phi[group_id[plot_id[i]]])
                              })
                   
                     Ysim = 
                       sapply(X = 1:N,
                              FUN = function(i){
                                y = 0
                                while (y<=cwddiams_data$lb |
                                       y>=cwddiams_data$ub){
                                  
                                  y = 
                                    invgamma::rinvgamma(n = 1,
                                                        shape = shape[i],
                                                        rate = rate[i])
                                }
                                return(y)
                              })
                   
                     results = 
                       data.frame(plot_id = plot_id,
                                  group_id = group_id[plot_id],
                                  Ysim = Ysim,
                                  draw = draw) %>%
                       rowid_to_column('obs_id')
                         
                     
                     return(results)
                   }))

cwddiams_sim %>%
  ggplot(aes(x = Ysim))+
  geom_histogram()+
  theme_minimal()

#### veg p/a ###################################################################

# just to get the data structure, not actually looking at the real data yet
vegpa_data = 
  list(
    G = 1,
    P = 10,
    L = readRDS(here::here('02-data','05-for_analysis','vegpa_training_data.rds'))$L,
    coords = readRDS(here::here('02-data','05-for_analysis','vegpa_training_data.rds'))$coords

  )

vegpa_data$N = vegpa_data$P*vegpa_data$L
vegpa_data$location_id = rep(1:vegpa_data$L, times = vegpa_data$P)
vegpa_data$plot_id = rep(1:vegpa_data$P, each = vegpa_data$L)
vegpa_data$group_id = rep(1, times = vegpa_data$P)

set.seed(110819)
vegpa_prior = 
  data.frame(
    'beta' = rnorm(n = 4000, mean = 0, sd = 1.5),
    'alpha' = truncnorm::rtruncnorm(n = 4000, mean = 0, sd = 1.5, a = 0),
    'rho' = invgamma::rinvgamma(n = 4000, shape = 5, rate = 40),
    'sigmaPlot' = truncnorm::rtruncnorm(n = 4000, mean = 0, sd = 1.5, a = 0),
    'kappa' = truncnorm::rtruncnorm(n = 4000, mean = 0, sd = 1.5, a = 0),
    draw = 1:4000
  )

vegpa_sim = 
  do.call('bind_rows',
          lapply(X = 1:400,
                 FUN = function(draw){
                   
                   # extract parameters for this draw
                   beta = 
                     vegpa_prior %>%
                     select(contains('beta')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   alpha = 
                     vegpa_prior %>%
                     select(contains('alpha')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   rho = 
                     vegpa_prior %>%
                     select(contains('rho')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   sigmaPlot = 
                     vegpa_prior %>%
                     select(contains('sigmaPlot')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   plotEffects = 
                     rnorm(n = vegpa_data$P,
                           mean = 0,
                           sd = sigmaPlot)
                     
                   zGP = 
                     matrix(nrow = vegpa_data$L,
                            ncol = vegpa_data$P,
                            byrow = FALSE,
                            data = 
                              rnorm(n = 
                                      vegpa_data$L*
                                      vegpa_data$P,
                                    mean = 0,
                                    sd = 1))
                   
                   
                   dists = 
                     vegpa_data$coords %>%
                     dist() %>%
                     as.matrix()
                   
                   # construct GP effects
                   K = 
                     array(dim = 
                             c(vegpa_data$L,
                               vegpa_data$L,
                               vegpa_data$G),
                           dimnames = 
                             list('l' = 1:vegpa_data$L,
                                  'lprime' = 1:vegpa_data$L,
                                  'g' = 1:vegpa_data$G),
                           data = 
                             sapply(X = 1:vegpa_data$G,
                                    FUN = function(g){
                                      sapply(X = 1:vegpa_data$L,
                                             FUN = function(l){
                                               sapply(X = 1:vegpa_data$L,
                                                      FUN = function(lprime){
                                                        ((alpha[g]**2)*
                                                          exp(-(dists[lprime,l]**2)/
                                                                (2*rho[g]**2)))
                                                      })
                                             })
                                    }))
                   
                   for (g in 1:vegpa_data$G){
                     K[,,g] = K[,,g]+diag(nrow = nrow(K[,,g]),
                                          ncol = ncol(K[,,g]),
                                          x = 1e-9)
                   }
                   
                   KL = 
                     array(dim = 
                             c(vegpa_data$L,
                               vegpa_data$L,
                               vegpa_data$G),
                           dimnames = 
                             list('l' = 1:vegpa_data$L,
                                  'lprime' = 1:vegpa_data$L,
                                  'g' = 1:vegpa_data$G),
                           data = 
                             sapply(X = 1:vegpa_data$G,
                                    FUN = function(g){
                                      t(chol(K[,,g]))
                                    }))
                 
                   # realized random effects
                   GP = 
                     matrix(nrow = vegpa_data$L,
                            ncol = vegpa_data$P,
                            data = 
                              sapply(X = 1:vegpa_data$P,
                                     FUN = function(p){
                                       KL[,,vegpa_data$group_id[p]] %*%
                                         zGP[,p]
                                     }))
                   
                   # linear predictors
                   logiteta = 
                     sapply(X = 1:vegpa_data$N,
                            FUN = function(i){
                              
                               beta[vegpa_data$group_id[vegpa_data$plot_id[i]]] +
                                plotEffects[vegpa_data$plot_id[i]]+
                                as.numeric(GP[vegpa_data$location_id[i],
                                              vegpa_data$plot_id][i])
                              
                            })
                    
                   # realizations
                   Y = 
                     sapply(X = 1:vegpa_data$N,
                            FUN = function(i){
                              rbinom(n = 1,
                                     size = 1,
                                     prob = boot::inv.logit(logiteta[i]))
                            })
                   
                   
                   results = 
                     data.frame(Y_sim = Y,
                                group_id = 
                                  vegpa_data$group_id[vegpa_data$plot_id],
                                plot_id = vegpa_data$plot_id,
                                location_id = vegpa_data$location_id,
                                logiteta = logiteta,
                                draw = draw) %>%
                     as_tibble() %>%
                     rowid_to_column('obs_id')
                   
                   return(results)
                 }))

# plot mean as a function of beta
vegpa_sim %>%
  mutate(p = boot::inv.logit(logiteta)) %>%
  group_by(draw, group_id, plot_id) %>%
  summarise(Y_sim = sum(Y_sim)/n()) %>%
  ungroup() %>%
  left_join(vegpa_prior) %>%
  ggplot(aes(x = Y_sim))+
  geom_histogram()+
  theme_minimal()

#### trees #####################################################################


# just to get the data structure, not actually looking at the real data yet
trees_data = 
  list(
    G = 1,
    P = 10,
    L = readRDS(here::here('02-data','05-for_analysis','trees_training_data.rds'))$L,
    coords = readRDS(here::here('02-data','05-for_analysis','trees_training_data.rds'))$coords
  )

trees_data$N = trees_data$P*trees_data$L
trees_data$location_id = rep(1:trees_data$L, times = trees_data$P)
trees_data$plot_id = rep(1:trees_data$P, each = trees_data$L)
trees_data$group_id = rep(1, times = trees_data$P)

trees_prior = 
  data.frame(
    'beta' = rnorm(n = 4000, mean = -1, sd = 0.5),
    'alpha' = truncnorm::rtruncnorm(n = 4000, mean = 0, sd = 0.5, a = 0),
    'rho' = invgamma::rinvgamma(n = 4000, shape = 5, rate = 40),
    'sigmaPlot' = truncnorm::rtruncnorm(n = 4000, mean = 0, sd = 0.5, a = 0),
    draw = 1:4000
  )

set.seed(110819)
trees_sim = 
  do.call('bind_rows',
          lapply(X = 1:100,
                 FUN = function(draw){
                   
                   # extract parameters for this draw
                   beta = 
                     trees_prior %>%
                     select(contains('beta')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   alpha = 
                     trees_prior %>%
                     select(contains('alpha')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   rho = 
                     trees_prior %>%
                     select(contains('rho')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   sigmaPlot = 
                     trees_prior %>%
                     select(contains('sigmaPlot')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   plotEffects = 
                     rnorm(n = trees_data$P,
                           mean = 0,
                           sd = sigmaPlot)
                     
                   zGP = 
                     matrix(nrow = trees_data$L,
                            ncol = trees_data$P,
                            byrow = FALSE,
                            data = 
                              rnorm(n = 
                                      trees_data$L*
                                      trees_data$P,
                                    mean = 0,
                                    sd = 1))

                   dists = 
                     trees_data$coords %>%
                     dist() %>%
                     as.matrix()
                   
                   # construct GP effects
                   K = 
                     array(dim = 
                             c(trees_data$L,
                               trees_data$L,
                               trees_data$G),
                           dimnames = 
                             list('l' = 1:trees_data$L,
                                  'lprime' = 1:trees_data$L,
                                  'g' = 1:trees_data$G),
                           data = 
                             sapply(X = 1:trees_data$G,
                                    FUN = function(g){
                                      sapply(X = 1:trees_data$L,
                                             FUN = function(l){
                                               sapply(X = 1:trees_data$L,
                                                      FUN = function(lprime){
                                                        ((alpha[g]**2)*
                                                          exp(-(dists[lprime,l]**2)/
                                                                (2*rho[g]**2)))
                                                      })
                                             })
                                    }))
                   
                   for (g in 1:trees_data$G){
                     K[,,g] = K[,,g]+diag(nrow = nrow(K[,,g]),
                                          ncol = ncol(K[,,g]),
                                          x = 1e-9)
                   }
                   
                   KL = 
                     array(dim = 
                             c(trees_data$L,
                               trees_data$L,
                               trees_data$G),
                           dimnames = 
                             list('l' = 1:trees_data$L,
                                  'lprime' = 1:trees_data$L,
                                  'g' = 1:trees_data$G),
                           data = 
                             sapply(X = 1:trees_data$G,
                                    FUN = function(g){
                                      t(chol(K[,,g]))
                                    }))
                 
                   # realized random effects
                   GP = 
                     matrix(nrow = trees_data$L,
                            ncol = trees_data$P,
                            data = 
                              sapply(X = 1:trees_data$P,
                                     FUN = function(p){
                                       KL[,,trees_data$group_id[p]] %*%
                                         zGP[,p]
                                     }))
                   
                   # linear predictors
                   logmu = 
                     sapply(X = 1:trees_data$N,
                            FUN = function(i){
                              
                               beta[trees_data$group_id[trees_data$plot_id[i]]] +
                                plotEffects[trees_data$plot_id[i]]+
                                as.numeric(GP[trees_data$location_id[i],
                                              trees_data$plot_id][i])
                              
                            })
                    
                   # realizations
                   Y = 
                     sapply(X = 1:trees_data$N,
                            FUN = function(i){
                              rpois(n = 1,
                                    lambda = exp(logmu[i]))})
                   
                   
                   results = 
                     data.frame(Y_sim = Y,
                                group_id = 
                                  trees_data$group_id[trees_data$plot_id],
                                plot_id = trees_data$plot_id,
                                location_id = trees_data$location_id,
                                logmu = logmu,
                                draw = draw) %>%
                     as_tibble() %>%
                     rowid_to_column('obs_id')
                   
                   return(results)
                 }))

# plot mean as a function of beta
trees_sim %>%
  left_join(trees_prior) %>%
  ggplot(aes(x = Y_sim))+
  geom_histogram()+
  theme_minimal()


