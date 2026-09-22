# todos:
# remove or change main titles
# 

#### setup #####################################################################

library(here)
library(tidyverse)
library(posterior)


#### litter ####################################################################



litter_validation_data = 
  readRDS(here::here('02-data', 
                     '05-for_analysis',
                     'litter_validation_data.rds'))

litter_posterior = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'litter_fit.rds'))$draws() %>%
  as_draws_df()

set.seed(110819)
litter_sim = 
  do.call('bind_rows',
          lapply(X = 1:4000,
                 FUN = function(draw){
                   
                   print(paste0('Working on draw ', draw))
                   
                   # extract parameters for this draw
                   beta = 
                     litter_posterior %>%
                     select(contains('beta')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   alpha = 
                     litter_posterior %>%
                     select(contains('alpha')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   rho = 
                     litter_posterior %>%
                     select(contains('rho')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   plotEffects = 
                     litter_posterior %>%
                     select(contains('plot_effect')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   zGP = 
                     matrix(nrow = litter_validation_data$L,
                            ncol = litter_validation_data$P,
                            byrow = FALSE,
                            data = 
                              litter_posterior %>%
                              select(contains('zGP')) %>%
                              slice(draw) %>%
                              as.data.frame() %>%
                              as.numeric())
                   
                   kappa = 
                     litter_posterior %>%
                     select(contains('kappa')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   dists = 
                     litter_validation_data$coords %>%
                     dist() %>%
                     as.matrix()
                   
                   # construct GP effects
                   K = 
                     array(dim = 
                             c(litter_validation_data$L,
                               litter_validation_data$L,
                               litter_validation_data$G),
                           dimnames = 
                             list('l' = 1:litter_validation_data$L,
                                  'lprime' = 1:litter_validation_data$L,
                                  'g' = 1:litter_validation_data$G),
                           data = 
                             sapply(X = 1:litter_validation_data$G,
                                    FUN = function(g){
                                      sapply(X = 1:litter_validation_data$L,
                                             FUN = function(l){
                                               sapply(X = 1:litter_validation_data$L,
                                                      FUN = function(lprime){
                                                        ((alpha[g]**2)*
                                                           exp(-(dists[lprime,l]**2)/
                                                                 (2*rho[g]**2)))
                                                      })
                                             })
                                    }))
                   
                   for (g in 1:litter_validation_data$G){
                     K[,,g] = K[,,g]+diag(nrow = nrow(K[,,g]),
                                          ncol = ncol(K[,,g]),
                                          x = 1e-9)
                   }
                   
                   KL = 
                     array(dim = 
                             c(litter_validation_data$L,
                               litter_validation_data$L,
                               litter_validation_data$G),
                           dimnames = 
                             list('l' = 1:litter_validation_data$L,
                                  'lprime' = 1:litter_validation_data$L,
                                  'g' = 1:litter_validation_data$G),
                           data = 
                             sapply(X = 1:litter_validation_data$G,
                                    FUN = function(g){
                                      t(chol(K[,,g]))
                                    }))
                   
                   # realized random effects
                   GP = 
                     matrix(nrow = litter_validation_data$L,
                            ncol = litter_validation_data$P,
                            data = 
                              sapply(X = 1:litter_validation_data$P,
                                     FUN = function(p){
                                       KL[,,litter_validation_data$group_id[p]] %*%
                                         zGP[,p]
                                     }))
                   
                   # linear predictors
                   logmu = 
                     sapply(X = 1:litter_validation_data$N,
                            FUN = function(i){
                              
                              beta[litter_validation_data$group_id[litter_validation_data$plot_id[i]]] +
                                plotEffects[litter_validation_data$plot_id[i]]+
                                as.numeric(GP[litter_validation_data$location_id[i],
                                              litter_validation_data$plot_id][i])
                              
                            })
                   
                   # realizations
                   Y = 
                     sapply(X = 1:litter_validation_data$N,
                            FUN = function(i){
                              rnbinom(n = 1,
                                      mu = exp(logmu[i]),
                                      size = 
                                        kappa[litter_validation_data$group_id[litter_validation_data$plot_id[i]]])
                            })
                   
                   
                   results = 
                     data.frame(Y_obs = litter_validation_data$Y,
                                Y_sim = Y,
                                group_id = 
                                  litter_validation_data$group_id[litter_validation_data$plot_id],
                                plot_id = litter_validation_data$plot_id,
                                location_id = litter_validation_data$location_id,
                                logmu = logmu,
                                draw = draw) %>%
                     as_tibble() %>%
                     rowid_to_column('obs_id')
                   
                   return(results)
                 }))

saveRDS(litter_sim,
        here::here('02-data', '06-results', 'posterior_predictions',
                   'litter_sim.rds'))

litter_sim = 
  readRDS(here::here('02-data', '06-results', 'posterior_predictions',
                     'litter_sim.rds'))

litter_sim %>%
  group_by(obs_id, Y_obs, group_id) %>%
  summarise(Y_sim = mean(Y_sim)) %>%
  ungroup() %>%
  mutate(ae = abs(Y_obs-Y_sim)) %>%
  pull(ae) %>%
  mean() #1.71

1.71 / mean(litter_validation_data$Y) # 0.356

litter_plot_1 = 
  litter_sim %>%
  group_by(obs_id, Y_obs, group_id) %>%
  summarise(Y_sim = mean(Y_sim)) %>%
  ungroup() %>%
  ggplot(aes(x = Y_sim, y = Y_obs))+
  geom_abline(intercept=0,slope=1, color = 'red')+
  geom_smooth(method = 'lm')+
  geom_point(size = 0)+
  theme_minimal()+
  coord_fixed()+
  facet_wrap(~group_id)+
  labs(x = 'Mean simulated litter depth (cm)',
       y = 'Observed (validation data) litter depth (cm)')

litter_plot_1

ggsave(litter_plot_1,
       filename = here::here('04-communication','figures', 'manuscript',
                             'litter_prediction_1.png'),
       height = 4, width = 6.5, units = 'in')

litter_plot_2 = 
  litter_sim %>%
  group_by(obs_id, Y_obs, group_id) %>%
  summarise() %>%
  ungroup() %>%
  mutate(source = 'obs') %>%
  group_by(Y_obs, group_id, source) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  left_join(
    litter_sim %>%
      group_by(obs_id, Y_obs, group_id) %>%
      summarise() %>%
      ungroup() %>%
      mutate(source = 'obs') %>%
      group_by(group_id, source) %>%
      summarise(total = n()) %>%
      ungroup() 
  ) %>%
  mutate(proportion = count / total) %>%
  select(group_id, Y = Y_obs, proportion, source) %>%
  bind_rows(
    litter_sim %>%
      group_by(Y_sim, group_id) %>%
      summarise(count = n()) %>%
      ungroup() %>%
      mutate(source = 'sim') %>%
      left_join(
        litter_sim %>%
          group_by(group_id) %>%
          summarise(total = n()) %>%
          ungroup() %>%
          mutate(source = 'sim')
      ) %>%
      mutate(proportion = count / total) %>%
      select(group_id, Y = Y_sim, proportion, source)
  ) %>%
  ggplot()+
  geom_col(aes(x = Y, y = proportion, fill = source),
           position = position_dodge())+
  facet_grid(group_id~., scales = 'free')+
  theme_minimal()+
  labs(title = 'litter prediction')

litter_plot_2

ggsave(litter_plot_2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'litter_prediction_2.png'),
       height = 6.5, width = 6.5, units = 'in')

litter_plot_3 = 
  litter_sim %>%
  group_by(group_id, obs_id, Y_obs) %>%
  summarise(logmu.mean = mean(logmu),
            Y_sim.mean = mean(Y_sim),
            Y_sim.05 = quantile(Y_sim, 0.05),
            Y_sim.95 = quantile(Y_sim, 0.95)) %>%
  ungroup() %>%
  ggplot()+
  geom_errorbar(aes(x = logmu.mean, ymin = Y_sim.05, ymax = Y_sim.95),
                alpha = 0.75)+
  geom_point(aes(x = logmu.mean, y = Y_obs), color = 'red')+
  theme_minimal()+
  labs(title = 'litter prediction')+
  facet_grid(group_id~., scales = 'free')

litter_plot_3

ggsave(litter_plot_3,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'litter_prediction_3.png'),
       height = 6.5, width = 6.5, units = 'in')


litter_plot_4 = 
  ggplot()+
  geom_errorbar(
    data = 
      litter_sim %>%
      group_by(group_id, Y_sim, draw) %>%
      summarise(N_obs = n()) %>%
      ungroup() %>%
      right_join(litter_sim %>%
                   expand(nesting(group_id, draw),
                          Y_sim = 0:max(.$Y_sim))) %>%
      mutate(N_obs = ifelse(is.na(N_obs), 0, N_obs)) %>%
      group_by(group_id, Y_sim) %>%
      summarise(N_obs.025 = quantile(N_obs, 0.025),
                N_obs.975 = quantile(N_obs, 0.975)) %>%
      ungroup(),
    aes(x = Y_sim, ymin = N_obs.025, ymax = N_obs.975, color = 'Simulated')
  )+
  facet_grid(group_id~., scales = 'free')+
  theme_minimal()+
  labs(x = 'Litter depth (cm)', y = 'N observations')+
  geom_col(
    data = 
      litter_sim %>%
      group_by(group_id, obs_id, Y_obs) %>%
      summarise() %>%
      ungroup() %>%
      group_by(group_id, Y_obs) %>%
      summarise(count = n()) %>%
      ungroup(),
    aes(x = Y_obs, y = count, fill = 'Observed (validation data)'),
    width = 0.5,
    alpha = 0.75
  )+
  scale_fill_manual(values = c('Observed (validation data)' = 'red'))+
  scale_color_manual(values = c('Simulated' = 'blue'))+
  theme(legend.position = 'bottom',
        legend.title = element_blank())


litter_plot_4

ggsave(litter_plot_4,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'litter_prediction_4.png'),
       height = 6.5, width = 6.5, units = 'in')

head(litter_sim)

litter_plot_5 = 
  ggplot()+
  stat_ecdf(
    data = litter_sim,
    aes(x = Y_sim, group = draw, color = 'Simulated'),
    geom = 'step', 
    pad = FALSE,
    #color = 'blue',
    lwd = 1,
    alpha = 0.1)+
  stat_ecdf(
    data = litter_sim %>% group_by(group_id, Y_obs, obs_id) %>% summarise() %>% ungroup(),
    aes(x = Y_obs, color = 'Observed (validation data)'),
    geom = 'step',
    pad = 'false',
    #color = 'red',
    lwd = 1
  )+
  theme_minimal()+
  scale_color_manual(values = c('Observed (validation data)' = 'red', 'Simulated' = 'blue'))+
  facet_grid(group_id~.)+
  labs(x = 'Litter depth (cm)', y = 'Cumulative density')+
  theme(legend.title = element_blank(),
        legend.position = 'bottom')

litter_plot_5

ggsave(litter_plot_5,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'litter_prediction_5.png'),
       height = 6.5, width = 6.5, units = 'in')


rm(litter_sim)

#### duff ######################################################################

duff_validation_data = 
  readRDS(here::here('02-data', 
                     '05-for_analysis',
                     'duff_validation_data.rds'))

duff_posterior = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'duff_fit.rds'))$draws() %>%
  as_draws_df()

set.seed(110819)
duff_sim = 
  do.call('bind_rows',
          lapply(X = 1:4000,
                 FUN = function(draw){
                   
                   print(paste0('Working on draw ', draw))
                   
                   # extract parameters for this draw
                   beta = 
                     duff_posterior %>%
                     select(contains('beta')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   alpha = 
                     duff_posterior %>%
                     select(contains('alpha')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   rho = 
                     duff_posterior %>%
                     select(contains('rho')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   plotEffects = 
                     duff_posterior %>%
                     select(contains('plot_effect')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   zGP = 
                     matrix(nrow = duff_validation_data$L,
                            ncol = duff_validation_data$P,
                            byrow = FALSE,
                            data = 
                              duff_posterior %>%
                              select(contains('zGP')) %>%
                              slice(draw) %>%
                              as.data.frame() %>%
                              as.numeric())
                   
                   kappa = 
                     duff_posterior %>%
                     select(contains('kappa')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   dists = 
                     duff_validation_data$coords %>%
                     dist() %>%
                     as.matrix()
                   
                   # construct GP effects
                   K = 
                     array(dim = 
                             c(duff_validation_data$L,
                               duff_validation_data$L,
                               duff_validation_data$G),
                           dimnames = 
                             list('l' = 1:duff_validation_data$L,
                                  'lprime' = 1:duff_validation_data$L,
                                  'g' = 1:duff_validation_data$G),
                           data = 
                             sapply(X = 1:duff_validation_data$G,
                                    FUN = function(g){
                                      sapply(X = 1:duff_validation_data$L,
                                             FUN = function(l){
                                               sapply(X = 1:duff_validation_data$L,
                                                      FUN = function(lprime){
                                                        ((alpha[g]**2)*
                                                           exp(-(dists[lprime,l]**2)/
                                                                 (2*rho[g]**2)))
                                                      })
                                             })
                                    }))
                   
                   for (g in 1:duff_validation_data$G){
                     K[,,g] = K[,,g]+diag(nrow = nrow(K[,,g]),
                                          ncol = ncol(K[,,g]),
                                          x = 1e-9)
                   }
                   
                   KL = 
                     array(dim = 
                             c(duff_validation_data$L,
                               duff_validation_data$L,
                               duff_validation_data$G),
                           dimnames = 
                             list('l' = 1:duff_validation_data$L,
                                  'lprime' = 1:duff_validation_data$L,
                                  'g' = 1:duff_validation_data$G),
                           data = 
                             sapply(X = 1:duff_validation_data$G,
                                    FUN = function(g){
                                      t(chol(K[,,g]))
                                    }))
                   
                   # realized random effects
                   GP = 
                     matrix(nrow = duff_validation_data$L,
                            ncol = duff_validation_data$P,
                            data = 
                              sapply(X = 1:duff_validation_data$P,
                                     FUN = function(p){
                                       KL[,,duff_validation_data$group_id[p]] %*%
                                         zGP[,p]
                                     }))
                   
                   # linear predictors
                   logmu = 
                     sapply(X = 1:duff_validation_data$N,
                            FUN = function(i){
                              
                              beta[duff_validation_data$group_id[duff_validation_data$plot_id[i]]] +
                                plotEffects[duff_validation_data$plot_id[i]]+
                                as.numeric(GP[duff_validation_data$location_id[i],
                                              duff_validation_data$plot_id][i])
                              
                            })
                   
                   # realizations
                   Y = 
                     sapply(X = 1:duff_validation_data$N,
                            FUN = function(i){
                              rnbinom(n = 1,
                                      mu = exp(logmu[i]),
                                      size = 
                                        kappa[duff_validation_data$group_id[duff_validation_data$plot_id[i]]])
                            })
                   
                   
                   results = 
                     data.frame(Y_obs = duff_validation_data$Y,
                                Y_sim = Y,
                                group_id = 
                                  duff_validation_data$group_id[duff_validation_data$plot_id],
                                plot_id = duff_validation_data$plot_id,
                                location_id = duff_validation_data$location_id,
                                logmu = logmu,
                                draw = draw) %>%
                     as_tibble() %>%
                     rowid_to_column('obs_id')
                   
                   return(results)
                 }))

saveRDS(duff_sim,
        here::here('02-data',
                   '06-results',
                   'posterior_predictions',
                   'duff_sim.rds'))

duff_sim = readRDS(here::here('02-data',
                              '06-results',
                              'posterior_predictions',
                              'duff_sim.rds'))

duff_sim %>%
  group_by(obs_id, Y_obs, group_id) %>%
  summarise(Y_sim = mean(Y_sim)) %>%
  ungroup() %>%
  mutate(ae = abs(Y_obs-Y_sim)) %>%
  pull(ae) %>%
  mean() # 1.52

1.52 / mean(duff_validation_data$Y)

duff_plot_1 = 
  duff_sim %>%
  group_by(obs_id, Y_obs, group_id) %>%
  summarise(Y_sim = mean(Y_sim)) %>%
  ungroup() %>%
  ggplot(aes(x = Y_sim, y = Y_obs))+
  geom_abline(intercept=0,slope=1, color = 'red')+
  geom_smooth(method = 'lm')+
  geom_point(size = 0)+
  theme_minimal()+
  coord_fixed()+
  facet_wrap(~group_id)+
  labs(x = 'Mean simulated duff depth (cm)',
       y = 'Observed (validation data) duff depth (cm)')

duff_plot_1

ggsave(duff_plot_1,
       filename = here::here('04-communication','figures', 'manuscript',
                             'duff_prediction_1.png'),
       height = 4, width = 6.5, units = 'in')

duff_plot_2 = 
  duff_sim %>%
  group_by(obs_id, Y_obs, group_id) %>%
  summarise() %>%
  ungroup() %>%
  mutate(source = 'obs') %>%
  group_by(Y_obs, group_id, source) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  left_join(
    duff_sim %>%
      group_by(obs_id, Y_obs, group_id) %>%
      summarise() %>%
      ungroup() %>%
      mutate(source = 'obs') %>%
      group_by(group_id, source) %>%
      summarise(total = n()) %>%
      ungroup() 
  ) %>%
  mutate(proportion = count / total) %>%
  select(group_id, Y = Y_obs, proportion, source) %>%
  bind_rows(
    duff_sim %>%
      group_by(Y_sim, group_id) %>%
      summarise(count = n()) %>%
      ungroup() %>%
      mutate(source = 'sim') %>%
      left_join(
        duff_sim %>%
          group_by(group_id) %>%
          summarise(total = n()) %>%
          ungroup() %>%
          mutate(source = 'sim')
      ) %>%
      mutate(proportion = count / total) %>%
      select(group_id, Y = Y_sim, proportion, source)
  ) %>%
  ggplot()+
  geom_col(aes(x = Y, y = proportion, fill = source),
           position = position_dodge())+
  facet_grid(group_id~., scales = 'free')+
  theme_minimal()+
  labs(title = 'duff prediction')+
  scale_x_continuous(limits = c(0, 50))

duff_plot_2

ggsave(duff_plot_2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'duff_prediction_2.png'),
       height = 6.5, width = 6.5, units = 'in')

duff_plot_3 = 
  duff_sim %>%
  group_by(group_id, obs_id, Y_obs) %>%
  summarise(logmu.mean = mean(logmu),
            Y_sim.mean = mean(Y_sim),
            Y_sim.05 = quantile(Y_sim, 0.05),
            Y_sim.95 = quantile(Y_sim, 0.95)) %>%
  ungroup() %>%
  ggplot()+
  geom_errorbar(aes(x = logmu.mean, ymin = Y_sim.05, ymax = Y_sim.95),
                alpha = 0.75)+
  geom_point(aes(x = logmu.mean, y = Y_obs), color = 'red')+
  theme_minimal()+
  labs(title = 'duff prediction')+
  facet_grid(group_id~., scales = 'free')

duff_plot_3

ggsave(duff_plot_3,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'duff_prediction_3.png'),
       height = 6.5, width = 6.5, units = 'in')


duff_plot_4 = 
  ggplot()+
  geom_errorbar(
    data = 
      duff_sim %>%
      group_by(group_id, Y_sim, draw) %>%
      summarise(N_obs = n()) %>%
      ungroup() %>%
      right_join(duff_sim %>%
                   expand(nesting(group_id, draw),
                          Y_sim = 0:max(.$Y_sim))) %>%
      mutate(N_obs = ifelse(is.na(N_obs), 0, N_obs)) %>%
      group_by(group_id, Y_sim) %>%
      summarise(N_obs.025 = quantile(N_obs, 0.025),
                N_obs.975 = quantile(N_obs, 0.975)) %>%
      ungroup(),
    aes(x = Y_sim, ymin = N_obs.025, ymax = N_obs.975, color = 'Simulated')
  )+
  facet_grid(group_id~., scales = 'free')+
  theme_minimal()+
  labs(x = 'Duff depth (cm)', y = 'N observations')+
  geom_col(
    data = 
      duff_sim %>%
      group_by(group_id, obs_id, Y_obs) %>%
      summarise() %>%
      ungroup() %>%
      group_by(group_id, Y_obs) %>%
      summarise(count = n()) %>%
      ungroup(),
    aes(x = Y_obs, y = count, fill = 'Observed (validation data)'),
    width = 0.5,
    alpha = 0.75
  )+
  scale_fill_manual(values = c('Observed (validation data)' = 'red'))+
  scale_color_manual(values = c('Simulated' = 'blue'))+
  theme(legend.position = 'bottom',
        legend.title = element_blank())+
  scale_x_continuous(limits = c(-1, 50))


duff_plot_4

ggsave(duff_plot_4,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'duff_prediction_4.png'),
       height = 6.5, width = 6.5, units = 'in')

duff_plot_5 = 
  ggplot()+
  stat_ecdf(
    data = duff_sim,
    aes(x = Y_sim, group = draw, color = 'Simulated'),
    geom = 'step', 
    pad = FALSE,
    #color = 'blue',
    lwd = 1,
    alpha = 0.01)+
  stat_ecdf(
    data = duff_sim %>% group_by(group_id, Y_obs, obs_id) %>% summarise() %>% ungroup(),
    aes(x = Y_obs, color = 'Observed (validation data)'),
    geom = 'step',
    pad = 'false',
    #color = 'red',
    lwd = 1
  )+
  theme_minimal()+
  scale_color_manual(values = c('Observed (validation data)' = 'red', 'Simulated' = 'blue'))+
  facet_grid(group_id~.)+
  labs(x = 'Duff depth (cm)', y = 'Cumulative density')+
  scale_x_continuous(limits = c(-1, 50))+
  theme(legend.position = 'bottom',
        legend.title = element_blank())

duff_plot_5

ggsave(duff_plot_5,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'duff_prediction_5.png'),
       height = 6.5, width = 6.5, units = 'in')

rm(duff_sim)

#### fwdtallies ######################################################################

fwdtallies_validation_data = 
  readRDS(here::here('02-data', 
                     '05-for_analysis',
                     'fwdtallies_validation_data.rds'))

fwdtallies_posterior = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'fwdtallies_fit.rds'))$draws() %>%
  as_draws_df()

set.seed(110819)
fwdtallies_sim = 
  do.call('bind_rows',
          lapply(X = 1:4000,
                 FUN = function(draw){
                   
                   # extract parameters for this draw
                   beta = 
                     fwdtallies_posterior %>%
                     select(contains('beta')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   alpha = 
                     fwdtallies_posterior %>%
                     select(contains('alpha')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   rho = 
                     fwdtallies_posterior %>%
                     select(contains('rho')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   plotEffects = 
                     fwdtallies_posterior %>%
                     select(contains('plot_effect')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   zGP = 
                     matrix(nrow = fwdtallies_validation_data$L,
                            ncol = fwdtallies_validation_data$P,
                            byrow = FALSE,
                            data = 
                              fwdtallies_posterior %>%
                              select(contains('zGP')) %>%
                              slice(draw) %>%
                              as.data.frame() %>%
                              as.numeric())
                   
                   tau = 
                     fwdtallies_posterior %>%
                     select(contains('tau')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   kappa = 
                     fwdtallies_posterior %>%
                     select(contains('kappa')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   dists = 
                     fwdtallies_validation_data$coords %>%
                     dist() %>%
                     as.matrix()
                   
                   # construct GP effects
                   K = 
                     array(dim = 
                             c(fwdtallies_validation_data$L,
                               fwdtallies_validation_data$L,
                               fwdtallies_validation_data$G),
                           dimnames = 
                             list('l' = 1:fwdtallies_validation_data$L,
                                  'lprime' = 1:fwdtallies_validation_data$L,
                                  'g' = 1:fwdtallies_validation_data$G),
                           data = 
                             sapply(X = 1:fwdtallies_validation_data$G,
                                    FUN = function(g){
                                      sapply(X = 1:fwdtallies_validation_data$L,
                                             FUN = function(l){
                                               sapply(X = 1:fwdtallies_validation_data$L,
                                                      FUN = function(lprime){
                                                        ((alpha[g]**2)*
                                                           exp(-(dists[lprime,l]**2)/
                                                                 (2*rho[g]**2)))
                                                      })
                                             })
                                    }))
                   
                   for (g in 1:fwdtallies_validation_data$G){
                     K[,,g] = K[,,g]+diag(nrow = nrow(K[,,g]),
                                          ncol = ncol(K[,,g]),
                                          x = tau[g])
                   }
                   
                   KL = 
                     array(dim = 
                             c(fwdtallies_validation_data$L,
                               fwdtallies_validation_data$L,
                               fwdtallies_validation_data$G),
                           dimnames = 
                             list('l' = 1:fwdtallies_validation_data$L,
                                  'lprime' = 1:fwdtallies_validation_data$L,
                                  'g' = 1:fwdtallies_validation_data$G),
                           data = 
                             sapply(X = 1:fwdtallies_validation_data$G,
                                    FUN = function(g){
                                      t(chol(K[,,g]))
                                    }))
                   
                   # realized random effects
                   GP = 
                     matrix(nrow = fwdtallies_validation_data$L,
                            ncol = fwdtallies_validation_data$P,
                            data = 
                              sapply(X = 1:fwdtallies_validation_data$P,
                                     FUN = function(p){
                                       KL[,,fwdtallies_validation_data$group_id[p]] %*%
                                         zGP[,p]
                                     }))
                   
                   # linear predictors
                   logmu = 
                     sapply(X = 1:fwdtallies_validation_data$N,
                            FUN = function(i){
                              
                              beta[fwdtallies_validation_data$group_id[fwdtallies_validation_data$plot_id[i]]] +
                                plotEffects[fwdtallies_validation_data$plot_id[i]]+
                                as.numeric(GP[fwdtallies_validation_data$location_id[i],
                                              fwdtallies_validation_data$plot_id][i])
                              
                            })
                   
                   # realizations
                   Y = 
                     sapply(X = 1:fwdtallies_validation_data$N,
                            FUN = function(i){
                              rnbinom(n = 1,
                                      mu = exp(logmu[i]),
                                      size = 
                                        kappa[fwdtallies_validation_data$group_id[fwdtallies_validation_data$plot_id[i]]])
                            })
                   
                   
                   results = 
                     data.frame(Y_obs = fwdtallies_validation_data$Y,
                                Y_sim = Y,
                                group_id = 
                                  fwdtallies_validation_data$group_id[fwdtallies_validation_data$plot_id],
                                plot_id = fwdtallies_validation_data$plot_id,
                                location_id = fwdtallies_validation_data$location_id,
                                logmu = logmu,
                                draw = draw) %>%
                     as_tibble() %>%
                     rowid_to_column('obs_id')
                   
                   return(results)
                 }))

saveRDS(fwdtallies_sim,
        here::here('02-data',
                   '06-results',
                   'posterior_predictions',
                   'fwdtallies_sim.rds'))

fwdtallies_sim = 
  readRDS(here::here('02-data',
                     '06-results',
                     'posterior_predictions',
                     'fwdtallies_sim.rds'))



fwdtallies_sim %>%
  group_by(obs_id, Y_obs, group_id) %>%
  summarise(Y_sim = mean(Y_sim)) %>%
  ungroup() %>%
  mutate(ae = abs(Y_obs-Y_sim)) %>%
  pull(ae) %>%
  mean() # 3.78

3.78/mean(fwdtallies_validation_data$Y)

fwdtallies_plot_1 = 
  fwdtallies_sim %>%
  group_by(obs_id, Y_obs, group_id) %>%
  summarise(Y_sim = mean(Y_sim)) %>%
  ungroup() %>%
  ggplot(aes(x = Y_sim, y = Y_obs))+
  geom_abline(intercept=0,slope=1, color = 'red')+
  geom_smooth(method = 'lm')+
  geom_point(size = 0)+
  theme_minimal()+
  coord_fixed()+
  facet_wrap(~group_id)+
  labs(x = 'Mean simulated FWD tally',
       y = 'Observed (validation data) FWD tally')

fwdtallies_plot_1

ggsave(fwdtallies_plot_1,
       filename = here::here('04-communication','figures', 'manuscript',
                             'fwdtallies_prediction_1.png'),
       height = 4, width = 6.5, units = 'in')

fwdtallies_plot_2 = 
  fwdtallies_sim %>%
  group_by(obs_id, Y_obs, group_id) %>%
  summarise() %>%
  ungroup() %>%
  mutate(source = 'obs') %>%
  group_by(Y_obs, group_id, source) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  left_join(
    fwdtallies_sim %>%
      group_by(obs_id, Y_obs, group_id) %>%
      summarise() %>%
      ungroup() %>%
      mutate(source = 'obs') %>%
      group_by(group_id, source) %>%
      summarise(total = n()) %>%
      ungroup() 
  ) %>%
  mutate(proportion = count / total) %>%
  select(group_id, Y = Y_obs, proportion, source) %>%
  bind_rows(
    fwdtallies_sim %>%
      group_by(Y_sim, group_id) %>%
      summarise(count = n()) %>%
      ungroup() %>%
      mutate(source = 'sim') %>%
      left_join(
        fwdtallies_sim %>%
          group_by(group_id) %>%
          summarise(total = n()) %>%
          ungroup() %>%
          mutate(source = 'sim')
      ) %>%
      mutate(proportion = count / total) %>%
      select(group_id, Y = Y_sim, proportion, source)
  ) %>%
  ggplot()+
  geom_col(aes(x = Y, y = proportion, fill = source),
           position = position_dodge())+
  facet_grid(group_id~., scales = 'free')+
  theme_minimal()+
  labs(title = 'fwdtallies prediction')+
  scale_x_continuous(limits = c(0, 150))

fwdtallies_plot_2

ggsave(fwdtallies_plot_2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwdtallies_prediction_2.png'),
       height = 6.5, width = 6.5, units = 'in')

fwdtallies_plot_3 = 
  fwdtallies_sim %>%
  group_by(group_id, obs_id, Y_obs) %>%
  summarise(logmu.mean = mean(logmu),
            Y_sim.mean = mean(Y_sim),
            Y_sim.05 = quantile(Y_sim, 0.05),
            Y_sim.95 = quantile(Y_sim, 0.95)) %>%
  ungroup() %>%
  ggplot()+
  geom_errorbar(aes(x = logmu.mean, ymin = Y_sim.05, ymax = Y_sim.95),
                alpha = 0.75)+
  geom_point(aes(x = logmu.mean, y = Y_obs), color = 'red')+
  theme_minimal()+
  labs(title = 'fwdtallies prediction')+
  facet_grid(group_id~., scales = 'free')

fwdtallies_plot_3

ggsave(fwdtallies_plot_3,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwdtallies_prediction_3.png'),
       height = 6.5, width = 6.5, units = 'in')


fwdtallies_plot_4 = 
  ggplot()+
  geom_errorbar(
    data = 
      fwdtallies_sim %>%
      group_by(group_id, Y_sim, draw) %>%
      summarise(N_obs = n()) %>%
      ungroup() %>%
      right_join(fwdtallies_sim %>%
                   expand(nesting(group_id, draw),
                          Y_sim = 0:max(.$Y_sim))) %>%
      mutate(N_obs = ifelse(is.na(N_obs), 0, N_obs)) %>%
      group_by(group_id, Y_sim) %>%
      summarise(N_obs.025 = quantile(N_obs, 0.025),
                N_obs.975 = quantile(N_obs, 0.975)) %>%
      ungroup(),
    aes(x = Y_sim, ymin = N_obs.025, ymax = N_obs.975, color = 'Simulated')
  )+
  facet_grid(group_id~., scales = 'free')+
  theme_minimal()+
  labs(x = 'FWD tally', y = 'N observations')+
  geom_col(
    data = 
      fwdtallies_sim %>%
      group_by(group_id, obs_id, Y_obs) %>%
      summarise() %>%
      ungroup() %>%
      group_by(group_id, Y_obs) %>%
      summarise(count = n()) %>%
      ungroup(),
    aes(x = Y_obs, y = count, fill = 'Observed (validation data)'),
    width = 0.5,
    alpha = 0.75
  )+
  scale_fill_manual(values = c('Observed (validation data)' = 'red'))+
  scale_color_manual(values = c('Simulated' = 'blue'))+
  theme(legend.position = 'bottom',
        legend.title = element_blank())+
  scale_x_continuous(limits = c(-1, 100))

fwdtallies_plot_4

ggsave(fwdtallies_plot_4,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwdtallies_prediction_4.png'),
       height = 6.5, width = 6.5, units = 'in')


fwdtallies_plot_5 = 
  ggplot()+
  stat_ecdf(
    data = fwdtallies_sim,
    aes(x = Y_sim, group = draw, color = 'Simulated'),
    geom = 'step', 
    pad = FALSE,
    #color = 'blue',
    lwd = 1,
    alpha = 0.01)+
  stat_ecdf(
    data = fwdtallies_sim %>% group_by(group_id, Y_obs, obs_id) %>% summarise() %>% ungroup(),
    aes(x = Y_obs, color = 'Observed (validation data)'),
    geom = 'step',
    pad = 'false',
    #color = 'red',
    lwd = 1
  )+
  theme_minimal()+
  scale_color_manual(values = c('Observed (validation data)' = 'red', 'Simulated' = 'blue'))+
  facet_grid(group_id~.)+
  labs(x = 'FWD tally', y = 'Cumulative density')+
  scale_x_continuous(limits = c(-1, 100))+
  theme(legend.position = 'bottom',
        legend.title = element_blank())

fwdtallies_plot_5

ggsave(fwdtallies_plot_5,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwdtallies_prediction_5.png'),
       height = 6.5, width = 6.5, units = 'in')

rm(fwdtallies_sim)
#### fwddiams ##################################################################

fwddiams_validation_data = 
  readRDS(here::here('02-data',
                     '05-for_analysis',
                     'fwddiams_validation_data.rds'))

fwddiams_posterior = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'fwddiams_fit.rds'))$draws() %>%
  as_draws_df()

set.seed(110819)
fwddiams_sim = 
  do.call('rbind',
          lapply(X = 1:4000,
                 FUN = 
                   function(draw){
                     print(paste0('Working on draw', draw))
                     
                     intercept = 
                       fwddiams_posterior %>%
                       select(contains('intercept')) %>%
                       slice(draw) %>%
                       as.data.frame() %>%
                       as.numeric()
                     
                     phi = 
                       fwddiams_posterior %>%
                       select(contains('phi')) %>%
                       slice(draw) %>%
                       as.data.frame() %>%
                       as.numeric()
                     
                     effect_plot = 
                       fwddiams_posterior %>%
                       select(contains('effect_plot')) %>%
                       slice(draw) %>%
                       as.data.frame() %>%
                       as.numeric()
                     
                     effect_trans = 
                       fwddiams_posterior %>%
                       select(contains('effect_trans')) %>%
                       slice(draw) %>%
                       as.data.frame() %>%
                       as.numeric()
                     
                     N = fwddiams_validation_data$N
                     group_id = fwddiams_validation_data$group_id
                     plot_id = fwddiams_validation_data$plot_id
                     transect_id = fwddiams_validation_data$transect_id
                     
                     logmu = 
                       sapply(X = 1:fwddiams_validation_data$N,
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
                     
                     scale = 
                       sapply(X = 1:N,
                              FUN = function(i){
                                exp(logmu[i])*(1+phi[group_id[plot_id[transect_id[i]]]])
                              })
                     
                     Ysim = 
                       sapply(X = 1:N,
                              FUN = function(i){
                                y = 0
                                while (y<=fwddiams_validation_data$lb |
                                       y>=fwddiams_validation_data$ub){
                                  
                                  # this is really confusing, but it looks like 
                                  # the rinvgamma function calls the beta 
                                  # parameter the rate, and 1/beta the scale.
                                  # comparing the PDF formulae given in 
                                  # the rinvgamma documentation, bourgingnon 2020,
                                  # and the stan invgamma documentation, its 
                                  # clear that what stan and bourgingnon call the
                                  # scale is called the rate in rinvgamma
                                  y = 
                                    invgamma::rinvgamma(n = 1,
                                                        shape = shape[i],
                                                        rate = scale[i])
                                }
                                return(y)
                              })
                     
                     results = 
                       data.frame(Yobs = fwddiams_validation_data$Y,
                                  transect_id = transect_id,
                                  plot_id = plot_id[transect_id],
                                  group_id = group_id[plot_id[transect_id]],
                                  Ysim = Ysim,
                                  draw = draw) %>%
                       rowid_to_column('obs_id')
                     
                     
                     return(results)
                   }))

saveRDS(fwddiams_sim,
        here::here('02-data',
                   '06-results',
                   'posterior_predictions',
                   'fwddiams_sim.rds'))

fwddiams_sim = 
  readRDS(here::here('02-data',
                     '06-results',
                     'posterior_predictions',
                     'fwddiams_sim.rds'))

fwddiams_plot_1 = 
  fwddiams_sim %>%
  filter(draw<=4000) %>%
  ggplot()+
  geom_line(aes(color = 'Simulated', group = draw, x = Ysim),
            stat = 'density',
            alpha = 0.1)+
  geom_line(
    stat = 'density',
    data = 
      fwddiams_sim %>%
      group_by(group_id, obs_id, Yobs) %>%
      summarise() %>%
      ungroup(),
    aes(x = Yobs,
        color = 'Observed (validation data)'),
    lwd = 1
  )+
  facet_grid(group_id~.)+
  scale_color_manual(values = c('Observed (validation data)' = 'red', 'Simulated' = 'blue'))+
  theme_minimal()+
  labs(x = 'FWD diameter (cm)', y = 'Probability density')+
  theme(legend.position = 'bottom',
        legend.title = element_blank())

fwddiams_plot_1

ggsave(fwddiams_plot_1,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwddiams_prediction_1.png'),
       height = 6.5, width = 6.5, units = 'in')


fwddiams_plot_2 = 
  fwddiams_sim %>%
  filter(draw<=4000) %>%
  ggplot()+
  stat_ecdf(aes(x = Ysim, group = draw, color = 'Simulated'),
            alpha = 0.1,
            pad = FALSE)+
  stat_ecdf(
    data = fwddiams_sim %>%
      group_by(group_id, obs_id, Yobs) %>%
      summarise() %>%
      ungroup(),
    aes(x = Yobs, color = 'Observed (validation data)'),
    lwd = 1,
    pad = FALSE
  )+
  facet_grid(group_id~.)+
  theme_minimal()+
  scale_color_manual(values = c('Observed (validation data)' = 'red', 'Simulated' = 'blue'))+
  labs(x = 'FWD diameter (cm)', y = 'Cumulative density')+
  theme(legend.position = 'bottom',
        legend.title = element_blank())

fwddiams_plot_2

ggsave(fwddiams_plot_2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwddiams_prediction_2.png'),
       height = 6.5, width = 6.5, units = 'in')

rm(fwddiams_sim)

#### cwdtallies ################################################################


cwdtallies_validation_data = 
  readRDS(here::here('02-data', 
                     '05-for_analysis',
                     'cwdtallies_validation_data.rds'))

cwdtallies_posterior = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'cwdtallies_fit.rds'))$draws() %>%
  as_draws_df()

set.seed(110819)
cwdtallies_sim = 
  do.call('bind_rows',
          lapply(X = 1:4000,
                 FUN = function(draw){
                   
                   print(paste0('Working on draw ', draw))
                   
                   # extract parameters for this draw
                   beta = 
                     cwdtallies_posterior %>%
                     select(contains('beta')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   alpha = 
                     cwdtallies_posterior %>%
                     select(contains('alpha')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   rho = 
                     cwdtallies_posterior %>%
                     select(contains('rho')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   plotEffects = 
                     cwdtallies_posterior %>%
                     select(contains('plot_effect')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   transectEffects = 
                     cwdtallies_posterior %>%
                     select(contains('transect_effect')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   zGP = 
                     matrix(nrow = cwdtallies_validation_data$L,
                            ncol = cwdtallies_validation_data$P,
                            byrow = FALSE,
                            data = 
                              cwdtallies_posterior %>%
                              select(contains('zGP')) %>%
                              slice(draw) %>%
                              as.data.frame() %>%
                              as.numeric())
                   
                   
                   dists = 
                     cwdtallies_validation_data$coords %>%
                     dist() %>%
                     as.matrix()
                   
                   # construct GP effects
                   K = 
                     array(dim = 
                             c(cwdtallies_validation_data$L,
                               cwdtallies_validation_data$L,
                               cwdtallies_validation_data$G),
                           dimnames = 
                             list('l' = 1:cwdtallies_validation_data$L,
                                  'lprime' = 1:cwdtallies_validation_data$L,
                                  'g' = 1:cwdtallies_validation_data$G),
                           data = 
                             sapply(X = 1:cwdtallies_validation_data$G,
                                    FUN = function(g){
                                      sapply(X = 1:cwdtallies_validation_data$L,
                                             FUN = function(l){
                                               sapply(X = 1:cwdtallies_validation_data$L,
                                                      FUN = function(lprime){
                                                        ((alpha[g]**2)*
                                                           exp(-(dists[lprime,l]**2)/
                                                                 (2*rho[g]**2)))
                                                      })
                                             })
                                    }))
                   
                   for (g in 1:cwdtallies_validation_data$G){
                     K[,,g] = K[,,g]+diag(nrow = nrow(K[,,g]),
                                          ncol = ncol(K[,,g]),
                                          x = 1e-9)
                   }
                   
                   KL = 
                     array(dim = 
                             c(cwdtallies_validation_data$L,
                               cwdtallies_validation_data$L,
                               cwdtallies_validation_data$G),
                           dimnames = 
                             list('l' = 1:cwdtallies_validation_data$L,
                                  'lprime' = 1:cwdtallies_validation_data$L,
                                  'g' = 1:cwdtallies_validation_data$G),
                           data = 
                             sapply(X = 1:cwdtallies_validation_data$G,
                                    FUN = function(g){
                                      t(chol(K[,,g]))
                                    }))
                   
                   # realized random effects
                   GP = 
                     matrix(nrow = cwdtallies_validation_data$L,
                            ncol = cwdtallies_validation_data$P,
                            data = 
                              sapply(X = 1:cwdtallies_validation_data$P,
                                     FUN = function(p){
                                       KL[,,cwdtallies_validation_data$group_id[p]] %*%
                                         zGP[,p]
                                     }))
                   
                   # linear predictors
                   logmu = 
                     sapply(X = 1:cwdtallies_validation_data$N,
                            FUN = function(i){
                              
                              beta[
                                cwdtallies_validation_data$group_id[
                                  cwdtallies_validation_data$plot_id[
                                    cwdtallies_validation_data$transect_id[i]]]] +
                                plotEffects[
                                  cwdtallies_validation_data$plot_id[
                                    cwdtallies_validation_data$transect_id[i]]]+
                                transectEffects[
                                  cwdtallies_validation_data$transect_id[i]]+
                                as.numeric(GP[cwdtallies_validation_data$location_id[i],
                                              cwdtallies_validation_data$plot_id][
                                                cwdtallies_validation_data$transect_id[i]
                                              ])
                              
                            })
                   
                   # realizations
                   Y = 
                     sapply(X = 1:cwdtallies_validation_data$N,
                            FUN = function(i){
                              rpois(n = 1,
                                    lambda = exp(logmu[i]))})
                   
                   
                   results = 
                     data.frame(Y_obs = cwdtallies_validation_data$Y,
                                Y_sim = Y,
                                group_id = 
                                  cwdtallies_validation_data$group_id[
                                    cwdtallies_validation_data$plot_id[
                                      cwdtallies_validation_data$transect_id]],
                                plot_id = 
                                  cwdtallies_validation_data$plot_id[
                                    cwdtallies_validation_data$transect_id],
                                transect_id = 
                                  cwdtallies_validation_data$transect_id,
                                location_id = cwdtallies_validation_data$location_id,
                                logmu = logmu,
                                draw = draw) %>%
                     as_tibble() %>%
                     rowid_to_column('obs_id')
                   
                   return(results)
                 }))

saveRDS(cwdtallies_sim,
        here::here('02-data',
                   '06-results',
                   'posterior_predictions',
                   'cwdtallies_sim.rds'))

cwdtallies_sim = 
  readRDS(here::here('02-data',
                     '06-results',
                     'posterior_predictions',
                     'cwdtallies_sim.rds'))


cwdtallies_mae = 
  cwdtallies_sim %>%
  group_by(obs_id, Y_obs, group_id) %>%
  summarise(Y_sim = mean(Y_sim)) %>%
  ungroup() %>%
  mutate(ae = abs(Y_obs-Y_sim)) %>%
  pull(ae) %>%
  mean()

cwdtallies_plot_1 = 
  cwdtallies_sim %>%
  group_by(obs_id, Y_obs, group_id) %>%
  summarise(Y_sim = mean(Y_sim)) %>%
  ungroup() %>%
  ggplot(aes(x = Y_sim, y = Y_obs))+
  geom_abline(intercept=0,slope=1, color = 'red')+
  geom_smooth(method = 'lm')+
  geom_point(size = 0)+
  theme_minimal()+
  coord_fixed()+
  facet_wrap(~group_id)+
  labs(x = 'Mean simulated CWD tally',
       y = 'Observed (validation data) CWD tally')

cwdtallies_plot_1

ggsave(cwdtallies_plot_1,
       filename = here::here('04-communication','figures', 'manuscript',
                             'cwdtallies_prediction_1.png'),
       height = 6.5, width = 6.5, units = 'in')

cwdtallies_plot_2 = 
  cwdtallies_sim %>%
  group_by(obs_id, Y_obs, group_id) %>%
  summarise() %>%
  ungroup() %>%
  mutate(source = 'obs') %>%
  group_by(Y_obs, group_id, source) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  left_join(
    cwdtallies_sim %>%
      group_by(obs_id, Y_obs, group_id) %>%
      summarise() %>%
      ungroup() %>%
      mutate(source = 'obs') %>%
      group_by(group_id, source) %>%
      summarise(total = n()) %>%
      ungroup() 
  ) %>%
  mutate(proportion = count / total) %>%
  select(group_id, Y = Y_obs, proportion, source) %>%
  bind_rows(
    cwdtallies_sim %>%
      group_by(Y_sim, group_id) %>%
      summarise(count = n()) %>%
      ungroup() %>%
      mutate(source = 'sim') %>%
      left_join(
        cwdtallies_sim %>%
          group_by(group_id) %>%
          summarise(total = n()) %>%
          ungroup() %>%
          mutate(source = 'sim')
      ) %>%
      mutate(proportion = count / total) %>%
      select(group_id, Y = Y_sim, proportion, source)
  ) %>%
  ggplot()+
  geom_col(aes(x = Y, y = proportion, fill = source),
           position = position_dodge())+
  facet_grid(group_id~., scales = 'free')+
  theme_minimal()+
  labs(title = 'cwdtallies prediction')

cwdtallies_plot_2

ggsave(cwdtallies_plot_2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'cwdtallies_prediction_2.png'),
       height = 6.5, width = 6.5, units = 'in')

cwdtallies_plot_3 = 
  cwdtallies_sim %>%
  group_by(group_id, obs_id, Y_obs) %>%
  summarise(logmu.mean = mean(logmu),
            Y_sim.mean = mean(Y_sim),
            Y_sim.05 = quantile(Y_sim, 0.05),
            Y_sim.95 = quantile(Y_sim, 0.95)) %>%
  ungroup() %>%
  ggplot()+
  geom_errorbar(aes(x = logmu.mean, ymin = Y_sim.05, ymax = Y_sim.95),
                alpha = 0.75)+
  geom_point(aes(x = logmu.mean, y = Y_obs), color = 'red')+
  theme_minimal()+
  labs(title = 'cwdtallies prediction')+
  facet_grid(group_id~., scales = 'free')

cwdtallies_plot_3

ggsave(cwdtallies_plot_3,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'cwdtallies_prediction_3.png'),
       height = 6.5, width = 6.5, units = 'in')


cwdtallies_plot_4 = 
  ggplot()+
  geom_errorbar(
    data = 
      cwdtallies_sim %>%
      group_by(group_id, Y_sim, draw) %>%
      summarise(N_obs = n()) %>%
      ungroup() %>%
      right_join(cwdtallies_sim %>%
                   expand(nesting(group_id, draw),
                          Y_sim = 0:max(.$Y_sim))) %>%
      mutate(N_obs = ifelse(is.na(N_obs), 0, N_obs)) %>%
      group_by(group_id, Y_sim) %>%
      summarise(N_obs.025 = quantile(N_obs, 0.025),
                N_obs.975 = quantile(N_obs, 0.975)) %>%
      ungroup(),
    aes(x = Y_sim, ymin = N_obs.025, ymax = N_obs.975, color = 'Simulated')
  )+
  facet_grid(group_id~., scales = 'free')+
  theme_minimal()+
  labs(x = 'CWD tally', y = 'N observations')+
  geom_col(
    data = 
      cwdtallies_sim %>%
      group_by(group_id, obs_id, Y_obs) %>%
      summarise() %>%
      ungroup() %>%
      group_by(group_id, Y_obs) %>%
      summarise(count = n()) %>%
      ungroup(),
    aes(x = Y_obs, y = count, fill = 'Observed (validation data)'),
    width = 0.5,
    alpha = 0.75
  )+
  scale_fill_manual(values = c('Observed (validation data)' = 'red'))+
  scale_color_manual(values = c('Simulated' = 'blue'))+
  theme(legend.position = 'bottom',
        legend.title = element_blank())

cwdtallies_plot_4
ggsave(cwdtallies_plot_4,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'cwdtallies_prediction_4.png'),
       height = 6.5, width = 6.5, units = 'in')

head(cwdtallies_sim)

cwdtallies_plot_5 = 
  ggplot()+
  stat_ecdf(
    data = cwdtallies_sim,
    aes(x = Y_sim, group = draw, color = 'Simulated'),
    geom = 'step', 
    pad = FALSE,
    #color = 'blue',
    lwd = 1,
    alpha = 0.01)+
  stat_ecdf(
    data = cwdtallies_sim %>% group_by(group_id, Y_obs, obs_id) %>% summarise() %>% ungroup(),
    aes(x = Y_obs, color = 'Observed (validation data)'),
    geom = 'step',
    pad = 'false',
    #color = 'red',
    lwd = 1
  )+
  theme_minimal()+
  scale_color_manual(values = c('Observed (validation data)' = 'red', 'Simulated' = 'blue'))+
  facet_grid(group_id~.)+
  labs(x = 'CWD tally', y = 'Cumulative density')+
  theme(legend.position = 'bottom',
        legend.title = element_blank())

cwdtallies_plot_5

ggsave(cwdtallies_plot_5,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'cwdtallies_prediction_5.png'),
       height = 6.5, width = 6.5, units = 'in')


rm(cwdtallies_sim)

#### cwddiams ##################################################################

cwddiams_validation_data = 
  readRDS(here::here('02-data',
                     '05-for_analysis',
                     'cwddiams_validation_data.rds'))

cwddiams_posterior = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'cwddiams_fit.rds'))$draws() %>%
  as_draws_df()

set.seed(110819)
cwddiams_sim = 
  do.call('rbind',
          lapply(X = 1:4000,
                 FUN = 
                   function(draw){
                     
                     print(paste0('Working on draw ', draw))
                     
                     mu = 
                       cwddiams_posterior %>%
                       select(contains('mu')) %>%
                       slice(draw) %>%
                       as.data.frame() %>%
                       as.numeric()
                     
                     phi = 
                       cwddiams_posterior %>%
                       select(contains('phi')) %>%
                       slice(draw) %>%
                       as.data.frame() %>%
                       as.numeric()
                     
                     N = cwddiams_validation_data$N
                     group_id = cwddiams_validation_data$group_id
                     plot_id = cwddiams_validation_data$plot_id
                     
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
                                while (y<=cwddiams_validation_data$lb |
                                       y>=cwddiams_validation_data$ub){
                                  
                                  y = 
                                    invgamma::rinvgamma(n = 1,
                                                        shape = shape[i],
                                                        rate = rate[i])
                                }
                                return(y)
                              })
                     
                     results = 
                       data.frame(Yobs = cwddiams_validation_data$Y,
                                  plot_id = plot_id,
                                  group_id = group_id[plot_id],
                                  Ysim = Ysim,
                                  draw = draw) %>%
                       rowid_to_column('obs_id')
                     
                     
                     return(results)
                   }))

saveRDS(cwddiams_sim,
        here::here('02-data',
                   '06-results',
                   'posterior_predictions',
                   'cwddiams_sim.rds'))

cwddiams_sim = 
  readRDS(here::here('02-data',
                     '06-results',
                     'posterior_predictions',
                     'cwddiams_sim.rds'))

cwddiams_plot_1 = 
  cwddiams_sim %>%
  filter(draw<=4000) %>%
  ggplot()+
  geom_line(aes(color = 'Simulated', group = draw, x = Ysim),
            stat = 'density',
            alpha = 0.1)+
  geom_line(
    stat = 'density',
    data = 
      cwddiams_sim %>%
      group_by(group_id, obs_id, Yobs) %>%
      summarise() %>%
      ungroup(),
    aes(x = Yobs,
        color = 'Observed (validation data)'),
    lwd = 1
  )+
  facet_grid(group_id~.)+
  scale_color_manual(values = c('Observed (validation data)' = 'red', 'Simulated' = 'blue'))+
  theme_minimal()+
  labs(x = 'CWD diam (cm)', y = 'Probability density')+
  theme(legend.position = 'bottom',
        legend.title = element_blank())

cwddiams_plot_1

ggsave(cwddiams_plot_1,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'cwddiams_prediction_1.png'),
       height = 6.5, width = 6.5, units = 'in')


cwddiams_plot_2 = 
  cwddiams_sim %>%
  filter(draw<=4000) %>%
  ggplot()+
  stat_ecdf(aes(x = Ysim, group = draw, color = 'Simulated'),
            alpha = 0.1,
            pad = FALSE)+
  stat_ecdf(
    data = cwddiams_sim %>%
      group_by(group_id, obs_id, Yobs) %>%
      summarise() %>%
      ungroup(),
    aes(x = Yobs, color = 'Observed (validation data)'),
    lwd = 1,
    pad = FALSE
  )+
  facet_grid(group_id~.)+
  theme_minimal()+
  scale_color_manual(values = c('Observed (validation data)' = 'red', 'Simulated' = 'blue'))+
  labs(x = 'CWD diam (cm)', y = 'Cumulative density function')+
  theme(legend.position = 'bottom',
        legend.title = element_blank())

cwddiams_plot_2

ggsave(cwddiams_plot_2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'cwddiams_prediction_2.png'),
       height = 6.5, width = 6.5, units = 'in')

rm(cwddiams_sim)

#### vegpa #####################################################################

vegpa_validation_data = 
  readRDS(here::here('02-data', 
                     '05-for_analysis',
                     'vegpa_validation_data.rds'))

vegpa_posterior = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'vegpa_fit.rds'))$draws() %>%
  as_draws_df()

set.seed(110819)
vegpa_sim = 
  do.call('bind_rows',
          lapply(X = 1:4000,
                 FUN = function(draw){
                   
                   print(paste0('Working on draw ', draw))
                   
                   # extract parameters for this draw
                   beta = 
                     vegpa_posterior %>%
                     select(contains('beta')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   alpha = 
                     vegpa_posterior %>%
                     select(contains('alpha')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   rho = 
                     vegpa_posterior %>%
                     select(contains('rho')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   plotEffects = 
                     vegpa_posterior %>%
                     select(contains('plot_effect')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   zGP = 
                     matrix(nrow = vegpa_validation_data$L,
                            ncol = vegpa_validation_data$P,
                            byrow = FALSE,
                            data = 
                              vegpa_posterior %>%
                              select(contains('zGP')) %>%
                              slice(draw) %>%
                              as.data.frame() %>%
                              as.numeric())
                   
                   
                   dists = 
                     vegpa_validation_data$coords %>%
                     dist() %>%
                     as.matrix()
                   
                   # construct GP effects
                   K = 
                     array(dim = 
                             c(vegpa_validation_data$L,
                               vegpa_validation_data$L,
                               vegpa_validation_data$G),
                           dimnames = 
                             list('l' = 1:vegpa_validation_data$L,
                                  'lprime' = 1:vegpa_validation_data$L,
                                  'g' = 1:vegpa_validation_data$G),
                           data = 
                             sapply(X = 1:vegpa_validation_data$G,
                                    FUN = function(g){
                                      sapply(X = 1:vegpa_validation_data$L,
                                             FUN = function(l){
                                               sapply(X = 1:vegpa_validation_data$L,
                                                      FUN = function(lprime){
                                                        ((alpha[g]**2)*
                                                           exp(-(dists[lprime,l]**2)/
                                                                 (2*rho[g]**2)))
                                                      })
                                             })
                                    }))
                   
                   for (g in 1:vegpa_validation_data$G){
                     K[,,g] = K[,,g]+diag(nrow = nrow(K[,,g]),
                                          ncol = ncol(K[,,g]),
                                          x = 1e-9)
                   }
                   
                   KL = 
                     array(dim = 
                             c(vegpa_validation_data$L,
                               vegpa_validation_data$L,
                               vegpa_validation_data$G),
                           dimnames = 
                             list('l' = 1:vegpa_validation_data$L,
                                  'lprime' = 1:vegpa_validation_data$L,
                                  'g' = 1:vegpa_validation_data$G),
                           data = 
                             sapply(X = 1:vegpa_validation_data$G,
                                    FUN = function(g){
                                      t(chol(K[,,g]))
                                    }))
                   
                   # realized random effects
                   GP = 
                     matrix(nrow = vegpa_validation_data$L,
                            ncol = vegpa_validation_data$P,
                            data = 
                              sapply(X = 1:vegpa_validation_data$P,
                                     FUN = function(p){
                                       KL[,,vegpa_validation_data$group_id[p]] %*%
                                         zGP[,p]
                                     }))
                   
                   # linear predictors
                   logiteta = 
                     sapply(X = 1:vegpa_validation_data$N,
                            FUN = function(i){
                              
                              beta[vegpa_validation_data$group_id[vegpa_validation_data$plot_id[i]]] +
                                plotEffects[vegpa_validation_data$plot_id[i]]+
                                as.numeric(GP[vegpa_validation_data$location_id[i],
                                              vegpa_validation_data$plot_id][i])
                              
                            })
                   
                   # realizations
                   Y = 
                     sapply(X = 1:vegpa_validation_data$N,
                            FUN = function(i){
                              rbinom(n = 1,
                                     size = 1,
                                     prob = boot::inv.logit(logiteta[i]))})
                   
                   
                   results = 
                     data.frame(Y_obs = vegpa_validation_data$Y,
                                Y_sim = Y,
                                group_id = 
                                  vegpa_validation_data$group_id[vegpa_validation_data$plot_id],
                                plot_id = vegpa_validation_data$plot_id,
                                location_id = vegpa_validation_data$location_id,
                                logiteta = logiteta,
                                draw = draw) %>%
                     as_tibble() %>%
                     rowid_to_column('obs_id')
                   
                   return(results)
                 }))

saveRDS(vegpa_sim,
        here::here('02-data',
                   '06-results',
                   'posterior_predictions',
                   'vegpa_sim.rds'))

vegpa_sim = readRDS(here::here('02-data',
                               '06-results',
                               'posterior_predictions',
                               'vegpa_sim.rds'))



vegpa_prediction_1 = 
  vegpa_sim %>%
  mutate(p = boot::inv.logit(logiteta)) %>%
  group_by(group_id, obs_id, Y_obs) %>%
  summarise(p.mean = mean(p),
            p.50 = median(p),
            p.05 = quantile(p, 0.05),
            p.95 = quantile(p, 0.95)) %>%
  ungroup() %>%
  group_by(group_id) %>%
  mutate(r = dense_rank(p.mean),
         r_bin = cut(r, breaks = seq(from = 0, 
                                     to = n()+1, length.out = 10))) %>%
  ungroup() %>%
  ggplot()+
  geom_jitter(aes(x = r, y = Y_obs, color = 'Observed (validation data)'), 
              size = 0, width = 0, height = 0.1)+
  geom_ribbon(aes(x = r, ymin = p.05, ymax = p.95), alpha = 0.1)+
  geom_point(aes(x = r, y = p.50, color = 'Simulated'), size = 0)+
  geom_point(data = 
               vegpa_sim %>%
               mutate(p = boot::inv.logit(logiteta)) %>%
               group_by(group_id, obs_id, Y_obs) %>%
               summarise(p.mean = mean(p),
                         p.50 = median(p),
                         p.05 = quantile(p, 0.05),
                         p.95 = quantile(p, 0.95)) %>%
               ungroup() %>%
               group_by(group_id) %>%
               mutate(r = dense_rank(p.mean),
                      r_bin = cut(r, breaks =
                                    seq(from = 0, to = n()+1, length.out = 10))) %>%
               ungroup() %>%
               group_by(group_id, r_bin) %>%
               summarise(r_loc = mean(r),
                         p.actual = mean(Y_obs)) %>%
               ungroup(),
             aes(x = r_loc, y = p.actual, color = 'Observed (validation data)'),
             pch = 18, size = 3)+
  
  facet_grid(group_id~.)+
  theme_minimal()+
  labs(x = 'Rank (by probability of presence)', y = 'Presence of understory vegetation (validation data)')+
  scale_color_manual(values = c('Observed (validation data)' = 'red', 'Simulated' = 'blue'))+
  theme(legend.position = 'bottom',
        legend.title = element_blank())

vegpa_prediction_1

ggsave(vegpa_prediction_1,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'vegpa_prediction_1.png'),
       height = 6.5, width = 6.5, units = 'in')

ggplot()+
  geom_jitter(
    data = 
      vegpa_sim %>%
      mutate(p = boot::inv.logit(logiteta)) %>%
      group_by(group_id, obs_id, Y_obs) %>%
      summarise(p.mean = mean(p),
                p.50 = median(p),
                p.05 = quantile(p, 0.05),
                p.95 = quantile(p, 0.95)) %>%
      ungroup() %>%
      mutate(p.bin = cut(p.50, seq(0, 1, 0.1))),
    aes(x = p.50, y = Y_obs),
    color = 'red', size = 0, height = 0.1, width = 0
  )+
  geom_point(
    data = 
      vegpa_sim %>%
      mutate(p = boot::inv.logit(logiteta)) %>%
      group_by(group_id, obs_id, Y_obs) %>%
      summarise(p.50 = median(p),
                p.05 = quantile(p, 0.05),
                p.95 = quantile(p, 0.95)) %>%
      ungroup() %>%
      mutate(p.bin = round(p.50, 1)) %>%
      group_by(group_id, p.bin) %>%
      summarise(p.50 = median(p.50),
                p.05 = median(p.05),
                p.95 = median(p.95)) %>%
      ungroup(),
    aes(x = p.bin, y = p.50),
    size = 3
  )+
  geom_errorbar(
    data = 
      vegpa_sim %>%
      mutate(p = boot::inv.logit(logiteta)) %>%
      group_by(group_id, obs_id, Y_obs) %>%
      summarise(p.50 = median(p),
                p.05 = quantile(p, 0.05),
                p.95 = quantile(p, 0.95)) %>%
      ungroup() %>%
      mutate(p.bin = round(p.50, 1)) %>%
      group_by(group_id, p.bin) %>%
      summarise(p.50 = median(p.50),
                p.05 = median(p.05),
                p.95 = median(p.95)) %>%
      ungroup(),
    aes(x = p.bin, ymin = p.05, ymax = p.95),
    width = 0.05
  )+
  geom_point(
    data = 
      vegpa_sim %>%
      mutate(p = boot::inv.logit(logiteta)) %>%
      group_by(group_id, obs_id, Y_obs) %>%
      summarise(p.mean = mean(p),
                p.50 = median(p),
                p.05 = quantile(p, 0.05),
                p.95 = quantile(p, 0.95)) %>%
      ungroup() %>%
      mutate(p.bin = round(p.50, 1)) %>%
      group_by(group_id, p.bin) %>%
      summarise(Y_obs.mean = mean(Y_obs)) %>%
      ungroup(),
    aes(x = p.bin, y = Y_obs.mean),
    pch = 18, color = 'blue', size = 3
  )+
  facet_grid(group_id~.)+
  theme_minimal()

head(vegpa_sim)

ggplot()+
  stat_ecdf(
    data = vegpa_sim %>%
      group_by(group_id, draw, plot_id) %>%
      summarise(percent_cover = mean(Y_sim)) %>%
      ungroup(),
    aes(x = percent_cover, group = draw, color = 'Simulated'),
    geom = 'step', 
    pad = FALSE,
    #color = 'blue',
    lwd = 1,
    alpha = 0.01)+
  stat_ecdf(
    data = 
      vegpa_sim %>%
      group_by(group_id, plot_id, obs_id, Y_obs) %>%
      summarise() %>%
      ungroup() %>%
      group_by(group_id, plot_id) %>%
      summarise(percent_cover = mean(Y_obs)) %>%
      ungroup(),
    aes(x = percent_cover, color = 'Observed (validation data)'),
    geom = 'step',
    pad = 'false',
    #color = 'red',
    lwd = 1
  )+
  theme_minimal()+
  scale_color_manual(values = c('Observed (validation data)' = 'red', 'Simulated' = 'blue'))+
  facet_grid(group_id~.)+
  labs(x = 'Plot-level % cover', y = 'Cumulative density', title = 'Veg pres/abs prediction')

rm(vegpa_sim)

#### trees #####################################################################


trees_validation_data = 
  readRDS(here::here('02-data', 
                     '05-for_analysis',
                     'trees_validation_data.rds'))

trees_posterior = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'trees_fit.rds'))$draws() %>%
  as_draws_df()

set.seed(110819)
trees_sim = 
  do.call('bind_rows',
          lapply(X = 1:4000,
                 FUN = function(draw){
                   
                   print(paste0('Working on draw', draw))
                   
                   # extract parameters for this draw
                   beta = 
                     trees_posterior %>%
                     select(contains('beta')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   alpha = 
                     trees_posterior %>%
                     select(contains('alpha')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   rho = 
                     trees_posterior %>%
                     select(contains('rho')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   plotEffects = 
                     trees_posterior %>%
                     select(contains('plot_effect')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   zGP = 
                     matrix(nrow = trees_validation_data$L,
                            ncol = trees_validation_data$P,
                            byrow = FALSE,
                            data = 
                              trees_posterior %>%
                              select(contains('zGP')) %>%
                              slice(draw) %>%
                              as.data.frame() %>%
                              as.numeric())
                   
                   
                   dists = 
                     trees_validation_data$coords %>%
                     dist() %>%
                     as.matrix()
                   
                   # construct GP effects
                   K = 
                     array(dim = 
                             c(trees_validation_data$L,
                               trees_validation_data$L,
                               trees_validation_data$G),
                           dimnames = 
                             list('l' = 1:trees_validation_data$L,
                                  'lprime' = 1:trees_validation_data$L,
                                  'g' = 1:trees_validation_data$G),
                           data = 
                             sapply(X = 1:trees_validation_data$G,
                                    FUN = function(g){
                                      sapply(X = 1:trees_validation_data$L,
                                             FUN = function(l){
                                               sapply(X = 1:trees_validation_data$L,
                                                      FUN = function(lprime){
                                                        ((alpha[g]**2)*
                                                           exp(-(dists[lprime,l]**2)/
                                                                 (2*rho[g]**2)))
                                                      })
                                             })
                                    }))
                   
                   for (g in 1:trees_validation_data$G){
                     K[,,g] = K[,,g]+diag(nrow = nrow(K[,,g]),
                                          ncol = ncol(K[,,g]),
                                          x = 1e-9)
                   }
                   
                   KL = 
                     array(dim = 
                             c(trees_validation_data$L,
                               trees_validation_data$L,
                               trees_validation_data$G),
                           dimnames = 
                             list('l' = 1:trees_validation_data$L,
                                  'lprime' = 1:trees_validation_data$L,
                                  'g' = 1:trees_validation_data$G),
                           data = 
                             sapply(X = 1:trees_validation_data$G,
                                    FUN = function(g){
                                      t(chol(K[,,g]))
                                    }))
                   
                   # realized random effects
                   GP = 
                     matrix(nrow = trees_validation_data$L,
                            ncol = trees_validation_data$P,
                            data = 
                              sapply(X = 1:trees_validation_data$P,
                                     FUN = function(p){
                                       KL[,,trees_validation_data$group_id[p]] %*%
                                         zGP[,p]
                                     }))
                   
                   # linear predictors
                   logmu = 
                     sapply(X = 1:trees_validation_data$N,
                            FUN = function(i){
                              
                              beta[trees_validation_data$group_id[trees_validation_data$plot_id[i]]] +
                                plotEffects[trees_validation_data$plot_id[i]]+
                                as.numeric(GP[trees_validation_data$location_id[i],
                                              trees_validation_data$plot_id][i])
                              
                            })
                   
                   # realizations
                   Y = 
                     sapply(X = 1:trees_validation_data$N,
                            FUN = function(i){
                              rpois(n = 1,
                                    lambda = exp(logmu[i]))
                            })
                   
                   
                   results = 
                     data.frame(Y_obs = trees_validation_data$Y,
                                Y_sim = Y,
                                group_id = 
                                  trees_validation_data$group_id[trees_validation_data$plot_id],
                                plot_id = trees_validation_data$plot_id,
                                location_id = trees_validation_data$location_id,
                                logmu = logmu,
                                draw = draw) %>%
                     as_tibble() %>%
                     rowid_to_column('obs_id')
                   
                   return(results)
                 }))

saveRDS(trees_sim,
        here::here('02-data',
                   '06-results',
                   'posterior_predictions',
                   'trees_sim.rds'))

trees_sim = 
  readRDS(here::here('02-data',
                     '06-results',
                     'posterior_predictions',
                     'trees_sim.rds'))

trees_mae = 
  trees_sim %>%
  group_by(obs_id, Y_obs, group_id) %>%
  summarise(Y_sim = mean(Y_sim)) %>%
  ungroup() %>%
  mutate(ae = abs(Y_obs-Y_sim)) %>%
  pull(ae) %>%
  mean()

trees_plot_1 = 
  trees_sim %>%
  group_by(obs_id, Y_obs, group_id) %>%
  summarise(Y_sim = mean(Y_sim)) %>%
  ungroup() %>%
  ggplot(aes(x = Y_sim, y = Y_obs))+
  geom_abline(intercept=0,slope=1, color = 'red')+
  geom_smooth(method = 'lm')+
  geom_point(size = 0)+
  theme_minimal()+
  #coord_fixed()+
  facet_wrap(~group_id)+
  labs(x = 'Mean simulated tree tally',
       y = 'Observed (validation data) tree tally')

trees_plot_1

ggsave(trees_plot_1,
       filename = here::here('04-communication','figures', 'manuscript',
                             'trees_prediction_1.png'),
       height = 4, width = 6.5, units = 'in')

trees_plot_2 = 
  trees_sim %>%
  group_by(obs_id, Y_obs, group_id) %>%
  summarise() %>%
  ungroup() %>%
  mutate(source = 'obs') %>%
  group_by(Y_obs, group_id, source) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  left_join(
    trees_sim %>%
      group_by(obs_id, Y_obs, group_id) %>%
      summarise() %>%
      ungroup() %>%
      mutate(source = 'obs') %>%
      group_by(group_id, source) %>%
      summarise(total = n()) %>%
      ungroup() 
  ) %>%
  mutate(proportion = count / total) %>%
  select(group_id, Y = Y_obs, proportion, source) %>%
  bind_rows(
    trees_sim %>%
      group_by(Y_sim, group_id) %>%
      summarise(count = n()) %>%
      ungroup() %>%
      mutate(source = 'sim') %>%
      left_join(
        trees_sim %>%
          group_by(group_id) %>%
          summarise(total = n()) %>%
          ungroup() %>%
          mutate(source = 'sim')
      ) %>%
      mutate(proportion = count / total) %>%
      select(group_id, Y = Y_sim, proportion, source)
  ) %>%
  ggplot()+
  geom_col(aes(x = Y, y = proportion, fill = source),
           position = position_dodge())+
  facet_grid(group_id~., scales = 'free')+
  theme_minimal()

trees_plot_2

ggsave(trees_plot_2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'trees_prediction_2.png'),
       height = 6.5, width = 6.5, units = 'in')

trees_plot_3 = 
  trees_sim %>%
  group_by(group_id, obs_id, Y_obs) %>%
  summarise(logmu.mean = mean(logmu),
            Y_sim.mean = mean(Y_sim),
            Y_sim.05 = quantile(Y_sim, 0.05),
            Y_sim.95 = quantile(Y_sim, 0.95)) %>%
  ungroup() %>%
  ggplot()+
  geom_errorbar(aes(x = logmu.mean, ymin = Y_sim.05, ymax = Y_sim.95),
                alpha = 0.75)+
  geom_point(aes(x = logmu.mean, y = Y_obs), color = 'red')+
  theme_minimal()+
  labs(title = 'trees prediction')+
  facet_grid(group_id~., scales = 'free')

trees_plot_3

ggsave(trees_plot_3,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'trees_prediction_3.png'),
       height = 6.5, width = 6.5, units = 'in')



trees_plot_4 = 
  ggplot()+
  geom_errorbar(
    data = 
      trees_sim %>%
      group_by(group_id, Y_sim, draw) %>%
      summarise(N_obs = n()) %>%
      ungroup() %>%
      right_join(trees_sim %>%
                   expand(nesting(group_id, draw),
                          Y_sim = 0:max(.$Y_sim))) %>%
      mutate(N_obs = ifelse(is.na(N_obs), 0, N_obs)) %>%
      group_by(group_id, Y_sim) %>%
      summarise(N_obs.025 = quantile(N_obs, 0.025),
                N_obs.975 = quantile(N_obs, 0.975)) %>%
      ungroup(),
    aes(x = Y_sim, ymin = N_obs.025, ymax = N_obs.975, color = 'Simulated')
  )+
  facet_grid(group_id~., scales = 'free')+
  theme_minimal()+
  labs(x = 'Trees tally', y = 'N observations')+
  geom_col(
    data = 
      trees_sim %>%
      group_by(group_id, obs_id, Y_obs) %>%
      summarise() %>%
      ungroup() %>%
      group_by(group_id, Y_obs) %>%
      summarise(count = n()) %>%
      ungroup(),
    aes(x = Y_obs, y = count, fill = 'Observed (validation data)'),
    width = 0.5,
    alpha = 0.75
  )+
  scale_fill_manual(values = c('Observed (validation data)' = 'red'))+
  scale_color_manual(values = c('Simulated' = 'blue'))+
  theme(legend.position = 'bottom',
        legend.title = element_blank())

trees_plot_4

ggsave(trees_plot_4,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'trees_prediction_4.png'),
       height = 6.5, width = 6.5, units = 'in')

head(trees_sim)

trees_plot_5 = 
  ggplot()+
  stat_ecdf(
    data = trees_sim,
    aes(x = Y_sim, group = draw, color = 'Simulated'),
    geom = 'step', 
    pad = FALSE,
    #color = 'blue',
    lwd = 1,
    alpha = 0.1)+
  stat_ecdf(
    data = trees_sim %>% group_by(group_id, Y_obs, obs_id) %>% summarise() %>% ungroup(),
    aes(x = Y_obs, color = 'Observed (validation data)'),
    geom = 'step',
    pad = 'false',
    #color = 'red',
    lwd = 1
  )+
  theme_minimal()+
  scale_color_manual(values = c('Observed (validation data)' = 'red', 'Simulated' = 'blue'))+
  facet_grid(group_id~.)+
  labs(x = 'Tree tally', y = 'Cumulative density')+
  theme(legend.position = 'bottom',
        legend.title = element_blank())

trees_plot_5

ggsave(trees_plot_5,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'trees_prediction_5.png'),
       height = 6.5, width = 6.5, units = 'in')

head(trees_sim)

trees_plot_6 = 
  trees_sim %>%
  group_by(draw, group_id, plot_id) %>%
  summarise(tph = sum(Y_sim)/0.05) %>%
  ungroup() %>%
  ggplot()+
  stat_ecdf(
    aes(x = tph, group = draw, color = 'Simulated'),
    alpha = 0.1,
    geom = 'step'
  )+
  stat_ecdf(
    data = 
      trees_sim %>%
      group_by(group_id, plot_id, obs_id, Y_obs) %>%
      summarise() %>%
      ungroup() %>%
      group_by(group_id, plot_id) %>%
      summarise(tph = sum(Y_obs)/0.05) %>%
      ungroup(),
    aes(x = tph, color = 'Observed (validation data)'),
    geom = 'step', 
    lwd = 1
  )+
  scale_color_manual(values = c('Observed (validation data)' = 'red','Simulated' ='blue'))+
  theme_minimal()

ggsave(trees_plot_6,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'trees_prediction_6.png'),
       height = 4.5, width = 6.5, units = 'in')


rm(trees_sim)

#### saplings #####################################################################


saplings_validation_data = 
  readRDS(here::here('02-data', 
                     '05-for_analysis',
                     'saplings_validation_data.rds'))

saplings_posterior = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'saplings_fit.rds'))$draws() %>%
  as_draws_df()

set.seed(110819)
saplings_sim = 
  do.call('bind_rows',
          lapply(X = 1:4000,
                 FUN = function(draw){
                   
                   print(paste0('Working on draw', draw))
                   
                   # extract parameters for this draw
                   beta = 
                     saplings_posterior %>%
                     select(contains('beta')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   alpha = 
                     saplings_posterior %>%
                     select(contains('alpha')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   rho = 
                     saplings_posterior %>%
                     select(contains('rho')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   plotEffects = 
                     saplings_posterior %>%
                     select(contains('plot_effect')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   zGP = 
                     matrix(nrow = saplings_validation_data$L,
                            ncol = saplings_validation_data$P,
                            byrow = FALSE,
                            data = 
                              saplings_posterior %>%
                              select(contains('zGP')) %>%
                              slice(draw) %>%
                              as.data.frame() %>%
                              as.numeric())
                   
                   
                   dists = 
                     saplings_validation_data$coords %>%
                     dist() %>%
                     as.matrix()
                   
                   # construct GP effects
                   K = 
                     array(dim = 
                             c(saplings_validation_data$L,
                               saplings_validation_data$L,
                               saplings_validation_data$G),
                           dimnames = 
                             list('l' = 1:saplings_validation_data$L,
                                  'lprime' = 1:saplings_validation_data$L,
                                  'g' = 1:saplings_validation_data$G),
                           data = 
                             sapply(X = 1:saplings_validation_data$G,
                                    FUN = function(g){
                                      sapply(X = 1:saplings_validation_data$L,
                                             FUN = function(l){
                                               sapply(X = 1:saplings_validation_data$L,
                                                      FUN = function(lprime){
                                                        ((alpha[g]**2)*
                                                           exp(-(dists[lprime,l]**2)/
                                                                 (2*rho[g]**2)))
                                                      })
                                             })
                                    }))
                   
                   for (g in 1:saplings_validation_data$G){
                     K[,,g] = K[,,g]+diag(nrow = nrow(K[,,g]),
                                          ncol = ncol(K[,,g]),
                                          x = 1e-9)
                   }
                   
                   KL = 
                     array(dim = 
                             c(saplings_validation_data$L,
                               saplings_validation_data$L,
                               saplings_validation_data$G),
                           dimnames = 
                             list('l' = 1:saplings_validation_data$L,
                                  'lprime' = 1:saplings_validation_data$L,
                                  'g' = 1:saplings_validation_data$G),
                           data = 
                             sapply(X = 1:saplings_validation_data$G,
                                    FUN = function(g){
                                      t(chol(K[,,g]))
                                    }))
                   
                   # realized random effects
                   GP = 
                     matrix(nrow = saplings_validation_data$L,
                            ncol = saplings_validation_data$P,
                            data = 
                              sapply(X = 1:saplings_validation_data$P,
                                     FUN = function(p){
                                       KL[,,saplings_validation_data$group_id[p]] %*%
                                         zGP[,p]
                                     }))
                   
                   # linear predictors
                   logmu = 
                     sapply(X = 1:saplings_validation_data$N,
                            FUN = function(i){
                              
                              beta[saplings_validation_data$group_id[saplings_validation_data$plot_id[i]]] +
                                plotEffects[saplings_validation_data$plot_id[i]]+
                                as.numeric(GP[saplings_validation_data$location_id[i],
                                              saplings_validation_data$plot_id][i])
                              
                            })
                   
                   # realizations
                   Y = 
                     sapply(X = 1:saplings_validation_data$N,
                            FUN = function(i){
                              rpois(n = 1,
                                    lambda = exp(logmu[i]))
                            })
                   
                   
                   results = 
                     data.frame(Y_obs = saplings_validation_data$Y,
                                Y_sim = Y,
                                group_id = 
                                  saplings_validation_data$group_id[saplings_validation_data$plot_id],
                                plot_id = saplings_validation_data$plot_id,
                                location_id = saplings_validation_data$location_id,
                                logmu = logmu,
                                draw = draw) %>%
                     as_tibble() %>%
                     rowid_to_column('obs_id')
                   
                   return(results)
                 }))

saveRDS(saplings_sim,
        here::here('02-data',
                   '06-results',
                   'posterior_predictions',
                   'saplings_sim.rds'))

saplings_sim = readRDS(here::here('02-data',
                                  '06-results',
                                  'posterior_predictions',
                                  'saplings_sim.rds'))

saplings_mae = 
  saplings_sim %>%
  group_by(obs_id, Y_obs, group_id) %>%
  summarise(Y_sim = mean(Y_sim)) %>%
  ungroup() %>%
  mutate(ae = abs(Y_obs-Y_sim)) %>%
  pull(ae) %>%
  mean()

saplings_plot_1 = 
  saplings_sim %>%
  group_by(obs_id, Y_obs, group_id) %>%
  summarise(Y_sim = mean(Y_sim)) %>%
  ungroup() %>%
  ggplot(aes(x = Y_sim, y = Y_obs))+
  geom_abline(intercept=0,slope=1, color = 'red')+
  geom_smooth(method = 'lm')+
  geom_point(size = 0)+
  theme_minimal()+
  #coord_fixed()+
  facet_wrap(~group_id)+
  labs(x = 'Mean simulated sapling tally',
       y = 'Observed (validation data) sapling tally')

saplings_plot_1

ggsave(saplings_plot_1,
       filename = here::here('04-communication','figures', 'manuscript',
                             'saplings_prediction_1.png'),
       height = 4, width = 6.5, units = 'in')

saplings_plot_2 = 
  saplings_sim %>%
  group_by(obs_id, Y_obs, group_id) %>%
  summarise() %>%
  ungroup() %>%
  mutate(source = 'obs') %>%
  group_by(Y_obs, group_id, source) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  left_join(
    saplings_sim %>%
      group_by(obs_id, Y_obs, group_id) %>%
      summarise() %>%
      ungroup() %>%
      mutate(source = 'obs') %>%
      group_by(group_id, source) %>%
      summarise(total = n()) %>%
      ungroup() 
  ) %>%
  mutate(proportion = count / total) %>%
  select(group_id, Y = Y_obs, proportion, source) %>%
  bind_rows(
    saplings_sim %>%
      group_by(Y_sim, group_id) %>%
      summarise(count = n()) %>%
      ungroup() %>%
      mutate(source = 'sim') %>%
      left_join(
        saplings_sim %>%
          group_by(group_id) %>%
          summarise(total = n()) %>%
          ungroup() %>%
          mutate(source = 'sim')
      ) %>%
      mutate(proportion = count / total) %>%
      select(group_id, Y = Y_sim, proportion, source)
  ) %>%
  ggplot()+
  geom_col(aes(x = Y, y = proportion, fill = source),
           position = position_dodge())+
  facet_grid(group_id~., scales = 'free')+
  theme_minimal()+
  labs(title = 'saplings prediction')

saplings_plot_2

ggsave(saplings_plot_2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'saplings_prediction_2.png'),
       height = 6.5, width = 6.5, units = 'in')

saplings_plot_3 = 
  saplings_sim %>%
  group_by(group_id, obs_id, Y_obs) %>%
  summarise(logmu.mean = mean(logmu),
            Y_sim.mean = mean(Y_sim),
            Y_sim.05 = quantile(Y_sim, 0.05),
            Y_sim.95 = quantile(Y_sim, 0.95)) %>%
  ungroup() %>%
  ggplot()+
  geom_errorbar(aes(x = logmu.mean, ymin = Y_sim.05, ymax = Y_sim.95),
                alpha = 0.75)+
  geom_point(aes(x = logmu.mean, y = Y_obs), color = 'red')+
  theme_minimal()+
  labs(title = 'saplings prediction')+
  facet_grid(group_id~., scales = 'free')

saplings_plot_3

ggsave(saplings_plot_3,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'saplings_prediction_3.png'),
       height = 6.5, width = 6.5, units = 'in')


saplings_plot_4 = 
  ggplot()+
  geom_errorbar(
    data = 
      saplings_sim %>%
      group_by(group_id, Y_sim, draw) %>%
      summarise(N_obs = n()) %>%
      ungroup() %>%
      right_join(saplings_sim %>%
                   expand(nesting(group_id, draw),
                          Y_sim = 0:max(.$Y_sim))) %>%
      mutate(N_obs = ifelse(is.na(N_obs), 0, N_obs)) %>%
      group_by(group_id, Y_sim) %>%
      summarise(N_obs.025 = quantile(N_obs, 0.025),
                N_obs.975 = quantile(N_obs, 0.975)) %>%
      ungroup(),
    aes(x = Y_sim, ymin = N_obs.025, ymax = N_obs.975, color = 'Simulated')
  )+
  facet_grid(group_id~., scales = 'free')+
  theme_minimal()+
  labs(x = 'Saplings tally', y = 'N observations')+
  geom_col(
    data = 
      saplings_sim %>%
      group_by(group_id, obs_id, Y_obs) %>%
      summarise() %>%
      ungroup() %>%
      group_by(group_id, Y_obs) %>%
      summarise(count = n()) %>%
      ungroup(),
    aes(x = Y_obs, y = count, fill = 'Observed (validation data)'),
    width = 0.5,
    alpha = 0.75
  )+
  scale_fill_manual(values = c('Observed (validation data)' = 'red'))+
  scale_color_manual(values = c('Simulated' = 'blue'))+
  theme(legend.position = 'bottom',
        legend.title = element_blank())

saplings_plot_4


saplings_plot_4
ggsave(saplings_plot_4,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'saplings_prediction_4.png'),
       height = 6.5, width = 6.5, units = 'in')

head(saplings_sim)

saplings_plot_5 = 
  ggplot()+
  stat_ecdf(
    data = saplings_sim,
    aes(x = Y_sim, group = draw, color = 'Simulated'),
    geom = 'step', 
    pad = FALSE,
    #color = 'blue',
    lwd = 1,
    alpha = 0.1)+
  stat_ecdf(
    data = saplings_sim %>% group_by(group_id, Y_obs, obs_id) %>% summarise() %>% ungroup(),
    aes(x = Y_obs, color = 'Observed (validation data)'),
    geom = 'step',
    pad = 'false',
    #color = 'red',
    lwd = 1
  )+
  theme_minimal()+
  scale_color_manual(values = c('Observed (validation data)' = 'red', 'Simulated' = 'blue'))+
  facet_grid(group_id~.)+
  labs(x = 'Sapling tally', y = 'Cumulative density')+
  theme(legend.position = 'bottom',
        legend.title = element_blank())

saplings_plot_5

ggsave(saplings_plot_5,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'saplings_prediction_5.png'),
       height = 6.5, width = 6.5, units = 'in')

head(saplings_sim)

saplings_plot_6 = 
  saplings_sim %>%
  group_by(draw, group_id, plot_id) %>%
  summarise(tph = sum(Y_sim)/0.05) %>%
  ungroup() %>%
  ggplot()+
  stat_ecdf(
    aes(x = tph, group = draw, color = 'Simulated'),
    alpha = 0.1,
    geom = 'step'
  )+
  stat_ecdf(
    data = 
      saplings_sim %>%
      group_by(group_id, plot_id, obs_id, Y_obs) %>%
      summarise() %>%
      ungroup() %>%
      group_by(group_id, plot_id) %>%
      summarise(tph = sum(Y_obs)/0.05) %>%
      ungroup(),
    aes(x = tph, color = 'Observed (validation data)'),
    geom = 'step', 
    lwd = 1
  )+
  scale_color_manual(values = c('Observed (validation data)' = 'red','Simulated' ='blue'))+
  theme_minimal()

ggsave(saplings_plot_6,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'saplings_prediction_6.png'),
       height = 4.5, width = 6.5, units = 'in')


rm(saplings_sim)

# cwdtallies mae 0.263
0.263/mean(cwdtallies_validation_data$Y)

# saplings mae 0.112
0.112/mean(saplings_validation_data$Y)
# trees mae 0.096
0.096/mean(trees_validation_data$Y)
