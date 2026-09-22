library(here)
library(tidyverse)
library(posterior)
library(cmdstanr)

#### setup surface fuels biomass ###############################################

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


# first aggregate to plot level by summing individual trees
# note that we lump live and dead trees together for the purpose of estimating 
# fuels, because both contribute (or have contributed) to the existing fuel bed
overstory_composition = 
  
  
  readRDS(here::here('02-data',
                     '04-geolocated',
                     'seki_trees.rds')) %>%
  
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
  
  group_by(mort, plot_id, spp) %>%
  summarise(ba_m2ha = sum(ba_m2ha)) %>%
  ungroup() %>%
  
  # fill in 0s for species which werent present on a plot
  complete(nesting(mort, plot_id), spp) %>%
  mutate(ba_m2ha = ifelse(is.na(ba_m2ha), 0, ba_m2ha)) %>%
  
  # get the total BA/ha on each plot
  left_join(x = .,
            y = 
              group_by(., mort, plot_id) %>% 
              summarise(total_ba_m2ha = sum(ba_m2ha)) %>%
              ungroup()) %>%
  
  # get the proportion of plot total basal area occoupied by each species
  mutate(pBA = ba_m2ha / total_ba_m2ha) %>%
  
  # aggregate to mortality-class average composition
  group_by(mort, spp) %>%
  summarise(pBA = mean(pBA)) %>%
  ungroup()


overstory_composition # very interesting the different mortality classes 
# have very different spp composition

# get BA-weighted average coefficient for fuel load (kg/m2) as a function of 
# litterduff_depth (cm)
litter_coeffs = 
  overstory_composition %>%
  left_join(litterduff_coeffs) %>%
  mutate(weighted = pBA * litter_coeff) %>%
  group_by(mort) %>%
  summarise(weighted_coeff = sum(weighted)) %>%
  ungroup()

duff_coeffs = 
  overstory_composition %>%
  left_join(litterduff_coeffs) %>%
  mutate(weighted = pBA * duff_coeff) %>%
  group_by(mort) %>%
  summarise(weighted_coeff = sum(weighted)) %>%
  ungroup()

litter_coeffs

# get BA-weighted SEC (secant of acute angle) for each timelag class
SEC = 
  overstory_composition %>%
  left_join(SEC) %>%
  select(mort, spp, pBA, x1h, x10h, x100h, x1000s = x1000h, x1000r = x1000h) %>%
  pivot_longer(c(x1h, x10h, x100h, x1000s, x1000r),
               names_to = 'timelag_class',
               values_to = 'sec',
               names_prefix = 'x') %>%
  mutate(weighted = sec*pBA) %>%
  group_by(mort, timelag_class) %>%
  summarise(weighted_sec = sum(weighted)) %>%
  ungroup()


# get the BA-weighted average specific gravity for each timelag class
SG = 
  overstory_composition %>%
  left_join(SG) %>%
  select(mort, spp, pBA, x1h, x10h, x100h, x1000s) %>%
  pivot_longer(c(x1h, x10h, x100h, x1000s),
               names_to = 'timelag_class',
               values_to = 'sg',
               names_prefix = 'x') %>%
  mutate(weighted = sg*pBA) %>%
  group_by(mort, timelag_class) %>%
  summarise(weighted_sg = sum(weighted)) %>%
  ungroup() %>%
  
  # and add an entry for 1000h rotten, which doesn't vary by species
  rbind(.,
        data.frame(mort = c('low', 'mid', 'high'),
                   timelag_class = c('1000r', '1000r', '1000r'),
                   weighted_sg = rep(sg_1000r, times = 3)) %>%
          mutate(mort = as.character(mort),
                 timelag_class = as.character(timelag_class)))

# get BA-weighted average QMD (cm2) for each FWD timelag class
QMDcm = 
  overstory_composition %>%
  left_join(QMDcm) %>%
  select(mort, spp, pBA, x1h, x10h, x100h) %>%
  pivot_longer(c(x1h, x10h, x100h),
               names_to = 'timelag_class',
               values_to = 'qmd_cm',
               names_prefix = 'x') %>%
  mutate(weighted = qmd_cm*pBA) %>%
  group_by(mort, timelag_class) %>%
  summarise(weighted_qmd = sum(weighted)) %>%
  ungroup()


#### litter ####################################################################
# load the observed data

# load the posterior distribution
litter_posterior = 
  readRDS(here::here('02-data',
                                    '06-results',
                                    'real_fits',
                                    'litter_fit.rds'))$draws() %>%
  as_draws_df()

litter_training_data = 
  readRDS(here::here('02-data', 
                     '05-for_analysis',
                     'litter_training_data.rds'))

# simulate 4k datasets, including new draws of random effect realizations
litter_sim = 
  do.call('bind_rows',
          lapply(X = 1:4000,
                 FUN = function(draw){
                   
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
                   
                   sigmaPlot = 
                     litter_posterior %>%
                     select(contains('sigmaPlot')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   kappa = 
                     litter_posterior %>%
                     select(contains('kappa')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   dists = 
                     litter_training_data$coords %>%
                     dist() %>%
                     as.matrix()
                   
                   # draw new random effect realizations
                   zPlot = rnorm(n = litter_training_data$P,
                                 mean = 0,
                                 sd = 1)
                   
                   zGP = 
                     matrix(nrow = litter_training_data$L,
                            ncol = litter_training_data$P,
                            byrow = FALSE,
                            data = 
                              rnorm(n = litter_training_data$L*litter_training_data$P,
                                    mean = 0,
                                    sd = 1))
                   
                   # construct plot effects
                   plotEffects = zPlot * sigmaPlot[litter_training_data$group_id]
                     
                   # construct GP effects
                   K = 
                     array(dim = 
                             c(litter_training_data$L,
                               litter_training_data$L,
                               litter_training_data$G),
                           dimnames = 
                             list('l' = 1:litter_training_data$L,
                                  'lprime' = 1:litter_training_data$L,
                                  'g' = 1:litter_training_data$G),
                           data = 
                             sapply(X = 1:litter_training_data$G,
                                    FUN = function(g){
                                      sapply(X = 1:litter_training_data$L,
                                             FUN = function(l){
                                               sapply(X = 1:litter_training_data$L,
                                                      FUN = function(lprime){
                                                        ((alpha[g]**2)*
                                                          exp(-(dists[lprime,l]**2)/
                                                                (2*rho[g]**2)))
                                                      })
                                             })
                                    }))
                   
                   for (g in 1:litter_training_data$G){
                     K[,,g] = K[,,g]+diag(nrow = nrow(K[,,g]),
                                          ncol = ncol(K[,,g]),
                                          x = 1e-9)
                   }
                   
                   KL = 
                     array(dim = 
                             c(litter_training_data$L,
                               litter_training_data$L,
                               litter_training_data$G),
                           dimnames = 
                             list('l' = 1:litter_training_data$L,
                                  'lprime' = 1:litter_training_data$L,
                                  'g' = 1:litter_training_data$G),
                           data = 
                             sapply(X = 1:litter_training_data$G,
                                    FUN = function(g){
                                      t(chol(K[,,g]))
                                    }))
                 
                   # realized random effects
                   GP = 
                     matrix(nrow = litter_training_data$L,
                            ncol = litter_training_data$P,
                            data = 
                              sapply(X = 1:litter_training_data$P,
                                     FUN = function(p){
                                       KL[,,litter_training_data$group_id[p]] %*%
                                         zGP[,p]
                                     }))
                   
                   # linear predictors
                   logmu = 
                     sapply(X = 1:litter_training_data$N,
                            FUN = function(i){
                              
                               beta[litter_training_data$group_id[litter_training_data$plot_id[i]]] +
                                plotEffects[litter_training_data$plot_id[i]]+
                                as.numeric(GP[litter_training_data$location_id[i],
                                              litter_training_data$plot_id][i])
                              
                            })
                    
                   # realizations
                   Y = 
                     sapply(X = 1:litter_training_data$N,
                            FUN = function(i){
                              rnbinom(n = 1,
                                      mu = exp(logmu[i]),
                                      size = 
                                        kappa[litter_training_data$group_id[litter_training_data$plot_id[i]]])
                            })
                   
                   
                   results = 
                     data.frame(Y = Y,
                                group_id = 
                                  litter_training_data$group_id[litter_training_data$plot_id],
                                plot_id = litter_training_data$plot_id,
                                location_id = litter_training_data$location_id,
                                logmu = logmu,
                                draw = draw,
                                source = 'sim') %>%
                     as_tibble() %>%
                     rowid_to_column('obs_id')
                   
                   return(results)
                 }))



# combine observed and simulated datasets
litter_sim = 
  litter_sim %>%
  bind_rows(
    data.frame(Y = litter_training_data$Y,
               group_id = litter_training_data$group_id[litter_training_data$plot_id],
               plot_id = litter_training_data$plot_id,
               location_id = litter_training_data$location_id,
               logmu = NA,
               draw = 4001,
               source = 'observed')
  )


# convert to fuel loads
litter_sim = 
  litter_sim %>%
  left_join(
    data.frame(group_id = c(1, 2, 3),
               mort = c('low', 'mid', 'high'))
  ) %>%
  left_join(litter_coeffs) %>%
  # coefficient goes from depth in cm to biomass in kg/m2, so multiply by 10 
  # to get biomass in Mg/ha
  mutate(litter_mgha = Y * weighted_coeff * 10)

# aggregate to plot level
head(litter_sim)

litter_sim_plots = 
  litter_sim %>%
  group_by(group_id, plot_id, draw, source) %>%
  summarise(litter_mgha = mean(litter_mgha)) %>%
  ungroup() 

litter_sim_plots %>%
  ggplot(aes(x = source, y = litter_mgha))+
  geom_boxplot(
    data = litter_sim_plots %>%
      filter(source == 'sim'),
    width = 0.25,
    position = position_nudge(x = -0.35))+
  geom_jitter(
    data = litter_sim_plots %>%
      filter(source == 'observed'),
    height = 0, width = 0.15)+
  geom_violin(
    data = litter_sim_plots %>%
      filter(source == 'sim'),
    width = 0.25,
    position = position_nudge(x = 0.35))+
  facet_grid(group_id ~.)+
  scale_y_log10()+
  theme_minimal()

litter_sim_plots %>%
  filter(source=='sim') %>%
  pull(litter_mgha) %>%
  summary()

ggplot(data = litter_sim_plots,
       aes(x = litter_mgha, color = source, group = draw))+
  stat_ecdf(geom = 'step',
            lwd = 1,
            pad = TRUE,
            alpha = 0.1)+
  scale_x_log10()+
  theme_minimal()

plotlevel_litter_cdf = 
  ggplot()+
  stat_ecdf(
    data = litter_sim_plots %>%
      filter(source == 'sim'),
    aes(x = litter_mgha, group = draw),
    geom = 'step',
    alpha = 0.1,
    pad = FALSE,
    color = 'blue')+
  stat_ecdf(
    data = litter_sim_plots %>%
      filter(source=='observed'),
    aes(x = litter_mgha),
    geom = 'step',
    lwd = 1,
    pad = FALSE,
    color = 'red'
  )+
  scale_x_log10()+
  theme_minimal()

plotlevel_litter_cdf

ggsave(plotlevel_litter_cdf,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'litter_plotlevel.png'),
       height = 4.5, width = 6.5, units = 'in')


head(litter_sim_plots)

# all of the observed data are plausible assuming the model is true, but the 
# model predicts (rare) behaviors that seem unlikely

# plot results

rm(list = ls('litter_sim'))

#### duff ######################################################################

# load the observed data

# load the posterior distribution
duff_posterior = 
  readRDS(here::here('02-data',
                                    '06-results',
                                    'real_fits',
                                    'duff_fit.rds'))$draws() %>%
  as_draws_df()

duff_training_data = 
  readRDS(here::here('02-data', 
                     '05-for_analysis',
                     'duff_training_data.rds'))

# simulate 4k datasets, including new draws of random effect realizations
duff_sim = 
  do.call('bind_rows',
          lapply(X = 1:4000,
                 FUN = function(draw){
                   
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
                   
                   sigmaPlot = 
                     duff_posterior %>%
                     select(contains('sigmaPlot')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   kappa = 
                     duff_posterior %>%
                     select(contains('kappa')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   dists = 
                     duff_training_data$coords %>%
                     dist() %>%
                     as.matrix()
                   
                   # draw new random effect realizations
                   zPlot = rnorm(n = duff_training_data$P,
                                 mean = 0,
                                 sd = 1)
                   
                   zGP = 
                     matrix(nrow = duff_training_data$L,
                            ncol = duff_training_data$P,
                            byrow = FALSE,
                            data = 
                              rnorm(n = duff_training_data$L*duff_training_data$P,
                                    mean = 0,
                                    sd = 1))
                   
                   # construct plot effects
                   plotEffects = zPlot * sigmaPlot[duff_training_data$group_id]
                     
                   # construct GP effects
                   K = 
                     array(dim = 
                             c(duff_training_data$L,
                               duff_training_data$L,
                               duff_training_data$G),
                           dimnames = 
                             list('l' = 1:duff_training_data$L,
                                  'lprime' = 1:duff_training_data$L,
                                  'g' = 1:duff_training_data$G),
                           data = 
                             sapply(X = 1:duff_training_data$G,
                                    FUN = function(g){
                                      sapply(X = 1:duff_training_data$L,
                                             FUN = function(l){
                                               sapply(X = 1:duff_training_data$L,
                                                      FUN = function(lprime){
                                                        ((alpha[g]**2)*
                                                          exp(-(dists[lprime,l]**2)/
                                                                (2*rho[g]**2)))
                                                      })
                                             })
                                    }))
                   
                   for (g in 1:duff_training_data$G){
                     K[,,g] = K[,,g]+diag(nrow = nrow(K[,,g]),
                                          ncol = ncol(K[,,g]),
                                          x = 1e-9)
                   }
                   
                   KL = 
                     array(dim = 
                             c(duff_training_data$L,
                               duff_training_data$L,
                               duff_training_data$G),
                           dimnames = 
                             list('l' = 1:duff_training_data$L,
                                  'lprime' = 1:duff_training_data$L,
                                  'g' = 1:duff_training_data$G),
                           data = 
                             sapply(X = 1:duff_training_data$G,
                                    FUN = function(g){
                                      t(chol(K[,,g]))
                                    }))
                 
                   # realized random effects
                   GP = 
                     matrix(nrow = duff_training_data$L,
                            ncol = duff_training_data$P,
                            data = 
                              sapply(X = 1:duff_training_data$P,
                                     FUN = function(p){
                                       KL[,,duff_training_data$group_id[p]] %*%
                                         zGP[,p]
                                     }))
                   
                   # linear predictors
                   logmu = 
                     sapply(X = 1:duff_training_data$N,
                            FUN = function(i){
                              
                               beta[duff_training_data$group_id[duff_training_data$plot_id[i]]] +
                                plotEffects[duff_training_data$plot_id[i]]+
                                as.numeric(GP[duff_training_data$location_id[i],
                                              duff_training_data$plot_id][i])
                              
                            })
                    
                   # realizations
                   Y = 
                     sapply(X = 1:duff_training_data$N,
                            FUN = function(i){
                              rnbinom(n = 1,
                                      mu = exp(logmu[i]),
                                      size = 
                                        kappa[duff_training_data$group_id[duff_training_data$plot_id[i]]])
                            })
                   
                   
                   results = 
                     data.frame(Y = Y,
                                group_id = 
                                  duff_training_data$group_id[duff_training_data$plot_id],
                                plot_id = duff_training_data$plot_id,
                                location_id = duff_training_data$location_id,
                                logmu = logmu,
                                draw = draw,
                                source = 'sim') %>%
                     as_tibble() %>%
                     rowid_to_column('obs_id')
                   
                   return(results)
                 }))



# combine observed and simulated datasets
duff_sim = 
  duff_sim %>%
  bind_rows(
    data.frame(Y = duff_training_data$Y,
               group_id = duff_training_data$group_id[duff_training_data$plot_id],
               plot_id = duff_training_data$plot_id,
               location_id = duff_training_data$location_id,
               logmu = NA,
               draw = 4001,
               source = 'observed')
  )


# convert to fuel loads
duff_sim = 
  duff_sim %>%
  left_join(
    data.frame(group_id = c(1, 2, 3),
               mort = c('low', 'mid', 'high'))
  ) %>%
  left_join(duff_coeffs) %>%
  # coefficient goes from depth in cm to biomass in kg/m2, so multiply by 10 
  # to get biomass in Mg/ha
  mutate(duff_mgha = Y * weighted_coeff * 10)

# aggregate to plot level
head(duff_sim)

 duff_sim_plots = 
  duff_sim %>%
  group_by(group_id, plot_id, draw, source) %>%
  summarise(duff_mgha = mean(duff_mgha)) %>%
  ungroup() 

duff_sim_plots %>%
  ggplot(aes(x = source, y = duff_mgha))+
  geom_boxplot(
    data = duff_sim_plots %>%
      filter(source == 'sim'),
    width = 0.25,
    position = position_nudge(x = -0.35))+
  geom_jitter(
    data = duff_sim_plots %>%
      filter(source == 'observed'),
    height = 0, width = 0.15)+
  geom_violin(
    data = duff_sim_plots %>%
      filter(source == 'sim'),
    width = 0.25,
    position = position_nudge(x = 0.35))+
  facet_grid(group_id ~.)+
  scale_y_log10()+
  theme_minimal()

duff_sim_plots %>%
  filter(source=='sim') %>%
  pull(duff_mgha) %>%
  summary()

ggplot(data = duff_sim_plots,
       aes(x = duff_mgha, color = source, group = draw))+
  stat_ecdf(geom = 'step',
            lwd = 1,
            pad = TRUE,
            alpha = 0.1)+
  scale_x_log10()+
  theme_minimal()

plotlevel_duff_cdf = 
  ggplot()+
  stat_ecdf(
    data = duff_sim_plots %>%
      filter(source == 'sim'),
    aes(x = duff_mgha, group = draw),
    geom = 'step',
    alpha = 0.1,
    pad = FALSE,
    color = 'blue')+
  stat_ecdf(
    data = duff_sim_plots %>%
      filter(source=='observed'),
    aes(x = duff_mgha),
    geom = 'step',
    lwd = 1,
    pad = FALSE,
    color = 'red'
  )+
  scale_x_log10()+
  theme_minimal()

plotlevel_duff_cdf

head(duff_sim_plots)



ggsave(plotlevel_duff_cdf,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'duff_plotlevel.png'),
       height = 4.5, width = 6.5, units = 'in')


head(duff_sim_plots)

