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
set.seed(110819)
litter_sim = 
  do.call('bind_rows',
          lapply(X = 1:4000,
                 FUN = function(draw){
                   
                   print(paste0('Working on draw', draw))
                   
                   
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
                                              litter_training_data$plot_id[i]])
                              
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
  theme_minimal()+
  scale_color_manual(values = c('observed' = 'red', 'simulated' = 'blue'))+
  labs(x = 'Plot mean litter load (Mg / ha)', y = 'Cumulative density function')

plotlevel_litter_cdf

ggsave(plotlevel_litter_cdf,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'litter_plotlevel.png'),
       height = 6.5, width = 6.5, units = 'in')


head(litter_sim_plots)

plotlevel_litter_cdf2 = 
  litter_sim_plots %>%
  filter(source=='sim') %>%
  group_by(source, draw) %>%
  mutate(plot_rank = row_number(litter_mgha)) %>%
  ungroup() %>%
  group_by(source, plot_rank) %>%
  summarise(litter_mgha.01 = quantile(litter_mgha, 0.01),
            litter_mgha.10 = quantile(litter_mgha, 0.1),
            litter_mgha.25 = quantile(litter_mgha, 0.25),
            litter_mgha.75 = quantile(litter_mgha, 0.75),
            litter_mgha.90 = quantile(litter_mgha, 0.90),
            litter_mgha.99 = quantile(litter_mgha, 0.99)) %>%
  ungroup() %>%
  mutate(y_min = (plot_rank/21)-(1/21),
         y_max = (plot_rank/21)) %>%
  ggplot()+
  geom_rect(aes(xmin = litter_mgha.01, xmax = litter_mgha.99, ymin = y_min, ymax = y_max,
                fill = 'simulated'),
            alpha = 0.25)+
  geom_rect(aes(xmin = litter_mgha.10, xmax = litter_mgha.90, ymin = y_min, ymax = y_max,
                fill = 'simulated'),
            alpha = 0.5)+
  geom_rect(aes(xmin = litter_mgha.25, xmax = litter_mgha.75, ymin = y_min, ymax = y_max,
                fill = 'simulated'),
            alpha = 1)+
  stat_ecdf(
    data = litter_sim_plots %>%
      filter(source=='observed'),
    aes(x = litter_mgha, color = 'observed'),
    geom = 'step',
    lwd = 1)+
  scale_x_log10()+
  theme_minimal()+
  scale_fill_manual(values = c('simulated' = 'blue'))+
  scale_color_manual(values = c('observed' = 'red'))+
  labs(x = 'Plot mean litter load (Mg / ha)', y = 'Cumulative density function')+
  theme(legend.title = element_blank())

ggsave(plotlevel_litter_cdf2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'litter_plotlevel_2.png'),
       height = 6.5, width = 6.5, units = 'in')


# all of the observed data are plausible assuming the model is true, but the 
# model predicts (rare) behaviors that seem unlikely

# plot results

rm(litter_sim)

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
set.seed(110819)
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
                                              duff_training_data$plot_id[i]])
                              
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
  theme_minimal()+
  scale_color_manual(values = c('observed' = 'red', 'simulated' = 'blue'))+
  labs(x = 'Plot mean duff load (Mg / ha)', y = 'Cumulative density function')

plotlevel_duff_cdf

ggsave(plotlevel_duff_cdf,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'duff_plotlevel.png'),
       height = 6.5, width = 6.5, units = 'in')


head(duff_sim_plots)

plotlevel_duff_cdf2 = 
  duff_sim_plots %>%
  filter(source=='sim') %>%
  group_by(source, draw) %>%
  mutate(plot_rank = row_number(duff_mgha)) %>%
  ungroup() %>%
  group_by(source, plot_rank) %>%
  summarise(duff_mgha.01 = quantile(duff_mgha, 0.01),
            duff_mgha.10 = quantile(duff_mgha, 0.1),
            duff_mgha.25 = quantile(duff_mgha, 0.25),
            duff_mgha.75 = quantile(duff_mgha, 0.75),
            duff_mgha.90 = quantile(duff_mgha, 0.90),
            duff_mgha.99 = quantile(duff_mgha, 0.99)) %>%
  ungroup() %>%
  mutate(y_min = (plot_rank/21)-(1/21),
         y_max = (plot_rank/21)) %>%
  ggplot()+
  geom_rect(aes(xmin = duff_mgha.01, xmax = duff_mgha.99, ymin = y_min, ymax = y_max,
                fill = 'simulated'),
            alpha = 0.25)+
  geom_rect(aes(xmin = duff_mgha.10, xmax = duff_mgha.90, ymin = y_min, ymax = y_max,
                fill = 'simulated'),
            alpha = 0.5)+
  geom_rect(aes(xmin = duff_mgha.25, xmax = duff_mgha.75, ymin = y_min, ymax = y_max,
                fill = 'simulated'),
            alpha = 1)+
  stat_ecdf(
    data = duff_sim_plots %>%
      filter(source=='observed'),
    aes(x = duff_mgha, color = 'observed'),
    geom = 'step',
    lwd = 1)+
  scale_x_log10()+
  theme_minimal()+
  scale_fill_manual(values = c('simulated' = 'blue'))+
  scale_color_manual(values = c('observed' = 'red'))+
  labs(x = 'Plot mean duff load (Mg / ha)', y = 'Cumulative density function')+
  theme(legend.title = element_blank())

plotlevel_duff_cdf2

ggsave(plotlevel_duff_cdf2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'duff_plotlevel_2.png'),
       height = 6.5, width = 6.5, units = 'in')

rm(duff_sim)

#### FWD #######################################################################

# first, simulate tallies
fwdtallies_training_data = 
  readRDS(here::here('02-data', 
                     '05-for_analysis',
                     'fwdtallies_training_data.rds'))

fwddiams_training_data = 
  readRDS(here::here('02-data',
                     '05-for_analysis',
                     'fwddiams_training_data.rds'))

fwdtallies_posterior = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'fwdtallies_fit.rds'))$draws() %>%
  as_draws_df()

fwddiams_posterior = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'fwddiams_fit.rds'))$draws() %>%
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
                   
                   sigmaPlot = 
                     fwdtallies_posterior %>%
                     select(contains('sigmaPlot')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   zPlot = 
                     rnorm(n = fwdtallies_training_data$P,
                           mean = 0,
                           sd = 1)
                   plotEffects = 
                     zPlot * sigmaPlot[fwdtallies_training_data$group_id]
                   
                   zGP = 
                     matrix(nrow = fwdtallies_training_data$L,
                            ncol = fwdtallies_training_data$P,
                            byrow = FALSE,
                            data = 
                              rnorm(n = fwdtallies_training_data$L*fwdtallies_training_data$P,
                                    mean = 0,
                                    sd = 1)) 
                   
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
                     fwdtallies_training_data$coords %>%
                     dist() %>%
                     as.matrix()
                   
                   # construct GP effects
                   K = 
                     array(dim = 
                             c(fwdtallies_training_data$L,
                               fwdtallies_training_data$L,
                               fwdtallies_training_data$G),
                           dimnames = 
                             list('l' = 1:fwdtallies_training_data$L,
                                  'lprime' = 1:fwdtallies_training_data$L,
                                  'g' = 1:fwdtallies_training_data$G),
                           data = 
                             sapply(X = 1:fwdtallies_training_data$G,
                                    FUN = function(g){
                                      sapply(X = 1:fwdtallies_training_data$L,
                                             FUN = function(l){
                                               sapply(X = 1:fwdtallies_training_data$L,
                                                      FUN = function(lprime){
                                                        ((alpha[g]**2)*
                                                          exp(-(dists[lprime,l]**2)/
                                                                (2*rho[g]**2)))
                                                      })
                                             })
                                    }))
                   
                   for (g in 1:fwdtallies_training_data$G){
                     K[,,g] = K[,,g]+diag(nrow = nrow(K[,,g]),
                                          ncol = ncol(K[,,g]),
                                          x = tau[g])
                   }
                   
                   KL = 
                     array(dim = 
                             c(fwdtallies_training_data$L,
                               fwdtallies_training_data$L,
                               fwdtallies_training_data$G),
                           dimnames = 
                             list('l' = 1:fwdtallies_training_data$L,
                                  'lprime' = 1:fwdtallies_training_data$L,
                                  'g' = 1:fwdtallies_training_data$G),
                           data = 
                             sapply(X = 1:fwdtallies_training_data$G,
                                    FUN = function(g){
                                      t(chol(K[,,g]))
                                    }))
                 
                   # realized random effects
                   GP = 
                     matrix(nrow = fwdtallies_training_data$L,
                            ncol = fwdtallies_training_data$P,
                            data = 
                              sapply(X = 1:fwdtallies_training_data$P,
                                     FUN = function(p){
                                       KL[,,fwdtallies_training_data$group_id[p]] %*%
                                         zGP[,p]
                                     }))
                   
                   # linear predictors
                   logmu = 
                     sapply(X = 1:fwdtallies_training_data$N,
                            FUN = function(i){
                              
                               beta[fwdtallies_training_data$group_id[fwdtallies_training_data$plot_id[i]]] +
                                plotEffects[fwdtallies_training_data$plot_id[i]]+
                                as.numeric(GP[fwdtallies_training_data$location_id[i],
                                              fwdtallies_training_data$plot_id[i]])
                              
                            })
                    
                   # realizations
                   Y = 
                     sapply(X = 1:fwdtallies_training_data$N,
                            FUN = function(i){
                              rnbinom(n = 1,
                                      mu = exp(logmu[i]),
                                      size = 
                                        kappa[fwdtallies_training_data$group_id[fwdtallies_training_data$plot_id[i]]])
                            })
                   
                   
                   results = 
                     data.frame(Y_obs = fwdtallies_training_data$Y,
                                Y_sim = Y,
                                group_id = 
                                  fwdtallies_training_data$group_id[fwdtallies_training_data$plot_id],
                                plot_id = fwdtallies_training_data$plot_id,
                                location_id = fwdtallies_training_data$location_id,
                                logmu = logmu,
                                draw = draw) %>%
                     as_tibble() %>%
                     rowid_to_column('obs_id')
                   
                   return(results)
                 }))

set.seed(110819)
# next, reformat the simulated tallies to be input data for the diameters model, 
# and simulate a diameter for each particle
fwddiams_simulatedtallies_data = 
  do.call('rbind',
          lapply(X = 1:4000,#sample(1:4000, 300),#1:100,
             FUN = function(d){
               
               print(paste0('Simulating draw ', d))
               
               # reformat the simulated tallies to be input data for the diameters model
               
               sim_tally_data = fwdtallies_sim %>% filter(draw == d)
               
               #sim_tally_long = 
              #   do.call('bind_rows',
               #          lapply(X = 1:nrow(sim_tally_data),
                #                FUN = function(i){
                 #                 if (sim_tally_data$Y_sim[i]==0){
                #                    return(NULL)
                 #                 }
                #                  longwise = 
                 #                   sim_tally_data[i,] %>%
                 #                   select(obs_id, group_id, plot_id, location_id, draw) %>%
                  #                  right_join(data.frame(particle_id = 
                   #                                         1:sim_tally_data$Y_sim[i]),
                    #                           by = character())
                     #             return(longwise)
                      #          }))
               
               sim_tally_long = 
                 sim_tally_data %>%
                 filter(Y_sim > 0) %>%
                 select(obs_id, group_id, plot_id, location_id, draw, Y_sim) %>%
                 group_by(obs_id, group_id, plot_id, location_id, draw) %>%
                 mutate(particle_id = list(seq(from = 1, to = Y_sim, by = 1))) %>%
                 unnest(cols = particle_id) %>%
                 ungroup() %>%
                 select(-Y_sim)
               
               sim_tally_long = 
                 sim_tally_long %>%
                 mutate(transect_id = row_number(obs_id))
               
               sim_diam_data = 
                 list(N = nrow(sim_tally_long),
                      P = fwdtallies_training_data$P,
                      G = fwdtallies_training_data$G,
                      TS = max(sim_tally_long$transect_id),
                      transect_id = sim_tally_long$transect_id,
                      plot_id = 
                        sim_tally_long %>%
                        group_by(plot_id, transect_id) %>%
                        summarise() %>%
                        ungroup() %>%
                        arrange(transect_id) %>%
                        pull(plot_id) %>%
                        as.integer(),
                      group_id = 
                        fwdtallies_training_data$group_id,
                      lb = fwddiams_training_data$lb,
                      ub = fwddiams_training_data$ub)
               
               # simulate a diameter for each particle
               intercept = 
                 fwddiams_posterior %>%
                 select(contains('intercept')) %>%
                 slice(d) %>%
                 as.data.frame() %>%
                 as.numeric()
               
               phi = 
                 fwddiams_posterior %>%
                 select(contains('phi')) %>%
                 slice(d) %>%
                 as.data.frame() %>%
                 as.numeric()
                       
               SD_plot = 
                 fwddiams_posterior %>%
                 select(contains('SD_plot')) %>%
                 slice(d) %>%
                 as.data.frame() %>%
                 as.numeric()
               
               SD_trans = 
                 fwddiams_posterior %>%
                 select(contains('SD_trans')) %>%
                 slice(d) %>%
                 as.data.frame() %>%
                 as.numeric()
               
               zPlot = rnorm(n = fwdtallies_training_data$P,
                             mean = 0,
                             sd = 1)
               
               zTrans = rnorm(n = sim_diam_data$TS,
                              mean = 0,
                              sd = 1)
               
               effectPlot = zPlot * SD_plot[sim_diam_data$group_id]
               effectTrans = zTrans * SD_trans[sim_diam_data$group_id[sim_diam_data$plot_id]]
               
               logmu = 
                 sapply(X = 1:sim_diam_data$N,
                        FUN = function(i){
                          intercept[sim_diam_data$group_id[sim_diam_data$plot_id[sim_diam_data$transect_id[i]]]]+
                            effectPlot[sim_diam_data$plot_id[sim_diam_data$transect_id[i]]]+
                            effectTrans[sim_diam_data$transect_id[i]]
                        })
               
               phi = 
                 sapply(X = 1:sim_diam_data$N,
                        FUN = function(i){
                          phi[sim_diam_data$group_id[sim_diam_data$plot_id[sim_diam_data$transect_id[i]]]]
                        })
               
               shape = phi+2
             
               scale = exp(logmu)*(1+phi)
    
               diam_cm = 
                 sapply(X = 1:sim_diam_data$N,
                        FUN = function(i){
                          y = 0
                          while (y <= fwddiams_training_data$lb | 
                                 y >= fwddiams_training_data$ub){
                            # this is really confusing, but it looks like 
                            # the rinvgamma function calls the beta 
                            # parameter the rate, and 1/beta the scale.
                            # comparing the PDF formulae given in 
                            # the rinvgamma documentation, bourgingnon 2020,
                            # and the stan invgamma documentation, its 
                            # clear that what stan and bourgingnon call the
                            # scale is called the rate in rinvgamma
                            y = invgamma::rinvgamma(n = 1, shape = shape[i], rate = scale[i])
                          }
                          return(y)
                        })
               
               sim_tally_long$diam_cm = diam_cm
               
               #sim_tally_long$draw = d
               #return(sim_tally_long)
               
               # aggregate back up to 1h, 10h, 100h tallies
               sim_fwd = 
                 sim_tally_long %>%
                 mutate(timelag_class = 
                          ifelse(diam_cm > 0 & diam_cm <= 0.635,
                                 '1h',
                                 ifelse(diam_cm > 0.635 & diam_cm <= 2.54,
                                        '10h',
                                        ifelse(diam_cm > 2.54 & diam_cm <= 7.63,
                                               '100h',
                                               NA)))) %>%
                 group_by(draw, group_id, plot_id, location_id, obs_id, timelag_class) %>%
                 summarise(count = n()) %>%
                 ungroup()
              
               # add back in the 0 counts
               results = 
                 sim_tally_data  %>%
                 select(draw, group_id, plot_id, location_id, obs_id) %>%
                 expand(nesting(draw, group_id, plot_id, location_id, obs_id),
                          timelag_class = c('1h', '10h', '100h')) %>%
                 left_join(sim_fwd,
                           by = c('draw' = 'draw',
                                  'group_id' = 'group_id',
                                  'plot_id' = 'plot_id',
                                  'location_id' = 'location_id',
                                  'obs_id' = 'obs_id',
                                  'timelag_class' = 'timelag_class')) %>%
                 mutate(count = ifelse(is.na(count),0,count))
               
               return(results)
               
             })
          )
  

# load the observed data, estimate fuel load for each sample, aggregate to plot level
fwd_observed = 
  
  readRDS(here::here('02-data', '03-clean', 'seki_fwdtallies.rds'))%>%
  left_join(readRDS(here::here('02-data', '03-clean', 'seki_metadata.rds')) %>%
              left_join(readRDS(here::here('02-data', '02-intermediate', 'plot_locations_true.rds')) %>%
                          select(mort, plot_id, x_center = x_true, y_center = y_true),
                        by = c('plot_id')) %>%
              
              select(plot_id, transect_id, mort)) %>%
  
  # reformat tallies to longwise
  pivot_longer(cols = c(a1h, a10h, a100h, b1h, b10h, b100h),
               names_to = 'subsample_id',
               values_to = 'count') %>%
  mutate(timelag_class = gsub(x = subsample_id,
                              pattern = 'a|b',
                              replacement = ''),
         subsample = gsub(x = subsample_id,
                          pattern = '1h|10h|100h',
                          replacement = ''))  %>%
  
  # add columns for slope, QMD, SEC, SG, and transect length, and the conversion 
  # constant k
  left_join(., 
            y = QMDcm) %>%
  left_join(.,
            y = SEC) %>%
  left_join(.,
            y = SG) %>%
  mutate(slp_c = 1) %>%
  mutate(transect_length_m = 1,
         k = 1.234) %>%
  
  # use brown's equations to estimate fuel load
  mutate(mgha = 
           (k*weighted_qmd*weighted_sec*slp_c*weighted_sg*count)/
           transect_length_m)  %>%
  select(mort, plot_id, transect_id, location_m, subsample, timelag_class,count, mgha)

fwd_counts_plot = 
  ggplot()+
  stat_ecdf(
    data = fwddiams_simulatedtallies_data,
    aes(x = count, group = draw),
    color = 'blue',
    alpha = 0.1
  )+
  stat_ecdf(
    data = fwd_observed,
    aes(x = count),
    color = 'red'
  )+
  facet_grid(timelag_class~group_id)+
  theme_minimal()+
  scale_x_log10()

ggsave(fwd_counts_plot,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwd_counts_sim.png'),
       height = 6.5, width = 6.5, units = 'in')


fwddiams_simulatedtallies_data = 
  fwddiams_simulatedtallies_data %>%
  mutate(mort = ifelse(group_id==1,'low',
                       ifelse(group_id==2,'mid',
                              'high'))) %>%
  left_join(QMDcm) %>%
  left_join(SG) %>%
  left_join(SEC) %>%
  mutate(slp_c = 1,
         transect_length_m = 1,
         k = 1.234) %>%
  mutate(mgha = (k*weighted_qmd*weighted_sg*weighted_sec*slp_c*count)/transect_length_m) %>%
  select(draw, group_id, plot_id, location_id, obs_id, timelag_class,count, mgha)

# fifth, aggregate to plot level
fwd_sim_plots = 
  fwddiams_simulatedtallies_data %>%
  group_by(draw, group_id, plot_id, location_id, obs_id) %>%
  summarise(mgha = sum(mgha)) %>%
  ungroup() %>%
  group_by(draw, group_id, plot_id) %>%
  summarise(mgha = mean(mgha)) %>%
  ungroup() %>%
  mutate(source = 'simulated')

# seventh, combine the observed and simulated data
fwd_sim_plots = 
  fwd_sim_plots %>%
  bind_rows(fwd_observed %>%
              mutate(group_id = ifelse(mort=='low',1,
                                       ifelse(mort=='mid',2,
                                              3)),
                     source = 'observed',
                     draw = 4001,
                     plot_id = as.integer(factor(plot_id))) %>%
              group_by(source, draw, group_id, plot_id, 
                       transect_id, location_m, subsample) %>%
              summarise(mgha = sum(mgha,na.rm = TRUE)) %>%
              ungroup() %>%
              group_by(source, draw, group_id, plot_id) %>%
              summarise(mgha = mean(mgha, na.rm = TRUE)) %>%
              ungroup())




# finally, plot them
plotlevel_fwd_cdf = 
  ggplot()+
  stat_ecdf(
    data = fwd_sim_plots %>%
      filter(source == 'simulated'),
    aes(x = mgha, group = draw),
    geom = 'step',
    alpha = 0.1,
    pad = FALSE,
    color = 'blue')+
  stat_ecdf(
    data = fwd_sim_plots %>%
      filter(source=='observed'),
    aes(x = mgha),
    geom = 'step',
    lwd = 1,
    pad = FALSE,
    color = 'red'
  )+
  scale_x_log10()+
  theme_minimal()+
  scale_color_manual(values = c('observed' = 'red', 'simulated' = 'blue'))+
  labs(x = 'Plot mean FWD load (Mg / ha)', y = 'Cumulative density function')

plotlevel_fwd_cdf


ggsave(plotlevel_fwd_cdf,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwd_plotlevel.png'),
       height = 6.5, width = 6.5, units = 'in')

plotlevel_fwd_cdf2 = 
  fwd_sim_plots %>%
  filter(source=='simulated') %>%
  group_by(source, draw) %>%
  mutate(plot_rank = row_number(mgha)) %>%
  ungroup() %>%
  group_by(source, plot_rank) %>%
  summarise(fwd_mgha.01 = quantile(mgha, 0.01),
            fwd_mgha.10 = quantile(mgha, 0.1),
            fwd_mgha.25 = quantile(mgha, 0.25),
            fwd_mgha.75 = quantile(mgha, 0.75),
            fwd_mgha.90 = quantile(mgha, 0.90),
            fwd_mgha.99 = quantile(mgha, 0.99)) %>%
  ungroup() %>%
  mutate(y_min = (plot_rank/21)-(1/21),
         y_max = (plot_rank/21)) %>%
  ggplot()+
  geom_rect(aes(xmin = fwd_mgha.01, xmax = fwd_mgha.99, ymin = y_min, ymax = y_max,
                fill = 'simulated'),
            alpha = 0.25)+
  geom_rect(aes(xmin = fwd_mgha.10, xmax = fwd_mgha.90, ymin = y_min, ymax = y_max,
                fill = 'simulated'),
            alpha = 0.5)+
  geom_rect(aes(xmin = fwd_mgha.25, xmax = fwd_mgha.75, ymin = y_min, ymax = y_max,
                fill = 'simulated'),
            alpha = 1)+
  stat_ecdf(
    data = fwd_sim_plots %>%
      filter(source=='observed'),
    aes(x = mgha, color = 'observed'),
    geom = 'step',
    lwd = 1)+
  scale_x_log10()+
  theme_minimal()+
  scale_fill_manual(values = c('simulated' = 'blue'))+
  scale_color_manual(values = c('observed' = 'red'))+
  labs(x = 'Plot mean fwd load (Mg / ha)', y = 'Cumulative density function')+
  theme(legend.title = element_blank())

plotlevel_fwd_cdf2

ggsave(plotlevel_fwd_cdf2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwd_plotlevel_2.png'),
       height = 6.5, width = 6.5, units = 'in')


rm(fwddiams_simulatedtallies_data)

#### coarse woody debris #######################################################


# first, simulate tallies
cwdtallies_training_data = 
  readRDS(here::here('02-data', 
                     '05-for_analysis',
                     'cwdtallies_training_data.rds'))

cwddiams_training_data = 
  readRDS(here::here('02-data',
                     '05-for_analysis',
                     'cwddiams_training_data.rds'))

cwdtallies_posterior = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'cwdtallies_fit.rds'))$draws() %>%
  as_draws_df()

cwddiams_posterior = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'cwddiams_fit.rds'))$draws() %>%
  as_draws_df()

set.seed(110819)
cwdtallies_sim = 
  do.call('bind_rows',
          lapply(X = 1:4000,
                 FUN = function(draw){
                   
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
                   
                   sigmaPlot = 
                     cwdtallies_posterior %>%
                     select(contains('sigmaPlot')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   sigmaTrans = 
                     cwdtallies_posterior %>%
                     select(contains('sigmaTrans')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   zPlot = 
                     rnorm(n = cwdtallies_training_data$P,
                           mean = 0,
                           sd = 1)
                   
                   zTrans = 
                     rnorm(n = cwdtallies_training_data$TR,
                           mean = 0, sd = 1)
                   
                   plotEffects = 
                     zPlot * sigmaPlot[cwdtallies_training_data$group_id]
                   
                   transEffects = 
                     zPlot * 
                     sigmaTrans[cwdtallies_training_data$group_id[cwdtallies_training_data$plot_id]]
                   
                   zGP = 
                     matrix(nrow = cwdtallies_training_data$L,
                            ncol = cwdtallies_training_data$P,
                            byrow = FALSE,
                            data = 
                              rnorm(n = cwdtallies_training_data$L*cwdtallies_training_data$P,
                                    mean = 0,
                                    sd = 1)) 
                   
                   dists = 
                     cwdtallies_training_data$coords %>%
                     dist() %>%
                     as.matrix()
                   
                   # construct GP effects
                   K = 
                     array(dim = 
                             c(cwdtallies_training_data$L,
                               cwdtallies_training_data$L,
                               cwdtallies_training_data$G),
                           dimnames = 
                             list('l' = 1:cwdtallies_training_data$L,
                                  'lprime' = 1:cwdtallies_training_data$L,
                                  'g' = 1:cwdtallies_training_data$G),
                           data = 
                             sapply(X = 1:cwdtallies_training_data$G,
                                    FUN = function(g){
                                      sapply(X = 1:cwdtallies_training_data$L,
                                             FUN = function(l){
                                               sapply(X = 1:cwdtallies_training_data$L,
                                                      FUN = function(lprime){
                                                        ((alpha[g]**2)*
                                                          exp(-(dists[lprime,l]**2)/
                                                                (2*rho[g]**2)))
                                                      })
                                             })
                                    }))
                   
                   for (g in 1:cwdtallies_training_data$G){
                     K[,,g] = K[,,g]+diag(nrow = nrow(K[,,g]),
                                          ncol = ncol(K[,,g]),
                                          x = 1e-9)
                   }
                   
                   KL = 
                     array(dim = 
                             c(cwdtallies_training_data$L,
                               cwdtallies_training_data$L,
                               cwdtallies_training_data$G),
                           dimnames = 
                             list('l' = 1:cwdtallies_training_data$L,
                                  'lprime' = 1:cwdtallies_training_data$L,
                                  'g' = 1:cwdtallies_training_data$G),
                           data = 
                             sapply(X = 1:cwdtallies_training_data$G,
                                    FUN = function(g){
                                      t(chol(K[,,g]))
                                    }))
                 
                   # realized random effects
                   GP = 
                     matrix(nrow = cwdtallies_training_data$L,
                            ncol = cwdtallies_training_data$P,
                            data = 
                              sapply(X = 1:cwdtallies_training_data$P,
                                     FUN = function(p){
                                       KL[,,cwdtallies_training_data$group_id[p]] %*%
                                         zGP[,p]
                                     }))
                   
                   # linear predictors
                   logmu = 
                     sapply(X = 1:cwdtallies_training_data$N,
                            FUN = function(i){
                              
                               beta[cwdtallies_training_data$group_id[
                                 cwdtallies_training_data$plot_id[
                                   cwdtallies_training_data$transect_id[i]]]] +
                                plotEffects[
                                  cwdtallies_training_data$plot_id[
                                    cwdtallies_training_data$transect_id[i]]]+
                                as.numeric(GP[cwdtallies_training_data$location_id[i],
                                              cwdtallies_training_data$plot_id[
                                                cwdtallies_training_data$transect_id[i]]]) +
                                transEffects[cwdtallies_training_data$transect_id[i]]
                              
                            })
                    
                   # realizations
                   Y = 
                     sapply(X = 1:cwdtallies_training_data$N,
                            FUN = function(i){
                              rpois(n = 1,
                                    lambda = exp(logmu[i]))})
                   
                   
                   results = 
                     data.frame(Y_obs = cwdtallies_training_data$Y,
                                Y_sim = Y,
                                group_id = 
                                  cwdtallies_training_data$group_id[
                                    cwdtallies_training_data$plot_id[
                                      cwdtallies_training_data$transect_id
                                    ]],
                                plot_id = cwdtallies_training_data$plot_id[
                                  cwdtallies_training_data$transect_id
                                ],
                                transect_id = 
                                  cwdtallies_training_data$transect_id,
                                location_id = cwdtallies_training_data$location_id,
                                logmu = logmu,
                                draw = draw) %>%
                     as_tibble() %>%
                     rowid_to_column('obs_id')
                   
                   return(results)
                 }))

set.seed(110819)
# next, reformat the simulated tallies to be input data for the diameters model, 
# and simulate a diameter for each particle
cwddiams_simulatedtallies_data = 
  do.call('rbind',
          lapply(X = 1:4000,#sample(1:4000, 300),#1:100,
             FUN = function(d){
               
               print(paste0('Simulating draw ', d))
               
               # reformat the simulated tallies to be input data for the diameters model
               
               sim_tally_data = cwdtallies_sim %>% filter(draw == d)
               
               #sim_tally_long = 
              #   do.call('bind_rows',
               #          lapply(X = 1:nrow(sim_tally_data),
                #                FUN = function(i){
                 #                 if (sim_tally_data$Y_sim[i]==0){
                #                    return(NULL)
                 #                 }
                #                  longwise = 
                 #                   sim_tally_data[i,] %>%
                 #                   select(obs_id, group_id, plot_id, location_id, draw) %>%
                  #                  right_join(data.frame(particle_id = 
                   #                                         1:sim_tally_data$Y_sim[i]),
                    #                           by = character())
                     #             return(longwise)
                      #          }))
               
               sim_tally_long = 
                 sim_tally_data %>%
                 filter(Y_sim > 0) %>%
                 select(obs_id, group_id, plot_id, location_id, transect_id, draw, Y_sim) %>%
                 group_by(obs_id, group_id, plot_id, transect_id, location_id, draw) %>%
                 mutate(particle_id = list(seq(from = 1, to = Y_sim, by = 1))) %>%
                 unnest(cols = particle_id) %>%
                 ungroup() %>%
                 select(-Y_sim)
               
               sim_diam_data = 
                 list(N = nrow(sim_tally_long),
                      P = cwdtallies_training_data$P,
                      G = cwdtallies_training_data$G,
                      plot_id = sim_tally_long$plot_id,
                      group_id = 
                        cwdtallies_training_data$group_id,
                      lb = cwddiams_training_data$lb,
                      ub = cwddiams_training_data$ub)
               
               # simulate a diameter for each particle
               mu = 
                 cwddiams_posterior %>%
                 select(contains('mu')) %>%
                 slice(d) %>%
                 as.data.frame() %>%
                 as.numeric()
               
               phi = 
                 cwddiams_posterior %>%
                 select(contains('phi')) %>%
                 slice(d) %>%
                 as.data.frame() %>%
                 as.numeric()
               
               mu_i = 
                 sapply(X = 1:sim_diam_data$N,
                        FUN = function(i){
                          mu[sim_diam_data$group_id[
                            sim_diam_data$plot_id
                          ]]
                        })
               
               phi_i = 
                 sapply(X = 1:sim_diam_data$N,
                        FUN = function(i){
                          phi[sim_diam_data$group_id[
                            sim_diam_data$plot_id
                          ]]
                        })
                 
               shape = phi_i+2
             
               scale = mu_i*(1+phi_i)
    
               diam_cm = 
                 sapply(X = 1:sim_diam_data$N,
                        FUN = function(i){
                          y = 0
                          while (y <= cwddiams_training_data$lb | 
                                 y >= cwddiams_training_data$ub){
                            # this is really confusing, but it looks like 
                            # the rinvgamma function calls the beta 
                            # parameter the rate, and 1/beta the scale.
                            # comparing the PDF formulae given in 
                            # the rinvgamma documentation, bourgingnon 2020,
                            # and the stan invgamma documentation, its 
                            # clear that what stan and bourgingnon call the
                            # scale is called the rate in rinvgamma
                            y = invgamma::rinvgamma(n = 1, shape = shape[i], rate = scale[i])
                          }
                          return(y)
                        })
               
               sim_tally_long$diam_cm = diam_cm
               
               #sim_tally_long$draw = d
               #return(sim_tally_long)
               
               # aggregate back up to sum of squared diams
               sim_cwd = 
                 sim_tally_long %>%
                 group_by(draw, group_id, plot_id, transect_id, 
                          location_id) %>%
                 summarise(ssd_cm2 = sum(diam_cm**2)) %>%
                 ungroup()
              
               # add back in the 0 counts
               results = 
                 sim_tally_data  %>%
                 select(draw, group_id, plot_id, transect_id, location_id) %>%
                 left_join(sim_cwd,
                           by = c('draw' = 'draw',
                                  'group_id' = 'group_id',
                                  'plot_id' = 'plot_id',
                                  'transect_id' = 'transect_id',
                                  'location_id' = 'location_id')) %>%
                 mutate(ssd_cm2 = ifelse(is.na(ssd_cm2),0,ssd_cm2))
               
               return(results)
               
             })
          )
  

# load the observed data, estimate fuel load for each sample, aggregate to plot level
cwd_observed = 
  # load the data
  readRDS(here::here('02-data',
                     '03-clean',
                     'seki_cwd.rds')) %>%
  
  # each piece counts for 1 piece, and bin locations into 1m subtransects
  mutate(count = 1,
         location_m = (floor(location_m/1)*1)+0.5) %>%

  select(plot_id, transect_id, az, location_m, count, diam_cm) %>%
  
  # add in rows for 0 intersections on each meter
  bind_rows(
    readRDS(here::here('02-data', '02-intermediate', 'plot_locations_true.rds')) %>%
      select(plot_id) %>%
      expand(plot_id,
             az = c(0, 90, 180, 270),
             location_m = seq(from = 0.5, to = 29.5, by = 1)) %>%
      mutate(diam_cm = 0)
  ) %>%
  
  # get the total number of intersections on each meter
  group_by(plot_id, az, location_m) %>%
  summarise(ssd_cm2 = sum(diam_cm**2)) %>%
  ungroup() %>%
  
  # get mortality level and transect id columns
  left_join(readRDS(here::here('02-data', '02-intermediate', 'plot_locations_true.rds')) %>%
              select(mort, plot_id),
            by = c('plot_id' = 'plot_id')) %>%  
  # want to treat the N and S transects as a single transect, same with the E/W
  mutate(direction = ifelse(is.element(az, c(0, 180)),'NS', 'EW'),
         transect_id = paste0(plot_id,'-',direction)) %>%
  mutate(plot_id.i = as.integer(factor(plot_id)),
         transect_id.i = as.integer(factor(transect_id)),
         group_id = as.integer(factor(mort, levels = c('low', 'mid', 'high')))) %>%
  mutate(draw = 4001) %>%
  select(draw, mort, group_id, plot_id = plot_id.i, transect_id = transect_id.i,
         location_m, ssd_cm2) %>%
  mutate(timelag_class = '1000s') %>%
  left_join(SEC) %>%
  left_join(SG) %>%
  mutate(slp_c = 1) %>%
  mutate(transect_length_m = 1,
         k = 1.234) %>%
  
  mutate(mgha = 
           (k*ssd_cm2*weighted_sec*slp_c*weighted_sg)/
           transect_length_m) %>%
  select(group_id, plot_id, transect_id, location_m, ssd_cm2, mgha)

head(cwddiams_simulatedtallies_data)

cwd_simulated = 
  cwddiams_simulatedtallies_data %>% 
  mutate(timelag_class = '1000s') %>%
  left_join(SEC) %>%
  left_join(SG) %>%
  mutate(slp_c = 1) %>%
  mutate(transect_length_m = 1,
         k = 1.234) %>%
  mutate(mgha = 
           (k*ssd_cm2*weighted_sec*slp_c*weighted_sg)/
           transect_length_m) %>%
  select(draw, group_id, plot_id, transect_id, location_id, ssd_cm2, mgha)

head(cwd_simulated)





cwd_simulated = 
  cwd_simulated %>%
  mutate(mort = ifelse(group_id==1,'low',
                       ifelse(group_id==2,'mid',
                              'high')),
         timelag_class = '1000s') %>%
  left_join(SG) %>%
  left_join(SEC) %>%
  mutate(slp_c = 1,
         transect_length_m = 1,
         k = 1.234) %>%
  mutate(mgha = (k*weighted_sg*weighted_sec*slp_c*ssd_cm2)/transect_length_m) %>%
  select(draw, group_id, plot_id, transect_id, location_id, timelag_class, mgha)

# fifth, aggregate to plot level
cwd_sim_plots = 
  cwd_simulated %>%
  group_by(draw, group_id, plot_id) %>%
  summarise(mgha = mean(mgha)) %>%
  ungroup() %>%
  mutate(source = 'simulated') %>%
  bind_rows(cwd_observed %>%
              mutate(group_id = ifelse(mort=='low',1,
                                       ifelse(mort=='mid',2,
                                              3)),
                     source = 'observed',
                     draw = 4001,
                     plot_id = as.integer(factor(plot_id))) %>%
              group_by(source, draw, group_id, plot_id) %>%
              summarise(mgha = mean(mgha, na.rm = TRUE)) %>%
              ungroup())

# finally, plot them
plotlevel_cwd_cdf = 
  ggplot()+
  stat_ecdf(
    data = cwd_sim_plots %>%
      filter(source == 'simulated'),
    aes(x = mgha, group = draw),
    geom = 'step',
    alpha = 0.1,
    pad = FALSE,
    color = 'blue')+
  stat_ecdf(
    data = cwd_sim_plots %>%
      filter(source=='observed'),
    aes(x = mgha),
    geom = 'step',
    lwd = 1,
    pad = FALSE,
    color = 'red'
  )+
  scale_x_log10()+
  theme_minimal()+
  scale_color_manual(values = c('observed' = 'red', 'simulated' = 'blue'))+
  labs(x = 'Plot mean cwd load (Mg / ha)', y = 'Cumulative density function')

plotlevel_cwd_cdf


ggsave(plotlevel_cwd_cdf,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'cwd_plotlevel.png'),
       height = 6.5, width = 6.5, units = 'in')

plotlevel_cwd_cdf2 = 
  cwd_sim_plots %>%
  filter(source=='simulated') %>%
  group_by(source, draw) %>%
  mutate(plot_rank = row_number(mgha)) %>%
  ungroup() %>%
  group_by(source, plot_rank) %>%
  summarise(cwd_mgha.01 = quantile(mgha, 0.01),
            cwd_mgha.10 = quantile(mgha, 0.1),
            cwd_mgha.25 = quantile(mgha, 0.25),
            cwd_mgha.75 = quantile(mgha, 0.75),
            cwd_mgha.90 = quantile(mgha, 0.90),
            cwd_mgha.99 = quantile(mgha, 0.99)) %>%
  ungroup() %>%
  mutate(y_min = (plot_rank/21)-(1/21),
         y_max = (plot_rank/21)) %>%
  ggplot()+
  geom_rect(aes(xmin = cwd_mgha.01, xmax = cwd_mgha.99, ymin = y_min, ymax = y_max,
                fill = 'simulated'),
            alpha = 0.25)+
  geom_rect(aes(xmin = cwd_mgha.10, xmax = cwd_mgha.90, ymin = y_min, ymax = y_max,
                fill = 'simulated'),
            alpha = 0.5)+
  geom_rect(aes(xmin = cwd_mgha.25, xmax = cwd_mgha.75, ymin = y_min, ymax = y_max,
                fill = 'simulated'),
            alpha = 1)+
  stat_ecdf(
    data = cwd_sim_plots %>%
      filter(source=='observed'),
    aes(x = mgha, color = 'observed'),
    geom = 'step',
    lwd = 1)+
  scale_x_log10()+
  theme_minimal()+
  scale_fill_manual(values = c('simulated' = 'blue'))+
  scale_color_manual(values = c('observed' = 'red'))+
  labs(x = 'Plot mean cwd load (Mg / ha)', y = 'Cumulative density function')+
  theme(legend.title = element_blank())

plotlevel_cwd_cdf2

ggsave(plotlevel_cwd_cdf2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'cwd_plotlevel_2.png'),
       height = 6.5, width = 6.5, units = 'in')


rm(cwddiams_simulatedtallies_data, cwd_simulated)


#### understory vegetation presences ###########################################


vegpa_training_data = 
  readRDS(here::here('02-data', 
                     '05-for_analysis',
                     'vegpa_training_data.rds'))

vegpa_posterior = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'vegpa_fit.rds'))$draws() %>%
  as_draws_df()

set.seed(110819)
veg_presences = 
  do.call('bind_rows',
          lapply(X = 1:4000,
                 FUN = function(draw){
                   
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
                   
                   sigmaPlot = 
                     vegpa_posterior %>%
                     select(contains('sigmaPlot')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   zPlot = 
                     rnorm(n = vegpa_training_data$P,
                           mean = 0,
                           sd = 1)
                   plotEffects = 
                     zPlot * sigmaPlot[vegpa_training_data$group_id]
                   
                   zGP = 
                     matrix(nrow = vegpa_training_data$L,
                            ncol = vegpa_training_data$P,
                            byrow = FALSE,
                            data = 
                              rnorm(n = vegpa_training_data$L*
                                      vegpa_training_data$P,
                                    mean = 0,
                                    sd = 1))
                   
                   
                   dists = 
                     vegpa_training_data$coords %>%
                     dist() %>%
                     as.matrix()
                   
                   # construct GP effects
                   K = 
                     array(dim = 
                             c(vegpa_training_data$L,
                               vegpa_training_data$L,
                               vegpa_training_data$G),
                           dimnames = 
                             list('l' = 1:vegpa_training_data$L,
                                  'lprime' = 1:vegpa_training_data$L,
                                  'g' = 1:vegpa_training_data$G),
                           data = 
                             sapply(X = 1:vegpa_training_data$G,
                                    FUN = function(g){
                                      sapply(X = 1:vegpa_training_data$L,
                                             FUN = function(l){
                                               sapply(X = 1:vegpa_training_data$L,
                                                      FUN = function(lprime){
                                                        ((alpha[g]**2)*
                                                          exp(-(dists[lprime,l]**2)/
                                                                (2*rho[g]**2)))
                                                      })
                                             })
                                    }))
                   
                   for (g in 1:vegpa_training_data$G){
                     K[,,g] = K[,,g]+diag(nrow = nrow(K[,,g]),
                                          ncol = ncol(K[,,g]),
                                          x = 1e-9)
                   }
                   
                   KL = 
                     array(dim = 
                             c(vegpa_training_data$L,
                               vegpa_training_data$L,
                               vegpa_training_data$G),
                           dimnames = 
                             list('l' = 1:vegpa_training_data$L,
                                  'lprime' = 1:vegpa_training_data$L,
                                  'g' = 1:vegpa_training_data$G),
                           data = 
                             sapply(X = 1:vegpa_training_data$G,
                                    FUN = function(g){
                                      t(chol(K[,,g]))
                                    }))
                 
                   # realized random effects
                   GP = 
                     matrix(nrow = vegpa_training_data$L,
                            ncol = vegpa_training_data$P,
                            data = 
                              sapply(X = 1:vegpa_training_data$P,
                                     FUN = function(p){
                                       KL[,,vegpa_training_data$group_id[p]] %*%
                                         zGP[,p]
                                     }))
                   
                   # linear predictors
                   logiteta = 
                     sapply(X = 1:vegpa_training_data$N,
                            FUN = function(i){
                              
                               beta[vegpa_training_data$group_id[vegpa_training_data$plot_id[i]]] +
                                plotEffects[vegpa_training_data$plot_id[i]]+
                                as.numeric(GP[vegpa_training_data$location_id[i],
                                              vegpa_training_data$plot_id][i])
                              
                            })
                   
                   # realizations
                   Y = 
                     sapply(X = 1:vegpa_training_data$N,
                            FUN = function(i){
                              rbinom(n = 1,
                                     size = 1,
                                     prob = boot::inv.logit(logiteta[i]))})
                   
                   
                   results = 
                     data.frame(Y_sim = Y,
                                group_id = 
                                  vegpa_training_data$group_id[vegpa_training_data$plot_id],
                                plot_id = vegpa_training_data$plot_id,
                                location_id = vegpa_training_data$location_id,
                                logiteta = logiteta,
                                draw = draw) %>%
                     as_tibble() %>%
                     rowid_to_column('obs_id')
                   
                   return(results)
                 }))

# add in the observed data
veg_presences = 
  veg_presences %>%
  mutate(source = 'simulated') %>%
  bind_rows(
    data.frame(
      draw = 4001,
      source = 'observed',
      Y_sim = vegpa_training_data$Y,
      obs_id = 1:vegpa_training_data$N,
      group_id = vegpa_training_data$group_id[vegpa_training_data$plot_id],
      plot_id = vegpa_training_data$plot_id,
      location_id = vegpa_training_data$location_id,
      logiteta = NA
    )
  )

# using the transect lengths as weights, randomly pull a spp and height for 
# each pixel
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
  ungroup()

head(veg_observed)  

veg_presences = 
  veg_presences %>%
  filter(group_id==1) %>%
  bind_cols(
    veg_observed %>%
      filter(group_id==1) %>%
      slice_sample(n = nrow(veg_presences %>% filter(group_id==1)),
                   weight_by = transect_length_m,
                   replace = TRUE) %>%
      select(spp, status, height_m)
                   
  ) %>%
  bind_rows(
    veg_presences %>%
      filter(group_id==2) %>%
      bind_cols(
        veg_observed %>%
          filter(group_id==2) %>%
          slice_sample(n = nrow(veg_presences %>% filter(group_id==2)),
                       weight_by = transect_length_m,
                       replace = TRUE) %>%
          select(spp, status, height_m)
      )
  ) %>%
  bind_rows(
    veg_presences %>%
      filter(group_id==3) %>%
      bind_cols(
        veg_observed %>%
          filter(group_id==3) %>%
          slice_sample(n = nrow(veg_presences %>% filter(group_id==3)),
                       weight_by = transect_length_m,
                       replace = TRUE) %>%
          select(spp, status, height_m)
      )
  )


head(veg_presences)

veg_cover_plot_1 = 
  ggplot()+
  stat_ecdf(
    data = 
      veg_presences %>%
      filter(source=='simulated' & is.element(draw, sample(1:4000,100))) %>%
      group_by(draw, source, group_id, plot_id) %>%
      summarise(pcover = sum(Y_sim)/n()) %>%
      ungroup(),
    aes(x = pcover, group = draw), 
    color = 'blue',
    alpha = 0.25
  )+
  stat_ecdf(
    data = 
      veg_presences %>%
      filter(source == 'observed') %>%
      group_by(draw, source, group_id, plot_id) %>%
      summarise(pcover = sum(Y_sim)/n()) %>%
      ungroup(),
    aes(x = pcover), 
    color = 'red',
    lwd = 1
  )+
  scale_color_manual(values = c('observed' = 'red', 'simulated' = 'blue'))+
  theme_minimal()

veg_cover_plot_1
  
ggsave(veg_cover_plot_1,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'veg_plotlevel_1.png'),
       width = 6.5, height = 4.5, units = 'in')



veg_cover_plot_2 = 
  veg_presences %>%
  filter(source=='simulated') %>%
  group_by(source, draw, group_id, plot_id) %>%
  summarise(pcover = sum(Y_sim)/n()) %>%
  ungroup() %>%
  group_by(source, draw) %>%
  mutate(plot_rank = row_number(pcover)) %>%
  ungroup() %>%
  group_by(source, plot_rank) %>%
  summarise(pcover.01 = quantile(pcover, 0.01),
            pcover.10 = quantile(pcover, 0.1),
            pcover.25 = quantile(pcover, 0.25),
            pcover.75 = quantile(pcover, 0.75),
            pcover.90 = quantile(pcover, 0.90),
            pcover.99 = quantile(pcover, 0.99)) %>%
  ungroup() %>%
  mutate(y_min = (plot_rank/21)-(1/21),
         y_max = (plot_rank/21)) %>%
  ggplot()+
  geom_rect(aes(xmin = pcover.01, xmax = pcover.99, ymin = y_min, ymax = y_max,
                fill = 'simulated'),
            alpha = 0.25)+
  geom_rect(aes(xmin = pcover.10, xmax = pcover.90, ymin = y_min, ymax = y_max,
                fill = 'simulated'),
            alpha = 0.5)+
  geom_rect(aes(xmin = pcover.25, xmax = pcover.75, ymin = y_min, ymax = y_max,
                fill = 'simulated'),
            alpha = 1)+
  stat_ecdf(
    data = 
      veg_presences %>%
      filter(source == 'observed') %>%
      group_by(draw, source, group_id, plot_id) %>%
      summarise(pcover = sum(Y_sim)/n()) %>%
      ungroup(),
    aes(x = pcover), 
    color = 'red',
    lwd = 1
  )+
  theme_minimal()+
  scale_fill_manual(values = c('simulated' = 'blue'))+
  scale_color_manual(values = c('observed' = 'red'))+
  labs(x = 'Understory vegetation cover proportion', 
       y = 'Cumulative density function')+
  theme(legend.title = element_blank())

veg_cover_plot_2

ggsave(veg_cover_plot_2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'veg_plotlevel_2.png'),
       width = 6.5, height = 4.5, units = 'in')


head(veg_presences)


#### trees #####################################################################

# load the observed data

# load the posterior distribution
trees_posterior = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'trees_fit.rds'))$draws() %>%
  as_draws_df()

trees_training_data = 
  readRDS(here::here('02-data', 
                     '05-for_analysis',
                     'trees_training_data.rds'))

# simulate 4k datasets, including new draws of random effect realizations
set.seed(110819)
trees_sim = 
  do.call('bind_rows',
          lapply(X = 1:4000,
                 FUN = function(draw){
                   
                   print(paste0('Working on draw ', draw))
                   
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
                   
                   sigmaPlot = 
                     trees_posterior %>%
                     select(contains('sigmaPlot')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   dists = 
                     trees_training_data$coords %>%
                     dist() %>%
                     as.matrix()
                   
                   # draw new random effect realizations
                   zPlot = rnorm(n = trees_training_data$P,
                                 mean = 0,
                                 sd = 1)
                   
                   zGP = 
                     matrix(nrow = trees_training_data$L,
                            ncol = trees_training_data$P,
                            byrow = FALSE,
                            data = 
                              rnorm(n = trees_training_data$L*trees_training_data$P,
                                    mean = 0,
                                    sd = 1))
                   
                   # construct plot effects
                   plotEffects = zPlot * sigmaPlot[trees_training_data$group_id]
                   
                   # construct GP effects
                   K = 
                     array(dim = 
                             c(trees_training_data$L,
                               trees_training_data$L,
                               trees_training_data$G),
                           dimnames = 
                             list('l' = 1:trees_training_data$L,
                                  'lprime' = 1:trees_training_data$L,
                                  'g' = 1:trees_training_data$G),
                           data = 
                             sapply(X = 1:trees_training_data$G,
                                    FUN = function(g){
                                      sapply(X = 1:trees_training_data$L,
                                             FUN = function(l){
                                               sapply(X = 1:trees_training_data$L,
                                                      FUN = function(lprime){
                                                        ((alpha[g]**2)*
                                                           exp(-(dists[lprime,l]**2)/
                                                                 (2*rho[g]**2)))
                                                      })
                                             })
                                    }))
                   
                   for (g in 1:trees_training_data$G){
                     K[,,g] = K[,,g]+diag(nrow = nrow(K[,,g]),
                                          ncol = ncol(K[,,g]),
                                          x = 1e-9)
                   }
                   
                   KL = 
                     array(dim = 
                             c(trees_training_data$L,
                               trees_training_data$L,
                               trees_training_data$G),
                           dimnames = 
                             list('l' = 1:trees_training_data$L,
                                  'lprime' = 1:trees_training_data$L,
                                  'g' = 1:trees_training_data$G),
                           data = 
                             sapply(X = 1:trees_training_data$G,
                                    FUN = function(g){
                                      t(chol(K[,,g]))
                                    }))
                   
                   # realized random effects
                   GP = 
                     matrix(nrow = trees_training_data$L,
                            ncol = trees_training_data$P,
                            data = 
                              sapply(X = 1:trees_training_data$P,
                                     FUN = function(p){
                                       KL[,,trees_training_data$group_id[p]] %*%
                                         zGP[,p]
                                     }))
                   
                   # linear predictors
                   logmu = 
                     sapply(X = 1:trees_training_data$N,
                            FUN = function(i){
                              
                              beta[trees_training_data$group_id[trees_training_data$plot_id[i]]] +
                                plotEffects[trees_training_data$plot_id[i]]+
                                as.numeric(GP[trees_training_data$location_id[i],
                                              trees_training_data$plot_id[i]])
                              
                            })
                   
                   # realizations
                   Y = 
                     sapply(X = 1:trees_training_data$N,
                            FUN = function(i){
                              rpois(n = 1,
                                    lambda = exp(logmu[i]))
                            })
                   
                   
                   results = 
                     data.frame(Y = Y,
                                group_id = 
                                  trees_training_data$group_id[trees_training_data$plot_id],
                                plot_id = trees_training_data$plot_id,
                                location_id = trees_training_data$location_id,
                                logmu = logmu,
                                draw = draw,
                                source = 'sim') %>%
                     as_tibble() %>%
                     rowid_to_column('obs_id')
                   
                   return(results)
                 }))



# combine observed and simulated datasets
trees_sim = 
  trees_sim %>%
  bind_rows(
    data.frame(Y = trees_training_data$Y,
               group_id = trees_training_data$group_id[trees_training_data$plot_id],
               plot_id = trees_training_data$plot_id,
               location_id = trees_training_data$location_id,
               logmu = NA,
               draw = 4001,
               source = 'observed')
  )


head(trees_sim)

trees_sim_plots = 
  trees_sim %>%
  mutate(tph = Y/0.05) %>%
  group_by(group_id, plot_id, draw, source) %>%
  summarise(tph = sum(tph)) %>%
  ungroup() 

trees_sim_plots %>%
  ggplot(aes(x = source, y = tph))+
  geom_boxplot(
    data = trees_sim_plots %>%
      filter(source == 'sim'),
    width = 0.25,
    position = position_nudge(x = -0.35))+
  geom_jitter(
    data = trees_sim_plots %>%
      filter(source == 'observed'),
    height = 0, width = 0.15)+
  geom_violin(
    data = trees_sim_plots %>%
      filter(source == 'sim'),
    width = 0.25,
    position = position_nudge(x = 0.35))+
  facet_grid(group_id ~.)+
  scale_y_log10()+
  theme_minimal()

trees_sim_plots %>%
  filter(source=='sim') %>%
  pull(tph) %>%
  summary()

ggplot(data = trees_sim_plots,
       aes(x = tph, color = source, group = draw))+
  stat_ecdf(geom = 'step',
            lwd = 1,
            pad = TRUE,
            alpha = 0.1)+
  scale_x_log10()+
  theme_minimal()

plotlevel_trees_cdf = 
  ggplot()+
  stat_ecdf(
    data = trees_sim_plots %>%
      filter(source == 'sim'),
    aes(x = tph, group = draw),
    geom = 'step',
    alpha = 0.1,
    pad = FALSE,
    color = 'blue')+
  stat_ecdf(
    data = trees_sim_plots %>%
      filter(source=='observed'),
    aes(x = tph),
    geom = 'step',
    lwd = 1,
    pad = FALSE,
    color = 'red'
  )+
  scale_x_log10()+
  theme_minimal()+
  scale_color_manual(values = c('observed' = 'red', 'simulated' = 'blue'))+
  labs(x = 'Plot TPH >11.4 cm DBH', y = 'Cumulative density function')

plotlevel_trees_cdf

ggsave(plotlevel_trees_cdf,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'trees_plotlevel.png'),
       height = 6.5, width = 6.5, units = 'in')


head(trees_sim_plots)

plotlevel_trees_cdf2 = 
  trees_sim_plots %>%
  filter(source=='sim') %>%
  group_by(source, draw) %>%
  mutate(plot_rank = row_number(tph)) %>%
  ungroup() %>%
  group_by(source, plot_rank) %>%
  summarise(tph.01 = quantile(tph, 0.01),
            tph.10 = quantile(tph, 0.1),
            tph.25 = quantile(tph, 0.25),
            tph.75 = quantile(tph, 0.75),
            tph.90 = quantile(tph, 0.90),
            tph.99 = quantile(tph, 0.99)) %>%
  ungroup() %>%
  mutate(y_min = (plot_rank/21)-(1/21),
         y_max = (plot_rank/21)) %>%
  ggplot()+
  geom_rect(aes(xmin = tph.01, xmax = tph.99, ymin = y_min, ymax = y_max,
                fill = 'simulated'),
            alpha = 0.25)+
  geom_rect(aes(xmin = tph.10, xmax = tph.90, ymin = y_min, ymax = y_max,
                fill = 'simulated'),
            alpha = 0.5)+
  geom_rect(aes(xmin = tph.25, xmax = tph.75, ymin = y_min, ymax = y_max,
                fill = 'simulated'),
            alpha = 1)+
  stat_ecdf(
    data = trees_sim_plots %>%
      filter(source=='observed'),
    aes(x = tph, color = 'observed'),
    geom = 'step',
    lwd = 1)+
  scale_x_log10()+
  theme_minimal()+
  scale_fill_manual(values = c('simulated' = 'blue'))+
  scale_color_manual(values = c('observed' = 'red'))+
  labs(x = 'Plot TPH >11.4 cm DBH', y = 'Cumulative density function')+
  theme(legend.title = element_blank())

plotlevel_trees_cdf2

ggsave(plotlevel_trees_cdf2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'trees_plotlevel_2.png'),
       height = 6.5, width = 6.5, units = 'in')



# plot results

rm(trees_sim)

#### saplings #####################################################################

# load the observed data

# load the posterior distribution
saplings_posterior = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'saplings_fit.rds'))$draws() %>%
  as_draws_df()

saplings_training_data = 
  readRDS(here::here('02-data', 
                     '05-for_analysis',
                     'saplings_training_data.rds'))

# simulate 4k datasets, including new draws of random effect realizations
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
                   
                   sigmaPlot = 
                     saplings_posterior %>%
                     select(contains('sigmaPlot')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   dists = 
                     saplings_training_data$coords %>%
                     dist() %>%
                     as.matrix()
                   
                   # draw new random effect realizations
                   zPlot = rnorm(n = saplings_training_data$P,
                                 mean = 0,
                                 sd = 1)
                   
                   zGP = 
                     matrix(nrow = saplings_training_data$L,
                            ncol = saplings_training_data$P,
                            byrow = FALSE,
                            data = 
                              rnorm(n = saplings_training_data$L*saplings_training_data$P,
                                    mean = 0,
                                    sd = 1))
                   
                   # construct plot effects
                   plotEffects = zPlot * sigmaPlot[saplings_training_data$group_id]
                   
                   # construct GP effects
                   K = 
                     array(dim = 
                             c(saplings_training_data$L,
                               saplings_training_data$L,
                               saplings_training_data$G),
                           dimnames = 
                             list('l' = 1:saplings_training_data$L,
                                  'lprime' = 1:saplings_training_data$L,
                                  'g' = 1:saplings_training_data$G),
                           data = 
                             sapply(X = 1:saplings_training_data$G,
                                    FUN = function(g){
                                      sapply(X = 1:saplings_training_data$L,
                                             FUN = function(l){
                                               sapply(X = 1:saplings_training_data$L,
                                                      FUN = function(lprime){
                                                        ((alpha[g]**2)*
                                                           exp(-(dists[lprime,l]**2)/
                                                                 (2*rho[g]**2)))
                                                      })
                                             })
                                    }))
                   
                   for (g in 1:saplings_training_data$G){
                     K[,,g] = K[,,g]+diag(nrow = nrow(K[,,g]),
                                          ncol = ncol(K[,,g]),
                                          x = 1e-9)
                   }
                   
                   KL = 
                     array(dim = 
                             c(saplings_training_data$L,
                               saplings_training_data$L,
                               saplings_training_data$G),
                           dimnames = 
                             list('l' = 1:saplings_training_data$L,
                                  'lprime' = 1:saplings_training_data$L,
                                  'g' = 1:saplings_training_data$G),
                           data = 
                             sapply(X = 1:saplings_training_data$G,
                                    FUN = function(g){
                                      t(chol(K[,,g]))
                                    }))
                   
                   # realized random effects
                   GP = 
                     matrix(nrow = saplings_training_data$L,
                            ncol = saplings_training_data$P,
                            data = 
                              sapply(X = 1:saplings_training_data$P,
                                     FUN = function(p){
                                       KL[,,saplings_training_data$group_id[p]] %*%
                                         zGP[,p]
                                     }))
                   
                   # linear predictors
                   logmu = 
                     sapply(X = 1:saplings_training_data$N,
                            FUN = function(i){
                              
                              beta[saplings_training_data$group_id[saplings_training_data$plot_id[i]]] +
                                plotEffects[saplings_training_data$plot_id[i]]+
                                as.numeric(GP[saplings_training_data$location_id[i],
                                              saplings_training_data$plot_id[i]])
                              
                            })
                   
                   # realizations
                   Y = 
                     sapply(X = 1:saplings_training_data$N,
                            FUN = function(i){
                              rpois(n = 1,
                                    lambda = exp(logmu[i]))
                            })
                   
                   
                   results = 
                     data.frame(Y = Y,
                                group_id = 
                                  saplings_training_data$group_id[saplings_training_data$plot_id],
                                plot_id = saplings_training_data$plot_id,
                                location_id = saplings_training_data$location_id,
                                logmu = logmu,
                                draw = draw,
                                source = 'sim') %>%
                     as_tibble() %>%
                     rowid_to_column('obs_id')
                   
                   return(results)
                 }))



# combine observed and simulated datasets
saplings_sim = 
  saplings_sim %>%
  bind_rows(
    data.frame(Y = saplings_training_data$Y,
               group_id = saplings_training_data$group_id[saplings_training_data$plot_id],
               plot_id = saplings_training_data$plot_id,
               location_id = saplings_training_data$location_id,
               logmu = NA,
               draw = 4001,
               source = 'observed')
  )


head(saplings_sim)

saplings_sim_plots = 
  saplings_sim %>%
  mutate(tph = Y/0.0116) %>%
  group_by(group_id, plot_id, draw, source) %>%
  summarise(tph = sum(tph)) %>%
  ungroup() 

saplings_sim_plots %>%
  group_by(source, draw) %>%
  mutate(rank = row_number(tph)) %>%
  ungroup() %>%
  filter(rank >= 20) %>%
  arrange(desc(tph))

saplings_sim_plots %>%
  group_by(source, draw) %>%
  mutate(rank = row_number(tph)) %>%
  ungroup() %>%
  filter(draw == 1170) %>%
  arrange(desc(rank)) %>%
  print(n = Inf)

saplings_sim %>% filter(draw == 1170) %>% pull(plot_id) %>% unique() %>% length()

saplings_sim %>% filter(draw==1170) %>%
  group_by(plot_id) %>%
  summarise(tph = sum(Y)/0.0116) %>%
  ungroup() %>%
  pull(tph) %>%
  unique() %>%
  length()


saplings_sim_plots %>%
  ggplot(aes(x = source, y = tph))+
  geom_boxplot(
    data = saplings_sim_plots %>%
      filter(source == 'sim'),
    width = 0.25,
    position = position_nudge(x = -0.35))+
  geom_jitter(
    data = saplings_sim_plots %>%
      filter(source == 'observed'),
    height = 0, width = 0.15)+
  geom_violin(
    data = saplings_sim_plots %>%
      filter(source == 'sim'),
    width = 0.25,
    position = position_nudge(x = 0.35))+
  facet_grid(group_id ~.)+
  scale_y_log10()+
  theme_minimal()

saplings_sim_plots %>%
  filter(source=='sim') %>%
  pull(tph) %>%
  summary()

ggplot(data = saplings_sim_plots,
       aes(x = tph, color = source, group = draw))+
  stat_ecdf(geom = 'step',
            lwd = 1,
            pad = TRUE,
            alpha = 0.1)+
  scale_x_log10()+
  theme_minimal()

plotlevel_saplings_cdf = 
  ggplot()+
  stat_ecdf(
    data = saplings_sim_plots %>%
      filter(source == 'sim'),
    aes(x = tph, group = draw),
    geom = 'step',
    alpha = 0.1,
    pad = FALSE,
    color = 'blue')+
  stat_ecdf(
    data = saplings_sim_plots %>%
      filter(source=='observed'),
    aes(x = tph),
    geom = 'step',
    lwd = 1,
    pad = FALSE,
    color = 'red'
  )+
  scale_x_log10()+
  theme_minimal()+
  scale_color_manual(values = c('observed' = 'red', 'simulated' = 'blue'))+
  labs(x = 'Plot TPH <11.4 cm DBH', y = 'Cumulative density function')

plotlevel_saplings_cdf

ggsave(plotlevel_saplings_cdf,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'saplings_plotlevel.png'),
       height = 6.5, width = 6.5, units = 'in')


head(saplings_sim_plots)

plotlevel_saplings_cdf2 = 
  saplings_sim_plots %>%
  filter(source=='sim') %>%
  group_by(source, draw) %>%
  mutate(plot_rank = row_number(tph)) %>%
  ungroup() %>%
  group_by(source, plot_rank) %>%
  summarise(tph.01 = quantile(tph+1, 0.01),
            tph.10 = quantile(tph+1, 0.1),
            tph.25 = quantile(tph+1, 0.25),
            tph.75 = quantile(tph+1, 0.75),
            tph.90 = quantile(tph+1, 0.90),
            tph.99 = quantile(tph+1, 0.99)) %>%
  ungroup() %>%
  mutate(y_min = (plot_rank/21)-(1/21),
         y_max = (plot_rank/21)) %>%
  ggplot()+
  geom_rect(aes(xmin = tph.01, xmax = tph.99, ymin = y_min, ymax = y_max,
                fill = 'simulated'),
            alpha = 0.25)+
  geom_rect(aes(xmin = tph.10, xmax = tph.90, ymin = y_min, ymax = y_max,
                fill = 'simulated'),
            alpha = 0.5)+
  geom_rect(aes(xmin = tph.25, xmax = tph.75, ymin = y_min, ymax = y_max,
                fill = 'simulated'),
            alpha = 1)+
  stat_ecdf(
    data = saplings_sim_plots %>%
      filter(source=='observed'),
    aes(x = tph, color = 'observed'),
    geom = 'step',
    lwd = 1)+
  scale_x_log10()+
  theme_minimal()+
  scale_fill_manual(values = c('simulated' = 'blue'))+
  scale_color_manual(values = c('observed' = 'red'))+
  labs(x = 'Plot TPH <11.4 cm DBH', y = 'Cumulative density function')+
  theme(legend.title = element_blank())

plotlevel_saplings_cdf2

ggsave(plotlevel_saplings_cdf2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'saplings_plotlevel_2.png'),
       height = 6.5, width = 6.5, units = 'in')


# all of the observed data are plausible assuming the model is true, but the 
# model predicts (rare) behaviors that seem unlikely

# plot results

rm(saplings_sim)



#### scratch ###################################################################

fwddiams_simulatedtallies_data %>%
  ggplot()+
  geom_line(aes(color = 'simulated', group = draw, x = diam_cm),
            stat = 'density',
            alpha = 0.1)+
  geom_density(
    data = 
      data.frame(
        group_id = 
          readRDS(here::here('02-data', '05-for_analysis', 'fwddiams_training_data.rds'))$group_id[
            readRDS(here::here('02-data', '05-for_analysis', 'fwddiams_training_data.rds'))$plot_id[
              readRDS(here::here('02-data', '05-for_analysis', 'fwddiams_training_data.rds'))$transect_id]],
        plot_id = 
          readRDS(here::here('02-data', '05-for_analysis', 'fwddiams_training_data.rds'))$plot_id[
              readRDS(here::here('02-data', '05-for_analysis', 'fwddiams_training_data.rds'))$transect_id],
        diam_cm = 
          readRDS(here::here('02-data', '05-for_analysis', 'fwddiams_training_data.rds'))$Y),
    aes(x = diam_cm,
        color = 'observed'),
    lwd = 1
  )+
  facet_grid(group_id~.)+
  geom_vline(xintercept = 0.635, lty = 2)+
  geom_vline(xintercept = 2.54, lty = 2)+
  scale_color_manual(values = c('observed' = 'red', 'simulated' = 'blue'))+
  theme_minimal()+
  labs(title = 'FWD diameters retrodiction 1',
       x = 'Diam (cm)', y = 'Probability density')

head(fwddiams_simulatedtallies_data)
# estimate fuel load for each sample


head(fwddiams_simulatedtallies_data)





fwddiams_simulatedtallies_data %>%
  #filter(is.element(draw, sample(x = 1:100, size = 21))) %>%
  group_by(draw, group_id, plot_id, location_id, timelag_class) %>%
  summarise(count = mean(count),
            mgha = mean(mgha)) %>%
  ungroup() %>%
  group_by(draw, group_id, plot_id, timelag_class) %>%
  summarise(count = mean(count),
            mgha = mean(mgha)) %>%
  ungroup() %>%
  mutate(source = 'simulated') %>% 
  select(plot_id, source, timelag_class, count) %>%
  bind_rows(
    fwd_observed %>% 
      group_by(mort, plot_id, location_m, timelag_class) %>%
      summarise(count = mean(count),
                mgha = mean(mgha)) %>%
      ungroup() %>%
      group_by(mort, plot_id, timelag_class) %>%
      summarise(count = mean(count),
                mgha = mean(mgha)) %>%
      ungroup() %>%
      mutate(source = 'observed',
             plot_id = as.integer(factor(plot_id))) %>%
      select(plot_id, source, timelag_class, count)
  ) %>%
  ggplot()+
  geom_histogram(aes(x = count))+
  facet_grid(source~timelag_class, scales = 'free_y')




head(fwd_plots)  



head(fwd_observed)

head(fwddiams_simulatedtallies_data)





head(fwd_observed)

fwd_observed %>%
  mutate(source = 'observed') %>%
  mutate(draw = 4001,
         group_id = ifelse(mort=='low',1,ifelse(mort=='mid',2,3)),
         plot_id = as.integer(factor(plot_id))) %>%
  select(draw, group_id, plot_id, timelag_class, count, mgha, source) %>%
  bind_rows(
    fwddiams_simulatedtallies_data %>%
      select(draw, group_id, plot_id, timelag_class, count, mgha) %>%
      mutate(source = 'simulated')
  ) %>%
  #filter(draw == 1) %>%


  group_by(group_id, plot_id, transect_id, location_m, subsample) %>%
  summarise(mgha = sum(mgha, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(mort, plot_id) %>%
  summarise(mgha = mean(mgha, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(group_id = ifelse(mort=='low',1,
                           ifelse(mort=='mid',2,3)),
         source = 'observed',
         draw = 4001,
         plot_id = as.integer(factor(plot_id))) %>%
  select(draw, group_id, plot_id, mgha, source)




head(fwd_sim_plots)




