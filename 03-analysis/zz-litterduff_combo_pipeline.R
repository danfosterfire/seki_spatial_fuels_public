
library(here)
library(tidyverse)
library(sf)
library(cmdstanr)

#### prepare model data ########################################################

metadata = 
  readRDS(here::here('02-data', '03-clean', 'seki_metadata.rds')) %>%
  left_join(
    readRDS(here::here('02-data', '02-intermediate', 'plot_locations_true.rds')) %>%
      select(mort, plot_id, x_center = x_true, y_center = y_true),
    by = c('plot_id')
  )

litterduff = 
  readRDS(here::here('02-data', '03-clean', 'seki_litterduff.rds')) %>%
  left_join(metadata %>%
              select(plot_id, transect_id, mort)) %>%
  mutate(x_rel = (location_m*sin(az*(pi/180))),
         y_rel = (location_m*cos(az*(pi/180)))) %>%
  filter(!is.na(litter_cm) & !is.na(duff_cm))

coords_litterduff = 
  litterduff %>%
  group_by(x_rel, y_rel) %>%
  summarise() %>%
  ungroup() %>%
  rowid_to_column('location_id')

litterduff = 
  litterduff %>%
  left_join(coords_litterduff) %>%
  mutate(plot_id.i = as.integer(factor(plot_id))) %>%
  mutate(group_id = as.integer(factor(mort, levels = c('low', 'mid', 'high')))) %>%
  rowid_to_column('obs_id')

set.seed(110819)
litterduff.train = 
  litterduff %>%
  group_by(mort) %>%
  slice_sample(prop = 0.9) %>%
  ungroup()
litterduff.valid = 
  litterduff %>%
  filter(!is.element(obs_id, litterduff.train$obs_id))

litterduff.train = 
  litterduff.train %>%
  pivot_longer(cols = c('litter_cm', 'duff_cm'),
               names_to = 'component',
               values_to = 'depth_cm') %>%
  mutate(depth_cm = round(depth_cm, 0),
         component_id = as.integer(factor(component,
                                          levels = c('duff_cm', 'litter_cm'))))

litterduff.valid = 
  litterduff.valid %>%
  pivot_longer(cols = c('litter_cm', 'duff_cm'),
               names_to = 'component',
               values_to = 'depth_cm') %>%
  mutate(depth_cm = round(depth_cm, 0),
         component_id = as.integer(factor(componenet,
                                          levels = c('duff_cm', 'litter_cm'))))



litterduff_training_data = 
  list(N = nrow(litterduff.train),
       C = 2,
       P = length(unique(litterduff$plot_id)),
       L = nrow(coords_litterduff),
       coords = 
         coords_litterduff[,c('x_rel', 'y_rel')] %>%
         as.matrix(),
       component_id = litterduff.train$component_id,
       plot_id = litterduff.train$plot_id.i,
       location_id = litterduff.train$location_id,
       Y = litterduff.train$depth_cm)







#### estimate parameters #######################################################

library(cmdstanr)
library(posterior)
library(bayesplot)

litterduff_model = 
  cmdstanr::cmdstan_model(here::here('03-analysis', 'litterduff_combo.stan'))

litterduff_fit = 
  litterduff_model$sample(
    data = litterduff_training_data,
    seed = 110819,
    chains = 4,
    parallel_chains = 4,
    output_dir = here::here('02-data',
                            '06-results',
                            'real_fits'),
    output_basename = 'litterduff_combo'
  )

litterduff_fit$cmdstan_diagnose()

litterduff_fit$summary()

mcmc_pairs(litterduff_fit$draws(),
           pars = c('beta[1]', 'beta[2]', 'kappa[1]', 'kappa[2]',
                    'sigmaPlot_common', 'sigmaPlot_component[1]','sigmaPlot_component[2]'))


#### posterior retrodictions ###################################################


#### posterior predictions #####################################################


#### prior vs posterior graphs #################################################


#### posterior graphs ##########################################################