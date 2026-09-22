
#### setup #####################################################################

library(here)
library(cmdstanr)

library(tidyverse)
library(posterior)
library(bayesplot)

litter_training_data = readRDS(here::here('02-data', 
                                          '05-for_analysis',
                                          'litter_training_data.rds'))

duff_training_data = readRDS(here::here('02-data',
                                        '05-for_analysis',
                                        'duff_training_data.rds'))


fwdtallies_data = 
  readRDS(here::here('02-data',
                     '05-for_analysis',
                     'fwdtallies_training_data.rds'))


fwddiams_data = 
  readRDS(here::here('02-data',
                     '05-for_analysis',
                     'fwddiams_training_data.rds'))

cwdtallies_data = 
  readRDS(here::here('02-data',
                     '05-for_analysis',
                     'cwdtallies_training_data.rds'))

cwddiams_data = 
  readRDS(here::here('02-data',
                     '05-for_analysis',
                     'cwddiams_training_data.rds'))

vegpa_data = 
  readRDS(here::here('02-data',
                     '05-for_analysis',
                     'vegpa_training_data.rds'))

trees_training_data = 
  readRDS(here::here('02-data',
                     '05-for_analysis',
                     'trees_training_data.rds'))

saplings_training_data = 
  readRDS(here::here('02-data',
                     '05-for_analysis',
                     'saplings_training_data.rds'))

stan_model.litterduff = cmdstan_model(here::here('03-analysis',
                                      'litterduff_depths.stan'))

stan_model.fwddiams = 
  cmdstan_model(here::here('03-analysis',
                           'fwddiams.stan'))

stan_model.cwddiams = 
  cmdstan_model(here::here('03-analysis',
                           'cwddiams.stan'))

stan_model.fwd_tallies = 
  cmdstan_model(here::here('03-analysis',
                           'fwd_tallies.stan'))

stan_model.cwd_tallies = 
  cmdstan_model(here::here('03-analysis', 
                           'cwd_tallies.stan'))

stan_model.veg_presence = 
  cmdstan_model(here::here('03-analysis',
                           'veg_presence.stan'))

stan_model.veg_height = 
  cmdstan_model(here::here('03-analysis',
                           'veg_heights.stan'))

stan_model.stems = 
  cmdstan_model(here::here('03-analysis',
                           'trees_saplings2.stan'))



#### litter ####################################################################

litter_fit = stan_model.litterduff$sample(data = litter_training_data,
                                   parallel_chains = 4,
                                   output_dir = 
                                     here::here('02-data', '06-results', 'real_fits'),
                                   output_basename = 'litter',
                                   seed = 110819,
                                   adapt_delta = 0.99)

litter_fit$save_object(here::here('02-data', '06-results', 'real_fits', 
                                      'litter_fit.rds'))

litter_fit = readRDS(here::here('02-data', '06-results', 'real_fit', 'litter_fit.rds'))

litter_fit$cmdstan_diagnose()

litter_fit$summary(c('beta', 'alpha', 'rho',
                         'sigmaPlot', 'kappa'))

litter_pairs_1 = 
  mcmc_pairs(litter_fit$draws(),
           pars = c('beta[1]', 'alpha[1]', 'rho[1]', 
                    'sigmaPlot[1]', 'kappa[1]'),
           grid_args = list(top = 'Litter'))

litter_pairs_1

ggsave(litter_pairs_1,
       filename = 
         here::here('04-communication', 
                    'figures',
                    'manuscript',
                    'litter_pairs_1.png'),
       height = 6.5, width = 6.5, units = 'in')

litter_pairs_2 = 
  mcmc_pairs(litter_fit$draws(),
           pars = c('beta[2]', 'alpha[2]', 'rho[2]', 
                    'sigmaPlot[2]', 'kappa[2]'),
           grid_args = list(top = 'Litter'))

litter_pairs_2

ggsave(litter_pairs_2,
       filename = 
         here::here('04-communication', 
                    'figures',
                    'manuscript',
                    'litter_pairs_2.png'),
       height = 6.5, width = 6.5, units = 'in')


litter_pairs_3 = 
  mcmc_pairs(litter_fit$draws(),
           pars = c('beta[3]', 'alpha[3]', 'rho[3]', 
                    'sigmaPlot[3]', 'kappa[3]'),
           grid_args = list(top = 'Litter'))

litter_pairs_3

ggsave(litter_pairs_3,
       filename = 
         here::here('04-communication', 
                    'figures',
                    'manuscript',
                    'litter_pairs_3.png'),
       height = 6.5, width = 6.5, units = 'in')

mcmc_pairs(litter_fit$draws(),
           pars = c('beta[1]', 'beta[2]', 'beta[3]'))

mcmc_dens_overlay(litter_fit$draws(variables = c('alpha[1]', 'rho[1]','beta[1]',
                                                     'alpha[2]', 'rho[2]','beta[2]',
                                                     'alpha[3]', 'rho[3]','beta[3]',
                                                     'sigmaPlot[1]',
                                                     'sigmaPlot[2]', 
                                                     'sigmaPlot[3]', 'kappa[1]',
                                                     'kappa[2]', 'kappa[3]')))






#### duff ######################################################################

duff_fit = stan_model.litterduff$sample(data = duff_training_data,
                                   parallel_chains = 4,
                                   output_dir = 
                                     here::here('02-data', '06-results', 'real_fits'),
                                   output_basename = 'duff',
                                   seed = 110819,
                                   adapt_delta = 0.99)


duff_fit$save_object(here::here('02-data', '06-results', 'real_fits', 
                                      'duff_fit.rds'))

duff_fit = readRDS(here::here('02-data', '06-results', 'real_fits','duff_fit.rds'))

duff_fit$cmdstan_diagnose()

duff_fit$summary(c('beta', 'alpha', 'rho',
                         'sigmaPlot', 'kappa')) %>%
  print(n = Inf)

duff_pairs_1 = 
  mcmc_pairs(duff_fit$draws(),
           pars = c('beta[1]', 'alpha[1]', 'rho[1]', 
                    'sigmaPlot[1]', 'kappa[1]'),
           grid_args = list(top = 'Duff'))

duff_pairs_1

ggsave(duff_pairs_1,
       filename = 
         here::here('04-communication', 
                    'figures',
                    'manuscript',
                    'duff_pairs_1.png'),
       height = 6.5, width = 6.5, units = 'in')

duff_pairs_2 = 
  mcmc_pairs(duff_fit$draws(),
           pars = c('beta[2]', 'alpha[2]', 'rho[2]', 
                    'sigmaPlot[2]', 'kappa[2]'),
           grid_args = list(top = 'Duff'))

duff_pairs_2

ggsave(duff_pairs_2,
       filename = 
         here::here('04-communication', 
                    'figures',
                    'manuscript',
                    'duff_pairs_2.png'),
       height = 6.5, width = 6.5, units = 'in')


duff_pairs_3 = 
  mcmc_pairs(duff_fit$draws(),
           pars = c('beta[3]', 'alpha[3]', 'rho[3]', 
                    'sigmaPlot[3]', 'kappa[3]'),
           grid_args = list(top = 'Duff'))

duff_pairs_3

ggsave(duff_pairs_3,
       filename = 
         here::here('04-communication', 
                    'figures',
                    'manuscript',
                    'duff_pairs_3.png'),
       height = 6.5, width = 6.5, units = 'in')

mcmc_pairs(duff_fit$draws(),
           pars = c('beta[1]', 'beta[2]', 'beta[3]'))

mcmc_dens_overlay(duff_fit$draws(variables = c('alpha[1]', 'rho[1]','beta[1]',
                                                     'alpha[2]', 'rho[2]','beta[2]',
                                                     'alpha[3]', 'rho[3]','beta[3]',
                                                     'sigmaPlot[1]',
                                                     'sigmaPlot[2]', 
                                                     'sigmaPlot[3]', 'kappa[1]',
                                                     'kappa[2]', 'kappa[3]')))








#### fwd tallies ###############################################################

fwdtallies_fit = stan_model.fwd_tallies$sample(data = fwdtallies_data,
                                   parallel_chains = 4,
                                   output_dir = 
                                     here::here('02-data', '06-results', 'real_fits'),
                                   output_basename = 'fwdtallies',
                                   seed = 110819,
                                   adapt_delta = 0.99)


fwdtallies_fit$save_object(here::here('02-data', '06-results', 'real_fits', 
                                      'fwdtallies_fit.rds'))

fwdtallies_fit = readRDS(here::here('02-data', '06-results', 'real_fits', 'fwdtallies_fit.rds'))

fwdtallies_fit$cmdstan_diagnose()

fwdtallies_fit$summary(c('beta', 'alpha', 'rho', 'tau',
                         'sigmaPlot', 'kappa')) %>%
  print(n = Inf)

fwdtallies_pairs_1 = 
  mcmc_pairs(fwdtallies_fit$draws(),
           pars = c('beta[1]', 'alpha[1]', 'rho[1]','tau[1]', 
                    'sigmaPlot[1]', 'kappa[1]'))

fwdtallies_pairs_1

ggsave(fwdtallies_pairs_1,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwdtallies_pairs_1.png'),
       height = 6.5, width = 6.5, units = 'in')

fwdtallies_pairs_2 = 
  mcmc_pairs(fwdtallies_fit$draws(),
           pars = c('beta[2]', 'alpha[2]', 'rho[2]','tau[2]', 
                    'sigmaPlot[2]', 'kappa[2]'))

fwdtallies_pairs_2

ggsave(fwdtallies_pairs_2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwdtallies_pairs_2.png'),
       height = 6.5, width = 6.5, units = 'in')

fwdtallies_pairs_3 = 
  mcmc_pairs(fwdtallies_fit$draws(),
           pars = c('beta[3]', 'alpha[3]', 'rho[3]','tau[3]', 
                    'sigmaPlot[3]', 'kappa[3]'))

fwdtallies_pairs_3

ggsave(fwdtallies_pairs_3,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwdtallies_pairs_3.png'),
       height = 6.5, width = 6.5, units = 'in')



mcmc_pairs(fwdtallies_fit$draws(),
           pars = c('beta[1]', 'beta[2]', 'beta[3]'))

mcmc_dens_overlay(fwdtallies_fit$draws(variables = 
                                    c('alpha[1]', 'rho[1]','beta[1]','tau[1]',
                                      'alpha[2]', 'rho[2]','beta[2]','tau[2]',
                                       'alpha[3]', 'rho[3]','beta[3]','tau[3]',
                                      'sigmaPlot[1]', 'sigmaPlot[2]', 'sigmaPlot[3]', 
                                      'kappa[1]', 'kappa[2]', 'kappa[3]')))



#### fwd diameters #############################################################

fwddiams_fit = stan_model.fwddiams$sample(data = fwddiams_data,
                                    parallel_chains = 4,
                                    output_dir = 
                                      here::here('02-data', '06-results', 'real_fits'),
                                    output_basename = 'fwddiams',
                                    seed = 110819,
                                    adapt_delta = 0.9)

fwddiams_fit$save_object(here::here('02-data', '06-results', 'real_fits',
                                 'fwddiams_fit.rds'))

fwddiams_fit = readRDS(here::here('02-data', '06-results', 'real_fits', 'fwddiams_fit.rds'))

fwddiams_fit$cmdstan_diagnose()

fwddiams_fit$summary(c('intercept', 'phi', 'SD_plot', 'SD_trans'))

fwddiams_pairs_1 = 
  mcmc_pairs(fwddiams_fit$draws(),
           np = nuts_params(fwddiams_fit),
           pars = c('intercept[1]', 'SD_plot[1]', 'phi[1]', 'SD_trans[1]'))

ggsave(fwddiams_pairs_1,
       filename = here::here('04-communication',
                  'figures',
                  'fwddiams_pairs_1.png'),
       height = 6.5, width = 6.5, units = 'in')

fwddiams_pairs_2 = 
  mcmc_pairs(fwddiams_fit$draws(),
           np = nuts_params(fwddiams_fit),
           pars = c('intercept[2]', 'SD_plot[2]', 'phi[2]', 'SD_trans[2]'))


ggsave(fwddiams_pairs_2,
       filename = here::here('04-communication',
                  'figures',
                  'fwddiams_pairs_2.png'),
       height = 6.5, width = 6.5, units = 'in')

fwddiams_pairs_3 = 
  mcmc_pairs(fwddiams_fit$draws(),
           np = nuts_params(fwddiams_fit),
           pars = c('intercept[3]', 'SD_plot[3]', 'phi[3]', 'SD_trans[3]'))


ggsave(fwddiams_pairs_3,
       filename = here::here('04-communication',
                  'figures',
                  'fwddiams_pairs_3.png'),
       height = 6.5, width = 6.5, units = 'in')

fwddiams_pairs_1
fwddiams_pairs_2
fwddiams_pairs_3



#### cwd tallies ###############################################################

cwdtallies_fit = 
  stan_model.cwd_tallies$sample(data = cwdtallies_data,
                                parallel_chains = 4,
                                output_dir = here::here('02-data', '06-results', 'real_fits'),
                                output_basename = 'cwdtallies',
                                seed = 112188,
                                adapt_delta = 0.8)

cwdtallies_fit$save_object(here::here('02-data', '06-results', 'real_fits', 
                                      'cwdtallies_fit.rds'))

cwdtallies_fit = readRDS(here::here('02-data', '06-results', 'real_fits', 'cwdtallies_fit.rds'))


cwdtallies_fit$cmdstan_diagnose()

cwdtallies_fit$summary(c('beta', 'alpha', 'rho', 
                         'sigmaPlot', 'sigmaTrans')) %>%
  print(n = Inf)

cwdtallies_pairs_1 = 
  mcmc_pairs(cwdtallies_fit$draws(),
           pars = c('beta[1]', 'alpha[1]', 'rho[1]',
                    'sigmaPlot[1]', 'sigmaTrans[1]'))

cwdtallies_pairs_1

ggsave(cwdtallies_pairs_1,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'cwdtallies_pairs_1.png'),
       height = 6.5, width = 6.5, units = 'in')

cwdtallies_pairs_2 = 
  mcmc_pairs(cwdtallies_fit$draws(),
           pars = c('beta[2]', 'alpha[2]', 'rho[2]',
                    'sigmaPlot[2]', 'sigmaTrans[2]'))

cwdtallies_pairs_2

ggsave(cwdtallies_pairs_2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'cwdtallies_pairs_2.png'),
       height = 6.5, width = 6.5, units = 'in')

cwdtallies_pairs_3 = 
  mcmc_pairs(cwdtallies_fit$draws(),
           pars = c('beta[3]', 'alpha[3]', 'rho[3]',
                    'sigmaPlot[3]', 'sigmaTrans[3]'))

cwdtallies_pairs_3

ggsave(cwdtallies_pairs_3,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'cwdtallies_pairs_3.png'),
       height = 6.5, width = 6.5, units = 'in')




#### cwd diams #################################################################


cwddiams_fit = stan_model.cwddiams$sample(data = cwddiams_data,
                                    parallel_chains = 4,
                                    output_dir = 
                                      here::here('02-data', '06-results', 'real_fits'),
                                    output_basename = 'cwddiams',
                                    seed = 110819,
                                    adapt_delta = 0.8)


cwddiams_fit$save_object(here::here('02-data', '06-results', 'real_fits',
                                 'cwddiams_fit.rds'))

cwddiams_fit = readRDS(here::here('02-data', '06-results', 'real_fits', 'cwddiams_fit.rds'))

cwddiams_fit$cmdstan_diagnose()

cwddiams_fit$summary(c('mu', 'phi'))

cwddiams_pairs_1 = 
  mcmc_pairs(cwddiams_fit$draws(),
           np = nuts_params(cwddiams_fit),
           pars = c('mu[1]', 'phi[1]', 'mu[2]', 'phi[2]', 'mu[3]', 'phi[3]'))

cwddiams_pairs_1

ggsave(cwddiams_pairs_1,
       filename = here::here('04-communication',
                  'figures',
                  'cwddiams_pairs_1.png'),
       height = 6.5, width = 6.5, units = 'in')





#### veg presence ##############################################################


vegpa_fit = stan_model.veg_presence$sample(data = vegpa_data,
                                    parallel_chains = 4,
                                    output_dir = 
                                      here::here('02-data', '06-results', 'real_fits'),
                                    output_basename = 'vegpa',
                                    seed = 110819,
                                    adapt_delta = 0.8)

vegpa_fit$save_object(here::here('02-data', '06-results', 'real_fits',
                                 'vegpa_fit.rds'))

vegpa_fit = readRDS(here::here('02-data', '06-results', 'real_fits', 'vegpa_fit.rds'))

vegpa_fit$cmdstan_diagnose()

vegpa_fit$summary(variables = c('beta', 'alpha', 'rho', 'sigmaPlot')) %>%
  print(n = Inf)

vegpa_pairs_1 = 
  mcmc_pairs(vegpa_fit$draws(
  variables = c('beta[1]', 'alpha[1]', 'rho[1]', 'sigmaPlot[1]')
  ))


vegpa_pairs_2 = 
  mcmc_pairs(vegpa_fit$draws(
  variables = c('beta[2]', 'alpha[2]', 'rho[2]', 'sigmaPlot[2]')
  ))


vegpa_pairs_3 = 
  mcmc_pairs(vegpa_fit$draws(
  variables = c('beta[3]', 'alpha[3]', 'rho[3]', 'sigmaPlot[3]')
  ))

mcmc_dens_overlay(vegpa_fit$draws(
  variables = c('beta', 'alpha', 'rho', 'sigmaPlot')
))


ggsave(vegpa_pairs_1,
       filename = here::here('04-communication',
                  'figures',
                  'vegpa_pairs_1.png'),
       height = 6.5, width = 6.5, units = 'in')


ggsave(vegpa_pairs_2,
       filename = here::here('04-communication',
                  'figures',
                  'vegpa_pairs_2.png'),
       height = 6.5, width = 6.5, units = 'in')


ggsave(vegpa_pairs_3,
       filename = here::here('04-communication',
                  'figures',
                  'vegpa_pairs_3.png'),
       height = 6.5, width = 6.5, units = 'in')




#### trees and snags ###########################################################


trees_fit = stan_model.stems$sample(data = trees_training_data,
                                   parallel_chains = 4,
                                   output_dir = 
                                     here::here('02-data', '06-results', 'real_fits'),
                                   output_basename = 'trees',
                                   seed = 110819,
                                   adapt_delta = 0.8)


trees_fit$save_object(here::here('02-data', '06-results', 'real_fits', 
                                      'trees_fit.rds'))

trees_fit = readRDS(here::here('02-data', '06-results', 'real_fits', 'trees_fit.rds'))

trees_fit$cmdstan_diagnose()

trees_fit$summary(c('beta', 'alpha', 'rho',
                         'sigmaPlot')) %>%
  print(n = Inf)

trees_pairs_1 = 
  mcmc_pairs(trees_fit$draws(),
           pars = c('beta[1]', 'alpha[1]', 'rho[1]', 
                    'sigmaPlot[1]'),
           grid_args = list(top = 'trees'))

trees_pairs_1

ggsave(trees_pairs_1,
       filename = 
         here::here('04-communication', 
                    'figures',
                    'manuscript',
                    'trees_pairs_1.png'),
       height = 6.5, width = 6.5, units = 'in')

trees_pairs_2 = 
  mcmc_pairs(trees_fit$draws(),
           pars = c('beta[2]', 'alpha[2]', 'rho[2]', 
                    'sigmaPlot[2]'),
           grid_args = list(top = 'trees'))

trees_pairs_2

ggsave(trees_pairs_2,
       filename = 
         here::here('04-communication', 
                    'figures',
                    'manuscript',
                    'trees_pairs_2.png'),
       height = 6.5, width = 6.5, units = 'in')


trees_pairs_3 = 
  mcmc_pairs(trees_fit$draws(),
           pars = c('beta[3]', 'alpha[3]', 'rho[3]', 
                    'sigmaPlot[3]'),
           grid_args = list(top = 'trees'))

trees_pairs_3

ggsave(trees_pairs_3,
       filename = 
         here::here('04-communication', 
                    'figures',
                    'manuscript',
                    'trees_pairs_3.png'),
       height = 6.5, width = 6.5, units = 'in')

mcmc_pairs(trees_fit$draws(),
           pars = c('beta[1]', 'beta[2]', 'beta[3]'))

mcmc_dens_overlay(trees_fit$draws(variables = c('alpha[1]', 'rho[1]','beta[1]',
                                                     'alpha[2]', 'rho[2]','beta[2]',
                                                     'alpha[3]', 'rho[3]','beta[3]',
                                                     'sigmaPlot[1]',
                                                     'sigmaPlot[2]', 
                                                     'sigmaPlot[3]')))



#### saplings (incl dead) ######################################################

saplings_fit = stan_model.stems$sample(data = saplings_training_data,
                                   parallel_chains = 4,
                                   output_dir = 
                                     here::here('02-data', '06-results', 'real_fits'),
                                   output_basename = 'saplings',
                                   seed = 110819,
                                   adapt_delta = 0.8)

saplings_fit$save_object(here::here('02-data', '06-results', 'real_fits', 
                                      'saplings_fit.rds'))

saplings_fit = readRDS(here::here('02-data', '06-results', 'real_fits', 'saplings_fit.rds'))

saplings_fit$cmdstan_diagnose()

saplings_fit$summary(c('beta', 'alpha', 'rho',
                         'sigmaPlot')) %>%
  print(n = Inf)

saplings_pairs_1 = 
  mcmc_pairs(saplings_fit$draws(),
           pars = c('beta[1]', 'alpha[1]', 'rho[1]', 
                    'sigmaPlot[1]'),
           grid_args = list(top = 'saplings'))

saplings_pairs_1

ggsave(saplings_pairs_1,
       filename = 
         here::here('04-communication', 
                    'figures',
                    'manuscript',
                    'saplings_pairs_1.png'),
       height = 6.5, width = 6.5, units = 'in')

saplings_pairs_2 = 
  mcmc_pairs(saplings_fit$draws(),
           pars = c('beta[2]', 'alpha[2]', 'rho[2]', 
                    'sigmaPlot[2]'),
           grid_args = list(top = 'saplings'))

saplings_pairs_2

ggsave(saplings_pairs_2,
       filename = 
         here::here('04-communication', 
                    'figures',
                    'manuscript',
                    'saplings_pairs_2.png'),
       height = 6.5, width = 6.5, units = 'in')


saplings_pairs_3 = 
  mcmc_pairs(saplings_fit$draws(),
           pars = c('beta[3]', 'alpha[3]', 'rho[3]', 
                    'sigmaPlot[3]'),
           grid_args = list(top = 'saplings'))

saplings_pairs_3

ggsave(saplings_pairs_3,
       filename = 
         here::here('04-communication', 
                    'figures',
                    'manuscript',
                    'saplings_pairs_3.png'),
       height = 6.5, width = 6.5, units = 'in')

mcmc_pairs(saplings_fit$draws(),
           pars = c('beta[1]', 'beta[2]', 'beta[3]'))

mcmc_dens_overlay(saplings_fit$draws(variables = c('alpha[1]', 'rho[1]','beta[1]',
                                                     'alpha[2]', 'rho[2]','beta[2]',
                                                     'alpha[3]', 'rho[3]','beta[3]',
                                                     'sigmaPlot[1]',
                                                     'sigmaPlot[2]', 
                                                     'sigmaPlot[3]')))


