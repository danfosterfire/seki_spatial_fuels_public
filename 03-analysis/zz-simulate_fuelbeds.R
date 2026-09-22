
#### setup #####################################################################

library(here)
library(posterior)
library(cmdstanr)
library(tidyverse)

# load fitted models and extract median parameter values,
litterduff_params = 
  readRDS(here::here('02-data',
                      '06-results',
                      'real_fits',
                      'litterduff_fit.rds'))$draws(
                        variables = 
                          c('beta[1]', 'beta[2]', 'beta[3]',
                            'alpha[1]', 'alpha[2]', 'alpha[3]',
                            'rho[1]', 'rho[2]', 'rho[3]',
                            'sigmaPlot[1]','sigmaPlot[2]','sigmaPlot[3]', 
                            'kappa')
                      ) %>%
  as_draws_df() %>%
  summarise_all(median)


fwddiams_params = 
  readRDS(here::here('02-data', 
                     '06-results', 
                     'real_fits',
                     'fwddiams_fit.rds'))$draws(
                    variables = c('beta[1]', 'beta[2]', 'beta[3]',
                                  'alpha[1]', 'alpha[2]', 'alpha[3]')
                     ) %>%
  as_draws_df() %>%
  summarise_all(median)

fwd1to100h_params = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'fwd1to100h_fit.rds'))$draws(
                       variables = c('beta[1]', 'beta[2]', 'beta[3]',
                                     'alpha[1]', 'alpha[2]', 'alpha[3]',
                                     'rho[1]', 'rho[2]', 'rho[3]',
                                     'tau[1]', 'tau[2]', 'tau[3]',
                                     'sigmaPlot[1]','sigmaPlot[2]','sigmaPlot[3]', 
                                     'kappa')
                     ) %>%
  as_draws_df() %>%
  summarise_all(median)
  


veg_pa_params = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'veg_pa_fit.rds'))$draws(
                       variables = c('beta[1]', 'beta[2]', 'beta[3]',
                                     'alpha[1]', 'alpha[2]', 'alpha[3]',
                                     'rho[1]', 'rho[2]', 'rho[3]',
                                     'sigmaPlot[1]', 'sigmaPlot[2]',
                                     'sigmaPlot[3]')
                     ) %>%
  as_draws_df() %>%
  summarise_all(median)

veg_height_params = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'veg_height_fit.rds'))$draws(
                       variables = c('beta[1]', 'beta[2]', 'beta[3]',
                                     'alpha[1]', 'alpha[2]', 'alpha[3]',
                                     'rho[1]', 'rho[2]', 'rho[3]',
                                     'sigmaPlot[1]', 'sigmaPlot[2]',
                                     'sigmaPlot[3]')
                     ) %>%
  as_draws_df() %>%
  summarise_all(median)



cwdtallies_params = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'cwd_fit.rds'))$draws(
                       variables = c('beta[1]', 'beta[2]', 'beta[3]',
                                     'alpha[1]', 'alpha[2]', 'alpha[3]',
                                     'rho[1]', 'rho[2]', 'rho[3]',
                                     'sigmaPlot[1]', 'sigmaPlot[2]', 
                                     'sigmaPlot[3]', 'sigmaTransect')
                     ) %>%
  as_draws_df() %>%
  summarise_all(median)

cwddiams_params = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'cwddiams_fit.rds'))$draws(
                       variables = c('alpha[1]', 'alpha[2]', 'alpha[3]',
                                     'beta[1]', 'beta[2]', 'beta[3]')
                     ) %>%
  as_draws_df() %>%
  summarise_all(median)

trees_params = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'trees_fit.rds'))$draws(
                       variables = c('beta[1]', 'beta[2]', 'beta[3]',
                                     'kappa')
                     ) %>%
  as_draws_df() %>%
  summarise_all(median)

sapling_params = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'saplings_fit.rds'))$draws(
                       variables = c('beta[1]', 'beta[2]', 'beta[3]',
                                     'kappa')
                     ) %>%
  as_draws_df() %>%
  summarise_all(median)



#### simulate observations: litterduff  ########################################
source(here::here('00-R', 'simulate_observations.R'))

# NOTE: WANT TO HAVE DIFFERENT SEEDS FOR EACH SIMULATION, OTHERWISE WE MAY 
# BE INDUCING CORRELATION BETWEEN THEM BY TAKING IDENTICAL RANDOM EFFECTS 
# DRAWS (E.G. THE SAME PLOT EFFECTS VECTOR) ACROSS 'INDEPENDENT' SIMULATIONS
set.seed(110819)

litterduff_obs = 
  simulate_obs_litterduff(
    seed = NULL,
    G = 3,
    L = 3721,
    P = 21,
    locations = 
      expand.grid(x_coord = seq(from = 0, to = 30, by = 0.5),
                  y_coord = seq(from = 0, to = 30, by = 0.5)),
    beta = 
      c(litterduff_params$`beta[1]`,
        litterduff_params$`beta[2]`,
        litterduff_params$`beta[3]`),
    
    sigmaPlot = 
      c(litterduff_params$`sigmaPlot[1]`,
        litterduff_params$`sigmaPlot[2]`,
        litterduff_params$`sigmaPlot[3]`),
    
    kappa = litterduff_params$kappa,
    
    alpha = 
      c(litterduff_params$`alpha[1]`,
        litterduff_params$`alpha[2]`,
        litterduff_params$`alpha[3]`),
    
    rho = 
      c(litterduff_params$`rho[1]`,
        litterduff_params$`rho[2]`,
        litterduff_params$`rho[3]`)
    
  )

ggplot(data = 
         litterduff_obs,
       aes(x = x_coord, y = y_coord, fill = litterduff_cm))+
  geom_tile()+
  scale_fill_viridis_c()+
  facet_wrap(~plot_id, labeller = 'label_both')+
  theme_minimal()

#### simulate observations: FWD  ###############################################

fwd1to100h_obs = 
  simulate_obs_fwd_tallies(
    seed = NULL,
    G = 3,
    L = 961,
    P = 21,
    locations = 
      expand.grid(x_coord = seq(from = 0, to = 30, by = 1),
                  y_coord = seq(from = 0, to = 30, by = 1)),
    beta = 
      c(fwd1to100h_params$`beta[1]`,
        fwd1to100h_params$`beta[2]`,
        fwd1to100h_params$`beta[3]`),
    
    sigmaPlot = 
      c(fwd1to100h_params$`sigmaPlot[1]`,
        fwd1to100h_params$`sigmaPlot[2]`,
        fwd1to100h_params$`sigmaPlot[3]`),
    
    kappa = fwd1to100h_params$kappa,
    
    alpha = 
      c(fwd1to100h_params$`alpha[1]`,
        fwd1to100h_params$`alpha[2]`,
        fwd1to100h_params$`alpha[3]`),
    
    rho = 
      c(fwd1to100h_params$`rho[1]`,
        fwd1to100h_params$`rho[2]`,
        fwd1to100h_params$`rho[3]`),
    
    tau = 
      c(fwd1to100h_params$`tau[1]`,
        fwd1to100h_params$`tau[2]`,
        fwd1to100h_params$`tau[3]`)
    
  )


ggplot(data = 
         fwd1to100h_obs,
       aes(x = x_coord, y = y_coord, fill = woody_tally))+
  geom_tile()+
  scale_fill_viridis_c()+
  facet_wrap(~plot_id, labeller = 'label_both')+
  theme_minimal()+
  coord_fixed()


fwd_diams_obs = 
  simulate_obs_woody_diams(seed = NULL,
                           woody_tallies = fwd1to100h_obs,
                           beta = 
                             c(fwddiams_params$`beta[1]`,
                               fwddiams_params$`beta[2]`,
                               fwddiams_params$`beta[3]`),
                           alpha = 
                             c(fwddiams_params$`alpha[1]`,
                               fwddiams_params$`alpha[2]`,
                               fwddiams_params$`alpha[3]`),
                         lb = 0,
                         ub = 7.62)

summary(fwd_diams_obs)

#### simulate observations: CWD ################################################

cwd_obs = 
  simulate_obs_cwd_tallies(
    seed = NULL,
    G = 3,
    L = 961,
    #L = 120,
    P = 21,
    locations = 
      expand.grid(x_coord = seq(from = 0, to = 30, by = 1),
                  y_coord = seq(from = 0, to = 30, by = 1)),
      #data.frame(x_coord = c(seq(-30, 29, 1), rep(0, 60)),
      #           y_coord = c(rep(0,60), seq(-30, 29, 1))),
    beta = 
      c(cwdtallies_params$`beta[1]`,
        cwdtallies_params$`beta[2]`,
        cwdtallies_params$`beta[3]`),
    
    sigmaPlot = 
      c(cwdtallies_params$`sigmaPlot[1]`,
        cwdtallies_params$`sigmaPlot[2]`,
        cwdtallies_params$`sigmaPlot[3]`),
    
    alpha = 
      c(cwdtallies_params$`alpha[1]`,
        cwdtallies_params$`alpha[2]`,
        cwdtallies_params$`alpha[3]`),
    
    rho = 
      c(cwdtallies_params$`rho[1]`,
        cwdtallies_params$`rho[2]`,
        cwdtallies_params$`rho[3]`)
  )

ggplot(data = 
         cwd_obs,
       aes(x = x_coord, y = y_coord, fill = woody_tally))+
  geom_tile()+
  scale_fill_viridis_c()+
  facet_wrap(~plot_id, labeller = 'label_both')+
  theme_minimal()+
  coord_fixed()

cwd_diams_obs = 
  simulate_obs_woody_diams(seed = NULL,
                           woody_tallies = cwd_obs,
                           beta = 
                             c(cwddiams_params$`beta[1]`,
                               cwddiams_params$`beta[2]`,
                               cwddiams_params$`beta[3]`),
                           alpha = 
                             c(cwddiams_params$`alpha[1]`,
                               cwddiams_params$`alpha[2]`,
                               cwddiams_params$`alpha[3]`),
                         lb = 7.62,
                         ub = 110)

summary(cwd_diams_obs)

#### simulate observations: veg ################################################

# resolution of the grid such that each pixel is 1 individual plant of about 
# the middle (not really mean) size,
# with average ln(crown_diameter) = 4.127 -> 0.62m across all samples (all species), 
# so average crown area of ~0.302m2, or a pixel resolution of ~0.549m =~ 0.5m
# note that in the mcginnis 2010 paper I'm relying on, they didn't take a 
# random sample of individuals, but resampled to get wider coverage of 
# the existing size classes, so this dimension is really like a midpoint
veg_pixel_size = 0.5
#veg_pixel_size = 0.62

veg_pa_obs = 
  simulate_obs_veg_presence(
    seed = NULL,
    G = 3,
    M = 14,
    L = 3721,
    #L = 2401,
    P = 21,
    locations = 
 
      # height rather than an average
      expand.grid(x_coord = seq(from = 0, to = 30, by = veg_pixel_size),
                  y_coord = seq(from = 0, to = 30, by = veg_pixel_size)),
    
    beta = 
      c(veg_pa_params$`beta[1]`,
        veg_pa_params$`beta[2]`,
        veg_pa_params$`beta[3]`),
    sigmaPlot = 
      c(veg_pa_params$`sigmaPlot[1]`,
        veg_pa_params$`sigmaPlot[2]`,
        veg_pa_params$`sigmaPlot[3]`),
    alpha = 
      c(veg_pa_params$`alpha[1]`,
        veg_pa_params$`alpha[2]`,
        veg_pa_params$`alpha[3]`),
    rho = 
      c(veg_pa_params$`rho[1]`,
        veg_pa_params$`rho[2]`,
        veg_pa_params$`rho[3]`)
  )

ggplot(data = 
         veg_pa_obs,
       aes(x = x_coord, y = y_coord, fill = veg_present))+
  geom_tile()+
  scale_fill_viridis_d(begin = 0.05, end = 0.85)+
  facet_wrap(~plot_id, labeller = 'label_both')+
  theme_minimal()+
  coord_fixed()

veg_height_obs = 
  simulate_obs_veg_height(
    seed = 110819,
    G = 3,
    M = 14,
    L = 3721,
    P = 21,
    locations = 
 
      # height rather than an average
      expand.grid(x_coord = seq(from = 0, to = 30, by = veg_pixel_size),
                  y_coord = seq(from = 0, to = 30, by = veg_pixel_size)),
    
    beta = 
      c(veg_height_params$`beta[1]`,
        veg_height_params$`beta[2]`,
        veg_height_params$`beta[3]`),
    sigmaPlot = 
      c(veg_height_params$`sigmaPlot[1]`,
        veg_height_params$`sigmaPlot[2]`,
        veg_height_params$`sigmaPlot[3]`),
    alpha = 
      c(veg_height_params$`alpha[1]`,
        veg_height_params$`alpha[2]`,
        veg_height_params$`alpha[3]`),
    rho = 
      c(veg_height_params$`rho[1]`,
        veg_height_params$`rho[2]`,
        veg_height_params$`rho[3]`)
  )

# set the heights for absences to NA
veg_obs = 
  veg_pa_obs %>%
  left_join(veg_height_obs) %>%
  mutate(veg_height = ifelse(veg_present, veg_height, NA),
         veg_N = ifelse(veg_present, veg_N, NA))

head(veg_obs)

ggplot(data = 
         veg_obs,
       aes(x = x_coord, y = y_coord, fill = veg_height))+
  geom_tile()+
  scale_fill_viridis_c()+
  facet_wrap(~plot_id, labeller = 'label_both')+
  theme_minimal()+
  coord_fixed()



#### simulate observations: trees ##############################################

tree_counts = 
  simulate_stem_counts(
    seed = NULL,
    P = 21,
    G = 3,
    beta = c(trees_params$`beta[1]`,trees_params$`beta[2]`,trees_params$`beta[3]`),
    kappa = trees_params$kappa,
    plot_area_ha = (30*30)/10000
  )

tree_obs = 
  simulate_stem_maps(
    seed = NULL,
    stem_counts = tree_counts,
    treelist = 
      readRDS(here::here('02-data',
                                  '03-clean',
                                  'seki_trees.rds')) %>%
      left_join(readRDS(here::here('02-data',
                                   '04-geolocated',
                                   'seki_metadata.rds')) %>%
                  select(plot_id, mort)) %>%
      mutate(group_id = 
               ifelse(mort == 'low',
                      1,
                      ifelse(mort == 'mid',
                             2,
                             3))) %>%
      select(group_id, status, spp, dbh_cm, height_m, htcb_m, decay_class) %>%
      filter(dbh_cm >= 11.4),
    plot_window = list(x_min = 0, x_max = 30, y_min = 0, y_max = 30)
  )


tree_obs

#### simulate observations: saplings ###########################################


sapling_counts = 
  simulate_stem_counts(
    seed = NULL,
    P = 21,
    G = 3,
    beta = c(sapling_params$`beta[1]`,sapling_params$`beta[2]`,sapling_params$`beta[3]`),
    kappa = sapling_params$kappa,
    plot_area_ha = (30*30)/10000
  )

sapling_obs = 
  simulate_stem_maps(
    seed = NULL,
    stem_counts = sapling_counts,
    treelist = 
      readRDS(here::here('02-data',
                                  '03-clean',
                                  'seki_trees.rds')) %>%
      left_join(readRDS(here::here('02-data',
                                   '04-geolocated',
                                   'seki_metadata.rds')) %>%
                  select(plot_id, mort)) %>%
      mutate(group_id = 
               ifelse(mort == 'low',
                      1,
                      ifelse(mort == 'mid',
                             2,
                             3))) %>%
      select(group_id, status, spp, dbh_cm, height_m, htcb_m, decay_class) %>%
      filter(dbh_cm < 11.4),
    plot_window = list(x_min = 0, x_max = 30, y_min = 0, y_max = 30)
  )


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
litterduff_coeffs = 
  overstory_composition %>%
  left_join(litterduff_coeffs) %>%
  mutate(weighted = pBA * litterduff_coeff) %>%
  group_by(mort) %>%
  summarise(weighted_coeff = sum(weighted)) %>%
  ungroup()

litterduff_coeffs

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

#### biomass conversion: litterduff ############################################


litterduff_loads = 
  litterduff_obs %>%
  select(group_id, plot_id, x_coord, y_coord, litterduff_cm) %>%
  
  left_join(data.frame(group_id = c(1, 2, 3),
                       mort = c('low', 'mid', 'high'))) %>%
  
  
  # get the stimated litter+duff load for each observation
  left_join(litterduff_coeffs) %>%
  mutate(litterduff_kgm2 = litterduff_cm * weighted_coeff,
         litterduff_mgha = litterduff_kgm2 * 10) %>%
  
  # standardize columns
  select(mort, group_id, plot_id, x_coord, y_coord, 
         litterduff_cm, litterduff_mgha)



ggplot(data = 
         litterduff_loads %>%
         mutate(plot_id = floor((plot_id-1)/3)+1,
                mort = factor(mort, levels = c('low', 'mid', 'high'))),
       aes(x = x_coord, y = y_coord, fill = litterduff_mgha))+
  geom_tile()+
  scale_fill_viridis_c()+
  facet_grid(mort~plot_id, labeller = 'label_both')+
  theme_minimal()+
  coord_fixed()


ggplot(data = litterduff_loads,
       aes(x = litterduff_mgha))+
  geom_histogram()





#### biomass conversion: FWD ###################################################

fwd_loading_diams = 
  fwd_diams_obs %>%
  mutate(timelag_class = 
           ifelse(diam_cm > 0 & diam_cm < 0.64,
                  '1h',
                  ifelse(diam_cm >= 0.64 & diam_cm < 2.54,
                         '10h',
                         '100h')),
         count = 1) %>%
  group_by(group_id, plot_id, x_coord, y_coord, timelag_class) %>%
  summarise(count = sum(count)) %>%
  ungroup() %>%
  complete(nesting(group_id, plot_id), 
           x_coord, y_coord, timelag_class) %>%
  mutate(count = ifelse(is.na(count),0,count)) %>%
  
  left_join(data.frame(group_id = c(1, 2, 3),
                       mort = c('low', 'mid', 'high'))) %>%
  left_join(QMDcm) %>%
  left_join(SEC) %>%
  left_join(SG) %>%
  mutate(slp_c = 1,
         transect_length_m = 1,
         k = 1.234) %>%
  
  mutate(mgha = (k*weighted_qmd*weighted_sec*slp_c*weighted_sg*count)/
                    transect_length_m) %>%
  select(mort, plot_id, x_coord, y_coord, timelag_class, count, mgha)




#### biomass conversion: cwd ###################################################

length(seq(0,30,1))**2*length(unique(cwd_diams_obs$plot_id))


cwd_loading = 
  cwd_diams_obs %>%
  mutate(timelag_class = '1000s') %>%
  group_by(group_id, plot_id, timelag_class, x_coord, y_coord) %>%
  summarise(ssd_cm2 = sum(diam_cm**2)) %>%
  ungroup() %>%
  
  # fill in the zeroes
  bind_rows(
    data.frame(x_coord = rep(c(seq(-30, 29, 1), rep(0, 60)), times = 21),
               y_coord = rep(c(rep(0,60), seq(-30, 29, 1)), times = 21),
               group_id = rep(x = c(1, 2, 3), each = 120*7),
               plot_id = rep(x = 1:21, each = 120)) %>%
      mutate(ssd_cm2 = 0,
             timelag_class = '1000s')
  ) %>%
  group_by(group_id, plot_id, timelag_class, x_coord, y_coord) %>%
  summarise(ssd_cm2 = sum(ssd_cm2)) %>%
  ungroup() %>%

# and join in the coefficients 
  left_join(.,
            y = SEC) %>%
  left_join(.,
            y = SG) %>%
  mutate(transect_length_m = 1, # using 1m subtransects
         k = 1.234,
         slp_c = 1) %>%
  
  # and use brown's equation (with the new k const) to estimate fuel load
  mutate(mgha = 
           (k*ssd_cm2*weighted_sec*slp_c*weighted_sg)/
           transect_length_m)

ggplot(cwd_loading,
       aes(x = x_coord, y = y_coord, fill = mgha))+
  geom_tile()+
  facet_wrap(~plot_id)+
  coord_fixed()+
  scale_fill_viridis_c()






#### biomass conversion: understory veg ########################################

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
head(veg_obs)

veg_obs = 
  veg_obs %>%
  mutate(
    crown_vol_m3 = (veg_pixel_size**2)*veg_height,
    total_biomass_g_vol = 
      ifelse(veg_present,
             exp(7.278+(0.812*log(crown_vol_m3)))*1.27,
             0),
    foliage_biomass_g_vol = 
      ifelse(veg_present,
             exp(5.297+(0.662*log(crown_vol_m3)))*1.65,
             0),
    x1h_biomass_g_vol = 
      ifelse(veg_present,
             exp(6.362+(0.739*log(crown_vol_m3)))*1.25,
             0),
    ln_diameter_cm = log(56),
    ln_height_cm = log(veg_height*100),
    total_biomass_g_dh = 
      ifelse(veg_present,
             exp(-4.658 + (2.078*ln_diameter_cm) + (0.336*ln_height_cm))*1.28,
             0),
    foliage_biomass_g_dh = 
      ifelse(veg_present,
             exp(-4.363 + (1.754*ln_diameter_cm)+(0.193*ln_height_cm))*1.69,
             0),
    x1h_biomass_g_dh = 
      ifelse(veg_present,
             exp(-4.250 + (2.124*ln_diameter_cm)+(-0.004*ln_height_cm))*1.25,
             0)
  )
  

head(veg_obs)





#### biomass conversion: trees and saplings ####################################

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

tree_sapling_biomass = 
  tree_obs %>%
  bind_rows(sapling_obs) %>%
  
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
                 1/0.05,
                 1/((2*30+(2*28))/10000)),
    ba_m2 = pi*((dbh_cm/2/100)**2),
    ba_m2ha = ba_m2 * tph,
    agbiomass_mgha = 
      (agbiomass_kg/1000)*tph
  )

head(tree_sapling_biomass)

#### voxelize understory #######################################################

# convert from pixels to voxels
veg_voxels = 
  veg_obs %>%
  select(group_id, plot_id, x_coord, y_coord, veg_height_m = veg_height,
         total_biomass_g = total_biomass_g_dh, 
         foliage_biomass_g = foliage_biomass_g_dh, 
         x1h_biomass_g = x1h_biomass_g_dh) %>%
  expand(nesting(group_id, plot_id, x_coord, y_coord,
                 veg_height_m, total_biomass_g, foliage_biomass_g, x1h_biomass_g),
         z_coord = seq(from = 0.125, to = 3.325, by = 0.25)) %>%
  
  # divide the total pixel biomass amongst the vertical layers
  mutate(total_biomass_g = 
           # if the vegetation height is above this z coord then biomass is 
           # the total biomass divided by the number of occupied voxels, otherwise 0
           ifelse(veg_height_m > z_coord & !is.na(veg_height_m),
                  total_biomass_g / (veg_height_m/0.25),
                  0),
         
         foliage_biomass_g = 
           ifelse(veg_height_m > z_coord & !is.na(veg_height_m),
                  foliage_biomass_g / (veg_height_m/0.25),
                  0),
         
         x1h_biomass_g = 
           ifelse(veg_height_m > z_coord & !is.na(veg_height_m),
                  x1h_biomass_g / (veg_height_m/0.25),
                  0)
         ) %>%
  
    
  # convert to bulk density (g/m3)
  mutate(total_biomass_gm3 = total_biomass_g / (veg_pixel_size*veg_pixel_size*0.25),
         foliage_biomass_gm3 = foliage_biomass_g / (veg_pixel_size*veg_pixel_size*0.25),
         x1h_biomass_gm3 = x1h_biomass_g / (veg_pixel_size*veg_pixel_size*0.25)) %>%
  
  # clean up intermediate columns for clarity
  select(group_id, plot_id, x_coord, y_coord, z_coord, total_biomass_gm3, 
         foliage_biomass_gm3, x1h_biomass_gm3)


#### voxelize trees and saplings ###############################################

head(tree_sapling_biomass)

summary(tree_sapling_biomass)

#### reality comparison: litterduff ############################################

readRDS(here::here('02-data',
                   '05-for_analysis',
                   'biomass',
                   'seki_litterduff.rds')) %>%
  select(mort, litterduff_mgha) %>%
  mutate(source = 'observed') %>%
  bind_rows(litterduff_loads %>%
              select(mort,litterduff_mgha) %>%
              mutate(source = 'simulated')) %>%
  ggplot(aes(x = litterduff_mgha))+
  geom_histogram()+
  facet_grid(source~mort, scales = 'free_y')


readRDS(here::here('02-data',
                   '05-for_analysis',
                   'biomass',
                   'seki_litterduff.rds')) %>%
  select(mort, litterduff_mgha) %>%
  mutate(source = 'observed') %>%
  bind_rows(litterduff_loads %>%
              select(mort,litterduff_mgha) %>%
              mutate(source = 'simulated')) %>%
  ggplot(aes(x = litterduff_mgha, color = source))+
  geom_density()+
  facet_grid(.~mort, scales = 'free_y')+
  theme_minimal()

#### reality comparison: FWD ###################################################

fwd_loading_diams %>%
  group_by(plot_id, timelag_class) %>%
  summarise(mgha = mean(mgha)) %>%
  ungroup() %>%
  pivot_wider(names_from = 'timelag_class', values_from = 'mgha') %>%
  summary()

fwd_loading_diams %>%
  group_by(plot_id, timelag_class) %>%
  summarise(mgha = mean(mgha)) %>%
  ungroup() %>%
  mutate(source = 'simulated',
         plot_id = as.character(plot_id)) %>%
  bind_rows(
    readRDS(here::here('02-data', 
                   '05-for_analysis', 
                   'biomass',
                   'seki_fwd.rds')) %>%
  pivot_longer(cols = c(count_a1h, count_a10h, count_a100h,
                        count_b1h, count_b10h, count_b100h,
                        mgha_a1h, mgha_a10h, mgha_a100h,
                        mgha_b1h, mgha_b10h, mgha_b100h),
               names_sep = '_',
               names_to = c('response', 'subsampleclass')) %>%
  mutate(subsample = gsub(x = subsampleclass, pattern = '1h|10h|100h',
                          replacement = ''),
         timelag_class = gsub(x = subsampleclass, pattern = 'a|b',
                              replacement = '')) %>%
  pivot_wider(names_from = 'response', values_from = 'value') %>%
    group_by(plot_id, timelag_class) %>%
    summarise(mgha = mean(mgha, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(source = 'real')
  ) %>%
  ggplot(aes(x = mgha, color = source))+
  geom_density()+
  facet_wrap(~timelag_class, scales = 'free')+
  theme_minimal()

fwd_loading_diams %>%
  group_by(plot_id, timelag_class) %>%
  summarise(mgha = mean(mgha)) %>%
  ungroup() %>%
  mutate(source = 'simulated',
         plot_id = as.character(plot_id)) %>%
  bind_rows(
    readRDS(here::here('02-data', 
                   '05-for_analysis', 
                   'biomass',
                   'seki_fwd.rds')) %>%
  pivot_longer(cols = c(count_a1h, count_a10h, count_a100h,
                        count_b1h, count_b10h, count_b100h,
                        mgha_a1h, mgha_a10h, mgha_a100h,
                        mgha_b1h, mgha_b10h, mgha_b100h),
               names_sep = '_',
               names_to = c('response', 'subsampleclass')) %>%
  mutate(subsample = gsub(x = subsampleclass, pattern = '1h|10h|100h',
                          replacement = ''),
         timelag_class = gsub(x = subsampleclass, pattern = 'a|b',
                              replacement = '')) %>%
  pivot_wider(names_from = 'response', values_from = 'value') %>%
    group_by(plot_id, timelag_class) %>%
    summarise(mgha = mean(mgha, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(source = 'real')
  ) %>%
  ggplot(aes(x = mgha, fill = source))+
  geom_histogram()+
  theme_bw()+
  facet_grid(source~timelag_class, scales = 'free')




readRDS(here::here('02-data',
                   '05-for_analysis',
                   'biomass',
                   'seki_fwd.rds')) %>%
  group_by(plot_id) %>%
  summarise(mgha_a1h = mean(mgha_a1h, na.rm = TRUE),
            mgha_b1h = mean(mgha_b1h, na.rm = TRUE),
            mgha_a10h = mean(mgha_a10h, na.rm = TRUE),
            mgha_b10h = mean(mgha_b10h, na.rm = TRUE),
            mgha_a100h = mean(mgha_a100h, na.rm = TRUE),
            mgha_b100h = mean(mgha_b100h, na.rm = TRUE)) %>%
  ungroup() %>%
  summary()

readRDS(here::here('02-data', 
                   '05-for_analysis', 
                   'biomass',
                   'seki_fwd.rds')) %>%
  pivot_longer(cols = c(count_a1h, count_a10h, count_a100h,
                        count_b1h, count_b10h, count_b100h,
                        mgha_a1h, mgha_a10h, mgha_a100h,
                        mgha_b1h, mgha_b10h, mgha_b100h),
               names_sep = '_',
               names_to = c('response', 'subsampleclass')) %>%
  mutate(subsample = gsub(x = subsampleclass, pattern = '1h|10h|100h',
                          replacement = ''),
         timelag_class = gsub(x = subsampleclass, pattern = 'a|b',
                              replacement = '')) %>%
  pivot_wider(names_from = 'response', values_from = 'value') %>%
  select(mort, plot_id, x_coord, y_coord, timelag_class, count, mgha) %>%
  mutate(source = 'real') %>%
  bind_rows(
    fwd_loading_diams %>%
      select(mort, plot_id, x_coord, y_coord, timelag_class, count, mgha) %>%
      mutate(source = 'simulated',
             plot_id = as.character(plot_id))
  ) %>%
  ggplot(aes(x = mgha, fill = source))+
  geom_histogram(position = position_dodge())+
  facet_grid(source~mort+timelag_class, scales = 'free')+
  theme_minimal()

#### reality comparison: CWD ###################################################

cwd_loading %>%
  filter(x_coord == 15 | y_coord == 15) %>%
  group_by(plot_id, timelag_class) %>%
  summarise(mgha = mean(mgha), ssd = mean(ssd_cm2)) %>%
  ungroup() %>%
  summary()


# looks good now
readRDS(here::here('02-data',
                   '05-for_analysis',
                   'biomass',
                   'seki_cwd.rds')) %>%
  group_by(plot_id, timelag_class) %>%
  summarise(mgha = mean(mgha), ssd = mean(ssd_cm2)) %>%
  ungroup() %>%
  summary()

readRDS(here::here('02-data',
                   '05-for_analysis',
                   'biomass',
                   'seki_cwd.rds')) %>%
  group_by(plot_id, timelag_class) %>%
  summarise(mgha = mean(mgha), ssd = mean(ssd_cm2)) %>%
  ungroup() %>%
  mutate(source = 'real') %>%
  bind_rows(
    cwd_loading %>%
      group_by(plot_id, timelag_class) %>%
      summarise(mgha = mean(mgha)) %>%
      ungroup() %>%
      mutate(plot_id = as.character(plot_id),
             source = 'simulated')
  ) %>%
  ggplot(aes(x = mgha, fill = source))+
  geom_histogram()+
  facet_grid(source~.)+
  theme_minimal()


# biomass predictions are too low; is it because of the diameters or the counts?

# check the counts; they look OK
sum(cwd_obs$woody_tally) # 428 intersections on 120m of transect time 21 plots
sum(cwd_obs$woody_tally)/nrow(cwd_obs) # 0.17 intersections per meter of transect

nrow(seki.cwd) # 430 intersections on 120m of transect times 21 plots
nrow(seki.cwd)/(120*21) # 0.17 intersections per meter


test = readRDS(here::here('02-data',
                       '04-geolocated',
                       'seki_cwd.rds')) %>%
      select(plot_id, az, location_m) %>%
      mutate(woody_tally = 1, location_m = floor(location_m)) %>%
      bind_rows(
        readRDS(here::here('02-data',
                           '04-geolocated',
                           'seki_metadata.rds')) %>%
          select(plot_id) %>%
          tidyr::expand(plot_id,
                 az = c(0, 90, 180, 270),
                 location_m = 0:29) %>%
          mutate(woody_tally = 0)
      ) %>%
      group_by(plot_id, az, location_m) %>%
      summarise(woody_tally = sum(woody_tally)) %>%
      ungroup() %>%
      select(plot_id, woody_tally)


cwd_obs %>%
  select(plot_id, woody_tally) %>%
  mutate(source = 'simulated',
         plot_id = paste0(plot_id, 's')) %>%
  bind_rows(test %>%
              mutate(source = 'real')) %>%
  ggplot(aes(x = woody_tally, fill = source))+
  geom_histogram(position = position_dodge())


# check the diameters
head(cwd_diams_obs) 

cwd_diams_obs %>%
  select(group_id, diam_cm) %>%
  mutate(source = 'sim') %>%
  bind_rows(
    data.frame(
      diam_cm = readRDS(here::here('02-data',
                                   '05-for_analysis',
                                   'cwd_diams_data.rds'))$Y,
      group_id = readRDS(here::here('02-data',
                                    '05-for_analysis',
                                    'cwd_diams_data.rds'))$group_id[
                                      readRDS(here::here('02-data',
                                                         '05-for_analysis',
                                                         'cwd_diams_data.rds'))$plot_id
                                    ]
    ) %>%
      mutate(source = 'real')
  ) %>%
  ggplot(aes(x = diam_cm, color = source))+
  geom_density()+
  facet_wrap(~group_id)

cwd_diams_obs %>%
  select(plot_id, diam_cm) %>%
  mutate(source = 'sim',
         plot_id = paste0(as.character(plot_id),'s')) %>%
  bind_rows(
      readRDS(here::here('02-data',
                     '04-geolocated',
                     'seki_cwd.rds')) %>%
      select(plot_id, diam_cm)%>%
      mutate(source = 'real')
  ) %>%
  ggplot(aes(x = diam_cm, color = source))+
  geom_density()

# both diameters and counts look fine what's going on?
seki.cwd = 
  readRDS(here::here('02-data',
                     '04-geolocated',
                     'seki_cwd.rds'))

seki.metadata = readRDS(here::here('02-data',
                                    '04-geolocated',
                                    'seki_metadata.rds'))


# full set of subtransects
subtransects = 
  seki.metadata %>%
               filter(cwd=='1') %>%
               mutate(
                 slp_pcent = 
                   tan(slp_degrees*(pi/180))*100,
                 slp_c = sqrt(1+(slp_pcent/100)^2)) %>%
               select(mort, plot_id, transect_id, az, x_center, y_center, slp_c) %>%
                tidyr::expand(nesting(mort, plot_id, transect_id, az, x_center,
                                      y_center, slp_c),
                              location_m = seq(0, 29, 1)) %>%
               tidyr::expand(nesting(mort, plot_id, transect_id, az, x_center,
                                     y_center, slp_c, location_m),
                      timelag_class = c('1000s', '1000r'))

observed.cwd = 
  seki.cwd %>%
  # tidy up 
  mutate(timelag_class = ifelse(as.numeric(decay_class) <= 2, '1000s', '1000r')) %>%
  select(mort, plot_id, transect_id, az, location_m, 
         timelag_class, diam_cm, comments, x_center, y_center) %>%
  
  # get the sum of squared diameters for each subtransect and timelag class
  # round the location to the nearest 5m midpoint (2.5, 7.5,...27.5m)
  mutate(location_m = floor(location_m)) %>%
  group_by(mort, plot_id, transect_id, az, x_center, y_center,  location_m,
           timelag_class)  %>%
  summarize(ssd_cm2 = sum(diam_cm^2)) %>%
  ungroup() %>%
  right_join(subtransects) %>%
  mutate(ssd_cm2 = ifelse(is.na(ssd_cm2), 0, ssd_cm2)) %>%
  
  # and join in the coefficients 
  left_join(.,
            y = SEC) %>%
  left_join(.,
            y = SG) %>%
  mutate(transect_length_m = 1, # using 5m subtransects
         k = 1.234) %>%
  
  # and use brown's equation (with the new k const) to estimate fuel load
  mutate(mgha = 
           (k*ssd_cm2*weighted_sec*slp_c*weighted_sg)/
           transect_length_m) %>%
  
  # drop some intermediate columns
  select(mort, plot_id, transect_id, az, x_center, y_center, 
         timelag_class, location_m, ssd_cm2, mgha) %>%
  
  # going to lump sound and rotten together to reduce the number of 0s and the 
  # dimensionality of the response, mark will @ me if he wants it changed
  mutate(timelag_class = 'cwd') %>%
  group_by(mort, plot_id, transect_id, az, x_center, y_center, 
           location_m, timelag_class) %>%
  summarise(ssd_cm2 = sum(ssd_cm2),
            mgha = sum(mgha)) %>%
  ungroup()

observed.cwd %>%
  group_by(plot_id) %>%
  summarise(ssd_cm2 = mean(ssd_cm2),
            mgha = mean(mgha)) %>%
  ungroup() %>%
  summary()



#### reality comparison: understory veg ########################################

head(veg_obs)

veg_obs %>%
  group_by(group_id, plot_id) %>%
  mutate(total_biomass_gm2_vol = total_biomass_g_vol / (veg_pixel_size**2),
         total_biomass_mgha_vol = total_biomass_gm2_vol * 0.01,
         total_biomass_gm2_dh = total_biomass_g_dh / (veg_pixel_size**2),
         total_biomass_mgha_dh = total_biomass_gm2_dh * 0.01) %>%
  summarise(total_biomass_mgha_vol = mean(total_biomass_mgha_vol),
            total_biomass_mgha_dh = mean(total_biomass_mgha_dh)) %>%
  ungroup() %>%
  mutate(source = 'simulated') %>%
  bind_rows(
   readRDS(here::here('02-data',
                   '04-geolocated',
                   'seki_metadata.rds')) %>%
  select(mort, plot_id, transect_id) %>%
  expand(nesting(mort, plot_id, transect_id),
         location_m = seq(0.5, 14.5, 1)) %>%
  left_join(readRDS(here::here('02-data',
                               '03-clean',
                               'seki_veg.rds')) %>%
              #filter(!is.element(spp, c('G', 'F'))) %>%
              select(plot_id, transect_id, start_m, end_m, height_m)) %>%
  group_by(mort, plot_id, transect_id, location_m) %>%
  summarise(veg_present = any(location_m>=start_m & location_m<=end_m),
            height_m = max(height_m, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(veg_present = ifelse(is.na(veg_present), FALSE, veg_present),
         height_m = ifelse(is.finite(height_m),height_m, NA)) %>%
  group_by(mort, plot_id) %>%
  summarise(p_cover = mean(veg_present),
            height_m = mean(height_m, na.rm = TRUE)) %>%
  ungroup() %>%
  # individuals_per_ha = (1 ind / 0.302 m2) * 10,000m2 / ha * p_cover
  mutate(ind_crown_area_m2 = pi*(exp(4.127)/200)**2,
         individuals_per_ha = (1/ind_crown_area_m2)*10000*p_cover,
         # uses crown diameter and height, assuming crown diameter
         biomass_g_per_ind_1 = exp(-4.658+2.079*4.127+0.336*log(height_m*100))*1.28,
         ind_crown_vol_m3 = ind_crown_area_m2*height_m,
         
         # uses crown volume
         biomass_g_per_ind_2 = exp(7.278+0.812*log(ind_crown_vol_m3))*1.27,
         total_biomass_mgha_dh = (biomass_g_per_ind_1/1000000)*individuals_per_ha,
         total_biomass_mgha_vol = (biomass_g_per_ind_2/1000000)*individuals_per_ha) %>%
    mutate(source = 'observed',
           plot_id = as.numeric(plot_id),
           group_id = as.integer(factor(mort, levels = c('low', 'mid', 'high')))) %>%
    select(group_id, plot_id, source, total_biomass_mgha_dh, total_biomass_mgha_vol)
  ) %>%
  
  ggplot(aes(x = total_biomass_mgha_vol, fill = source))+
  geom_histogram()+
  facet_grid(source+group_id~.)+
  theme_minimal()


veg_obs %>%
  group_by(group_id, plot_id) %>%
  mutate(total_biomass_gm2_vol = total_biomass_g_vol / (veg_pixel_size**2),
         total_biomass_mgha_vol = total_biomass_gm2_vol * 0.01,
         total_biomass_gm2_dh = total_biomass_g_dh / (veg_pixel_size**2),
         total_biomass_mgha_dh = total_biomass_gm2_dh * 0.01) %>%
  summarise(total_biomass_mgha_vol = mean(total_biomass_mgha_vol),
            total_biomass_mgha_dh = mean(total_biomass_mgha_dh)) %>%
  ungroup() %>%
  mutate(source = 'simulated') %>%
  bind_rows(
   readRDS(here::here('02-data',
                   '04-geolocated',
                   'seki_metadata.rds')) %>%
  select(mort, plot_id, transect_id) %>%
  expand(nesting(mort, plot_id, transect_id),
         location_m = seq(0.5, 14.5, 1)) %>%
  left_join(readRDS(here::here('02-data',
                               '03-clean',
                               'seki_veg.rds')) %>%
              #filter(!is.element(spp, c('G', 'F'))) %>%
              select(plot_id, transect_id, start_m, end_m, height_m)) %>%
  group_by(mort, plot_id, transect_id, location_m) %>%
  summarise(veg_present = any(location_m>=start_m & location_m<=end_m),
            height_m = max(height_m, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(veg_present = ifelse(is.na(veg_present), FALSE, veg_present),
         height_m = ifelse(is.finite(height_m),height_m, NA)) %>%
  group_by(mort, plot_id) %>%
  summarise(p_cover = mean(veg_present),
            height_m = mean(height_m, na.rm = TRUE)) %>%
  ungroup() %>%
  # individuals_per_ha = (1 ind / 0.302 m2) * 10,000m2 / ha * p_cover
  mutate(ind_crown_area_m2 = pi*(exp(4.127)/200)**2,
         individuals_per_ha = (1/ind_crown_area_m2)*10000*p_cover,
         # uses crown diameter and height, assuming crown diameter
         biomass_g_per_ind_1 = exp(-4.658+2.079*4.127+0.336*log(height_m*100))*1.28,
         ind_crown_vol_m3 = ind_crown_area_m2*height_m,
         
         # uses crown volume
         biomass_g_per_ind_2 = exp(7.278+0.812*log(ind_crown_vol_m3))*1.27,
         total_biomass_mgha_dh = (biomass_g_per_ind_1/1000000)*individuals_per_ha,
         total_biomass_mgha_vol = (biomass_g_per_ind_2/1000000)*individuals_per_ha) %>%
    mutate(source = 'observed',
           plot_id = as.numeric(plot_id),
           group_id = as.integer(factor(mort, levels = c('low', 'mid', 'high')))) %>%
    select(group_id, plot_id, source, total_biomass_mgha_dh, total_biomass_mgha_vol)
  ) %>%
  
  ggplot(aes(x = total_biomass_mgha_dh, fill = source))+
  geom_histogram()+
  facet_grid(source+group_id~.)+
  theme_minimal()



#### START HERE AND SHOW BATTLES: BIOMASS UNCERTAINTY ##########################
readRDS(here::here('02-data',
                   '03-clean',
                   'seki_metadata.rds')) %>%
  select(plot_id, transect_id) %>%
  expand(nesting(plot_id, transect_id),
         location_m = seq(0.5, 14.5, 1)) %>%
  left_join(readRDS(here::here('02-data',
                               '03-clean',
                               'seki_veg.rds')) %>%
              #filter(!is.element(spp, c('G', 'F'))) %>%
              select(plot_id, transect_id, start_m, end_m, height_m)) %>%
  group_by(plot_id, transect_id, location_m) %>%
  summarise(veg_present = any(location_m>=start_m & location_m<=end_m),
            height_m = max(height_m, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(veg_present = ifelse(is.na(veg_present), FALSE, veg_present),
         height_m = ifelse(is.finite(height_m),height_m, NA)) %>%
  group_by(plot_id) %>%
  summarise(p_cover = mean(veg_present),
            height_m = mean(height_m, na.rm = TRUE)) %>%
  ungroup() %>%
  # individuals_per_ha = (1 ind / 0.302 m2) * 10,000m2 / ha * p_cover
  mutate(ind_crown_area_m2 = pi*(exp(4.127)/200)**2,
         individuals_per_ha = (1/ind_crown_area_m2)*10000*p_cover,
         # uses crown diameter and height, assuming crown diameter
         biomass_g_per_ind_1 = exp(-4.658+2.079*4.127+0.336*log(height_m*100))*1.28,
         ind_crown_vol_m3 = ind_crown_area_m2*height_m,
         
         # uses crown volume
         biomass_g_per_ind_2 = exp(7.278+0.812*log(ind_crown_vol_m3))*1.27,
         biomass_mgha_1 = (biomass_g_per_ind_1/1000000)*individuals_per_ha,
         biomass_mgha_2 = (biomass_g_per_ind_2/1000000)*individuals_per_ha) %>%
  #select(biomass_mgha_1, biomass_mgha_2) %>%
  summary()
  #pull(biomass_mgha_1) %>%
  #hist(main = 'Mgha biomass (crown diam + heigeht)') 
  #pull(biomass_mgha_2) %>%
  #hist(main = 'Mgha biomass (crown vol.)')


# assuming that each individual has a crown diameter of exp(4.127) cm and some 
# plot-average height, it 
# matters a LOT whether you're estimating the biomass using the diameter+height 
# equation or the volume equation. It looks like the simulations are good, but 
# my sense of shrub biomass is calibrated to the diameter+height equation 
# rather than the crown volume equation. Probably the diameter+height equation 
# is more reliable for our hypothetical tall+skinny plants?


head(veg_voxels)

veg_voxels %>%
  group_by(group_id, plot_id, x_coord, y_coord) %>%
  summarise(total_biomass_gm2 = sum(total_biomass_gm3*0.25)) %>%
  ungroup() %>%
  group_by(group_id, plot_id) %>%
  summarise(total_biomass_gm2 = mean(total_biomass_gm2,na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(total_biomass_Mgha = total_biomass_gm2*0.01) %>%
  pull(total_biomass_Mgha) %>%
  #ggplot(aes(x = total_biomass_Mgha))+
  #geom_histogram()
  hist(main = 'Total shrub biomass (Mg/ha)')
  #summary()
# understory veg biomass looks pretty high to me

head(veg_obs)

veg_obs %>%
  group_by(group_id, plot_id) %>%
  summarise(p_cover = mean(veg_present),
            veg_height = mean(veg_height, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = p_cover, y = veg_height, color = factor(group_id)))+
  geom_point()

veg_obs %>%
  # get just 30m transects, 4 per plot
  filter(y_coord==15|x_coord==15) %>%
  group_by(group_id, plot_id) %>%
  summarise(p_cover = mean(veg_present),
            veg_height = mean(veg_height, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(source = 'simulated',
         plot_id = as.character(plot_id)) %>%
  bind_rows(
    readRDS(here::here('02-data',
                       '04-geolocated',
                       'seki_metadata.rds')) %>%
      mutate(group_id = as.integer(factor(mort, levels = c('low', 'mid', 'high')))) %>%
      select(group_id, plot_id, transect_id) %>%
      expand(nesting(group_id, plot_id, transect_id),
             location_m = seq(0.5, 14.5, 1)) %>%
      left_join(readRDS(here::here('02-data',
                                   '03-clean',
                                   'seki_veg.rds')) %>%
                  select(plot_id, transect_id, start_m, end_m)) %>%
      group_by(group_id, plot_id, transect_id, location_m) %>%
      summarise(veg_present = any(location_m>=start_m & location_m<=end_m)) %>%
      ungroup() %>%
      mutate(veg_present = ifelse(is.na(veg_present), FALSE, veg_present)) %>%
      group_by(group_id, plot_id) %>%
      summarise(p_cover = mean(veg_present)) %>%
      ungroup() %>%
      mutate(source = 'real')
  ) %>%
  ggplot(aes(x = p_cover, fill = source))+
  geom_histogram()+
  facet_grid(source+group_id~.)


veg_obs %>%
  # get just 30m transects, 4 per plot
  filter(y_coord==15|x_coord==15) %>%
  group_by(group_id, plot_id) %>%
  summarise(p_cover = mean(veg_present),
            veg_height = mean(veg_height, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(source = 'simulated',
         plot_id = as.character(plot_id)) %>%
  bind_rows(
    readRDS(here::here('02-data',
                       '04-geolocated',
                       'seki_metadata.rds')) %>%
      mutate(group_id = as.integer(factor(mort, levels = c('low', 'mid', 'high')))) %>%
      select(group_id, plot_id, transect_id) %>%
      expand(nesting(group_id, plot_id, transect_id),
             location_m = seq(0.5, 14.5, 1)) %>%
      left_join(readRDS(here::here('02-data',
                                   '03-clean',
                                   'seki_veg.rds')) %>%
                  select(plot_id, transect_id, start_m, end_m)) %>%
      group_by(group_id, plot_id, transect_id, location_m) %>%
      summarise(veg_present = any(location_m>=start_m & location_m<=end_m),
                veg_height = max(veg_height, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(veg_present = ifelse(is.na(veg_present), FALSE, veg_present),
             veg_height = ifelse(veg_present, veg_height, NA)) %>%
      group_by(group_id, plot_id) %>%
      summarise(p_cover = mean(veg_present),
                veg_height = mean(veg_height, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(source = 'real')
  ) %>%
  ggplot(aes(x = veg_height, fill = source))+
  geom_histogram()+
  facet_grid(source+group_id~.)



#### reality comparison: trees #################################################


#### reality comparison: saplings ##############################################

