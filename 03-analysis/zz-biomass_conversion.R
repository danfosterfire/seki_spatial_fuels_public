# Rfuels expects a pretty rigid data format which the spatially explicit 
# protocol doesn't work with. So we can't just call Rfuels::estimate_fuel_loads() 
# to convert the direct observations to load estimates. However, we will 
# follow the same conversion scheme (using species-specific coefficients, 
# weighted by basal area of the local overstory) as used by Rfuels 
# and described at https://github.com/danfosterfire/Rfuels.

# I will use a study-wide overstory composition, because we are interested 
# in analysing within-study spatial variation, and I don't want to be 
# introducing variation via the species coefficients. 

#### setup #####################################################################

# store the constants tables reported from the van wagtendonk papers and 
# dump them into some RDS files; this is a modified version of the RFuels 
# script (which dumps data into RDS files instead of constructing it 
# implicitly as part of the R package)
source(here::here('00-R', 'build_constants_tables.R'))
rm(list = ls())

# load packages
library(here)
library(tidyverse)

# build the constants tables and dump them into some .RDS files


# load data
seki.cwd = 
  readRDS(here::here('02-data',
                     '04-geolocated',
                     'seki_cwd.rds'))

seki.fwd = 
  readRDS(here::here('02-data',
                     '04-geolocated',
                     'seki_fwdtallies.rds'))


seki.litterduff = 
  readRDS(here::here('02-data',
                     '04-geolocated',
                     'seki_litterduff.rds'))

seki.veg = 
  readRDS(here::here('02-data',
                     '04-geolocated',
                     'seki_veg.rds'))

seki.trees = 
  readRDS(here::here('02-data',
                     '04-geolocated',
                     'seki_trees.rds'))

seki.metadata = readRDS(here::here('02-data',
                                    '04-geolocated',
                                    'seki_metadata.rds'))


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


#### get overstory composition #################################################

# need species composition for each mortality class expressed as propotion of total BA
head(seki.trees)

# first aggregate to plot level by summing individual trees
# note that we lump live and dead trees together for the purpose of estimating 
# fuels, because both contribute (or have contributed) to the existing fuel bed
overstory_composition = 
  seki.trees %>%
  
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

#### calculate composition-weighted coefficients ###############################

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


# could compare to the observed QMD for each size class here....

#### estimate litter + duff loads ##############################################

#' Litter and duff are measured as depths at specific points along a sampling
#' transect. Van WAgtendonk et al. (1998) developed regressions for litter,
#' duff, and combined-litter-and-duff loading (kg / m^2) as a function of
#' depth (cm) for 19 different Sierra Nevada conifer species. See vignette
#' for details

# we're going to use the combined loading

seki.litterduff = 
  seki.litterduff %>%
  mutate(litterduff_cm = round(litter_cm + duff_cm,0)) %>%
  select(mort, plot_id, transect_id, az, location_m, 
         x_coord, y_coord, litter_cm, duff_cm, fuel_cm, litterduff_cm, comments) %>%
  
  # get the stimated litter+duff load for each observation
  left_join(litterduff_coeffs) %>%
  mutate(litterduff_kgm2 = litterduff_cm * weighted_coeff,
         litterduff_mgha = litterduff_kgm2 * 10) %>%
  
  # standardize columns
  select(mort, plot_id, transect_id, az, location_m,
         x_coord, y_coord, litter_cm, duff_cm, fuel_cm, litterduff_cm, 
         litterduff_mgha, comments)

seki.litterduff %>% print(width = Inf)

ggplot(seki.litterduff,
       aes(x = mort, y = litterduff_mgha))+
  geom_boxplot()

ggplot(seki.litterduff,
       aes(x = mort, y = litter_cm))+
  geom_boxplot()

ggplot(seki.litterduff,
       aes(x = mort, y = duff_cm))+
  geom_boxplot()


# compare with FFS data; only works on UC laptop
ffs.fuels = 
  readRDS(here::here('02-data', '00-source', 'ffs_fuels.rds'))

ffs.fuels %>%
  filter(timestep=='pre_treatment'|
           treatment=='control') %>%
  mutate(litterduff_mgha = litter_mgha+duff_mgha) %>%
  ggplot(aes(x =treatment, y = litterduff_mgha))+
  geom_boxplot()

ffs.fuels %>%
  filter(timestep=='pre_treatment'|
           treatment=='control') %>%
  select(treatment, timestep, duff_a_cm, duff_b_cm, litter_a_cm, litter_b_cm, 
         litter_mgha, duff_mgha) %>%
  mutate(litterduff_mgha = litter_mgha+duff_mgha) %>%
  pivot_longer(cols = c(duff_a_cm, duff_b_cm, litter_a_cm, litter_b_cm),
               names_to = c('class', 'position', 'units'),
               names_sep = '_',
               values_to = 'value') %>%
  ggplot(aes(x = class, y = value))+
  geom_boxplot()

ffs.fuels %>%
  filter(timestep=='pre_treatment'|treatment=='control') %>%
  mutate(mort = 'ffs_pretreatment') %>%
  select(mort, plot_id, timestep, azimuth, duff_a_cm, duff_b_cm, litter_a_cm, litter_b_cm) %>%
  pivot_longer(cols = c(duff_a_cm, duff_b_cm, litter_a_cm, litter_b_cm),
               names_to = c('class', 'position', 'units'),
               names_sep = '_',
               values_to = c('depth_cm')) %>%
  pivot_wider(id_cols = c('mort','timestep', 'plot_id', 'azimuth', 'position'),
              names_from = 'class', values_from = 'depth_cm') %>%
  select(mort, plot_id, duff_cm = duff, litter_cm = litter) %>%
  bind_rows(seki.litterduff %>%
              select(mort, plot_id, duff_cm, litter_cm)) %>%
  pivot_longer(cols = c(duff_cm, litter_cm),
               names_to = 'class', values_to = 'depth_cm') %>%
  ggplot(aes(x = mort, y = depth_cm))+
  geom_boxplot(position = position_nudge(x = -0.25),
               width = 0.2)+
  geom_jitter(height = 0, width = 0.1, alpha = 0.1)+
  geom_violin(width = 0.2, position = position_nudge(x = 0.25))+
  facet_grid(.~class)+
  theme_bw()

ffs.fuels %>%
  filter(timestep=='pre_treatment'|treatment=='control') %>%
  mutate(mort = 'ffs',
         source = 'ffs') %>%
  select(source, mort, timestep, plot_id, azimuth, duff_a_cm, duff_b_cm, litter_a_cm, litter_b_cm) %>%
  pivot_longer(cols = c(duff_a_cm, duff_b_cm, litter_a_cm, litter_b_cm),
               names_to = c('class', 'position', 'units'),
               names_sep = '_',
               values_to = c('depth_cm')) %>%
  pivot_wider(id_cols = c('source', 'mort', 'timestep', 'plot_id', 'azimuth', 'position'),
              names_from = 'class', values_from = 'depth_cm') %>%
  select(source, mort, plot_id, duff_cm = duff, litter_cm = litter) %>%
  bind_rows(seki.litterduff %>%
              mutate(source = 'seki') %>%
              select(source, mort, plot_id, duff_cm, litter_cm)) %>%
  pivot_longer(cols = c(duff_cm, litter_cm),
               names_to = 'class', values_to = 'depth_cm') %>%
  ggplot(aes(x = depth_cm))+
  geom_histogram()+
  facet_grid(source~class)




# oooooofff... i don't like this.... I see two modes in the SEKI data 
# which makes "they were sometimes measuring in inches and sometimes cm" 
# seem more likely...
# OK, but if they were sometimes incorrectly measuring in inches, instead of 
# cm, then the observed depth values for SEKI should be LOWER than the ffs 
# data, not higher. That's not it. 


ffs.fuels %>%
  filter(timestep=='pre_treatment'|treatment=='control') %>%
  mutate(mort = 'ffs',
         litterduff_mgha = litter_mgha+duff_mgha) %>%
  select(mort, plot_id, litterduff_mgha) %>%
  bind_rows(seki.litterduff %>%
              select(mort, plot_id, litterduff_mgha)) %>%
  ggplot(aes(x = mort, y = litterduff_mgha))+
  geom_boxplot(width = 0.25, position = position_nudge(x = -.25))+
  geom_violin(width = 0.1, position = position_nudge(x = 0.25))+
  geom_jitter(height = 0, width = 0.1, alpha = 0.5)+
  geom_point(data = 
               ffs.fuels %>%
               filter(timestep=='pre_treatment'|treatment=='control') %>%
               mutate(mort = 'ffs',
                      litterduff_mgha = litter_mgha+duff_mgha) %>%
               select(mort, plot_id, litterduff_mgha) %>%
               bind_rows(seki.litterduff %>%
                           select(mort, plot_id, litterduff_mgha)) %>%
               group_by(mort) %>%
               summarise(mgha.mean = mean(litterduff_mgha, na.rm = TRUE),
                         mgha.se = sd(litterduff_mgha, na.rm = TRUE)/sqrt(n())) %>%
               ungroup(),
             aes(x = mort, y = mgha.mean), size = 5, color = 'red')+
  geom_errorbar(data = 
                  ffs.fuels %>%
                  filter(timestep=='pre_treatment'|treatment=='control') %>%
                  mutate(mort = 'ffs',
                         litterduff_mgha = litter_mgha+duff_mgha) %>%
                  select(mort, plot_id, litterduff_mgha) %>%
                  bind_rows(seki.litterduff %>%
                              select(mort, plot_id, litterduff_mgha)) %>%
                  group_by(mort) %>%
                  
                  summarise(mgha.mean = mean(litterduff_mgha, na.rm = TRUE),
                            mgha.se = sd(litterduff_mgha, na.rm = TRUE)/sqrt(n())) %>%
                  ungroup(),
                aes(x = mort, ymin = mgha.mean-(2*mgha.se), ymax = mgha.mean+(2*mgha.se),
                    y = mgha.mean), 
                width = 0, color = 'red', lwd = 1)+
  theme_bw()

ffs.fuels
summary(ffs.fuels %>% filter(timestep=='pre_treatment'|treatment=='control') %>%
          mutate(litterduff_mgha = litter_mgha+duff_mgha) %>%
          select(litter_mgha, duff_mgha, litterduff_mgha))

summary(ffs.fuels %>%
          filter(timestep=='post_1'&treatment=='control') %>%
          select(litter_mgha, duff_mgha))

summary(seki.litterduff %>% select(litter_cm, duff_cm, litterduff_cm, litterduff_mgha))
summary(ffs.fuels %>% 
          filter(timestep=='post_1'|treatment=='control') %>%
          mutate(litterduff_a_cm = litter_a_cm+duff_a_cm,
                 litterduff_b_cm = litter_b_cm+duff_b_cm,
                 litterduff_mgha = litter_mgha+duff_mgha) %>%
          select(litter_a_cm, litter_b_cm, duff_a_cm, duff_b_cm, 
                 litterduff_a_cm, litterduff_b_cm, litterduff_mgha))
# i think the bimodal thing is an artifact, I'm not seeing it on jitter or 
# violin plots


seki.litterduff %>%
  ggplot(aes(x = litterduff_mgha))+
  #geom_histogram()+
  geom_density(color = 'red', lwd = 1)+
  facet_grid(mort~.)

# i think this is real; the recorded litter depths in SEKI are about 2x 
# the recorded depths in the pre-treatment FFS data, but the duff depths 
# are comparable in SEKI and BFRS (which kind of rules out a "inches instead of cm" 
# sort of issue with the SEKI data; if the crew were using the wrong units it 
# shoult affect both the litter and duff depths); it also rules out a 
# "danny messed up the coding for biomass here" issue because the litter 
# depths are so different (not just the litterduff_mgha values)

# the bimodal depth distribution is weird; though maybe an artifact of the 
# historgams? doesn't show up in a density graph

#### FWD loading ###############################################################

seki.fwd = 
  seki.fwd %>%
  
  # reformat tallies to longwise
  pivot_longer(cols = c(a1h, a10h, a100h, b1h, b10h, b100h),
               names_to = 'subsample_id',
               values_to = 'count') %>%
  mutate(timelag_class = gsub(x = subsample_id,
                              pattern = 'a|b',
                              replacement = ''),
         subsample = gsub(x = subsample_id,
                          pattern = '1h|10h|100h',
                          replacement = '')) %>%
  
  # add columns for slope, QMD, SEC, SG, and transect length, and the conversion 
  # constant k
  left_join(., 
            y = QMDcm) %>%
  left_join(.,
            y = SEC) %>%
  left_join(.,
            y = SG) %>%
  left_join(.,
            y = 
              seki.metadata %>%
              mutate(
                slp_pcent = 
                  tan(slp_degrees*(pi/180))*100,
                slp_c = sqrt(1+(slp_pcent/100)^2)) %>%
              select(transect_id, slp_c)) %>%
  mutate(transect_length_m = 1,
         k = 1.234) %>%
  
  # use brown's equations to estimate fuel load
  mutate(mgha = 
           (k*weighted_qmd*weighted_sec*slp_c*weighted_sg*count)/
           transect_length_m) %>%
  
  # pivot back wider
  pivot_wider(id_cols = c('mort', 'plot_id', 'transect_id', 'az', 'location_m',
                          'x_coord', 'y_coord'),
              names_from = 'subsample_id',
              values_from = c('count', 'mgha'))

ffs.fuels %>%
  filter(timestep=='pre_treatment') %>%
  select(x1h_mgha, x10h_mgha, x100h_mgha) %>% 
  mutate(source = 'ffs', mort = 'ffs',
         fwd_mgha = x1h_mgha+x10h_mgha+x100h_mgha) %>%
  bind_rows(seki.fwd %>%
              mutate(source = 'seki') %>%
              select(source, mort, transect_id, location_m,
                     mgha_a1h, mgha_a10h, mgha_a100h,
                     mgha_b1h, mgha_b10h, mgha_b100h) %>%
              pivot_longer(cols = c('mgha_a1h', 'mgha_b1h',
                                    'mgha_a10h', 'mgha_b10h',
                                    'mgha_a100h', 'mgha_b100h'),
                           names_to = c('units', 'group'),
                           names_sep = '_',
                           values_to = 'mgha') %>%
              mutate(subsample = 
                       gsub(x = group, pattern = '1h|10h|100h', replacement = ''),
                     timelag = 
                       paste0('x',
                            gsub(x = group, pattern = 'a|b', replacement = ''))) %>%
              pivot_wider(id_cols = c('source', 'mort','transect_id', 'location_m', 'subsample'),
                          names_from = 'timelag', values_from = 'mgha') %>%
              select(source, mort, x1h_mgha = x1h, x10h_mgha = x10h, x100h_mgha = x100h) %>%
              mutate(fwd_mgha = x1h_mgha+x10h_mgha+x100h_mgha)) %>%
  pivot_longer(cols = c('x1h_mgha', 'x10h_mgha', 'x100h_mgha', 'fwd_mgha'),
               names_to = 'variable', values_to = 'value') %>%
  ggplot(aes(x = source, y = value))+
  geom_boxplot()+
  facet_wrap(~variable, scales = 'free_y')


# fwd loadings look very similar between seki and FFS pretreatment... 
# so duff and fwd are similar, but seki has 2x the litter? does that
# make sense?


#### CWD loading ###############################################################

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
                              location_m = 
                                c(5, 15, 25)) %>%
               tidyr::expand(nesting(mort, plot_id, transect_id, az, x_center,
                                     y_center, slp_c, location_m),
                      timelag_class = c('1000s', '1000r'))

seki.cwd = 
  seki.cwd %>%
  # tidy up 
  mutate(timelag_class = ifelse(as.numeric(decay_class) <= 2, '1000s', '1000r')) %>%
  select(mort, plot_id, transect_id, az, location_m, 
         timelag_class, diam_cm, comments, x_center, y_center) %>%
  
  # get the sum of squared diameters for each subtransect and timelag class
  # round the location to the nearest 5m midpoint (2.5, 7.5,...27.5m)
  mutate(location_m = floor(location_m/10)*10+5) %>%
  group_by(mort, plot_id, transect_id, az, x_center, y_center, location_m, 
           timelag_class) %>%
  summarize(ssd_cm2 = sum(diam_cm^2)) %>%
  ungroup() %>%
  right_join(subtransects) %>%
  mutate(ssd_cm2 = ifelse(is.na(ssd_cm2), 0, ssd_cm2))# %>%
  
  # and join in the coefficients 
  left_join(.,
            y = SEC) %>%
  left_join(.,
            y = SG) %>%
  mutate(transect_length_m = 10, # using 5m subtransects
         k = 1.234) %>%
  
  # and use brown's equation (with the new k const) to estimate fuel load
  mutate(mgha = 
           (k*ssd_cm2*weighted_sec*slp_c*weighted_sg)/
           transect_length_m) %>%
  
  # drop some intermediate columns
  select(mort, plot_id, transect_id, az, x_center, y_center, 
         timelag_class, location_m, ssd_cm2, mgha)
  
  # going to lump sound and rotten together to reduce the number of 0s and the 
  # dimensionality of the response, mark will @ me if he wants it changed
  mutate(timelag_class = 'cwd') %>%
  group_by(mort, plot_id, transect_id, az, x_center, y_center, 
           location_m, timelag_class) %>%
  summarise(ssd_cm2 = sum(ssd_cm2),
            mgha = sum(mgha)) %>%
  ungroup()

head(seki.cwd)
summary(seki.cwd)
ggplot(data = seki.cwd,
       aes(x = mgha))+
  geom_histogram()

ggplot(data = seki.cwd,
       aes(x = mort, y = mgha))+
  geom_boxplot()




#### trees #####################################################################

seki.trees



table9 = read.csv(here::here('02-data',
                             '00-source',
                             'gill',
                             'table9.csv'))

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
  
  # map
unique(seki.trees$spp)[order(unique(seki.trees$spp))]

ref_species %>% print(width = Inf)
test = 
  ref_species %>%
  filter(str_detect(SPECIES_SYMBOL,
                    pattern = 
                      paste(unique(seki.trees$spp),
                            collapse = '|'))) %>%
  arrange(GENUS, SPECIES, SPECIES_SYMBOL)

test = 
  ref_species %>%
  select(GENUS, SPECIES, SPECIES_SYMBOL) %>%
  arrange(GENUS, SPECIES, SPECIES_SYMBOL)


seki.trees = 
  seki.trees %>%
  
  # estimate the crown radius of each tree using the results of gill et al. 2000
  # first join in the regression coefficients
  left_join(table9, by = c('spp' = 'Species')) %>%
  
  # fill in unmatched values with the 'other conifer' value (not a great fit for 
  # the hardwoods in the dataset, but OK for a first pass) 
  mutate(b0 = ifelse(is.na(b0), table9[table9$Species=='OTHER','b0'], b0),
         b1 = ifelse(is.na(b1), table9[table9$Species=='OTHER','b1'], b1)) %>%
  
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
      agbiomass_lb-stembiomass_lb-barkbiomass_lb-foliagebiomass_lb
    
    # this was spitting out negative bole biomass for some saplings
    #saplingbiomass_lb = 
    #  (agbiomass_lb-foliagebiomass_lb)*JENKINS_SAPLING_ADJUSTMENT
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
    
    # the sapling adjustment lumps the branches in with the rest of the bole,
    # we're going to do the same for saplings
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
    mort, plot_id, tree_id, x_coord, y_coord, ts, location_m, distance_m,
    status, spp, dbh_cm, height_m, htcb_m, decay_class, comments,
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


ggplot(data = 
         seki.trees,
       aes(x = dbh_cm, y = agbiomass_kg))+
  geom_point()

ggplot(data = 
         seki.trees %>%
         group_by(mort, plot_id, status, spp) %>%
         summarise(tph = sum(tph),
                   ba_m2ha = sum(ba_m2ha)) %>%
         ungroup() %>%
         complete(nesting(mort, plot_id), status, spp) %>%
         mutate(tph = ifelse(is.na(tph), 0, tph),
                ba_m2ha = ifelse(is.na(ba_m2ha), 0, ba_m2ha)),
       aes(x = plot_id, y = ba_m2ha, fill = status))+
  geom_col(position = position_stack())+
  facet_grid(.~mort, scales = 'free_x')+
  theme_bw()

ggplot(data = 
         seki.trees %>%
         group_by(mort, plot_id, status, spp) %>%
         summarise(tph = sum(tph),
                   agbiomass_mgha = sum(agbiomass_mgha)) %>%
         ungroup() %>%
         complete(nesting(mort, plot_id), status, spp) %>%
         mutate(tph = ifelse(is.na(tph), 0, tph),
                agbiomass_mgha = ifelse(is.na(agbiomass_mgha), 0, agbiomass_mgha)),
       aes(x = plot_id, y = agbiomass_mgha, fill = status))+
  geom_col(position = position_stack())+
  facet_grid(.~mort, scales = 'free_x')+
  theme_bw()


ggplot(data = 
         seki.trees %>%
         group_by(mort, plot_id, status, spp) %>%
         summarise(tph = sum(tph),
                   ba_m2ha = sum(ba_m2ha)) %>%
         ungroup() %>%
         complete(nesting(mort, plot_id), status, spp) %>%
         mutate(tph = ifelse(is.na(tph), 0, tph),
                ba_m2ha = ifelse(is.na(ba_m2ha), 0, ba_m2ha)),
       aes(x = plot_id, y = tph, fill = status))+
  geom_col(position = position_stack())+
  facet_grid(.~mort, scales = 'free_x')+
  theme_bw()

ggplot(data = 
         seki.trees %>%
         group_by(mort, plot_id, status) %>%
         summarise(ba_m2ha = sum(ba_m2ha)) %>%
         ungroup() %>%
         mutate(mort = factor(mort, levels = c('low', 'mid', 'high'))) %>%
         complete(nesting(mort, plot_id), status) %>%
         mutate(ba_m2ha = ifelse(is.na(ba_m2ha), 0, ba_m2ha)) %>%
         group_by(mort, status) %>%
         summarise(ba_m2ha.mean = mean(ba_m2ha),
                   ba_m2ha.se = sd(ba_m2ha)/sqrt(n())) %>%
         ungroup(),
       aes(x = mort))+
  geom_errorbar(aes(ymin = ba_m2ha.mean-(2*ba_m2ha.se),
                    ymax = ba_m2ha.mean+(2*ba_m2ha.se)),
                width = 0.1, lwd = 1)+
  geom_point(aes(y = ba_m2ha.mean, color = status), size = 5)+
  facet_grid(status~.)+
  theme_bw()

ggplot(data = 
         seki.trees %>%
         group_by(mort, plot_id, status) %>%
         summarise(ba_m2ha = sum(ba_m2ha)) %>%
         ungroup() %>%
         complete(nesting(mort, plot_id), status) %>%
         mutate(ba_m2ha = ifelse(is.na(ba_m2ha), 0, ba_m2ha)) %>%
         left_join(seki.trees %>%
                     group_by(mort, plot_id) %>%
                     summarise(totalba_m2ha = sum(ba_m2ha)) %>%
                     ungroup()) %>%
         mutate(mort = factor(mort, levels = c('low', 'mid', 'high')),
                pBA = ba_m2ha / totalba_m2ha) %>%
         filter(status=='L') %>%
         group_by(mort) %>%
         summarise(pLive.mean = mean(pBA),
                   pLive.se = sd(pBA)/sqrt(n())) %>%
         ungroup(),
       aes(x = mort))+
  geom_errorbar(aes(ymin = pLive.mean-(2*pLive.se),
                    ymax = pLive.mean+(2*pLive.se)),
                width = 0.1, lwd = 1)+
  geom_point(aes(y = pLive.mean, color = mort), size = 5)+
  labs(x = 'Mortality Class',
       y = 'Proportion Live')+
  theme_bw()+
  scale_color_viridis_d(option = 'C',
                        begin = 0.05, end = 0.85)

# BA, TPH, and biomass all look high but very plausible at this site
# low/med mortality classes are pretty similar, and both distinct from high
# which makes sense, given that the classes were based on quantiles, so
# "high" mortality is a very wide range

#### shrubs ####################################################################

# tbd


#### scratch ###################################################################
head(seki.litterduff)

library(lme4)

# these fits are busted, they don't take into account the within-transect 
# spatial correlation. I'm just testing to see that even under these 
# pseudoreplicated models whether fuel loads vary significantly across mortality
# classes
litterduff_fit = 
  nlme::lme(data = seki.litterduff %>% filter(!is.na(litterduff_mgha)), 
              litterduff_mgha ~ mort,
            random =~ 1|plot_id/transect_id)

summary(litterduff_fit)

head(seki.fwd)
fwd_fit = 
  nlme::lme(data = 
              seki.fwd %>% 
              mutate(fwd_mgha = 
                       (mgha_a1h+mgha_b1h)*0.5 + 
                       (mgha_a10h+mgha_b10h)*0.5 +
                       (mgha_a100h+mgha_b100h)*0.5) %>%
              filter(!is.na(fwd_mgha)), 
              fwd_mgha ~ mort,
            random =~ 1|plot_id/transect_id)

summary(fwd_fit)


#### write results #############################################################

saveRDS(seki.cwd,
        here::here('02-data',
                   '05-for_analysis',
                   'biomass',
                   'seki_cwd.rds'))
write.csv(seki.cwd,
          here::here('02-data',
                     '05-for_analysis',
                     'biomass',
                     'seki_cwd.csv'),
          row.names = FALSE)


saveRDS(seki.fwd,
        here::here('02-data',
                   '05-for_analysis',
                   'biomass',
                   'seki_fwd.rds'))
write.csv(seki.fwd,
          here::here('02-data',
                     '05-for_analysis',
                     'biomass',
                     'seki_fwd.csv'),
          row.names = FALSE)


saveRDS(seki.litterduff,
        here::here('02-data',
                   '05-for_analysis',
                   'biomass',
                   'seki_litterduff.rds'))
write.csv(seki.litterduff,
          here::here('02-data',
                     '05-for_analysis',
                     'biomass',
                     'seki_litterduff.csv'),
          row.names = FALSE)


saveRDS(seki.trees,
        here::here('02-data',
                   '05-for_analysis',
                   'biomass',
                   'seki_trees.rds'))
write.csv(seki.trees,
          here::here('02-data',
                     '05-for_analysis',
                     'biomass',
                     'seki_trees.csv'),
          row.names = FALSE)







