

#### setup #####################################################################

library(here)
library(sf)

library(tidyverse)

metadata = 
  readRDS(here::here('02-data', '03-clean', 'seki_metadata.rds'))

head(metadata)


plot_locations.tibble = 
  readRDS(here::here('02-data', '02-intermediate', 'plot_locations_true.rds'))

head(plot_locations.tibble)
plot_locations.tibble[,c('x_true', 'y_true')] %>%
  as.matrix() %>%
  dist() %>%
  summary()

#### geolocate metadata ########################################################




metadata = 
  metadata %>%
  left_join(plot_locations.tibble %>%
              select(mort, plot_id, x_center = x_true, y_center = y_true),
            by = c('plot_id'))



#### litter and duff separately ################################################

litterduff = 
  readRDS(here::here('02-data', '03-clean', 'seki_litterduff.rds')) %>%
  mutate(litterduff_cm0 = round(litter_cm+duff_cm, 0)) %>%
  left_join(metadata %>%
              select(plot_id, transect_id, mort)) %>%
  filter(!is.na(litterduff_cm0)) %>%
  mutate(x_rel = (location_m*sin(az*(pi/180))),
         y_rel = (location_m*cos(az*(pi/180))))
  
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
litterduff_training = 
  litterduff %>%
  group_by(mort) %>%
  slice_sample(prop = 0.9) %>%
  ungroup()

litterduff_validation = 
  litterduff %>%
  filter(!is.element(obs_id,
                     litterduff_training %>%
                       pull(obs_id)))



litter_training_data = 
  list(N = nrow(litterduff_training),
       P = length(unique(litterduff$plot_id)),
       L = nrow(coords_litterduff),
       G = length(unique(litterduff$mort)),
       coords = coords_litterduff[,c('x_rel', 'y_rel')] %>% as.matrix(),
       plot_id = litterduff_training$plot_id.i,
       location_id = litterduff_training$location_id,
       group_id = 
         litterduff %>%
         group_by(plot_id.i, group_id) %>%
         summarise() %>%
         ungroup() %>%
         arrange(plot_id.i) %>%
         pull(group_id) %>%
         as.numeric(),
       Y = round(litterduff_training$litter_cm,0))


litter_validation_data = 
  list(N = nrow(litterduff_validation),
       P = length(unique(litterduff$plot_id)),
       L = nrow(coords_litterduff),
       G = length(unique(litterduff$mort)),
       coords = coords_litterduff[,c('x_rel', 'y_rel')] %>% as.matrix(),
       plot_id = litterduff_validation$plot_id.i,
       location_id = litterduff_validation$location_id,
       group_id = 
         litterduff %>%
         group_by(plot_id.i, group_id) %>%
         summarise() %>%
         ungroup() %>%
         arrange(plot_id.i) %>%
         pull(group_id) %>%
         as.numeric(),
       Y = round(litterduff_validation$litter_cm, 0))


duff_training_data = 
  list(N = nrow(litterduff_training),
       P = length(unique(litterduff$plot_id)),
       L = nrow(coords_litterduff),
       G = length(unique(litterduff$mort)),
       coords = coords_litterduff[,c('x_rel', 'y_rel')] %>% as.matrix(),
       plot_id = litterduff_training$plot_id.i,
       location_id = litterduff_training$location_id,
       group_id = 
         litterduff %>%
         group_by(plot_id.i, group_id) %>%
         summarise() %>%
         ungroup() %>%
         arrange(plot_id.i) %>%
         pull(group_id) %>%
         as.numeric(),
       Y = round(litterduff_training$duff_cm,0))


duff_validation_data = 
  list(N = nrow(litterduff_validation),
       P = length(unique(litterduff$plot_id)),
       L = nrow(coords_litterduff),
       G = length(unique(litterduff$mort)),
       coords = coords_litterduff[,c('x_rel', 'y_rel')] %>% as.matrix(),
       plot_id = litterduff_validation$plot_id.i,
       location_id = litterduff_validation$location_id,
       group_id = 
         litterduff %>%
         group_by(plot_id.i, group_id) %>%
         summarise() %>%
         ungroup() %>%
         arrange(plot_id.i) %>%
         pull(group_id) %>%
         as.numeric(),
       Y = round(litterduff_validation$duff_cm, 0))




saveRDS(litter_training_data,
        here::here('02-data', '05-for_analysis', 'litter_training_data.rds'))


saveRDS(litter_validation_data,
        here::here('02-data', '05-for_analysis', 'litter_validation_data.rds'))


saveRDS(duff_training_data,
        here::here('02-data', '05-for_analysis', 'duff_training_data.rds'))


saveRDS(duff_validation_data,
        here::here('02-data', '05-for_analysis', 'duff_validation_data.rds'))



#### fwd tallies ###############################################################


fwd = 
  readRDS(here::here('02-data', '03-clean', 'seki_fwdtallies.rds'))   %>%

  # get relative spatial coordinates
  mutate(x_rel = round(location_m*sin(az*pi/180),2),
       y_rel = round(location_m*cos(az*pi/180),2)) %>%
  
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
  left_join(metadata %>%
              select(plot_id, transect_id, mort))

head(fwd)

fwd_wide = 
  fwd %>%
  group_by(mort, plot_id, transect_id, az, location_m, timelag_class) %>%
  summarise(count = mean(count, na.rm = TRUE)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c('plot_id', 'transect_id', 'az', 'location_m'),
              names_from = c('timelag_class'),
              values_from = 'count',
              names_prefix = 'x')

ggplot(fwd_wide,
       aes(x = x1h, y = x10h))+
  geom_point()+
  geom_smooth(method = 'lm')

ggplot(fwd_wide,
       aes(x = x10h, y = x100h))+
  geom_jitter()+
  geom_smooth(method = 'lm')


coords_fwd = 
  fwd %>%
  group_by(x_rel, y_rel) %>%
  summarise() %>%
  ungroup() %>%
  rowid_to_column('location_id')

fwd = 
  fwd %>%
  left_join(coords_fwd)

# from notes after fitting separate models for 1h, 10h, and 100h:

# takeaway: 
#   - realistic simulations clearly require correlation amongst the FWD classes, 
#   which is clear in the raw data. 
#   - 100h fwd, and to a lesser extent 10h fwd, are rare enough that the posterior 
#   distributions for the GP parameters are mostly informed by the prior, rather than 
#   the data
#   - there is mixed evidence about the GP parameters varying by size class. 
#   Some parameters vary more by mortality class than size class (alpha), and 
#  the data are not informative about some parameters for 10h and 100h fuels (rho), 
#  and the data are suggestive about differences among size classes for tau.
#  - given the 3 above points, I think the most realistic simulations will be 
# for a generic FWD model, lumping the tallies to include 1-100h and then assigning 
# particle diameters from the FWD diameters model.

# so making a dataset that gets total tally of 1-100h fuels at each location
fwd1to100h = 
  fwd %>%
  group_by(mort, plot_id, transect_id, location_id, x_rel, y_rel, subsample) %>%
  summarise(count = sum(count, na.rm = FALSE)) %>%
  ungroup() %>%
  filter(!is.na(count)) %>%
  mutate(plot_id.i = as.integer(factor(plot_id)),
         group_id = as.integer(factor(mort, levels = c('low', 'mid', 'high')))) 

head(fwd1to100h)

ggplot(data = fwd1to100h,
       aes(x = count))+
  geom_density()+
  facet_grid(mort~.)

fwd1to100h = 
  fwd1to100h %>%
  rowid_to_column('obs_id')

set.seed(110819)
fwd1to100h.training = 
  fwd1to100h %>%
  group_by(mort) %>%
  slice_sample(prop = 0.9) %>%
  ungroup()

fwd1to100h.validation = 
  fwd1to100h %>%
  filter(!is.element(obs_id,
                     fwd1to100h.training$obs_id))

fwdtallies_training_data = 
  list(N = nrow(fwd1to100h.training),
       P = length(unique(fwd1to100h$plot_id)),
       L = nrow(coords_fwd),
       G = length(unique(fwd1to100h$mort)),
       coords = coords_fwd[,c('x_rel', 'y_rel')] %>% as.matrix(),
       plot_id = fwd1to100h.training$plot_id.i,
       location_id = fwd1to100h.training$location_id,
       group_id = 
         fwd1to100h %>%
         group_by(plot_id.i, group_id) %>%
         summarise() %>%
         ungroup() %>%
         arrange(plot_id.i) %>%
         pull(group_id) %>%
         as.numeric(),
       Y = fwd1to100h.training$count)


fwdtallies_validation_data = 
  list(N = nrow(fwd1to100h.validation),
       P = length(unique(fwd1to100h$plot_id)),
       L = nrow(coords_fwd),
       G = length(unique(fwd1to100h$mort)),
       coords = coords_fwd[,c('x_rel', 'y_rel')] %>% as.matrix(),
       plot_id = fwd1to100h.validation$plot_id.i,
       location_id = fwd1to100h.validation$location_id,
       group_id = 
         fwd1to100h %>%
         group_by(plot_id.i, group_id) %>%
         summarise() %>%
         ungroup() %>%
         arrange(plot_id.i) %>%
         pull(group_id) %>%
         as.numeric(),
       Y = fwd1to100h.validation$count)

saveRDS(fwdtallies_training_data,
        here::here('02-data',
                   '05-for_analysis',
                   'fwdtallies_training_data.rds'))


saveRDS(fwdtallies_validation_data,
        here::here('02-data',
                   '05-for_analysis',
                   'fwdtallies_validation_data.rds'))




#### fwd diameters #############################################################



fwd_diams = 
  
  # start with the FWD diams
  readRDS(here::here('02-data',
                     '03-clean',
                     'seki_fwddiams.rds')) %>%
  
  left_join(metadata %>%
              select(plot_id, transect_id, az, mort, x_center, y_center)) %>%
  
  mutate(direction = ifelse(is.element(az, c(0, 180)),'NS', 'EW'),
         transect_id = paste0(plot_id,'-',az)) %>%
  
  mutate(plot_id.i = as.integer(factor(plot_id)),
         transect_id.i = as.integer(factor(transect_id)),
         group_id = as.integer(factor(mort,
                                      levels = c('low', 'mid', 'high'))))

ggplot(data = fwd_diams,
       aes(x = diam_mm, fill = mort))+
  geom_histogram(position = position_dodge())+
  theme_minimal()+
  facet_grid(mort~.)

ggplot(data = fwd_diams,
       aes(x = diam_mm, fill = mort))+
  geom_histogram(position = position_dodge())+
  theme_minimal()

fwd_diams = 
  fwd_diams %>%
  rowid_to_column('obs_id')

set.seed(110819)
fwd_diams.training = 
  fwd_diams %>%
  group_by(mort) %>%
  slice_sample(prop = 0.9) %>%
  ungroup()

fwd_diams.validation = 
  fwd_diams %>%
  filter(!is.element(obs_id,
                     fwd_diams.training$obs_id))

fwddiams_training_data = 
  list(N = nrow(fwd_diams.training),
       P = length(unique(fwd_diams$plot_id)),
       TS = length(unique(fwd_diams$transect_id)),
       G = length(unique(fwd_diams$mort)),
       plot_id = 
         fwd_diams %>%
         group_by(transect_id.i, plot_id.i) %>%
         summarise() %>%
         ungroup() %>%
         arrange(transect_id.i) %>%
         pull(plot_id.i) %>%
         as.numeric(),
       transect_id = fwd_diams.training$transect_id.i,
       group_id = 
         fwd_diams %>%
         group_by(plot_id.i, group_id) %>%
         summarise() %>%
         ungroup() %>%
         arrange(plot_id.i) %>%
         pull(group_id) %>%
         as.numeric(),
       Y = fwd_diams.training$diam_mm/10,
       lb = 0,
       ub = 7.55)

fwddiams_validation_data = 
  list(N = nrow(fwd_diams.validation),
       P = length(unique(fwd_diams$plot_id)),
       TS = length(unique(fwd_diams$transect_id)),
       G = length(unique(fwd_diams$mort)),
       plot_id = 
         fwd_diams %>%
         group_by(transect_id.i, plot_id.i) %>%
         summarise() %>%
         ungroup() %>%
         arrange(transect_id.i) %>%
         pull(plot_id.i) %>%
         as.numeric(),
       transect_id = fwd_diams.validation$transect_id.i,
       group_id = 
         fwd_diams %>%
         group_by(plot_id.i, group_id) %>%
         summarise() %>%
         ungroup() %>%
         arrange(plot_id.i) %>%
         pull(group_id) %>%
         as.numeric(),
       Y = fwd_diams.validation$diam_mm/10,
       lb = 0,
       ub = 7.55)

saveRDS(fwddiams_training_data,
        here::here('02-data',
                   '05-for_analysis',
                   'fwddiams_training_data.rds'))

saveRDS(fwddiams_validation_data,
        here::here('02-data',
                   '05-for_analysis',
                   'fwddiams_validation_data.rds'))



#### CWD #######################################################################

cwd = 
  
  # load the data
  readRDS(here::here('02-data',
                     '03-clean',
                     'seki_cwd.rds')) %>%
  
  # each piece counts for 1 piece, and bin locations into 1m subtransects
  mutate(count = 1,
         location_m = (floor(location_m/1)*1)+0.5) %>%
  
  select(plot_id, transect_id, az, location_m, count) %>%
  
  # add in rows for 0 intersections on each meter
  bind_rows(
    plot_locations.tibble %>%
      select(plot_id) %>%
      expand(plot_id,
             az = c(0, 90, 180, 270),
             location_m = seq(from = 0.5, to = 29.5, by = 1)) %>%
      mutate(count = 0)
  ) %>%
  
  # get the total number of intersections on each meter
  group_by(plot_id, az, location_m) %>%
  summarise(count = sum(count)) %>%
  ungroup() %>%
  
  # get mortality level and transect id columns
  left_join(plot_locations.tibble %>%
              select(mort, plot_id),
            by = c('plot_id' = 'plot_id')) %>%  
  # want to treat the N and S transects as a single transect, same with the E/W
  mutate(direction = ifelse(is.element(az, c(0, 180)),'NS', 'EW'),
         transect_id = paste0(plot_id,'-',direction)) %>%
  mutate(plot_id.i = as.integer(factor(plot_id)),
         transect_id.i = as.integer(factor(transect_id)),
         group_id = as.integer(factor(mort, levels = c('low', 'mid', 'high')))) %>%
  mutate(
    x_rel = round(location_m*sin(az*pi/180),2),
    y_rel = round(location_m*cos(az*pi/180),2)
  )

head(cwd)


cwd_coords = 
  cwd  %>%
  group_by(x_rel, y_rel) %>%
  summarise() %>%
  ungroup() %>%
  rowid_to_column('location_id')

cwd = 
  cwd %>%
  left_join(cwd_coords,
            by = c('x_rel' = 'x_rel',
                   'y_rel' = 'y_rel'))

cwd = 
  cwd %>%
  rowid_to_column('obs_id')

set.seed(110819)
cwd.training = 
  cwd %>%
  group_by(mort) %>%
  slice_sample(prop = 0.9) %>%
  ungroup()

cwd.validation = 
  cwd %>%
  filter(!is.element(obs_id, cwd.training$obs_id))


cwdtallies_training_data = 
  list(N = nrow(cwd.training),
       P = length(unique(cwd$plot_id)),
       TR = length(unique(cwd$transect_id)),
       L = nrow(cwd_coords),
       G = length(unique(cwd$mort)),
       
       Y = cwd.training$count,
       location_id = cwd.training$location_id,
       transect_id = cwd.training$transect_id.i,
       plot_id = 
         cwd %>%
         group_by(transect_id.i, plot_id.i) %>%
         summarise() %>%
         ungroup() %>%
         arrange(transect_id.i) %>%
         pull(plot_id.i) %>%
         as.numeric(),
       group_id = 
         cwd %>%
         group_by(plot_id.i, group_id) %>%
         summarise() %>%
         ungroup() %>%
         arrange(plot_id.i) %>%
         pull(group_id) %>%
         as.numeric(),
       
       coords = cwd_coords[,c('x_rel', 'y_rel')] %>% as.matrix())

cwdtallies_validation_data = 
  list(N = nrow(cwd.validation),
       P = length(unique(cwd$plot_id)),
       TR = length(unique(cwd$transect_id)),
       L = nrow(cwd_coords),
       G = length(unique(cwd$mort)),
       
       Y = cwd.validation$count,
       location_id = cwd.validation$location_id,
       transect_id = cwd.validation$transect_id.i,
       plot_id = 
         cwd %>%
         group_by(transect_id.i, plot_id.i) %>%
         summarise() %>%
         ungroup() %>%
         arrange(transect_id.i) %>%
         pull(plot_id.i) %>%
         as.numeric(),
       group_id = 
         cwd %>%
         group_by(plot_id.i, group_id) %>%
         summarise() %>%
         ungroup() %>%
         arrange(plot_id.i) %>%
         pull(group_id) %>%
         as.numeric(),
       
       coords = cwd_coords[,c('x_rel', 'y_rel')] %>% as.matrix())

saveRDS(cwdtallies_training_data,
        here::here('02-data', 
                   '05-for_analysis',
                   'cwdtallies_training_data.rds'))

saveRDS(cwdtallies_validation_data,
        here::here('02-data',
                   '05-for_analysis',
                   'cwdtallies_validation_data.rds'))


#### CWD diameters #############################################################

cwd_diams = 
  # load the data
  readRDS(here::here('02-data',
                     '03-clean',
                     'seki_cwd.rds')) %>%
  
  # get mortality level and transect id columns
  left_join(plot_locations.tibble %>%
              select(mort, plot_id),
            by = c('plot_id' = 'plot_id')) %>%
  mutate(transect_id = 
           paste0(plot_id, '-',
                  stringr::str_pad(az, width = 3, side = 'left', pad = 0)))%>%
  mutate(plot_id.i = as.integer(factor(plot_id)),
         group_id = as.integer(factor(mort, levels = c('low', 'mid', 'high'))),
         transect_id.i = as.integer(factor(transect_id))) %>%
  rowid_to_column('obs_id')

set.seed(110819)
cwd_diams.training = 
  cwd_diams %>%
  group_by(mort) %>%
  slice_sample(prop = 0.9) %>%
  ungroup()

cwd_diams.validation = 
  cwd_diams %>%
  filter(!is.element(obs_id, cwd_diams.training$obs_id))

cwddiams_training_data = 
  list(N = nrow(cwd_diams.training),
       P = length(unique(cwd_diams$plot_id)),
       G = length(unique(cwd_diams$mort)),

       group_id = 
         cwd_diams %>%
         group_by(plot_id.i, group_id) %>%
         summarise() %>%
         ungroup() %>%
         arrange(plot_id.i) %>%
         pull(group_id) %>%
         as.numeric(),
       
       plot_id = 
         cwd_diams.training$plot_id.i,
       Y = cwd_diams.training$diam_cm,
       lb = 7.55, 
       ub = max(cwd_diams$diam_cm)+0.1) # could set this to max tree size or something?


cwddiams_validation_data = 
  list(N = nrow(cwd_diams.validation),
       P = length(unique(cwd_diams$plot_id)),
       G = length(unique(cwd_diams$mort)),

       group_id = 
         cwd_diams %>%
         group_by(plot_id.i, group_id) %>%
         summarise() %>%
         ungroup() %>%
         arrange(plot_id.i) %>%
         pull(group_id) %>%
         as.numeric(),
       
       plot_id = 
         cwd_diams.validation$plot_id.i,
       
       Y = cwd_diams.validation$diam_cm,
       lb = 7.55, 
       ub = max(cwd_diams$diam_cm)+0.1)

saveRDS(cwddiams_training_data,
        here::here('02-data',
                   '05-for_analysis',
                   'cwddiams_training_data.rds'))

saveRDS(cwddiams_validation_data,
        here::here('02-data',
                   '05-for_analysis',
                   'cwddiams_validation_data.rds'))




#### understory veg P/A ########################################################

veg = 
  readRDS(here::here('02-data',
                     '03-clean',
                     'seki_veg.rds'))
head(veg)  

# making quadrature points for understory veg presence or absence or 
# height
veg_quadrature = 
  metadata %>%
  select(mort, plot_id, transect_id, az) %>%
  left_join(veg) %>%
  expand(nesting(mort, plot_id, transect_id, az, start_m, end_m, spp, status, height_m, comments),
         sample_point = seq(from = 0.25, to = 14.75, by = 0.5)) %>%
  # veg present is NA at all points if there was no veg on the transect
  mutate(veg_present = start_m <= sample_point & end_m >= sample_point) %>%
  
  # convert NAs to FALSE and get height
  mutate(veg_present = ifelse(is.na(veg_present), FALSE, veg_present),
         veg_height = ifelse(veg_present, height_m, NA)) %>%
  group_by(mort, plot_id, transect_id, az, sample_point) %>%
  summarise(veg_present = any(veg_present),
            veg_height = max(veg_height, na.rm = TRUE)) %>%
  ungroup() %>%

  # get relative spatial coordinates
  mutate(
    veg_height = ifelse(is.finite(veg_height), 
                        veg_height, 
                        NA),
    x_rel = round(sample_point*sin(az*pi/180),2),
       y_rel = round(sample_point*cos(az*pi/180),2),
       plot_id.i = as.integer(factor(plot_id)),
       group_id = as.integer(factor(mort,
                                    levels = c('low', 'mid', 'high')))) %>%
  
  # get height as number of layers
  mutate(height_N = ceiling(veg_height / 0.25))


ggplot(data = veg_quadrature,
       aes(x = as.factor(veg_height)))+
  geom_bar()

ggplot(data = veg_quadrature,
       aes(x = height_N))+
  geom_histogram()+
  facet_grid(mort~.)


veg_coords = 
  veg_quadrature %>%
  group_by(x_rel, y_rel) %>%
  summarise() %>%
  ungroup() %>%
  rowid_to_column('location_id')


veg_quadrature = 
  veg_quadrature %>%
  left_join(veg_coords)

head(veg_quadrature)

ggplot(data = veg_quadrature %>% filter(plot_id=='01'),
       aes(x = x_rel, y = y_rel, color = veg_present))+
  geom_point()


ggplot(data = veg_quadrature %>% filter(plot_id=='01'),
       aes(x = x_rel, y = y_rel, color = veg_height))+
  geom_point()+
  scale_color_viridis_c()

veg_quadrature = 
  veg_quadrature %>%
  rowid_to_column('obs_id')

set.seed(110819)
veg_pa.training = 
  veg_quadrature %>%
  group_by(mort) %>%
  slice_sample(prop = 0.9) %>%
  ungroup()

veg_pa.validation = 
  veg_quadrature %>%
  filter(!is.element(obs_id, veg_pa.training$obs_id))

vegpa_training_data = 
  list(N = nrow(veg_pa.training),
       P = length(unique(veg_quadrature$plot_id)),
       L = nrow(veg_coords),
       G = length(unique(veg_quadrature$mort)),
       coords = veg_coords[,c('x_rel', 'y_rel')] %>% as.matrix(),
       plot_id = veg_pa.training$plot_id.i,
       location_id = veg_pa.training$location_id,
       group_id = 
         veg_quadrature %>%
         group_by(plot_id.i, group_id) %>%
         summarise() %>%
         ungroup() %>%
         arrange(plot_id.i) %>%
         pull(group_id) %>%
         as.numeric(),
       Y = as.integer(veg_pa.training$veg_present))

vegpa_validation_data = 
  list(N = nrow(veg_pa.validation),
       P = length(unique(veg_quadrature$plot_id)),
       L = nrow(veg_coords),
       G = length(unique(veg_quadrature$mort)),
       coords = veg_coords[,c('x_rel', 'y_rel')] %>% as.matrix(),
       plot_id = veg_pa.validation$plot_id.i,
       location_id = veg_pa.validation$location_id,
       group_id = 
         veg_quadrature %>%
         group_by(plot_id.i, group_id) %>%
         summarise() %>%
         ungroup() %>%
         arrange(plot_id.i) %>%
         pull(group_id) %>%
         as.numeric(),
       Y = as.integer(veg_pa.validation$veg_present))

saveRDS(
  vegpa_training_data,
  here::here('02-data',
                   '05-for_analysis',
                   'vegpa_training_data.rds'))

saveRDS(
  vegpa_validation_data,
  here::here('02-data',
                   '05-for_analysis',
                   'vegpa_validation_data.rds'))




#### veg height ################################################################

veg_height = 
  veg_quadrature %>%
  filter(!is.na(height_N)) %>%
  mutate(plot_id.i = as.integer(factor(plot_id)))

veg_height_data = 
    list(N = nrow(veg_height),
       M = as.integer(max(veg_height$veg_height,na.rm = TRUE)/0.25),
       P = length(unique(veg_height$plot_id)),
       L = nrow(veg_coords),
       G = length(unique(veg_height$mort)),
       J = length(unique(veg_height$mort)),
       
       coords = veg_coords[,c('x_rel', 'y_rel')] %>% as.matrix(),
       X = matrix(nrow = nrow(veg_height), 
                  ncol = length(unique(veg_height$mort)),
                  byrow = FALSE,
                  data = c(rep(1, times = nrow(veg_height)),
                           as.integer(veg_height$mort=='mid'),
                           as.integer(veg_height$mort=='high'))),
       plot_id = veg_height$plot_id.i,
       location_id = veg_height$location_id,
       group_id = 
         veg_height %>%
         group_by(plot_id.i, group_id) %>%
         summarise() %>%
         ungroup() %>%
         arrange(plot_id.i) %>%
         pull(group_id) %>%
         as.numeric(),
       Y = veg_height$height_N)

#### trees #####################################################################



# make a quadrature of 1m2 cells for each plot
trees_quadrature = 
  expand.grid(x_rel = seq(-14.5, 14.5, 1),
              y_rel = seq(-14.5, 14.5, 1)) %>%
  filter((x_rel >= -4.5 & x_rel <= 4.5) | (y_rel >= -4.5 & y_rel <= 4.5)) %>%
  expand(nesting(x_rel, y_rel),
           plot_locations.tibble %>%
               group_by(plot_id, mort) %>%
               summarise() %>%
               ungroup())  %>%
  left_join(
    readRDS(here::here('02-data', '03-clean', 'seki_trees.rds')) %>%
      filter(dbh_cm >= 11.4) %>%
      mutate(count = 1) %>%
        mutate(x_rel = ifelse(is.element(ts, c('NE', 'SE')),
                            distance_m,
                            ifelse(is.element(ts, c('NW', 'SW')),
                                   -distance_m,
                                   ifelse(is.element(ts, c('EN', 'ES')),
                                          location_m,
                                          ifelse(is.element(ts, c('WN', 'WS')),
                                                 -location_m,
                                                 NA)))),
                 y_rel = ifelse(is.element(ts, c('NE', 'NW')),
                            location_m,
                            ifelse(is.element(ts, c('SE', 'SW')),
                                   -location_m,
                                   ifelse(is.element(ts, c('EN', 'WN')),
                                          distance_m,
                                          ifelse(is.element(ts, c('ES', 'WS')),
                                                 -distance_m,
                                                 NA))))) %>%
      select(plot_id, x_rel, y_rel, count) %>%
      mutate(x_rel = floor(x_rel)+0.5,
             y_rel = floor(y_rel)+0.5) %>%
      filter((x_rel >= -14.5 & x_rel <= 14.5 & y_rel >= -4.5 & y_rel <= 4.5)|
               (x_rel >= -4.5 & x_rel <= 4.5 & y_rel >= -14.5 & y_rel <= 14.5)) %>%
      group_by(plot_id, x_rel, y_rel) %>%
      summarise(count = sum(count)) %>%
      ungroup()) %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  rowid_to_column('obs_id') %>%
  mutate(group_id = ifelse(mort=='low', 
                           1, 
                           ifelse(mort == 'mid', 2, 3)))

head(trees_quadrature)

ggplot(trees_quadrature,
       aes(x = x_rel, y = y_rel, fill = count))+
  geom_tile()+
  coord_fixed()+
  theme_minimal()+
  facet_wrap(~plot_id)


trees_coords = 
  trees_quadrature %>%
  group_by(x_rel, y_rel) %>%
  summarise() %>%
  ungroup() %>%
  rowid_to_column('location_id')

trees_quadrature = 
  trees_quadrature %>%
  left_join(trees_coords) %>%
  mutate(plot_id.i = as.integer(factor(plot_id)))


set.seed(110819)
trees_quadrature.training = 
  trees_quadrature %>%
  group_by(group_id) %>%
  slice_sample(prop = 0.9) %>%
  ungroup()

trees_quadrature.validation = 
  trees_quadrature %>%
  filter(!is.element(obs_id, trees_quadrature.training$obs_id))

trees_training_data = 
  list(N = nrow(trees_quadrature.training),
       P = length(unique(trees_quadrature$plot_id)),
       G = length(unique(trees_quadrature$group_id)),
       L = nrow(trees_coords),
       coords = trees_coords[,c('x_rel', 'y_rel')] %>% as.matrix,
       location_id = trees_quadrature.training$location_id,
       plot_id = trees_quadrature.training$plot_id.i,
       group_id = 
         trees_quadrature %>%
         select(group_id, plot_id.i) %>%
         group_by(group_id, plot_id.i) %>%
         summarise() %>%
         ungroup() %>%
         arrange(plot_id.i) %>%
         pull(group_id),
       Y = trees_quadrature.training$count)

trees_validation_data = 
  list(N = nrow(trees_quadrature.validation),
       P = length(unique(trees_quadrature$plot_id)),
       G = length(unique(trees_quadrature$group_id)),
       L = nrow(trees_coords),
       coords = trees_coords[,c('x_rel', 'y_rel')] %>% as.matrix,
       location_id = trees_quadrature.validation$location_id,
       plot_id = trees_quadrature.validation$plot_id.i,
       group_id = 
         trees_quadrature %>%
         select(group_id, plot_id.i) %>%
         group_by(group_id, plot_id.i) %>%
         summarise() %>%
         ungroup() %>%
         arrange(plot_id.i) %>%
         pull(group_id),
       Y = trees_quadrature.validation$count)

saveRDS(trees_training_data,
        here::here('02-data',
                   '05-for_analysis',
                   'trees_training_data.rds'))

saveRDS(trees_validation_data,
        here::here('02-data',
                   '05-for_analysis',
                   'trees_validation_data.rds'))


#### saplings ##################################################################



# make a quadrature of 1m2 cells for each plot
saplings_quadrature = 
  expand.grid(x_rel = seq(-14.5, 14.5, 1),
              y_rel = seq(-14.5, 14.5, 1)) %>%
  filter((x_rel >= -0.5 & x_rel <= 0.5) | (y_rel >= -0.5 & y_rel <= 0.5)) %>%
  expand(nesting(x_rel, y_rel),
           plot_locations.tibble %>%
               group_by(plot_id, mort) %>%
               summarise() %>%
               ungroup())  %>%
  left_join(
    readRDS(here::here('02-data', '03-clean', 'seki_trees.rds')) %>%
      filter(dbh_cm < 11.4) %>%
      mutate(count = 1) %>%
        mutate(x_rel = ifelse(is.element(ts, c('NE', 'SE')),
                            distance_m,
                            ifelse(is.element(ts, c('NW', 'SW')),
                                   -distance_m,
                                   ifelse(is.element(ts, c('EN', 'ES')),
                                          location_m,
                                          ifelse(is.element(ts, c('WN', 'WS')),
                                                 -location_m,
                                                 NA)))),
                 y_rel = ifelse(is.element(ts, c('NE', 'NW')),
                            location_m,
                            ifelse(is.element(ts, c('SE', 'SW')),
                                   -location_m,
                                   ifelse(is.element(ts, c('EN', 'WN')),
                                          distance_m,
                                          ifelse(is.element(ts, c('ES', 'WS')),
                                                 -distance_m,
                                                 NA))))) %>%
      select(plot_id, x_rel, y_rel, count) %>%
      mutate(x_rel = floor(x_rel)+0.5,
             y_rel = floor(y_rel)+0.5) %>%
      filter((x_rel >= -14.5 & x_rel <= 14.5 & y_rel >= -0.5 & y_rel <= 0.5)|
               (x_rel >= -0.5 & x_rel <= 0.5 & y_rel >= -14.5 & y_rel <= 14.5)) %>%
      group_by(plot_id, x_rel, y_rel) %>%
      summarise(count = sum(count)) %>%
      ungroup()) %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  rowid_to_column('obs_id') %>%
  mutate(group_id = ifelse(mort=='low', 
                           1, 
                           ifelse(mort == 'mid', 2, 3)))

head(saplings_quadrature)

ggplot(saplings_quadrature,
       aes(x = x_rel, y = y_rel, fill = count))+
  geom_tile()+
  coord_fixed()+
  theme_minimal()+
  facet_wrap(~plot_id)


saplings_coords = 
  saplings_quadrature %>%
  group_by(x_rel, y_rel) %>%
  summarise() %>%
  ungroup() %>%
  rowid_to_column('location_id')

saplings_quadrature = 
  saplings_quadrature %>%
  left_join(saplings_coords) %>%
  mutate(plot_id.i = as.integer(factor(plot_id)))


set.seed(110819)
saplings_quadrature.training = 
  saplings_quadrature %>%
  group_by(group_id) %>%
  slice_sample(prop = 0.9) %>%
  ungroup()

saplings_quadrature.validation = 
  saplings_quadrature %>%
  filter(!is.element(obs_id, saplings_quadrature.training$obs_id))

saplings_training_data = 
  list(N = nrow(saplings_quadrature.training),
       P = length(unique(saplings_quadrature$plot_id)),
       G = length(unique(saplings_quadrature$group_id)),
       L = nrow(saplings_coords),
       coords = saplings_coords[,c('x_rel', 'y_rel')] %>% as.matrix,
       location_id = saplings_quadrature.training$location_id,
       plot_id = saplings_quadrature.training$plot_id.i,
       group_id = 
         saplings_quadrature %>%
         select(group_id, plot_id.i) %>%
         group_by(group_id, plot_id.i) %>%
         summarise() %>%
         ungroup() %>%
         arrange(plot_id.i) %>%
         pull(group_id),
       Y = saplings_quadrature.training$count)

saplings_validation_data = 
  list(N = nrow(saplings_quadrature.validation),
       P = length(unique(saplings_quadrature$plot_id)),
       G = length(unique(saplings_quadrature$group_id)),
       L = nrow(saplings_coords),
       coords = saplings_coords[,c('x_rel', 'y_rel')] %>% as.matrix,
       location_id = saplings_quadrature.validation$location_id,
       plot_id = saplings_quadrature.validation$plot_id.i,
       group_id = 
         saplings_quadrature %>%
         select(group_id, plot_id.i) %>%
         group_by(group_id, plot_id.i) %>%
         summarise() %>%
         ungroup() %>%
         arrange(plot_id.i) %>%
         pull(group_id),
       Y = saplings_quadrature.validation$count)

saveRDS(saplings_training_data,
        here::here('02-data',
                   '05-for_analysis',
                   'saplings_training_data.rds'))

saveRDS(saplings_validation_data,
        here::here('02-data',
                   '05-for_analysis',
                   'saplings_validation_data.rds'))



