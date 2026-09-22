
#### thinking about CWD ########################################################

library(here)
library(sf)
library(tidyverse)

# dataframe with plot centers
plot_locations = 
  readRDS(here::here('02-data',
                     '04-geolocated',
                     'seki_plot_locations.rds')) %>%
  mutate(x_actual = st_coordinates(.)[,'X'],
         y_actual = st_coordinates(.)[,'Y']) %>%
  as.data.frame() %>%
  select(-geometry)


head(plot_locations)

# get the pairwise distances between plots
dist(plot_locations[,c('x_actual', 'y_actual')]) %>%
  as.matrix() %>%
  as.vector() %>%
  summary()


cwd = readRDS(here::here('02-data',
                         '04-geolocated',
                         'seki_cwd.rds'))

head(cwd)

cwd %>%
  ggplot(aes(x = location_m, size = diam_cm, y = transect_id))+
  geom_point()+
  scale_x_continuous(limits = c(0, 30))+
  facet_grid(mort~., scales = 'free_y')

cwd %>%
  ggplot(aes(x = diam_cm, color = mort))+
  geom_density()


cwd %>%
  mutate(location_m = floor(location_m)+0.5,
         count = 1) %>%
  select(mort, plot_id, az, location_m, count) %>%
  bind_rows(
    plot_locations %>%
      select(mort, plot_id) %>%
      expand(nesting(mort, plot_id),
             az = c(0, 90, 180, 270),
             location_m = 0:29+0.5) %>%
      mutate(count = 0)
  ) %>%
  group_by(mort, plot_id, az, location_m) %>%
  summarise(count = sum(count)) %>%
  ungroup() %>%
  mutate(transect_id = paste0(plot_id, str_pad(az, 3, 'left', 0))) %>%
  ggplot()+
  geom_segment(aes(x = location_m-0.5, y = transect_id, 
                   xend = location_m+0.5, yend = transect_id,
                   color = count),
               lwd = 3)+
  scale_color_viridis_c()+
  theme_minimal()


#### comparing the DX fuel loads against MMI ###################################


library(here)
library(sf)
library(tidyverse)
library(raster)
dx_fuels_2021 = 
  read.csv(here::here('02-data',
                      '00-source',
                      'dx',
                      'SEKI_fuels_processed_2021.csv'))

head(dx_fuels_2021)

dx_plots = 
  st_read(here::here('02-data',
                     '00-source',
                     'seki_dx',
                     'DX_plots_29orig.shp'))

length(unique(dx_fuels_2021$plot_id))

head(dx_plots)

head(dx_fuels_2021)

dx_plots = 
  dx_plots %>%
  left_join(dx_fuels_2021 %>%
              group_by(plot_id) %>%
              summarise(duff_depth_cm = mean(duff_depth_cm, na.rm = TRUE),
                        litter_depth_cm = mean(litter_depth_cm, na.rm = TRUE),
                        count_x1h = mean(count_x1h, na.rm = TRUE),
                        count_x10h = mean(count_x10h, na.rm = TRUE),
                        count_x100h = mean(count_x100h, na.rm = TRUE),
                        sum_d2_1000h = mean(sum_d2_1000r_cm2+sum_d2_1000s_cm2, na.rm = TRUE)) %>%
              ungroup() %>%
              mutate(plot_id = as.numeric(plot_id)),
            by = c('Plot' = 'plot_id'))


mmi = raster::raster(here::here('02-data',
                                '00-source',
                                'edart',
                                'mmi_sum5_sc408ss_2016.bsq'))

dx_plots$mmi = raster::extract(mmi, dx_plots)

dx_plots %>%
  as.data.frame() %>%
  dplyr::select(-geometry) %>%
  pivot_longer(cols = c(duff_depth_cm,
                        litter_depth_cm,
                        count_x1h,
                        count_x10h,
                        count_x100h,
                        sum_d2_1000h),
               names_to = 'component',
               values_to = 'value') %>%
  ggplot(aes(x = mmi, y = value))+
  geom_point()+
  geom_smooth(method = 'lm')+
  facet_wrap(~component, scales = 'free')+
  theme_minimal()

dx_overstory_2017 = 
  read.csv(here::here('02-data',
                      '00-source',
                      'dx',
                      'SEKI_treelist_processed_2017.csv'))

head(dx_overstory_2017)

head(seki_trees_2021)

seki_trees_2021 = 
  readRDS(here::here('02-data',
                     '05-for_analysis',
                     'biomass',
                     'seki_trees.rds')) %>%
  filter(status=='L' | decay_class=='1') %>%
  group_by(plot_id, status) %>%
  summarise(ba_m2ha = sum(ba_m2ha)) %>%
  ungroup() %>%
  pivot_wider(names_from = c('status'),
              values_from = c('ba_m2ha')) %>%
  mutate(proportion_dead = D/(D+L)) %>%
  left_join(
    readRDS(here::here('02-data',
                       '04-geolocated',
                       'seki_plot_locations.rds')) %>%
      mutate(x_center = st_coordinates(.)[,'X'],
             y_center = st_coordinates(.)[,'Y']) %>%
      as.data.frame() %>%
      dplyr::select(-geometry)
  ) %>%
  st_as_sf(coords = c('x_center', 'y_center'),
           crs = readRDS(here::here('02-data',
                                    '04-geolocated',
                                    'seki_plot_locations.rds')) %>%
             st_crs())

head(seki_trees_2021)

seki_trees_2021$mmi = raster::extract(mmi, seki_trees_2021)

ggplot(seki_trees_2021,
       aes(x = mmi/100, y = proportion_dead))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0, color = 'red')+
  geom_smooth(method = 'lm')+
  theme_minimal()+
  coord_fixed()


# still weird. It looks like MMI does a good job predicting recently dead, but 
# does not track with fuel loads in either the spatial fuels or the DX data. 


#### linear point pattern analysis of cwd ######################################

library(here)
library(tidyverse)
library(spatstat)


cwd = readRDS(here::here('02-data',
                         '04-geolocated',
                         'seki_cwd.rds'))

head(cwd)
help("linnet")
help("as.ppp")

transect_linear_network = 
  linnet(vertices = 
           ppp(x = c(0, 30), y= c(0, 0),
               window = owin(xrange = c(-1, 31), yrange = c(-1, 1))),
         m = matrix(nrow = 2, ncol = 2, data = c(FALSE, TRUE, TRUE, FALSE), byrow = FALSE))

help("lpp")

plot(transect_linear_network)

cwd_lpp = 
  lapply(X = unique(cwd$transect_id),
         FUN = function(transect){
           on_transect = cwd %>% filter(transect_id==transect)
           coords = 
             matrix(nrow = nrow(on_transect),
                    ncol = 2,
                    data = 
                      c(on_transect$location_m,
                        rep(0, times = nrow(on_transect))),
                    byrow = FALSE)
           
           linear_pattern = 
             lpp(X = coords, L = transect_linear_network)
           
           return(linear_pattern)
         })

head(cwd_lpp)

help("density.lpp")

density(x = cwd_lpp[[1]])
summary(cwd_lpp[[1]])
plot(linearpcf(X = cwd_lpp[[1]], r = 0:15))

cwd_pcf = 
  lapply(X = 1:length(cwd_lpp),
         FUN = function(transect){
           
           tryCatch({
             linearpcf(X = cwd_lpp[[transect]], r = 0:15)
           },
           error = function(cond){
             print(cond)
             return(NA)
           })
         })

cwd_pcf = cwd_pcf[!is.na(cwd_pcf)]


cwd_envelopes = 
  lapply(X = 1:length(cwd_lpp),
         FUN = function(transect){
           envelope(cwd_lpp[[transect]], nsim = 199, nrank = 5)
         })

lapply(cwd_envelopes, plot)



#### union of DX and spatial plots #############################################


union_plots = 
  readRDS(here::here('02-data',
                     '04-geolocated',
                     'seki_plot_locations.rds')) %>%
  mutate(x_actual = st_coordinates(.)[,'X'],
         y_actual = st_coordinates(.)[,'Y']) %>%
  as.data.frame() %>%
  mutate(type = 'spatial') %>%
  select(-geometry) %>%
  select(plot_id, type, x_coord = x_actual, y_coord = y_actual)  %>%
  bind_rows(
    st_read(here::here('02-data',
                     '00-source',
                     'seki_dx',
                     'DX_plots_29orig.shp')) %>%
      st_transform(crs = st_crs(readRDS(here::here('02-data',
                                                   '04-geolocated',
                                                   'seki_plot_locations.rds')))) %>%
      mutate(x_coord = st_coordinates(.)[,'X'],
             y_coord = st_coordinates(.)[,'Y']) %>%
      as.data.frame() %>%
      mutate(type = 'dx',
             Plot = as.character(Plot)) %>%
      select(-geometry) %>%
      select(plot_id = Plot, type, x_coord, y_coord)
  )

head(union_plots)

data.frame(distance_m = 
             dist(union_plots) %>%
             as.matrix() %>%
             as.vector()) %>%
  filter(distance_m >0) %>%
  ggplot(aes(x = distance_m))+
  geom_histogram()
  

