
#### setup #####################################################################

library(here)
library(sf)
library(tidyverse)

#### clean up plot centers #####################################################

# the true plot center locations were recorded into several 
# different files depending on the chronology of when the 
# plot was visited; just keeping notes here while i peer at them 
# in arcgis, excel, and datasheets...
# the chain of events is:
#  a) late may plot_locations_arc created and loaded into arcgis
#  b) week of june 1: record true locations for plots 23 and 2, 
#     create plot_locations_redo_arc (true locatinos of plots 2 and 23, 
#     and new target locations for plots 1-17 respecting a distance-from
#     road mask)
#  c) approx june 7-14: record true locatinos for plots 
#    1, 3, 4, 5, 7, 8, 9, 10, 14, 15, 16
#  d) approx june 18: create target locations for plots 
#    25, 30, and 35 (plus other unvisited plots)
#  e) approx june 18: create file plot_loc_061821 (true locations for 
#    plots 1, 2, 3, 4, 5, 7, 8, 9, 10, 14, 15, 16, 23
#  f) approx june 22-july1: record true locations for plots
#     6 (in avenza/excel), 11 (backup_ptge), 12 (in excel), 13 (in excel and backup_ptge), 17 (backup_ptge), 25, 30, 35
#  g) aug 8: create files seki_fuels_080921, seki_fuels_backup_plotloc_080921,
#     seki_fuels_backup_ptge_080921, and seki_fuels_additional_080921 
#     pulling from the trimble
#  h) aug 8: this file written

# PLOT | target location file | true location file(s)
# 1|plot_locatinos_redo_arc, seki_fuels_backup_plotloc_080821 | plot_loc_061821, seki_fuels_080921
# 2 | plot_locatinos_arc | plot_locations_redo_arc, plot_loc_061821, seki_fuels_080921, seki_fuels_backup_plotloc_080921
# 3| plot_locatinos_redo_arc, seki_fuels_backup_plotloc_080821 | plot_loc_061821, seki_fuels_080921
# 4| plot_locatinos_redo_arc, seki_fuels_backup_plotloc_080821 | plot_loc_061821, seki_fuels_080921
# 5| plot_locatinos_redo_arc, seki_fuels_backup_plotloc_080821 | plot_loc_061821, seki_fuels_080921
# 6| plot_locatinos_redo_arc, plot_loc_061821, seki_fuels_080921, seki_fuels_backup_plotloc_080921  | datasheet_digital_seki_7-18
# 7| plot_locatinos_redo_arc, seki_fuels_backup_plotloc_080821 | plot_loc_061821, seki_fuels_080921
# 8| plot_locatinos_redo_arc, seki_fuels_backup_plotloc_080821 | plot_loc_061821, seki_fuels_08092
# 9| plot_locatinos_redo_arc, seki_fuels_backup_plotloc_080821 | plot_loc_061821, seki_fuels_08092
# 10| plot_locatinos_redo_arc, seki_fuels_backup_plotloc_080821 | plot_loc_061821, seki_fuels_08092
# 11| plot_locatinos_redo_arc, plot_loc_061821, seki_fuels_080921, seki_fuels_backup_plotloc_080921 | seki_Fuels_backup_ptge_080921
# 12| plot_locatinos_redo_arc, plot_loc_061821, seki_fuels_080921, seki_fuels_backup_plotloc_080921 | datasheet_digital_seki_7-18_completed
# 13| plot_locatinos_redo_arc, plot_loc_061821, seki_fuels_080921, seki_fuels_backup_plotloc_080921 | seki_fuels_backup_ptgen_080921 
# 14| plot_locatinos_redo_arc, seki_fuels_backup_plotloc_080821 | plot_loc_061821, seki_fuels_080921
# 15| plot_locatinos_redo_arc, seki_fuels_backup_plotloc_080821 | plot_loc_061821, seki_fuels_080921
# 16| plot_locatinos_redo_arc, seki_fuels_backup_plotloc_080821 | plot_loc_061821, seki_fuels_080921
# 17| plot_locatinos_redo_arc, plot_loc_061821, seki_fuels_080921, seki_fuels_backup_plotloc_080921 | seki_fuels_backup_ptgen_080921 
# 23 | plot_locations_arc | plot_locations_redo_arc, plot_loc_061821, seki_fuels_080921, seki_fuels_backup_plotloc_080921
# 25 | plot_locations_additional_arc | seki_fuels_added_plotloc_080921
# 30 | plot_locations_additional_arc |  seki_fuels_added_plotloc_080921
# 35 | plot_locations_additional_arc |  seki_fuels_added_plotloc_080921

# first load all these different shapefiles...


# target locations for plots 1-17, 23 and 25, 30, 35, respectively
# the base versions of these files were created in R and not recording the 
# projection correctly when writing, so projection files were added in arc
plot_locations_redo_arc = 
  st_read(here::here('02-data',
                     '02-intermediate',
                     'plot_locations_redo_arc.shp')) %>%
  mutate(mort = as.character(mort),
         plot_id = as.character(plot_id))

plot_locations_additional_arc = 
  st_read(here::here('02-data',
                     '02-intermediate',
                     'plot_locations_additional_arc.shp'))%>%
  mutate(mort = as.character(mort),
         plot_id = as.character(plot_id))

# true locations:
# plot_locations_redo_arc with some plot locations modified in the field using 
# the trimble to match the true location, rather than the target location
seki_fuels_080921 = 
  st_read(here::here('02-data',
                     '02-intermediate',
                     'recorded_plots_080921',
                     'seki-fuels',
                     'plot_loc.shp')) %>%
  st_set_crs(st_crs(plot_locations_redo_arc))%>%
  mutate(mort = as.character(mort),
         plot_id = as.character(plot_id))

seki_fuels_backup_ptge_080921 = 
  st_read(here::here('02-data',
                     '02-intermediate',
                     'recorded_plots_080921',
                     'seki-fuels-backup',
                     'Point_ge.shp')) %>%
  st_set_crs(st_crs(plot_locations_redo_arc))%>%
  mutate(Comment = as.character(Comment))

seki_fuels_added_plotloc_080921 = 
  st_read(here::here('02-data',
                     '02-intermediate',
                     'recorded_plots_080921',
                     'seki-fuels-added',
                     'plot_loc.shp')) %>%
  st_set_crs(st_crs(plot_locations_redo_arc))%>%
  mutate(mort = as.character(mort),
         plot_id = as.character(plot_id))


target_locations.tibble = 
  bind_rows(plot_locations_redo_arc,
        plot_locations_additional_arc) %>%
  filter(is.element(plot_id,
                    c(1:17, 23, 25, 30, 35))) %>%
  mutate(x_target = st_coordinates(.)[,'X'],
         y_target = st_coordinates(.)[,'Y']) %>%
  as_tibble() %>%
  select(-geometry)


true_locations.tibble = 
  seki_fuels_080921 %>%
  filter(is.element(plot_id, 
                    c(1:5, 7:10, 14:16, 23))) %>%
  select(-mort) %>%
  bind_rows(seki_fuels_added_plotloc_080921 %>%
              filter(is.element(plot_id,
                                c(25, 30, 35))) %>%
              select(-mort)) %>%
  bind_rows(seki_fuels_backup_ptge_080921 %>%
              mutate(plot_id = gsub(x = Comment,
                                    pattern = 'plot | actual',
                                    replacement='')) %>%
              select(-Comment)) %>%
  mutate(x_true = st_coordinates(.)[,'X'],
         y_true = st_coordinates(.)[,'Y']) %>%
  as_tibble() %>%
  select(-geometry) %>%
  bind_rows(data.frame(plot_id = c('6', '12'),
                       x_true = c(337482, 337744),
                       y_true = c(4047635, 4047769)) %>%
              mutate(plot_id = as.character(plot_id)) %>%
              as_tibble())

plot_locations.tibble = 
  left_join(target_locations.tibble,
            true_locations.tibble,
            by = c('plot_id')) %>%
  mutate(plot_id = str_pad(plot_id, width = 2, side = 'left', pad = '0'))

plot_locations.tibble

plot_locations.sf = 
  plot_locations.tibble %>%
  st_as_sf(coords = c('x_true', 'y_true'),
           crs = st_crs(seki_fuels_080921))

aoi = 
  st_read(here::here('02-data',
                     '00-source',
                     'seki_dx',
                     'dxstudyarea.shp')) %>%
  st_transform(crs = st_crs(plot_locations.sf))

ggplot()+
  geom_sf(data = aoi, fill = NA)+
  geom_sf(data = plot_locations.sf,
          aes(color = mort))+
  theme_bw()

#### geolocate metadata ########################################################

seki.metadata = readRDS(here::here('02-data', '03-clean', 'seki_metadata.rds'))

seki.metadata = 
  seki.metadata %>%
  left_join(plot_locations.tibble %>%
              select(mort, plot_id, x_center = x_true, y_center = y_true),
            by = c('plot_id'))

#### geolocate litter and duff #################################################

seki.litterduff = 
  readRDS(here::here('02-data',
                     '03-clean',
                     'seki_litterduff.rds'))

plot_locations.tibble

check_translation = 
  data.frame(compass_az = c(0, 90, 180, 270)) %>%
  mutate(shifted_az = compass_az+90,
         shifted_radian = shifted_az*(pi/180),
         shifted_cos = round(cos(shifted_radian), 3),
         shifted_sin = round(sin(shifted_radian), 3),
         old_x = round(sin(compass_az*(pi/180)), 3),
         old_y = round(cos(compass_az*(pi/180)), 3))

# need to just use sin for x and cos for y on the original angle, 
# rather than try to translate the angle from compass-degres (which increase
# in the clockwise direction) to radians (which increase in the counterclockwise
# direction); not smart enough to translate on the frontend and for some reason
# the hack works
check_translation

seki.litterduff = 
  seki.litterduff %>%
  left_join(.,
            plot_locations.tibble %>%
              select(plot_id, mort, x_center = x_true, y_center = y_true)) %>%
  mutate(
    # not really sure why this works but it does translate between compass-angle
    # and unit-circle angle correctly to just use sin for the x and cos for the 
    # y...
    x_coord = x_center + (location_m*sin(az*(pi/180))),
    y_coord = y_center + (location_m*cos(az*(pi/180)))
  ) %>%
  select(-x_center, -y_center)



#### geolocate fwd #############################################################

seki.fwd_tallies = 
  readRDS(here::here('02-data',
                     '03-clean',
                     'seki_fwdtallies.rds'))

seki.fwd_tallies = 
  seki.fwd_tallies %>%
  left_join(.,
            plot_locations.tibble %>%
              select(plot_id, mort, x_center = x_true, y_center = y_true)) %>%
  mutate(x_coord = x_center + location_m*sin(az*(pi/180)),
         y_coord = y_center + location_m*cos(az*(pi/180))) %>%
  select(-x_center, -y_center)

ggplot(data = seki.fwd_tallies %>% filter(plot_id=='01'),
       aes(x = x_coord, y = y_coord, color = transect_id))+
  geom_point()

#### geolocate cwd #############################################################

seki.cwd = 
  readRDS(here::here('02-data',
                     '03-clean',
                     'seki_cwd.rds')) %>%
  left_join(.,
            plot_locations.tibble %>%
              select(plot_id, mort, x_center = x_true, y_center = y_true))
  


head(seki.cwd)

seki.cwd %>%
  group_by(plot_id, mort) %>%
  summarise(ssd_cm2 = sum(diam_cm**2)) %>%
  ungroup() %>%
  ggplot(aes(x = mort, y = ssd_cm2))+
  geom_point()





#### geolocate veg #############################################################

seki.veg = 
  readRDS(here::here('02-data',
                     '03-clean',
                     'seki_veg.rds'))

seki.veg = 
  seki.veg %>%
  left_join(.,
            plot_locations.tibble %>%
              select(plot_id, mort, x_center = x_true, y_center = y_true))

#### geolocate trees ###########################################################

seki.trees = 
  readRDS(here::here('02-data',
                     '03-clean',
                     'seki_trees.rds'))

seki.trees = 
  seki.trees %>%
  left_join(.,
            plot_locations.tibble %>%
              select(plot_id, mort, x_center = x_true, y_center = y_true)) %>%
  # translate each tree's location according to its location along the 
  # transect, distance from the transect, and quadrant
  mutate(x_coord = ifelse(is.element(ts, c('NE', 'SE')),
                    x_center+distance_m,
                    ifelse(is.element(ts, c('NW', 'SW')),
                           x_center-distance_m,
                           ifelse(is.element(ts, c('EN', 'ES')),
                                  x_center+location_m,
                                  ifelse(is.element(ts, c('WN', 'WS')),
                                         x_center-location_m,
                                         NA)))),
         y_coord = ifelse(is.element(ts, c('NE', 'NW')),
                    y_center+location_m,
                    ifelse(is.element(ts, c('SE', 'SW')),
                           y_center-location_m,
                           ifelse(is.element(ts, c('EN', 'WN')),
                                  y_center+distance_m,
                                  ifelse(is.element(ts, c('ES', 'WS')),
                                         y_center-distance_m,
                                         NA))))) %>%
  select(-x_center, -y_center)


#### write results #############################################################

st_write(plot_locations.sf,
         here::here('02-data',
                    '04-geolocated',
                    'seki_plot_locations.shp'),
         delete_dsn=TRUE)
write.csv(plot_locations.tibble,
          here::here('02-data',
                     '04-geolocated',
                     'seki_plot_locations.csv'),
          row.names = FALSE)
saveRDS(plot_locations.sf,
        here::here('02-data',
                   '04-geolocated',
                   'seki_plot_locations.rds'))

write.csv(seki.cwd,
          here::here('02-data',
                     '04-geolocated',
                     'seki_cwd.csv'),
          row.names = FALSE)
saveRDS(seki.cwd,
        here::here('02-data',
                   '04-geolocated',
                   'seki_cwd.rds'))

write.csv(seki.fwd_tallies,
          here::here('02-data',
                     '04-geolocated',
                     'seki_fwdtallies.csv'),
          row.names = FALSE)
saveRDS(seki.fwd_tallies,
        here::here('02-data',
                   '04-geolocated',
                   'seki_fwdtallies.rds'))

write.csv(seki.litterduff,
          here::here('02-data',
                     '04-geolocated',
                     'seki_litterduff.csv'),
          row.names = FALSE)
saveRDS(seki.litterduff,
        here::here('02-data',
                   '04-geolocated',
                   'seki_litterduff.rds'))

write.csv(seki.veg,
          here::here('02-data',
                     '04-geolocated',
                     'seki_veg.csv'),
          row.names = FALSE)
saveRDS(seki.veg,
        here::here('02-data',
                   '04-geolocated',
                   'seki_veg.rds'))

write.csv(seki.trees,
          here::here('02-data',
                     '04-geolocated',
                     'seki_trees.csv'),
          row.names = FALSE)
saveRDS(seki.trees,
        here::here('02-data',
                   '04-geolocated',
                   'seki_trees.rds'))

write.csv(seki.metadata,
          here::here('02-data',
                     '04-geolocated',
                     'seki_metadata.csv'),
          row.names = FALSE)
saveRDS(seki.metadata,
        here::here('02-data',
                   '04-geolocated',
                   'seki_metadata.rds'))
