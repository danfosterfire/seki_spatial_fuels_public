# Decided to limit plot placement to within 600m of a road or trail because of 
# difficult access during the first week. We did a 600m bushwack to get to a 
# plot and it took about an hour each way, having 6 plots which are > 1km from 
# a trail (and over a ridge) was too problematic in terms of crew fatigue and 
# being able to do a complete plot in a single workday. This script generates 
# combines the 2 plots already completed from the initial draw with a new set 
# of plot locations, generated with the following criteria:
#  - less than 25 degree slope at plot center
#  - within 600m of a road or trail
#  - at least 50m from a road or trail
#  - in mixed-conifer (pine) forest type from the calveg layer
#  - within the DX crystal cave AOI
#  - stratified by mortality class
# The script also includes a comparison of the AOI defined by these criteria 
# (with the road/trail restriction) with the an AOI without the road/trail 
# restrictions, to check whether the restriction causes a bias in terms of the 
# slope/aspect/elevation. 


library(here)
library(raster)
library(sf)
library(tidyverse)
library(elevatr)
library(ggplot2)
library(tmap)


#### load data #################################################################

# crystal cave DX AOI, defined using an elevation band around crystal cave road
# and environs
dx_aoi = st_read(here::here('02-data',
                         '00-source',
                         'seki_dx',
                         'dxstudyarea.shp')) 

# write this out for use in arcgis
st_write(dx_aoi %>%
           st_buffer(100) %>%
           st_bbox() %>%
           st_as_sfc(),
         here::here('02-data',
                    '02-intermediate',
                    'intermediate_plot_location_files',
                    'buffered_aoi.shp'),
         delete_dsn = TRUE)

# polyline tracing crystal cave road and mapped trails coming off the road
# (hand traced from the ESRI topo map of the AOI). I've excluded the trail
# running east from crystal cave road up to sunset rock because the terrain 
# beneath the trail is not passable, but the canyon is so steep that the 
# far side of it is within 600m as the crow flies and it was getting included 
# in the "accessible area" even though it's not accessible.
access = st_read(here::here('02-data', 
                            '00-source',
                            'seki_access.shp'))


# DEM (30.7x30.6m resolution)
dem = 
  elevatr::get_elev_raster(locations = dx_aoi %>% st_buffer(dist = 100), 
                           z = 11,
                           clip = 'bbox', 
                           prj = st_crs(dx_aoi)$wkt)


mmi = 
  raster::raster(here::here('02-data',
                            '00-source',
                            'edart',
                            'mmi_sum5_sc408ss_2016.bsq')) %>%
  raster::crop(dx_aoi %>%
                 sf::st_buffer(100) %>%
                 sf::st_transform(crs(.))) %>%
  raster::mask(dx_aoi %>%
                 sf::st_buffer(100) %>%
                 sf::st_transform(crs(.)))

# https://www.fs.usda.gov/detailfull/r5/landmanagement/gis/?cid=fsbdev3_048029&width=full
calveg = 
  
  # load the calveg shapefile
  #st_read(here::here('02-data',
  #                   '00-from_external',
  #                   'r5',
  #                   'S_USA.EVMid_R05_SouthSierra.gdb'))
  st_read(here::here('02-data',
                     '02-intermediate',
                     'calveg_clipped.shp'))

#### site description for manuscript ###########################################

plot(calveg)

names(calveg)

calveg %>%
  st_intersection(dx_aoi %>% st_transform(st_crs(calveg))) %>%
  group_by(REGIONAL_D) %>%
  summarise() %>%
  ungroup() %>%
  mutate(area_m2 = as.numeric(st_area(.)),
         area_ha = area_m2 / 10000) %>%
  mutate(total_area_ha = sum(.$area_ha)) %>%
  mutate(prop_area = area_ha / total_area_ha) %>%
  arrange(desc(prop_area)) %>%
  dplyr::select(REGIONAL_D, area_ha, prop_area) %>%
  print(n = Inf)

#### define AOI ################################################################

# want to make inference about: mixed conifer - pine vegetation in the DX 
# crystal cave AOI
# calveg cropped / clipped to the AOI, reprojected, filtered to only the mixed 
# pine regional dominance type, and transformed to the crs of the edart MMI
aoi = 
  calveg %>%
  st_crop(dx_aoi %>%
            sf::st_buffer(100) %>%
            sf::st_transform(st_crs(calveg))) %>%
  st_intersection(dx_aoi %>%
                    sf::st_transform(st_crs(calveg))) %>%
  dplyr::select(REGIONAL_DOMINANCE_TYPE = REGIONAL_D) %>%
  filter(REGIONAL_DOMINANCE_TYPE=='MP') %>%
  group_by(REGIONAL_DOMINANCE_TYPE) %>%
  summarise() %>%
  ungroup() %>%
  sf::st_transform(crs = crs(mmi))

tmap_mode('view')
tm_shape(dx_aoi)+tm_borders('black')+tm_shape(aoi)+tm_polygons('green')

# make a raster version for masking the MMI layer
aoi_raster = 
  rasterize(aoi, mmi)


tm_shape(dx_aoi)+tm_borders('black')+tm_shape(aoi_raster)+tm_raster()+tm_shape(aoi)+tm_borders('red')


#### build drought mortality strata ############################################


plot(mmi)


mmi[mmi > 100] = NA

# Create the elevation mask to get rid of that island within the AOI


mort_smoothed =
  raster::focal(mmi,
                w = matrix(1, nrow = 3, ncol = 3),
                fun = mean) %>%
  raster::mask(aoi_raster)

tm_shape(dx_aoi)+tm_borders('black')+tm_shape(mort_smoothed)+tm_raster(palette = 'viridis')+tm_shape(aoi)+tm_borders('red')


ggplot()+
  geom_raster(data = as.data.frame(mort_smoothed, xy = TRUE) %>%
                filter(!is.na(layer)),
              aes(x = x, y = y, fill = layer))+
  theme_bw()+
  geom_sf(data = aoi, fill = NA, color = 'black')+
  scale_fill_viridis_c(option = 'B')


quantile(values(mort_smoothed), na.rm = TRUE)

low_mort_smoothed = 
  mort_smoothed < as.numeric(quantile(values(mort_smoothed), 0.25, na.rm = TRUE))

low_mort_smoothed[low_mort_smoothed == 0] = NA

plot(low_mort_smoothed)



high_mort_smoothed = 
  mort_smoothed >= as.numeric(quantile(values(mort_smoothed), 0.75, na.rm = TRUE))

high_mort_smoothed[high_mort_smoothed == 0] = NA

plot(high_mort_smoothed)

mid_mort_smoothed = 
  mort_smoothed >= as.numeric(quantile(values(mort_smoothed), 0.25, na.rm = TRUE)) &
  mort_smoothed < as.numeric(quantile(values(mort_smoothed), 0.75, na.rm = TRUE))

mid_mort_smoothed[mid_mort_smoothed == 0] = NA

mid_mort_smoothed.sf = 
  rasterToPolygons(mid_mort_smoothed, dissolve = TRUE) %>%
  sf::st_as_sf() %>%
  #sf::st_buffer(-30) %>%
  #sf::st_make_valid() %>%
  #sf::st_buffer(0) %>%
  mutate(mort = 'mid')

high_mort_smoothed.sf = 
  rasterToPolygons(high_mort_smoothed, dissolve = TRUE)%>%
  sf::st_as_sf() %>%
  #sf::st_buffer(-30) %>%
  #sf::st_make_valid() %>%
  #sf::st_buffer(0) %>%
  mutate(mort = 'high')

low_mort_smoothed.sf = 
  rasterToPolygons(low_mort_smoothed, dissolve = TRUE) %>%
  sf::st_as_sf() %>%
  #sf::st_buffer(-30) %>%
  #sf::st_make_valid() %>%
  #sf::st_buffer(0) %>%
  mutate(mort = 'low')

mort_categories = 
  rbind(mid_mort_smoothed.sf,
        high_mort_smoothed.sf,
        low_mort_smoothed.sf) %>%
  mutate(mort = factor(mort,
                       levels = c('high', 'mid', 'low')))

saveRDS(mort_categories,
        here::here('02-data',
                   '02-intermediate',
                   'intermediate_plot_location_files',
                   'mort_categories.rds'))

mort_categories_buffered = 
  mort_categories %>%
  sf::st_buffer(-30) %>%
  sf::st_make_valid() %>%
  sf::st_buffer(0)

ggplot()+
  geom_sf(data = mort_categories, aes(fill = mort))+
  geom_sf(data = aoi, fill = NA, color = 'black')+
  theme_bw()+
  scale_fill_brewer(palette = 'RdYlGn')

ggplot()+
  geom_sf(data = mort_categories_buffered, aes(fill = mort))+
  geom_sf(data = aoi, fill = NA, color = 'black')+
  theme_bw()+
  scale_fill_brewer(palette = 'RdYlGn')

#### build masks ###############################################################

# building an "accessible" mask (distance to road / trail > 50 and < 600, slope < 25 degrees)


# slope raster built form the DEM
slope = 
  raster::terrain(dem, opt = 'slope', unit = 'degrees')

# mask - polygons where slope is less than 25 degrees
slope25 = 
  rasterToPolygons(slope > 25, dissolve = TRUE) %>%
  sf::st_as_sf() %>%
  group_by(layer) %>%
  summarise() %>%
  filter(layer == 0)

# distance to road
access_buffer = 
  access %>%
  sf::st_buffer(600) %>%
  summarise() %>%
  sf::st_difference(st_combine(access %>%
                                 sf::st_buffer(100) %>%
                                 summarise()))


access_mask = 
  access_buffer %>%
  st_intersection(slope25 %>%
                    st_transform(crs = st_crs(access_buffer))) %>%
  st_intersection(aoi) %>%
  summarise()

saveRDS(access_mask,
        here::here('02-data',
                   '02-intermediate',
                   'intermediate_plot_location_files',
                   'access_mask.rds'))

tm_shape(access_buffer)+tm_polygons()
names(dem) = 'layer'
ggplot()+
  geom_raster(data = as.data.frame(dem, xy = TRUE),
              aes(x = x, y = y, fill = layer))+
  scale_fill_viridis_c()+
  geom_sf(data = aoi, fill = NA, color = 'black')+
  geom_sf(data = access_mask, fill = 'green', color = NA, alpha = 0.5)+
  theme_bw()+
  labs(fill = 'Elev (m)')


tm_shape(dx_aoi)+tm_borders('black')+tm_shape(aoi)+tm_borders('blue')+tm_shape(access_mask)+tm_polygons('green')

sample_area = 
  mort_categories_buffered %>%
  st_intersection(access_mask)

# as of 06-04-21, the crew has sampled two plots (23 and 2) and uploaded the true updated 
# coordinates for them
existing_plots = 
  st_read(here::here('02-data',
                     '00-from_external',
                     'trimble',
                     'seki_060421.shp')) %>%
  st_set_crs(crs(mmi)) %>%
  filter(is.element(plot_id, c(2, 23)))

# plot 23 actually falls slightly outside (~15m) the core area of the low-mortality patch 
# it's within, but it's still within a low mortality zone and for one plot we've 
# already collected I'm willing to slightly fudge the "30m inward buffer" rule.
tm_shape(dx_aoi)+tm_borders('black')+tm_shape(aoi)+tm_borders('blue')+tm_shape(sample_area)+tm_polygons('mort')+tm_shape(existing_plots)+tm_dots()+tm_scale_bar()


#### comparing mask versions ###################################################

# fasterize handles holes more consistently
access_mask_raster = 
  fasterize::fasterize(access_mask, dem) %>%
  projectRaster(dem)



aoi_raster = 
  rasterize(aoi, dem) %>%
  projectRaster(dem)

aoi_raster

accessible_elev = 
  mask(dem, access_mask_raster)
plot(accessible_elev)
plot(access_mask, add = TRUE, color = NA)


aoi_elev = 
  mask(dem, aoi_raster)

plot(aoi_elev)
plot(aoi, add = TRUE, color = NA)


accessible_slope = 
  mask(slope, access_mask_raster)

aoi_slope = 
  mask(slope, aoi_raster)

southwestness = 
  terrain(dem, opt = 'aspect', unit = 'degrees')
southwestness = 1-(abs(southwestness - 225)/225)

accessible_sw = mask(southwestness, access_mask_raster)

aoi_sw = mask(southwestness, aoi_raster)


terrain_compare = 
  as.data.frame(aoi_elev) %>%
  filter(!is.na(layer)) %>%
  select(value = layer) %>%
  mutate(param = 'elev',
         mask = 'aoi') %>%
  
  bind_rows(as.data.frame(accessible_elev) %>%
              filter(!is.na(layer)) %>%
              select(value = layer) %>%
              mutate(param = 'elev',
                     mask = 'accessible')) %>%
  
  bind_rows(as.data.frame(aoi_slope) %>%
              filter(!is.na(slope)) %>%
              select(value = slope) %>%
              mutate(param = 'slope',
                     mask = 'aoi')) %>%
  
  bind_rows(as.data.frame(accessible_slope) %>%
              filter(!is.na(slope)) %>%
              select(value = slope) %>%
              mutate(param = 'slope',
                     mask = 'accessible')) %>%

  bind_rows(as.data.frame(aoi_sw) %>%
              filter(!is.na(layer)) %>%
              select(value = layer) %>%
              mutate(param = 'sw',
                     mask = 'aoi')) %>%
  
  bind_rows(as.data.frame(accessible_sw) %>%
              filter(!is.na(layer)) %>%
              select(value = layer) %>%
              mutate(param = 'sw',
                     mask = 'accessible')) 

terrain_compare_plot = 
  ggplot(data = terrain_compare,
       aes(x = mask, y = value)) +
  geom_boxplot()+
  facet_wrap(~param, scales = 'free_y')

terrain_compare_plot

ggsave(terrain_compare_plot,
       filename = 
         here::here('04-communication',
                    'figures',
                    'manuscript',
                    'aoi_access_comparison.png'),
       height = 4, width = 6.5, units = 'in')

tm_shape(slope)+
  tm_raster(palette = 'viridis', style = 'cont')+
  tm_shape(slope25)+
  tm_borders('white')+
  tm_shape(access_mask)+
  tm_polygons()+
  tm_shape(access_mask_raster)+
  tm_raster()+
  tm_shape(accessible_slope)+
  tm_raster(style = 'cont', palette = 'inferno')

#### draw plots ################################################################

set.seed(112188)
plot_locations = 
  sf::st_sample(sample_area, size = 1000) %>%
  sf::st_as_sf() %>%
  sf::st_join(sample_area) %>%
  select(-layer) %>%
  group_by(mort) %>%
  sample_n(size = 6) %>%
  ungroup() %>%
  rowid_to_column('plot_id')


tm_shape(dx_aoi)+tm_borders('black')+
  tm_shape(access_mask)+tm_borders('green')+
  tm_shape(sample_area)+
  tm_polygons('mort')+
  tm_shape(plot_locations)+
  tm_dots('red')+
  tm_shape(existing_plots)+
  tm_dots('blue')+
  tm_scale_bar()

tm_shape(dx_aoi)+tm_borders('black')+
  tm_shape(mort_smoothed)+
  tm_raster(palette = 'inferno')+
  #tm_shape(access_mask)+tm_borders('green')+
  tm_shape(plot_locations)+
  tm_dots('green')+
  tm_shape(existing_plots)+
  tm_dots('blue')+
  tm_scale_bar()
#### combine with old plots ####################################################

# we already have one low and one high mort plot, so drop plots 6 and 18 
plot_locations = 
  plot_locations %>%
  filter(!is.element(plot_id, c(6, 18)))

# plot_id 2 is already taken in the existing data, so relabel it plot 6 
plot_locations[plot_locations$plot_id==2,'plot_id'] = 6

# add back in the existing plots (2 and 23)
plot_locations = 
  plot_locations %>%
  bind_rows(st_zm(existing_plots, what = 'ZM'))

tm_shape(dx_aoi)+tm_borders('black')+
  tm_shape(access_mask)+tm_borders('green')+
  tm_shape(sample_area)+
  tm_polygons('mort')+
  tm_shape(plot_locations)+
  tm_dots('red')+
  tm_scale_bar()


plot_locations.df = 
  plot_locations %>%
  mutate(x_coord = sf::st_coordinates(.)[,'X'],
         y_coord = sf::st_coordinates(.)[,'Y']) %>%
  as.data.frame() %>%
  select(-geometry)


#### write results #############################################################

sf::st_write(plot_locations,
             here::here('02-data', '02-intermediate', 'plot_locations_redo.shp'),
             delete_dsn = TRUE)

write.csv(plot_locations.df,
          here::here('02-data', '02-intermediate', 'plot_locations_redo.csv'))

#### map for paper #############################################################

library(ggspatial)
library(cowplot)

plot_locations_true = 
  st_read(here::here('02-data', '02-intermediate', 'plot_locations_true.shp'))


overview_map = 
  ggplot(data = 
         spData::us_states %>%
           filter(NAME=='California'))+
  geom_sf(fill = NA, lwd = 1)+
  geom_sf(data = 
               aoi %>%
               st_centroid,
             color = 'red',
          size = 5,
          pch = 8)+
  coord_sf(crs = "EPSG:26910")+
  theme_minimal()+
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = 'white'),
        plot.margin = unit(c(0, 0, 0, 0), 'cm'))

aoi_map = 
  ggplot(data = aoi)+
  geom_sf(data = mort_categories %>%
            mutate(`Mortality Class` = 
                     ifelse(mort=='low',
                            'Low (<7%)',
                            ifelse(mort == 'mid',
                                   'Medium (7-18%)',
                                   'High (>18%)'))) %>%
            mutate(`Mortality Class` = 
                     factor(`Mortality Class`,
                            levels = c('High (>18%)',
                                       'Medium (7-18%)',
                                       'Low (<7%)'))),
          aes(fill = `Mortality Class`),
          color = NA,
          alpha = 0.5)+
  #geom_sf(fill = NA, lwd = 1)+
  geom_sf(
    data = plot_locations_true,
    show.legend = TRUE
  )+
  coord_sf(crs = 'EPSG:26910')+
  theme_minimal()+
  annotation_scale()+
  scale_fill_viridis_d(begin = 0.05, 
                       end = 0.85, 
                       option = 'C',
                       direction = -1,
                       guide = guide_legend(override.aes = list(shape = NA)))+
  scale_color_manual(values = c('Inventory plots' = 'black'),
                     name = NULL,
                     guide = guide_legend(override.aes = list(shape = 16)))+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank())+
  annotation_north_arrow(pad_y = unit(1, 'cm'))

aoi_map

map_figure = 
  ggdraw()+
  draw_plot(aoi_map)+
  draw_plot(overview_map,
            x = 0.75, y = 0.05, width = 0.25, height = 0.25)

map_figure

ggsave(map_figure,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'overview_map.png'),
       height = 4.5, width = 6.5, units = 'in')

