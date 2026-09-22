library(here)
library(raster)
library(sf)
library(tidyverse)
library(elevatr)
library(ggplot2)
set.seed(110819)

#### aoi and masks #############################################################
aoi = st_read(here::here('02-data',
                         '00-source',
                         'seki_dx',
                         'dxstudyarea.shp')) 
st_crs(aoi)
dem = 
  elevatr::get_elev_raster(locations = aoi %>% st_buffer(dist = 100), 
                           z = 11,
                           clip = 'bbox', 
                           prj = 'EPSG:26911')

test = raster::extract(dem, aoi, df = TRUE)

head(test)

summary(test) # min 1517 max 1844

worldclim = getData(name = 'worldclim',
                    download = TRUE,
                    path = here::here('02-data', '00-source', 'worldclim'),
                    res = 0.5,
                    lat = 36.76,
                    lon = -118.80,
                    var = 'bio')
plot(worldclim)

aoi_clim = raster::extract(worldclim, st_transform(aoi, 'EPSG:4326'), df = TRUE)

summary(aoi_clim)

head(aoi_clim)


names(dem) = 'layer'

slope = 
  raster::terrain(dem, opt = 'slope', unit = 'degrees')

slope50 = 
  rasterToPolygons(slope > 25, dissolve = TRUE) %>%
  sf::st_as_sf() %>%
  group_by(layer) %>%
  summarise() %>%
  filter(layer == 0)





ggplot()+
  geom_raster(data = as.data.frame(dem, xy = TRUE),
              aes(x = x, y = y, fill = layer))+
  scale_fill_viridis_c()+
  geom_sf(data = aoi, fill = NA, color = 'black')+
  geom_sf(data = slope50, fill = 'green', color = NA, alpha = 0.5)+
  theme_bw()+
  labs(fill = 'Elev (m)')



#### calveg ####################################################################


calveg = 
  
  # load the calveg shapefile
  st_read(here::here('02-data',
                     '00-source',
                     'calveg',
                     'S_USA.EVMid_R05_SouthSierra.gdb'))

# this is throwing errors as of 8/16/2021; i did a rough crop in arcGIS instead
#calveg_cropped = 
#  calveg %>%
  
#  # rough crop the shapefile by the buffered aoi
#  st_crop(aoi %>%
#            sf::st_buffer(100) %>%
#            sf::st_transform(st_crs(calveg)))

calveg_cropped = 
  st_read(here::here('02-data',
                     '02-intermediate',
                     'calveg_clipped.shp'))

calveg_clipped = 
  calveg_cropped %>%
  
  # select the regional dominance type data
  dplyr::select(REGIONAL_DOMINANCE_TYPE = REGIONAL_D) %>%
  
  # transform to the aoi CRS and clip by the aoi and the slope mask
  st_transform(st_crs(aoi)) %>%
  st_intersection(aoi) %>%
  st_intersection(slope50) %>%
  
  # calculate the area of each polygon
  mutate(area_ha = as.numeric(st_area(.)/10000)) %>%
  
  # group and aggregate by the veg type
  group_by(REGIONAL_DOMINANCE_TYPE) %>%
  summarise(area_ha = sum(area_ha)) %>%
  ungroup() %>%
  
  # calculate the proportion of the total aarea occupied by each veg type
  mutate(total_area = sum(pull(., area_ha)),
         prop_area = area_ha / total_area) %>%
  arrange(desc(prop_area)) %>%
  # filter to only the mixed pine cover type (using CWHR==CON was too broad,
  # included some bare areas which were obvious on imagery)
  filter(is.element(REGIONAL_DOMINANCE_TYPE,
                    c('MP')))

calveg_lifeform = 
  calveg_cropped %>%
  
  # select the regional dominance type data
  dplyr::select(CWHR_LIFEFORM = CWHR_LIFEF) %>%
  
  # transform to the aoi CRS and clip by the aoi and the slope mask
  st_transform(st_crs(aoi)) %>%
  st_intersection(aoi) %>%
  st_intersection(slope50) %>%
  
  # calculate the area of each polygon
  mutate(area_ha = as.numeric(st_area(.)/10000)) %>%
  
  # group and aggregate by the veg type
  group_by(CWHR_LIFEFORM) %>%
  summarise(area_ha = sum(area_ha)) %>%
  ungroup() %>%
  
  # calculate the proportion of the total aarea occupied by each veg type
  mutate(total_area = sum(pull(., area_ha)),
         prop_area = area_ha / total_area) %>%
  arrange(desc(prop_area))

ggplot()+
  geom_sf(data = aoi, fill = NA, color = 'black')+
  geom_sf(data = calveg_clipped, aes(fill = REGIONAL_DOMINANCE_TYPE))+
  theme_bw()

ggplot()+
  geom_sf(data = aoi, fill = NA, color = 'black')+
  geom_sf(data = calveg_lifeform, aes(fill = CWHR_LIFEFORM))+
  theme_bw()

library(tmap)
tmap_mode('view')
tm_shape(calveg_clipped)+tm_polygons('REGIONAL_DOMINANCE_TYPE')

#### drought mortality #########################################################


mmi = 
  raster::raster(here::here('02-data',
                            '00-source',
                            'edart',
                            'mmi_sum5_sc408ss_2016.bsq')) %>%
  raster::crop(aoi %>%
                 sf::st_buffer(100) %>%
                 sf::st_transform(crs(.))) %>%
  raster::mask(aoi %>%
                 sf::st_transform(crs(.)))

# https://www.fs.usda.gov/detailfull/r5/landmanagement/gis/?cid=fsbdev3_048029&width=full
calveg_mask = 
#  calveg %>%
  calveg_cropped %>%
  st_crop(aoi %>%
            sf::st_buffer(100) %>%
            sf::st_transform(st_crs(calveg_cropped))) %>%
  dplyr::select(REGIONAL_DOMINANCE_TYPE = REGIONAL_D) %>%
  filter(REGIONAL_DOMINANCE_TYPE=='MP') %>%
  sf::st_transform(crs = crs(mmi))


plot(mmi)


mmi[mmi > 100] = NA

# Create the elevation mask to get rid of that island within the AOI


mort_smoothed =
  raster::focal(mmi,
                w = matrix(1, nrow = 3, ncol = 3),
                fun = mean) %>%
  raster::mask(slope50) %>%
  raster::mask(calveg_mask)


elev_mask = 
  raster::projectRaster(from = dem, to = mort_smoothed)

elev_mask = 
  (elev_mask < (6000*0.3048) & elev_mask > (5000*0.3048))

elev_mask[elev_mask==0] = NA

mort_smoothed = 
  raster::mask(mort_smoothed, elev_mask)

plot(elev_mask)
plot(mort_smoothed)

ggplot()+
  geom_raster(data = as.data.frame(mort_smoothed, xy = TRUE) %>%
                filter(!is.na(layer)),
              aes(x = x, y = y, fill = layer))+
  theme_bw()+
  geom_sf(data = aoi, fill = NA, color = 'black')+
  scale_fill_viridis_c(option = 'B')


quantile(values(mort_smoothed), na.rm = TRUE)
# 25%: 7.222 50%: 12.222    75%: 18.5556  100%: 54.555

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


mort_categories_buffered = 
  mort_categories %>%
  sf::st_buffer(-30) %>%
  sf::st_make_valid() %>%
  sf::st_buffer(0)


# write out some intermediate layers for making maps
st_write(mort_categories,
         here::here('02-data',
                    '02-intermediate',
                    'intermediate_plot_location_files',
                    'mort_categories.shp'))


st_write(calveg_clipped,
         here::here('02-data',
                    '02-intermediate',
                    'intermediate_plot_location_files',
                    'calveg_clipped.shp'),
         delete_dsn = TRUE)

saveRDS(calveg_clipped,
        here::here('02-data',
                   '02-intermediate',
                   'intermediate_plot_location_files',
                   'calveg_clipped.rds'))

saveRDS(calveg_mask,
        here::here('02-data',
                   '02-intermediate',
                   'intermediate_plot_location_files',
                   'calveg_mask.rds'))

saveRDS(slope50,
        here::here('02-data',
                   '02-intermediate',
                   'intermediate_plot_location_files',
                   'slope_mask.rds'))

saveRDS(mort_smoothed,
        here::here('02-data',
                   '02-intermediate',
                   'intermediate_plot_location_files',
                   'mort_smoothed.rds'))

saveRDS(mort_categories,
        here::here('02-data',
                   '02-intermediate',
                   'intermediate_plot_location_files',
                   'mort_categories.rds'))

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

library(tmap)
tmap_mode('view')

tm_shape(mort_smoothed)+
  tm_raster(palette = '-inferno')+
  tm_shape(low_mort_smoothed.sf)+tm_borders('green', lwd = 2)+
  tm_shape(high_mort_smoothed.sf)+tm_borders('red', lwd = 2)+
  tm_shape(mid_mort_smoothed.sf)+tm_borders('yellow', lwd = 2)

plot_locations = 
  sf::st_sample(mort_categories_buffered, size = 1000) %>%
  sf::st_as_sf() %>%
  sf::st_join(mort_categories_buffered) %>%
  select(-layer) %>%
  group_by(mort) %>%
  sample_n(size = 10) %>%
  ungroup() %>%
  rowid_to_column('plot_id')

# some plots are too close (within ~100m) of a road, ditch these
# everything remaining is AT LEAST 50m from a road
# did this manually because finding a road shapefile for a national park
# was too much hassle
plot_locations_final = 
  plot_locations %>%
  filter(!is.element(plot_id,
                     c(22, 30, 26, 17, 15, 25))) %>%
  
  # then select the 6 (5 plots + 1 alternate) first ID values in each 
  # drought mortality class
  filter(is.element(plot_id,
                    c(1, 2, 3, 4, 5, 6,
                      11, 12, 13, 14, 16, 18,
                      21, 23, 24, 27, 28, 29)))


plot_locations.df = 
  plot_locations_final %>%
  mutate(x_coord = sf::st_coordinates(.)[,'X'],
         y_coord = sf::st_coordinates(.)[,'Y']) %>%
  as.data.frame() %>%
  select(-geometry)

tm_shape(aoi)+
  tm_borders('black')+
  tm_shape(plot_locations_final)+
  tm_dots('mort')+
  tm_scale_bar()

#### write results #############################################################

sf::st_write(plot_locations_final,
             here::here('02-data', '02-intermediate', 'plot_locations.shp'),
             delete_dsn = TRUE)

write.csv(plot_locations.df,
          here::here('02-data', '02-intermediate', 'plot_locations.csv'))
