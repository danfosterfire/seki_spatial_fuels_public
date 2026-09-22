library(here)
library(raster)
library(sf)
library(tidyverse)
library(elevatr)
library(ggplot2)
library(tmap)


#### load files ################################################################

# dx plot locations
dx_plots = 
  read_csv(here::here('02-data',
                      '00-source',
                      'seki_dx_plots.csv')) %>%
  st_as_sf(coords = c('Longitude.DD','Latitude.DD'),
           crs = 4326)

head(dx_plots)
# crystal cave DX AOI, defined using an elevation band around crystal cave road
# and environs
dx_aoi = st_read(here::here('02-data',
                            '00-source',
                            'seki_dx',
                            'dxstudyarea.shp')) 

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

#### clip the AOI to mixed pine regional dominance type ########################

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
tm_shape(dx_aoi)+tm_borders('black')+tm_shape(aoi)+tm_polygons('green')+
  tm_shape(dx_plots)+tm_dots()



#### bin the mortality index values into quantiles and assemble shapefile ######

# make a raster version for masking the MMI layer
aoi_raster = 
  rasterize(aoi, mmi)

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

# combine the separate layers; this is the shape we'll extract the classes from
mort_categories = 
  rbind(mid_mort_smoothed.sf,
        high_mort_smoothed.sf,
        low_mort_smoothed.sf) %>%
  mutate(mort = factor(mort,
                       levels = c('high', 'mid', 'low')))


tm_shape(mort_categories)+tm_polygons('mort')

#### prepare products for john #################################################

# total area by mort class (mixed pine dominance type only)
mort_categories %>%
  mutate(area_ha = as.numeric(st_area(.))/10000) %>%
  st_write(here::here('02-data',
                      '05-for_analysis',
                      'aoi_mort_categories.shp'))

mort_categories %>%
  mutate(area_ha = as.numeric(st_area(.))/10000) %>%
  left_join(as_data_frame(.) %>%
              group_by(layer) %>%
              summarise(total_area_ha = sum(area_ha)) %>%
              ungroup()) %>%
  mutate(proportion_of_total = area_ha/total_area_ha) %>%
  as_data_frame() %>%
  select(-geometry)

dx_plots %>%
  st_transform(crs = st_crs(mort_categories)) %>%
  st_intersection(mort_categories) %>%
  mutate(mort = as.character(mort)) %>%
  mutate(mort = ifelse(is.na(mort),'NOT IN MIXED PINE AOI',mort)) %>%
  st_transform(crs = 4326) %>%
  mutate(lat = st_coordinates(.)[,'Y'],
         lon = st_coordinates(.)[,'X']) %>%
  as_data_frame() %>%
  select(-geometry) %>%
  write_csv(here::here('02-data',
                       '05-for_analysis',
                       'dx_plots_mort_classes.csv'))

