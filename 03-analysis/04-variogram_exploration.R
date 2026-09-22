

#### setup #####################################################################

# load packages
library(here)
library(tidyverse)
library(gstat)




# load data
metadata = readRDS(here::here('02-data', '03-clean', 'seki_metadata.rds'))

plot_locations.tibble = 
  readRDS(here::here('02-data', '02-intermediate', 'plot_locations_true.rds'))

# geolocate metadata
metadata = 
  metadata %>%
  left_join(plot_locations.tibble %>%
              select(mort, plot_id, x_center = x_true, y_center = y_true),
              by = c('plot_id'))


#### prepare fuel components data ##############################################

litterduff = 
  readRDS(here::here('02-data', '03-clean', 'seki_litterduff.rds')) %>%
  left_join(metadata %>%
              select(plot_id, transect_id, mort, x_center, y_center)) %>%
  filter(!is.na(litter_cm) & !is.na(duff_cm)) %>%
  mutate(x_rel = (location_m * sin(az*(pi/180))),
         y_rel = (location_m * cos(az*(pi/180))),
         x_abs = x_rel + x_center,
         y_abs = y_rel + y_center)

head(litterduff)

fwd = 
  
  readRDS(here::here('02-data', '03-clean', 'seki_fwdtallies.rds'))  %>%
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
              select(plot_id, transect_id, mort, x_center, y_center)) %>%
  mutate(x_abs = x_rel + x_center,
         y_abs = y_rel + y_center) %>%
  filter(!is.na(count))


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
  left_join(metadata %>%
              select(plot_id, mort, x_center, y_center),
            by = c('plot_id' = 'plot_id')) %>%  
  # want to treat the N and S transects as a single transect, same with the E/W
  mutate(direction = ifelse(is.element(az, c(0, 180)),'NS', 'EW'),
         transect_id = paste0(plot_id,'-',direction)) %>%
  mutate(plot_id.i = as.integer(factor(plot_id)),
         transect_id.i = as.integer(factor(transect_id)),
         group_id = as.integer(factor(mort, levels = c('low', 'mid', 'high')))) %>%
  mutate(
    x_rel = round(location_m*sin(az*pi/180),2),
    y_rel = round(location_m*cos(az*pi/180),2),
    x_abs = x_rel + x_center,
    y_abs = y_rel + y_center
  ) %>%
  filter(!is.na(count))

head(cwd)

# making quadrature points for understory veg presence or absence or 
# height
veg = 
  metadata %>%
  select(mort, plot_id, transect_id, az, x_center, y_center) %>%
  left_join(readRDS(here::here('02-data',
                               '03-clean',
                               'seki_veg.rds'))) %>%
  expand(nesting(mort, plot_id, transect_id, az, start_m, end_m, spp, status, height_m, comments,
                 x_center, y_center),
         sample_point = seq(from = 0.25, to = 14.75, by = 0.5)) %>%
  # veg present is NA at all points if there was no veg on the transect
  mutate(veg_present = start_m <= sample_point & end_m >= sample_point) %>%
  
  # convert NAs to FALSE and get height
  mutate(veg_present = ifelse(is.na(veg_present), FALSE, veg_present),
         veg_height = ifelse(veg_present, height_m, NA)) %>%
  group_by(mort, plot_id, transect_id, az, sample_point, x_center, y_center) %>%
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
                                    levels = c('low', 'mid', 'high'))),
    x_abs = x_rel + x_center,
    y_abs = y_rel + y_center) %>%
  
  # get height as number of layers
  mutate(height_N = ceiling(veg_height / 0.25)) %>%
  filter(!is.na(veg_height))

head(veg)


# make a quadrature of 1m2 cells for each plot
trees_quadrature = 
  expand.grid(x_rel = seq(-14.5, 14.5, 1),
              y_rel = seq(-14.5, 14.5, 1)) %>%
  filter((x_rel >= -4.5 & x_rel <= 4.5) | (y_rel >= -4.5 & y_rel <= 4.5)) %>%
  expand(nesting(x_rel, y_rel),
         plot_locations.tibble %>%
           group_by(plot_id, mort, x_true, y_true) %>%
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
                           ifelse(mort == 'mid', 2, 3))) %>%
  mutate(x_abs = x_true + x_rel,
         y_abs = y_true + y_rel)


# make a quadrature of 1m2 cells for each plot
saplings_quadrature = 
  expand.grid(x_rel = seq(-14.5, 14.5, 1),
              y_rel = seq(-14.5, 14.5, 1)) %>%
  filter((x_rel >= -0.5 & x_rel <= 0.5) | (y_rel >= -0.5 & y_rel <= 0.5)) %>%
  expand(nesting(x_rel, y_rel),
         plot_locations.tibble %>%
           group_by(plot_id, mort, x_true, y_true) %>%
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
                           ifelse(mort == 'mid', 2, 3))) %>%
  mutate(x_abs = x_true + x_rel,
         y_abs = y_true + y_rel)

#### make plot diagram #########################################################

fwd_tallies_transects = 
  data.frame(x = c(-30, -20, -10, -8, -6, -5, -4, -3, 2, 3, 4, 5, 7, 9, 19, 29,
                   rep(0, 16),
                   rep(-0.5, 16),
                   -29.5, -19.5, -9.5, -7.5, -5.5, -4.5, -3.5, -2.5, 2.5, 3.5, 4.5, 5.5, 7.5, 9.5, 19.5, 29.5),
             xend = c(-29, -19, -9, -7, -5, -4, -3, -2, 3, 4, 5, 6, 8, 10, 20, 30,
                      rep(0, 16),
                      rep(0.5, 16),
                      -29.5, -19.5, -9.5, -7.5, -5.5, -4.5, -3.5, -2.5, 2.5, 3.5, 4.5, 5.5, 7.5, 9.5, 19.5, 29.5),
             y = c(rep(0, 16), 
                   -30, -20, -10, -8, -6, -5, -4, -3, 2, 3, 4, 5, 7, 9, 19, 29,
                   -29.5, -19.5, -9.5, -7.5, -5.5, -4.5, -3.5, -2.5, 2.5, 3.5, 4.5, 5.5, 7.5, 9.5, 19.5, 29.5,
                   rep(-0.5, 16)),
             yend = c(rep(0, 16),
                      -29, -19, -9, -7, -5, -4, -3, -2, 3, 4, 5, 6, 8, 10, 20, 30,
                      -29.5, -19.5, -9.5, -7.5, -5.5, -4.5, -3.5, -2.5, 2.5, 3.5, 4.5, 5.5, 7.5, 9.5, 19.5, 29.5,
                      rep(0.5, 16)),
             component = 'fwd')

cwd_tallies_transects = 
  data.frame(x = c(-30, 0, 0, 0),
             xend = c(0, 30, 0, 0),
             y = c(0, 0, -30, 0),
             yend = c(0, 0, 0, 30),
             component = 'main')

trees_bounds = 
  data.frame(x = c(-15, -15, -15, -5, -5, -5, -5, 5, 5, 5, 5, 15),
             xend = c(-15, -5, -5, -5, -5, 5, 5, 5, 5, 15, 15, 15),
             y = c(-5, 5, -5, 5, -5, 15, -15, 15, -15, 5, -5, 5),
             yend = c(5, 5, -5, 15, -15, 15, -15, 5, -5, 5, -5, -5),
             component = 'trees')

saplings_bounds = 
  data.frame(x = c(-15, -15, -15, -1, -1, -1, -1, 1, 1, 1, 1, 15),
             xend = c(-15, -1, -1, -1, -1, 1, 1, 1, 1, 15, 15, 15),
             y = c(-1, 1, -1, 1, -1, 15, -15, 15, -15, 1, -1, 1),
             yend = c(1, 1, -1, 15, -15, 15, -15, 1, -1, 1, -1, -1),
             component = 'saplings')

plot_diagram = 
  ggplot()+  
  scale_x_continuous(limits = c(-31, 31))+
  scale_y_continuous(limits = c(-31, 31))+
  geom_segment(data = cwd_tallies_transects,
               aes(x = x, y = y, xend = xend, yend = yend, color = component),
               lwd = 1)+
  geom_segment(data = fwd_tallies_transects,
               aes(color = component,
                   x = x, y = y, xend = xend, yend = yend),
               lwd = 1)+
  geom_segment(data = trees_bounds,
               aes(x = x, y = y, xend = xend, yend = yend, color = component),
               lwd = 1)+
  geom_segment(data = saplings_bounds,
               aes(x = x, y = y, xend = xend, yend = yend, color = component),
               lwd = 1)+
  geom_point(data = data.frame(x = 0, y = 0, component = 'center'),
             aes(x = x, y = y, pch = 'Plot center'), color = 'black', size = 4)+
  coord_fixed()+
  theme_minimal()+
  scale_color_brewer(labels = c(fwd = '1-meter FWD subtransects',
                                main = '30-meter main transects',
                                trees = expression(paste('500 ', m^2, ' trees plot')),
                                saplings = expression(paste('116 ', m^2, ' saplings plot'))),
                     palette = 'Set1')+
  scale_shape_manual(values = c('Plot center' = 4))+
  theme(legend.title = element_blank(),
        text = element_text(size = 16))+
  labs(x = 'Easting (m)',
       y = 'Northing (m)')


plot_diagram

ggsave(plot_diagram,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'plot_diagram.png'),
       height = 9, width = 13, units = 'in')


htrees_sample = 
  ggplot(data = 
         seki.trees %>%
         filter(plot_id=='14') %>%
         left_join(seki.metadata %>%
                     dplyr::select(plot_id, x_center, y_center)) %>%
         mutate(x_relative = x_coord-x_center, y_relative = y_coord-y_center),
       aes(x = x_relative, y = y_relative, size = dbh_cm, color = height_m))+
  geom_point(alpha = 0.75)+
  geom_segment(x = -15, y = -5, xend = -15, yend = 5, col = 'black', lwd = 1)+
  geom_segment(x = 15, y = -5, xend = 15, yend = 5, col = 'black', lwd = 1)+
  geom_segment(x = -5, y = 15, xend = 5, yend = 15, col = 'black', lwd = 1)+
  geom_segment(x = -5, y = -15, xend = 5, yend = -15, col = 'black', lwd = 1)+
  geom_segment(x = -15, y = -5, xend = -5, yend = -5, col = 'black', lwd = 1)+
  geom_segment(x = -15, y = 5, xend = -5, yend = 5, col = 'black', lwd = 1)+
  geom_segment(x = -5, y = -15, xend = -5, yend = -5, col = 'black', lwd = 1)+
  geom_segment(x = -5, y = 15, xend = -5, yend = 5, col = 'black', lwd = 1)+
  geom_segment(x = 5, y= 15, xend = 5, yend = 5, col = 'black', lwd = 1)+
  geom_segment(x = 5, y = -15, xend = 5, yend = -5, col = 'black', lwd = 1)+
  geom_segment(x = 15, y = 5, xend = 5, yend = 5, col = 'black', lwd = 1)+
  geom_segment(x = 15, y = -5, xend = 5, yend = -5, col = 'black', lwd = 1)+
  geom_segment(x = -15, y = -1, xend = -15, yend = 1, col = 'black', lwd = 1, lty = 2)+
  geom_segment(x = 15, y = -1, xend = 15, yend = 1, col = 'black', lwd = 1, lty = 2)+
  geom_segment(x = -1, y = 15, xend = 1, yend = 15, col = 'black', lwd = 1, lty = 2)+
  geom_segment(x = -1, y = -15, xend = 1, yend = -15, col = 'black', lwd = 1, lty = 2)+
  geom_segment(x = -15, y = -1, xend = -1, yend = -1, col = 'black', lwd = 1, lty = 2)+
  geom_segment(x = -15, y = 1, xend = -1, yend = 1, col = 'black', lwd = 1, lty = 2)+
  geom_segment(x = -1, y = -15, xend = -1, yend = -1, col = 'black', lwd = 1, lty = 2)+
  geom_segment(x = -1, y = 15, xend = -1, yend = 1, col = 'black', lwd = 1, lty = 2)+
  geom_segment(x = 1, y= 15, xend = 1, yend = 1, col = 'black', lwd = 1, lty = 2)+
  geom_segment(x = 1, y = -15, xend = 1, yend = -1, col = 'black', lwd = 1, lty = 2)+
  geom_segment(x = 15, y = 1, xend = 1, yend = 1, col = 'black', lwd = 1, lty = 2)+
  geom_segment(x = 15, y = -1, xend = 1, yend = -1, col = 'black', lwd = 1, lty = 2)+
  scale_x_continuous(limits = c(-16, 16))+
  scale_y_continuous(limits = c(-16, 16))+
  coord_fixed()+
  scale_color_viridis_c()+
  theme_bw()+
  labs(x = 'Easting (m)', y = 'Northing (m)', color = 'Height (m)', 
       size = 'DBH (cm)', title = 'Trees and saplings')+
  ppt_theme+
  theme(legend.box = 'vertical',
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.spacing = unit(0.1, units = 'inches'),
        legend.margin = margin(unit(0.1, units = 'inches')))

trees_sample

#### combine fuel components ###################################################

head(litterduff)
head(fwd)

# start with duff
g.low = 
  gstat(id = 'duff',
        formula = duff_cm ~ 1,
        locations = ~x_abs+y_abs,
        data = litterduff %>% filter(mort == 'low'))

# add litter 
g.low = 
  gstat(g = g.low,
        id = 'litter',
        formula = litter_cm ~ 1,
        locations = ~x_abs+y_abs,
        data = litterduff %>% filter(mort == 'low'))

# add 1h, 10h, 100h
g.low = 
  gstat(g = g.low,
        id = 'fwd1h',
        formula = count ~ 1,
        locations = ~x_abs+y_abs,
        data = fwd %>% filter(mort == 'low' & timelag_class=='1h'))
g.low = 
  gstat(g = g.low,
        id = 'fwd10h',
        formula = count ~ 1,
        locations = ~x_abs+y_abs,
        data = fwd %>% filter(mort == 'low' & timelag_class=='10h'))
g.low = 
  gstat(g = g.low,
        id = 'fwd100h',
        formula = count ~ 1,
        locations = ~x_abs+y_abs,
        data = fwd %>% filter(mort == 'low' & timelag_class=='100h'))

# add in 1000h
g.low = 
  gstat(g = g.low,
        id = 'cwd',
        formula = count ~ 1,
        locations = ~x_abs+y_abs,
        data = cwd %>% filter(mort == 'low'))

# add in veg heights
g.low = 
  gstat(g = g.low,
        id = 'vegpres',
        formula = veg_present ~ 1,
        locations = ~x_abs+y_abs,
        data = veg %>% filter(mort == 'low'))

g.low = 
  gstat(g = g.low,
        id = 'trees',
        formula = count ~ 1,
        locations = ~x_abs+y_abs,
        data = trees_quadrature %>% filter(mort=='low'))

g.low = 
  gstat(g = g.low,
        id = 'saplings',
        formula = count ~ 1,
        locations = ~x_abs+y_abs,
        data = saplings_quadrature %>% filter(mort=='low'))



# start with duff
g.mid = 
  gstat(id = 'duff',
        formula = duff_cm ~ 1,
        locations = ~x_abs+y_abs,
        data = litterduff %>% filter(mort == 'mid'))

# add litter 
g.mid = 
  gstat(g = g.mid,
        id = 'litter',
        formula = litter_cm ~ 1,
        locations = ~x_abs+y_abs,
        data = litterduff %>% filter(mort == 'mid'))

# add 1h, 10h, 100h
g.mid = 
  gstat(g = g.mid,
        id = 'fwd1h',
        formula = count ~ 1,
        locations = ~x_abs+y_abs,
        data = fwd %>% filter(mort == 'mid' & timelag_class=='1h'))
g.mid = 
  gstat(g = g.mid,
        id = 'fwd10h',
        formula = count ~ 1,
        locations = ~x_abs+y_abs,
        data = fwd %>% filter(mort == 'mid' & timelag_class=='10h'))
g.mid = 
  gstat(g = g.mid,
        id = 'fwd100h',
        formula = count ~ 1,
        locations = ~x_abs+y_abs,
        data = fwd %>% filter(mort == 'mid' & timelag_class=='100h'))

# add in 1000h
g.mid = 
  gstat(g = g.mid,
        id = 'cwd',
        formula = count ~ 1,
        locations = ~x_abs+y_abs,
        data = cwd %>% filter(mort == 'mid'))

g.mid = 
  gstat(g = g.mid,
        id = 'vegpres',
        formula = veg_present ~ 1,
        locations = ~x_abs+y_abs,
        data = veg %>% filter(mort == 'mid'))

g.mid = 
  gstat(g = g.mid,
        id = 'trees',
        formula = count ~ 1,
        locations = ~x_abs+y_abs,
        data = trees_quadrature %>% filter(mort=='mid'))

g.mid = 
  gstat(g = g.mid,
        id = 'saplings',
        formula = count ~ 1,
        locations = ~x_abs+y_abs,
        data = saplings_quadrature %>% filter(mort=='mid'))

# start with duff
g.high = 
  gstat(id = 'duff',
        formula = duff_cm ~ 1,
        locations = ~x_abs+y_abs,
        data = litterduff %>% filter(mort == 'high'))

# add litter 
g.high = 
  gstat(g = g.high,
        id = 'litter',
        formula = litter_cm ~ 1,
        locations = ~x_abs+y_abs,
        data = litterduff %>% filter(mort == 'high'))

# add 1h, 10h, 100h
g.high = 
  gstat(g = g.high,
        id = 'fwd1h',
        formula = count ~ 1,
        locations = ~x_abs+y_abs,
        data = fwd %>% filter(mort == 'high' & timelag_class=='1h'))
g.high = 
  gstat(g = g.high,
        id = 'fwd10h',
        formula = count ~ 1,
        locations = ~x_abs+y_abs,
        data = fwd %>% filter(mort == 'high' & timelag_class=='10h'))
g.high = 
  gstat(g = g.high,
        id = 'fwd100h',
        formula = count ~ 1,
        locations = ~x_abs+y_abs,
        data = fwd %>% filter(mort == 'high' & timelag_class=='100h'))

# add in 1000h
g.high = 
  gstat(g = g.high,
        id = 'cwd',
        formula = count ~ 1,
        locations = ~x_abs+y_abs,
        data = cwd %>% filter(mort == 'high'))

g.high = 
  gstat(g = g.high,
        id = 'vegpres',
        formula = veg_present ~ 1,
        locations = ~x_abs+y_abs,
        data = veg %>% filter(mort == 'high'))

g.high = 
  gstat(g = g.high,
        id = 'trees',
        formula = count ~ 1,
        locations = ~x_abs+y_abs,
        data = trees_quadrature %>% filter(mort=='high'))

g.high = 
  gstat(g = g.high,
        id = 'saplings',
        formula = count ~ 1,
        locations = ~x_abs+y_abs,
        data = saplings_quadrature %>% filter(mort=='high'))

#### build variograms and cross-variograms #####################################


# 60m: low mortality
vm_low_30 = variogram(g.low, cutoff = 30, width = 1)

plot(vm_low_30)

png(filename = here::here('04-communication',
                          'figures',
                          'manuscript',
                          'variograms_low.png'),
    width = 18, height = 13, units = 'in',
    res = 300)
plot(vm_low_30)
dev.off()


# duff: very fine scale correlation
# litter: very fine cale correlation or not spatially correlated
# fwd1h: correlated
# fwd10h: correlated
# 100h: correlated
# cwd:  maybe correlated
# vegheight: correlated
# ** duff.litter: noisily cross-correlated **
# duff.1h: no cross corr
# duff.10h: no cross-correlation
# duff.100h: no cross.correlation
# duff.cwd: no cross correlation
# duff.vegheight: no cross correlation
# * litter.1h: weak cross correlation *
# litter.10h: no cross corr
# litter.100h: no cross cor
# litter.cwd: no cross cor
# litter.veg no cross cor
# ** 1h.10h: cross corr **
# * 1h.100h noisy cross corr *
# 1h.cwd: no cross corr
# * 1h.vegheight: weird cross cor*
# * 10h.100h: noisy cross cor *
# 10h.cwd: no cross cor
# 10h.veg: no cross cor
# 100h.cwd: no cross cor 
# 100h.veg: no cross cor
# cwd.veg: no cross cor


# 30m: mid mortality
vm_mid_30 = variogram(g.mid, cutoff = 30, width = 1)

plot(vm_mid_30)

png(filename = here::here('04-communication',
                          'figures',
                          'manuscript',
                          'variograms_mid.png'),
    width = 18, height = 13, units = 'in',
    res = 300)
plot(vm_mid_30)
dev.off()

# duff: very fine scale correlation
# litter: very fine cale correlation or not spatially correlated
# fwd1h: correlated
# fwd10h: correlated
# 100h: weakly correlated
# cwd:  maybe correlated
# vegheight: correlated
# * duff.litter: no cross corr *
# duff.1h: no cross-correlation 
# duff.10h: no (reverse?) cross-correlation
# duff.100h: no (reverse?) cross.correlation
# duff.cwd: no (reverse?) cross correlation
# duff.vegheight: no cross correlation
# litter.1h: no cross correlation
# litter.10h: no cross corr
# litter.100h: no cross cor
# litter.cwd: no cross cor
# litter.veg no cross cor
# 1h.10h: no cross corr
# 1h.100h no cross corr 
# 1h.cwd: no cross cor 
# * v1h.vegheight: cross corr *
# ** 10h.100h: cross cor **
# 10h.cwd: no cross cor
# 10h.veg: no cross cor
# 100h.cwd: no cross cor
# 100h.veg: no cross cor
# cwd.veg: no cross cor


# 60m: low mortality
vm_high_30 = variogram(g.high, cutoff = 30, width = 1)

plot(vm_high_30)

png(filename = here::here('04-communication',
                          'figures',
                          'manuscript',
                          'variograms_high.png'),
    width = 18, height = 13, units = 'in',
    res = 300)
plot(vm_high_30)
dev.off()

# duff: very fine scale correlation
# litter: very fine cale correlation or not spatially correlated
# fwd1h: correlated
# fwd10h: correlated
# 100h: correlated
# cwd:  maybe correlated
# vegheight: correlated
# ** duff.litter: noisily cross-correlated **
# duff.1h: no cross corr
# duff.10h: no cross-correlation
# duff.100h: no cross.correlation
# duff.cwd: no cross correlation
# duff.vegheight: no cross correlation
# litter.1h: no cross correlation 
# litter.10h: no cross corr
# litter.100h: no cross cor
# litter.cwd: no cross cor
# litter.veg no cross cor
# ** 1h.10h: cross corr **
# 1h.100h no cross corr
# 1h.cwd: no cross corr
# 1h.vegheight: no cross cor
# * 10h.100h: noisy cross cor *
# 10h.cwd: no cross cor
# *10h.veg: weak cross cor*
# 100h.cwd: no cross cor 
# *100h.veg: weak cross cor* 
# cwd.veg: no cross cor

cor(litterduff$litter_cm, litterduff$duff_cm) # weak correlation between litter and duff

# moderate correlation between 1h and 10h:
cor(
  fwd %>%
    pivot_wider(id_cols = c('plot_id', 'transect_id', 'x_rel', 'y_rel', 'subsample'),
              names_from = 'timelag_class',
              values_from = 'count') %>%
    pull(`1h`),
  fwd %>%
    pivot_wider(id_cols = c('plot_id', 'transect_id', 'x_rel', 'y_rel', 'subsample'),
              names_from = 'timelag_class',
              values_from = 'count') %>%
    pull(`10h`))

# weak correlation between 10h and 100h
cor(
  fwd %>%
    pivot_wider(id_cols = c('plot_id', 'transect_id', 'x_rel', 'y_rel', 'subsample'),
              names_from = 'timelag_class',
              values_from = 'count') %>%
    pull(`10h`),
  fwd %>%
    pivot_wider(id_cols = c('plot_id', 'transect_id', 'x_rel', 'y_rel', 'subsample'),
              names_from = 'timelag_class',
              values_from = 'count') %>%
    pull(`100h`))


# weak correlation between 1h and 10h
cor(
  fwd %>%
    pivot_wider(id_cols = c('plot_id', 'transect_id', 'x_rel', 'y_rel', 'subsample'),
              names_from = 'timelag_class',
              values_from = 'count') %>%
    pull(`1h`),
  fwd %>%
    pivot_wider(id_cols = c('plot_id', 'transect_id', 'x_rel', 'y_rel', 'subsample'),
              names_from = 'timelag_class',
              values_from = 'count') %>%
    pull(`100h`))


# duff: consistent spatial correlation at different scales
# litter: inconsistent spatial correlation at fine scales
# 1h: 2/3 spatially correlated
# 10h: 3/3 spatially correlated
# 100h: weak spatial correlation
# cwd: very fine scale correlation
# veg:  3/3/ spatially correlated
# duff.litter: noisy/weak/present corr
# duff.1h: no/no/no
# duff.10h: no/no/no
# duff.100h: no/weak/no
# duff.cwd: no/no/no
# duff.veg: no/no/no
# litter.1h: weak-course/no/no
# litter.10h: no/no/no
# litter.100h: no/no/no
# litter.cwd: no/no/no
# litter.veg: no/no/no
# 1h.10h: yes/noisy/noisy
# 1h.100h: yes/noisy/no
# 1h.cwd: no/no/no
# 1h.veg: yes/yes/no
# 10h.100h: noisy/yes/noisy
# 10h.cwd: no/no/no
# 10h.veg: no/no/weak
# 100h.cwd: no/no/no
# 100h.veg: no/no/maybe
# cwd.veg: no/no/no

# thoughts:
# variograms suggest some cross correlation between litter and duff, 
#  show pretty clear cross correlation between 1h and 10h fuels, and suggest 
#  cross correlation between 1h and 100h fuels. They also suggest correlation 
#  between 1h fuels and veg heights. A bivariate response litter and duff model 
# with cross correlation between the two (and separate within-component GPs) 
# doesn't converge, and the correlation between litter and duff is weak (~0.18). 
# there's stronger correlation between 1h and 10h fuels, but similarly weak 
# correlation between 1h and 100h and between 1h and 10h. (does a cross component 
# correlation fit for FWD?). Given the difficulty of estimating the parameters for 
# the cross-correlation type models, I'm going to ignore cross correlation 
# for litter and duff (where its relatively weak in any case) 
# again, these are limitations with the generative statistical model approach, 
# where subtle behaviors and dependencies can't always be included in the model. 
# A complete census data collection approach would result in fewer fuelbeds to 
# run fire models on, but the fuelbeds themselves would have all of these subtle 
# behaviors baked in in a realistic way. 

# I'm probably going to include the cross correlation between FWD classes by 
# modelling the total count of FWD particles, rather than a per-timelag approach. 
# Pending a cross correlation model also not working for FWD. 

# for now, I'm not planning on trying to do a fwd:veg interaction model because 
# the response distributions have to be different for those and I'm not really 
# sure how to represent cross correlation in that case (I can give them a common 
# spatial random effect, but would it be meaningful or realistic with the 
# different responses on different scales?).
