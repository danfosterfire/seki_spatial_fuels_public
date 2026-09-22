library(tidyverse)

seki.litterduff = 
  readRDS(here::here('02-data',
                     '03-clean',
                     'seki_litterduff.rds'))

head(seki.litterduff)

seki.litterduff = 
  seki.litterduff %>%
  mutate(x_center = 0, y_center = 0) %>%
  filter(plot_id == '01') %>%
  mutate(
    # not really sure why this works but it does translate between compass-angle
    # and unit-circle angle correctly to just use sin for the x and cos for the 
    # y...
    x_coord = x_center + (location_m*sin(az*(pi/180))),
    y_coord = y_center + (location_m*cos(az*(pi/180)))
  ) %>%
  select(-x_center, -y_center)


seki.fwd_tallies = 
  readRDS(here::here('02-data',
                     '03-clean',
                     'seki_fwdtallies.rds'))

seki.fwd_tallies = 
  seki.fwd_tallies %>%
  mutate(x_center = 0, y_center = 0) %>%
  filter(plot_id == '01') %>%
  mutate(x_coord = x_center + location_m*sin(az*(pi/180)),
         y_coord = y_center + location_m*cos(az*(pi/180))) %>%
  select(-x_center, -y_center)

seki.fwd_diams = 
  seki.fwd_tallies %>%
  filter(location_m == 9.5)

litterduff_pairwise = 
  seki.litterduff %>%
  select(x_coord, y_coord) %>%
  as.matrix() %>%
  dist()

litterduff_pairwise = 
  tibble(id = 1:length(litterduff_pairwise),
         distance_m = litterduff_pairwise)

fwd_pairwise = 
  seki.fwd_tallies %>%
  select(x_coord, y_coord) %>%
  as.matrix() %>%
  dist()

fwd_pairwise = 
  tibble(id = 1:length(fwd_pairwise),
         distance_m = fwd_pairwise)

combined_locations = 
  seki.litterduff %>%
  select(x_coord, y_coord) %>%
  mutate(component = 'Litter+Duff depths') %>%
  bind_rows(
    seki.fwd_tallies %>%
      select(x_coord, y_coord) %>%
      mutate(component = 'FWD tallies')
  ) %>%
  bind_rows(
    seki.fwd_diams %>%
      select(x_coord, y_coord) %>%
      mutate(component = 'FWD diameters')
  )

plot_overview_gg = 
  ggplot()+
  geom_point(
    data = combined_locations,
    aes(x = x_coord, y = y_coord, pch = component),
    size = 2,
    alpha = 0.75
  )+
  coord_fixed()+
  theme_bw()+
  labs(x = 'Meters east from plot center',
       y = 'Meters north from plot center',
       pch = 'Fuel loading observations')+
  scale_shape_manual(values = c(19, 3, 1))+
  theme(legend.position = 'bottom',
        legend.direction = 'vertical')


plot_overview_gg

litterduff_histogram = 
  ggplot(data = litterduff_pairwise,
         aes(x = distance_m))+
  geom_histogram()+
  theme_bw()+
  labs(
    title = 'Litter+duff pairwise distances',
    x = 'Distance between observations (m)',
    y = 'Number of pairs'
  )

litterduff_histogram

fwd_histogram = 
  ggplot(data = fwd_pairwise,
         aes(x = distance_m))+
  geom_histogram()+
  theme_bw()+
  labs(
    title = 'FWD tallies pairwise distances',
    x = 'Distance between observations (m)',
    y = 'Number of pairs'
  )

fwd_histogram

library(cowplot)

histograms = 
  cowplot::plot_grid(litterduff_histogram, fwd_histogram,
                     ncol = 1)


overall_figure = 
  cowplot::plot_grid(plot_overview_gg, histograms,
                     rel_widths = c(1, 0.75), rel_heights = c(1, 1),
                     nrow = 1)

overall_figure

ggsave(overall_figure,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'plot_overview_litterduff_fwd.png'),
       height = 5, width = 9, units = 'in')


