library(here)
library(viridis)
library(tidyverse)
library(rgl)


duff = 
  read_csv(here::here('02-data',
                      '06-results',
                      'simulated_fuelbeds',
                      'duff.csv'))

litter = 
  read_csv(here::here('02-data',
                      '06-results',
                      'simulated_fuelbeds',
                      'litter.csv'))

fwd = 
  read_csv(here::here('02-data',
                      '06-results',
                      'simulated_fuelbeds',
                      'fwd.csv'))

cwd = 
  read_csv(here::here('02-data',
                      '06-results',
                      'simulated_fuelbeds',
                      'cwd.csv'))

understory = 
  read_csv(here::here('02-data',
                      '06-results',
                      'simulated_fuelbeds',
                      'understory.csv'))

trees_saplings = 
  read_csv(here::here('02-data',
                      '06-results',
                      'simulated_fuelbeds',
                      'trees_saplings.csv'))

#### build 2d images ###########################################################


#### stack to 3d ###############################################################



#### scratch ###################################################################

library(data.table)
library(viridis)
library(rgl)


bulkdensity_array = 
  bulkdensity_trees+bulkdensity_shrubs

library(viridis)
library(rgl)

voxels.bulkdensity = as.data.table(bulkdensity_array)
voxels.bulkdensity.thinned = voxels.bulkdensity[value>0,]
voxels.bulkdensity.thinned$x = as.numeric(voxels.bulkdensity.thinned$x)
voxels.bulkdensity.thinned$y = as.numeric(voxels.bulkdensity.thinned$y)
voxels.bulkdensity.thinned$z = as.numeric(voxels.bulkdensity.thinned$z)


head(voxels.bulkdensity.thinned)



voxels.bulkdensity.thinned$vcolor = 
  map_viridis(x = voxels.bulkdensity.thinned$value,
              res = 2, begin = 0.8, end = 0)

install.packages('VoxR')



VoxR::plot_voxels(data = voxels.bulkdensity.thinned,
                  alpha = 0.5)

library(rgl)

# THIS IS HOW TO DO IT!!!!!!!!!!!!!!!!!!!#############################



# NEXT STEPS:
# set a seed for each fuel component individually to get nice looking ones
# run the trees and save it 
# write out the long dataframes or data.tables for the surface fuels and the voxels
# do a plotting script, based off of 'make_voxel_plot.R' but using the above approach 
# to get beautiful transparent viridisfied cubes for the shrubs and the trees/saplings
# (separately)


#### simulate grids of surface fuels ###########################################

# make nice viridis color pallette
map_viridis = 
  function(x, res, begin = 0, end = 1, opt = 'D', direction = 1){
    
    # define the sequence
    minval = min(x)
    maxval = max(x)
    
    scaled_values = 
      as.integer(round((x-minval) / (maxval-minval), res) * 10 * res)
    
    viridis_colors = viridis(n = 10*res+1, begin = begin, end = end, option = opt,
                             direction = direction)
    
    return(viridis_colors[scaled_values+1])
    
  }

litter_surface = 
  ggplot(data = 
           sim_results %>%
           filter(component=='litterduff'),
         aes(x = x, y = y, fill = capped_mgha))+
  geom_tile()+
  scale_fill_viridis_c(option = 'B', direction = -1)+
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
  theme_minimal()+
  coord_fixed()+
  theme(text = element_text(size = 14),
        legend.direction = 'vertical',
        legend.text = element_text(size = 12))+
  labs(fill = 'Mg/ha')
litter_surface

fwd_surface = 
  
  ggplot(data = 
           sim_results %>%
           filter(component=='fwd'),
         aes(x = x, y = y, fill = capped_mgha))+
  geom_tile()+
  scale_fill_viridis_c(option = 'B', direction = -1)+
  theme_minimal()+
  coord_fixed()+  
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
  theme(text = element_text(size = 14),
        legend.direction = 'vertical',
        legend.text = element_text(size = 12))+
  labs(fill = 'Mg/ha')

fwd_surface

seki.metadata = 
  readRDS(here::here('02-data',
                     '04-geolocated',
                     'seki_metadata.rds'))

trees_pp = 
  ggplot(data = 
           treepoints %>%
           left_join(seki.metadata %>%
                       dplyr::select(plot_id, x_center, y_center)) %>%
           mutate(x_relative = x_coord-x_center, y_relative = y_coord-y_center),
         aes(x = x_relative, y = y_relative, size = dbh_cm))+
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
  theme_bw()+
  labs(x = 'Easting (m)', y = 'Northing (m)', color = 'Height (m)', 
       size = 'DBH (cm)', title = 'Trees and saplings')
trees_pp


#### plot litterduff ###########################################################

duffz = -0.5
fwdz = -.25
treesppz = -0.05




# plot litterduff+fwd


# plot litterduff+fwd+tree point pattern


# plot litterduff+fwd+voxels


# save movie 





# top outline
rgl.lines(x = c(-15, -5), z = c(-5, -5), y = c(30,30), color = 'black')
rgl.lines(x = c(-15, -5), z = c(5, 5), y = c(30, 30), color = 'black')
rgl.lines(x = c(-15, -15), z = c(5, -5), y = c(30, 30), color = 'black')
rgl.lines(x = c(15, 5), z = c(-5, -5), y = c(30,30), color = 'black')
rgl.lines(x = c(15, 5), z = c(5, 5), y = c(30, 30), color = 'black')
rgl.lines(x = c(15, 15), z = c(5, -5), y = c(30, 30), color = 'black')
rgl.lines(z = c(-15, -5), x = c(-5, -5), y = c(30,30), color = 'black')
rgl.lines(z = c(-15, -5), x = c(5, 5), y = c(30, 30), color = 'black')
rgl.lines(z = c(-15, -15), x = c(5, -5), y = c(30, 30), color = 'black')
rgl.lines(z = c(15, 5), x = c(-5, -5), y = c(30,30), color = 'black')
rgl.lines(z = c(15, 5), x = c(5, 5), y = c(30, 30), color = 'black')
rgl.lines(z = c(15, 15), x = c(5, -5), y = c(30, 30), color = 'black')

# vertical joins 
rgl.lines(x = c(-15, -15), z = c(-5, -5), y = c(0.1,30), color = 'black')
rgl.lines(x = c(-15, -15), z = c(5, 5), y = c(0.1,30), color = 'black')
rgl.lines(x = c(15, 15), z = c(-5, -5), y = c(0.1,30), color = 'black')
rgl.lines(x = c(15, 15), z = c(5, 5), y = c(0.1,30), color = 'black')
rgl.lines(z = c(-15, -15), x = c(-5, -5), y = c(0.1,30), color = 'black')
rgl.lines(z = c(-15, -15), x = c(5, 5), y = c(0.1,30), color = 'black')
rgl.lines(z = c(15, 15), x = c(-5, -5), y = c(0.1,30), color = 'black')
rgl.lines(z = c(15, 15), x = c(5, 5), y = c(0.1,30), color = 'black')
rgl.lines(x = c(-5, -5), z = c(-5, -5), y = c(0.1,30), color = 'black')
rgl.lines(x = c(-5, -5), z = c(5, 5), y = c(0.1,30), color = 'black')
rgl.lines(x = c(5, 5), z = c(-5, -5), y = c(0.1,30), color = 'black')
rgl.lines(x = c(5, 5), z = c(5, 5), y = c(0.1,30), color = 'black')
rgl.lines(z = c(-5, -5), x = c(-5, -5), y = c(0.1,30), color = 'black')
rgl.lines(z = c(-5, -5), x = c(5, 5), y = c(0.1,30), color = 'black')
rgl.lines(z = c(5, 5), x = c(-5, -5), y = c(0.1,30), color = 'black')
rgl.lines(z = c(5, 5), x = c(5, 5), y = c(0.1,30), color = 'black')



                  
