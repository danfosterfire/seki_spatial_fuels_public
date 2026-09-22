library(here)
library(tidyverse)


#### litter ################################################################

litter_training_data = 
  readRDS(here::here('02-data',
                     '05-for_analysis',
                     'litter_training_data.rds'))

litter_df = 
  data.frame(Y = litter_training_data$Y,
             plot_id = litter_training_data$plot_id,
             location_id = litter_training_data$location_id,
             group_id = litter_training_data$group_id[litter_training_data$plot_id])


# y distribution
ggplot(litter_df,
       aes(x = Y))+
  geom_bar()

# x distribution
ggplot(litter_df,
       aes(x = group_id))+
  geom_bar()

ggplot(litter_df,
       aes(x = plot_id))+
  geom_bar()

ggplot(litter_df,
       aes(x = location_id))+
  geom_bar()

# XY distribution
ggplot(litter_df,
       aes(x = Y))+
  geom_bar()+
  facet_grid(group_id~.)



#### duff ################################################################

duff_training_data = 
  readRDS(here::here('02-data',
                     '05-for_analysis',
                     'duff_training_data.rds'))

duff_df = 
  data.frame(Y = duff_training_data$Y,
             plot_id = duff_training_data$plot_id,
             location_id = duff_training_data$location_id,
             group_id = duff_training_data$group_id[duff_training_data$plot_id])


# y distribution
ggplot(duff_df,
       aes(x = Y))+
  geom_bar()

# x distribution
ggplot(duff_df,
       aes(x = group_id))+
  geom_bar()

ggplot(duff_df,
       aes(x = plot_id))+
  geom_bar()

ggplot(duff_df,
       aes(x = location_id))+
  geom_bar()

# XY distribution
ggplot(duff_df,
       aes(x = Y))+
  geom_bar()+
  facet_grid(group_id~.)



#### fwdtallies ################################################################

fwdtallies_training_data = 
  readRDS(here::here('02-data',
                     '05-for_analysis',
                     'fwdtallies_training_data.rds'))

fwdtallies_df = 
  data.frame(Y = fwdtallies_training_data$Y,
             plot_id = fwdtallies_training_data$plot_id,
             location_id = fwdtallies_training_data$location_id,
             group_id = fwdtallies_training_data$group_id[fwdtallies_training_data$plot_id])


# y distribution
ggplot(fwdtallies_df,
       aes(x = Y))+
  geom_bar()

# x distribution
ggplot(fwdtallies_df,
       aes(x = group_id))+
  geom_bar()

ggplot(fwdtallies_df,
       aes(x = plot_id))+
  geom_bar()

ggplot(fwdtallies_df,
       aes(x = location_id))+
  geom_bar()

# XY distribution
ggplot(fwdtallies_df,
       aes(x = Y))+
  geom_bar()+
  facet_grid(group_id~.)

#### fwd diameters #############################################################

fwddiams_training_data = 
  readRDS(here::here('02-data',
                     '05-for_analysis',
                     'fwddiams_training_data.rds'))

fwddiams_df = 
  data.frame(Y = 
               fwddiams_training_data$Y,
             group_id = 
               fwddiams_training_data$group_id[
                 fwddiams_training_data$plot_id[
                   fwddiams_training_data$transect_id]],
             plot_id = 
               fwddiams_training_data$plot_id[
                 fwddiams_training_data$transect_id],
             transect_id = 
               fwddiams_training_data$transect_id)

# y distribution
ggplot(data = fwddiams_df,
       aes(x = Y))+
  geom_histogram()

# x distributions
ggplot(data = fwddiams_df,
       aes(x = group_id))+
  geom_bar()

ggplot(data = fwddiams_df,
       aes(x = plot_id))+
  geom_bar()

ggplot(data = fwddiams_df,
       aes(x = transect_id))+
  geom_bar()

# XY distributions
ggplot(data = fwddiams_df,
       aes(x = Y, color = as.factor(group_id)))+
  geom_density(lwd = 1)+
  theme_minimal()

ggplot(data = fwddiams_df,
       aes(x = Y, color = as.factor(plot_id)))+
  geom_density(lwd = 1, alpha = 0.5)+
  theme_minimal()

ggplot(data = fwddiams_df,
       aes(x = Y, color = as.factor(transect_id)))+
  geom_density(lwd = 1, alpha = 0.5)+
  theme_minimal()

### CWD tallies ################################################################

cwdtallies_training_data = 
  readRDS(here::here('02-data',
                     '05-for_analysis',
                     'cwdtallies_training_data.rds'))

names(cwdtallies_training_data)

cwdtallies_training_data$N

cwdtallies_training_data$P

cwdtallies_training_data$TR

cwdtallies_training_data$L

cwdtallies_training_data$G

cwdtallies.df = 
  data.frame(Y = cwdtallies_training_data$Y,
             location_id = cwdtallies_training_data$location_id,
             transect_id = cwdtallies_training_data$transect_id,
             plot_id = cwdtallies_training_data$plot_id[cwdtallies_training_data$transect_id],
             group_id = 
               cwdtallies_training_data$group_id[
                 cwdtallies_training_data$plot_id[
                   cwdtallies_training_data$transect_id]])

# y distribution
ggplot(data = cwdtallies.df,
       aes(x = Y))+
  geom_bar()

# X distributions
ggplot(data = cwdtallies.df,
       aes(x = location_id))+
  geom_bar()

ggplot(data = cwdtallies.df,
       aes(x = transect_id))+
  geom_bar()

ggplot(data = cwdtallies.df,
       aes(x = plot_id))+
  geom_bar()

ggplot(data = cwdtallies.df,
       aes(x = group_id))+
  geom_bar()

# XY distributions
ggplot(data = cwdtallies.df,
       aes(x = Y))+
  geom_bar()+
  facet_grid(group_id~.)

#### cwd diameters #############################################################

cwddiams_training_data = 
  readRDS(here::here('02-data',
                     '05-for_analysis',
                     'cwddiams_training_data.rds'))

cwddiams_df = 
  data.frame(Y = 
               cwddiams_training_data$Y,
             group_id = 
               cwddiams_training_data$group_id[
                 cwddiams_training_data$plot_id],
             plot_id = 
               cwddiams_training_data$plot_id)

# y distribution
ggplot(data = cwddiams_df,
       aes(x = Y))+
  geom_histogram()

# x distributions
ggplot(data = cwddiams_df,
       aes(x = group_id))+
  geom_bar()

ggplot(data = cwddiams_df,
       aes(x = plot_id))+
  geom_bar()


# XY distributions
ggplot(data = cwddiams_df,
       aes(x = Y, color = as.factor(group_id)))+
  geom_density(lwd = 1)+
  theme_minimal()

ggplot(data = cwddiams_df,
       aes(x = Y, color = as.factor(plot_id)))+
  geom_density(lwd = 1, alpha = 0.5)+
  theme_minimal()


#### veg presence/absence ######################################################

vegpa_training_data = 
  readRDS(here::here('02-data',
                     '05-for_analysis',
                     'vegpa_training_data.rds'))

vegpa_training_data$N
vegpa_training_data$P
vegpa_training_data$L
vegpa_training_data$G
vegpa_training_data$coords

vegpa.df = 
  data.frame(
    Y = vegpa_training_data$Y,
    plot_id = vegpa_training_data$plot_id,
    location_id = vegpa_training_data$location_id,
    group_id = vegpa_training_data$group_id[vegpa_training_data$plot_id]
  )

# y distribution
ggplot(vegpa.df,
       aes(x = Y))+
  geom_bar()

# x distributions
ggplot(vegpa.df,
       aes(x = plot_id))+
  geom_bar()
ggplot(vegpa.df,
       aes(x = location_id))+
  geom_bar()
ggplot(vegpa.df,
       aes(x = group_id))+
  geom_bar()

# XY distributions
ggplot(vegpa.df,
       aes(x = plot_id, fill = as.factor(Y)))+
  geom_bar(position = position_fill())
ggplot(vegpa.df,
       aes(x = location_id, fill = as.factor(Y)))+
  geom_bar(position = position_fill())
ggplot(vegpa.df,
       aes(x = group_id, fill = as.factor(Y)))+
  geom_bar(position = position_fill())


#### trees #####################################################################


trees_training_data = 
  readRDS(here::here('02-data',
                     '05-for_analysis',
                     'trees_training_data.rds'))

trees_df = 
  data.frame(Y = trees_training_data$Y,
             plot_id = trees_training_data$plot_id,
             location_id = trees_training_data$location_id,
             group_id = trees_training_data$group_id[trees_training_data$plot_id])


# y distribution
ggplot(trees_df,
       aes(x = Y))+
  geom_bar()

# x distribution
ggplot(trees_df,
       aes(x = group_id))+
  geom_bar()

ggplot(trees_df,
       aes(x = plot_id))+
  geom_bar()

ggplot(trees_df,
       aes(x = location_id))+
  geom_bar()

# XY distribution
ggplot(trees_df,
       aes(x = Y))+
  geom_bar()+
  facet_grid(group_id~.)

#### saplings ##################################################################


saplings_training_data = 
  readRDS(here::here('02-data',
                     '05-for_analysis',
                     'saplings_training_data.rds'))

saplings_df = 
  data.frame(Y = saplings_training_data$Y,
             plot_id = saplings_training_data$plot_id,
             location_id = saplings_training_data$location_id,
             group_id = saplings_training_data$group_id[saplings_training_data$plot_id])


# y distribution
ggplot(saplings_df,
       aes(x = Y))+
  geom_bar()

# x distribution
ggplot(saplings_df,
       aes(x = group_id))+
  geom_bar()

ggplot(saplings_df,
       aes(x = plot_id))+
  geom_bar()

ggplot(saplings_df,
       aes(x = location_id))+
  geom_bar()

# XY distribution
ggplot(saplings_df,
       aes(x = Y))+
  geom_bar()+
  facet_grid(group_id~.)




