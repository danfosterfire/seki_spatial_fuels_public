
#### setup #####################################################################

library(here)
library(tidyverse)
library(readxl)

# load data
seki.metadata = 
  readxl::read_excel(here::here('02-data',
                                '01-raw_field_sheets',
                                'datasheet_digital_seki_080122_completed.xlsx'),
                     sheet = 'plot_metadata',
                     col_names = 
                       c('plot', 'az', 'slp_degrees',
                         'cwd', 'fwd', 'litter_duff', 'veg', 'trees',
                         'inv_date', 'comment', 'comment2', 'comment3'),
                     col_types = 
                       c('text', 'text', 'numeric', 
                         'text', 'text', 'text', 'text', 'text', 
                         'date', 'text', 'text', 'text'),
                     skip = 1) %>%
  mutate(inv_date = as.Date(inv_date, format = '%m/%d/%Y'))

seki.fwd_tallies = 
  readxl::read_excel(here::here('02-data',
                              '01-raw_field_sheets',
                              'datasheet_digital_seki_080122_completed.xlsx'),
                   sheet = 'fwd_tallies',
                   col_names = 
                     c('plot', 'az', 'location_m',
                       'a1h', 'a10h', 'a100h', 'b1h', 'b10h', 'b100h', 'comments'),
                   col_types = 
                     c('text', 'text', 'numeric', 
                       'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric',
                       'text'),
                   skip = 1)


seki.fwd_diams = 
  
  readxl::read_excel(here::here('02-data',
                                '01-raw_field_sheets',
                                'datasheet_digital_seki_080122_completed.xlsx'),
                     sheet = 'fwd_diams',
                     col_names = 
                       c('plot', 'az', 'diam_mm'),
                     col_types = 
                       c('text', 'text', 'numeric'),
                     skip = 1)

seki.litterduff = 
  
  readxl::read_excel(here::here('02-data',
                                '01-raw_field_sheets',
                                'datasheet_digital_seki_080122_completed.xlsx'),
                     sheet = 'litterduff',
                     col_names = 
                       c('plot', 'az', 'location_m', 'litter_cm',
                         'duff_cm', 'fuel_cm', 'comments'),
                     col_types = 
                       c('text', 'text',
                         'numeric',
                         'numeric', 'numeric', 'numeric',
                         'text'),
                     skip = 1)

seki.cwd = 
  readxl::read_excel(here::here('02-data',
                                '01-raw_field_sheets',
                                'datasheet_digital_seki_080122_completed.xlsx'),
                     sheet = 'cwd',
                     col_names = 
                       c('plot', 'az', 'location_m',
                         'diam_cm', 'decay', 'comments'),
                     col_types = 
                       c('text', 'text', 'numeric',
                         'numeric', 'text', 'text'),
                     skip = 1)

seki.veg = 
  readxl::read_excel(here::here('02-data',
                                '01-raw_field_sheets',
                                'datasheet_digital_seki_080122_completed.xlsx'),
                     sheet = 'veg',
                     col_names = 
                       c('plot', 'az', 'start_m', 'end_m',
                         'spp', 'status', 'height_m', 'comments'),
                     col_types = 
                       c('text', 'text', 'numeric', 'numeric',
                         'text', 'text', 'numeric', 'text'),
                     skip = 1)

seki.trees = 
  readxl::read_excel(here::here('02-data',
                                '01-raw_field_sheets',
                                'datasheet_digital_seki_080122_completed.xlsx'),
                     sheet = 'trees',
                     col_names = 
                       c('plot', 'ts', 'status', 'spp', 'location_m',
                         'distance_m', 'dbh_cm', 'height_m', 'htcb_m', 'decay',
                         'comments'),
                     col_types = 
                       c('text', 'text', 'text', 'text', 'numeric',
                         'numeric', 'numeric', 'numeric', 'numeric', 'text',
                         'text'),
                     skip = 1)


#### check comments ############################################################

seki.cwd %>%
  filter(!is.na(comments))

seki.fwd_tallies %>%
  filter(!is.na(comments)) %>%
  print(n = Inf, width = Inf)
# the crew added comments where there was some obstruction on a fwd subtransect;
# set these measurements to NA
seki.fwd_tallies %>%
  filter((!is.na(plot)& !is.na(az) & !is.na(location_m))&
           ((plot=='1'&az=='270'&location_m==4.5)|
           (plot=='14'&az=='0'&location_m==2.5)|
           (plot=='8'&az=='0'&location_m==9.5)|
           (plot=='5'&az=='270'&location_m==5.5)|
           (plot=='4'&az=='270'&location_m==19.5)|
           (plot=='15'&az=='90'&location_m==19.5)|
           (plot=='11'&az=='180'&location_m==7.5)|
           (plot=='11'&az=='180'&location_m==9.5)|
           (plot=='11'&az=='270'&location_m==5.5)|
           (plot=='13'&az=='180'&location_m==19.5)))

seki.fwd_tallies = 
  seki.fwd_tallies %>%
  mutate(
    a1h = 
      ifelse((!is.na(plot)& !is.na(az) & !is.na(location_m))&
               ((plot=='5'&az=='270'&location_m==5.5)|
                  (plot=='4'&az=='270'&location_m==19.5)),
             NA,
             a1h),
    a10h = 
      ifelse((!is.na(plot)& !is.na(az) & !is.na(location_m))&
               ((plot=='5'&az=='270'&location_m==5.5)|
                  (plot=='4'&az=='270'&location_m==19.5)),
             NA,
             a10h),
    
    a100h = 
      ifelse((!is.na(plot)& !is.na(az) & !is.na(location_m))&
               ((plot=='5'&az=='270'&location_m==5.5)|
                  (plot=='4'&az=='270'&location_m==19.5)),
             NA,
             a100h),
    
    b1h = 
      ifelse((!is.na(plot)& !is.na(az) & !is.na(location_m))&
               ((plot=='1'&az=='270'&location_m==4.5)|
                  (plot=='14'&az=='0'&location_m==2.5)|
                  (plot=='8'&az=='0'&location_m==9.5)|
                  (plot=='5'&az=='270'&location_m==5.5)|
                  (plot=='4'&az=='270'&location_m==19.5)|
                  (plot=='15'&az=='90'&location_m==19.5)|
                  (plot=='11'&az=='180'&location_m==7.5)|
                  #(plot=='11'&az=='180'&location_m==9.5)|
                  (plot=='11'&az=='270'&location_m==5.5)|
                  (plot=='13'&az=='180'&location_m==19.5)),
             NA,
             b1h),
    
    b10h = 
      ifelse((!is.na(plot)& !is.na(az) & !is.na(location_m))&
               ((plot=='1'&az=='270'&location_m==4.5)|
                  (plot=='14'&az=='0'&location_m==2.5)|
                  (plot=='8'&az=='0'&location_m==9.5)|
                  (plot=='5'&az=='270'&location_m==5.5)|
                  (plot=='4'&az=='270'&location_m==19.5)|
                  (plot=='15'&az=='90'&location_m==19.5)|
                  (plot=='11'&az=='180'&location_m==7.5)|
                  #(plot=='11'&az=='180'&location_m==9.5)|
                  (plot=='11'&az=='270'&location_m==5.5)|
                  (plot=='13'&az=='180'&location_m==19.5)),
             NA,
             b10h),
    
    b100h = 
      ifelse((!is.na(plot)& !is.na(az) & !is.na(location_m))&
               ((plot=='1'&az=='270'&location_m==4.5)|
                  (plot=='14'&az=='0'&location_m==2.5)|
                  (plot=='8'&az=='0'&location_m==9.5)|
                  (plot=='5'&az=='270'&location_m==5.5)|
                  (plot=='4'&az=='270'&location_m==19.5)|
                  (plot=='15'&az=='90'&location_m==19.5)|
                  (plot=='11'&az=='180'&location_m==7.5)|
                  #(plot=='11'&az=='180'&location_m==9.5)|
                  (plot=='11'&az=='270'&location_m==5.5)|
                  (plot=='13'&az=='180'&location_m==19.5)),
             NA,
             b100h))

seki.litterduff %>%
  filter(!is.na(comments)) %>%
  print(n = Inf)

# I mostly agree with the crew's decisions on interpretation here, except for 
# interpreting "trace litter" as NA. I'm going to intrepret those as 0s:
seki.litterduff %>%
  filter(str_detect(comments, 'TRACE LITTER'))
seki.litterduff = 
  seki.litterduff %>%
  mutate(litter_cm = 
           ifelse(!is.na(comments)&str_detect(comments, 'TRACE LITTER'),
                  0,
                  litter_cm))

seki.trees %>%
  filter(!is.na(comments)) %>%
  print(n = Inf, width = Inf)

# if it's a fill-in row for "no trees present", just ditch it
seki.trees = 
  seki.trees %>%
  filter(!(!is.na(comments)&str_detect(comments, 'NO TREES')))

# note missing decay for 25-SE-9.7-3.6, which is coded as "NA" text
# not worried about broken/dead tops
# i'm trusting the dist/location flip
# going to leave 35-se-9.6-2.8 as-is

seki.veg %>%
  filter(!is.na(comments))

seki.veg = 
  seki.veg %>%
  filter(!(!is.na(comments) & 
           str_detect(comments, 'NONE|NO SHRUBS')))

#### check for missing ID info #################################################

# all observations should have complete ID info (plot, azimuth, potentially 
# location and distance)
seki.cwd %>%
  filter(is.na(plot)|is.na(az))

seki.fwd_diams %>%
  filter(is.na(plot)|is.na(az))

seki.fwd_tallies %>%
  filter(is.na(plot)|is.na(az)|is.na(location_m))
unique(seki.fwd_tallies$location_m)[order(unique(seki.fwd_tallies$location_m))]

seki.litterduff %>%
  filter(is.na(plot)|is.na(az)|is.na(location_m))
unique(seki.litterduff$location_m)[order(unique(seki.litterduff$location_m))]

seki.metadata %>%
  filter(is.na(plot)|is.na(az))

seki.trees %>%
  filter(is.na(plot)|is.na(ts)|is.na(location_m)|is.na(distance_m)|
           ts=='NA') %>%
  print(n = Inf)




#### check for missing observation info ########################################

# all cwd observations should have decay and diam_cm
seki.cwd %>%
  filter(is.na(decay)|is.na(diam_cm)|decay=='NA'|is.na(location_m))
# drop the fill-in rows with NA location and diam
seki.cwd = 
  seki.cwd %>%
  filter(!is.na(location_m))
# interpolate the missing decay class by pulling from the observed distribution
# of the classes for the rest of the data
decay_probs = 
  seki.cwd %>%
  filter(!is.na(decay)&decay!='NA') %>%
  group_by(decay) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(id = '1') %>%
  left_join(.,
            seki.cwd %>%
              summarise(total = n()) %>%
              ungroup() %>%
              mutate(id = '1')) %>%
  mutate(p_class = count / total)

# Comment by LMR: This seems like a questionable method. It's not going to affect 
# your results whatsoever with only one observation, but there are better ways to 
# deal with missing data. I'll probably have much more missing data than you, so 
# it'll be good for me to do it the correct way. Essentially, I'd include another 
# line in the model to account for the variable with missing data. The line 
# would function as a likelihood for observed data, and as a prior for missing 
# data. McElreath's book is good at explaining it.
# Edit: I see other missing data was filled in based on the most likely value. I 
# doubt it affects things much at all, and you could make arguments about human
# error that we never account for, but it's not the most honest to treat 
# guessed values as true observations when there are statistical methods to deal 
# with them. 
set.seed(110819)
seki.cwd[seki.cwd$decay=='NA','decay'] = 
  sample(c('1', '2', '3', '4', '5'), 
       prob = decay_probs$p_class,
       size = 1)[1]

# all fwd_diams should have a diam_mm
seki.fwd_diams %>%
  filter(is.na(diam_mm))
seki.fwd_diams %>%
  filter(plot=='1'&az=='180')
# drop the fill-in
seki.fwd_diams = 
  seki.fwd_diams %>%
  filter(!is.na(diam_mm))

# fwd_tallies can have NAs

# litterduff can have NAs

# all metadata should have slp_degrees
seki.metadata %>%
  filter(is.na(slp_degrees)) 
# ditch the undone transect and fill in the missing slope with the mean
# slope
seki.metadata[seki.metadata$plot=='17'&seki.metadata$az=='0','slp_degrees'] = 
  mean(seki.metadata$slp_degrees, na.rm = TRUE)
seki.metadata = 
  seki.metadata %>%
  filter(!is.na(slp_degrees))

# all trees shoudl have spp, status, dbh_cm, height_m
seki.trees %>%
  filter(is.na(spp)|spp=='NA'|is.na(status)|status=='NA'|
           is.na(dbh_cm)|is.na(height_m))

# all live trees should have htcb_m
seki.trees %>%
  filter(status=='L'&is.na(htcb_m))
# interpolate the missing values using a generic (all spp) equation
ggplot(seki.trees,
       aes(x = height_m, y = htcb_m))+
  geom_point()+
  geom_smooth(method = 'lm', formula = y~I(x)+I(x**2)+I(x**3))
set.seed(110819)
htcb_fit = 
  lm(data = seki.trees,
     htcb_m ~ I(height_m)+I(height_m**2)+I(height_m**3))
seki.trees[seki.trees$status=='L'&
             is.na(seki.trees$htcb_m),'htcb_m'] = 
  predict(htcb_fit, 
          newdata = seki.trees %>% 
            filter(status=='L'&is.na(htcb_m)))

# all snags should have decay
seki.trees %>%
  filter(status=='D'&(is.na(decay)|decay=='NA'))
# fill in the blank by pulling from the observed distribution
decay_probs_trees = 
  seki.trees %>%
  filter(!(status=='D'&(is.na(decay)|decay=='NA'))) %>%
  filter(status=='D') %>%
  group_by(decay) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(id = '1') %>%
  left_join(.,
            seki.trees %>%
              filter(status=='D'&!(is.na(decay)|decay=='NA')) %>%
              summarise(total = n()) %>%
              ungroup() %>%
              mutate(id = '1')) %>%
  mutate(p_class = count / total)
set.seed(110819)
seki.trees[seki.trees$status=='D'&
             seki.trees$decay=='NA', 'decay'] = 
  sample(x = decay_probs_trees$decay,
         prob = decay_probs_trees$p_class,
         size = 1)

# all veg should have start_m, end_m, status, spp, height_m
seki.veg %>%
  filter(is.na(start_m)|is.na(end_m)|is.na(status)|status=='NA'|
           is.na(spp)|spp=='NA'|is.na(height_m))


#### initial cleaning ##########################################################

# drop NA rows in fwd_diams; these were inserted for samples with 0 intersections
seki.fwd_diams = 
  seki.fwd_diams %>%
  filter(!is.na(diam_mm))

# drop NA location rows in CWD; also inserted for samples with 0 intersections
seki.cwd = 
  seki.cwd %>%
  filter(!is.na(location_m))

# drop NA start rows for veg; inserted for transects with 0 intersections
seki.veg = 
  seki.veg %>%
  filter(!is.na(start_m))

# drop NA status rows for trees; inserted for transects with 0 trees
seki.trees = 
  seki.trees %>%
  filter(!is.na(status)&status!='NA')


summary(seki.fwd_diams)

summary(seki.fwd_tallies)

summary(seki.litterduff)

summary(seki.trees)

summary(seki.veg)

seki.veg %>% filter(start_m > 15 | end_m > 15)


# 11-180 14.8 20.0 is a typo; fix it
seki.veg %>% filter(plot == '11' & az == '180')
seki.veg[seki.veg$plot=='11'&
           seki.veg$az=='180'&
           seki.veg$start_m==14.8&
           seki.veg$end_m==20.0,'end_m'] = 14.95

# 35-0 11.5-17.2 is a typo; fix it
seki.veg[seki.veg$plot=='35'&
           seki.veg$az=='0'&
           seki.veg$start_m==11.5&
           seki.veg$end_m==17.2,'end_m'] = 13.2

# 16-90 13.3 15.1 not a typo; truncate
# 9-270 14.1 15.5 not a typo; truncate
# 3-270 start 12.9 end 15.2 not a typo; truncate
# 9-270 14.4 15.1 not a typo; truncate
seki.veg[seki.veg$end_m>15,'end_m'] = 15

# 1-180 15.5 15.7 not a typo; remove
seki.veg = seki.veg[!(seki.veg$start_m>=15&seki.veg$start_m>=15),]


summary(seki.cwd)

seki.cwd %>% filter(diam_cm < 7.6)

# the 7.5cm diam was transcribed correctly from the field sheet; crew probably 
# just recorded a 7.5cm piece on accident and it should be dropped

# the 5 cm diam is actually a 15 on the datasheet
seki.cwd[seki.cwd$plot=='12'&seki.cwd$location_m==9.95,'diam_cm'] = 15

seki.cwd = seki.cwd %>% filter(diam_cm >= 7.6)

#### spot checking #############################################################

# checking the quality of data transcription by pulling 10 transects 
# at random...

set.seed(110819)
seki.metadata %>%
  sample_n(10) %>%
  arrange(plot, az)

print_transect = 
  function(p, a){
    seki.cwd %>% 
      filter(plot == p & az == a) %>%
      print(n = Inf, width = Inf)
    
    seki.fwd_diams %>%
      filter(plot==p & az == a) %>%
      print(n = Inf, width = Inf)
    
    seki.fwd_tallies%>%
      filter(plot==p & az == a) %>%
      print(n = Inf, width = Inf)
    
    seki.litterduff %>%
      filter(plot==p & az == a) %>%
      print(n = Inf, width = Inf)
    
    seki.veg %>%
      filter(plot==p & az == a) %>%
      print(n = Inf, width = Inf)
    
    seki.trees %>%
      filter(plot==p) %>%
      print(n = Inf, width = Inf)
      
  }


print_transect('2', '180') # looks good

print_transect('3', '270') # unchecked from here

print_transect('9', '180')

print_transect('9', '270')

print_transect('10', '180')

print_transect('12', '90')

print_transect('14', '90')

print_transect('16', '270')

print_transect('17', '0')

print_transect('17', '270')


#### standardize codes #########################################################

# tree spp
unique(seki.trees$spp)[order(unique(seki.trees$spp))]
# TOCA is correct (nutmeg), CAD presumably is CADE
seki.trees[seki.trees$spp=='CAD','spp'] = 'CADE'

# tree status
unique(seki.trees$status)

# tree ts
unique(seki.trees$ts)[order(unique(seki.trees$ts))]

# tree decay
unique(seki.trees$decay)[order(unique(seki.trees$decay))]

# cwd decay
unique(seki.cwd$decay)[order(unique(seki.cwd$decay))]

# shrub spp
unique(seki.veg$spp)[order(unique(seki.veg$spp))]
# i'm betting CUCH is QUCH
seki.veg = 
  seki.veg %>%
  left_join(.,
            data.frame(spp_recorded = 
                         c('ABCO', 'ALRH', 'ARPA', 'CADE', 'CEIN', 'CHFO',
                           'COCA', 'COCO', 'CONU', 'CUCH', 'F', 'G', 'PILA', 
                           'PIPO', 'QUCH', 'QUKE', 'QUVA', 'RIBES SP', 
                           'RIRO', 'ROGY', 'RUBUS', 'SYMO'),
                       spp_actual = 
                         c('ABCO', 'ALRH', 'ARPA', 'CADE', 'CEIN', 'CHFO',
                           'COCO', 'COCO', 'CONU', 'QUCH', 'F', 'G', 'PILA',
                           'PIPO', 'QUCH', 'QUKE', 'QUVA', 'RIB-', 'RIRO',
                           'ROGY', 'RUB-', 'SYMO')) %>%
              mutate(spp_recorded = as.character(spp_recorded),
                     spp_actual = as.character(spp_actual)),
            by = c('spp' = 'spp_recorded')) %>%
  select(plot, az, start_m, end_m, spp = spp_actual, status, height_m, comments)


# shrub status
unique(seki.veg$status)

# shrub height (should be in increments of 0.25)
unique(seki.veg$height_m)[order(unique(seki.veg$height_m))]

# plot, az for all tables
seki.cwd = 
  seki.cwd %>%
  mutate(plot = str_pad(plot,
                        width = 2,
                        side = 'left',
                        pad = '0'),
         transect_id = 
           paste0(plot, '-', str_pad(az, width = 3, side = 'left', pad = '0')),
         az = 
           as.numeric(az)) %>%
  select(plot_id = plot, 
         transect_id,
         az,
         location_m, diam_cm, decay_class = decay, comments)

seki.fwd_diams = 
  seki.fwd_diams %>%
  mutate(plot = str_pad(plot,
                        width = 2,
                        side = 'left',
                        pad = '0'),
         transect_id = 
           paste0(plot, '-', str_pad(az, width = 3, side = 'left', pad = '0')),
         az = 
           as.numeric(az)) %>%
  select(plot_id = plot, 
         transect_id,
         az,
         diam_mm)

seki.fwd_tallies = 
  seki.fwd_tallies %>%
  mutate(plot = str_pad(plot,
                        width = 2,
                        side = 'left',
                        pad = '0'),
         transect_id = 
           paste0(plot, '-', str_pad(az, width = 3, side = 'left', pad = '0')),
         az = 
           as.numeric(az)) %>%
  select(plot_id = plot, 
         transect_id,
         az,
         location_m, 
         a1h, a10h, a100h, b1h, b10h, b100h, comments)

seki.litterduff = 
  seki.litterduff %>%
  mutate(plot = str_pad(plot,
                        width = 2,
                        side = 'left',
                        pad = '0'),
         transect_id = 
           paste0(plot, '-', str_pad(az, width = 3, side = 'left', pad = '0')),
         az = 
           as.numeric(az)) %>%
  select(plot_id = plot, 
         transect_id,
         az,
         location_m, 
         litter_cm, duff_cm, fuel_cm, comments)

seki.trees = 
  seki.trees %>%
  mutate(plot = str_pad(plot,
                        width = 2,
                        side = 'left',
                        pad = '0'),
         tree_id = 
           paste0(plot,
                  '-',
                  ts,
                  '-',
                  location_m,
                  '-',
                  distance_m)) %>%
  select(plot_id = plot, 
         tree_id,
         ts, location_m, distance_m,
         status, spp, dbh_cm, height_m, htcb_m, decay_class = decay, comments)

seki.veg = 
  seki.veg %>%
  mutate(plot = str_pad(plot,
                        width = 2,
                        side = 'left',
                        pad = '0'),
         transect_id = 
           paste0(plot, '-', str_pad(az, width = 3, side = 'left', pad = '0')),
         az = 
           as.numeric(az)) %>%
  select(plot_id = plot, 
         transect_id,
         az,
         start_m, end_m, spp, status, height_m, comments)

seki.metadata

seki.metadata = 
  seki.metadata %>%
  mutate(plot_id = 
           str_pad(plot, width = 2, side = 'left', pad = '0'),
         transect_id = 
           paste0(plot_id,
                  '-',
                  str_pad(az, width = 3, side = 'left', pad = '0')),
         az = as.numeric(az)) %>%
  select(plot_id, transect_id, az, slp_degrees, cwd, fwd, litter_duff, veg, trees,
         inv_date, comment, comment2, comment3)

#### round all the shrub location codes ########################################

# data were inconsistently rounded to the nearest 0.1m in both the field datasheets
# and on digitization, so just round everything to nearest 0.1
seki.veg = 
  seki.veg %>%
  mutate(start_m = round(start_m, digits = 1),
         end_m = round(end_m, digits = 1))

#### check location:distance:dbh combos ########################################

# there shouldnt be any small (dbh < 11.4dm) trees more than 1m distance or 15m 
# location
seki.trees %>%
  filter(dbh_cm < 11.4 & (distance_m>1 | location_m > 15))
# info is all correct; should just have been not included so drop it
seki.trees = 
  seki.trees %>%
  filter(tree_id!='12-ES-6-2')

# there shouldnt be any big trees more than 15m location
seki.trees %>%
  filter(location_m > 15)
# typo; should be 10.6
seki.trees[seki.trees$tree_id=='09-ES-19.6-2.3','location_m'] = 
  10.6
seki.trees[seki.trees$tree_id=='09-ES-19.6-2.6', 'tree_id'] = 
  '09-ES-10.6-2.6'

# there can be big trees with dist > 5m, but they must have dist < 15m and 
# location < 5m
seki.trees %>%
  filter(distance_m > 15 | 
           (distance_m > 5 & location_m > 5))

#### check DBH:height and height:htcb ##########################################

ggplot(data = seki.trees,
       aes(x = dbh_cm, y = height_m, col = status))+
  geom_point()

ggplot(data = seki.trees %>% filter(status=='L'),
       aes(x = htcb_m, y = height_m))+
  geom_point()+
  coord_fixed()+
  geom_abline(intercept = 0, slope = 1, color = 'red')

seki.trees %>%
  filter(!(htcb_m<height_m))
# matches the datasheet, but looks like a mistake to me. set to NA and 
# interpolate

seki.trees[seki.trees$tree_id=='14-NE-5.7-4.1',
           'htcb_m'] = 
  predict(htcb_fit, newdata = seki.trees %>% filter(tree_id=='14-NE-5.7-4.1'))


#### check outliers ############################################################

ggplot(data = seki.cwd,
       aes(y = diam_cm))+
  geom_boxplot()

ggplot(data = seki.fwd_diams,
       aes(x = diam_mm))+
  geom_histogram()

ggplot(data = 
         seki.fwd_tallies %>%
         pivot_longer(cols = c(a1h, b1h, a10h, b10h, a100h, b100h),
                      names_to = 'timelag', values_to = 'count'),
       aes(x = count))+
  geom_histogram()+
  facet_wrap(~timelag, scales = 'free_y')

seki.litterduff %>%
  pivot_longer(cols = c(litter_cm, duff_cm, fuel_cm),
               names_to = 'class', values_to = 'depth_cm') %>%
  ggplot(aes(x = depth_cm))+
  geom_histogram()+
  facet_wrap(~class, scales = 'free')

ggplot(seki.metadata,
       aes(x = slp_degrees))+
  geom_histogram()

ggplot(seki.trees,
       aes(x = dbh_cm))+
  geom_histogram()+
  facet_grid(spp~status)

ggplot(seki.trees,
       aes(x = height_m))+
  geom_histogram()+
  facet_grid(spp~status)

ggplot(seki.trees %>% filter(status=='L'),
       aes(x = dbh_cm, y= height_m, color = spp))+
  geom_point()

ggplot(seki.veg,
       aes(x = start_m))+
  geom_histogram()

ggplot(seki.veg,
       aes(x = end_m))+
  geom_histogram()

ggplot(seki.veg,
       aes(x = height_m))+
  geom_histogram()


ggplot(seki.veg %>%
         mutate(cover_m = abs(start_m - end_m)) %>%
         group_by(plot_id, transect_id, spp) %>%
         summarise(cover_m = sum(cover_m)) %>%
         ungroup() %>%
         mutate(cover_perc = cover_m / 30),
       aes(x = spp, y = cover_perc))+
  geom_boxplot()



# looks good; i really need to know what coco is
# Corylus cornuta! https://calscape.org/Corylus-cornuta-ssp.-californica-()

#### check metadata ############################################################

# check for any weirdness noted in the metadata
seki.metadata %>%
  print(n = Inf, width = Inf)

seki.veg %>% filter(transect_id=='08-180')

# went back and got that shrubs data; set it to 1
seki.metadata[seki.metadata$plot_id=='08'&
                seki.metadata$az==180,
              'veg'] = '1'


#### data display figure #######################################################

library(cowplot)
litterduff_sample = 
  ggplot(data = 
         seki.litterduff %>%
         filter(plot_id=='14'&
                  is.element(transect_id, c('14-090', '14-270'))) %>%
         pivot_longer(cols = c(litter_cm, duff_cm),
                      names_to = c('component', 'units'),
                      names_sep = '_', 
                      values_to = 'depth_cm') %>%
         mutate(component = factor(str_to_title(as.character(component)),
                                   levels = c('Litter', 'Duff')),
                x_position = sin(az*(pi/180))*location_m),
       aes(x = x_position, y = depth_cm, fill = component))+
  geom_col(position = position_stack())+
  scale_x_continuous(limits = c(-30, 30))+
  theme_bw()+
  labs(x = 'Along transect (m)',
       y = 'Depth (cm)',
       fill = 'Component',
       title = 'Litter and Duff')
litterduff_sample

cwd_sample = 
  ggplot(data = seki.cwd %>%
         filter(is.element(transect_id, c('14-090', '14-270'))) %>%
         mutate(x_relative = sin(az*(pi/180))*location_m,
                decay_class = factor(decay_class,
                                     levels = c('1', '2', '3', '4', '5'))),
       aes(y = plot_id, x = x_relative, color = decay_class, size = diam_cm))+
  scale_x_continuous(limits = c(-30, 30))+
  geom_point()+
  #scale_size_identity()+
  theme_bw()+
  scale_color_viridis_d(option = 'C',
                        begin = 0.05, end = 0.85)+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank())+
  labs(x = 'Along transect (m)',
       size = 'Diameter (cm)',
       color = 'Decay class',
       title = 'Coarse woody debris')+
  theme(legend.direction = 'horizontal',
        legend.position = 'bottom',
        legend.box = 'vertical',
        legend.spacing = unit(0.1, units = 'inches'),
        legend.margin = margin(unit(0.1, units = 'inches')))

cwd_sample

fwddiams_sample = 
  ggplot(data = 
         seki.fwd_diams %>%
         filter(is.element(transect_id, c('14-090', '14-270'))),
       aes(x = diam_mm))+
  geom_histogram()+
  theme_bw()+
  labs(x = 'Diameter (mm)', y = 'Frequency', 
       title = 'FWD Diameters')

fwddiams_sample

fwdtallies_sample = 
  ggplot(data = seki.fwd_tallies %>%
         filter(is.element(transect_id, c('14-090', '14-270'))) %>%
         dplyr::select(az, location_m, a1h, a10h, a100h, 
                b1h, b10h, b100h) %>%
         pivot_longer(cols = c(a1h, a10h, a100h, b1h, b10h, b100h),
                      names_to = 'measure', values_to = 'count') %>%
         mutate(
                timelag_class = gsub(x = measure, pattern = 'a|b', replacement = ''),
                subsample = gsub(x = measure, pattern = '1h|10h|100h', replacement = ''),
                x_relative = sin(az*(pi/180))*location_m,
                timelag_class = factor(timelag_class,
                                       levels = c('100h', '10h', '1h'))) %>%
         filter(subsample=='a'),
       aes(x = x_relative, y = count, fill = timelag_class))+
  geom_col()+
  theme_bw()+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'A')+
  labs(x = 'Along transect (m)', y = 'Count (per meter)', fill = 'Timelag Class',
       title = 'FWD Tallies')

fwdtallies_sample




trees_sample = 
  ggplot(data = 
         seki.trees %>%
         filter(plot_id=='14') %>%
        mutate(x_coord = ifelse(is.element(ts, c('NE', 'SE')),
                         distance_m,
                          ifelse(is.element(ts, c('NW', 'SW')),
                                 -1*distance_m,
                                 ifelse(is.element(ts, c('EN', 'ES')),
                                        location_m,
                                        ifelse(is.element(ts, c('WN', 'WS')),
                                               -1*location_m,
                                               NA)))),
               y_coord = ifelse(is.element(ts, c('NE', 'NW')),
                          location_m,
                          ifelse(is.element(ts, c('SE', 'SW')),
                                 -1*location_m,
                                 ifelse(is.element(ts, c('EN', 'WN')),
                                        distance_m,
                                        ifelse(is.element(ts, c('ES', 'WS')),
                                               -1*distance_m,
                                               NA))))),
       aes(x = x_coord, y = y_coord, size = dbh_cm, color = height_m))+
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
  theme(legend.box = 'vertical',
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.spacing = unit(0.1, units = 'inches'),
        legend.margin = margin(unit(0.1, units = 'inches')))

trees_sample

veg_sample = 
  ggplot(data = 
         seki.veg %>%
         filter(is.element(transect_id, c('14-090', '14-270'))) %>%
         mutate(start_relative = sin(az*(pi/180))*start_m,
                end_relative = sin(az*(pi/180))*end_m),
       aes(xmin = start_relative, xmax = end_relative, ymax = height_m, fill = spp))+
  geom_rect(aes(xmin = start_relative, xmax = end_relative,
                ymax = height_m, fill = spp), ymin = 0, alpha = 0.75)+
  theme_bw()+
  scale_x_continuous(limits = c(-15, 15))+
  scale_y_continuous(limits = c(0, 2))+
  scale_fill_brewer(type = 'qual', palette = 'Set1')+
  labs(x = 'Along transect (m)', y = 'Above ground (m)', fill = 'Species',
       title = 'Understory Veg')
veg_sample

all_together = 
    
    plot_grid(
    
      plot_grid(trees_sample, veg_sample,
              nrow = 2, ncol = 1,
              rel_heights = c(1.75, 1)),
      plot_grid(litterduff_sample, 
                plot_grid(fwddiams_sample,
                          fwdtallies_sample,
                          nrow = 1, ncol = 2,
                          rel_widths = c(1, 2)), 
                cwd_sample,
              nrow = 3, ncol = 1,
              rel_heights = c(1, 1, 1)),
      ncol = 2, nrow = 1,
      rel_widths = c(1, 1.5))

all_together

ggsave(all_together,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'data_sample.png'),
       height = 6.5, width = 9, units = 'in')

#### write results #############################################################

write.csv(seki.cwd,
          here::here('02-data',
                     '03-clean',
                     'seki_cwd.csv'),
          row.names = FALSE)

saveRDS(seki.cwd,
        here::here('02-data',
                   '03-clean',
                   'seki_cwd.rds'))

write.csv(seki.fwd_diams,
          here::here('02-data',
                     '03-clean',
                     'seki_fwddiams.csv'),
          row.names = FALSE)

saveRDS(seki.fwd_diams,
        here::here('02-data',
                   '03-clean',
                   'seki_fwddiams.rds'))

write.csv(seki.fwd_tallies,
          here::here('02-data',
                     '03-clean',
                     'seki_fwdtallies.csv'),
          row.names = FALSE)

saveRDS(seki.fwd_tallies,
        here::here('02-data',
                   '03-clean',
                   'seki_fwdtallies.rds'))

write.csv(seki.litterduff,
          here::here('02-data',
                     '03-clean',
                     'seki_litterduff.csv'),
          row.names = FALSE)

saveRDS(seki.litterduff,
        here::here('02-data',
                   '03-clean',
                   'seki_litterduff.rds'))

write.csv(seki.metadata,
          here::here('02-data',
                     '03-clean',
                     'seki_metadata.csv'),
          row.names = FALSE)

saveRDS(seki.metadata,
        here::here('02-data',
                   '03-clean',
                   'seki_metadata.rds'))

write.csv(seki.trees,
          here::here('02-data',
                     '03-clean',
                     'seki_trees.csv'),
          row.names = FALSE)

saveRDS(seki.trees,
        here::here('02-data',
                   '03-clean',
                   'seki_trees.rds'))

write.csv(seki.veg,
          here::here('02-data',
                     '03-clean',
                     'seki_veg.csv'),
          row.names = FALSE)

saveRDS(seki.veg,
        here::here('02-data',
                   '03-clean',
                   'seki_veg.rds'))

