# Goal: Figure out the optimal number of subsamples in each plot to estimate 
# parameters in the litter/duff and FWD models.

# OPTION 2: ALTER NUMBER OF SUBSAMPELES 

# Approach:
# - Fit FWD models with Danny's SEKI data.
# - alter the number of transects or subsamples within a plot
# - Fit models for each of those iterations, see how it affects the parameters

rm(list = ls())
library(tidyverse)
library(rethinking)
library(cowplot)
theme_set(theme_classic())


# Setup -------------------------------------------------------------------

# load data
FWDdf <- read_rds('02-data/03-clean/seki_fwdtallies_alternative.rds')


# get rid of set b tallies and NAs, make plot_id consecutive
FWDdf <- FWDdf %>% filter(!is.na(h1), set != 'b') %>% 
  mutate(plot_id = as.integer(as.factor(plot_id)))


# MODEL
# Y_xp ~ NB(mu_xp, phi)
# log(mu_xp) = a0[groupID] + ap + GP_xp
# ap ~ N(0, sigma_p)
# GP ~ MVN(0, K(d))
# [K(d)]_ij = eta^2 * exp(-.5/rho^2 * dist_ij^2)




# CHOOSE HOW TO SUBSET LOCATIONS ------------------------------------------

# Create coordinates of fuel samples, different designs
getLocations <- function(locations, printPlots = T){
  x = locations
  x <- c(rev(x*-1), x)
  y <- x
  coords <- rbind(
    cbind(x, y = 0),
    cbind(x = 0, y = y)
  )
  coordPlot <- ggplot(as.data.frame(coords), aes(x, y)) + geom_point() + coord_equal()
  
  # look at pairwise distances
  distmatLower <- fields::rdist(coords, compact = T)
  histPlot <- ggplot(data.frame(x = as.vector(distmatLower)), aes(x)) +
    geom_histogram(fill = 'grey', color = 'black', bins = 30)
  plots <- plot_grid(coordPlot, histPlot)
  
  if(printPlots == F) return(coords)
  else return(list(plots = plots, coords = coords))
  
}

# Locations. the number refers to number of samples.
x8 <- c(2.5, 3.5, 4.5, 5.5, 7.5, 9.5, 19.5, 29.5) #current design
x7 <- c(2.5, 3.5, 4.5, 5.5, 7.5, 9.5, 19.5)
x6 <- c(2.5, 3.5, 5.5, 7.5, 9.5, 19.5)
x5 <- c(2.5, 3.5, 5.5, 9.5, 19.5)
getLocations(x8)
getLocations(x7)
getLocations(x6)
getLocations(x5)



# Define function to subset the datasets with a smaller set of locations
subsetDat <- function(locations){
  locations2 <- c(rev(locations*-1), locations)
  datNew <- FWDdf %>% 
    filter(x_rel %in% locations2 | y_rel %in% locations2) %>% 
    mutate(location_id = as.integer(as.factor(location_id)))
  coordsdf <- datNew %>% 
    distinct(x_rel, y_rel, location_id) %>% 
    arrange(location_id)
  coords <- cbind(coordsdf$x_rel, coordsdf$y_rel)
  return(list(datNew, coords))
}



# DEFINE DATA -------------------------------------------------------------


# generate datalist of potentially subsetted data
getDataList_Loc <- function(keep_loc, response_var){
  
  # susample the data
  datlist <- subsetDat(keep_loc)
  dat <- datlist[[1]]
  coords <- datlist[[2]]
  
  
  # define datalist for Stan
  response = pull(dat, all_of(response_var))
  list(
    N = nrow(dat),
    L = nrow(coords),
    P = max(dat$plot_id),
    G = max(dat$group_id),
    coords = coords,
    plotID = dat$plot_id,
    locID = dat$location_id,
    groupID = dat$group_id,
    y = as.integer(response)
  )
}



# RUN STAN MODELS ---------------------------------------------------------

stan_model.fwd = cmdstan_model('04-sensitivity_analyses/GP_NB.stan')

#fit model
fitmodel <- function(keep_loc, responseVar, nchains = 3, iter = 2000){
  df_sample <- getDataList_Loc(keep_loc, responseVar)
  fit <- stan_model.fwd$sample(data = df_sample, 
                               adapt_delta = .95, 
                               #init = list(inits1, inits2, inits3),
                               iter_warmup = iter/2, iter_sampling = iter/2,
                               parallel_chains = nchains, chains = nchains)
  return(fit)
}

# locations list
locList <- list(x8, x7, x6, x5)
nsubsamples <- c(8, 7, 6, 5)

# fit models
fits10h <- lapply(locList, function(x) fitmodel(x, 'h10'))
fits1h <- lapply(locList, function(x) fitmodel(x, 'h1'))

#saveRDS(fits10h, '02-data/06-results/sensitivity/fits/fits10h_nsubsamples.rds')
#saveRDS(fits1h, '02-data/06-results/sensitivity/fits/fits1h_nsubsamples.rds')

# SUMMARIZE POSTERIORS ----------------------------------------------------


# 10-hr fuels ================================

# combine draws 
draws10h <- list(NULL)
for(i in 1:length(fits10h)){
  draws <- fits10h[[i]]$draws(
    variables = c('a0', 'alpha', 'rho', 'sigma_ap', 'kappa'), 
    format = 'draws_df')
  draws10h[[i]] <- draws %>% 
    mutate(nSubsamples = nsubsamples[i])
}
draws10h <- bind_rows(draws10h)

# summarize posteriors, grouped by nSubsamples
summary10h <- draws10h %>% 
  select(-c(.chain, .iteration, .draw)) %>% 
  pivot_longer(-c(nSubsamples)) %>% 
  group_by(name, nSubsamples) %>% 
  summarize(median = median(value),
            SD = sd(value),
            lower90 = HPDI(value, .9)[1],
            upper90 = HPDI(value, .9)[2],
            lower50 = HPDI(value, .5)[1],
            upper50 = HPDI(value, .5)[2]) %>% 
  mutate(fuelSize = '10hr')


# 1-hr fuels ================================


# combine draws 
draws1h <- list(NULL)
for(i in 1:length(fits1h)){
  draws <- fits1h[[i]]$draws(
    variables = c('a0', 'alpha', 'rho', 'sigma_ap', 'kappa'), 
    format = 'draws_df')
  draws1h[[i]] <- draws %>% 
    mutate(nSubsamples = nsubsamples[i])
}
draws1h <- bind_rows(draws1h)

# summarize posteriors, grouped by number of transects
summary1h <- draws1h %>% 
  select(-c(.chain, .iteration, .draw)) %>% 
  pivot_longer(-c(nSubsamples)) %>% 
  group_by(name, nSubsamples) %>% 
  summarize(median = median(value),
            SD = sd(value),
            lower90 = HPDI(value, .9)[1],
            upper90 = HPDI(value, .9)[2],
            lower50 = HPDI(value, .5)[1],
            upper50 = HPDI(value, .5)[2]) %>% 
  mutate(fuelSize = '1hr')



# VISUALIZE PARAM COMPARISON ----------------------------------------------

# viz posteriors 10-hr fuels
ggplot(summary10h, aes(nSubsamples, median)) +
  geom_linerange(aes(ymin = lower50, ymax = upper50), 
                 alpha = 1, size = 1.5, color = '#dd9933') +
  geom_pointrange(aes(ymin = lower90, ymax = upper90), 
                  alpha = .4, size = 1, fatten = 2) +
  facet_wrap(~name, scales = 'free_y', nrow = 2) +
  labs(x = 'Number of Subsamples', y = 'median, 90% and 50% HDPI',
       title = 'Tallies of 10-hr fuels, Dannys Crystal Cave data', 
       subtitle = 'altered number of subsamples.')
ggsave('02-data/06-results/sensitivity/FWD10h_nsubsamples_realdat.png', width = 7, height = 3.5)


# viz posteriors 1-hr fuels
ggplot(summary1h, aes(nSubsamples, median)) +
  geom_linerange(aes(ymin = lower50, ymax = upper50), 
                 alpha = 1, size = 1.5, color = '#05878a') +
  geom_pointrange(aes(ymin = lower90, ymax = upper90), 
                  alpha = .4, size = 1, fatten = 2) +
  facet_wrap(~name, scales = 'free_y', nrow = 2) +
  labs(x = 'Number of Subsamples', y = 'median, 90% and 50% HDPI',
       title = 'Tallies of 1-hr fuels, Dannys Crystal Cave data', 
       subtitle = 'altered number of subsamples.')
ggsave('02-data/06-results/sensitivity/FWD1h_nsubsamples_realdat.png', width = 7, height = 3.5)

# viz 1 and 10 together
bind_rows(summary10h, summary1h) %>% 
  ggplot(., aes(nSubsamples, median, group = fuelSize)) +
  geom_linerange(aes(ymin = lower50, ymax = upper50, color = fuelSize), 
                 position = position_dodge2(width = .5),
                 alpha = 1, size = 1.5) +
  geom_pointrange(aes(ymin = lower90, ymax = upper90), 
                  position = position_dodge2(width = .5),
                  alpha = .4, size = 1, fatten = 1.5) +
  scale_colour_manual(values = c('#dd9933', '#05878a')) +
  facet_wrap(~name, scales = 'free_y', nrow = 2) +
  labs(x = 'Number of transects', y = 'median, 90% and 50% HDPI',
       title = 'Tallies of 1 & 10-hr fuels, Dannys Crystal Cave data')
ggsave('02-data/06-results/sensitivity/FWD1_10h_nsubsamples_realdat.png', width = 7, height = 3.5)
