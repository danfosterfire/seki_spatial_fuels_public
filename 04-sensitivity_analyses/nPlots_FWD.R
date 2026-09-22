# Goal: figure out how the number of plots changes the precision of the parameter
# estimates. Do for the FWD and litter/duff since those models are written already.

# Approach: 
# - fit the models, monitor the distribution of the parameters
# - for mortality groups, include a fixed effect but don't bother varying alpha or rho
# - repeat many times with a differently sized dataset, derived by resampling



rm(list = ls())
library(tidyverse)
library(rethinking)
library(cowplot)
library(cmdstanr)
theme_set(theme_classic())


# load data
dataList_fwd1h <- read_rds('02-data/05-for_analysis/fwd1h_data.rds')
str(dataList_fwd1h)
coords <- dataList_fwd1h$coords
FWDdf_raw <- read_rds('02-data/04-geolocated/seki_fwdtallies.rds')

# manipulate the data a little
FWDdf <- FWDdf_raw %>% 
  #column for a and b subsamples
  pivot_longer(
    cols = a1h:b100h,
    names_to = c('set', ".value"),
    names_pattern = "(.)(.*)", names_prefix = 'h'
  ) %>% 
  rename(h1 = `1h`, h10 = `10h`, h100 = `100h`) %>% 
  # remove NAs 
  filter(!is.na(h1)) %>% 
  # add relative locations
  mutate(
    x_rel = case_when(
      az == 90 ~ location_m,
      az == 270 ~ -location_m,
      T ~ 0),
    y_rel = case_when(
      az == 0 ~ location_m,
      az == 180 ~ -location_m,
      T ~ 0
    )) %>% 
  mutate(mort = factor(mort, levels = c('low', 'mid', 'high'))) %>% 
  # add ID codes
  left_join(as.data.frame(coords) %>% mutate(location_id = 1:nrow(.))) %>% 
  mutate(plot_id = as.integer(plot_id),
         group_id = as.integer(mort)) 
print(FWDdf, width = Inf)
saveRDS(FWDdf, '02-data/03-clean/seki_fwdtallies_alternative.rds')

# actually, I don't want set b. remove that too.
FWDdf <- filter(FWDdf, set == 'a')


# DEFINE RESAMPLING AND DATALIST FUNCTIONS --------------------------------


# resample the dataset with replacement. Balance by mortality class?
mortClasses = distinct(FWDdf, plot_id, mort)
plotIDs <- split(mortClasses$plot_id, mortClasses$mort)

resampleData <- function(nPerMort){
  choosePlots <- sapply(plotIDs, function(x) sample(x, nPerMort, replace = T)) %>% 
    as.vector()
  df <- lapply(1:length(choosePlots), function(i) 
    FWDdf[FWDdf$plot_id == choosePlots[i],] %>% 
      mutate(plot_id = i)) %>% #same plots can now have different plot_ids
    bind_rows()
  # make sure plotid and groupid's are consecutive, make mort a character string
  #df <- df %>% 
  #  mutate(#plot_id = sample, #as.integer(as.factor(plot_id)),
  #         group_id = as.integer(as.factor(group_id)))
  return(df)
}


# prep datalist
getDataList <- function(dat, response_var){
  response = pull(dat, all_of(response_var))
  list(
    N = nrow(dat),
    L = max(dat$location_id),
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


PosteriorFWD <- function(plotsPerMort, nsims, responseVar){
  
  # define list to dump the posteriors into and start values for kappa and alpha
  post <- list(NULL)
  inits1 <- list(kappa = .1, alpha = 1)
  inits2 <- list(kappa = 4, alpha = .1)
  inits3 <- list(kappa = .5, alpha = .5)
  
  # fit resampled data several times and get summary of posteriors
  for(s in 1:nsims){
    #fit
    df_sample <- getDataList(resampleData(plotsPerMort), responseVar)
    fit <- stan_model.fwd$sample(data = df_sample, 
                                 adapt_delta = .95, 
                                 init = list(inits1, inits2, inits3),
                                 iter_warmup = 1000, iter_sampling = 1000,
                                 parallel_chains = 3, chains = 3)
    
    #extract post & summarize
    draws <- fit$draws(
      variables = c('a0', 'alpha', 'rho', 'sigma_ap', 'kappa'), 
      format = 'draws_df')
    post[[s]] <- draws %>% 
      select(-c(.chain, .iteration, .draw)) %>% 
      pivot_longer(everything()) %>% 
      group_by(name) %>% 
      summarize(median = median(value),
                SD = sd(value),
                lower90 = HPDI(value, .9)[1],
                upper90 = HPDI(value, .9)[2],
                lower50 = HPDI(value, .5)[1],
                upper50 = HPDI(value, .5)[2]) %>% 
      mutate(nPlots = plotsPerMort*3,
             simulation = s, 
             fuelSize = responseVar)
  }
  
  # bind posteriors together
  post <- bind_rows(post)
  return(post)
}

# fit lots of resampled data across a gradient of sample sizes
plotsTotal <- c(3,6, 9, 12, 15, 21, 24, 30)
plotsPerMort_vector <- plotsTotal / 3

# 1- hr fuels
FWD_10h <- lapply(1:length(plotsPerMort_vector), 
                  function(x) PosteriorFWD(plotsPerMort_vector[x], 5, 'h10'))
FWD_10h <- bind_rows(FWD_10h)

# 1-hr fuels
FWD_1h <- lapply(1:length(plotsPerMort_vector), 
                 function(x) PosteriorFWD(plotsPerMort_vector[x], 5, 'h1'))
FWD_1h <- bind_rows(FWD_1h)


# VISUALIZE COMPARISON OF POSTERIORS --------------------------------------

# viz posteriors 10-hr fuels
ggplot(FWD_10h, aes(nPlots, median)) +
  geom_pointrange(aes(ymin = lower90, ymax = upper90), 
                  position = position_dodge2(width = 1),
                  alpha = .7, size = .2) +
  facet_wrap(~name, scales = 'free_y') +
  labs(x = 'Total number of plots', y = 'median, 90% HDPI',
       title = 'Tallies of 10-hr fuels, Dannys Crystal Cave data', 
       subtitle = '5 sets of resampled data per sample size')
ggsave('02-data/06-results/sensitivity/FWD10h_nplots.png', width = 7, height = 5)


# viz posteriors 1-hr fuels
ggplot(FWD_1h, aes(nPlots, median)) +
  geom_pointrange(aes(ymin = lower90, ymax = upper90), 
                  position = position_dodge2(width = 1),
                  alpha = .7, size = .2) +
  facet_wrap(~name, scales = 'free_y') +
  labs(x = 'Total number of plots', y = 'median, 90% HDPI',
       title = 'Tallies of 1-hr fuels, Dannys Crystal Cave data', 
       subtitle = '5 sets of resampled data per sample size')
ggsave('02-data/06-results/sensitivity/FWD1h_nplots.png', width = 7, height = 5)


# viz them togehter
bind_rows(FWD_1h, FWD_10h) %>% 
  ggplot(., aes(nPlots, median, color = fuelSize)) +
  geom_pointrange(aes(ymin = lower90, ymax = upper90), 
                  position = position_dodge2(width = 1),
                  alpha = .7, size = .2) +
  facet_wrap(~name, scales = 'free_y') +
  scale_colour_manual(values = c('#05878a', '#dd9933')) +
  labs(x = 'Total number of plots', y = 'median, 90% HDPI',
       title = 'Tallies of 1- & 10-hr fuels, Dannys Crystal Cave data', 
       subtitle = '5 sets of resampled data per sample size')
ggsave('02-data/06-results/sensitivity/FWD1_10h_nplots.png', width = 8, height = 5)


# SINGLE TEST RUN AND EVALUATION ------------------------------------------

# 1 hr fuels ==============================
df1hr <- getDataList(FWDdf, 'h1')
fit_1hr = stan_model.fwd$sample(data = df1hr, parallel_chains = 2, chains = 2)
#fit_test$cmdstan_diagnose()
#fit_test$summary(c('a0', 'alpha', 'rho', 'sigma_ap', 'kappa')) %>% print(n = Inf)

# goodness of fit
y_rep <- fit_1hr$draws(variables = c('y_rep'), format = 'draws_matrix')
bayesplot::ppc_dens_overlay(df1hr$y, y_rep[1:50,])

# look at covariance/correlation function
simCov1hr <- fit_1hr$draws(variables = 'K', format = 'draws_matrix')
simCovList1hr <- lapply(1:50, function(x) matrix(simCov1hr[x,], nrow = sqrt(dim(simCov1hr)[2]), byrow = T))
simCorList1hr <- lapply(simCovList1hr, cov2cor)#convert cov to cor matrix

# put into dataframe to plot correlation vs distance
plotCorCov <- function(CorCovList, ylab, Title = 'GP posterior', xmax = 20){
  dmat <- as.matrix(dist(coords, diag = T, upper = T))
  tmp <- lapply(1:length(CorCovList), function(z) 
    cbind(draw = z, 
          cor = as.vector(CorCovList[[z]]), 
          d = as.vector(dmat))
  )
  tmp2 <- as.data.frame(do.call(rbind, tmp))
  p <- tmp2 %>%
    ggplot(., aes(d, cor)) +
    geom_line(lwd = .2, aes(group = draw), alpha = .2) +
    coord_cartesian(xlim = c(0,xmax)) +
    labs(x = 'Distance (m)', y = ylab, title = Title)
  return(list(tmp2, p))
}

#
cov1h <- plotCorCov(simCovList1hr, 'Covariance', '1-hr GP posterior')
cor1h <- plotCorCov(simCorList1hr, 'Correlation')



# 10 hr fuels ============================

df10hr <- getDataList(FWDdf, 'h10')
fit_10hr = stan_model.fwd$sample(data = df10hr, parallel_chains = 2, chains = 2)

# look at covariance/correlation function
simCov10hr <- fit_10hr$draws(variables = 'K', format = 'draws_matrix')
simCovList10hr <- lapply(1:50, function(x) matrix(simCov10hr[x,], nrow = sqrt(dim(simCov10hr)[2]), byrow = T))
simCorList10hr <- lapply(simCovList10hr, cov2cor)#convert cov to cor matrix

# put into dataframe to plot correlation vs distance
cov10h <- plotCorCov(simCovList10hr, 'Covariance', '10-hr GP posterior')
cor10h <- plotCorCov(simCorList10hr, 'Correlation', '10-hr GP posterior')


# plot them together?
binded_dat <- bind_rows(
  cov1h[[1]] %>% mutate(fuelSize = '1hr'), 
  cov10h[[1]] %>% mutate(fuelSize = '10hr'))

ggplot(filter(binded_dat, draw %in% 1:50), aes(d, cor, color = fuelSize)) +
  geom_line(lwd = .3, aes(group = paste0(draw, fuelSize)), alpha = .5) +
  coord_cartesian(xlim = c(0,30)) +
  scale_colour_manual(values = c('#dd9933', '#05878a')) +
  labs(x = 'Distance (m)', y = 'covariance', title = '1 and 10 hr fuels GP posterior',
       subtitle = 'mortality included as fixed effect, alpha or rho dont vary by it')
ggsave('02-data/06-results/sensitivity/FWD1_10h_covariance.png', width = 6, height = 3.5)
