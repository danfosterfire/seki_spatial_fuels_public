# Goal: Figure out the optimal number of subsamples in each plot to estimate 
# parameters in the litter/duff and FWD models.

# Approach:
# - Fit FWD models with Danny's SEKI data.
# - alter the number of transects or subsamples within a plot
# - Fit models for each of those iterations, see how it affects the parameters

rm(list = ls())
library(tidyverse)
library(rethinking)
library(cmdstanr)
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



# OPTION 1: ALTER NUMBER OF TRANSECTS -------------------------------------



# Define function to subset the datasets with a smaller set of locations
subsetDat <- function(keep_az){
  
  # susample the data
  dat <- FWDdf %>% 
    filter(az %in% keep_az) %>% 
    mutate(location_id = as.integer(as.factor(location_id)))
  
  # get coords
  coordsdf <- dat %>% 
    distinct(x_rel, y_rel, location_id) %>% 
    arrange(location_id)
  coords <- cbind(coordsdf$x_rel, coordsdf$y_rel)
  
  return(list(dat, coords))
}


# generate datalist of potentially subsetted data
getDataList_Az <- function(keep_az, response_var){
  
  # susample the data
  datlist <- subsetDat(keep_az)
  dat <- datlist[[1]]
  coords = datlist[[2]]
  
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
fitmodel <- function(keep_az, responseVar, nchains = 3, iter = 2000){
  df_sample <- getDataList_Az(keep_az, responseVar)
  fit <- stan_model.fwd$sample(data = df_sample, 
                               adapt_delta = .95, 
                               #init = list(inits1, inits2, inits3),
                               iter_warmup = iter/2, iter_sampling = iter/2,
                               parallel_chains = nchains, chains = nchains)
  return(fit)
}

# include possible azimuth combos. designs with 2 and 3 transect will be joined 
# as ensembles, which creates a weighted avg
azList <- list(
  c(0, 90, 180, 270), 
  c(0, 90, 180),
  c(0, 90, 270),
  c(0, 180, 270), 
  c(90, 180, 270), 
  c(0, 180),
  c(90, 270)
)
azMetadat <- data.frame(
  nTransects = c(4, 3, 3, 3, 3, 2, 2),
  set = c(1, 1, 2, 3, 4, 1, 2)
)


fits10h <- lapply(azList, function(x) fitmodel(x, 'h10'))
fits1h <- lapply(azList, function(x) fitmodel(x, 'h1'))

#saveRDS(fits10h, '02-data/06-results/sensitivity/fits/fits10h_ntransects.rds')
#saveRDS(fits1h, '02-data/06-results/sensitivity/fits/fits1h_ntransects.rds')
#fits10h <- read_rds('02-data/06-results/sensitivity/fits/fits10h_ntransects.rds')
#fits1h <- read_rds('02-data/06-results/sensitivity/fits/fits1h_ntransects.rds')

# SUMMARIZE POSTERIORS ----------------------------------------------------


# 10-hr fuels ================================

# combine draws 
draws10h <- list(NULL)
for(i in 1:length(fits10h)){
  draws <- fits10h[[i]]$draws(
    variables = c('a0', 'alpha', 'rho', 'sigma_ap', 'kappa'), 
    format = 'draws_df')
  draws10h[[i]] <- draws %>% 
    mutate(nTransects = azMetadat$nTransects[i],
           set = azMetadat$set[i])
}
draws10h <- bind_rows(draws10h)

# summarize posteriors, grouped by number of transects
summary10h <- draws10h %>% 
  select(-c(.chain, .iteration, .draw)) %>% 
  pivot_longer(-c(nTransects, set)) %>% 
  group_by(name, nTransects) %>% 
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
    mutate(nTransects = azMetadat$nTransects[i],
           set = azMetadat$set[i])
}
draws1h <- bind_rows(draws1h)

# summarize posteriors, grouped by number of transects
summary1h <- draws1h %>% 
  select(-c(.chain, .iteration, .draw)) %>% 
  pivot_longer(-c(nTransects, set)) %>% 
  group_by(name, nTransects) %>% 
  summarize(median = median(value),
            SD = sd(value),
            lower90 = HPDI(value, .9)[1],
            upper90 = HPDI(value, .9)[2],
            lower50 = HPDI(value, .5)[1],
            upper50 = HPDI(value, .5)[2]) %>% 
  mutate(fuelSize = '1hr')



# VISUALIZE PARAM COMPARISON ----------------------------------------------

# viz posteriors 10-hr fuels
ggplot(summary10h, aes(nTransects, median)) +
  geom_linerange(aes(ymin = lower50, ymax = upper50), 
                  alpha = 1, size = 1.5, color = '#dd9933') +
  geom_pointrange(aes(ymin = lower90, ymax = upper90), 
                  alpha = .4, size = 1, fatten = 2) +
  facet_wrap(~name, scales = 'free_y', nrow = 2) +
  labs(x = 'Number of transects', y = 'median, 90% and 50% HDPI',
       title = 'Tallies of 10-hr fuels, Dannys Crystal Cave data', 
       subtitle = 'altered number of transects. data with 2 and 3 trans. ensembles of each combo.')
ggsave('02-data/06-results/sensitivity/FWD10h_ntransects.png', width = 7, height = 3.5)


# viz posteriors 1-hr fuels
ggplot(summary1h, aes(nTransects, median)) +
  geom_linerange(aes(ymin = lower50, ymax = upper50), 
                 alpha = 1, size = 1.5, color = '#05878a') +
  geom_pointrange(aes(ymin = lower90, ymax = upper90), 
                  alpha = .4, size = 1, fatten = 2) +
  facet_wrap(~name, scales = 'free_y', nrow = 2) +
  labs(x = 'Number of transects', y = 'median, 90% and 50% HDPI',
       title = 'Tallies of 1-hr fuels, Dannys Crystal Cave data', 
       subtitle = 'altered number of transects. data with 2 and 3 trans. ensembles of each combo.')
ggsave('02-data/06-results/sensitivity/FWD1h_ntransects.png', width = 7, height = 3.5)

# viz 1 and 10 together
bind_rows(summary10h, summary1h) %>% 
  ggplot(., aes(nTransects, median, group = fuelSize)) +
  geom_linerange(aes(ymin = lower50, ymax = upper50, color = fuelSize), 
                 position = position_dodge2(width = .5),
                 alpha = 1, size = 1.5) +
  geom_pointrange(aes(ymin = lower90, ymax = upper90), 
                  position = position_dodge2(width = .5),
                  alpha = .4, size = 1, fatten = 1.5) +
  scale_colour_manual(values = c('#dd9933', '#05878a')) +
  facet_wrap(~name, scales = 'free_y', nrow = 2) +
  labs(x = 'Number of transects', y = 'median, 90% and 50% HDPI',
       title = 'Tallies of 10-hr fuels, Dannys Crystal Cave data', 
       subtitle = 'altered number of transects. data with 2 and 3 trans. ensembles of each combo.')
ggsave('02-data/06-results/sensitivity/FWD1_10h_ntransects.png', width = 7, height = 3.5)


# CHECK OUT COVAR FUNCTION ------------------------------------------------

# need list of coords. each dataset has a unique coords list. 
coordsList <- lapply(azList, function(x) getDataList_Az(x, 'h1')$coords)


# extract Covariance or Correlatin matrix
extractCovCor <- function(fit, draws){
  
  # look at covariance/correlation function
  Cov <- fit$draws(variables = 'K', format = 'draws_matrix')
  totaldraws = dim(Cov)[1]
  selectdraws <- sample(totaldraws, draws)
  CovList <- lapply(1:length(selectdraws), 
                    function(x) matrix(Cov[selectdraws[x],], nrow = sqrt(dim(Cov)[2]), byrow = T))
  CorList <- lapply(CovList, cov2cor)#convert cov to cor matrix
  return(list(CovList = CovList, CorList = CorList))
}

# put into dataframe to plot correlation vs distance
dist_vs_covcor <- function(CorCovList, coords, ylab = 'cor/cov', Title = 'GP posterior', xmax = 20){
  dmat <- as.matrix(dist(coords, diag = T, upper = T))
  tmp <- lapply(1:length(CorCovList), function(z) 
    cbind(draw = z, 
          cor = as.vector(CorCovList[[z]]), 
          d = as.vector(dmat))
  )
  tmp2 <- as.data.frame(do.call(rbind, tmp))
  return(tmp2)
}


# APOLOGIES THAT THIS CODE IS REALLY REPETITIVE!!

# 10 hr fuels =============================================

# get cov/cor matrix for each fit
fit10h_4trans <- extractCovCor(fits10h[[1]], 100)
fit10h_3trans1 <- extractCovCor(fits10h[[2]], 25)
fit10h_3trans2 <- extractCovCor(fits10h[[3]], 25)
fit10h_3trans3 <- extractCovCor(fits10h[[4]], 25)
fit10h_3trans4 <- extractCovCor(fits10h[[5]], 25)
fit10h_2trans1 <- extractCovCor(fits10h[[6]], 50)
fit10h_2trans2 <- extractCovCor(fits10h[[7]], 50)

# get distance vs cov/cor for each fit
cov10h_4trans <- dist_vs_covcor(fit10h_4trans[[1]], coordsList[[1]])
cov10h_3trans1 <- dist_vs_covcor(fit10h_3trans1[[1]], coordsList[[2]])
cov10h_3trans2 <- dist_vs_covcor(fit10h_3trans2[[1]], coordsList[[3]])
cov10h_3trans3 <- dist_vs_covcor(fit10h_3trans3[[1]], coordsList[[4]])
cov10h_3trans4 <- dist_vs_covcor(fit10h_3trans4[[1]], coordsList[[5]])
cov10h_2trans1 <- dist_vs_covcor(fit10h_2trans1[[1]], coordsList[[6]])
cov10h_2trans2 <- dist_vs_covcor(fit10h_2trans2[[1]], coordsList[[7]])
cor10h_4trans <- dist_vs_covcor(fit10h_4trans[[2]], coordsList[[1]])
cor10h_3trans1 <- dist_vs_covcor(fit10h_3trans1[[2]], coordsList[[2]])
cor10h_3trans2 <- dist_vs_covcor(fit10h_3trans2[[2]], coordsList[[3]])
cor10h_3trans3 <- dist_vs_covcor(fit10h_3trans3[[2]], coordsList[[4]])
cor10h_3trans4 <- dist_vs_covcor(fit10h_3trans4[[2]], coordsList[[5]])
cor10h_2trans1 <- dist_vs_covcor(fit10h_2trans1[[2]], coordsList[[6]])
cor10h_2trans2 <- dist_vs_covcor(fit10h_2trans2[[2]], coordsList[[7]])

# bind posteriors together for 3 and 2 transects
cov10h_3trans <- bind_rows(cov10h_3trans1 %>% mutate(set = 1), 
                           cov10h_3trans2 %>% mutate(set = 2), 
                           cov10h_3trans3 %>% mutate(set = 3), 
                           cov10h_3trans4 %>% mutate(set = 4)) %>% 
  mutate(draw = as.integer(paste0(draw, set)))
cov10h_2trans <- bind_rows(cov10h_2trans1 %>% mutate(set = 1), 
                           cov10h_2trans2 %>% mutate(set = 2)) %>% 
  mutate(draw = as.integer(paste0(draw, set)))
cor10h_3trans <- bind_rows(cor10h_3trans1 %>% mutate(set = 1), 
                           cor10h_3trans2 %>% mutate(set = 2), 
                           cor10h_3trans3 %>% mutate(set = 3), 
                           cor10h_3trans4 %>% mutate(set = 4)) %>% 
  mutate(draw = as.integer(paste0(draw, set)))
cor10h_2trans <- bind_rows(cor10h_2trans1 %>% mutate(set = 1), 
                           cor10h_2trans2 %>% mutate(set = 2)) %>% 
  mutate(draw = as.integer(paste0(draw, set)))



# plot all together
binded_cov10h <- bind_rows(
  cov10h_4trans %>% mutate(ntransects = 4), 
  cov10h_3trans %>% mutate(ntransects = 3),
  cov10h_2trans %>% mutate(ntransects = 2))


p10cov <- ggplot(binded_cov10h, aes(d, cor, color = as.factor(ntransects))) +
  geom_line(lwd = .3, aes(group = paste0(draw, ntransects)), alpha = .5, show.legend = F) +
  coord_cartesian(xlim = c(0,30)) +
  scale_colour_manual(values = c('#dd9933', '#05878a', 'black')) +
  labs(x = 'Distance (m)', y = 'covariance', title = '10 hr fuels GP posterior',
       subtitle = 'number of transects varied') +
  facet_wrap(~ntransects)

# plot all together
binded_cor10h <- bind_rows(
  cor10h_4trans %>% mutate(ntransects = 4), 
  cor10h_3trans %>% mutate(ntransects = 3),
  cor10h_2trans %>% mutate(ntransects = 2))

p10cor <- ggplot(binded_cor10h, aes(d, cor, color = as.factor(ntransects))) +
  geom_line(lwd = .3, aes(group = paste0(draw, ntransects)), alpha = .5, show.legend = F) +
  coord_cartesian(xlim = c(0,30)) +
  scale_colour_manual(values = c('#dd9933', '#05878a', 'black')) +
  labs(x = 'Distance (m)', y = 'correlation', title = '10 hr fuels GP posterior',
       subtitle = 'number of transects varied') +
  facet_wrap(~ntransects)

cowplot::plot_grid(p10cov, p10cor, nrow = 2, rel_heights = c(1, .85))
ggsave('02-data/06-results/sensitivity/FWD10h_ntransects_covcor.png', width = 7, height = 5)




# APOLOGIES THAT THIS CODE IS REALLY REPETITIVE!!

# 1 hr fuels =============================================

# get cov/cor matrix for each fit
fit1h_4trans <- extractCovCor(fits1h[[1]], 100)
fit1h_3trans1 <- extractCovCor(fits1h[[2]], 25)
fit1h_3trans2 <- extractCovCor(fits1h[[3]], 25)
fit1h_3trans3 <- extractCovCor(fits1h[[4]], 25)
fit1h_3trans4 <- extractCovCor(fits1h[[5]], 25)
fit1h_2trans1 <- extractCovCor(fits1h[[6]], 50)
fit1h_2trans2 <- extractCovCor(fits1h[[7]], 50)

# get distance vs cov/cor for each fit
cov1h_4trans <- dist_vs_covcor(fit1h_4trans[[1]], coordsList[[1]])
cov1h_3trans1 <- dist_vs_covcor(fit1h_3trans1[[1]], coordsList[[2]])
cov1h_3trans2 <- dist_vs_covcor(fit1h_3trans2[[1]], coordsList[[3]])
cov1h_3trans3 <- dist_vs_covcor(fit1h_3trans3[[1]], coordsList[[4]])
cov1h_3trans4 <- dist_vs_covcor(fit1h_3trans4[[1]], coordsList[[5]])
cov1h_2trans1 <- dist_vs_covcor(fit1h_2trans1[[1]], coordsList[[6]])
cov1h_2trans2 <- dist_vs_covcor(fit1h_2trans2[[1]], coordsList[[7]])
cor1h_4trans <- dist_vs_covcor(fit1h_4trans[[2]], coordsList[[1]])
cor1h_3trans1 <- dist_vs_covcor(fit1h_3trans1[[2]], coordsList[[2]])
cor1h_3trans2 <- dist_vs_covcor(fit1h_3trans2[[2]], coordsList[[3]])
cor1h_3trans3 <- dist_vs_covcor(fit1h_3trans3[[2]], coordsList[[4]])
cor1h_3trans4 <- dist_vs_covcor(fit1h_3trans4[[2]], coordsList[[5]])
cor1h_2trans1 <- dist_vs_covcor(fit1h_2trans1[[2]], coordsList[[6]])
cor1h_2trans2 <- dist_vs_covcor(fit1h_2trans2[[2]], coordsList[[7]])

# bind posteriors together for 3 and 2 transects
cov1h_3trans <- bind_rows(cov1h_3trans1 %>% mutate(set = 1), 
                           cov1h_3trans2 %>% mutate(set = 2), 
                           cov1h_3trans3 %>% mutate(set = 3), 
                           cov1h_3trans4 %>% mutate(set = 4)) %>% 
  mutate(draw = as.integer(paste0(draw, set)))
cov1h_2trans <- bind_rows(cov1h_2trans1 %>% mutate(set = 1), 
                           cov1h_2trans2 %>% mutate(set = 2)) %>% 
  mutate(draw = as.integer(paste0(draw, set)))
cor1h_3trans <- bind_rows(cor1h_3trans1 %>% mutate(set = 1), 
                           cor1h_3trans2 %>% mutate(set = 2), 
                           cor1h_3trans3 %>% mutate(set = 3), 
                           cor1h_3trans4 %>% mutate(set = 4)) %>% 
  mutate(draw = as.integer(paste0(draw, set)))
cor1h_2trans <- bind_rows(cor1h_2trans1 %>% mutate(set = 1), 
                           cor1h_2trans2 %>% mutate(set = 2)) %>% 
  mutate(draw = as.integer(paste0(draw, set)))



# plot all together
binded_cov1h <- bind_rows(
  cov1h_4trans %>% mutate(ntransects = 4), 
  cov1h_3trans %>% mutate(ntransects = 3),
  cov1h_2trans %>% mutate(ntransects = 2))


p1cov <- ggplot(binded_cov1h, aes(d, cor, color = as.factor(ntransects))) +
  geom_line(lwd = .3, aes(group = paste0(draw, ntransects)), alpha = .5, show.legend = F) +
  coord_cartesian(xlim = c(0,30)) +
  scale_colour_manual(values = c('#dd9933', '#05878a', 'black')) +
  labs(x = 'Distance (m)', y = 'covariance', title = '1 hr fuels GP posterior',
       subtitle = 'number of transects varied') +
  facet_wrap(~ntransects)

# plot all together
binded_cor1h <- bind_rows(
  cor1h_4trans %>% mutate(ntransects = 4), 
  cor1h_3trans %>% mutate(ntransects = 3),
  cor1h_2trans %>% mutate(ntransects = 2))

p1cor <- ggplot(binded_cor1h, aes(d, cor, color = as.factor(ntransects))) +
  geom_line(lwd = .3, aes(group = paste0(draw, ntransects)), alpha = .5, show.legend = F) +
  coord_cartesian(xlim = c(0,30)) +
  scale_colour_manual(values = c('#dd9933', '#05878a', 'black')) +
  labs(x = 'Distance (m)', y = 'correlation', title = '1 hr fuels GP posterior',
       subtitle = 'number of transects varied') +
  facet_wrap(~ntransects)

cowplot::plot_grid(p1cov, p1cor, nrow = 2, rel_heights = c(1, .85))
ggsave('02-data/06-results/sensitivity/FWD1h_ntransects_covcor.png', width = 7, height = 5)
