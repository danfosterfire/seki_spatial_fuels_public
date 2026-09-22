# Goal: Figure out the optimal number of subsamples in each plot to estimate 
# parameters in the litter/duff and FWD models.

# Approach:
# - Simulate litter/duff and FWD data with known parameters. Assume a fixed number
# of plots 
# - alter the number of subsamples within a plot
# - Fit models for each of those iterations, see how it affects the mean and sd
# of the parameters

# Improvements:
# - use more realistic parameters
# - try out 3 transects instead of 4
# - instead of increasing # of plots to 50, keep at 10 and resimulate data several times?

rm(list = ls())
library(tidyverse)
library(rethinking)
library(cowplot)
theme_set(theme_classic())


# Setup -------------------------------------------------------------------


# MODEL
# Y_xp ~ NB(mu_xp, phi)
# log(mu_xp) = a0 + ap + GP_xp
# ap ~ N(0, sigma_p)
# GP ~ MVN(0, K(d))
# [K(d)]_ij = eta^2 * exp(-.5/rho^2 * dist_ij^2)


# look at kernel function
distv <- seq(0, 30, by = 1) #distance vector
fkernel <- function(Alpha, rho, d){
  return(Alpha^2 * exp(-.5/rho^2 * d^2))
}
covv <- fkernel(.5, 2, distv) # covariance vector
plot(distv, covv, type = 'l')


# get coordinates at every meter along the transects. You'll simulate a master
# dataset and then subsample from it.
x <- seq(1.5, 29.5, by = 1)
x <- c(rev(x*-1), x)
y <- x
coords <- rbind(
  cbind(x, y = 0),
  cbind(x = 0, y = y)
)


# SIMULATE MASTER DATASET -----------------------------------------------------------

# define parameters
a0 = 1.5
sigma_p = .5
Alpha = .5
rho = 2
Kappa = .5
nplots = 50
nPerPlot = nrow(coords)
truePars <- c(a0, Kappa, sigma_p, Alpha, rho)
names(truePars) <- c('a0', 'kappa', 'sigma_ap', 'alpha', 'rho' )


# create structure of dataframe
plotID = rep(1:nplots, each = nPerPlot)
locID = rep(1:nPerPlot, times = nplots)
obs = 1:length(plotID)
X = rep(coords[,'x'], times = nplots)
Y = rep(coords[,'y'], times = nplots)
dat <-  tibble(obs, plotID, locID, X, Y)
N <- nrow(dat)

# simulate data
ap <- rnorm(nplots, 0, sigma_p) # plot-level intercepts
distMat <- fields::rdist(coords) # GP for observations in each plot
K <- fkernel(Alpha, rho, distMat)
diag(K) <- diag(K) + 1e-9 #make it positive definite
zGP <- matrix(rnorm(nPerPlot * nplots), nrow = nPerPlot, ncol = nplots)
GP <- t(chol(K)) %*% zGP

# mean model
mu <- rep(NA, N)
for(i in 1:N){
  mu[i] = exp(with(dat, a0 + ap[plotID[i]] + GP[locID[i], plotID[i]]))
}
dat$count <- rnbinom(N, size = Kappa, mu = mu)



# visualize a little bit
# variation across plots
ggplot(dat, aes(plotID, count, group = plotID)) +
  geom_boxplot() 
# variation within plots
ggplot(dat %>% filter(plotID %in% 1:4), aes(X, Y, color = log(count))) +
  geom_point(size = 2) +
  scale_color_viridis_c() +
  coord_equal() +
  facet_wrap(~plotID)



# SUBSAMPLE FROM MASTER DATASET -------------------------------------------

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
x10 <- c(2.5, 3.5, 4.5, 5.5, 7.5, 9.5, 12.5, 15.5, 18.5, 29.5)
x8 <- c(2.5, 3.5, 4.5, 5.5, 7.5, 9.5, 18.5, 29.5)
x6 <- c(2.5, 4.5, 5.5, 9.5, 18.5, 29.5)
x5 <- c(2.5, 5.5, 9.5, 18.5, 29.5)
x4 <- c(2.5, 5.5, 18.5, 29.5)
x3 <- c(2.5, 9.5, 29.5)
getLocations(x10)
getLocations(x8)
getLocations(x6)
getLocations(x5)
getLocations(x4)
getLocations(x3)

# Define function to subset the datasets with a smaller set of locations
subsetDat <- function(locations){
  locations2 <- c(rev(locations*-1), locations)
  datNew <- dat %>% 
    filter(X %in% locations2 | Y %in% locations2) %>%
    mutate(locID = as.integer(as.factor(locID)))
  return(datNew)
}


# PREP DATA FOR STAN MODEL ------------------------------------------------

prepData <- function(dataset, coordinates){
  dataList <- list(
    N = nrow(dataset),
    L = max(dataset$locID),
    P = max(dataset$plotID),
    locID = dataset$locID,
    plotID = dataset$plotID,
    coords = coordinates,
    y = dataset$count
  )
  return(dataList)
}

# get datalists for each of the different designs
locationVectors <- list(x3, x4, x5, x6, x8, x10)
dataLists <- lapply(locationVectors, 
                    function(i) prepData(subsetDat(i), getLocations(i, F)))
str(dataLists) #verify variables are correctly defined


# figure out good priors for Alpha and rho
nsim = 1000
Alphad <- 1/rgamma(nsim, 5, 2)
dens(Alphad); mean(Alphad); HPDI(Alphad)
rhod <- 1/rgamma(nsim, 3, 8)
dens(rhod); mean(rhod); HPDI(rhod)
covv <- fkernel(Alpha, rho, distv) # covariance vector
kernelsim <- sapply(1:100, function(x) fkernel(Alphad[x], rhod[x], distv))
matplot(kernelsim, type = 'l', col = scales::alpha('black', .3), lty = 1)




# RUN STAN MODELS ---------------------------------------------------------

stanNB <- stan_model(file = '04-sensitivity_analyses/GP_NB.stan')

# # do a single run
# fit_test <- sampling(stanNB, dataLists[[5]], iter = 2000, chains = 3, cores = 3)
# 
# # check out fits
# truePars
# precis(fit_test)
# precis(fit_test, depth = 3, pars = c('GP'))
# pairs(fit_test, pars = c('a0', 'rho', 'alpha', 'sigma_ap'))

# fit each of the datasets several times and summarize the posteriors
nsims = 1 #I ran this 10 times to see if the sampler added noise and it doesnt.
paramDF <- list(NULL)
for(s in 1:nsims){
  fits <- lapply(dataLists, 
                 function(i) sampling(stanNB, i, iter = 2000, chains = 3, cores = 3))
  
  summaryList <- list(NULL)
  for(i in 1:length(fits)){
    # get summaries of the parameters
    posts <- extract.samples(fits[[i]])
    
    summaryList[[i]] <- with(posts, cbind(a0, kappa, sigma_ap, alpha, rho)) %>% 
      data.frame() %>% 
      pivot_longer(everything()) %>% 
      group_by(name) %>% 
      summarise(median = median(value),
                SD = sd(value),
                lower90 = HPDI(value, .9)[1],
                upper90 = HPDI(value, .9)[2],
                lower50 = HPDI(value, .5)[1],
                upper50 = HPDI(value, .5)[2]) %>% 
      mutate(nPerTransect = dataLists[[i]]$L / 4,
             simulation = s)
  }
  
  paramDF[[s]] <- bind_rows(summaryList)

}

paramDF <- bind_rows(paramDF)
saveRDS(paramDF, '02-data/06-results/sensitivity/FWD_subsamples.RDS')



# COMPARE PARAMETERS ------------------------------------------------------

trueParsDF <- data.frame(name = names(truePars), value = truePars)

# viz posteriors
ggplot(paramDF, aes(nPerTransect, median)) +
  geom_linerange(aes(ymin = lower50, ymax = upper50), 
                 #position = position_nudge(x = .01),
                 alpha = .5, size = 1) +
  geom_pointrange(aes(ymin = lower90, ymax = upper90), 
                  #position = position_dodge2(width = .6),
                  alpha = .5) +
  geom_hline(data = trueParsDF, aes(yintercept = value), lty = 2) +
  facet_wrap(~name, scales = 'free_y') +
  labs(x = 'samples per transect', y = 'median, 90% & 50% HDPI',
       title = 'FWD tallies, 50 simulated plots')

ggsave('02-data/06-results/sensitivity/FWD_subsamples.pdf', width = 7, height = 4)
