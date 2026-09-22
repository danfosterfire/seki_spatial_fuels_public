# predicting fuels with a GP intercept
rm(list = ls())

library(tidyverse)
library(rstan)
library(rethinking)
#library(cmdstanr)
theme_set(theme_bw())

# MODEL
# Y_xp ~ Poisson(mu_xp)
# log(mu_xp) = a0 + ap + GP_xp
# ap ~ N(0, sigma_p)
# GP ~ MVN(0, K(d))
# [K(d)]_ij = eta^2 * exp(-.5/rho^2 * dist_ij^2)

# Create coordinates of fuel samples
x <- c(2.5, 3.5, 4.5, 5.5, 7.5, 9.5, 19.5, 29.5)
x <- c(rev(x*-1), x)
y <- x
coords <- rbind(
  cbind(x, y = 0),
  cbind(x = 0, y = y)
)
plot(coords)

# look at pairwise distances
distmatLower <- fields::rdist(coords, compact = T)
hist(as.vector(distmatLower), breaks = 30)

# look at kernel function
distv <- seq(0, 30, by = 1) #distance vector
fkernel <- function(eta, rho, d){
  return(eta^2 * exp(-.5/rho^2 * d^2))
}
covv <- fkernel(.5, 5, distv) # covariance vector
plot(distv, covv, type = 'l')





# SIMULATE DATA -----------------------------------------------------------

# define parameters
a0 = 1.5
sigma_p = .5
eta = .5
rho = 5
nplots = 10
nPerPlot = nrow(coords)

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
K <- fkernel(eta, rho, distMat)
zGP <- matrix(rnorm(nPerPlot * nplots), nrow = nPerPlot, ncol = nplots)
GP <- t(chol(K)) %*% zGP
# mean model
mu <- rep(NA, N)
for(i in 1:N){
  mu[i] = exp(with(dat, a0 + ap[plotID[i]] + GP[locID[i], plotID[i]]))
}
dat$count <- rpois(N, mu)

# remove some data to test how the model handles missing data
rmv <- sample(nrow(dat), round(.05 * nrow(dat)))
dat <- dat[-rmv,]

# visualize a little bit
# variation across plots
ggplot(dat, aes(plotID, count, group = plotID)) +
  geom_boxplot() 
# variation within plots
ggplot(dat %>% filter(plotID %in% 1:4), aes(X, Y, color = count)) +
  geom_point(size = 2) +
  scale_color_viridis_c() +
  coord_equal() +
  facet_wrap(~plotID)




# PREP DATA FOR STAN MODEL ------------------------------------------------
head(dat)

dataList <- list(
  N = nrow(dat),
  L = max(dat$locID),
  P = max(dat$plotID),
  locID = dat$locID,
  plotID = dat$plotID,
  coords = coords,
  y = dat$count
  )
str(dataList) #verify variables are correctly defined


# figure out good priors for eta and rho
nsim = 1000
etad <- 1/rgamma(nsim, 2, 1)
rhod <- 1/rgamma(nsim, 2, 1)
dens(etad); mean(etad)
dens(rhod); mean(rhod)
covv <- fkernel(eta, rho, distv) # covariance vector
kernelsim <- sapply(1:nsim, function(x) fkernel(etad[x], rhod[x], distv))
matplot(kernelsim, type = 'l', col = scales::alpha('black', .3), lty = 1)





# RUN STAN MODEL ----------------------------------------------------------

stan1 <- stan_model(file = 'stanmods/GPcountsv1.stan')

fit1 <- sampling(stan1, dataList, iter = 2000, chains = 3, cores = 3)
precis(fit1)
precis(fit1, depth = 3, pars = c('GP'))
pairs(fit1, pars = c('a0', 'rho', 'eta', 'sigma_ap'))

