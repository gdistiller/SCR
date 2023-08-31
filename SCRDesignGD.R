#SCR design, Feb 2023
#this script firstly explores how to use secrdesign
#it then plays with generating different designs
#sims were then run on the cluster, and the results are explored here too

#Useful packages: elevatr (generates elevation from (x,y)), 
#secrdesign (for lots of things, including optimized grid), 
#"terrain" function from raster package

tri <- terrain(elev, opt = "TRI", filename = "TRI.tiff")

#takes an elevation raster and makes a terrain ruggedness raster

#also 'extract' function, extract value of a raster at a point

#get proposed design from GAoptim for Enrm, plus other designs
#then can use std secrdesign tools to run the simulations

#hint: first just generate raw data summaries and check 

############################
#exploring secrdesign
############################

#code to test expected numbers

library(secr) ; library(secrdesign) ; library(kofnGA)

# generate scenario dataframe
scen <- make.scenarios (trapsindex = 1:2, detectfn = 'HHN', D = c(20),
                        lambda0 = c(0.05,0.2), sigma = 1:2, noccasions = 1)

# set sigma = k/sqrt(D) for k = 50, 100
scen$sigma <- c(50,100)[scen$sigma] / sqrt(scen$D)

# detector layouts
traplist <- list(
#  make.grid(6,6, detector = 'multi', spacing = 20),
  make.grid(6,6, detector = 'proximity', spacing = 20),
#  make.grid(6,6, detector = 'count', spacing = 20),
#  make.grid(10,10, detector = 'multi', spacing = 20),
  make.grid(10,10, detector = 'proximity', spacing = 20)
#  make.grid(10,10, detector = 'count', spacing = 20)
)

# deterministic summary: expected counts
nrm <- scenarioSummary(scen, traplist)

# stochastic simulations
sumnrm <- function(CH) {
  c(
    n = nrow(CH),
    r = sum(CH) - nrow(CH),
    # moves(CH, names = TRUE) to dodge bug in secr < 4.5.2
    m = sum(unlist(secr::moves(CH, names = TRUE))>0, na.rm = TRUE),
    n2 = sum( apply(apply(CH, c(1,3), sum)>0, 1, sum) >1 )
  )
}

sim <- run.scenarios(scen, trapset = traplist, nrepl = 1000,
                     extractfn = sumnrm, fit = FALSE)

get.est <- function(df) {
 est <- df[4,2]
}

vec <- lapply(sims2$output[[1]],get.est)
mean(as.numeric(vec))

# collate deterministic and stochastic results
out <- data.frame(nrm, simout)
plotc <- function (v = 'n') {
  dtype <- c('multi','proximity','count')
  detector <- rep(dtype,2)[out$trapsindex]
  out2 <- matrix(nrow=3, ncol = 2, dimnames = list(dtype, c('mean','sd')))
  for (d in 1:3) {
    OK <- (detector == dtype[d])
    RB <- (out[OK,v] - out[OK, paste0('E',v)]) / out[OK, paste0('E', v)]
    # use log scales to spread values
    plot(out[OK, paste0('E',v)], out[OK,v], log='xy',
         xlab = 'Expected', ylab = 'simulated')
    abline(0,1) # y = x line
    # return mean and sd of estimated relative bias
    out2[d,] <- round(c(mean = mean(RB), sd = sd(RB)),5)
    mtext (side=3, line=0.2, paste(v, " ", dtype[d]), cex = 0.9)
  }
  out2
}

par(mfrow=c(4,3), mgp=c(2.3,0.6,0), mar=c(4,4,2,1), pty='s')
cat("Relative discrepancy between expected and simulated counts\n")
cat("Number of individuals\n")
plotc('n')
cat("Number of recaptures\n")
plotc('r')
cat("Number of movements\n")
plotc('m')
cat("Number of individuals at two or more detectors\n")
plotc('n2')

####################################################################
##scenarios with a regular box as mask, count detectors, 2 levels of D, 2 levels of lambda0
####################################################################

# generate scenario dataframe
scen <- make.scenarios (trapsindex = 1, detectfn = 'HHN', D = c(0.001),
                        lambda0 = c(0.5), sigma = c(3000), noccasions = 1)

# create trap array
traps <- make.grid(10,10, detector = 'count', spacing = 4000)

# deterministic summary: expected counts ; note that this function excludes grouped scenarios
nrm <- scenarioSummary(scen, traps)

# stochastic simulations
raw <- run.scenarios(scen, trapset = traps, nrepl = 50,
                     extractfn = identity, fit = FALSE)
summary(raw)
raw$output[[1]]
str(raw$output[[1]])  #output from 1st scenario, CHs for each of the 5 reps

#fit models as a 2nd step
sims1 <- fit.models(raw, fit.args = list(detectfn = 'HHN'), fit = TRUE, ncores = 4, byscenario = FALSE)
summary(sims1)

D.ests <- select.stats(sims1, parameter = "D", statistics = c("estimate","SE.estimate","RB", "RSE", "COV","lcl","ucl"))
L0.ests <- select.stats(sims1, parameter = "lambda0", statistics = c("estimate","SE.estimate","RB", "RSE", "COV","lcl","ucl"))
Sig.ests <- select.stats(sims1, parameter = "sigma", statistics = c("estimate","SE.estimate","RB", "RSE", "COV","lcl","ucl"))

summary(D.ests, dec = 6, fields = c("n", "mean", "se"), alpha = 0.05, type = "list")
summary(L0.ests, dec = 3, fields = c("n", "mean", "se"), alpha = 0.05, type = "list")
summary(Sig.ests, dec = 2, fields = c("n", "mean", "se"), alpha = 0.05, type = "list")

##fit model manually, to 1st CH from 7th scen
CH1 <- (raw$output[[7]][[1]])
mod1 <- secr.fit(CH1, detectfn = 14, buffer = 4000)
summary(mod1)
exampleOutput <- predict(mod1)

#sim data and fit model
sims2 <- run.scenarios(scen, trapset = traps, nrepl = 5, fit = TRUE)

#Note that the output to run.scenarios includes the result of extractfn for each rep
#when fit = T, the default output from each rep is the result of applying predict 
#i.e. a data frame of 'real' estimates with se's

sims3 <- run.scenarios(scen, trapset = traps, nrepl = 5, fit = TRUE, extractfn = predict)
find.param(sims3)
Dests <- select.stats(sims3, parameter = "sigma", statistics = c("estimate","SE.estimate","RB", "RSE", "COV","lcl","ucl"))
#note that validate() can be used to remove rouge values
apply(as.matrix(Dests$output[[1]]),2,mean)
#summary can be used on the 'selectedstatistics' object
summary(Dests, dec = 2, fields = c("n", "mean", "se"), alpha = 0.05, type = "list")

#plot estimates from 'selectedstatistics' object
par(mfrow = c(2,2))
plot(Dests, type = "hist", statistic = "estimate")
plot(Dests, type = "CI")

############################
#generating optimal designs
############################
#use GAoptim to get proposed designs
#regular 2 sigma and optimised grid follow

# user parameters
lambda0 <- 2 
dens_per_100km2 <- 2 # mean animal density per 100km2, SLs are ~1
D <- dens_per_100km2 / 10000
sigma <- 3000
b_ac <- 1

nT <- 20
buffer <- 4 * sigma
dt <- "count"

# an artificial example, note the bottom-left corner is at origin
msk <- make.mask(type = 'rectangular', spacing = 1000, nx = 30, ny = 20, buffer = 0)
alltrps <- make.grid(nx = 29, ny = 19, origin = c(1000,1000), spacing = 1000, detector = "count")

# 50 generations for demonstration, use more in practice
opt <- GAoptim(mask = msk, alltraps = alltrps, ntraps = 20, detectpar = list(lambda0 = lambda0, sigma = sigma), 
               detectfn = 'HHN', D = D, noccasions = 1, ngen = 5, verbose = 1)

plot(msk)
plot(opt$optimaltraps, add = TRUE)
minnrRSE(opt, distribution = 'binomial')

# Using a criterion function
# En2 is unsuitable as a criterion function as it returns 2 values
# This function selects the second as the (unique) criterion
fn <- function(...) En2(...)[2]
opt2 <- GAoptim(msk, alltrps, ntraps = 20, detectpar = list(lambda0 = lambda0, sigma = sigma), 
                detectfn = 'HHN', D = D, noccasions = 1, ngen = 3, verbose = 1, criterion = fn)

plot(msk)
plot(opt2$optimaltraps, add = TRUE)
minnrRSE(opt2, distribution = 'binomial')

# grid design under constant D
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(stringr)
library(purrr)
library(secr)
library(secrdesign)
library(oSCR)
library(raster)
library(gdistance)
library(kofnGA)

# create mask and then a buffered mask (the mask used for the designs)
mask <- secr::raster(msk)
mask_df <- data.frame(coordinates(mask))
mask <- read.mask(data = mask_df)
plot(mask, dots=F)

buffer_mult <- ceiling(buffer / cellsize)
newmask_df <- expand.grid(x = seq(from = min(mask_df$x) - buffer_mult * cellsize, 
                                  to = max(mask_df$x) + buffer_mult * cellsize, 
                                  by = cellsize),
                          y = seq(from = min(mask_df$y) - buffer_mult * cellsize, 
                                  to = max(mask_df$y) + buffer_mult * cellsize, 
                                  by = cellsize)) 
newmask_df <- newmask_df %>% 
  mutate(keep = (apply(e2dist(as.matrix(newmask_df), as.matrix(mask_df)), 1, min) <= buffer)) %>%
  filter(keep == TRUE) %>% dplyr::select(-keep)
newmask <- read.mask(data = newmask_df)
plot(newmask,axes=T)

# create grid of possible trap locations
# reduce resolution of mesh so have fewer possible camera locations
cellsize <- 2/3 * sigma
red_factor <- cellsize[1] / attr(msk, "spacing")
if ((trunc(red_factor) - red_factor) != 0) stop("Check spacing, causing non-integer reduction factor for mesh")

alltraps_df <- data.frame(coordinates(data.frame(alltrps)))
alltraps <- as.matrix(alltraps_df)
plot(read.mask(data = alltraps_df), col = 'blue', add = TRUE, axes =T)  

#####################################
grid_traps_all <- data.frame(nT = as.integer(), sigma = as.numeric(), beta0 = as.numeric(), buffer = as.numeric(), dt = as.character(),
                             x = as.numeric(), y = as.numeric(), trap_id = as.integer())
opt_grid_traps_all <- data.frame(nT = as.integer(), sigma = as.numeric(), beta0 = as.numeric(), buffer = as.numeric(), dt = as.character(),
                                 x = as.numeric(), y = as.numeric(), trap_id = as.integer())
######################
#regular 2 sigma grid
######################

for(n in c(5,10,15)){

# hacky way to make a polygon just bigger than boundary of mask points
# used later to decide which randomly generated detectors are on the mask and so allowed
sap_1 <- mask_df %>% st_as_sf(coords = c("x", "y")) %>% st_buffer(dist = cellsize, endCapStyle = "SQUARE") %>% st_union() 
sap_1 <- sap_1 %>% st_buffer(dist = -cellsize * 0.4, endCapStyle = "SQUARE")

# place a grid over the area, with cells 2 * sigma apart
my_grid <- st_make_grid(sap_1, cellsize = c(2 * sigma, 2 * sigma), what = "centers") %>%
  st_intersection(sap_1)
grid_traps_full <- st_coordinates(my_grid) %>% as.data.frame() %>% rename(x = X, y = Y)

# choose a subset of nT of these, starting from a random trap and then choosing nT nearest neighbours
all_grid_traps <- list()

grid_traps <- grid_traps_full
set.seed(700)
xy_rand <- grid_traps[sample(1:nrow(grid_traps), 1), ]
grid_traps$dist2pt <- (grid_traps$x - xy_rand$x)^2 + (grid_traps$y - xy_rand$y)^2
grid_traps <- grid_traps %>% filter(rank(dist2pt, ties.method = "random") <= n) %>% mutate(trap_id = 1)
grid_traps_all <- rbind(grid_traps_all, cbind(nT = n, lambda0 = lambda0, sigma = sigma, buffer = buffer, dt = dt, grid_traps))

######################
#optimised grid
######################
grid_traps <- read.traps(data = grid_traps, detector = dt)

# compute optimal spacing for previous grid 
dd <- optimalSpacing(D = D, traps = grid_traps, detectpar = list(lambda0 = lambda0, sigma = sigma), noccasions = 1)$rotRSE$optimum.spacing

# sometimes optimal spacing is too big to fit desired number of traps in; if so, reduce dd until it is
n_opt_grid <- 0
red_dd <- 0.9

while(n_opt_grid < n){
  
  dd <- dd * red_dd
  
  # place a grid over the area, with cells optimally spaced
  my_grid <- st_make_grid(sap_1, cellsize = c(dd, dd), what = "centers") %>% st_intersection(sap_1)
  opt_grid_traps_full <- st_coordinates(my_grid) %>% as.data.frame() %>% rename(x = X, y = Y)
  
  # choose a subset of nT of these points, starting from a random trap and then choosing nT nearest neighbours
  opt_grid_traps <- opt_grid_traps_full
  set.seed(700)
  xy_rand <- opt_grid_traps[sample(1:nrow(opt_grid_traps), 1), ]
  opt_grid_traps$dist2pt <- (opt_grid_traps$x - xy_rand$x)^2 + (opt_grid_traps$y - xy_rand$y)^2
  opt_grid_traps <- opt_grid_traps %>% filter(rank(dist2pt, ties.method = "random") <= n) %>% mutate(trap_id = 1)
  
  n_opt_grid <- nrow(opt_grid_traps)
  
}

opt_grid_traps_all <- rbind(opt_grid_traps_all, cbind(nT = n, lambda0 = lambda0, sigma = sigma, buffer = buffer, dt = dt, opt_grid_traps))

}


