#Aug 2023
#settling on base setup for realistic leopard data
#using 1.5 sigma for trap spacing
library(secr) ; library(secrdesign)

#SCR parameters
MsigmaFactor = 2 #male sigma double female sigma
MDFactor = 0.25  #male D 1/4 female D

sigmaF = 3000 ; sigmaM = MsigmaFactor * sigmaF
DF = 0.00015 ; DM = MDFactor * DF
lambda0F = 5 ; lambda0M = 4.5

trapsF <- make.grid(7,7, detector = 'count', spacing = 1.5*sigmaF)
trapsM <- make.grid(7,7, detector = 'count', spacing = 1.5*sigmaM)

scen.M <- make.scenarios (trapsindex = 1, detectfn = 'HHN', D = DM,
                          lambda0 = lambda0M, sigma = sigmaM, noccasions = 1)
scen.F <- make.scenarios (trapsindex = 1, detectfn = 'HHN', D = DF,
                          lambda0 = lambda0F, sigma = sigmaF, noccasions = 1)
scenarioSummary(scen.F, trapsF)
scenarioSummary(scen.M, trapsM)

####################
#functions
####################

# Objective function with weighted average included
# Fn determines strata from the length of the sigma input
# still uses same criteria which determines what gets weighted
# max argument used to get max of En / Er / Em rather than weighted avg 
Enr.2pops <- function (v, alltraps, mask, detectpar, detectfn, noccasions, 
                 detector, D, crit, weights = c(0.5,0.5), max = FALSE) {
  
  traps <- subset(alltraps, v)
  
  strata = FALSE
  if (length(detectpar[[2]])>1) strata = TRUE
  
  if (is.function(crit)) {
    -crit(D = D, traps = traps, mask = mask, noccasions = noccasions, 
          detectpar = detectpar, detectfn = detectfn)[1]
  }
  if (strata==FALSE){
    if (crit<5) {
      enrm <- Enrm(D = D, traps = traps, mask = mask, noccasions = noccasions, 
                   detectpar = detectpar, detectfn = detectfn)
      c(-enrm[1], -enrm[2], -enrm[3], -(min(enrm[1],enrm[2])))[crit]    
    }
    else {
      en2 <- En2(D = D, traps = traps, mask = mask, noccasions = noccasions, 
                 detectpar = detectpar, detectfn = detectfn)
      -c(en2[2], sum(en2))[crit-5]
    } 
  } else {  #for 2 populations
    if (crit<5) {
      enrm1 <- Enrm(D = D[1], traps = traps, mask = mask, noccasions = noccasions, 
                    detectpar = list(lambda0 = detectpar[[1]][1], sigma = detectpar[[2]][1]), detectfn = detectfn)
      enrm2 <- Enrm(D = D[2], traps = traps, mask = mask, noccasions = noccasions, 
                    detectpar = list(lambda0 = detectpar[[1]][2], sigma = detectpar[[2]][2]), detectfn = detectfn)
      
      if (max==F){
        obj.n <- -(weights[1]*enrm1[1] + weights[2]*enrm2[1])
        obj.r <- -(weights[1]*enrm1[2] + weights[2]*enrm2[2])
        obj.m <- -(weights[1]*enrm1[3] + weights[2]*enrm2[3])
        obj.min <- -(weights[1]*(min(enrm1[1],enrm1[2]) + weights[2]*min(enrm2[1],enrm2[2])))
        
        c(obj.n, obj.r, obj.m, obj.min)[crit]  
      } else {
        obj.n <- -(max(enrm1[1], enrm2[1]))
        obj.r <- -(max(enrm1[2], enrm2[2]))
        obj.m <- -(max(enrm1[3], enrm2[3]))
        obj.min <- -(max(min(enrm1[1],enrm1[2]), min(enrm2[1],enrm2[2])))
        
        c(obj.n, obj.r, obj.m, obj.min)[crit]  
      }
    }
    else {
      en21 <- En2(D = D[1], traps = traps, mask = mask, noccasions = noccasions, 
                  detectpar = list(lambda0 = detectpar[[1]][1], sigma = detectpar[[2]][1]), detectfn = detectfn)
      en22 <- En2(D = D[2], traps = traps, mask = mask, noccasions = noccasions, 
                  detectpar = list(lambda0 = detectpar[[1]][2], sigma = detectpar[[2]][2]), detectfn = detectfn)
      
      obj.n2 <- (weights[1]*en21[2]+weights[2]*en22[2])
      obj.pb <- (weights[1]*sum(en21)+weights[2]*sum(en22))
      
      -c(obj.n2, obj.pb)[crit-5]
    }
  }
}

#version of GAoptim that uses the OF above
GAoptimTest <- function(mask, alltraps, ntraps, detectpar, noccasions, detectfn = c("HHN", "HHR", "HEX", "HAN", "HCG"),
                        D = NULL, criterion = 4, weights = c(0.5, 0.5), seed = NULL, ...){
  
  detectfn <- match.arg(detectfn)
  
  ## criterion (1 = En, 2 = Er, 3 = Em, 4 = min(En,Er), 5 = En2, 6 = En+En2)
  if (!is.function(criterion) && (criterion<1 || criterion>6)) stop ("invalid criterion code")
  
  if(missing(mask)) stop("Must supply a 'mask' object (coords of the study area)")
  if(missing(alltraps))   stop("Must supply a 'traps' object (all possible trap locations)")
  
  if (!inherits(mask, "mask")) stop ("mask should be a mask object")
  if (!inherits(alltraps, "traps")) stop ("alltraps should be a traps object")    
  
  detector <- match.arg(detector(alltraps), choices = c("count", "proximity", "multi"))
  if (noccasions == 1 && detector == "multi") stop ("multi detector requires > 1 occasion")
  if(!is.null(seed)) set.seed(seed)
  
  if (ms(mask) || ms(traps)) stop ("mask and traps should be single-session")
  
  #---------------------------------------------------------------------------
  
  des <- kofnGA::kofnGA(n = nrow(alltraps), 
                        k  = ntraps, 
                        OF = Enr.2pops,
                        ...,
                        alltraps    = alltraps,
                        mask        = mask,
                        detectpar   = detectpar,
                        noccasions  = noccasions,
                        detectfn    = detectfn,
                        detector    = detector,
                        D           = if (is.null(D)) 1 else D,
                        crit        = criterion)
  
  optimaltraps <- subset(alltraps, des$bestsol)
  
  if (!is.null(D)) {
    if (length(D)>1){
      optimalenrm1 <- Enrm(D = D[1], traps = optimaltraps, mask = mask, 
                          noccasions = noccasions, detectpar = list("lambda0" = detectpar[[1]][1], "sigma" = detectpar[[2]][1]), detectfn = detectfn)
      optimalenrm2 <- Enrm(D = D[2], traps = optimaltraps, mask = mask, 
                          noccasions = noccasions, detectpar = list("lambda0" = detectpar[[1]][2], "sigma" = detectpar[[2]][2]), detectfn = detectfn)
      optimalenrm <- c(optimalenrm1,optimalenrm2)
    } else {
      optimalenrm <- Enrm(D = D, traps = optimaltraps, mask = mask, 
                          noccasions = noccasions, detectpar = detectpar, detectfn = detectfn)
    }
  }
  else {
    optimalenrm <- NULL
  }
  
  out <- list(
    mask         = mask, 
    alltraps     = alltraps, 
    detectpar    = detectpar, 
    noccasions   = noccasions,
    detectfn     = detectfn,
    D            = D,
    criterion    = criterion,
    des          = des, 
    optimaltraps = optimaltraps,
    optimalenrm  = optimalenrm
    ## do not include minnrRSE - it depends on extra arguments CF, distribution
  )
  
  class(out) <- "GAoptim"
  out
  
}

# from oSCR package
e2dist <- function (x, y) {
  if (!is.matrix(x)) 
    x <- as.matrix(x)
  if (!is.matrix(y)) 
    y <- as.matrix(y)
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

#fn to generate a grid design
#randomly picks a grid point and selects nearest neighbours
grid.design <- function(mask, lambda0, sigma, nT){
  mask <- secr::raster(mask)
  mask_df <- data.frame(coordinates(mask))
  cellsize <- 2/3 * sigma
  buffer <- 4 * sigma
  
  #df to store results
  grid_traps_all <- data.frame(nT = as.integer(), sigma = as.numeric(), beta0 = as.numeric(), buffer = as.numeric(), dt = as.character(),
                               x = as.numeric(), y = as.numeric(), trap_id = as.integer())
  
  # hacky way to make a polygon just bigger than boundary of mask points
  # used later to decide which randomly generated detectors are on the mask and so allowed
  sap_1 <- mask_df %>% st_as_sf(coords = c("x", "y")) %>% st_buffer(dist = 0.6*cellsize, endCapStyle = "SQUARE") %>% st_union() 
  sap_1 <- sap_1 %>% st_buffer(dist = -cellsize * 0.4, endCapStyle = "SQUARE")
  
  # place a grid over the area, with cells 2 * sigma apart
  my_grid <- st_make_grid(sap_1, cellsize = c(2 * sigma, 2 * sigma), what = "centers") %>%
    st_intersection(sap_1)
  grid_traps_full <- st_coordinates(my_grid) %>% as.data.frame() %>% rename(x = X, y = Y)
  
  # choose a subset of nT of these, starting from a random trap and then choosing nT nearest neighbours
  grid_traps <- grid_traps_full
  xy_rand <- grid_traps[sample(1:nrow(grid_traps), 1), ]
  grid_traps$dist2pt <- (grid_traps$x - xy_rand$x)^2 + (grid_traps$y - xy_rand$y)^2
  grid_traps <- grid_traps %>% filter(rank(dist2pt, ties.method = "random") <= nT) %>% mutate(trap_id = 1)
  grid_traps_all <- rbind(grid_traps_all, cbind(nT = nT, lambda0 = lambda0, sigma = sigma, buffer = buffer, dt = dt, grid_traps))
  return(list("Selected grid" = grid_traps, "Selected grid df" = grid_traps_all))
}

#fn to create the buffered mask from df of coordinates
Create.buffmask <- function(sigma, mask){
  
  #get mask coords
  mask <- secr::raster(mask)
  mask_df <- data.frame(coordinates(mask))
  
  # create a buffered mask (the mask used for the designs) from the mask above
  buffer <- 4 * sigma
  cellsize <- 2/3 * sigma
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
  
  return(list("Buffered mask" = newmask ))
}

################################################################
#libraries
library(secr)
library(secrdesign)
library(kofnGA)
library(raster)
library(dplyr)
library(tidyr)
library(sf)
#library(oSCR) #for e2dist()

library(ggplot2)

library(stringr)
library(purrr)
library(gdistance)

#basic sim, to be expanded
#uses SCR pars as above, with a regular grid of 49 traps
#OFs compared are: WtEnr, Enr, and regular 1.5 sigma grid  

#initially 1 rep
##################################
#first propose different designs
###################################
#SCR parameters
MsigmaFactor = 2 #male sigma double female sigma
MDFactor = 0.25  #male D 1/4 female D

sigmaF = 3000 ; sigmaM = MsigmaFactor * sigmaF ; sigma <- c(sigmaF, sigmaM)
DF = 0.00015 ; DM = MDFactor * DF ; D <- c(DF, DM)
lambda0F = 5 ; lambda0M = 4.5 ; lambda0 <- c(lambda0F,lambda0M) 

nT <- 20
bufferF <- 4 * sigmaF ; bufferM <- 4 * sigmaM
dt <- "count"

# an artificial example, note the bottom-left corner is at origin
msk <- make.mask(type = 'rectangular', spacing = 1000, nx = 48, ny = 48, buffer = 0)
alltrps <- make.grid(nx = 29, ny = 19, origin = c(1000,1000), spacing = 1000, detector = "count")

mask.buffF <- Create.buffmask(sigma = sigmaF, mask = msk)
mask.buffM <- Create.buffmask(sigma = sigmaM, mask = msk)

############################################################
#get proposed designs, for Enr (F/M/both) and grid (F/M/avg)
############################################################
Proposed.Traps.Enr <- vector(mode = "list", length = 3) #F/M/both
Proposed.Traps.Grid <- vector(mode = "list", length = 3)#F/M/avg

  #first use GA to get Enr designs
  # 5 generations to check, increase for actual sim
  optF <- GAoptimTest(mask = msk, alltraps = alltrps, ntraps = nT, detectpar = list(lambda0 = lambda0F, sigma = sigmaF), 
                     criterion = 4, detectfn = 'HHN', D = DF, noccasions = 1, ngen = 5, verbose = 1)
  optM <- GAoptimTest(mask = msk, alltraps = alltrps, ntraps = nT, detectpar = list(lambda0 = lambda0M, sigma = sigmaM), 
                      criterion = 4, detectfn = 'HHN', D = DM, noccasions = 1, ngen = 5, verbose = 1)
  optBoth <- GAoptimTest(mask = msk, alltraps = alltrps, ntraps = nT, detectpar = list(lambda0 = c(lambda0F,lambda0M), sigma = c(sigmaF,sigmaM)), 
                      criterion = 4, detectfn = 'HHN', D = c(DF,DM), noccasions = 1, ngen = 5, verbose = 1)
  
  Proposed.Traps.Enr[[1]] <- optF$optimaltraps
  Proposed.Traps.Enr[[2]] <- optM$optimaltraps
  Proposed.Traps.Enr[[3]] <- optBoth$optimaltraps
  
  #now a regular grid design for each sigma
  gridF <- grid.design(msk,lambda0F, sigmaF,16)
  gridM <- grid.design(msk,lambda0M, sigmaM,16)
  gridBoth <- grid.design(msk,mean(c(lambda0F,lambda0M)), mean(c(sigmaF,sigmaM)),16)

  trapsF <- read.traps(data = gridF[[1]], detector = dt)
  trapsM <- read.traps(data = gridM[[1]], detector = dt)
  trapsBoth <- read.traps(data = gridBoth[[1]], detector = dt)
  
  Proposed.Traps.Grid[[1]] <- trapsF
  Proposed.Traps.Grid[[2]] <- trapsM
  Proposed.Traps.Grid[[3]] <- trapsBoth
  
  #use secrdesign to run simulations comparing the designs
  scen <- make.scenarios (trapsindex = 1:6, detectfn = 'HHN', D = DF, lambda0 = lambda0F, 
                          sigma = sigmaF, noccasions = 1, groups = c('F','M'))
  
  male <- scen$group == 'M'
  scen$D[male] <- DM
  scen$lambda0[male] <- lambda0M
  scen$sigma[male] <- sigmaM
  
  traps.list <- list("Enr (F)" = Proposed.Traps.Enr[[1]], "Enr (M)" = Proposed.Traps.Enr[[2]], "Enr (both)" = Proposed.Traps.Enr[[3]],"Grid (F)" = Proposed.Traps.Grid[[1]],"Grid (M)" = Proposed.Traps.Grid[[2]], "Grid (avg)" = Proposed.Traps.Grid[[3]])
  mask.list <- list("Mask (F)" = mask.buffF[[1]], "Mask (M)" = mask.buffM[[1]], "Mask (M)" = mask.buffM[[1]], "Mask (F)" = mask.buffF[[1]], "Mask (M)" = mask.buffM[[1]], "Mask (M)" = mask.buffM[[1]])
  
  sims.group.both <- run.scenarios(2, scen, trapset = traps.list, fit = TRUE, extractfn = predict,
                                   fit.args = list(detectfn = 'HHN', model = list(D~g, sigma~g, lambda0~g), groups = "group"),
                                   maskset = mask.list, ncores = 4, byscenario = FALSE)
  
  
  #as function, creates mask from trap array
  #raw = T keeps raw datasets
  #if mask is NULL then alltraps used for both
  Sims1 <- function(mask = NULL, alltraps, nT, sigmaF, sigmaM, lambda0F, lambda0M, DF, DM, raw = FALSE, nreps = 10, ngen = 50){
    
    if (is.null(mask)) mask <- make.mask(alltraps, buffer = 0)
      
    mask.buffF <- Create.buffmask(sigma = sigmaF, mask = mask)  #creates a buffered mask with spacing 2/3 sigma
    mask.buffM <- Create.buffmask(sigma = sigmaM, mask = mask)
    
    Proposed.Traps.Enr <- vector(mode = "list", length = 3) #F/M/both
    Proposed.Traps.Grid <- vector(mode = "list", length = 3)#F/M/avg
    
    #first use GA to get Enr designs, both mask and alltrps are 2/3 sigma spacing
    optF <- GAoptimTest(mask = mask.buffF[[1]], alltraps = alltraps, ntraps = nT, detectpar = list(lambda0 = lambda0F, sigma = sigmaF), 
                        criterion = 4, detectfn = 'HHN', D = DF, noccasions = 1, ngen = ngen, verbose = 0)
    optM <- GAoptimTest(mask = mask.buffM[[1]], alltraps = alltraps, ntraps = nT, detectpar = list(lambda0 = lambda0M, sigma = sigmaM), 
                        criterion = 4, detectfn = 'HHN', D = DM, noccasions = 1, ngen = ngen, verbose = 0)
    optBoth <- GAoptimTest(mask = mask.buffM[[1]], alltraps = alltraps, ntraps = nT, detectpar = list(lambda0 = c(lambda0F,lambda0M), sigma = c(sigmaF,sigmaM)), 
                           criterion = 4, detectfn = 'HHN', D = c(DF,DM), noccasions = 1, ngen = ngen, verbose = 0)
    
    Proposed.Traps.Enr[[1]] <- optF$optimaltraps
    Proposed.Traps.Enr[[2]] <- optM$optimaltraps
    Proposed.Traps.Enr[[3]] <- optBoth$optimaltraps
    
    #now a regular grid design for each sigma
    #grid of 2 sigma is used to generate locations
    gridF <- grid.design(mask,lambda0F, sigmaF, nT)
    gridM <- grid.design(mask,lambda0M, sigmaM, nT)
    gridBoth <- grid.design(mask,mean(c(lambda0F,lambda0M)), mean(c(sigmaF,sigmaM)), nT)
    
    trapsF <- read.traps(data = gridF[[1]], detector = "count")
    trapsM <- read.traps(data = gridM[[1]], detector = "count")
    trapsBoth <- read.traps(data = gridBoth[[1]], detector = "count")
    
    Proposed.Traps.Grid[[1]] <- trapsF
    Proposed.Traps.Grid[[2]] <- trapsM
    Proposed.Traps.Grid[[3]] <- trapsBoth
    
    #use secrdesign to run simulations comparing the designs
    scen <- make.scenarios (trapsindex = 1:6, detectfn = 'HHN', D = DF, lambda0 = lambda0F, 
                            sigma = sigmaF, noccasions = 1, groups = c('F','M'))
    
    male <- scen$group == 'M'
    scen$D[male] <- DM
    scen$lambda0[male] <- lambda0M
    scen$sigma[male] <- sigmaM
    
    traps.list <- list("Enr (F)" = Proposed.Traps.Enr[[1]], "Enr (M)" = Proposed.Traps.Enr[[2]], "Enr (both)" = Proposed.Traps.Enr[[3]],"Grid (F)" = Proposed.Traps.Grid[[1]],"Grid (M)" = Proposed.Traps.Grid[[2]], "Grid (avg)" = Proposed.Traps.Grid[[3]])
    mask.list <- list("Mask (F)" = mask.buffF[[1]], "Mask (M)" = mask.buffM[[1]], "Mask (M)" = mask.buffM[[1]], "Mask (F)" = mask.buffF[[1]], "Mask (M)" = mask.buffM[[1]], "Mask (M)" = mask.buffM[[1]])
    
    #if raw = T will keep raw datasets ; extractfn = NULL defaults to reporting 4 summary stats
    data.stats <- NULL ; raw.data <- NULL
    if (raw==TRUE){
      data.stats <- run.scenarios(scen, trapset = traps.list, maskset = mask.list, nrepl = nreps, extractfn = NULL, fit = FALSE)
      raw.data <- run.scenarios(scen, trapset = traps.list, maskset = mask.list, nrepl = nreps, extractfn = identity, fit = FALSE)
      sims.group.both <- fit.models(raw.data, fit.args = list(detectfn = 'HHN', model = list(D~g, sigma~g, lambda0~g), groups = "group"),
                                    maskset = mask.list, fit = TRUE, ncores = 4, byscenario = FALSE)

    } else {
      sims.group.both <- run.scenarios(nreps, scen, trapset = traps.list, fit = TRUE, extractfn = predict,
                                       fit.args = list(detectfn = 'HHN', model = list(D~g, sigma~g, lambda0~g), groups = "group"),
                                       maskset = mask.list, ncores = 4, byscenario = FALSE)
    }
    return(list("Data summary stats" = data.stats, "CHs" = raw.data, "Sim results" = sims.group.both))
  }
  
  #test
  alltrpsF <- make.grid(nx = 119, ny = 119, origin = c(1000,1000), spacing = 2/3 *sigmaF, detector = "count")
  alltrpsM <- make.grid(nx = 119, ny = 119, origin = c(1000,1000), spacing = 2/3 *sigmaM, detector = "count")
  alltrpsBoth <- make.grid(nx = 119, ny = 119, origin = c(1000,1000), spacing = mean(c(2/3 *sigmaF, 2/3 * sigmaM)), detector = "count")
  
  test <- Sims1(alltraps = alltrpsF, nT = 100, sigmaF = sigmaF, sigmaM = sigmaM, 
                lambda0F = lambda0F, lambda0M = lambda0M, DF = DF, DM = DM, nreps = 1, raw = FALSE, ngen = 3 )
  
  
  estimateSummary(test, 'D',  format = 'data.frame', cols = c(1,2,6:8))
  estimateSummary(test, 'lambda0',  format = 'data.frame',cols = c(1,2,6:8))
  estimateSummary(test, 'sigma',  format = 'data.frame',cols = c(1,2,6:8))
  
  #########################
  #for new function that returns data summaries, need to go:
  estimateSummary(sim.results1.M[[3]], 'D',  format = 'data.frame', cols = c(1,2,6:8))
  estimateSummary(sim.results1.M[[3]], 'lambda0',  format = 'data.frame', cols = c(1,2,6:8))
  estimateSummary(sim.results1.M[[3]], 'sigma',  format = 'data.frame', cols = c(1,2,6:8))
  
  ######################
  #optimised grid
  ######################
  grid_traps <- read.traps(data = gridF[[1]], detector = dt)
  
  # compute optimal spacing for previous grid 
  dd <- optimalSpacing(D = DF, traps = grid_traps, detectpar = list(lambda0 = lambda0F, sigma = sigmaF), noccasions = 1)$rotRSE$optimum.spacing
  
  
 


#####bits of code
red_factor <- cellsize[1] / attr(msk, "spacing")
if ((trunc(red_factor) - red_factor) != 0) stop("Check spacing, causing non-integer reduction factor for mesh")

#mask and msk both the same thing? mask objects
mask <- secr::raster(msk)
mask_df <- data.frame(coordinates(mask))
mask <- read.mask(data = mask_df)
plot(mask, dots=F)

plot(msk)
plot(opt$optimaltraps, add = TRUE)
minnrRSE(opt, distribution = 'binomial')

buffer_mult <- ceiling(bufferF / cellsizeF)
newmask_df <- expand.grid(x = seq(from = min(mask_df$x) - buffer_mult * cellsizeF, 
                                  to = max(mask_df$x) + buffer_mult * cellsizeF, 
                                  by = cellsizeF),
                          y = seq(from = min(mask_df$y) - buffer_mult * cellsizeF, 
                                  to = max(mask_df$y) + buffer_mult * cellsizeF, 
                                  by = cellsizeF)) 
newmask_df <- newmask_df %>% 
  mutate(keep = (apply(e2dist(as.matrix(newmask_df), as.matrix(mask_df)), 1, min) <= bufferF)) %>%
  filter(keep == TRUE) %>% dplyr::select(-keep)
newmask <- read.mask(data = newmask_df)
plot(newmask,axes=T)

# create grid of possible trap locations
cellsize <- 2/3 * sigma
red_factor <- cellsize[1] / attr(msk, "spacing")
if ((trunc(red_factor) - red_factor) != 0) stop("Check spacing, causing non-integer reduction factor for mesh")

alltraps_df <- data.frame(coordinates(data.frame(alltrps)))
alltraps <- as.matrix(alltraps_df)[,c(1,2)]
plot(read.mask(data = alltraps_df), col = 'blue', add = TRUE, axes =T)  


mask.raster <- secr::raster(mask)
mask_df <- data.frame(coordinates(mask.raster))
mask <- read.mask(data = mask_df)

# create grid of possible trap locations
alltraps_df <- data.frame(coordinates(data.frame(alltraps)))
alltraps <- as.matrix(alltraps_df)
