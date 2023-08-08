#modifying OF for design, Aug 2023
#note that GAoptim calls kofnGA which uses an OF that matches the 1:6 criteria or allows a fn
#removed penalty for now
#criterion for two pops is crit 5

# Objective function
OF2 <- function (v, alltraps, mask, detectpar, detectfn, noccasions, 
                detector, D, crit, weights = c(0.5,0.5)) {
  
  traps <- subset(alltraps, v)
  
  strata = FALSE
  if (length(detectpar[[2]])>1) strata = TRUE
  
  if (is.function(crit)) {
    -crit(D = D, traps = traps, mask = mask, noccasions = noccasions, 
          detectpar = detectpar, detectfn = detectfn)[1]
  }
  if (strata==F){
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
      
      obj.n <- -(weights[1]*enrm1[1] + weights[2]*enrm2[1])
      obj.r <- -(weights[1]*enrm1[2] + weights[2]*enrm2[2])
      obj.m <- -(weights[1]*enrm1[3] + weights[2]*enrm2[3])
      obj.min <- -(weights[1]*(min(enrm1[1],enrm1[2]) + weights[2]*min(enrm2[1],enrm2[2])))
      
      c(obj.n, obj.r, obj.m, obj.min)[crit]    
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


#-------------------------------------------------------------------------------

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
                        OF = OF2,
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
    optimalenrm <- Enrm(D = D, traps = optimaltraps, mask = mask, 
                        noccasions = noccasions, detectpar = detectpar, detectfn = detectfn)
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

# user parameters
lambda0 <- 2 
dens_per_100km2 <- 2 # mean animal density per 100km2, SLs are ~1
D <- dens_per_100km2 / 10000
sigma <- 3000

nT <- 20
buffer <- 4 * sigma
dt <- "count"

# an artificial example, note the bottom-left corner is at origin
msk <- make.mask(type = 'rectangular', spacing = 1000, nx = 30, ny = 20, buffer = 0)
alltrps <- make.grid(nx = 29, ny = 19, origin = c(1000,1000), spacing = 1000, detector = "count")

# 50 generations for demonstration, use more in practice
opt <- GAoptimTest(mask = msk, alltraps = alltrps, ntraps = 20, detectpar = list(lambda0 = lambda0, sigma = sigma), 
               detectfn = 'HHN', D = D, noccasions = 1, ngen = 5, verbose = 1)

plot(msk)
plot(opt$optimaltraps, add = TRUE)

#now for two pops
lambda0 <- c(2,1.75) 
dens_per_100km2 <- c(0.0002,0.0001) # mean animal density per 100km2, SLs are ~1
D <- dens_per_100km2 / 10000
sigma <- c(3000, 5000)
nT <- 20
buffer <- 4 * max(sigma)
dt <- "count"

opt2a <- GAoptimTest(mask = msk, alltraps = alltrps, ntraps = 20, detectpar = list(lambda0 = lambda0, sigma = sigma), 
                   detectfn = 'HHN', D = D, weights = c(0.95,0.05), noccasions = 1, ngen = 100, verbose = 1)

opt2b <- GAoptimTest(mask = msk, alltraps = alltrps, ntraps = 20, detectpar = list(lambda0 = lambda0, sigma = sigma), 
                    detectfn = 'HHN', D = D, weights = c(0.05,0.95), noccasions = 1, ngen = 100, verbose = 1)

plot(msk)
plot(opt2a$optimaltraps, add = TRUE)
plot(opt2b$optimaltraps, add = TRUE)

####################################################################################
#Now start exploring scenarios to see how simulations perform



