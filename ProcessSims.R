#Aug 2023, process results of sims run on cluster
library(secr) ; library(secrdesign)
load("Cluster/Sims1.RData")

#for new function that returns data summaries, need to go:
estimateSummary(sim.results1.M[[3]], 'D',  format = 'data.frame', cols = c(1,2,6:8))
estimateSummary(sim.results1.M[[3]], 'lambda0',  format = 'data.frame', cols = c(1,2,6:8))
estimateSummary(sim.results1.M[[3]], 'sigma',  format = 'data.frame', cols = c(1,2,6:8))

#domensions for array [parameters, stats, Group, Scen, rep]
array1 <- estimateArray(sim.results1.M[[3]])

D.est.F.1 <- array1[1,1,1,1,]
D.est.M.1 <- array1[1,1,2,1,]

(mean(D.est.F.1[D.est.F.1<0.001]) - 0.00015)/0.00015

##########################
load("Cluster/Sims1b.RData")

#explore raw datasets

summary(sim.results1.M[[1]]$output[[1]])
