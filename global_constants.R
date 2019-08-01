##------------------------------------------------------------------------------
# Code for: "Performance of variable and function selection methods for
#            estimating the non-linear health effects of correlated
#            chemical mixtures: a simulation study."
# Authors:   Nina Lazarevic, Luke D. Knibbs, Peter D. Sly, Adrian G. Barnett
#
# Written by Nina Lazarevic using R 3.4.3
#
##------------------------------------------------------------------------------
#
# This script contains constants/static variables required by all functions.
#
# ***internal function, see master.R for master script***
#
##------------------------------------------------------------------------------

N <- 250 #sample size
reps <- 100 #simulation replications

#for each corr_struct:
# 8 DGPs = 2 SNR settings (lo,hi) * 4 exp-resp fns (l,s,q,a)
# 16 scenarios = estimate 2 models (of 6,12 exposures) * 8 DGPs
N.m <- 2 #number of models
N.snr <- 2 #number of signal-to-noise ratios
N.erfn <- 4 #number of exp-resp functions
N.dgp <- N.erfn * N.snr #number of DGPs
N.scenarios <- N.dgp * N.m #number of scenarios (per corr_struct)
N.corr<-2 #number of correlation structures (observed, half)

#exposure, dgp, erfn names
if (exists("simX")) {
  x.names <- dimnames(simX)[[2]]
} else if (exists("x")) {
  x.names <- dimnames(x)[[2]]
} else {
  warning("No data loaded!")
}
dgp.names<-c("loSNR.l","loSNR.s","loSNR.q","loSNR.a",
             "hiSNR.l","hiSNR.s","hiSNR.q","hiSNR.a")
erfn.names <- c("linear", "loglogS", "quadinvU", "asymminvU")

K <- 12 #total number of exposures

#outcome-associated exposures in each model
K.true <- 4 #number of true exposures
x.true.names <- c("MPB", "BP3", "PPB", "BPA")

#constants by model (m1, m2)
K.m1 <- 6 #number of exposures in model 1
K.m2 <- 12 #number of exposures in model 2
x.false.names.m1 <- c("MEP", "MHH")
x.inc.m1 <- x.names %in% c(x.true.names,x.false.names.m1) #included in model 1
x.inc.m2 <- rep(TRUE,K) #included exposures in model 2
x.names.m1 <- x.names[x.inc.m1] #names of exposures in model 1
x.names.m2 <- x.names[x.inc.m2] #names of exposures in model 2
scenario.names.m1 <- paste("m1", dgp.names, sep = ".")
scenario.names.m2 <- paste("m2", dgp.names, sep = ".")
x.true.m1 <- x.names.m1 %in% x.true.names #true exposures in model 1
x.true.m2 <- x.names.m2 %in% x.true.names #true exposures in model 2

#constants by scenario (length N.scenario)
m.type <- c(rep("m1", N.dgp), rep("m2", N.dgp)) #model type
K.scenario <- c(rep(K.m1, N.dgp), rep(K.m2, N.dgp)) #number exposures
x.inc <- c(rep(list(x.inc.m1), N.dgp), rep(list(x.inc.m2), N.dgp)) #included exp
x.true <- c(rep(list(x.true.m1), N.dgp), rep(list(x.true.m2), N.dgp)) #true exp
scenario.names <- c(scenario.names.m1, scenario.names.m2) #scenario names
erfn.names.scenario <- rep(erfn.names, N.erfn) #exp-resp fn names by scenario
x.true.ind <- list() #column indices for true x
for (j in 1:N.scenarios) {
  x.true.ind[[j]] <- match(x.true.names, x.names[x.inc[[j]]])
}

#additional constants for simulating data:
R2.loSNR <- 0.1
R2.hiSNR <- 0.3
#linear association strength: stronger=2 and weaker=1
b <- list(MPB = 2, BP3 = 2, PPB = 1, BPA = 1) 

#MCMC constants
iter.burnin.bkmr <- 8000
iter.bkmr <- 10000
iter.burnin.bart <- 4000
iter.bart <- 6000
iter.burnin.star <- 8000
iter.star <- 10000

#random number generation seeds
seed.sim_exposures <- 4936
seed.sim_outcomes_start <- 49360  #starting seeds, incremented by 1 
seed.method_est_start <- 52652    #for each subsequent rep

##------------------------------------------------------------------------------