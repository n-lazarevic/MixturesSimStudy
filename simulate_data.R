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
# This function simulates exposure and outcome data.
#
# ***internal function, see master.R for master script***
#
##------------------------------------------------------------------------------

simulate_data<-function(x) {
  
  #load libraries if not loaded
  require(copula)       #version 0.999-18
  require(copulaedas)   #version 1.4.2
  require(gsl)          #version 1.9-10.3
  
  #load global constants into current environment
  source("global_constants.R",local=TRUE)
  
  #create log-transformed and standardised data
  x.lnstd <- x
  for (i in 1:K) {
    x.lnstd[, i] <- (log(x[, i]) - mean(log(x[, i]))) / sd(log(x[, i]))
  }
  
  #----------------------------------------------------------------------------
  #simulate exposure data using multivariate t copula and kernel-smoothed
  #truncated empirical marginal distributions
  
  dtrunckernel <<- ftrunckernel
  
  #create multivarate distribution based on t copula
  tCop <- tCopula(
    param = cor(x.lnstd, method = "spearman")[lower.tri(cor(x.lnstd))],
    dim = K,
    dispstr = "un"
  )
  mvd <- mvdc(copula = tCop,
              margins = rep("trunckernel", K),
              paramMargins = list(
                list(
                  a = min(x.lnstd[, "ECP"]),
                  b = max(x.lnstd[, "ECP"]),
                  X = x.lnstd[, "ECP"],
                  h = bw.nrd(x.lnstd[, "ECP"])
                ),
                list(
                  a = min(x.lnstd[, "MEP"]),
                  b = max(x.lnstd[, "MEP"]),
                  X = x.lnstd[, "MEP"],
                  h = bw.nrd(x.lnstd[, "MEP"])
                ),
                list(
                  a = min(x.lnstd[, "MPB"]),
                  b = max(x.lnstd[, "MPB"]),
                  X = x.lnstd[, "MPB"],
                  h = bw.nrd(x.lnstd[, "MPB"])
                ),
                list(
                  a = min(x.lnstd[, "MOH"]),
                  b = max(x.lnstd[, "MOH"]),
                  X = x.lnstd[, "MOH"],
                  h = bw.nrd(x.lnstd[, "MOH"])
                ),
                list(
                  a = min(x.lnstd[, "MBP"]),
                  b = max(x.lnstd[, "MBP"]),
                  X = x.lnstd[, "MBP"],
                  h = bw.nrd(x.lnstd[, "MBP"])
                ),
                list(
                  a = min(x.lnstd[, "MHH"]),
                  b = max(x.lnstd[, "MHH"]),
                  X = x.lnstd[, "MHH"],
                  h = bw.nrd(x.lnstd[, "MHH"])
                ),
                list(
                  a = min(x.lnstd[, "BP3"]),
                  b = max(x.lnstd[, "BP3"]),
                  X = x.lnstd[, "BP3"],
                  h = bw.nrd(x.lnstd[, "BP3"])
                ),
                list(
                  a = min(x.lnstd[, "MZP"]),
                  b = max(x.lnstd[, "MZP"]),
                  X = x.lnstd[, "MZP"],
                  h = bw.nrd(x.lnstd[, "MZP"])
                ),
                list(
                  a = min(x.lnstd[, "PPB"]),
                  b = max(x.lnstd[, "PPB"]),
                  X = x.lnstd[, "PPB"],
                  h = bw.nrd(x.lnstd[, "PPB"])
                ),
                list(
                  a = min(x.lnstd[, "COP"]),
                  b = max(x.lnstd[, "COP"]),
                  X = x.lnstd[, "COP"],
                  h = bw.nrd(x.lnstd[, "COP"])
                ),
                list(
                  a = min(x.lnstd[, "MIB"]),
                  b = max(x.lnstd[, "MIB"]),
                  X = x.lnstd[, "MIB"],
                  h = bw.nrd(x.lnstd[, "MIB"])
                ),
                list(
                  a = min(x.lnstd[, "BPA"]),
                  b = max(x.lnstd[, "BPA"]),
                  X = x.lnstd[, "BPA"],
                  h = bw.nrd(x.lnstd[, "BPA"])
                )
              )
  )
  
  #simulate exposure data
  set.seed(seed.sim_exposures)
  simX <- array(NaN, dim = c(N, K, reps))
  for (i in 1:reps) simX[, , i] <- rMvdc(N, mvd)
  dimnames(simX)[[2]] <- colnames(x.lnstd)
  simXmat <- aperm(simX, c(1, 3, 2))
  dim(simXmat) <- c(N * reps, K)
  
  #create multivarate distribution based on t copula for "low correlation" data
  halfcorrmat <- 0.5 * cor(x.lnstd, method = "spearman")
  diag(halfcorrmat) <- 1
  halfcorrmat.lowertri <- halfcorrmat[lower.tri(cor(x.lnstd))]
  tCop.halfcorr <- tCopula(param = halfcorrmat.lowertri,
                           dim = K,
                           dispstr = "un")
  mvd.halfcorr <- mvdc(copula = tCop.halfcorr,
                       margins = rep("trunckernel", K),
                       paramMargins = list(
                         list(
                           a = min(x.lnstd[, "ECP"]),
                           b = max(x.lnstd[, "ECP"]),
                           X = x.lnstd[, "ECP"],
                           h = bw.nrd(x.lnstd[, "ECP"])
                         ),
                         list(
                           a = min(x.lnstd[, "MEP"]),
                           b = max(x.lnstd[, "MEP"]),
                           X = x.lnstd[, "MEP"],
                           h = bw.nrd(x.lnstd[, "MEP"])
                         ),
                         list(
                           a = min(x.lnstd[, "MPB"]),
                           b = max(x.lnstd[, "MPB"]),
                           X = x.lnstd[, "MPB"],
                           h = bw.nrd(x.lnstd[, "MPB"])
                         ),
                         list(
                           a = min(x.lnstd[, "MOH"]),
                           b = max(x.lnstd[, "MOH"]),
                           X = x.lnstd[, "MOH"],
                           h = bw.nrd(x.lnstd[, "MOH"])
                         ),
                         list(
                           a = min(x.lnstd[, "MBP"]),
                           b = max(x.lnstd[, "MBP"]),
                           X = x.lnstd[, "MBP"],
                           h = bw.nrd(x.lnstd[, "MBP"])
                         ),
                         list(
                           a = min(x.lnstd[, "MHH"]),
                           b = max(x.lnstd[, "MHH"]),
                           X = x.lnstd[, "MHH"],
                           h = bw.nrd(x.lnstd[, "MHH"])
                         ),
                         list(
                           a = min(x.lnstd[, "BP3"]),
                           b = max(x.lnstd[, "BP3"]),
                           X = x.lnstd[, "BP3"],
                           h = bw.nrd(x.lnstd[, "BP3"])
                         ),
                         list(
                           a = min(x.lnstd[, "MZP"]),
                           b = max(x.lnstd[, "MZP"]),
                           X = x.lnstd[, "MZP"],
                           h = bw.nrd(x.lnstd[, "MZP"])
                         ),
                         list(
                           a = min(x.lnstd[, "PPB"]),
                           b = max(x.lnstd[, "PPB"]),
                           X = x.lnstd[, "PPB"],
                           h = bw.nrd(x.lnstd[, "PPB"])
                         ),
                         list(
                           a = min(x.lnstd[, "COP"]),
                           b = max(x.lnstd[, "COP"]),
                           X = x.lnstd[, "COP"],
                           h = bw.nrd(x.lnstd[, "COP"])
                         ),
                         list(
                           a = min(x.lnstd[, "MIB"]),
                           b = max(x.lnstd[, "MIB"]),
                           X = x.lnstd[, "MIB"],
                           h = bw.nrd(x.lnstd[, "MIB"])
                         ),
                         list(
                           a = min(x.lnstd[, "BPA"]),
                           b = max(x.lnstd[, "BPA"]),
                           X = x.lnstd[, "BPA"],
                           h = bw.nrd(x.lnstd[, "BPA"])
                         )
                       )
  )
  
  #simulate "low correlation" exposure data
  set.seed(seed.sim_exposures)
  simX.halfcorr <- array(NaN, dim = c(N, K, reps))
  for (i in 1:reps) simX.halfcorr[, , i] <- rMvdc(N, mvd.halfcorr)
  dimnames(simX.halfcorr)[[2]] <- colnames(x.lnstd)
  dim(simX.halfcorr)
  simXmat.halfcorr <- aperm(simX.halfcorr, c(1, 3, 2))
  dim(simXmat.halfcorr) <- c(N * reps, K)
  
  rm(dtrunckernel, envir = .GlobalEnv)
  
  #----------------------------------------------------------------------------
  #simulate outcome data
  
  #simY has 8 columns, composed of:
  #cols 1:8       models of k=4 exposures associated with Y
  #   cols 1:4	    low SNR, corresponding to R2=0.1
  #   cols 5:8	    high SNR, corresponding to R2=0.3
  #then each set of 4 corresponds to 4 exposure-response functions
  #cols 1,5   linearly increasing
  #cols 2,6   log-logistic, S-shaped
  #cols 3,7   quadratic, symmetric inv-U-shaped
  #cols 4,8   asymmetric inv-U-shaped

  #initialise arrays to store simulated outcome
  simY <- array(NaN, dim = c(N, N.dgp, reps))
  simY.halfcorr <- simY
  dimnames(simY)[[2]] <-
    dimnames(simY.halfcorr)[[2]] <- 
    c("Y.l.k4.loSNR", "Y.s.k4.loSNR", "Y.q.k4.loSNR", "Y.a.k4.loSNR",
      "Y.l.k4.hiSNR", "Y.s.k4.hiSNR", "Y.q.k4.hiSNR", "Y.a.k4.hiSNR")

  #exposure-response functions
  #   for each exposure:
  #     - curves are anchored so that their lower y-limit is the same as the
  #       linear function
  #     - the scale of each function is adjusted so that the AUC is equal to the
  #       AUC of the linear function
  
  #linear function
  #  parameter choices:
  #  b=2 for MPB and BP3, b=1 for PPB and BPA
  #  vs vertically shifts the curve, only used to obtain AUC
  linfn <- function(x, b, vs = 0) (b * x + vs)
  
  #S-shaped function (log-logistic CDF) 
  #  b represents steepness, c and d are lower and upper y-limits, respectively,
  #   and e+min(x) is the ED50 (Ritz et al. 2015)
  #  parameter choices:
  #  c is anchored at the lower y-limit of the linear function above
  #  b is set to b=4, so that curve has the required "S" shape between x-limits
  #  e is set so that ED50==0
  #  d is chosen so that the AUC is = to the AUC of the linear function above
  #  vs vertically shifts the curve, only used to obtain AUC
  #  minx is the minimum of x
  loglogSfn <- function(x, b, c, d, e, minx, vs = 0)
      ((d + (c - d) / (1 + exp(b * (log(x - minx) - log(e))))) + vs)
  
  #inverse U-shaped quadratic function, with vertex (vx,vy)
  #  parameter choices:
  #  vx is set as the midpoint of the x-limits
  #  vy is chosen so that the AUC is = to the AUC of the linear function above
  #  a is set so that curve is anchored at the lower y-limit of the linear 
  #   function above, i.e., a=-1*(required range of y)/(half the range of x)^2 
  #  vs vertically shifts the curve, only used to obtain AUC
  quadinvUfn <- function(x, vx, vy, a, vs = 0) (vy + a * (x - vx) ^ 2 + vs)
  
  #asymmetric inverse U-shaped function, using the Dawson function, where 
  #  mfs*(range of y linear fn)/(max(dawson(x-min(x)))-min(dawson(x-min(x)))) is
  #  the scale factor, i.e., mfs*(range of f_linear(x))/(range of f_dawson(x))
  #  parameter choices:
  #  range_y is the required range for y, set to the range of linear fn above
  #  mfs is a multiplicative factor on the scale, used to obtain a curve with 
  #   AUC equal to the AUC of the linear function above 
  #  vs vertically shifts the curve; as the dawson function is always positive, 
  #   this is set to the lower y-limit of the linear function above
  #  minx is the minimum of x
  asymminvUfn <- function(x, range_y, minx, mfs = 1, vs, xorig)
    ((dawson(x - minx)) * mfs * 
       (range_y / (max(dawson(xorig - minx)) - min(dawson(xorig - minx)))) + vs)
  
  #function that returns R2 for specific standard deviation values
  R2fn <- function(a, expY, seed) {
    set.seed(seed)
    e = rnorm(N, sd = a)
    Y = expY + e
    R2 = 1 - sum(e ^ 2) / sum((Y - mean(Y)) ^ 2)
    return(R2)
  }
  
  #initialise arrays for exp-resp fns and curve evaluated at specified x stats
  x.p10.names <- c("0%", "10%", "20%", "30%", "40%", "50%",
                        "60%", "70%", "80%", "90%", "100%")
  x.p25.names <- c("25%", "50%", "75%")
  x.stat.names <- c("mean", x.p10.names, x.p25.names)
  x.true.stats <- array(
    NaN, 
    dim = c(reps, length(x.stat.names), K.true), 
    dimnames = list(c(), x.stat.names, x.true.names))
  erfn.stats <- array(
    NaN, 
    dim = c(reps, length(x.stat.names), K.true, N.erfn), 
    dimnames = list(c(), x.stat.names, x.true.names, erfn.names))
  true.erc.10 <- array(
    NaN, 
    dim = c(reps, length(x.p10.names), K.true, N.erfn), 
    dimnames = list(c(), x.p10.names, x.true.names, erfn.names))
  true.erc.10.halfcorr <- array(
    NaN, 
    dim = c(reps, length(x.p10.names), K.true, N.erfn), 
    dimnames = list(c(), x.p10.names, x.true.names, erfn.names))
  true.erc.25 <- array(
    NaN, 
    dim = c(reps, length(x.p25.names), K.true, N.erfn), 
    dimnames = list(c(), x.p25.names, x.true.names, erfn.names))
  true.erc.25.halfcorr <- array(
    NaN, 
    dim = c(reps, length(x.p25.names), K.true, N.erfn), 
    dimnames = list(c(), x.p25.names, x.true.names, erfn.names))
  erfn <- array(
    NaN, 
    dim = c(reps, N, K.true, N.erfn), 
    dimnames = list(c(), NULL, x.true.names, erfn.names))
  true.erc <- array(
    NaN, 
    dim = c(reps, N, K.true, N.erfn), 
    dimnames = list(c(), NULL, x.true.names, erfn.names))
  true.erc.halfcorr <- array(
    NaN, 
    dim = c(reps, N, K.true, N.erfn), 
    dimnames = list(c(), NULL, x.true.names, erfn.names))
  
  #adjustable parameters in exp-resp fns
  d <- vy <- mfs <- matrix(NaN, reps, K.true)
  colnames(d) <- colnames(vy) <- colnames(mfs) <- x.true.names
  
  #simulate y
  for (corr_struct in c("obscorr", "halfcorr")) {
    
    for (i in 1:reps) {
      
      if (i %in% seq(0, reps, 10))
        message(paste0("For correlation structure ",
          toupper(corr_struct), ", simulating data for replication ", i))
      
      seed = seed.sim_outcomes_start + i
      
      set.seed(seed)
      e <- rnorm(N)
      
      #assume only phenols and parabens truly associated with y
      if (corr_struct == "obscorr") {
        MPB <- simX[, "MPB", i]
        BP3 <- simX[, "BP3", i]
        PPB <- simX[, "PPB", i]
        BPA <- simX[, "BPA", i]
      } else if (corr_struct == "halfcorr") {
        MPB <- simX.halfcorr[, "MPB", i]
        BP3 <- simX.halfcorr[, "BP3", i]
        PPB <- simX.halfcorr[, "PPB", i]
        BPA <- simX.halfcorr[, "BPA", i]
      }
      
      #expected y, linear relationship
      MPB.l <- linfn(MPB, b = b$MPB)
      BP3.l <- linfn(BP3, b = b$BP3)
      PPB.l <- linfn(PPB, b = b$PPB)
      BPA.l <- linfn(BPA, b = b$BPA)
      expY.l.k4 <- MPB.l + BP3.l + PPB.l + BPA.l
      
      #simulated y, linear relationship
      Y.l.k4.loSNR <-
        expY.l.k4 + e * uniroot(function(x)
          R2fn(x, expY = expY.l.k4, seed) - R2.loSNR, c(0.1, 30))$root
      Y.l.k4.hiSNR <-
        expY.l.k4 + e * uniroot(function(x)
          R2fn(x, expY = expY.l.k4, seed) - R2.hiSNR, c(0.1, 30))$root

      #find AUC of linear function
      AUC.MPB.l <- integrate(linfn, lower = min(MPB), upper = max(MPB), 
                             b = b$MPB, vs = -1*min(MPB.l))$value
      AUC.BP3.l <- integrate(linfn, lower = min(BP3), upper = max(BP3), 
                             b = b$BP3, vs = -1*min(BP3.l))$value
      AUC.PPB.l <- integrate(linfn, lower = min(PPB), upper = max(PPB), 
                             b = b$PPB, vs = -1*min(PPB.l))$value
      AUC.BPA.l <- integrate(linfn, lower = min(BPA), upper = max(BPA), 
                             b = b$BPA, vs = -1*min(BPA.l))$value

      #find value of d argument required for S-shaped fn,  
      #so that AUC equal to linear function AUC
      d[i, "MPB"] <- uniroot(
        function(x) integrate(loglogSfn, lower = min(MPB), upper = max(MPB), 
                              b = 4, c = min(MPB.l), d = x, 
                              e = mean(MPB)-1*min(MPB), minx = min(MPB), 
                              vs = -1*min(MPB.l))$value - AUC.MPB.l, 
        interval = c(max(MPB.l)*0, max(MPB.l)*2))$root
      d[i, "BP3"] <- uniroot(
        function(x) integrate(loglogSfn, lower = min(BP3), upper = max(BP3), 
                              b = 4, c = min(BP3.l), d = x, 
                              e = mean(BP3)-1*min(BP3), minx = min(BP3), 
                              vs = -1*min(BP3.l))$value - AUC.BP3.l, 
        interval = c(max(BP3.l)*0, max(BP3.l)*2))$root
      d[i, "PPB"] <- uniroot(
        function(x) integrate(loglogSfn, lower = min(PPB), upper = max(PPB), 
                              b = 4, c = min(PPB.l), d = x, 
                              e = mean(PPB)-1*min(PPB), minx = min(PPB), 
                              vs = -1*min(PPB.l))$value - AUC.PPB.l, 
        interval = c(max(PPB.l)*0, max(PPB.l)*2))$root
      d[i, "BPA"] <- uniroot(
        function(x) integrate(loglogSfn, lower = min(BPA), upper = max(BPA), 
                              b = 4, c = min(BPA.l), d = x, 
                              e = -1*min(BPA), minx = min(BPA), 
                              vs = -1*min(BPA.l))$value - AUC.BPA.l, 
        interval = c(max(BPA.l)*0, max(BPA.l)*2))$root
      
      #expected y,  S-shaped relationship
      MPB.s <- loglogSfn(MPB, b = 4, c = min(MPB.l), d = d[i, "MPB"], 
                         e = -1*min(MPB), minx = min(MPB))
      BP3.s <- loglogSfn(BP3, b = 4, c = min(BP3.l), d = d[i, "BP3"], 
                         e = -1*min(BP3), minx = min(BP3))
      PPB.s <- loglogSfn(PPB, b = 4, c = min(PPB.l), d = d[i, "PPB"], 
                         e = -1*min(PPB), minx = min(PPB))
      BPA.s <- loglogSfn(BPA, b = 4, c = min(BPA.l), d = d[i, "BPA"], 
                         e = -1*min(BPA), minx = min(BPA))
      expY.s.k4 <- MPB.s + BP3.s + PPB.s + BPA.s
      
      #simulated y, S-shaped relationship
      Y.s.k4.loSNR <-
        expY.s.k4 + e * uniroot(function(x)
          R2fn(x, expY = expY.s.k4, seed) - R2.loSNR, c(0.01, 30))$root
      Y.s.k4.hiSNR <-
        expY.s.k4 + e * uniroot(function(x)
          R2fn(x, expY = expY.s.k4, seed) - R2.hiSNR, c(0.01, 30))$root
      
      #find value of vy argument required for inverse U-shaped quadratic fn, 
      #so that AUC equal to linear function AUC
      vy[i, "MPB"] <- uniroot(
        function(x) integrate(quadinvUfn,  lower = min(MPB), upper = max(MPB), 
                              vx = mean(c(min(MPB), max(MPB))), vy = x, 
                              a = (min(MPB.l) - x)/((max(MPB) - min(MPB))/2)^2, 
                              vs = -1*min(MPB.l))$value - AUC.MPB.l, 
        interval = c(0,  max(MPB.l)))$root
      vy[i, "BP3"] <- uniroot(
        function(x) integrate(quadinvUfn,  lower = min(BP3), upper = max(BP3), 
                              vx = mean(c(min(BP3), max(BP3))), vy = x, 
                              a = (min(BP3.l) - x)/((max(BP3) - min(BP3))/2)^2, 
                              vs = -1*min(BP3.l))$value - AUC.BP3.l, 
        interval = c(0, max(BP3.l)))$root  
      vy[i, "PPB"] <- uniroot(
        function(x) integrate(quadinvUfn, lower = min(PPB), upper = max(PPB), 
                              vx = mean(c(min(PPB), max(PPB))), vy = x, 
                              a = (min(PPB.l) - x)/((max(PPB) - min(PPB))/2)^2, 
                              vs = -1*min(PPB.l))$value - AUC.PPB.l, 
        interval = c(0, max(PPB.l)))$root
      vy[i, "BPA"] <- uniroot(
        function(x) integrate(quadinvUfn, lower = min(BPA), upper = max(BPA), 
                              vx = mean(c(min(BPA), max(BPA))), vy = x, 
                              a = (min(BPA.l) - x)/((max(BPA) - min(BPA))/2)^2, 
                              vs = -1*min(BPA.l))$value - AUC.BPA.l, 
        interval = c(0, max(BPA.l)))$root  
      
      #expected y, inverse U-shaped quadratic relationship
      MPB.q <- quadinvUfn(
        MPB, vx = mean(c(min(MPB), max(MPB))), 
        vy = vy[i, "MPB"], 
        a = (min(MPB.l) - vy[i, "MPB"])/((max(MPB) - min(MPB))/2)^2)
      BP3.q <- quadinvUfn(
        BP3, vx = mean(c(min(BP3), max(BP3))), vy = vy[i, "BP3"], 
        a = (min(BP3.l) - vy[i, "BP3"])/((max(BP3) - min(BP3))/2)^2)
      PPB.q <- quadinvUfn(
        PPB, vx = mean(c(min(PPB), max(PPB))), vy = vy[i, "PPB"], 
        a = (min(PPB.l) - vy[i, "PPB"])/((max(PPB) - min(PPB))/2)^2)
      BPA.q <- quadinvUfn(
        BPA, vx = mean(c(min(BPA), max(BPA))), vy = vy[i, "BPA"], 
        a = (min(BPA.l) - vy[i, "BPA"])/((max(BPA) - min(BPA))/2)^2)
      expY.q.k4 <- MPB.q + BP3.q + PPB.q + BPA.q
      
      #simulated y,  inverse U-shaped quadratic relationship
      Y.q.k4.loSNR <-
        expY.q.k4 + e * uniroot(function(x)
          R2fn(x, expY = expY.q.k4, seed) - R2.loSNR, c(0.1, 30))$root
      Y.q.k4.hiSNR <-
        expY.q.k4 + e * uniroot(function(x)
          R2fn(x, expY = expY.q.k4, seed) - R2.hiSNR, c(0.1, 30))$root
      
      #find value of multiplying factor on scale of asymmetric inverse U-shaped 
      #function, required so that AUC equal to linear function AUC
      mfs[i, "MPB"] <- uniroot(
        function(x) integrate(
          asymminvUfn, lower = min(MPB), upper = max(MPB), 
          range_y = (-1*min(MPB.l) + max(MPB.l)), minx = min(MPB), mfs = x, 
          vs = 0, xorig = MPB)$value - AUC.MPB.l, interval = c(0.01, 3))$root
      mfs[i, "BP3"] <- uniroot(
        function(x) integrate(
          asymminvUfn, lower = min(BP3), upper = max(BP3), 
          range_y = (-1*min(BP3.l) + max(BP3.l)), minx = min(BP3), mfs = x, 
          vs = 0, xorig = BP3)$value - AUC.BP3.l, interval = c(0.01, 3))$root
      mfs[i, "PPB"] <- uniroot(
        function(x) integrate(
          asymminvUfn, lower = min(PPB), upper = max(PPB), 
          range_y = (-1*min(PPB.l) + max(PPB.l)), minx = min(PPB), mfs = x, 
          vs = 0, xorig = PPB)$value - AUC.PPB.l, interval = c(0.01, 3))$root
      mfs[i, "BPA"] <- uniroot(
        function(x) integrate(
          asymminvUfn, lower = min(BPA), upper = max(BPA), 
          range_y = (-1*min(BPA.l) + max(BPA.l)), minx = min(BPA), mfs = x, 
          vs = 0, xorig = BPA)$value - AUC.BPA.l, interval = c(0.01, 3))$root
      
      #expected y, asymmetric inverse U-shaped relationship
      MPB.a <- asymminvUfn(MPB, range_y = (-1*min(MPB.l) + max(MPB.l)), 
                           minx = min(MPB), mfs = mfs[i, "MPB"], 
                           vs = min(MPB.l), xorig = MPB)
      BP3.a <- asymminvUfn(BP3, range_y = (-1*min(BP3.l) + max(BP3.l)), 
                           minx = min(BP3), mfs = mfs[i, "BP3"], 
                           vs = min(BP3.l), xorig = BP3)
      PPB.a <- asymminvUfn(PPB, range_y = (-1*min(PPB.l) + max(PPB.l)), 
                           minx = min(PPB), mfs = mfs[i, "PPB"], 
                           vs = min(PPB.l), xorig = PPB)
      BPA.a <- asymminvUfn(BPA, range_y = (-1*min(BPA.l) + max(BPA.l)), 
                           minx = min(BPA), mfs = mfs[i, "BPA"], 
                           vs = min(BPA.l), xorig = BPA)
      expY.a.k4 <- MPB.a + BP3.a + PPB.a + BPA.a
      
      #simulated y, asymmetric inverse U-shaped relationship
      Y.a.k4.loSNR <-
        expY.a.k4 + e * uniroot(function(x)
          R2fn(x, expY = expY.a.k4, seed) - R2.loSNR, c(0.1, 30))$root
      Y.a.k4.hiSNR <-
        expY.a.k4 + e * uniroot(function(x)
          R2fn(x, expY = expY.a.k4, seed) - R2.hiSNR, c(0.1, 30))$root
      
      #add to simY matrix
      if (corr_struct == "obscorr") {
        simY[, 1, i] <- Y.l.k4.loSNR
        simY[, 2, i] <- Y.s.k4.loSNR
        simY[, 3, i] <- Y.q.k4.loSNR
        simY[, 4, i] <- Y.a.k4.loSNR
        simY[, 5, i] <- Y.l.k4.hiSNR
        simY[, 6, i] <- Y.s.k4.hiSNR
        simY[, 7, i] <- Y.q.k4.hiSNR
        simY[, 8, i] <- Y.a.k4.hiSNR
      } else if (corr_struct == "halfcorr") {
        simY.halfcorr[, 1, i] <- Y.l.k4.loSNR
        simY.halfcorr[, 2, i] <- Y.s.k4.loSNR
        simY.halfcorr[, 3, i] <- Y.q.k4.loSNR
        simY.halfcorr[, 4, i] <- Y.a.k4.loSNR
        simY.halfcorr[, 5, i] <- Y.l.k4.hiSNR
        simY.halfcorr[, 6, i] <- Y.s.k4.hiSNR
        simY.halfcorr[, 7, i] <- Y.q.k4.hiSNR
        simY.halfcorr[, 8, i] <- Y.a.k4.hiSNR
      }
    
      #exposures evaluated at mean and percentiles
      p10 <- seq(from = 0, to = 1, by = 0.1)
      p25 <- c(0.25, 0.5, 0.75)
      x.true.stats[i, , "MPB"] <- c(mean(MPB), 
                                    quantile(MPB, probs = p10), 
                                    quantile(MPB, probs = p25))
      x.true.stats[i, , "BP3"] <- c(mean(BP3), 
                                    quantile(BP3, probs = p10), 
                                    quantile(BP3, probs = p25))
      x.true.stats[i, , "PPB"] <- c(mean(PPB), 
                                    quantile(PPB, probs = p10), 
                                    quantile(PPB, probs = p25))
      x.true.stats[i, , "BPA"] <- c(mean(BPA), 
                                    quantile(BPA, probs = p10), 
                                    quantile(BPA, probs = p25))
      
      #exposure-response fns evaluated at specified x percentiles and mean
      erfn.stats[i, , "MPB", "linear"] <- linfn(x.true.stats[i, , "MPB"], 
                                                b = b$MPB)
      erfn.stats[i, , "BP3", "linear"] <- linfn(x.true.stats[i, , "BP3"], 
                                                b = b$BP3)
      erfn.stats[i, , "PPB", "linear"] <- linfn(x.true.stats[i, , "PPB"], 
                                                b = b$PPB)
      erfn.stats[i, , "BPA", "linear"] <- linfn(x.true.stats[i, , "BPA"], 
                                                b = b$BPA)
      
      erfn.stats[i, , "MPB", "loglogS"] <- loglogSfn(
        x.true.stats[i, , "MPB"], b = 4, c = min(MPB.l), d = d[i, "MPB"], 
        e = -1*min(MPB), minx = min(MPB))
      erfn.stats[i, , "BP3", "loglogS"] <- loglogSfn(
        x.true.stats[i, , "BP3"], b = 4, c = min(BP3.l), d = d[i, "BP3"], 
        e = -1*min(BP3), minx = min(BP3))
      erfn.stats[i, , "PPB", "loglogS"] <- loglogSfn(
        x.true.stats[i, , "PPB"], b = 4, c = min(PPB.l), d = d[i, "PPB"], 
        e = -1*min(PPB), minx = min(PPB))
      erfn.stats[i, , "BPA", "loglogS"] <- loglogSfn(
        x.true.stats[i, , "BPA"], b = 4, c = min(BPA.l), d = d[i, "BPA"], 
        e = -1*min(BPA), minx = min(BPA))
      
      erfn.stats[i, , "MPB", "quadinvU"] <- quadinvUfn(
        x.true.stats[i, , "MPB"], vx = mean(c(min(MPB), max(MPB))), 
        vy = vy[i, "MPB"], 
        a = (min(MPB.l)-vy[i, "MPB"])/((max(MPB)-min(MPB))/2)^2)
      erfn.stats[i, , "BP3", "quadinvU"] <- quadinvUfn(
        x.true.stats[i, , "BP3"], vx = mean(c(min(BP3), max(BP3))), 
        vy = vy[i, "BP3"], 
        a = (min(BP3.l)-vy[i, "BP3"])/((max(BP3)-min(BP3))/2)^2)
      erfn.stats[i, , "PPB", "quadinvU"] <- quadinvUfn(
        x.true.stats[i, , "PPB"], vx = mean(c(min(PPB), max(PPB))), 
        vy = vy[i, "PPB"], 
        a = (min(PPB.l)-vy[i, "PPB"])/((max(PPB)-min(PPB))/2)^2)
      erfn.stats[i, , "BPA", "quadinvU"] <- quadinvUfn(
        x.true.stats[i, , "BPA"], vx = mean(c(min(BPA), max(BPA))), 
        vy = vy[i, "BPA"], 
        a = (min(BPA.l)-vy[i, "BPA"])/((max(BPA)-min(BPA))/2)^2)
      
      erfn.stats[i, , "MPB", "asymminvU"] <- asymminvUfn(
        x.true.stats[i, , "MPB"], range_y = (-1*min(MPB.l)+max(MPB.l)), 
        minx = min(MPB), mfs = mfs[i, "MPB"], vs = min(MPB.l), xorig = MPB)
      erfn.stats[i, , "BP3", "asymminvU"] <- asymminvUfn(
        x.true.stats[i, , "BP3"], range_y = (-1*min(BP3.l)+max(BP3.l)), 
        minx = min(BP3), mfs = mfs[i, "BP3"], vs = min(BP3.l), xorig = BP3)
      erfn.stats[i, , "PPB", "asymminvU"] <- asymminvUfn(
        x.true.stats[i, , "PPB"], range_y = (-1*min(PPB.l)+max(PPB.l)), 
        minx = min(PPB), mfs = mfs[i, "PPB"], vs = min(PPB.l), xorig = PPB)
      erfn.stats[i, , "BPA", "asymminvU"] <- asymminvUfn(
        x.true.stats[i, , "BPA"], range_y = (-1*min(BPA.l)+max(BPA.l)), 
        minx = min(BPA), mfs = mfs[i, "BPA"], vs = min(BPA.l), xorig = BPA)
      
      #exposure-response fns, stored in erfn
      erfn[i, , "MPB", "linear"] <- MPB.l
      erfn[i, , "BP3", "linear"] <- BP3.l
      erfn[i, , "PPB", "linear"] <- PPB.l
      erfn[i, , "BPA", "linear"] <- BPA.l
      
      erfn[i, , "MPB", "loglogS"] <- MPB.s
      erfn[i, , "BP3", "loglogS"] <- BP3.s
      erfn[i, , "PPB", "loglogS"] <- PPB.s
      erfn[i, , "BPA", "loglogS"] <- BPA.s
      
      erfn[i, , "MPB", "quadinvU"] <- MPB.q
      erfn[i, , "BP3", "quadinvU"] <- BP3.q
      erfn[i, , "PPB", "quadinvU"] <- PPB.q
      erfn[i, , "BPA", "quadinvU"] <- BPA.q
      
      erfn[i, , "MPB", "asymminvU"] <- MPB.a
      erfn[i, , "BP3", "asymminvU"] <- BP3.a
      erfn[i, , "PPB", "asymminvU"] <- PPB.a
      erfn[i, , "BPA", "asymminvU"] <- BPA.a
      
      if (i == reps) {
        if (corr_struct == "obscorr") {
          erfn.stats.obscorr <- erfn.stats
          x.true.stats.obscorr <- x.true.stats
          erfn.obscorr <- erfn
        } else if (corr_struct == "halfcorr") {
          erfn.stats.halfcorr <- erfn.stats
          x.true.stats.halfcorr <- x.true.stats
          erfn.halfcorr <- erfn
        }
      }
      
      #true exposure-response curve with one x evaluated at specified 
      #percentiles and others at mean
      for (fn in erfn.names) {
        for (x_j in x.true.names) {
          x_notj <- setdiff(x.true.names, x_j)
          if (corr_struct == "obscorr") {
            true.erc[i, , x_j, fn] <- erfn[i, , x_j, fn] + 
              sum(erfn.stats[i, 1, x_notj, fn])
            true.erc.10[i, , x_j, fn] <- erfn.stats[i, 2:12, x_j, fn] + 
              sum(erfn.stats[i, 1, x_notj, fn])
            true.erc.25[i, , x_j, fn] <- erfn.stats[i, 13:15, x_j, fn] + 
              sum(erfn.stats[i, 1, x_notj, fn])
          } else if (corr_struct == "halfcorr") {
            true.erc.halfcorr[i, , x_j, fn] <- erfn[i, , x_j, fn] + 
              sum(erfn.stats[i, 1, x_notj, fn])
            true.erc.10.halfcorr[i, , x_j, fn] <- erfn.stats[i, 2:12, x_j, fn] + 
              sum(erfn.stats[i, 1, x_notj, fn])
            true.erc.25.halfcorr[i, , x_j, fn] <- 
              erfn.stats[i, 13:15, x_j, fn] + sum(erfn.stats[i, 1, x_notj, fn])
          }
        }
      }
      
    } #close reps loop
  } #close corr_struct loop
  
  sim.data <- list(
    simX = simX,
    simX.halfcorr = simX.halfcorr,
    simY = simY,
    simY.halfcorr = simY.halfcorr,
    x.true.stats.obscorr = x.true.stats.obscorr,
    x.true.stats.halfcorr = x.true.stats.halfcorr,
    true.erc.10 = true.erc.10,
    true.erc.10.halfcorr = true.erc.10.halfcorr,
    true.erc.25 = true.erc.25,
    true.erc.25.halfcorr = true.erc.25.halfcorr
  )
  
  return(sim.data)
  
}

##------------------------------------------------------------------------------