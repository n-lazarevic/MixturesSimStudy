##------------------------------------------------------------------------------
# Code for: "Performance of variable and function selection methods for
#            estimating the non-linear health effects of correlated
#            chemical mixtures: a simulation study."
# Authors:   Nina Lazarevic,  Luke D. Knibbs,  Peter D. Sly,  Adrian G. Barnett
#
# Written by Nina Lazarevic using R 3.4.3
#
##------------------------------------------------------------------------------
#
# This function runs the selected method on the simulated data.
#
# ***internal function,  see master.R for master script***
#
##------------------------------------------------------------------------------

method_estimation <- function(
  sim.data, select_method, save_model_obj = TRUE, num.cores = 2) {
  
  #load libraries if not loaded
  require(parallel)     #version 3.4.3
  require(doParallel)   #version 1.0.11
  if (select_method == "bkmr") require(bkmr)         #version 0.2.0
  if (select_method == "bart") require(bartMachine)  #version 1.2.3
  if (select_method == "star") require(spikeSlabGAM) #version 1.1 - 14
  if (select_method == "lasso") require(glmnet)      #version 2.0 - 13
  if (select_method == "gam") require(mgcv)          #version 1.8 - 23
  
  #simulated data
  list2env(sim.data, envir = environment())
  
  #load global constants into current environment
  source("global_constants.R", local = TRUE)

  #column indices for dataset d,  i.e.,  data for one rep as a data.frame
  d.y.ind <- c(1:N.scenarios) #column indices for y in dataset d
  d.x.ind <- c((N.scenarios + 1):(N.scenarios + K.m2)) #col indices for x in d
  
  ##two result data frames for each method,  named "res.corr_struct.method":
  #
  #res.obscorr.method (observed correlation datasets)
  #res.halfcorr.method (half observed correlation datasets)
  #
  #16 columns each,  composed of:
  #cols 1:8       models of 4 exposures associated with Y, out of K=6 exposures
  #   cols 1:4	    low SNR, corresponding to R2 = 0.1	
  #   cols 5:8	    high SNR, corresponding to R2 = 0.3
  #cols 9:16	    models of 4 exposures associated with Y, out of K=12 exposures
  #   cols 9:12	    low SNR, corresponding to R2 = 0.1	
  #   cols 13:16	  high SNR, corresponding to R2 = 0.3
  #then each set of 4 corresponds to 4 exposure-response functions
  #cols 1, 5, 9, 13    linearly increasing
  #cols 2, 6, 10, 14   log-logistic, S-shaped
  #cols 3, 7, 11, 15   quadratic, symmetric inv-U-shaped
  #cols 4, 8, 12, 16   asymmetric inv-U-shaped
  
  for (corr_struct in c("obscorr", "halfcorr")) {
    
    message(paste(toupper(select_method), ": Started working on ", 
                  toupper(corr_struct), " at ", Sys.time(), ".", sep = ""))
    
    if (select_method != "gam") {
      
      #initialise res,  results list
      init.mat <- matrix(NaN, reps, N.scenarios, 
                         dimnames = list(NULL, scenario.names))
      init.mat.logical <- matrix(NA, reps, N.scenarios, 
                                 dimnames = list(NULL, scenario.names))
      init.arr1 <- array(NaN, dim = c(reps, K.m1, N.dgp), 
                         dimnames = list(NULL, x.names.m1, scenario.names.m1))
      init.arr2 <- array(NaN, dim = c(reps, K.m2, N.dgp), 
                         dimnames = list(NULL, x.names.m2, scenario.names.m2))
      
      res <- list(
        m1 = list(pips = init.arr1, sel = init.arr1), 
        m2 = list(pips = init.arr2, sel = init.arr2), 
        rank_correct = init.mat.logical, 
        rank_correct_gr = init.mat.logical, 
        prop_rank_correct = init.mat, 
        prop_rank_correct_gr = init.mat, 
        sum.sel = init.mat, sens = init.mat, spec = init.mat, 
        FDR = init.mat, prec = init.mat, NPV = init.mat, FPR = init.mat, 
        TP = init.mat, FP = init.mat, TN = init.mat, FN = init.mat, 
        F1 = init.mat)
      
      #include method-specific variables in res$m1 and res$m2
      if (select_method == "bkmr") {
      } else if (select_method == "bart") {
        res$m1$sel.global_max <- init.arr1
        res$m1$sel.global_se <- init.arr1
        res$m2$sel.global_max <- init.arr2
        res$m2$sel.global_se <- init.arr2
      } else if (select_method == "star") {
        res$m1$p1.lin <- init.arr1
        res$m1$p1.sm <- init.arr1
        res$m2$p1.lin <- init.arr2
        res$m2$p1.sm <- init.arr2
      } else if (select_method == "lasso") {
        res$m1$coefs <- init.arr1
        res$m2$coefs <- init.arr2
      }
      
      rm(init.mat, init.mat.logical, init.arr1, init.arr2)
      
    }
    
    if (select_method != "lasso") {
      
      if (select_method == "gam") res <- list()
      
      #include result lists for est exp-resp curve at every 10th and 25th %iles
      res$erc.10 <- array(
        NaN, dim = c(reps, 11, K.true, N.scenarios, 3), 
        dimnames = list(NULL, 
                        c("min", "10%", "20%", "30%", "40%", "50%", 
                                "60%", "70%", "80%", "90%", "max"), 
                        x.true.names, 
                        scenario.names, 
                        c("point", "loCI", "hiCI")))
      res$erc.25 <- array(
        NaN, dim = c(reps, 3, K.true, N.scenarios, 3), 
        dimnames = list(NULL, 
                        c("25%", "50%", "75%"), 
                        x.true.names, 
                        scenario.names, 
                        c("point", "loCI", "hiCI")))
      
    }
    
    #loop over reps
    for (i in 1:reps) {
      
      message(paste(toupper(select_method), ": Started rep ", i, " of ", reps, 
                    " at ", Sys.time(), ".", sep = ""))
      
      seed = as.integer(seed.method_est_start + i)
      set.seed(seed)
      
      #data for one rep as data.frame
      if (corr_struct == "obscorr") d <- as.data.frame(
        cbind(simY[, , i], simY[, , i], simX[, , i]))
      if (corr_struct == "halfcorr") d <- as.data.frame(
        cbind(simY.halfcorr[, , i], simY.halfcorr[, , i], simX.halfcorr[, , i]))
      
      #new data points at specified quantiles,  for prediction
      probs = sort(c(seq(from = 0, to = 1, by = 0.1), 
                     0.25, 0.75, 0.025, 0.975, 0.05, 0.95))
      d.new <- list()
      for (j in 1:N.scenarios) {
        d.new[[j]] <- array(NaN, 
                            dim = c(length(probs), K.scenario[j], K.true), 
                            dimnames = list(paste(probs*100, "%", sep = ""), 
                                            x.names[x.inc[[j]]], x.true.names))
        for (name.x in x.true.names) {
          d.new[[j]][, , name.x] <- t(matrix(
            rep(colMeans(d[, d.x.ind[x.inc[[j]]]]), 
                length(probs)), K.scenario[j], length(probs)))
          d.new[[j]][, name.x, name.x] <- quantile(d[, name.x], probs = probs)
        }
      }
      
      #initiate results class
      resClass <- function(vs = NULL, model = NULL, erc = NULL) {
        comb.res.list <- list(vs = vs, model = model, erc = erc)
        class(comb.res.list) <- append(class(comb.res.list), "resClass")
        return(comb.res.list)
      }
      
      #run R in parallel for estimation 
      cl <- makeCluster(num.cores, type = "PSOCK")
      cl.last.used <<- cl
      registerDoParallel(cl)
      
      #estimation
      if (select_method == "bkmr") {
        
        #scenario - specific tuning parameters
        r.jump2 <- rep(0.2, N.scenarios)
        r.jump2[c(5:6, 13:14)] <- 0.1
        lambda.jump <- rep(1, N.scenarios)
        lambda.jump[c(1:4, 9:12)] <- 0.5
        
        #BKMR: estimate models for rep i
        comb.res <- foreach(j = 1:N.scenarios, .packages = "bkmr") %dopar% {
          
          set.seed(seed)
          
          m.bkmr <- tryCatch(kmbayes(
            y = d[, j], Z = d[, d.x.ind[x.inc[[j]]]], 
            iter = iter.bkmr, verbose = FALSE, varsel = TRUE, est.h = T, 
            starting.values = list(r = 0.1, lambda = 1), 
            control.params = list(
              a.p0 = 1, b.p0 = 1, a.sigsq = 0.001, b.sigsq = 0.001, 
              r.prior = "gamma", mu.r = 0.25, sigma.r = 0.25, 
              r.muprop = 0.1, r.jump1 = 0.1, r.jump2 = r.jump2[j], 
              mu.lambda = 1, sigma.lambda = 1, lambda.jump = lambda.jump[j])),  
            error = function(e) e)
          
          if (inherits(m.bkmr, "error")) {
            
            message(paste(
              toupper(select_method), ": Rep ", i, " of ", reps, ",  scenario ",
              j, ",  in ", toupper(corr_struct), 
              " produced the following error at ", Sys.time(), ":", sep = ""))
            message(m.bkmr$message)
            
            result <- resClass()
            result$model <- m.bkmr
            
          } else {
            
            result <- resClass()
            result$vs <- ExtractPIPs(m.bkmr)
            result$model <- m.bkmr
            
            tmperc <- array(NaN, dim = c(length(probs), 3, K.true),
                            dimnames = list(paste(probs * 100, "%", sep = ""),
                                      c("point", "loCI", "hiCI"), x.true.names))
            for (name.x in x.true.names) {
              tmppost <- SamplePred(m.bkmr, 
                                    Znew = d.new[[j]][, , name.x], 
                                    Xnew = cbind(0), 
                                    sel = c((iter.burnin.bkmr + 1):(iter.bkmr)))
              tmperc[, , name.x] <- cbind(
                colMeans(tmppost), 
                apply(tmppost, 2, quantile, probs = 0.05), 
                apply(tmppost, 2, quantile, probs = 0.95))
            }
            result$erc <- tmperc
            
          }
          
          return(result)
          
        }
        
      } else if (select_method == "bart") {
        
        #BART: estimate models for rep i
        comb.res <- foreach(j = 1:N.scenarios, .packages = "bartMachine") %dopar% {
          m.bart = bartMachine(
            X = d[, d.x.ind[x.inc[[j]]]], y = d[, j],
            seed = seed, num_burn_in = iter.burnin.bart,
            num_iterations_after_burn_in = iter.bart - iter.burnin.bart,
            num_trees = 50, alpha = 0.95, beta = 2, k = 2, q = 0.9, nu = 3,
            verbose = F, serialize = T, mh_prob_steps = c(0.2, 0.6, 0.2)
          )
          
          result <- resClass()
          result$vs <- var_selection_by_permute(
            m.bart, num_permute_samples = 100, plot = FALSE)
          result$model <- m.bart
          
          tmperc <- array(NaN, dim = c(length(probs), 3, K.true),
                          dimnames = list(paste(probs * 100, "%", sep = ""),
                                     c("point", "loCI", "hiCI"), x.true.names))
          for (name.x in x.true.names) {
            tmppost <- bart_machine_get_posterior(
              m.bart, 
              as.data.frame(d.new[[j]][, , name.x]))$y_hat_posterior_samples
            tmperc[, , name.x] <- cbind(
              rowMeans(tmppost),
              apply(tmppost, 1, quantile, probs = 0.05),
              apply(tmppost, 1, quantile, probs = 0.95)
            )
          }
          result$erc <- tmperc
          
          return(result)
          
        }
        
      } else if (select_method == "star") {
        
        #STAR: estimate models for rep i
        comb.res <- foreach(j = 1:N.scenarios, .packages = "spikeSlabGAM") %dopar% {
          
          set.seed(seed)
          
          fmla <- as.formula(paste(
            colnames(d)[j], "~", paste(x.names[x.inc[[j]]], collapse = " + ")))
          m.star <- spikeSlabGAM(
            formula = fmla, data = d, 
            start = list(seed = as.integer(seed)), 
            mcmc = list(burnin = iter.burnin.star, nChains = 5, 
                        chainLength = iter.star, thin = 1, reduceRet = TRUE), 
            hyperparameters = list(w = c(1, 1), tau2 = c(5, 40), 
                                   gamma = c(0.025), sigma2 = c(0.001, 0.001)))
          
          result <- resClass()
          result$vs <- m.star$postMeans$pV1
          result$model <- m.star
          
          tmperc <- array(NaN, dim = c(length(probs), 3, K.true), 
                          dimnames = list(paste(probs*100, "%", sep = ""), 
                                     c("point", "loCI", "hiCI"), x.true.names))
          for (name.x in x.true.names) tmperc[, , name.x] <- predict(
            m.star, newdata = d.new[[j]][, , name.x], quantiles = c(0.05, 0.95))
          result$erc <- tmperc
          
          return(result)
          
        }
        
      } else if (select_method == "lasso") {
        
        #LASSO: estimate models for rep i
        comb.res <- foreach(j = 1:N.scenarios, .packages = "glmnet") %dopar% {
          
          set.seed(seed)
          
          m.lasso <- cv.glmnet(x = as.matrix(d[, d.x.ind[x.inc[[j]]]]), 
                               y = d[, j], nfolds = 10)
          
          result <- resClass()
          result$model <- m.lasso
          result$vs <- coef(m.lasso, s = m.lasso$lambda.min)[-1] != 0
          
          return(result)
          
        }
        
      } else if (select_method == "gam") {
        
        #GAM: estimate models for rep i
        comb.res <- foreach(j = 1:N.scenarios, .packages = "mgcv") %dopar% {
          
          set.seed(seed)
          
          fmla <- as.formula(paste(
            colnames(d)[j], "~s(", 
            paste(x.true.names, collapse = ") + s("), ")", sep = ""))
          m.gam <- gam(formula = fmla, family = gaussian(link = identity), 
                       data = d, method = "REML")
          
          result <- resClass()
          result$model <- m.gam
          
          tmperc <- array(NaN, dim = c(length(probs), 3, K.true), 
                          dimnames = list(paste(probs*100, "%", sep = ""), 
                                     c("point", "loCI", "hiCI"), x.true.names))
          for (name.x in x.true.names) {
            tmppredict <- predict(
              m.gam, newdata = as.data.frame(d.new[[j]][, , name.x]), 
              se.fit = T)
            tmperc[, , name.x] <- cbind(
              tmppredict$fit, 
              tmppredict$fit + qnorm(0.05)*tmppredict$se.fit, 
              tmppredict$fit + qnorm(0.95)*tmppredict$se.fit)
          }
          result$erc <- tmperc
          
          return(result)
          
        }
        
      }
      
      #stop cluster
      stopCluster(cl)
      
      #record variable selection results
      for (j in 1:N.scenarios) {
        
        if (inherits(comb.res[[j]]$model, "error")) next
        
        if (m.type[j] == "m1") j.m <- j
        if (m.type[j] == "m2") j.m <- j - N.dgp
        
        if (select_method != "gam") {
          
          sel <- rep(FALSE, K.scenario[j])
          pips <- rep(NaN, K.scenario[j])
          
          if (select_method == "bkmr") {
            
            #BKMR: record variable selection results and 
            #      posterior inclusion probabilities
            sel <- comb.res[[j]]$vs$PIP > 0.5
            pips <- comb.res[[j]]$vs$PIP
            
          } else if (select_method == "bart") {
            
            #BART: record variable selection results and 
            #      variable inclusion proportions (as variable "pips")
            sel[comb.res[[j]]$vs$important_vars_local_col_nums] <- TRUE
            pips <- comb.res[[j]]$vs$var_true_props_avg[c(x.names)[x.inc[[j]]]]
            
            #keep selections using other cutoffs
            sel.global_max <- sel.global_se <- rep(FALSE, K.scenario[j])
            sel.global_max[
              comb.res[[j]]$vs$important_vars_global_max_col_nums] <- TRUE
            sel.global_se[
              comb.res[[j]]$vs$important_vars_global_se_col_nums] <- TRUE
            res[[m.type[j]]]$sel.global_max[i, , j.m] <- sel.global_max
            res[[m.type[j]]]$sel.global_se[i, , j.m] <- sel.global_se
            
          } else if (select_method == "star") {
            
            #STAR: record variable selection results and 
            #      posterior inclusion probabilities
            for (l in 1:K.scenario[j]) sel[l] <- any(
              comb.res[[j]]$model$postMeans$pV1[c((2*(l - 1) + 1):(l*2))] > 0.5)
            for (l in 1:K.scenario[j]) pips[l] <- sum(
              comb.res[[j]]$model$postMeans$pV1[c((2*(l - 1) + 1):(l*2))])
            res[[m.type[j]]]$p1.lin[i, , j.m] <- 
              comb.res[[j]]$model$postMeans$pV1[seq(1, K.scenario[j]*2, 2)]
            res[[m.type[j]]]$p1.sm[i, , j.m] <- 
              comb.res[[j]]$model$postMeans$pV1[seq(2, K.scenario[j]*2, 2)]
            
          } else if (select_method == "lasso") {
            
            #LASSO: record variable selection results and coefficients
            sel <- comb.res[[j]]$vs
            res[[m.type[j]]]$coefs[i, , j.m] <- coef(
              comb.res[[j]]$model, s = comb.res[[j]]$model$lambda.min)[-1]
            
          }
          
          res[[m.type[j]]]$sel[i, , j.m] <- sel
          
          if (select_method != "lasso") {
            res[[m.type[j]]]$pips[i, , j.m] <- pips
            res$rank_correct[i, j] <-
              min(pips[x.true[[j]]]) >= max(pips[!x.true[[j]]])
            res$rank_correct_gr[i, j] <-
              min(pips[x.true[[j]]]) > max(pips[!x.true[[j]]])
            res$prop_rank_correct[i, j] <-
              sum(pips[x.true[[j]]] >= max(pips[!x.true[[j]]])) / K.true
            res$prop_rank_correct_gr[i, j] <-
              sum(pips[x.true[[j]]] > max(pips[!x.true[[j]]])) / K.true
          }
          
          TP <- sel * x.true[[j]]
          FP <- sel * (!x.true[[j]])
          TN <- (!sel) * (!x.true[[j]])
          FN <- (!sel) * x.true[[j]]
          
          res$sum.sel[i, j] <- sum(sel)
          res$TP[i, j] <- sum(TP)
          res$FP[i, j] <- sum(FP)
          res$TN[i, j] <- sum(TN)
          res$FN[i, j] <- sum(FN)
          res$sens[i, j] <- res$TP[i, j] / sum(x.true[[j]])
          res$spec[i, j] <- res$TN[i, j] / sum(!x.true[[j]])
          res$FDR[i, j] <- res$FP[i, j] / sum(sel)
          res$prec[i, j] <- 1 - res$FDR[i, j]
          res$NPV[i, j] <- res$TN[i, j] / sum(!sel)
          res$F1[i, j] <-
            2*(res$prec[i, j]*res$sens[i, j]) / (res$prec[i, j] + res$sens[i, j])
          res$FPR[i, j] <-
            res$FP[i, j] / (res$FP[i, j] + res$TN[i, j])
          
        }
        
        if (select_method != "lasso") {
          
          #record exp-resp curve estimates at required quantiles
          for (name.x in x.true.names) {
            for (k in c("point", "loCI", "hiCI")) {
              res$erc.10[i, , name.x, j, k] <- comb.res[[j]]$erc[
                c("0%", "10%", "20%", "30%", "40%", "50%", 
                  "60%", "70%", "80%", "90%", "100%"), k, name.x]
              res$erc.25[i, , name.x, j, k] <- comb.res[[j]]$erc[
                c("25%", "50%", "75%"), k, name.x]
            }
          }
          
        }
        
      } #close j
      
      #save comb.res
      if (save_model_obj == TRUE) save(
        comb.res, file = paste("combres_", corr_struct, "_", 
                               select_method, "_", i, ".Rdata", sep = ""))
      rm(comb.res)
      
      
    } #close i
    
    
    assign(paste("res.", corr_struct, ".", select_method, sep = ""), 
           res, envir = .GlobalEnv)
    
    if (select_method != "lasso") {
      
      if (corr_struct == "obscorr") true.erc.25.tmp <- true.erc.25
      if (corr_struct == "halfcorr") true.erc.25.tmp <- true.erc.25.halfcorr
      
      #MSEs for one method and corrn structure
      mse.25 <- matrix(NaN, N.scenarios, 4*3)
      colnames(mse.25) <- c("MPB.25%", "MPB.50%", "MPB.75%", 
                            "BP3.25%", "BP3.50%", "BP3.75%", 
                            "PPB.25%", "PPB.50%", "PPB.75%", 
                            "BPA.25%", "BPA.50%", "BPA.75%")
      rownames(mse.25) <- scenario.names
      for (j in 1:N.scenarios) {
        mse.25[j, 1:3] <- colMeans(
          (res$erc.25[, , "MPB", j, "point"] - 
             true.erc.25.tmp[, , "MPB", erfn.names.scenario[j]])^2, na.rm = T)
        mse.25[j, 4:6] <- colMeans(
          (res$erc.25[, , "BP3", j, "point"] - 
             true.erc.25.tmp[, , "BP3", erfn.names.scenario[j]])^2, na.rm = T)
        mse.25[j, 7:9] <- colMeans(
          (res$erc.25[, , "PPB", j, "point"] - 
             true.erc.25.tmp[, , "PPB", erfn.names.scenario[j]])^2, na.rm = T)
        mse.25[j, 10:12] <- colMeans(
          (res$erc.25[, , "BPA", j, "point"] - 
             true.erc.25.tmp[, , "BPA", erfn.names.scenario[j]])^2, na.rm = T)
      }
      
      #CI coverage for one method and corrn structure,  as a proportion
      coverage.25 <- matrix(NaN, N.scenarios, 4*3)
      colnames(coverage.25) <- c("MPB.25%", "MPB.50%", "MPB.75%", 
                                 "BP3.25%", "BP3.50%", "BP3.75%", 
                                 "PPB.25%", "PPB.50%", "PPB.75%", 
                                 "BPA.25%", "BPA.50%", "BPA.75%")
      rownames(coverage.25) <- scenario.names
      for (j in 1:N.scenarios) {
        coverage.25[j, 1:3] <- colMeans(
          (res$erc.25[, , "MPB", j, "loCI"] <= 
             true.erc.25.tmp[, , "MPB", erfn.names.scenario[j]])&
            (res$erc.25[, , "MPB", j, "hiCI"] >= 
               true.erc.25.tmp[, , "MPB", erfn.names.scenario[j]]), na.rm = T)
        coverage.25[j, 4:6] <- colMeans(
          (res$erc.25[, , "BP3", j, "loCI"] <= 
             true.erc.25.tmp[, , "BP3", erfn.names.scenario[j]])&
            (res$erc.25[, , "BP3", j, "hiCI"] >= 
               true.erc.25.tmp[, , "BP3", erfn.names.scenario[j]]), na.rm = T)
        coverage.25[j, 7:9] <- colMeans(
          (res$erc.25[, , "PPB", j, "loCI"] <= 
             true.erc.25.tmp[, , "PPB", erfn.names.scenario[j]])&
            (res$erc.25[, , "PPB", j, "hiCI"] >= 
               true.erc.25.tmp[, , "PPB", erfn.names.scenario[j]]), na.rm = T)
        coverage.25[j, 10:12] <- colMeans(
          (res$erc.25[, , "BPA", j, "loCI"] <= 
             true.erc.25.tmp[, , "BPA", erfn.names.scenario[j]])&
            (res$erc.25[, , "BPA", j, "hiCI"] >= 
               true.erc.25.tmp[, , "BPA", erfn.names.scenario[j]]), na.rm = T)
      }  
      
      assign(paste("mse.25.", corr_struct, ".", select_method, sep = ""), 
             mse.25, envir = .GlobalEnv)
      assign(paste("coverage.25.", corr_struct, ".", select_method, sep = ""), 
             coverage.25, envir = .GlobalEnv)
      
    }
    
    #clean up
    rm(d, res)
    if (select_method != "gam") rm(TP, FP, TN, FN, sel)
    if (select_method != "lasso") rm(pips, mse.25, coverage.25, true.erc.25.tmp, k)
    if (select_method == "bkmr") rm(lambda.jump, r.jump2)
    if (select_method == "bart") rm(sel.global_max, sel.global_se)
    
  } #close corr_struct
  
  #clean up
  rm(i, j, j.m, seed, resClass)
  
  #save workspace for each method
  save.image(file = paste("workspace", "_", select_method, "_", 
                          format(Sys.Date(), "%Y%m%d"), ".RData", sep = ""))
  
}

##------------------------------------------------------------------------------