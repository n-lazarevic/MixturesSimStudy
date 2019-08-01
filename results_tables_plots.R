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
# This function makes the following results tables and plots:
#
# MSE25_OBSCORR_<date>.csv & MSE25_HALFCORR_<date>.csv
# Coverage25_OBSCORR_<date>.csv & Coverage25_HALFCORR_<date>.csv
#   - these tables show MSEs and 90% credible interval coverage proportions
#     evaluated at the 25th, 50th, and 75th percentiles of an exposure
#   - 12 columns of 3 percentiles x 4 methods (BKMR, BART, STAR, GAM)
#   - 64 rows of 16 scenarios for each of 4 true exposures (MPB, BP3, PPB, BPA)
#
# VarSelection_<date>.csv
#   - this table contains 18 columns (9 statistics x 2 corr structures) and
#     64 rows of 16 scenarios x 4 methods
#
# Figures 3 to 10 as PDFs
#
# ***internal function,  see master.R for master script***
#
##------------------------------------------------------------------------------

results_tables_plots <- function(
  sim.data,
  res.halfcorr.bkmr, mse.25.halfcorr.bkmr, coverage.25.halfcorr.bkmr,
  res.obscorr.bkmr, mse.25.obscorr.bkmr, coverage.25.obscorr.bkmr,
  res.halfcorr.bart, mse.25.halfcorr.bart, coverage.25.halfcorr.bart,
  res.obscorr.bart, mse.25.obscorr.bart, coverage.25.obscorr.bart,
  res.halfcorr.star, mse.25.halfcorr.star, coverage.25.halfcorr.star,
  res.obscorr.star, mse.25.obscorr.star, coverage.25.obscorr.star,
  res.halfcorr.gam, mse.25.halfcorr.gam, coverage.25.halfcorr.gam,
  res.obscorr.gam, mse.25.obscorr.gam, coverage.25.obscorr.gam,
  res.halfcorr.lasso, 
  res.obscorr.lasso) {
  
  #load libraries if not loaded
  require(ggplot2)      #version 2.2.1,  for plotting only
  require(scales)       #version 0.5.0,  for plotting only
  require(cowplot)      #version 0.9.4,  for plotting only
  
  #simulated data
  list2env(sim.data, envir = environment())
  
  #load global constants into current environment
  source("global_constants.R", local = TRUE)
  
  #constants
  N.methods<-4
  N.ptiles<-3
  
  ##----------------------------------------------------------------------------
  message("Preparing tables of results...")
  #table of MSEs
  #array dims: N.scenarios, 3 percentiles x 4 methods, 4 true exposures
  tmp.mse.25 <- array(NaN, 
                      dim = c(N.scenarios, N.ptiles*N.methods, K.true), 
                      dimnames = list(scenario.names, c(
                        "BKMR.25%", "BART.25%", "STAR.25%", "GAM.25%", 
                        "BKMR.50%", "BART.50%", "STAR.50%", "GAM.50%", 
                        "BKMR.75%", "BART.75%", "STAR.75%", "GAM.75%"), 
                        x.true.names))
  for (corr_struct in c("obscorr", "halfcorr")) {
    counter = 0
    tmp.bkmr <- get(paste("mse.25.", corr_struct, ".bkmr", sep = ""))
    tmp.bart <- get(paste("mse.25.", corr_struct, ".bart", sep = ""))
    tmp.star <- get(paste("mse.25.", corr_struct, ".star", sep = ""))
    tmp.gam <- get(paste("mse.25.", corr_struct, ".gam", sep = ""))
    for (name.x in x.true.names) {
      counter = counter + 1
      tmp.mse.25[, , name.x] <- cbind(
        tmp.bkmr[, counter*3-2], tmp.bart[, counter*3-2], 
        tmp.star[, counter*3-2], tmp.gam[, counter*3-2], 
        tmp.bkmr[, counter*3-1], tmp.bart[, counter*3-1], 
        tmp.star[, counter*3-1], tmp.gam[, counter*3-1], 
        tmp.bkmr[, counter*3], tmp.bart[, counter*3], 
        tmp.star[, counter*3], tmp.gam[, counter*3])
    }
    if (corr_struct == "obscorr") mse.25.obscorr <- tmp.mse.25
    if (corr_struct == "halfcorr") mse.25.halfcorr <- tmp.mse.25
  }
  tmp.order <- c(seq(1, N.scenarios, N.erfn), 
                 seq(2, N.scenarios, N.erfn), 
                 seq(3, N.scenarios, N.erfn), 
                 seq(4, N.scenarios, N.erfn)) #all linear,  then s-shape,  etc.
  for (name.x in x.true.names) {
    write.table(mse.25.obscorr[tmp.order, , name.x], file = paste(
      "MSE25_OBSCORR_", format(Sys.Date(), "%Y%m%d"), ".csv", sep = ""), 
      sep = ",", append = T, col.names = NA)
    write.table(mse.25.halfcorr[tmp.order, , name.x], file = paste(
      "MSE25_HALFCORR_", format(Sys.Date(), "%Y%m%d"), ".csv", sep = ""), 
      sep = ",", append = T, col.names = NA)
  }
  rm(counter, name.x, tmp.bkmr, tmp.bart, tmp.star, tmp.mse.25, tmp.order)
  
  ##----------------------------------------------------------------------------
  #table of coverage proportions
  #array dims: N.scenarios, 3 percentiles x 4 methods, 4 true exposures
  tmp.coverage.25 <- array(
    NaN, 
    dim = c(N.scenarios, N.ptiles*N.methods, K.true), 
    dimnames = list(scenario.names, 
                    c("BKMR.25%", "BART.25%", "STAR.25%", "GAM.25%", 
                      "BKMR.50%", "BART.50%", "STAR.50%", "GAM.50%", 
                      "BKMR.75%", "BART.75%", "STAR.75%", "GAM.75%"), 
                    x.true.names))
  for (corr_struct in c("obscorr", "halfcorr")) {
    counter = 0
    tmp.bkmr <- get(paste("coverage.25.", corr_struct, ".bkmr", sep = ""))
    tmp.bart <- get(paste("coverage.25.", corr_struct, ".bart", sep = ""))
    tmp.star <- get(paste("coverage.25.", corr_struct, ".star", sep = ""))
    tmp.gam <- get(paste("coverage.25.", corr_struct, ".gam", sep = ""))
    for (name.x in x.true.names) {
      counter = counter + 1
      tmp.coverage.25[, , name.x] <- cbind(
        tmp.bkmr[, counter*3-2], tmp.bart[, counter*3-2], 
        tmp.star[, counter*3-2], tmp.gam[, counter*3-2], 
        tmp.bkmr[, counter*3-1], tmp.bart[, counter*3-1], 
        tmp.star[, counter*3-1], tmp.gam[, counter*3-1], 
        tmp.bkmr[, counter*3], tmp.bart[, counter*3], 
        tmp.star[, counter*3], tmp.gam[, counter*3])
    }
    if (corr_struct == "obscorr") coverage.25.propn.obscorr <- tmp.coverage.25
    if (corr_struct == "halfcorr") coverage.25.propn.halfcorr <- tmp.coverage.25
  }
  tmp.order <- c(seq(1, N.scenarios, N.erfn), 
                 seq(2, N.scenarios, N.erfn), 
                 seq(3, N.scenarios, N.erfn), 
                 seq(4, N.scenarios, N.erfn)) #all linear,  then s-shape,  etc.
  for (name.x in x.true.names) {
    write.table(coverage.25.propn.obscorr[tmp.order, , name.x], file = paste(
      "Coverage25_OBSCORR_", format(Sys.Date(), "%Y%m%d"), ".csv", sep = ""), 
      sep = ",", append = T, col.names = NA)
    write.table(coverage.25.propn.halfcorr[tmp.order, , name.x], file = paste(
      "Coverage25_HALFCORR_", format(Sys.Date(), "%Y%m%d"), ".csv", sep = ""), 
      sep = ",", append = T, col.names = NA)
  }
  rm(counter, name.x, tmp.bkmr, tmp.bart, tmp.star, tmp.coverage.25, tmp.order)
  
  ##----------------------------------------------------------------------------
  #table of variable selection results
  #dim: rows = N.scenarios x 4 methods; cols = 9 statistics x 2 corr structures
  N.stats<-9
  vs.table <- matrix(NaN, N.scenarios*N.methods, N.stats*N.corr) 
  colnames(vs.table) <- c("sens.halfcorr", "sens.obscorr", 
                          "spec.halfcorr", "spec.obscorr", 
                          "FDR.halfcorr", "FDR.obscorr", 
                          "FPR.halfcorr", "FPR.obscorr", 
                          "NPV.halfcorr", "NPV.obscorr", 
                          "F1.halfcorr", "F1.obscorr", 
                          "rank_correct_gr.halfcorr", "rank_correct_gr.obscorr", 
                          "prop_rank_correct_gr.halfcorr", 
                          "prop_rank_correct_gr.obscorr", 
                          "sum_sel.halfcorr", "sum_sel.obscorr")
  rownames(vs.table) <- as.vector(
    sapply(scenario.names, 
           FUN = function(x) paste(x, 
                                   c("bkmr", "bart", "star", "lasso"), 
                                   sep = ".")))
  counter = 0
  for (vs.stat in c("sens", "spec", "FDR", "FPR", "NPV", "F1", 
                    "rank_correct_gr", "prop_rank_correct_gr", "sum.sel")) {
    counter = counter + 1
    vs.table[, seq(1, N.stats*N.corr, N.corr)[counter]] <- as.vector(
      rbind(colMeans(res.halfcorr.bkmr[[vs.stat]], na.rm = T), 
            colMeans(res.halfcorr.bart[[vs.stat]], na.rm = T), 
            colMeans(res.halfcorr.star[[vs.stat]], na.rm = T), 
            colMeans(res.halfcorr.lasso[[vs.stat]], na.rm = T)))
    vs.table[, seq(2, N.stats*N.corr, N.corr)[counter]] <- as.vector(
      rbind(colMeans(res.obscorr.bkmr[[vs.stat]], na.rm = T), 
            colMeans(res.obscorr.bart[[vs.stat]], na.rm = T), 
            colMeans(res.obscorr.star[[vs.stat]], na.rm = T), 
            colMeans(res.obscorr.lasso[[vs.stat]], na.rm = T)))
  }
  #re-arrange rows:
  tmp.order <- c(seq(1, N.scenarios, N.erfn), 
                 seq(2, N.scenarios, N.erfn), 
                 seq(3, N.scenarios, N.erfn), 
                 seq(4, N.scenarios, N.erfn)) #all linear,  then s-shape,  etc.
  vs.table <- vs.table[as.vector(
    t(matrix(c(1:(N.scenarios*N.methods)), 
             N.scenarios, N.methods, byrow = T)[tmp.order, ])), ]
  write.table(vs.table, 
              file = paste("VarSelection_", format(Sys.Date(), "%Y%m%d"), 
                           ".csv", sep = ""), 
              sep = ",", col.names = NA)
  rm(vs.stat, counter, tmp.order)
  
  ##----------------------------------------------------------------------------
  ##VAR SELECTION RESULTS PLOTS

  message("Preparing Figures 3 to 10...")
  #variable selection results plots data.frame,  means only
  N.stats<-8
  vs.plots.df <- rbind(vs.table[, seq(1, N.stats*N.corr, N.corr)], 
                       vs.table[, seq(2, N.stats*N.corr, N.corr)])
  colnames(vs.plots.df) <- c("Sensitivity", "Specificity", "FDR", "FPR", "NPV", 
                             "F1", "rank_correct_gr", "prop_rank_correct_gr")
  vs.plots.df <- as.data.frame(vs.plots.df)
  vs.plots.df$mtype <- factor(
    rep(rep(c(rep("J  =  6,  low sparsity", N.snr*N.methods), 
              rep("J  =  12,  high sparsity", N.snr*N.methods)), 
            N.erfn), N.corr), 
    levels = c("J  =  6,  low sparsity", "J  =  12,  high sparsity"))
  vs.plots.df$SNR <- factor(
    rep(rep(c(rep("Low", N.methods), rep("High", N.methods)), 
            N.m*N.erfn), N.corr), levels = c("High", "Low"))
  vs.plots.df$erfn <- factor(
    rep(c(rep("Linear", N.m*N.snr*N.methods), 
          rep("S-shaped", N.m*N.snr*N.methods), 
          rep("Inverse-U-shaped\n(quadratic)", N.m*N.snr*N.methods), 
          rep("Inverse-U-shaped\n(asymmetric)", N.m*N.snr*N.methods)), N.corr), 
    levels = c("Linear", "S-shaped", "Inverse-U-shaped\n(quadratic)", 
               "Inverse-U-shaped\n(asymmetric)"))
  vs.plots.df$Method <- factor(
    rep(rep(c("BKMR", "BART", "BSTARSS", "LASSO"), N.m*N.snr*N.erfn), N.corr), 
    levels = c("BKMR", "BART", "BSTARSS", "LASSO"))
  vs.plots.df$corr <- factor(c(rep("Half", N.m*N.snr*N.erfn*N.methods), 
                               rep("Observed", N.m*N.snr*N.erfn*N.methods)), 
                             levels = c("Observed", "Half"))
  vs.plots.df$Precision <- 1-vs.plots.df$FDR
  
  #variable selection results plots data.frame,  all results
  N.stats<-8
  vs.plots.df.long <- matrix(NaN, N.scenarios*N.methods*reps*N.corr, N.stats) 
  colnames(vs.plots.df.long) <- c("Sensitivity", "Specificity", 
                                  "FDR", "FPR", "NPV", "F1", 
                                  "rank_correct_gr", "prop_rank_correct_gr")
  counter = 0
  for (vs.stat in c("sens", "spec", "FDR", "FPR", "NPV", "F1", 
                    "rank_correct_gr", "prop_rank_correct_gr")) {
    counter = counter + 1
    vs.plots.df.long[, counter] <- as.vector(
      rbind(res.halfcorr.bkmr[[vs.stat]], res.halfcorr.bart[[vs.stat]], 
            res.halfcorr.star[[vs.stat]], res.halfcorr.lasso[[vs.stat]], 
            res.obscorr.bkmr[[vs.stat]], res.obscorr.bart[[vs.stat]], 
            res.obscorr.star[[vs.stat]], res.obscorr.lasso[[vs.stat]]))
  }
  vs.plots.df.long <- as.data.frame(vs.plots.df.long)
  vs.plots.df.long$Scenario <- rep(
    as.vector(sapply(
      rep(scenario.names, each = N.corr), 
      FUN = function(x) paste(x, c("bkmr", "bart", "star", "lasso"), 
                              sep = "."))), each = reps)
  vs.plots.df.long$corr <- factor(
    rep(c(rep("Half", reps*N.methods), 
          rep("Observed", reps*N.methods)), N.scenarios), 
    levels = c("Observed", "Half"))
  vs.plots.df.long$Method <- factor(
    rep(c(rep(c(rep("BKMR", reps), 
                rep("BART", reps), 
                rep("BSTARSS", reps), 
                rep("LASSO", reps)), N.corr)), N.scenarios), 
    levels = c("BKMR", "BART", "BSTARSS", "LASSO"))
  vs.plots.df.long$mtype <- ""
  vs.plots.df.long$mtype[grep("^m1", vs.plots.df.long$Scenario)] <- 
    "J  =  6,  low sparsity"
  vs.plots.df.long$mtype[grep("^m2", vs.plots.df.long$Scenario)] <- 
    "J  =  12,  high sparsity"
  vs.plots.df.long$mtype <- factor(
    vs.plots.df.long$mtype, 
    levels = c("J  =  6,  low sparsity", "J  =  12,  high sparsity"))
  vs.plots.df.long$SNR <- ""
  vs.plots.df.long$SNR[grep("*loSNR", vs.plots.df.long$Scenario)] <- "Low"
  vs.plots.df.long$SNR[grep("*hiSNR", vs.plots.df.long$Scenario)] <- "High"
  vs.plots.df.long$SNR <- factor(vs.plots.df.long$SNR, 
                                 levels = c("High", "Low"))
  vs.plots.df.long$erfn <- ""
  vs.plots.df.long$erfn[grep("*SNR.l", vs.plots.df.long$Scenario)] <- 
    "Linear"
  vs.plots.df.long$erfn[grep("*SNR.s", vs.plots.df.long$Scenario)] <- 
    "S-shaped"
  vs.plots.df.long$erfn[grep("*SNR.q", vs.plots.df.long$Scenario)] <- 
    "Inverse-U-shaped\n(quadratic)"
  vs.plots.df.long$erfn[grep("*SNR.a", vs.plots.df.long$Scenario)] <- 
    "Inverse-U-shaped\n(asymmetric)"
  vs.plots.df.long$erfn <- factor(
    vs.plots.df.long$erfn, 
    levels = c("Linear", "S-shaped", "Inverse-U-shaped\n(quadratic)", 
               "Inverse-U-shaped\n(asymmetric)"))
  
  ##----------------------------------------------------------------------------
  #FIGURE 3: plot of y-axis = sens vs. x-axis = spec
  vs.plot.ySens_xSpec <- ggplot(
    vs.plots.df, aes(Specificity, Sensitivity, color = factor(Method))) + 
    facet_grid(mtype ~ erfn) + 
    geom_point(aes(shape = SNR), size = 4, stroke = 1) + 
    geom_line(aes(linetype = corr), size = 1.1) + 
    scale_colour_discrete(name = "Method") + 
    scale_shape_manual(name = "Signal-to-noise \nratio",  values = c(16, 2)) + 
    scale_linetype_discrete(name = "Correlation \nstructure") + 
    scale_x_continuous(breaks = c(0.25, 0.5, 0.75, 1), limits = c(0.2, 1)) + 
    theme_bw() + 
    theme(text = element_text(size = 16), axis.text = element_text(size = 16), 
          strip.text.x = element_text(size = 16), 
          strip.text.y = element_text(size = 16))
  print(vs.plot.ySens_xSpec)
  ggsave("Figure_3.pdf", plot = vs.plot.ySens_xSpec, 
         height = 752/72, width = 1335/72)
  
  #FIGURE 4: plot of y-axis = prec vs. x-axis = NPV
  vs.plot.yPrec_xNPV <- ggplot(
    vs.plots.df, aes(NPV, Precision, color = factor(Method))) + 
    facet_grid(mtype ~ erfn) + 
    geom_point(aes(shape = SNR), size = 4, stroke = 1) + 
    geom_line(aes(linetype = corr), size = 1.1) + 
    scale_colour_discrete(name = "Method") + 
    scale_shape_manual(name = "Signal-to-noise \nratio",  values = c(16, 2)) + 
    scale_linetype_discrete(name = "Correlation \nstructure") + 
    scale_y_continuous(breaks = c(0.25, 0.5, 0.75, 1)) + 
    scale_x_continuous(breaks = c(0.25, 0.5, 0.75, 1), limits = c(0.2, 1)) + 
    xlab("Negative predictive value") + 
    theme_bw() + 
    theme(text = element_text(size = 16), axis.text = element_text(size = 16), 
          strip.text.x = element_text(size = 16), 
          strip.text.y = element_text(size = 16))
  print(vs.plot.yPrec_xNPV)
  ggsave("Figure_4.pdf", plot = vs.plot.yPrec_xNPV, 
         height = 752/72, width = 1335/72)
  
  #FIGURE 5: boxplots of F1
  F1.boxplot <- ggplot(
    vs.plots.df.long, aes(interaction(corr, SNR), F1, fill = Method)) + 
    facet_grid(mtype ~ erfn) + 
    geom_boxplot(aes(color = Method), 
                 position = position_dodge(width = 0.7), 
                 size = 0.5, width = 0.6) + 
    stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2, 
                 position = position_dodge(width = 0.7)) + 
    stat_summary(fun.y = median, geom = "point", shape = 95, size = 3, 
                 stroke = 3, position = position_dodge(width = 0.7)) + 
    ylab("F1-statistic") + 
    xlab("Correlation structure (observed,  half) and 
         signal-to-noise ratio (low,  high) categories") + 
    theme_bw() + 
    theme(text = element_text(size = 16), axis.text = element_text(size = 16), 
          strip.text.x = element_text(size = 16), 
          strip.text.y = element_text(size = 16), 
          axis.text.x = element_text(size = 13)) + 
    scale_x_discrete(labels = c("Obs., \nHigh\n", "Half, \nHigh\n", 
                                "Obs., \nLow\n", "Half, \nLow\n"))
  print(F1.boxplot)
  ggsave("Figure_5.pdf", plot = F1.boxplot, height = 752/72, width = 1335/72)
  
  #FIGURE 6: plot of rank_correct_gr
  vs.plot.rank_correct_gr <- ggplot(
    subset(vs.plots.df, ((Method %in% c("BKMR", "BART", "BSTARSS")))), 
    aes(interaction(corr, SNR), rank_correct_gr, color = factor(Method))) + 
    facet_grid(mtype ~ erfn) + 
    geom_point(size = 4, stroke = 1, position = position_dodge(width = .1)) + 
    xlab("Correlation structure (observed,  half) and 
         signal-to-noise ratio (low,  high) categories") + 
    ylab("Proportion of replications with all outcome-associated\nexposures 
         ranked above unassociated exposures") + 
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 0.9)) + 
    scale_color_manual(name = "Method", 
                       values = c("#F8766D", "#7CAE00", "#00BFC4")) + 
    theme_bw() + 
    theme(text = element_text(size = 16), axis.text = element_text(size = 16), 
          strip.text.x = element_text(size = 16), 
          strip.text.y = element_text(size = 16), 
          axis.text.x = element_text(size = 13)) + 
    scale_x_discrete(labels = c("Obs., \nHigh\n", "Half, \nHigh\n", 
                                "Obs., \nLow\n", "Half, \nLow\n"))
  print(vs.plot.rank_correct_gr)
  ggsave("Figure_6.pdf", plot = vs.plot.rank_correct_gr, 
         height = 752/72, width = 1335/72)
  
  #FIGURE 7: plot of prop_rank_correct_gr with 95% confidence intervals
  vs.plots.df$prop_rank_correct_gr_ci <- qnorm(0.975)*
    sqrt(vs.plots.df$prop_rank_correct_gr*
           (1-vs.plots.df$prop_rank_correct_gr)/(K.true*reps))
  vs.plot.prop_rank_correct_gr_ci <- ggplot(
    subset(vs.plots.df, ((Method %in% c("BKMR", "BART", "BSTARSS")))), 
    aes(interaction(corr, SNR), prop_rank_correct_gr, color = factor(Method))) + 
    facet_grid(mtype ~ erfn) + 
    geom_point(size = 4, stroke = 1, position = position_dodge(width = .3)) + 
    geom_errorbar(aes(ymin = prop_rank_correct_gr-prop_rank_correct_gr_ci, 
                      ymax = prop_rank_correct_gr + prop_rank_correct_gr_ci), 
                  width = .4, position = position_dodge(width = .3)) + 
    xlab("Correlation structure (observed,  half) and 
         signal-to-noise ratio (low,  high) categories") + 
    ylab("Mean proportion of outcome-associated\nexposures 
         ranked above unassociated exposures") + 
    scale_y_continuous(breaks = c(0.25, 0.5, 0.75, 1), limits = c(0.35, 1)) + 
    scale_color_manual(name = "Method", 
                       values = c("#F8766D", "#7CAE00", "#00BFC4")) + 
    theme_bw() + 
    theme(text = element_text(size = 16), axis.text = element_text(size = 16), 
          strip.text.x = element_text(size = 16), 
          strip.text.y = element_text(size = 16), 
          axis.text.x = element_text(size = 13)) + 
    scale_x_discrete(labels = c("Obs., \nHigh\n", "Half, \nHigh\n", 
                                "Obs., \nLow\n", "Half, \nLow\n"))
  print(vs.plot.prop_rank_correct_gr_ci)
  ggsave("Figure_7.pdf", plot = vs.plot.prop_rank_correct_gr_ci, 
         height = 752/72, width = 1335/72)
  
  ##----------------------------------------------------------------------------
  ##MSE AND COVERAGE PLOTS 
  
  #MSE/oracle_MSE dataframe for plotting
  mse_propgam.25.obscorr <- mse.25.obscorr
  for (name.x in x.true.names) {
    mse_propgam.25.obscorr[, c(1:4), name.x] <- 
      mse_propgam.25.obscorr[, c(1:4), name.x]/
      mse_propgam.25.obscorr[, "GAM.25%", name.x]
    mse_propgam.25.obscorr[, c(5:8), name.x] <- 
      mse_propgam.25.obscorr[, c(5:8), name.x]/
      mse_propgam.25.obscorr[, "GAM.50%", name.x]
    mse_propgam.25.obscorr[, c(9:12), name.x] <- 
      mse_propgam.25.obscorr[, c(9:12), name.x]/
      mse_propgam.25.obscorr[, "GAM.75%", name.x]
  }
  mse_propgam.25.obscorr <- mse_propgam.25.obscorr[, -c(4, 8, 12), ] #remove GAM
  mse_propgam.25.halfcorr <- mse.25.halfcorr
  for (name.x in x.true.names) {
    mse_propgam.25.halfcorr[, c(1:4), name.x] <- 
      mse_propgam.25.halfcorr[, c(1:4), name.x]/
      mse_propgam.25.halfcorr[, "GAM.25%", name.x]
    mse_propgam.25.halfcorr[, c(5:8), name.x] <- 
      mse_propgam.25.halfcorr[, c(5:8), name.x]/
      mse_propgam.25.halfcorr[, "GAM.50%", name.x]
    mse_propgam.25.halfcorr[, c(9:12), name.x] <- 
      mse_propgam.25.halfcorr[, c(9:12), name.x]/
      mse_propgam.25.halfcorr[, "GAM.75%", name.x]
  }
  mse_propgam.25.halfcorr <- mse_propgam.25.halfcorr[, -c(4, 8, 12), ]
  #dataframe
  mse_propgam.plots.df.long <- rbind(
    as.data.frame.table(mse_propgam.25.halfcorr), 
    as.data.frame.table(mse_propgam.25.obscorr))
  colnames(mse_propgam.plots.df.long) <- 
    c("Scenario", "Method.Ptile", "Exposure", "MSEpropGAM")
  mse_propgam.plots.df.long$Exposure <- factor(
    mse_propgam.plots.df.long$Exposure, levels = x.true.names)
  mse_propgam.plots.df.long$corr <- factor(
    c(rep("Half", N.ptiles*(N.methods-1)*K.true*N.scenarios), 
      rep("Observed", N.ptiles*(N.methods-1)*K.true*N.scenarios)), 
    levels = c("Observed", "Half"))
  mse_propgam.plots.df.long$Method <- ""
  mse_propgam.plots.df.long$Method[
    grep("^BKMR", mse_propgam.plots.df.long$Method.Ptile)] <- "BKMR"
  mse_propgam.plots.df.long$Method[
    grep("^BART", mse_propgam.plots.df.long$Method.Ptile)] <- "BART"
  mse_propgam.plots.df.long$Method[
    grep("^STAR", mse_propgam.plots.df.long$Method.Ptile)] <- "BSTARSS"
  mse_propgam.plots.df.long$Method <- factor(
    mse_propgam.plots.df.long$Method, 
    levels = c("BKMR", "BART", "BSTARSS", "GAM"))
  mse_propgam.plots.df.long$Percentile <- ""
  mse_propgam.plots.df.long$Percentile[
    grep("*25", mse_propgam.plots.df.long$Method.Ptile)] <- "25%"
  mse_propgam.plots.df.long$Percentile[
    grep("*50", mse_propgam.plots.df.long$Method.Ptile)] <- "50%"
  mse_propgam.plots.df.long$Percentile[
    grep("*75", mse_propgam.plots.df.long$Method.Ptile)] <- "75%"
  mse_propgam.plots.df.long$Percentile <- factor(
    mse_propgam.plots.df.long$Percentile, levels = c("25%", "50%", "75%"))
  mse_propgam.plots.df.long$mtype <- ""
  mse_propgam.plots.df.long$mtype[
    grep("^m1", mse_propgam.plots.df.long$Scenario)] <- "J  =  6"
  mse_propgam.plots.df.long$mtype[
    grep("^m2", mse_propgam.plots.df.long$Scenario)] <- "J  =  12"
  mse_propgam.plots.df.long$mtype <- factor(
    mse_propgam.plots.df.long$mtype, levels = c("J  =  6", "J  =  12"))
  mse_propgam.plots.df.long$SNR <- ""
  mse_propgam.plots.df.long$SNR[
    grep("*loSNR", mse_propgam.plots.df.long$Scenario)] <- "Low"
  mse_propgam.plots.df.long$SNR[
    grep("*hiSNR", mse_propgam.plots.df.long$Scenario)] <- "High"
  mse_propgam.plots.df.long$SNR <- factor(
    mse_propgam.plots.df.long$SNR, levels = c("High", "Low"))
  mse_propgam.plots.df.long$erfn <- ""
  mse_propgam.plots.df.long$erfn[
    grep("*SNR.l", mse_propgam.plots.df.long$Scenario)] <- "Linear"
  mse_propgam.plots.df.long$erfn[
    grep("*SNR.s", mse_propgam.plots.df.long$Scenario)] <- "S-shaped"
  mse_propgam.plots.df.long$erfn[
    grep("*SNR.q", mse_propgam.plots.df.long$Scenario)] <- 
    "Inverse-U-shaped\n(quadratic)"
  mse_propgam.plots.df.long$erfn[
    grep("*SNR.a", mse_propgam.plots.df.long$Scenario)] <- 
    "Inverse-U-shaped\n(asymmetric)"
  mse_propgam.plots.df.long$erfn <- factor(
    mse_propgam.plots.df.long$erfn, 
    levels = c("Linear", "S-shaped", "Inverse-U-shaped\n(quadratic)", 
               "Inverse-U-shaped\n(asymmetric)"))
  
  ##----------------------------------------------------------------------------
  #FIGURE 8: all percentiles on same plot,  mean of all exposures
  mse_propgam.plot <- ggplot(
    subset(mse_propgam.plots.df.long, 
           ((Exposure %in% c("MPB", "BP3", "PPB", "BPA")))), 
    aes(interaction(corr, SNR), MSEpropGAM, color = factor(Method))) + 
    facet_grid(mtype ~ erfn) + 
    stat_summary(
      aes(shape = "25%"), 
      data = subset(mse_propgam.plots.df.long, 
                    ((Exposure %in% c("MPB", "BP3", "PPB", "BPA"))&
                       (Percentile %in% c("25%")))), 
      fun.y = mean, geom = "point", size = 4, 
      position = position_dodge(width = 0.5)) + 
    stat_summary(
      aes(shape = "50%"), 
      data = subset(mse_propgam.plots.df.long, 
                    ((Exposure %in% c("MPB", "BP3", "PPB", "BPA"))&
                       (Percentile %in% c("50%")))), 
      fun.y = mean, geom = "point", size = 4, 
      position = position_dodge(width = 0.5)) + 
    stat_summary(
      aes(shape = "75%"), 
      data = subset(mse_propgam.plots.df.long, 
                    ((Exposure %in% c("MPB", "BP3", "PPB", "BPA"))&
                       (Percentile %in% c("75%")))), 
      fun.y = mean, geom = "point", size = 4, 
      position = position_dodge(width = 0.5)) + 
    xlab("Correlation structure (observed,  half) and 
         signal-to-noise ratio (low,  high) categories") + 
    ylab("Ratio of method to oracle mean-squared error") + 
    geom_hline(yintercept = 1, linetype = 2) + 
    theme_bw() + 
    theme(text = element_text(size = 16), axis.text = element_text(size = 16), 
          strip.text.x = element_text(size = 16), 
          strip.text.y = element_text(size = 16), 
          axis.text.x = element_text(size = 13)) + 
    scale_color_manual(name = "Method", 
                       values = c("#F8766D", "#7CAE00", "#00BFC4")) + 
    scale_x_discrete(labels = c("Obs., \nHigh\n", "Half, \nHigh\n", 
                                "Obs., \nLow\n", "Half, \nLow\n")) + 
    scale_shape_manual(name = "Percentile", c("25%", "50%", "75%"), 
                       values = c(2, 19, 3))
  print(mse_propgam.plot)
  ggsave("Figure_8.pdf", plot = mse_propgam.plot, 
         height = 752/72, width = 1335/72)
  
  ##----------------------------------------------------------------------------
  #coverage dataframe for plotting
  coverage.plots.df.long <- rbind(
    as.data.frame.table(coverage.25.propn.halfcorr), 
    as.data.frame.table(coverage.25.propn.obscorr))
  colnames(coverage.plots.df.long) <- 
    c("Scenario", "Method.Ptile", "Exposure", "Coverage")
  coverage.plots.df.long <- 
    coverage.plots.df.long[-grep("^GAM", coverage.plots.df.long$Method.Ptile), ]  #remove GAM entries
  coverage.plots.df.long$Exposure <- factor(
    coverage.plots.df.long$Exposure, levels = x.true.names)
  coverage.plots.df.long$corr <- factor(
    c(rep("Half", N.ptiles*(N.methods-1)*K.true*N.scenarios), 
      rep("Observed", N.ptiles*(N.methods-1)*K.true*N.scenarios)), 
    levels = c("Observed", "Half"))
  coverage.plots.df.long$Method <- ""
  coverage.plots.df.long$Method[
    grep("^BKMR", coverage.plots.df.long$Method.Ptile)] <- "BKMR"
  coverage.plots.df.long$Method[
    grep("^BART", coverage.plots.df.long$Method.Ptile)] <- "BART"
  coverage.plots.df.long$Method[
    grep("^STAR", coverage.plots.df.long$Method.Ptile)] <- "BSTARSS"
  coverage.plots.df.long$Method <- factor(
    coverage.plots.df.long$Method, levels = c("BKMR", "BART", "BSTARSS", "GAM"))
  coverage.plots.df.long$Percentile <- ""
  coverage.plots.df.long$Percentile[
    grep("*25", coverage.plots.df.long$Method.Ptile)] <- "25%"
  coverage.plots.df.long$Percentile[
    grep("*50", coverage.plots.df.long$Method.Ptile)] <- "50%"
  coverage.plots.df.long$Percentile[
    grep("*75", coverage.plots.df.long$Method.Ptile)] <- "75%"
  coverage.plots.df.long$Percentile <- factor(
    coverage.plots.df.long$Percentile, levels = c("25%", "50%", "75%"))
  coverage.plots.df.long$mtype <- ""
  coverage.plots.df.long$mtype[
    grep("^m1", coverage.plots.df.long$Scenario)] <- "J  =  6"
  coverage.plots.df.long$mtype[
    grep("^m2", coverage.plots.df.long$Scenario)] <- "J  =  12"
  coverage.plots.df.long$mtype <- factor(
    coverage.plots.df.long$mtype, levels = c("J  =  6", "J  =  12"))
  coverage.plots.df.long$SNR <- ""
  coverage.plots.df.long$SNR[
    grep("*loSNR", coverage.plots.df.long$Scenario)] <- "Low"
  coverage.plots.df.long$SNR[
    grep("*hiSNR", coverage.plots.df.long$Scenario)] <- "High"
  coverage.plots.df.long$SNR <- factor(
    coverage.plots.df.long$SNR, levels = c("High", "Low"))
  coverage.plots.df.long$erfn <- ""
  coverage.plots.df.long$erfn[
    grep("*SNR.l", coverage.plots.df.long$Scenario)] <- "Linear"
  coverage.plots.df.long$erfn[
    grep("*SNR.s", coverage.plots.df.long$Scenario)] <- "S-shaped"
  coverage.plots.df.long$erfn[
    grep("*SNR.q", coverage.plots.df.long$Scenario)] <- 
    "Inverse-U-shaped\n(quadratic)"
  coverage.plots.df.long$erfn[
    grep("*SNR.a", coverage.plots.df.long$Scenario)] <- 
    "Inverse-U-shaped\n(asymmetric)"
  coverage.plots.df.long$erfn <- factor(
    coverage.plots.df.long$erfn, 
    levels = c("Linear", "S-shaped", "Inverse-U-shaped\n(quadratic)", 
               "Inverse-U-shaped\n(asymmetric)"))
  
  ##----------------------------------------------------------------------------
  #FIGURE 9: all percentiles on same plot,  mean of all exposures
  coverage.plot <- ggplot(
    subset(coverage.plots.df.long, 
           ((Exposure %in% c("MPB", "BP3", "PPB", "BPA")))), 
    aes(interaction(corr, SNR), Coverage, color = factor(Method))) + 
    facet_grid(mtype ~ erfn) + 
    stat_summary(
      aes(shape = "25%"), 
      data = subset(coverage.plots.df.long, 
                    ((Exposure %in% c("MPB", "BP3", "PPB", "BPA"))&
                       (Percentile %in% c("25%")))), 
      fun.y = mean, geom = "point", size = 4, 
      position = position_dodge(width = 0.5)) + 
    stat_summary(
      aes(shape = "50%"), 
      data = subset(coverage.plots.df.long, 
                    ((Exposure %in% c("MPB", "BP3", "PPB", "BPA"))&
                       (Percentile %in% c("50%")))), 
      fun.y = mean, geom = "point", size = 4, 
      position = position_dodge(width = 0.5)) + 
    stat_summary(
      aes(shape = "75%"), 
      data = subset(coverage.plots.df.long, 
                    ((Exposure %in% c("MPB", "BP3", "PPB", "BPA"))&
                       (Percentile %in% c("75%")))), 
      fun.y = mean, geom = "point", size = 4, 
      position = position_dodge(width = 0.5)) + 
    xlab("Correlation structure (observed,  half) and 
         signal-to-noise ratio (low,  high) categories") + 
    ylab("Coverage") + 
    geom_hline(yintercept = 0.9, linetype = 2) + 
    theme_bw() + 
    theme(text = element_text(size = 16), axis.text = element_text(size = 16), 
          strip.text.x = element_text(size = 16), 
          strip.text.y = element_text(size = 16), 
          axis.text.x = element_text(size = 13)) + 
    scale_color_manual(name = "Method", 
                       values = c("#F8766D", "#7CAE00", "#00BFC4")) + 
    scale_x_discrete(labels = c("Obs., \nHigh\n", "Half, \nHigh\n", 
                                "Obs., \nLow\n", "Half, \nLow\n")) + 
    scale_shape_manual(name = "Percentile", c("25%", "50%", "75%"), 
                       values = c(2, 19, 3))
  print(coverage.plot)
  ggsave("Figure_9.pdf", plot = coverage.plot, height = 752/72, width = 1335/72)
  
  ##----------------------------------------------------------------------------
  ##PLOTS OF EXPOSURE-RESPONSE CURVES
  
  #exp-resp curve dataframes for plotting
  N.stats<-11
  true.erc.plots.df <- rbind(
    as.data.frame.table(true.erc.10), 
    as.data.frame.table(true.erc.10.halfcorr))
  colnames(true.erc.plots.df) <- 
    c("RepID", "Ptile", "Exposure", "erfn", "Estimate")
  true.erc.plots.df$Percentile <- NaN
  true.erc.plots.df$Percentile[grep("*0", true.erc.plots.df$Ptile)] <- 0
  true.erc.plots.df$Percentile[grep("*10", true.erc.plots.df$Ptile)] <- 10
  true.erc.plots.df$Percentile[grep("*20", true.erc.plots.df$Ptile)] <- 20
  true.erc.plots.df$Percentile[grep("*30", true.erc.plots.df$Ptile)] <- 30
  true.erc.plots.df$Percentile[grep("*40", true.erc.plots.df$Ptile)] <- 40
  true.erc.plots.df$Percentile[grep("*50", true.erc.plots.df$Ptile)] <- 50
  true.erc.plots.df$Percentile[grep("*60", true.erc.plots.df$Ptile)] <- 60
  true.erc.plots.df$Percentile[grep("*70", true.erc.plots.df$Ptile)] <- 70
  true.erc.plots.df$Percentile[grep("*80", true.erc.plots.df$Ptile)] <- 80
  true.erc.plots.df$Percentile[grep("*90", true.erc.plots.df$Ptile)] <- 90
  true.erc.plots.df$Percentile[grep("*100", true.erc.plots.df$Ptile)] <- 100
  true.erc.plots.df$corr <- factor(
    c(rep("Observed", N.stats*reps*N.scenarios), 
      rep("Half", N.stats*reps*N.scenarios)),
    levels = c("Observed", "Half"))
  levels(true.erc.plots.df$erfn) <- c(
    "Linear", "S-shaped", "Inverse-U-shaped\n(quadratic)", 
    "Inverse-U-shaped\n(asymmetric)")
  true.erc.plots.df <- rbind(true.erc.plots.df, true.erc.plots.df, 
                             true.erc.plots.df, true.erc.plots.df)
  true.erc.plots.df$mtype <- factor(
    c(rep("J  =  6,  low", N.stats*reps*N.corr*N.snr*N.erfn*N.methods), 
      rep("J  =  12,  high", N.stats*reps*N.corr*N.snr*N.erfn*N.methods)), 
    levels = c("J  =  6,  low", "J  =  12,  high"))
  true.erc.plots.df$SNR <- factor(
    c(rep("High", N.stats*reps*N.m*N.erfn*N.methods), 
      rep("Low", N.stats*reps*N.m*N.erfn*N.methods), 
      rep("High", N.stats*reps*N.m*N.erfn*N.methods), 
      rep("Low", N.stats*reps*N.m*N.erfn*N.methods)), 
    levels = c("High", "Low"))
  #add x-values at percentiles
  tmp.df <- rbind(
    do.call("rbind", replicate(
      N.scenarios, as.data.frame.table(x.true.stats.obscorr[, c(2:12), ]), 
      simplify = FALSE)), 
    do.call("rbind", replicate(
      N.scenarios, as.data.frame.table(x.true.stats.halfcorr[, c(2:12), ]), 
      simplify = FALSE)))
  colnames(tmp.df) <- c("RepID", "Ptile", "Exposure", "X")
  true.erc.plots.df$X <- tmp.df$X
  rm(tmp.df)
  #estimated erc data.frame
  for (name.x in x.true.names) {
    tmp.df <- rbind(
      as.data.frame.table(res.halfcorr.bkmr$erc.10[, , name.x, , "point"]), 
      as.data.frame.table(res.obscorr.bkmr$erc.10[, , name.x, , "point"]), 
      as.data.frame.table(res.halfcorr.bart$erc.10[, , name.x, , "point"]), 
      as.data.frame.table(res.obscorr.bart$erc.10[, , name.x, , "point"]), 
      as.data.frame.table(res.halfcorr.star$erc.10[, , name.x, , "point"]), 
      as.data.frame.table(res.obscorr.star$erc.10[, , name.x, , "point"]), 
      as.data.frame.table(res.halfcorr.gam$erc.10[, , name.x, , "point"]), 
      as.data.frame.table(res.obscorr.gam$erc.10[, , name.x, , "point"]))
    if (name.x == x.true.names[1]) {
      curve.plots.df <- tmp.df
      colnames(curve.plots.df) <- c(
        "RepID", "Ptile", "Scenario", paste0("Y.", name.x))
    } else {
      curve.plots.df[[paste0("Y.", name.x)]] <- tmp.df$Freq
    }
  }
  for (name.x in x.true.names) {
    tmp.df.loCI <- rbind(
      as.data.frame.table(res.halfcorr.bkmr$erc.10[, , name.x, , "loCI"]), 
      as.data.frame.table(res.obscorr.bkmr$erc.10[, , name.x, , "loCI"]), 
      as.data.frame.table(res.halfcorr.bart$erc.10[, , name.x, , "loCI"]), 
      as.data.frame.table(res.obscorr.bart$erc.10[, , name.x, , "loCI"]), 
      as.data.frame.table(res.halfcorr.star$erc.10[, , name.x, , "loCI"]), 
      as.data.frame.table(res.obscorr.star$erc.10[, , name.x, , "loCI"]), 
      as.data.frame.table(res.halfcorr.gam$erc.10[, , name.x, , "loCI"]), 
      as.data.frame.table(res.obscorr.gam$erc.10[, , name.x, , "loCI"]))
    tmp.df.hiCI <- rbind(
      as.data.frame.table(res.halfcorr.bkmr$erc.10[, , name.x, , "hiCI"]), 
      as.data.frame.table(res.obscorr.bkmr$erc.10[, , name.x, , "hiCI"]), 
      as.data.frame.table(res.halfcorr.bart$erc.10[, , name.x, , "hiCI"]), 
      as.data.frame.table(res.obscorr.bart$erc.10[, , name.x, , "hiCI"]), 
      as.data.frame.table(res.halfcorr.star$erc.10[, , name.x, , "hiCI"]), 
      as.data.frame.table(res.obscorr.star$erc.10[, , name.x, , "hiCI"]), 
      as.data.frame.table(res.halfcorr.gam$erc.10[, , name.x, , "hiCI"]), 
      as.data.frame.table(res.obscorr.gam$erc.10[, , name.x, , "hiCI"]))
    curve.plots.df[[paste0("Y.loCI.", name.x)]] <- tmp.df.loCI$Freq
    curve.plots.df[[paste0("Y.hiCI.", name.x)]] <- tmp.df.hiCI$Freq
  }
  rm(tmp.df, tmp.df.loCI, tmp.df.hiCI)
  curve.plots.df$corr <- factor(
    rep(c(rep("Half", N.stats*reps*N.scenarios), 
          rep("Observed", N.stats*reps*N.scenarios)), N.methods), 
    levels = c("Observed", "Half"))
  curve.plots.df$Method <- factor(
    c(rep("BKMR", N.stats*reps*N.scenarios*N.corr), 
      rep("BART", N.stats*reps*N.scenarios*N.corr), 
      rep("BSTARSS", N.stats*reps*N.scenarios*N.corr), 
      rep("GAM", N.stats*reps*N.scenarios*N.corr)), 
    levels = c("BKMR", "BART", "BSTARSS", "GAM"))
  curve.plots.df$Percentile <- NaN
  curve.plots.df$Percentile[grep("*min", curve.plots.df$Ptile)] <- 0
  curve.plots.df$Percentile[grep("*10", curve.plots.df$Ptile)] <- 10
  curve.plots.df$Percentile[grep("*20", curve.plots.df$Ptile)] <- 20
  curve.plots.df$Percentile[grep("*30", curve.plots.df$Ptile)] <- 30
  curve.plots.df$Percentile[grep("*40", curve.plots.df$Ptile)] <- 40
  curve.plots.df$Percentile[grep("*50", curve.plots.df$Ptile)] <- 50
  curve.plots.df$Percentile[grep("*60", curve.plots.df$Ptile)] <- 60
  curve.plots.df$Percentile[grep("*70", curve.plots.df$Ptile)] <- 70
  curve.plots.df$Percentile[grep("*80", curve.plots.df$Ptile)] <- 80
  curve.plots.df$Percentile[grep("*90", curve.plots.df$Ptile)] <- 90
  curve.plots.df$Percentile[grep("*max", curve.plots.df$Ptile)] <- 100
  curve.plots.df$mtype <- ""
  curve.plots.df$mtype[
    grep("^m1", curve.plots.df$Scenario)] <- "J  =  6,  low"
  curve.plots.df$mtype[
    grep("^m2", curve.plots.df$Scenario)] <- "J  =  12,  high"
  curve.plots.df$mtype <- factor(
    curve.plots.df$mtype, levels = c("J  =  6,  low", "J  =  12,  high"))
  curve.plots.df$SNR <- ""
  curve.plots.df$SNR[grep("*loSNR", curve.plots.df$Scenario)] <- "Low"
  curve.plots.df$SNR[grep("*hiSNR", curve.plots.df$Scenario)] <- "High"
  curve.plots.df$SNR <- factor(curve.plots.df$SNR, levels = c("High", "Low"))
  curve.plots.df$erfn <- ""
  curve.plots.df$erfn[grep("*SNR.l", curve.plots.df$Scenario)] <- "Linear"
  curve.plots.df$erfn[grep("*SNR.s", curve.plots.df$Scenario)] <- "S-shaped"
  curve.plots.df$erfn[grep("*SNR.q", curve.plots.df$Scenario)] <- 
    "Inverse-U-shaped\n(quadratic)"
  curve.plots.df$erfn[grep("*SNR.a", curve.plots.df$Scenario)] <- 
    "Inverse-U-shaped\n(asymmetric)"
  curve.plots.df$erfn <- factor(
    curve.plots.df$erfn, 
    levels = c("Linear", "S-shaped", "Inverse-U-shaped\n(quadratic)", 
               "Inverse-U-shaped\n(asymmetric)"))
  #add x-values at percentiles
  for (name.x in x.true.names) {
    if (name.x == x.true.names[1]) {
      tmp.df <- rbind(
        do.call("rbind", replicate(
          N.scenarios, 
          as.data.frame.table(x.true.stats.halfcorr[, c(2:12), name.x]), 
          simplify = FALSE)), 
        do.call("rbind", replicate(
          N.scenarios, 
          as.data.frame.table(x.true.stats.obscorr[, c(2:12), name.x]), 
          simplify = FALSE)))
      colnames(tmp.df) <- c("RepID", "Ptile", paste0("X.", name.x))
    } else {
      tmp.df[[paste0("X.", name.x)]] <- rbind(
        do.call("rbind", replicate(
          N.scenarios, 
          as.data.frame.table(x.true.stats.halfcorr[, c(2:12), name.x]), 
          simplify = FALSE)), 
        do.call("rbind", replicate(
          N.scenarios, 
          as.data.frame.table(x.true.stats.obscorr[, c(2:12), name.x]), 
          simplify = FALSE)))$Freq
    }
  }
  for (name.x in x.true.names) curve.plots.df[[paste0("X.", name.x)]] <- 
    tmp.df[[paste0("X.", name.x)]]
  #combine true and estimated data.frames,  adding true erc values as an 
  #additional "method" for plotting purposes
  for (name.x in x.true.names) {
    if (name.x == x.true.names[1]) {
      true.erc.plots.df.wide <- 
        true.erc.plots.df[true.erc.plots.df$Exposure == name.x, ]
      colnames(true.erc.plots.df.wide)[5] <- c(paste0("Y.", name.x))
      colnames(true.erc.plots.df.wide)[10] <- c(paste0("X.", name.x))
    } else {
      true.erc.plots.df.wide[[paste0("Y.", name.x)]] <- 
        true.erc.plots.df[true.erc.plots.df$Exposure == name.x, ]$Estimate
      true.erc.plots.df.wide[[paste0("X.", name.x)]] <- 
        true.erc.plots.df[true.erc.plots.df$Exposure == name.x, ]$X
    }
    true.erc.plots.df.wide[[paste0("Y.loCI.", name.x)]] <- NaN
    true.erc.plots.df.wide[[paste0("Y.hiCI.", name.x)]] <- NaN
  }
  true.erc.plots.df.wide <- true.erc.plots.df.wide[, -3]
  true.erc.plots.df.wide$Method <- factor("True", levels = "True")
  true.erc.plots.df.wide$Scenario <- NA
  erc.plots.df <- rbind(curve.plots.df, true.erc.plots.df.wide)
  rm(tmp.df)
  
  ##----------------------------------------------------------------------------
  #FIGURE 10: plot est curves for MPB,  for all methods in one rep
  curve.plot.MPB <- ggplot(
    subset(erc.plots.df, 
           (Method %in% c("BKMR", "BART", "BSTARSS", "GAM", "True"))), 
    aes(X.MPB, Y.MPB)) + 
    facet_grid(erfn ~ mtype + SNR + corr) + 
    geom_line(
      data = subset(erc.plots.df, 
                    (RepID %in% c("X2"))&
                    (Method %in% c("BKMR", "BART", "BSTARSS", "GAM", "True"))), 
      aes(color = factor(Method)), size = 0.75) + 
    scale_color_manual(
      name = "Method",
      values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", "black")) + 
    xlab("MPB") + 
    ylab("Estimated value") + 
    theme_bw() + 
    theme(text = element_text(size = 14), axis.text = element_text(size = 12), 
          strip.text.x = element_text(size = 13), 
          strip.text.y = element_text(size = 13), 
          panel.grid.minor = element_blank())
  print(curve.plot.MPB)
  curve.plot.MPB <- ggdraw(curve.plot.MPB) + 
    draw_label("Model,  sparsity", x = 0.905, y = 0.9755, hjust = 0) + 
    draw_label("SNR", x = 0.905, y = (0.9755 + 0.92375)/2, hjust = 0) + 
    draw_label("Correlation", x = 0.905, y = 0.92375, hjust = 0)
  ggsave("Figure_10.pdf", plot = curve.plot.MPB, 
         height = 752/72, width = 1335/72)
  
}

##------------------------------------------------------------------------------