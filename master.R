##------------------------------------------------------------------------------
# Code for: "Performance of variable and function selection methods for
#            estimating the non-linear health effects of correlated
#            chemical mixtures: a simulation study."
# Authors:   Nina Lazarevic, Luke D. Knibbs, Peter D. Sly, Adrian G. Barnett
# Preprint:  https://arxiv.org/abs/1908.01583
#
# Written by Nina Lazarevic using R 3.4.3
##------------------------------------------------------------------------------

options(java.parameters = "-Xmx20g") 

#load required packages
library(copula)       #version 0.999-18
library(copulaedas)   #version 1.4.2
library(gsl)          #version 1.9-10.3
library(parallel)     #version 3.4.3
library(doParallel)   #version 1.0.11
library(bkmr)         #version 0.2.0
library(bartMachine)  #version 1.2.3
library(spikeSlabGAM) #version 1.1-14
library(glmnet)       #version 2.0-13
library(mgcv)         #version 1.8-23
library(earth)        #version 4.6.0
library(ggplot2)      #version 2.2.1, for plotting only
library(scales)       #version 0.5.0, for plotting only
library(cowplot)      #version 0.9.4, for plotting only

#required functions
source("simulate_data.R")
source("method_estimation.R")
source("results_tables_plots.R")

#import NHANES data from CSV file from Stata
nhanes_data <- read.csv("NHANES_ExposureData_17.csv", header = T)
    
#simulate data
sim.data <- simulate_data(x = nhanes_data)

#run simulation for each method 
#  ***please note storage space requirement if saving model objects***
method_estimation(sim.data, select_method = "bkmr", num.cores = 18,
                  save_model_obj = TRUE) #if TRUE requires 62Gb space
method_estimation(sim.data, select_method = "bart", num.cores = 18,
                  save_model_obj = TRUE) #if TRUE requires 11Gb space
method_estimation(sim.data, select_method = "star", num.cores = 18,
                  save_model_obj = TRUE) #if TRUE requires 86Gb space
method_estimation(sim.data, select_method = "lasso", num.cores = 18,
                  save_model_obj = TRUE) #if TRUE requires 240Mb space
method_estimation(sim.data, select_method = "gamdp", num.cores = 18,
                  save_model_obj = TRUE) #if TRUE requires 8Gb space
method_estimation(sim.data, select_method = "gamts", num.cores = 18,
                  save_model_obj = TRUE) #if TRUE requires 8Gb space
method_estimation(sim.data, select_method = "mars", num.cores = 18,
                  save_model_obj = TRUE) #if TRUE requires 83Gb space
method_estimation(sim.data, select_method = "gam", num.cores = 18,
                  save_model_obj = TRUE) #if TRUE requires 1.6Gb space

#prepare tables and plots of results, including:
#  "VarSelection_<date>.csv"
#  "MSE25_OBSCORR_<date>.csv" & "MSE25_HALFCORR_<date>.csv"
#  "Coverage25_OBSCORR_<date>.csv" & "Coverage25_HALFCORR_<date>.csv"
#  "Figure_3.pdf" to "Figure_10.pdf"
results_tables_plots(
  sim.data,
  res.halfcorr.bkmr, mse.25.halfcorr.bkmr, coverage.25.halfcorr.bkmr,
  res.obscorr.bkmr, mse.25.obscorr.bkmr, coverage.25.obscorr.bkmr,
  res.halfcorr.bart, mse.25.halfcorr.bart, coverage.25.halfcorr.bart,
  res.obscorr.bart, mse.25.obscorr.bart, coverage.25.obscorr.bart,
  res.halfcorr.star, mse.25.halfcorr.star, coverage.25.halfcorr.star,
  res.obscorr.star, mse.25.obscorr.star, coverage.25.obscorr.star,
  res.halfcorr.gamdp, mse.25.halfcorr.gamdp, coverage.25.halfcorr.gamdp,
  res.obscorr.gamdp, mse.25.obscorr.gamdp, coverage.25.obscorr.gamdp,
  res.halfcorr.gamts, mse.25.halfcorr.gamts, coverage.25.halfcorr.gamts,
  res.obscorr.gamts, mse.25.obscorr.gamts, coverage.25.obscorr.gamts,
  res.halfcorr.mars, mse.25.halfcorr.mars, coverage.25.halfcorr.mars,
  res.obscorr.mars, mse.25.obscorr.mars, coverage.25.obscorr.mars,
  res.halfcorr.gam, mse.25.halfcorr.gam, coverage.25.halfcorr.gam,
  res.obscorr.gam, mse.25.obscorr.gam, coverage.25.obscorr.gam,
  res.halfcorr.lasso, 
  res.obscorr.lasso
)

##------------------------------------------------------------------------------
