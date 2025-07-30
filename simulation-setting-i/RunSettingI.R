################################ RUN SIMULATION SETTING I ################################

setwd("YOUR-DIRECTORY-HERE") # FIXME: set to appropriate directory

# all required packages
library("CoxBoost")
library("Coxnet")
library("doBy")
library("doParallel")
library("genlasso")
library("glmnet")
library("mutoss")
library("randomForestSRC")
library("survival")
library("WGCNA")
library("yardstick")

source("SimLoop.R")
source("Lasso.R")
source("Coxnet.R")
source("Coxboost.R")
source("RandomForest.R")
source("BHAdjPValue.R")
source("QValue.R")
source("Cars.R")

# NOTE: the total number of synthetic datasets created =
# length(s_vals)*length(a_vals)*length(sig_str_vals)*num_trials.
# be aware of this when changing simulation parameters.

# set simulation parameters below

# all available methods listed. remove an element to not include it in the simulation.
methods = c("LASSO", "ALASSO", "CB", "RF_md", "RF_md_scr", "BH_Proc", "QV", "CARS_med",
            "CARS_msr")
# the number of effect variables ("sparsity") to use. add, remove, or change as desired.
s_vals = c(20, 50, 100)
# the levels of predictor correlation to use. add, remove, or change as desired.
a_vals = c(0, 0.5)
# the levels of signal strength to use. add, remove, or change as desired.
sig_str_vals = c(0.5, 1, 2)
# the number of unique synthetic datasets to create for every unique combination of
# s_vals, a_vals, and sig_str_vals values.
num_trials = 200
# number of observations in each synthetic dataset
n = 300
# total number of features in each synthetic dataset (including effect variables)
p = 1000

# run the study and save the results

results = run_setting_i(methods, s_vals, a_vals, sig_str_vals, num_trials, n, p)
saveRDS(results, file = "setting-i-results.RDS")

# NOTE: `results` is a tree of nested lists. The first level determines the method, the
# second determines the sparsity level (s), the third determines the level of predictor
# correlation (alpha), the fourth determines the strength of signal (gamma), and the fifth
# determines metric (e.g. F1-Score). Under all these levels are vectors of length
# `num_trials`.
