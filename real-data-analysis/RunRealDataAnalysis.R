################################# RUN REAL DATA ANALYSIS #################################

setwd("YOUR-DIRECTORY-HERE") # FIXME: set to appropriate directory

# required packages
library("carSurv")
library("CoxBoost")
library("Coxnet")
library("doBy")
library("doParallel")
library("dplyr")
library("glmnet")
library("mutoss")
library("randomForestSRC")
library("survival")
library("tidyr")
library("WGCNA")
library("yardstick")

source("AnalysisControl.R")

# the methods to use in the analysis. remove an item to remove it from the study.
methods = c("LASSO", "ALASSO", "CB", "RF_md", "RF_md_screen", "BHP", "QV", "CARS_med",
            "CARS_msr")

results = run_RDA(methods = methods, files = "BLCA.RDS", prelim = "UncorMRNA.RDS")

saveRDS(results, file = "real-data-analysis-results.RDS")

# NOTES ON RESULTING OBJECT STRUCTURE:
# The results object is a tree structure of nested lists. The first level is negligible
# with one branch, the second level determines which method was evaluated, the third
# level determines which portion of the result is desired (e.g. CI, computation time, or
# the mRNA features the method selected).
