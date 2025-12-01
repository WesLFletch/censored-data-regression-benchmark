############################# RUN REVISED REAL DATA ANALYSIS #############################

# QUICK START GUIDE:
# 1. Get correct package versions (SEE README)
# 2. Get BLCA.RDS (time stamped version in repository, or run TCGA_BLCA.R)
# 3. Perform PFS (run PFS.R, or get `CARS3000.RDS` from repository)
# 4. Set parameters (main parameter for this analysis is `nfolds` on line 49)
# 5. Modify working directory file path (line 11)
# 6. Run this file top to bottom

setwd("/FILEPATH/") # MODIFY
library("carSurv")
library("CoxBoost")
library("Coxnet")
library("doBy")
library("doParallel")
library("dplyr")
library("glmnet")
library("mutoss")
library("pec")
library("randomForestSRC")
library("survival")
library("tidyr")
library("WGCNA")
library("yardstick")
source("LASSOENET.R")
source("ALASSO.R")
source("CB.R")
source("RSF.R")
source("BHQV.R")
source("CARS.R")

# set parameters
cohort_file = "BLCA.RDS"
PFS_file = "CARS3000.RDS"
pos_controls_char = c(
  "IGF2BP1", "MMP19", "MMP20", "MMP27", "MMP8", "FOXA3", "FOXF1", "FOXF2", "FOXG1",
  "FOXI1", "FOXK2", "FOXL2", "FOXN4", "FOXP3", "FOXR1", "FGF14", "FGF22", "FGF5",
  "FGFRL1", "TOP1", "TOP3B", "PSMD4", "PSMD13", "SQLE", "PLA2G3", "PLA2G4A", "PLA2G4D",
  "ATG5", "IGF1R", "DLEU1"
)
methods = list(
  LASSO = LASSO, ALASSO = ALASSO, ENET = ENET, CB = CB, RSF = RSF, sRSF = sRSF, BHP = BHP,
  QV = QV, CARS_MED = CARS_MED, CARS_MSR = CARS_MSR
)
metrics = c(
  "Selected mRNAs", "True positives", "CI", "Brier365", "Brier1000", "Computation time"
)
nfolds = 10

# read and prepare full data
standardize = function(full_df, static_idxs){
  full_df[,-static_idxs] = scale(full_df[,-static_idxs])
  full_df = full_df[,colSums(is.na(full_df))==0]
  return(full_df)
}
# load data and standardize covaraites to mean 0 variance 1
full_df = standardize(readRDS(cohort_file), 1:4)
# remove the columns that were not retained by PFS
full_df = full_df[,colnames(full_df) %in% c(colnames(full_df)[1:2], readRDS(PFS_file))]
# get binary vector for use in calculating true positives
pos_controls_bin = colnames(full_df)[-(1:2)] %in% pos_controls_char

# NOTE: This is the beginning of the parallel computing over multiple cores for
# efficiency. Expect runtimes exceeding an hour on high-performance hardware.
cluster = makeCluster(detectCores()-1)
registerDoParallel(cluster)
raw_list = foreach(m = 1:length(methods), .packages = c( # each core gets a method
  "carSurv", "CoxBoost", "Coxnet", "doBy", "doParallel", "dplyr", "glmnet", "mutoss",
  "pec", "randomForestSRC", "survival", "tidyr", "WGCNA", "yardstick"
)) %dopar% {
  source("LASSOENET.R")
  source("ALASSO.R")
  source("CB.R")
  source("RSF.R")
  source("BHQV.R")
  source("CARS.R")
  method = methods[[m]] # get individual method
  set.seed(123) # fixed seed assignments so all methods get the same data
  n = nrow(full_df)
  p = ncol(full_df)-2
  row_nums = sample(1:n)
  fold_idxs = lapply(1:nfolds, \(f){ # get fold assignments for outer loop
    sort(row_nums[(floor((f-1)*n/nfolds)+1):floor(f*n/nfolds)])
  })
  # evaluate the method on all folds
  inner_results = lapply(fold_idxs, \(oob)tryCatch({ # apply over folds
    # perform train-test split
    X_train = full_df[-oob,-(1:2)]; X_test = full_df[oob, -(1:2)]
    Y_train = full_df$surv_time[-oob]; Y_test = full_df$surv_time[oob]
    start_time = Sys.time() # begin counting computation time
    out = method(X_train, Y_train) # perform the method on training data
    comp_time = difftime(Sys.time(), start_time, units = "sec") # stop comp timer
    ci_df = data.frame(true = Y_test, inv_risk = out$inv_risk(X_test)) # needed for CI
    selected = out$selected # get binary vector of selected features
    # return individual fold results
    list(
      selected, # selected features
      concordance_survival(ci_df, true, inv_risk)$.estimate, # out-of-sample CI
      out$bscore(X_test, Y_test, c(365, 1000)), # Brier score at 365 and 1,000 days
      comp_time
    )
  }, error = \(e)rep(NA, 4)))
  lapply(inner_results, \(f)list(   # list of metrics vector and set of selected features
    c(
      sum(f[[1]]),                  # num selected features for fold
      sum(f[[1]]&pos_controls_bin), # num true positives for fold
      f[[2]],                       # out-of-sample concordance index for fold
      f[[3]][1],                    # out-of-sample brier score at t=365 for fold
      f[[3]][2],                    # out-of-sample brier score at t=1000 for fold
      f[[4]]                        # computation time for fold
    ),
    f[[1]] # set of selected features in fold
  ))
}
stopCluster(cluster)

# reformat the raw list to make a better organized output
metrics_list = lapply(metrics, \(i)rep(NA, nfolds))
names(metrics_list) = metrics
results_list = lapply(methods, \(i)metrics_list)
names(results_list) = names(methods)
for (i in 1:length(methods)){
  for (j in 1:length(metrics)){
    for (k in 1:nfolds){
      results_list[[i]][[j]][k] = raw_list[[i]][[k]][[1]][j]
} } }
saveRDS(results_list, file = "RDAResults.RDS")
# reformat to extract the sets of selected features
folds_list = lapply(1:nfolds, \(i)NA)
selected_features_list = lapply(methods, \(i)folds_list)
names(selected_features_list) = names(methods)
for (i in 1:length(methods)){
  for (j in 1:nfolds){
    selected_features_list[[i]][[j]] = raw_list[[i]][[j]][[2]]
} }
saveRDS(selected_features_list, file = "RDASelectedFeatures.RDS")
