################################ RUN SIMULATION SETTING 2 ################################

# QUICK START GUIDE:
# 1. Get correct package versions (SEE README)
# 2. Get BLCA.RDS (time stamped version in repository, or run TCGA_BLCA.R)
# 3. Perform PFS (run PFS.R, or get `CARS3000.RDS` from repository)
# 4. Set parameters (main parameter for this analysis is `num_trials` on line 34)
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
num_trials = 200
cohort_file = "BLCA.RDS"
PFS_file = "CARS3000.RDS"
methods = list(
  LASSO = LASSO, ALASSO = ALASSO, ENET = ENET, CB = CB, RSF = RSF, sRSF = sRSF, BHP = BHP,
  QV = QV, CARS_MED = CARS_MED, CARS_MSR = CARS_MSR
)
metrics = c("FDR", "F1-Score", "CI", "Brier score", "RMSE", "Computation time")
shape = 2 # used in synthetic dataset generation
scale = 300 # used in synthetic dataset generation

# read and prepare full data
standardize = function(full_df, static_idxs){
  full_df[,-static_idxs] = scale(full_df[,-static_idxs])
  full_df = full_df[,colSums(is.na(full_df))==0]
  return(full_df)
}
# get cohort data and normalize features to mean 0 variance 1
full_df = standardize(readRDS(cohort_file), 1:4)
# only keep the features retained by PFS
full_df = full_df[,colnames(full_df) %in% c(colnames(full_df)[1:2], readRDS(PFS_file))]

# set true beta to 10 * CoxBoost estimates
time = full_df$surv_time[,1]
status = full_df$surv_time[,2]
X = as.matrix(full_df[,-(1:2)])
n = nrow(X)
beta = 10 * coef(CoxBoost(time, status, X))
linpred = as.vector(X %*% beta)
effect_features = beta != 0 # used for calculating FS metrics

# perform simulation loop (LONG COMPUTATION TIME, OPTIMIZED FOR CLUSTER COMPUTING)
cluster = makeCluster(detectCores()-1)
registerDoParallel(cluster)
raw_list = foreach(trial = 1:num_trials, .packages = c( # each core gets a replication
  "carSurv", "CoxBoost", "Coxnet", "doBy", "doParallel", "dplyr", "glmnet", "mutoss",
  "pec", "randomForestSRC", "survival", "tidyr", "WGCNA", "yardstick"
)) %dopar% {
  source("LASSOENET.R")
  source("ALASSO.R")
  source("CB.R")
  source("RSF.R")
  source("BHQV.R")
  source("CARS.R")
  set.seed(trial) # seed increments with each replication
  oob_rows = sort(sample(1:n, ceiling(n/10))) # obs for out-of-sample data of replication
  ty = scale * (-log(runif(n)) / exp(linpred))^(1/shape) # weibull distributed event time
  ty_test = ty[oob_rows]
  cy = runif(n, min = quantile(ty, prob=0.5), max = quantile(ty, prob=0.9)) # censor time
  X_full = as.data.frame(X)
  Y_full = Surv(apply(cbind(ty,cy), 1, min), as.numeric(ty<cy))
  # partition data into train and test
  X_train = X_full[-oob_rows,]; X_test = X_full[oob_rows,]
  Y_train = Y_full[-oob_rows]; Y_test = Y_full[oob_rows]
  # run and evaluate all methods on the dataset
  lapply(methods, \(method)tryCatch({ # apply over methods
    start_time = Sys.time() # start computation timer
    out = method(X_train, Y_train) # perform the method
    comp_time = difftime(Sys.time(), start_time, units = "sec") # stop comp timer
    surv_preds = out$predict(X_test) # used for RMSE calculation
    ci_df = data.frame(true = Y_test, inv_risk = out$inv_risk(X_test)) # needed for CI
    selected = out$selected # the method's selected features
    # needed for FDR and F1-score
    TP = sum(selected & effect_features)
    FP = sum(selected) - TP
    FN = sum(effect_features) - TP
    # return FDR, F1-Score, CI, Brier at 250, RMSE, and Computation time
    c(
      replaceMissing(FP/(FP+TP), 0),
      TP/(TP+(FP+FN)/2),
      concordance_survival(ci_df, true, inv_risk)$.estimate,
      out$bscore(X_test, Y_test, 250),
      sqrt(mean(log(surv_preds/ty_test)^2, na.rm = TRUE)),
      comp_time
    )
  }, error = \(e)rep(NA, 5)))
}
stopCluster(cluster)

# reformat and save results
metrics_list = lapply(metrics, \(i)rep(NA, num_trials))
names(metrics_list) = metrics
results_list = lapply(methods, \(i)metrics_list)
names(results_list) = names(methods)
for (i in 1:length(methods)){
  for (j in 1:length(metrics)){
    for (k in 1:num_trials){
      results_list[[i]][[j]][k] = raw_list[[k]][[i]][j]
} } }
saveRDS(results_list, file = "Setting2Results.RDS")
