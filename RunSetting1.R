################################ RUN SIMULATION SETTING 1 ################################

# QUICK START GUIDE:
# 1. Get correct package versions (SEE README)
# 2. Set parameters (lines 32-38)
# 3. Modify working directory file path (line 9)
# 4. Run this file top to bottom

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
n = 300
n_train = 200
p = 1000
s_vals = c(20, 50, 100)
alpha_vals = c(0, 0.5)
sig_str_vals = c(0.5, 1, 2)
methods = list(
  LASSO = LASSO, ALASSO = ALASSO, ENET = ENET, CB = CB, RSF = RSF, sRSF = sRSF, BHP = BHP,
  QV = QV, CARS_MED = CARS_MED, CARS_MSR = CARS_MSR
)
metrics = c("FDR", "F1-Score", "CI", "Brier score", "RMSE", "Computation time")

# helper to perform all methods on a single dataset
run_one_dataset = function(seed, A, s, sig_str){
  # simulate and partition data into training and testing splits
  set.seed(seed) # seeds increment with each iteration
  beta = rep(0, p)
  beta[1:floor(s/2)] = sig_str
  beta[(floor(s/2)+1):s] = -sig_str
  z = matrix(rnorm(n*p), nrow = p) # N(0,I) (uncorrelated features)
  x = t(A %*% z) # correlate features with Cholesky decomposition
  ty = rexp(n, exp(-(x %*% as.matrix(beta)))) # time follows exponential distribution
  ty_test = ty[-(1:n_train)]
  cy = runif(n, quantile(ty, prob=0.5), quantile(ty, prob=0.90)) # censoring time
  X_full = as.data.frame(x)
  Y_full = Surv(apply(cbind(ty,cy), 1, min), as.numeric(ty<cy))
  # partition into train and test data
  X_train = X_full[1:n_train,]; X_test = X_full[-(1:n_train),]
  Y_train = Y_full[1:n_train]; Y_test = Y_full[-(1:n_train)]
  # run and evaluate all methods on the dataset
  lapply(methods, \(method)tryCatch({
    start_time = Sys.time() # start comp timer
    out = method(X_train, Y_train) # perform the method
    comp_time = difftime(Sys.time(), start_time, units = "sec") # stop comp timer
    surv_preds = out$predict(X_test) # used for calculating RMSE
    ci_df = data.frame(true = Y_test, inv_risk = out$inv_risk(X_test)) # needed for CI
    selected = which(out$selected) # method's selected features
    # need for FDR and F1-score calculations
    TP = sum(selected <= s)
    FP = length(selected) - TP
    FN = s - TP
    # return FDR, F1-Score, CI, Brier at 0.5, RMSE, and Computation time
    c(
      replaceMissing(FP/(FP+TP), 0),
      TP/(TP+(FP+FN)/2),
      concordance_survival(ci_df, true, inv_risk)$.estimate,
      out$bscore(X_test, Y_test, 0.5),
      sqrt(mean(log(surv_preds/ty_test)^2, na.rm = TRUE)),
      comp_time
    )
  }, error = \(e)rep(NA, 5)))
}

# helper to get Cholesky decomposition for correlating features
get_chol_mat = function(p, alpha){
  sigma = matrix(alpha, p, p)
  diag(sigma) = 1
  return(t(chol(sigma)))
}

# perform simulation loop (OPTIMIZED FOR CLUSTER COMPUTING, TAKES MANY HOURS)
cluster = makeCluster(detectCores()-1)
registerDoParallel(cluster)
raw_list = foreach(trial = 1:num_trials, .packages = c(
  "carSurv", "CoxBoost", "Coxnet", "doBy", "doParallel", "dplyr", "glmnet", "mutoss",
  "pec", "randomForestSRC", "survival", "tidyr", "WGCNA", "yardstick"
)) %dopar% {
  source("LASSOENET.R")
  source("ALASSO.R")
  source("CB.R")
  source("RSF.R")
  source("BHQV.R")
  source("CARS.R")
  lapply(alpha_vals, \(alpha){
    A = get_chol_mat(p, alpha)
    lapply(s_vals, \(s){
      lapply(sig_str_vals, \(sig_str){
        run_one_dataset(trial, A, s, sig_str)
  })})})
}
stopCluster(cluster)

# reformat and save results
metrics_list = lapply(metrics, \(i)rep(NA, num_trials))
names(metrics_list) = metrics
sig_str_list = lapply(sig_str_vals, \(i)metrics_list)
names(sig_str_list) = as.character(sig_str_vals)
alpha_list = lapply(alpha_vals, \(i)sig_str_list)
names(alpha_list) = as.character(alpha_vals)
s_list = lapply(s_vals, \(i)alpha_list)
names(s_list) = as.character(s_vals)
results_list = lapply(methods, \(i)s_list)
names(results_list) = names(methods)
for (i in 1:length(methods)){
  for (j in 1:length(s_vals)){
    for (k in 1:length(alpha_vals)){
      for (l in 1:length(sig_str_vals)){
        for (m in 1:length(metrics)){
          for (o in 1:num_trials){
            results_list[[i]][[j]][[k]][[l]][[m]][o] = raw_list[[o]][[k]][[j]][[l]][[i]][m]
} } } } } }
saveRDS(results_list, file = "Setting1Results.RDS")
