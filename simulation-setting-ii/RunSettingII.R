############################## RUN SETTING II OF SIMULATIONS #############################

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

source("LassoRDSim.R")
source("AlassoRDSim.R")
source("CoxBoostRDSim.R")
source("RandomForestRDSim.R")
source("BHProcedureRDSim.R")
source("QvalueRDSim.R")
source("CarsRDSim.R")

# SET SIMULATION PARAMETERS

# when modifying num_reps or quantiles, be mindful that the total number of data sets
# analyzed = num_reps*length(quantiles).

# file names. DO NOT MODIFY.
out_file_name = "setting-ii-results.RDS"
prelim = "UncorMRNA.RDS"
# the number of times to simulate survival times (= the number of data sets to analyze for
# each approx. censorship case)
num_reps = 200
# the methods to perform on each data set. DO NOT MODIFY.
methods = list(LASSO, ALASSO, CB, RF_md, RF_md_screen, BHP, QV, CARS_med, CARS_msr)
# the quantiles to use in the make_df() function to achieve the different approx.
# censorship cases, adding more elements increases the number of data sets to analyze in
# the study
quantiles = list("10%cens" = c(0.3, 0.95), "30%cens" = c(0.1, 0.8))
# the factor to scale beta by, 6.35 is used to achieve having all nonzero absolute beta
# values between 0.1 and 0.69.
beta_factor = 6.35
# the names of the metrics used. DO NOT MODIFY.
metrics = c("CI", "RMSE", "F1-Score", "FDR", "comp_time")

# FOR SURVIVAL TIME GENERATION

make_df = function(X, beta, seed, quantiles){
  set.seed(seed) # ensure unique shuffles of observations with each survival simulation
  shape = 2
  scale = 300
  n = nrow(X)
  # shuffle observations so that out-of-sample observations are varied
  X = X[sample(1:n, n),]
  linpred = as.vector(X %*% beta)
  # survival function is Weibull with shape 2 and scale 300
  event_time = scale * (-log(runif(n)) / exp(linpred))^(1/shape)
  # using an arbitrary distribution for censoring that is independent of X
  q_lower = quantile(event_time, quantiles[1])
  q_upper = quantile(event_time, quantiles[2])
  censoring_time = runif(n, min = q_lower, max = q_upper)
  # determine observed times and censoring indicator
  obs_time = pmin(event_time, censoring_time)
  cens = as.integer(event_time <= censoring_time)
  # construct new df using simulated times-to-event and censoring times
  surv_obj = Surv(obs_time, cens)
  new_df = data.frame(surv_time = surv_obj, X)
  return(list(new_df = new_df, true_tte = event_time))
}

# GET X AND "REAL" BETA VECTOR FOR USE IN SIMULATING SURVIVAL, APPROXIMATED BY COXBOOST

full_df = readRDS("BLCA.RDS")
full_df = full_df[,colnames(full_df) %in% c(colnames(full_df)[1:4], readRDS(prelim))]
full_df[,-(1:4)] = scale(full_df[,-(1:4)])
X = as.matrix(full_df[,-(1:4)])
beta = coef(CoxBoost(full_df$surv_time[,1], full_df$surv_time[,2], X)) * beta_factor
true_sig_preds = names(beta[beta != 0])

# PERFORM MAIN SIMULATION LOOP

cluster = makeCluster(detectCores() - 1)
registerDoParallel(cluster)
raw_results = foreach(rep = 1:num_reps, .packages = c(
  "doBy", "doParallel", "dplyr", "survival", "tidyr", "yardstick")) %dopar% {
    source("LassoRDSim.R");source("AlassoRDSim.R");source("CoxBoostRDSim.R")
    source("RandomForestRDSim.R");source("BHProcedureRDSim.R");source("QvalueRDSim.R")
    source("CarsRDSim.R")
    lapply(quantiles, \(q_vec){
      data_list = make_df(X, beta, seed = rep, quantiles = q_vec)
      lapply(methods, \(method){
        tryCatch({
          method(data_list$new_df, data_list$true_tte, true_sig_preds)
        }, error = \(e)rep(NA, 5))
      })
    })
  }
stopCluster(cluster)

# REFORMAT THE RAW RESULTS
# current structure: raw_results[[1:num_reps]][[1:length(quantiles)]][[1:length(methods)]][1:length(metrics)]
# turn it into: results_list[[1:length(methods)]][[1:length(quantiles)]][[1:length(metrics)]][1:num_reps]
reps_vec = rep(NA, num_reps)
metrics_list = lapply(1:length(metrics), \(i)reps_vec)
names(metrics_list) = metrics
quantiles_list = lapply(1:length(quantiles), \(i)metrics_list)
names(quantiles_list) = names(quantiles)
methods_list = lapply(1:length(methods), \(i)quantiles_list)
names(methods_list) = c("LASSO", "ALASSO", "CB", "RF_md",
                        "RF_md_screen", "BHP", "QV", "CARS_med", "CARS_msr")
results_list = methods_list
for (i in 1:length(methods)){
  for (j in 1:length(quantiles)){
    for (k in 1:length(metrics)){
      for (l in 1:num_reps){
        results_list[[i]][[j]][[k]][l] = raw_results[[l]][[j]][[i]][k]
      }
    }
  }
}
saveRDS(results_list, file = out_file_name)

# NOTES OF THE STRUCTURE OF THE RESULTS OBJECT:
# results_list is a tree of nested lists. The first level determines the method, the
# second determines the approx. censorship rate examined, and the third determines which
# evaluation metric was used. Each leaf is a vector of length num_reps.
