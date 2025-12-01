############################## PERFORM PFS USING CARS METHOD #############################

# QUICK START GUIDE:
# Run this file top to bottom. The output file `CARS3000.RDS` is also available in the
# repository.

setwd("/FILEPATH/") # MODIFY
library("carSurv")
library("doBy")
library("dplyr")
library("survival")

# prep data
standardize = function(full_df, static_idxs){
  full_df[,-static_idxs] = scale(full_df[,-static_idxs])
  full_df = full_df[,colSums(is.na(full_df))==0]
  return(full_df)
}
full_df = standardize(readRDS("BLCA.RDS"), 1:4)
obs_time = full_df$surv_time[,1]
obs_event = full_df$surv_time[,2]
X = full_df[,-(1:4)]

# get CARS scores
cars_scores = carSurvScore(obs_time, obs_event, as.matrix(X))
cars_abs = abs(cars_scores)

# get 3,000 best features
best_idxs = which.maxn(cars_abs, 3000)
X_red = X[,best_idxs]
chosen_mrnas = colnames(X_red)

# save results
saveRDS(chosen_mrnas, "CARS3000.RDS")
