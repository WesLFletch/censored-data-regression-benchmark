########################## PERFORM PRELIMINARY FEATURE SELECTION #########################

setwd("YOUR-DIRECTORY-HERE") # FIXME: set to appropriate directory

# required packages
library("caret")
library("doParallel")
library("dplyr")
library("glmnet")
library("survival")

set.seed(1)

# load in data and normalize mRNA features
standardize = function(full_df, static_idxs){
  full_df[,-static_idxs] = scale(full_df[,-static_idxs])
  full_df = full_df[,colSums(is.na(full_df))==0]
  return(full_df)
}
full_df = standardize(readRDS("BLCA.RDS"), 1:4)
surv_time = full_df$surv_time

# remove features with IQR == 0
mrna_mat = full_df[,c(
  rep(FALSE, 4), # remove all non-mrna features
  apply(full_df[,-(1:4)], 2, \(i)IQR(i)!=0) # remove features with IQR == 0
)] %>% as.matrix(.) # cast to matrix

# remove features with high pairwise correlations
cor_mat = cor(mrna_mat)
high_pairwise_cors = findCorrelation(cor_mat, cutoff = 0.8)
mrna_mat = mrna_mat[,-high_pairwise_cors]

# get goodness-of-fit scores using cluster computing
cluster = makeCluster(detectCores() - 1)
registerDoParallel(cluster)
results_list = foreach(i = 1:ncol(mrna_mat), .packages = c(
  "doParallel", "dplyr", "glmnet", "survival")) %dopar% {
    tryCatch({
      test = coxph(surv_time ~ mrna_mat[,i]) # if it throws warning, infinite coefficient
      cv_coeffs = coef(cv.glmnet(mrna_mat[,-i], mrna_mat[,i], alpha = 1), s = "lambda.min")
      selected_features = (rownames(cv_coeffs)[which(cv_coeffs != 0)])[-1]
      if (length(selected_features) == 0) stop("no features selected by LASSO")
      my_lm = summary(lm(mrna_mat[,i] ~ mrna_mat[,selected_features]))
      c(my_lm$r.squared, my_lm$adj.r.squared, length(selected_features))
    }, warning = function(w) rep(NA, 3), error = function(e) rep(NA, 3))
  }
stopCluster(cluster)

# reformat the results
results_df = as.data.frame(t(sapply(results_list, I)))
colnames(results_df) = c("r_sq", "adj_r_sq", "num_selected_features")
results_df$gene = colnames(mrna_mat)

# save the results
saveRDS(results_df, file = "PrelimFSResults.RDS")
