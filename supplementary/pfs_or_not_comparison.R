##################### COMPARE EMBEDDED METHODS WHEN USING PFS OR NOT #####################

setwd("/scratch/user/weslfletch/Supplementary") # modified to work on cluster
library("carSurv")
library("CoxBoost")
library("Coxnet")
library("doParallel")
library("doBy")
library("dplyr")
library("glmnet")
library("randomForestSRC")
library("survival")
library("WGCNA")
library("yardstick")

################################### DEFINE THE METHODS ###################################

LASSO = function(X_train, Y_train, X_test, Y_test){
  a = 1; X = as.matrix(X_train); Y = Y_train
  beta_hat =
    cv.glmnet(x = X, y = Y, family = "cox", alpha = a, type.measure = "C")$lambda.min %>%
    glmnet(x = X, y = Y, family = "cox", alpha = a, lambda = .) %>%
    coef(.) %>%
    as.numeric(.)
  ci_df = data.frame(
    true = Y_test,
    pred = as.vector(exp(-(as.matrix(X_test) %*% beta_hat)))
  )
  concordance_survival(ci_df, true, pred)$.estimate
}
ALASSO = function(X_train, Y_train, X_test, Y_test){
  beta_hat =
    Coxnet(as.matrix(X_train), Y_train, "Lasso", adaptive = TRUE, nfolds = 10)$Beta
  ci_df = data.frame(
    true = Y_test,
    pred = as.vector(exp(-(as.matrix(X_test) %*% beta_hat)))
  )
  concordance_survival(ci_df, true, pred)$.estimate
}
ENET = function(X_train, Y_train, X_test, Y_test){
  a = 0.5; X = as.matrix(X_train); Y = Y_train
  beta_hat =
    cv.glmnet(x = X, y = Y, family = "cox", alpha = a, type.measure = "C")$lambda.min %>%
    glmnet(x = X, y = Y, family = "cox", alpha = a, lambda = .) %>%
    coef(.) %>%
    as.numeric(.)
  ci_df = data.frame(
    true = Y_test,
    pred = as.vector(exp(-(as.matrix(X_test) %*% beta_hat)))
  )
  concordance_survival(ci_df, true, pred)$.estimate
}
CB = function(X_train, Y_train, X_test, Y_test){
  beta_hat = coef(CoxBoost(Y_train[,1], Y_train[,2], as.matrix(X_train)))
  ci_df = data.frame(
    true = Y_test,
    pred = as.vector(exp(-(as.matrix(X_test) %*% beta_hat)))
  )
  concordance_survival(ci_df, true, pred)$.estimate
}
RSF = function(X_train, Y_train, X_test, Y_test){
  best_nodesize = 8
  best_mtry = 403
  train_df = data.frame(X_train, time = Y_train[,1], status = Y_train[,2])
  colnames(train_df)[1:ncol(X_train)] = colnames(X_train)
  my_rsf = rfsrc(Surv(time, status) ~ ., train_df, mtry = best_mtry, nodesize = best_nodesize)
  p_obj = predict.rfsrc(my_rsf, newdata = X_test)
  step_times = p_obj$time.interest
  sf_mat = p_obj$survival
  get_median = function(k)step_times[which(sf_mat[k,]<0.5)[1]]
  ci_df = data.frame(
    true = Y_test,
    pred = sapply(1:nrow(X_test), get_median)
  )
  concordance_survival(ci_df, true, pred)$.estimate
}
sRSF = function(X_train, Y_train, X_test, Y_test){
  best_nodesize = 8
  best_mtry = 403
  screen_cols = sapply(1:ncol(X_train), \(i){
    replaceMissing(summary(coxph(Y_train~X_train[,i]))$coefficients[5], 1)<0.1
  })
  train_df = data.frame(X_train[,screen_cols], time = Y_train[,1], status = Y_train[,2])
  colnames(train_df)[1:sum(screen_cols)] = colnames(X_train)[screen_cols]
  my_rsf = rfsrc(Surv(time, status) ~ ., train_df, mtry = best_mtry, nodesize = best_nodesize)
  p_obj = predict.rfsrc(my_rsf, newdata = X_test[,screen_cols])
  step_times = p_obj$time.interest
  sf_mat = p_obj$survival
  get_median = function(k)step_times[which(sf_mat[k,]<0.5)[1]]
  ci_df = data.frame(
    true = Y_test,
    pred = sapply(1:nrow(X_test), get_median)
  )
  concordance_survival(ci_df, true, pred)$.estimate
}

########################## PERFORM REAL DATA ANALYSIS COMPONENT ##########################

# set parameters
cohort_file = "BLCA.RDS"
PFS_file = "CARS3000.RDS"
methods = list(
  LASSO = LASSO, ALASSO = ALASSO, ENET = ENET, CB = CB, RSF = RSF, sRSF = sRSF
)
metrics = c("CI", "Computation time")
nfolds = 10

standardize = function(full_df, static_idxs){
  full_df[,-static_idxs] = scale(full_df[,-static_idxs])
  full_df = full_df[,colSums(is.na(full_df))==0]
  return(full_df)
}
full_df = standardize(readRDS(cohort_file), 1:4)[,-c(3,4)]
pfs_df = full_df[,colnames(full_df) %in% c(colnames(full_df)[1:2], readRDS(PFS_file))]

set.seed(123)
n = nrow(full_df)
row_nums = sample(1:n)
fold_idxs = lapply(1:nfolds, \(f){
  sort(row_nums[(floor((f-1)*n/nfolds)+1):floor(f*n/nfolds)])
})

cluster = makeCluster(detectCores()-1)
registerDoParallel(cluster)
raw_list = foreach(i = 1:nfolds, .packages = c(
  "carSurv", "CoxBoost", "Coxnet", "doParallel", "doBy", "dplyr", "glmnet",
  "randomForestSRC", "survival", "WGCNA", "yardstick"
)) %dopar% {
  oob = fold_idxs[[i]]
  full_X_train = full_df[-oob, -(1:2)]; full_X_test = full_df[oob, -(1:2)]
  full_Y_train = full_df$surv_time[-oob]; full_Y_test = full_df$surv_time[oob]
  pfs_X_train = pfs_df[-oob, -(1:2)]; pfs_X_test = pfs_df[oob, -(1:2)]
  pfs_Y_train = pfs_df$surv_time[-oob]; pfs_Y_test = pfs_df$surv_time[oob]
  lapply(methods, \(method){
    start_time_full = Sys.time()
    ci_full = tryCatch({
      method(full_X_train, full_Y_train, full_X_test, full_Y_test)
    }, error = \(e)NA)
    comp_time_full = if (is.na(ci_full)) {
      NA
    } else {difftime(Sys.time(), start_time_full, units = "sec")}
    start_time_pfs = Sys.time()
    ci_pfs = tryCatch({
      method(pfs_X_train, pfs_Y_train, pfs_X_test, pfs_Y_test)
    }, error = \(e)NA)
    comp_time_pfs = if (is.na(ci_pfs)) {
      NA
    } else {difftime(Sys.time(), start_time_pfs, units = "sec")}
    list(c(ci_full, comp_time_full), c(ci_pfs, comp_time_pfs))
  })
}
stopCluster(cluster)

# reformat raw results
metrics_list = lapply(metrics, \(i)rep(NA, nfolds))
names(metrics_list) = metrics
methods_list = lapply(methods, \(i)metrics_list)
names(methods_list) = names(methods)
results_list = lapply(1:2, \(i)methods_list)
names(results_list) = c("Full", "PFS")
for (i in 1:2){
  for (j in 1:length(methods)){
    for (k in 1:length(metrics)){
      for (l in 1:nfolds){
        results_list[[i]][[j]][[k]][l] = raw_list[[l]][[j]][[i]][k]
} } } }
saveRDS(results_list, file = "PFSorNotCompResults.RDS")

############################## PERFORM SIMULATION COMPONENT ##############################

# set parameters
num_trials = 20
cohort_file = "BLCA.RDS"
shape = 2
scale = 300

# get full data
standardize = function(full_df, static_idxs){
  full_df[,-static_idxs] = scale(full_df[,-static_idxs])
  full_df = full_df[,colSums(is.na(full_df))==0]
  return(full_df)
}
full_df = standardize(readRDS(cohort_file), 1:4)[,-c(3,4)]

# set true beta to 10 * CoxBoost estimates
time = full_df$surv_time[,1]
status = full_df$surv_time[,2]
X = as.matrix(full_df[,-(1:2)])
n = nrow(X)
beta = 10 * coef(CoxBoost(time, status, X))
linpred = as.vector(X %*% beta)
effect_features = beta != 0

cars_3k = function(X, Y){
  cars_scores = carSurvScore(Y[,1], Y[,2], as.matrix(X))
  top_3k = which.maxn(abs(cars_scores), 3000)
  top_3k_logical = rep(FALSE, ncol(X))
  top_3k_logical[top_3k] = TRUE
  names(top_3k_logical) = names(effect_features)
  top_3k_logical[effect_features]
}

cluster = makeCluster(detectCores()-1)
registerDoParallel(cluster)
results_list = foreach(seed = 1:num_trials, .packages = c(
  "carSurv", "CoxBoost", "Coxnet", "doParallel", "doBy", "dplyr", "glmnet",
  "randomForestSRC", "survival", "WGCNA", "yardstick"
)) %dopar% {
  set.seed(seed)
  ty = scale * (-log(runif(n)) / exp(linpred))^(1/shape)
  cy = runif(n, min = quantile(ty, prob=0.5), max = quantile(ty, prob=0.9))
  X_full = as.data.frame(X)
  Y_full = Surv(apply(cbind(ty,cy), 1, min), as.numeric(ty<cy))
  cars_3k(X_full, Y_full)
}
stopCluster(cluster)

saveRDS(list(beta = beta, fold_results = results_list), file = "PFSTrueSignalsKept.RDS")

################################## CREATE RESULTS PLOTS ##################################

# first plot the proportion of each gene's retention on its signal strength
library("ggplot2")
library("latex2exp")
signals_retained = readRDS("PFSTrueSignalsKept.RDS")
beta = signals_retained$beta[signals_retained$beta!=0]
effect_features = names(beta)
fold_results = signals_retained$fold_results
gene_retained_list = lapply(1:length(effect_features), \(i){
  sapply(fold_results, \(fold)fold[i])
})
names(gene_retained_list) = effect_features
prop_retention = sapply(gene_retained_list, mean)
plot_df = data.frame(
  gene = effect_features,
  beta = beta,
  abs_beta = abs(beta),
  prop_retention = prop_retention
)
retention_plot = ggplot(plot_df) +
  geom_point(aes(abs_beta, prop_retention)) +
  ylab("Retention rate") +
  xlab(TeX("|\\beta|")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15)
  )
ggsave("RetentionBySignalStrength.jpeg", retention_plot,
       width = 2500, height = 1500, units = "px")

# next compare performances of embedded methods with and without pfs
library("ggplot2")
library("patchwork")
pfs_contrast_results = readRDS("PFSorNotCompResults.RDS")
datas = c("Without PFS", "With PFS")
methods = names(pfs_contrast_results[[1]])
metrics = names(pfs_contrast_results[[1]][[1]])
{
  ci_df = data.frame()
  ct_df = data.frame()
  for (i in 1:length(datas)){
    for (j in 1:length(methods)){
      ci_df = rbind(ci_df, data.frame(
        data = datas[i],
        method = methods[j],
        value = pfs_contrast_results[[i]][[j]][[1]]
      ))
      ct_df = rbind(ct_df, data.frame(
        data = datas[i],
        method = methods[j],
        value = pfs_contrast_results[[i]][[j]][[2]]
      ))
    } }
}

ci_plot = ggplot(ci_df) +
  geom_boxplot(aes(y = value, fill = data)) +
  facet_wrap(vars(method), scales = "free", ncol = 2) +
  ggtitle("CI") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "white", color = "black")
  )
ct_plot = ggplot(ct_df) +
  geom_boxplot(aes(y = value, fill = data)) +
  facet_wrap(vars(method), scales = "free", ncol = 2) +
  ggtitle("Computation time") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "white", color = "black")
  )
out_plot = ci_plot + ct_plot
ggsave("PFSvsNoPFSCompPlot.jpeg", out_plot, width = 2500, height = 1500, units = "px")
