################# COMPARE RSF WITH DEFAULT, A PRIORI, AND RUNTIME TUNING #################

setwd("/scratch/user/weslfletch/Supplementary") # modified to work on cluster
library("doParallel")
library("randomForestSRC")
library("survival")
library("yardstick")

rsf_default = function(X_train, Y_train, X_test, Y_test){
  train_df = data.frame(X_train, time = Y_train[,1], status = Y_train[,2])
  my_rsf = rfsrc(Surv(time, status) ~ ., train_df)
  p_obj = predict.rfsrc(my_rsf, newdata = X_test)
  step_times = p_obj$time.interest
  sf_mat = p_obj$survival
  get_median = function(k)step_times[which(sf_mat[k,]<0.5)[1]]
  ci_df = data.frame(
    true = Y_test,
    preds = sapply(1:nrow(X_test), get_median)
  )
  concordance_survival(ci_df, true, preds)$.estimate
}
rsf_a_priori = function(X_train, Y_train, X_test, Y_test){
  a_priori_results = readRDS("RSFTuningParams.RDS")
  a_priori_nodesize = ceiling(mean(a_priori_results$nodesize))
  a_priori_mtry = ceiling(mean(a_priori_results$mtry))
  train_df = data.frame(X_train, time = Y_train[,1], status = Y_train[,2])
  my_rsf = rfsrc(Surv(time, status) ~ ., train_df,
                 mtry = a_priori_mtry, nodesize = a_priori_nodesize)
  p_obj = predict.rfsrc(my_rsf, newdata = X_test)
  step_times = p_obj$time.interest
  sf_mat = p_obj$survival
  get_median = function(k)step_times[which(sf_mat[k,]<0.5)[1]]
  ci_df = data.frame(
    true = Y_test,
    preds = sapply(1:nrow(X_test), get_median)
  )
  concordance_survival(ci_df, true, preds)$.estimate
}
rsf_runtime = function(X_train, Y_train, X_test, Y_test){
  train_df = data.frame(X_train, time = Y_train[,1], status = Y_train[,2])
  tune_obj = tune(Surv(time, status) ~ ., train_df)$optimal
  opt_nodesize = unname(tune_obj[1])
  opt_mtry = unname(tune_obj[2])
  my_rsf = rfsrc(Surv(time, status) ~ ., train_df,
                 mtry = opt_mtry, nodesize = opt_nodesize)
  p_obj = predict.rfsrc(my_rsf, newdata = X_test)
  step_times = p_obj$time.interest
  sf_mat = p_obj$survival
  get_median = function(k)step_times[which(sf_mat[k,]<0.5)[1]]
  ci_df = data.frame(
    true = Y_test,
    preds = sapply(1:nrow(X_test), get_median)
  )
  concordance_survival(ci_df, true, preds)$.estimate
}
{
  num_trials = 20
  n = 300
  n_train = 200
  p = 1000
  s = 50
  alpha = 0.5
  sig_str = 1
  sigma = matrix(alpha, p, p)
  diag(sigma) = 1
  A = t(chol(sigma))
  methods = list(
    rsf_default = rsf_default, rsf_a_priori = rsf_a_priori, rsf_runtime = rsf_runtime
  )
  metrics = c("CI", "Computation time")
}
cluster = makeCluster(detectCores()-1)
registerDoParallel(cluster)
raw_list = foreach(trial = 1:num_trials, .packages = c(
  "doParallel", "randomForestSRC", "survival", "yardstick"
)) %dopar% {
  set.seed(trial)
  beta = rep(0, p)
  beta[1:floor(s/2)] = sig_str
  beta[(floor(s/2)+1):s] = -sig_str
  z = matrix(rnorm(n*p), nrow = p)
  x = t(A %*% z)
  ty = rexp(n, exp(-(x %*% as.matrix(beta))))
  ty_test = ty[-(1:n_train)]
  cy = runif(n, quantile(ty, prob=0.5), quantile(ty, prob=0.90))
  X_full = as.data.frame(x)
  Y_full = Surv(apply(cbind(ty,cy), 1, min), as.numeric(ty<cy))
  X_train = X_full[1:n_train,]; X_test = X_full[-(1:n_train),]
  Y_train = Y_full[1:n_train]; Y_test = Y_full[-(1:n_train)]
  lapply(methods, \(method){
    start_time = Sys.time()
    out = method(X_train, Y_train, X_test, Y_test)
    comp_time = difftime(Sys.time(), start_time, units = "sec")
    c(
      out,
      comp_time
    )
  })
}
stopCluster(cluster)
{
  metrics_list = lapply(metrics, \(i)rep(NA, num_trials))
  names(metrics_list) = metrics
  results_list = lapply(methods, \(i)metrics_list)
  names(results_list) = names(methods)
  for (i in 1:length(methods)){
    for (j in 1:length(metrics)){
      for (k in 1:num_trials){
        results_list[[i]][[j]][k] = raw_list[[k]][[i]][j]
  } } }
}
saveRDS(results_list, "RSFTuningComparisonResults.RDS")

################################# CREATE TABLE OF RESULTS ################################

# extract individual metric results
results = readRDS("RSFTuningComparisonResults.RDS")
ci_vals = lapply(results, \(r)r[[1]])
names(ci_vals) = names(results)
ct_vals = lapply(results, \(r)r[[2]])
names(ct_vals) = names(results)

# create metric dataframes
{
  ci_df = data.frame()
  for (i in 1:length(ci_vals)){
    ci_df = rbind(ci_df, data.frame(
      method = names(ci_vals)[i],
      value = ci_vals[[i]]
    ))
  }
}
{
  ct_df = data.frame()
  for (i in 1:length(ct_vals)){
    ct_df = rbind(ct_df, data.frame(
      method = names(ct_vals)[i],
      value = ct_vals[[i]]
    ))
  }
}
# merge results
{
  results_table = as.data.frame(matrix(NA, length(results), length(results[[1]])))
  for (i in 1:length(results)){
    for (j in 1:length(results[[1]])){
      results_table[i,j] = paste0(
        round(mean(results[[i]][[j]], na.rm = TRUE), 2),
        "(",
        round(sd(results[[i]][[j]], na.rm = TRUE), 2),
        ")"
      )
    } }
  rownames(results_table) = names(results)
  colnames(results_table) = names(results[[1]])
}

print(results_table)
