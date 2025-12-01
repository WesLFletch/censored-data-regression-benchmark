######################### COMPARE OOB AND HELD-OUT METRICS IN RSF ########################

library("doParallel")
library("randomForestSRC")
library("survival")
library("yardstick")

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
  a_priori_results = readRDS("RSFTuningParams.RDS")
  a_priori_nodesize = ceiling(mean(a_priori_results$nodesize))
  a_priori_mtry = ceiling(mean(a_priori_results$mtry))
  train_df = data.frame(X_train, time = Y_train[,1], status = Y_train[,2])
  my_rsf = rfsrc(Surv(time, status) ~ ., train_df,
                 mtry = a_priori_mtry, nodesize = a_priori_nodesize)
  oob_ci = 1 - my_rsf$err.rate[length(my_rsf$err.rate)]
  p_obj = predict.rfsrc(my_rsf, newdata = X_test)
  step_times = p_obj$time.interest
  sf_mat = p_obj$survival
  get_median = function(k)step_times[which(sf_mat[k,]<0.5)[1]]
  heldout_ci_df = data.frame(
    true = Y_test,
    preds = sapply(1:nrow(X_test), get_median)
  )
  heldout_ci = concordance_survival(heldout_ci_df, true, preds)$.estimate
  c(oob = oob_ci, heldout = heldout_ci)
}
stopCluster(cluster)

results_list = lapply(1:2, \(i)rep(NA, num_trials))
names(results_list) = c("OOB", "Held-out")
for (i in 1:num_trials){
  results_list[[1]][i] = raw_list[[i]][1]
  results_list[[2]][i] = raw_list[[i]][2]
}
saveRDS(results_list, file = "RSFOOBvsHeldoutResults.RDS")

###################################### READ RESULTS ######################################

results_list = readRDS("RSFOOBvsHeldoutResults.RDS")
out_matrix = t(sapply(results_list, summary))
write.csv(out_matrix, file = "RSFOOBvsHeldoutResults.csv")
