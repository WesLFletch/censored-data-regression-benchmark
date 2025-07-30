###################################### OUTER CONTROL #####################################

run_setting_i = function(methods = NULL, s_vals = NULL, a_vals = NULL,
                         sig_str_vals = NULL, num_trials = 200, n = 300, p = 1000){
  # make raw_list from the ground up
  num_metrics = 5
  sig_str_list = lapply(1:length(sig_str_vals), \(i) rep(NA, num_metrics))
  alpha_list = lapply(1:length(a_vals), \(i) sig_str_list)
  raw_list = lapply(1:length(s_vals), \(i) alpha_list)
  # run all requested jobs
  for (i in 1:length(s_vals)){
    for (j in 1:length(a_vals)){
      A = get_chol_mat(p, a_vals[j]) # do once for each alpha for efficiency
      for (k in 1:length(sig_str_vals)){
        cluster = makeCluster(detectCores() - 1) # create parallel computing cluster
        registerDoParallel(cluster) # allow parallel computing
        raw_list[[i]][[j]][[k]] =
          foreach(l = 1:num_trials, .packages = c("survival")) %dopar% {
          source("SimLoop.R")
          run_one_dataset(methods, l, A, n, p, s_vals[i], sig_str_vals[k])}
        stopCluster(cluster) # destroy parallel computing cluster
      }
    }
  }
  return(
    reformat_raw_list(raw_list, methods, s_vals, a_vals, sig_str_vals, num_trials)
  )
}

run_one_dataset = function(methods, seed, A, n, p, s, sig_str){
  # make data
  num_metrics = 5
  set.seed(seed = seed)
  beta = rep(0, p)
  beta[1:floor(s/2)] = sig_str
  beta[(floor(s/2)+1):s] = -sig_str
  z = matrix(rnorm(n*p), nrow = p)
  x = t(A %*% matrix(rnorm(n*p), nrow = p))
  ty = rexp(n, exp(-(x %*% as.matrix(beta))))
  cy = runif(n, quantile(ty, prob=0.5), quantile(ty, prob=0.90))
  df = as.data.frame(x)
  df$surv_time = Surv(apply(cbind(ty,cy), 1, min), as.numeric(ty<cy))
  data = list(df = df, beta = beta, true_tte = ty)
  # run methods
  results = lapply(methods, \(m_name){
    method = names_to_methods[[m_name]]
    result = rep(NA, num_metrics)
    tryCatch({
      result = method(data, n, p, s)
    }, error = function(msg){})
    return(result)
  })
  names(results) = methods
  return(results)
}

reformat_raw_list = function(raw_list, methods, s_vals, a_vals, sig_str_vals, num_trials){
  # make methods_list from the ground up
  num_metrics = 5
  metrics_list = lapply(1:num_metrics, \(i) rep(NA, num_trials))
  names(metrics_list) =
    c("CI", "F1", "FDR", "comp_time", "rmse1")
  sig_str_list = lapply(1:length(sig_str_vals), \(i) metrics_list)
  names(sig_str_list) = paste("sig_str_", as.character(sig_str_vals), sep = "")
  alpha_list = lapply(1:length(a_vals), \(i) sig_str_list)
  names(alpha_list) = paste("alpha_", as.character(a_vals), sep = "")
  s_list = lapply(1:length(s_vals), \(i) alpha_list)
  names(s_list) = paste("s_", as.character(s_vals), sep = "")
  methods_list = lapply(1:length(methods), \(i) s_list)
  names(methods_list) = methods
  # populate methods_list
  for (i in 1:length(methods_list)){
    for (j in 1:length(s_list)){
      for (k in 1:length(alpha_list)){
        for (l in 1:length(sig_str_list)){
          for (m in 1:num_metrics){
            for (n in 1:num_trials){
              methods_list[[i]][[j]][[k]][[l]][[m]][n] =
                raw_list[[j]][[k]][[l]][[n]][[i]][m]
            }
          }
        }
      }
    }
  }
  return(methods_list)
}

#################################### HELPER FUNCTIONS ####################################

get_chol_mat = function(p, alpha){
  sigma = matrix(alpha, p, p)
  diag(sigma) = 1
  return(t(chol(sigma)))
}

################################### STRING TO FUNCTION ###################################

source("Lasso.R"); source("Coxnet.R"); source("Coxboost.R"); source("RandomForest.R");
source("BHAdjPValue.R"); source("QValue.R"); source("Cars.R")

names_to_methods = list(
  "LASSO" = LASSO, "ALASSO" = ALASSO, "CB" = CB, "RF_md" = RF_md, "RF_md_scr" = RF_md_scr,
  "BH_Proc" = BH_Proc, "QV" = QV, "CARS_med" = CARS, "CARS_msr" = CARS_residuals
)
