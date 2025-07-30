######################################### COXNET #########################################

run_Coxnet = function(data, n, p, s){
  
  library("Coxnet")
  library("survival")
  library("WGCNA")
  library("yardstick")
  library("mutoss")
  
  # prep data
  train_bound = (2*n/3)
  y = data$df$surv_time[1:train_bound]
  x = as.matrix(data$df[1:train_bound, 1:p])
  test_data = data$df[(train_bound+1):n,]
  test_true_tte = data$true_tte[(train_bound+1):n]
  
  # start computation timer
  start_time = Sys.time()
  
  # perform ALASSO and extract parameter estimates
  beta_hat = Coxnet(x, y, penalty = "Lasso", adaptive = TRUE, nfolds = 10)$Beta0
  
  # stop computation timer
  computation_time = difftime(Sys.time(), start_time, units = "secs")
  
  # predict risk scores using parameter estimates
  test_data$preds = as.vector(
    exp(-(as.matrix(test_data[,1:p]) %*% as.matrix(beta_hat)))
  )
  
  # get survival times using median estimate from proportional hazards model
  train_data = cbind(as.data.frame(x[,which(beta_hat != 0)]), data.frame(surv_time = y))
  model = coxph(surv_time ~ ., data = train_data)
  get_median = function(k){
    surv_fit = survfit(model, newdata = test_data[k, which(beta_hat != 0)])
    summary(surv_fit)$table["median"]
  }
  median_preds = unname(sapply(1:nrow(test_data), get_median))
  
  # get number of TP, TN, FP, and FN
  TP = sum(beta_hat[1:s] != 0)
  FN = s - TP
  FP = sum(beta_hat[(s+1):p] != 0)
  TN = p - s - FP
  
  # get metrics
  CI = concordance_survival(test_data, surv_time, preds)$.estimate
  precision = replaceMissing(TP / (TP + FP)) # higher is better
  recall = TP / (TP + FN) # = power, higher is better
  FDR = 1 - precision # lower is better
  f1_score = replaceMissing(2/((1/recall)+(1/precision)))
  rmse1 = sqrt(mean(log(median_preds/test_true_tte)^2, na.rm = TRUE))
  
  return(c(CI, f1_score, FDR, computation_time, rmse1))
  
}

####################################### VARIATIONS #######################################

ALASSO = function(data, n, p, s){
  run_Coxnet(data, n, p, s)
}
