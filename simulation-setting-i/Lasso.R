########################################## LASSO #########################################

run_Lasso = function(data, n, p, s){
  
  library("dplyr")
  library("glmnet")
  library("survival")
  library("WGCNA")
  library("yardstick")
  
  # prep data
  train_bound = (2*n/3)
  X = as.matrix(data$df[1:train_bound, 1:p])
  Y = data$df$surv_time[1:train_bound]
  test_data = data$df[(train_bound+1):n,]
  test_true_tte = data$true_tte[(train_bound+1):n]
  
  # start computation timer
  start_time = Sys.time()
  
  # perform LASSO and extract parameter estimates
  beta_hat =
    cv.glmnet(x = X, y = Y, family = "cox", type.measure = "C")$lambda.min %>%
    glmnet(x = X, y = Y, family = "cox", lambda = .) %>%
    coef(.) %>%
    as.numeric(.)
  
  # stop computation timer
  computation_time = difftime(Sys.time(), start_time, units = "secs")
  
  # predict survival times using parameter estimates
  test_data$preds = as.vector(
    exp(-(as.matrix(test_data[,1:p]) %*% as.matrix(beta_hat)))
  )
  
  # get survival times using median estimate from proportional hazards model
  train_data = cbind(as.data.frame(X[,which(beta_hat != 0)]), data.frame(surv_time = Y))
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

LASSO = function(data, n, p, s){
  run_Lasso(data, n, p, s)
}
