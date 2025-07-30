######################################### Q VALUE ########################################

run_QValue = function(data, n, p, s){
  
  library("doBy")
  library("dplyr")
  library("survival")
  library("WGCNA")
  library("yardstick")
  
  # prep data
  train_bound = (2*n/3)
  train_data = data$df[1:train_bound,]
  test_data = data$df[(train_bound+1):n,]
  test_true_tte = data$true_tte[(train_bound+1):n]
  
  # create beta_hat (have to make it before removing predictors)
  beta_hat = rep(0, p)
  names(beta_hat) = names(train_data)[1:p]
  
  # start computation timer
  start_time = Sys.time()
  
  # vectorize getting q_vals, if possible
  q_value_success = FALSE
  tryCatch({
    q_vals =
      (sapply(1:p, function(r){
        summary(coxph(surv_time ~ train_data[,r], data = train_data,
                      robust = TRUE))$coefficients[6]
      }) %>%
        qvalue(.))$qvalues
    q_value_success = TRUE
  }, error = function(msg){q_value_success = FALSE})
  
  if (q_value_success){
    # trim train_data
    train_data =
      intersect(
        which(q_vals < 0.05),
        which.minn(q_vals, train_bound - 2)
      ) %>%
      union(., which.min(q_vals)) %>% # FIXME: THIS IS THE CULPRIT
      c(., p+1) %>%
      sort(.) %>%
      train_data[,.]
    
    # fit coxph model and extract parameter estimates
    coxph_obj = coxph(surv_time ~ ., data = train_data, x = TRUE)
    beta_hat[names(coef(coxph_obj))] = coef(coxph_obj)
  }else{
    coxph_obj = coxph(surv_time ~ 1, data = train_data)
  }
  
  # stop computation timer
  computation_time = difftime(Sys.time(), start_time, units = "secs")
  
  # predict risk scores using parameter estimates
  test_data$preds = as.vector(
    exp(-(as.matrix(test_data[,1:p]) %*% as.matrix(beta_hat)))
  )
  
  # get survival times using median estimate from proportional hazards model
  get_median = function(k){
    surv_fit = survfit(coxph_obj, newdata = test_data[k,])
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

QV = function(data, n, p, s){
  run_QValue(data, n, p, s)
}
