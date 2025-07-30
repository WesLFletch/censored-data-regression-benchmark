###################################### RANDOM FOREST #####################################

run_RandomForest = function(data, n, p, s, use_min_depth = FALSE, use_screen = FALSE){
  
  library("doBy")
  library("randomForestSRC")
  library("survival")
  library("WGCNA")
  library("yardstick")
  
  # prep data
  train_bound = (2*n/3)
  train_data = data$df[1:train_bound,]
  test_data = data$df[(train_bound+1):n,]
  test_true_tte = data$true_tte[(train_bound+1):n]
  
  # adjust number of features selected for VIMP feature selection (NO LONGER USED)
  prop_to_select = ifelse(s < 100, 0.01, 0.1)
  
  # start computation timer
  start_time = Sys.time()
  
  # branch on use_min_depth
  if (!use_min_depth){
    # perform Random Forest and extract significant predictors
    out_rfsrc = rfsrc(surv_time ~ ., data = train_data, nodesize = 3, importance = TRUE)
    sig_preds = which.maxn(out_rfsrc$importance, floor(prop_to_select*p))
  }else{
    if (use_screen){
      p_vals = sapply(1:p, \(i){
        summary(coxph(train_data$surv_time ~ train_data[,i]))$coefficients[5]
      })
      filter_cols = c(which(p_vals < 0.1), p+1)
      train_data = train_data[,filter_cols]
      test_data = test_data[,filter_cols]
    }
    out_rfsrc = rfsrc(surv_time ~ ., data = train_data, nodesize = 3, importance = TRUE)
    sig_preds = var.select(out_rfsrc, verbose = FALSE)$topvars
    sig_preds = as.numeric(gsub("V", "", sig_preds))
  }
  
  # stop computation timer
  computation_time = difftime(Sys.time(), start_time, units = "secs")
  
  # predict survival times using out_rfsrc
  test_data$preds = as.vector(
    predict(out_rfsrc, newdata = test_data[,-ncol(test_data)])$predicted
  )
  
  # get number of TP, TN, FP, and FN
  TP = sum(sig_preds <= s)
  FN = s - TP
  FP = sum(sig_preds > s)
  TN = p - s - FP
  
  # get metrics
  CI = concordance_survival(test_data, surv_time, preds)$.estimate
  precision = replaceMissing(TP / (TP + FP)) # higher is better
  recall = TP / (TP + FN) # = power, higher is better
  FDR = 1 - precision # lower is better
  f1_score = replaceMissing(2/((1/recall)+(1/precision)))
  rmse1 = sqrt(mean(log(test_data$preds/test_true_tte)^2, na.rm = TRUE))
  
  return(c(CI, f1_score, FDR, computation_time, rmse1))
  
}

####################################### VARIATIONS #######################################

RF = function(data, n, p, s){ # NOTE: This method was removed due to poor performance
  run_RandomForest(data, n, p, s, use_min_depth = FALSE, use_screen = FALSE)
}

RF_md = function(data, n, p, s){
  run_RandomForest(data, n, p, s, use_min_depth = TRUE, use_screen = FALSE)
}

RF_md_scr = function(data, n, p, s){
  run_RandomForest(data, n, p, s, use_min_depth = TRUE, use_screen = TRUE)
}
