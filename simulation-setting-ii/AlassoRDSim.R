####################################### ALASSO RDS #######################################

alasso_rds = function(new_df, true_tte, true_sig_preds){
  
  library("Coxnet")
  
  # split into train and test dataframes
  n_train = floor(9*nrow(new_df)/10)
  train_df = new_df[1:n_train,]
  test_df = new_df[-(1:n_train),]
  test_true_tte = true_tte[-(1:n_train)]
  
  # estimate beta vector, measure computation time
  start_time = Sys.time()
  beta_hat = Coxnet(as.matrix(train_df[,-1]), train_df$surv_time,
                    penalty = "Lasso", adaptive = TRUE, nfolds = 10)$Beta0
  names(beta_hat) = colnames(train_df)[-1]
  sig_preds = names(beta_hat[beta_hat != 0])
  comp_time = difftime(Sys.time(), start_time, units = "secs")
  
  # estimate median survival times on test dataframe
  ph_model = coxph(
    surv_time ~ .,
    data = train_df[,colnames(train_df) %in% c("surv_time", sig_preds)]
  )
  get_median = function(i){
    surv_fit = survfit(ph_model, newdata = test_df[i, sig_preds])
    summary(surv_fit)$table["median"]
  }
  median_preds = unname(sapply(1:nrow(test_df), get_median))
  
  # get risk scores on test dataframe
  test_df$risk_scores = as.vector(exp(-(as.matrix(test_df[,-1]) %*% beta_hat)))
  
  # get metric values
  CI = concordance_survival(test_df, surv_time, risk_scores)$.estimate
  RMSE = sqrt(mean(log(median_preds/test_true_tte)^2, na.rm = TRUE))
  TP = length(intersect(sig_preds, true_sig_preds))
  FP = length(sig_preds) - TP
  FN = length(true_sig_preds) - TP
  F1 = TP/(TP+0.5*(FP+FN))
  FDR = replace_na(FP/(TP+FP), 0)
  
  # return metric values
  return(c(
    "CI" = CI,
    "RMSE" = RMSE,
    "F1-Score" = F1,
    "FDR" = FDR,
    "comp_time" = comp_time
  ))
  
}

ALASSO = function(new_df, true_tte, true_sig_preds){
  alasso_rds(new_df, true_tte, true_sig_preds)
}
