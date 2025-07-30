#################################### RANDOM FOREST RDS ###################################

randomforest_rds = function(new_df, true_tte, true_sig_preds, use_screen = FALSE){
  
  library("randomForestSRC")
  
  # split into train and test dataframes
  n_train = floor(9*nrow(new_df)/10)
  train_df = new_df[1:n_train,]
  test_df = new_df[-(1:n_train),]
  test_true_tte = true_tte[-(1:n_train)]
  
  # create survival forest, measure computation time
  start_time = Sys.time()
  get_screen_cols = function(new_df){
    c(TRUE, (sapply(2:ncol(train_df), \(i){summary(coxph(
      surv_time ~ new_df[,i], data = new_df
    ))$coefficients[5]}) %>% replace_na(., 1)) < 0.1) %>% which(.)
  }
  if (use_screen) screen_cols = get_screen_cols(new_df) else screen_cols = 1:ncol(new_df)
  train_df = train_df[,screen_cols]; test_df = test_df[,screen_cols]
  out_rfsrc = rfsrc(surv_time ~ ., data = train_df, nodesize = 3)
  sig_preds = var.select(out_rfsrc, verbose = FALSE)$topvars
  comp_time = difftime(Sys.time(), start_time, units = "secs")
  
  # estimate survival times on test dataframe
  test_df$surv_preds = predict(out_rfsrc, newdata = test_df[,-1])$predicted
  
  # get metric values
  CI = concordance_survival(test_df, surv_time, surv_preds)$.estimate
  RMSE = sqrt(mean(log(test_df$surv_preds/test_true_tte)^2, na.rm = TRUE))
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

RF_md = function(new_df, true_tte, true_sig_preds){
  randomforest_rds(new_df, true_tte, true_sig_preds, use_screen = FALSE)
}

RF_md_screen = function(new_df, true_tte, true_sig_preds){
  randomforest_rds(new_df, true_tte, true_sig_preds, use_screen = TRUE)
}
