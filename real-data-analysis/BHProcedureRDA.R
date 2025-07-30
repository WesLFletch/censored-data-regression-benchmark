#################################### BH PROCEDURE RDA ####################################

bhprocedure_rda = function(full_df, fold_idxs, prelim_mrna){
  
  library("mutoss")
  
  ########## FEATURE SELECTION PORTION ##########
  
  # filter mRNA features using preliminary results
  full_df = full_df[,colnames(full_df) %in% c(colnames(full_df)[1:4], prelim_mrna)]
  
  # perform feature selection
  sig_preds = c(rep(FALSE, 4), (sapply(5:ncol(full_df), function(i){
    summary(coxph(surv_time ~ full_df[,i], data = full_df))$coefficients[5]
  }) %>% replace_na(., 1) %>% BH(., alpha = 0.05, silent = TRUE))[[3]]) %>%
    colnames(full_df)[.]
  
  ########## PREDICTIVE ABILITY PORTION ##########
  
  # shuffle observations, so that folds are random
  set.seed(1)
  full_df = full_df[sample(1:nrow(full_df), nrow(full_df)),]
  # runs a single fold, measuring CI and computation time
  run_fold = function(full_df, test_idxs){
    # split into train and test sets
    train_df = full_df[-test_idxs,]
    test_df = full_df[test_idxs,]
    # remove predictors that are constant in train fold
    nonconst_preds = sapply(1:ncol(train_df), \(i)length(unique(train_df[,i])) > 1)
    train_df = train_df[,nonconst_preds]
    test_df = test_df[,nonconst_preds]
    # perform feature selection and estimate beta vector, measuring computation time
    start_time = Sys.time()
    filter_cols = c(rep(TRUE, 4), (sapply(5:ncol(train_df), function(i){
      summary(coxph(surv_time ~ train_df[,i], data = train_df))$coefficients[5]
    }) %>% replace_na(., 1) %>% BH(., alpha = 0.05, silent = TRUE))[[3]])
    train_df = train_df[,filter_cols]
    test_df = test_df[,filter_cols]
    beta_hat = coef(coxph(surv_time ~ ., data = train_df[,-1], x = TRUE)) %>%
      replace_na(., 0)
    comp_time = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    # get risk scores
    test_df$risk_score = as.vector(exp(-(as.matrix(test_df[,-(1:2)]) %*% beta_hat)))
    # compute out-of-sample CI
    CI = concordance_survival(test_df, surv_time, risk_score)$.estimate
    # return fold results
    return(c(CI, comp_time))
  }
  # perform 10 fold evaluations
  fold_results = sapply(fold_idxs, \(test_idxs)run_fold(full_df, test_idxs))
  
  ########## RETURN RESULTS ##########
  
  return(list(
    CI = fold_results[1,],
    computation_time = fold_results[2,],
    selected_features = sig_preds
  ))
  
}

BHP = function(full_df, fold_idxs, prelim_mrna){
  bhprocedure_rda(full_df, fold_idxs, prelim_mrna)
}
