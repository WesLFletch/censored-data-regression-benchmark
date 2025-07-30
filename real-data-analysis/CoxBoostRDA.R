###################################### COXBOOST RDA ######################################

coxboost_rda = function(full_df, fold_idxs, prelim_mrna){
  
  library("CoxBoost")
  
  ########## FEATURE SELECTION PORTION ##########
  
  # filter mRNA features using preliminary results
  full_df = full_df[,colnames(full_df) %in% c(colnames(full_df)[1:4], prelim_mrna)]
  
  # estimate beta vector
  beta_hat = coef(CoxBoost(
    full_df$surv_time[,1], full_df$surv_time[,2], as.matrix(full_df[,-(1:4)])
  ))
  names(beta_hat) = colnames(full_df)[-(1:4)]
  # extract names of significant predictors
  sig_preds = names(which(beta_hat != 0))
  
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
    # estimate beta vector, measuring computation time
    start_time = Sys.time()
    beta_hat = coef(CoxBoost(
      train_df$surv_time[,1], train_df$surv_time[,2], as.matrix(train_df[,-(1:2)]),
      unpen.index = 1:2
    ))
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

CB = function(full_df, fold_idxs, prelim_mrna){
  coxboost_rda(full_df, fold_idxs, prelim_mrna)
}
