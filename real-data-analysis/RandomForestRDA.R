#################################### RANDOM FOREST RDA ###################################

randomforest_rda = function(full_df, fold_idxs, prelim_mrna, use_screen = FALSE){
  
  library("randomForestSRC")
  
  ########## FEATURE SELECTION PORTION ##########
  
  # filter mRNA features using preliminary results
  full_df = full_df[,colnames(full_df) %in% c(colnames(full_df)[1:4], prelim_mrna)]
  
  # perform screening if desired
  get_screen_cols = function(full_df){
    c(rep(TRUE, 4),
      (sapply(5:ncol(full_df), \(i){
        summary(
          coxph(
            surv_time ~ full_df[,i], data = full_df, robust = TRUE
          )
        )$coefficients[6]
      }) %>% replace_na(., 1)) < 0.1) %>% which(.)
  }
  if (use_screen) fs_df = full_df[,get_screen_cols(full_df)] else fs_df = full_df
  # create survival forest
  fs_rfsrc = rfsrc(surv_time ~ ., data = fs_df[,-c(1,3,4)], nodesize = 3)
  # extract names of significant predictors
  sig_preds = var.select(fs_rfsrc, conservative = "medium", verbose = FALSE)$topvars
  
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
    # perform screening if desired and create survival forest, measuring computation time
    start_time = Sys.time()
    if (use_screen){
      screen_cols = get_screen_cols(train_df)
      train_df = train_df[,screen_cols]
      test_df = test_df[,screen_cols]
    }
    pred_rfsrc = rfsrc(surv_time ~ ., data = train_df[,-1], nodesize = 3)
    comp_time = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    # get predicted survival times (using them is equivalent to risk scores)
    test_df$preds = as.vector(predict(pred_rfsrc, newdata = test_df[,-(1:2)])$predicted)
    # compute out-of-sample CI
    CI = concordance_survival(test_df, surv_time, preds)$.estimate
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

RF_md = function(full_df, fold_idxs, prelim_mrna){
  randomforest_rda(full_df, fold_idxs, prelim_mrna, use_screen = FALSE)
}

RF_md_screen = function(full_df, fold_idxs, prelim_mrna){
  randomforest_rda(full_df, fold_idxs, prelim_mrna, use_screen = TRUE)
}
