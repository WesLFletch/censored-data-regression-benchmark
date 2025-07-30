######################################## CARS RDA ########################################

cars_rda = function(full_df, fold_idxs, prelim_mrna, residuals_elbow){
  
  library("carSurv")
  
  ########## FEATURE SELECTION PORTION ##########
  
  # filter mRNA features using preliminary results
  full_df = full_df[,colnames(full_df) %in% c(colnames(full_df)[1:4], prelim_mrna)]
  
  # perform feature selection
  # sig_preds =
  #   (carSurvScore(full_df$surv_time[,1], full_df$surv_time[,2], full_df[,-(1:4)]) %>%
  #      carVarSelect(., method = "threshold", threshold = 0.01) + 4) %>%
  #   colnames(full_df)[.]
  cars_scores = carSurvScore(full_df$surv_time[,1], full_df$surv_time[,2], full_df[,-(1:4)])
  # sig_preds = (sort(which.maxn(abs(cars_scores), floor(0.01*length(prelim_mrna)))) + 4) %>%
  #   colnames(full_df)[.]
  cars_curve = sort(abs(cars_scores), decreasing = TRUE)
  find_elbow = ifelse(residuals_elbow, # branch on residuals_elbow
    function(z){
      get_cum_err = function(z, split){
        df1 = data.frame(y = z[1:split], x = 1:split)
        df2 = data.frame(y = z[(split+1):length(z)], x = (split+1):length(z))
        lm1 = lm(y ~ x, data = df1); lm2 = lm(y ~ x, data = df2)
        return(sum(c(lm1$residuals^6, lm2$residuals^6)))
      }
      cum_errs = c(Inf, sapply(2:(length(z)-2), \(i)get_cum_err(z, i)))
      return(which.min(cum_errs))
    },
    function(z){
      n = length(z); p1 = c(1, z[1]); p2 = c(n, z[n])
      # Distance from each point to the line p1-p2
      distances = sapply(1:n, \(i){
        p = c(i, z[i])
        cross = abs((p2[1] - p1[1]) * (p1[2] - p[2]) - (p1[1] - p[1]) * (p2[2] - p1[2]))
        norm = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
        return(cross / norm)
      })
      return(which.max(distances))
    }
  )
  threshold = cars_curve[find_elbow(cars_curve)]
  sig_preds = colnames(full_df)[which(abs(cars_scores) >= threshold) + 4]
  
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
    # filter_cols = c(
    #   1:4,
    #   carSurvScore(full_df$surv_time[,1], full_df$surv_time[,2], full_df[,-(1:4)]) %>%
    #     carVarSelect(., method = "threshold", threshold = 0.01) + 4
    # )
    cars_scores = carSurvScore(full_df$surv_time[,1], full_df$surv_time[,2], full_df[,-(1:4)])
    # filter_cols = c(1:4, sort(which.maxn(abs(cars_scores), floor(0.01*length(prelim_mrna)))) + 4)
    cars_curve = sort(abs(cars_scores), decreasing = TRUE)
    threshold = cars_curve[find_elbow(cars_curve)]
    filter_cols = c(1:4, which(abs(cars_scores) >= threshold) + 4)
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

CARS = function(full_df, fold_idxs, prelim_mrna){
  cars_rda(full_df, fold_idxs, prelim_mrna, residuals_elbow = FALSE)
}

CARS_residuals = function(full_df, fold_idxs, prelim_mrna, residuals_elbow){
  cars_rda(full_df, fold_idxs, prelim_mrna, residuals_elbow = TRUE)
}
