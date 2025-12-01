###################################### RUN BH OR QV ######################################

priv_BHQV = function(X_train, Y_train, FDR_adj_char){
  features = colnames(X_train)
  FDR_adj = function(pvals) if (FDR_adj_char == "BH"){
    BH(pvals, 0.05, silent = TRUE)$rejected
  } else qvalue(pvals, pi0.method = "smoother", fdr.level = 0.05)$significant
  screen_cols = apply(X_train, 2, \(col)summary(coxph(Y_train ~ col))$coefficients[5]) %>%
    replaceMissing(., 1) %>%
    FDR_adj(.)
  # prepare outputs
  predict = function(X_test){
    if (any(screen_cols)){
      train_df = data.frame(X_train[,screen_cols], surv_time = Y_train)
      colnames(train_df)[-ncol(train_df)] = features[screen_cols]
      model = coxph(surv_time ~ ., train_df)
      get_median = function(k){
        test_df = data.frame(X_test[k, screen_cols])
        colnames(test_df) = features[screen_cols]
        surv_fit = survfit(model, newdata = test_df)
        summary(surv_fit)$table["median"]
      }
      unname(sapply(1:nrow(X_test), get_median))
    } else unname(rep(summary(survfit(coxph(Y_train ~ 1)))$table["median"], nrow(X_test)))
  }
  inv_risk = function(X_test){
    if (any(screen_cols)){
      train_df = data.frame(X_train[,screen_cols], surv_time = Y_train)
      colnames(train_df)[-ncol(train_df)] = features[screen_cols]
      beta_hat = replaceMissing(coef(coxph(surv_time ~ ., train_df)), 0)
      as.vector(exp(-(as.matrix(X_test[,screen_cols]) %*% beta_hat)))
    } else rep(1, nrow(X_test)) # with all beta_hat=0, exp(-(X %*% beta_hat))=1 for any X
  }
  bscore = function(X_test, Y_test, time_points){
    if (any(screen_cols)){
      train_df = data.frame(X_train[,screen_cols], surv_time = Y_train)
      colnames(train_df)[-ncol(train_df)] = features[screen_cols]
      test_df = data.frame(X_test[,screen_cols], surv_time = Y_test)
      colnames(test_df)[-ncol(test_df)] = features[screen_cols]
      model = coxph(surv_time ~ ., data = train_df)
      brier(model, time_points, newdata = test_df)$brier
    }else{
      train_df = data.frame(surv_time = Y_train, dummy = 1)
      test_df = data.frame(surv_time = Y_test, dummy = 1)
      model = coxph(surv_time ~ dummy, train_df)
      brier(model, time_points, newdata = test_df)$brier
    }
  }
  selected = screen_cols
  list(predict = predict, inv_risk = inv_risk, bscore = bscore, selected = selected)
}

BHP = function(X_train, Y_train)priv_BHQV(X_train, Y_train, FDR_adj_char = "BH")
QV = function(X_train, Y_train)priv_BHQV(X_train, Y_train, FDR_adj_char = "QV")
