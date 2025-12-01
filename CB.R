###################################### RUN COXBOOST ######################################

priv_CB = function(X_train, Y_train){
  X = as.matrix(X_train); time = Y_train[,1]; status = Y_train[,2]
  beta_hat = coef(CoxBoost(time, status, X))
  # prepare outputs
  predict = function(X_test){
    train_df = data.frame(X[,beta_hat!=0], surv_time = Y_train)
    model = coxph(surv_time ~ ., data = train_df)
    get_median = function(k){
      surv_fit = survfit(model, newdata = data.frame(X_test[k, beta_hat!=0]))
      summary(surv_fit)$table["median"]
    }
    unname(sapply(1:nrow(X_test), get_median))
  }
  inv_risk = function(X_test){
    as.vector(exp(-(as.matrix(X_test) %*% beta_hat)))
  }
  bscore = function(X_test, Y_test, time_points){
    train_df = data.frame(X[,beta_hat!=0], surv_time = Y_train)
    test_df = data.frame(X_test[,beta_hat!=0], surv_time = Y_test)
    model = coxph(surv_time ~ ., data = train_df)
    brier(model, time_points, newdata = test_df)$brier
  }
  selected = beta_hat != 0
  list(predict = predict, inv_risk = inv_risk, bscore = bscore, selected = selected)
}

CB = priv_CB
