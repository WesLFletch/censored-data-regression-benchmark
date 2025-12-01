######################################## RUN LASSO #######################################

priv_LASSO = function(X_train, Y_train, a){
  X = as.matrix(X_train); Y = Y_train
  beta_hat =
    cv.glmnet(x = X, y = Y, family = "cox", alpha = a, type.measure = "C")$lambda.min %>%
    glmnet(x = X, y = Y, family = "cox", alpha = a, lambda = .) %>%
    coef(.) %>%
    as.numeric(.)
  # prepare outputs
  predict = function(X_test){
    train_df = data.frame(X[,beta_hat!=0], surv_time = Y)
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
    train_df = data.frame(X[,beta_hat!=0], surv_time = Y)
    test_df = data.frame(X_test[,beta_hat!=0], surv_time = Y_test)
    model = coxph(surv_time ~ ., data = train_df)
    brier(model, time_points, newdata = test_df)$brier
  }
  selected = beta_hat != 0
  list(predict = predict, inv_risk = inv_risk, bscore = bscore, selected = selected)
}

LASSO = function(X_train, Y_train)priv_LASSO(X_train, Y_train, a = 1)
ENET = function(X_train, Y_train)priv_LASSO(X_train, Y_train, a = 0.5)
