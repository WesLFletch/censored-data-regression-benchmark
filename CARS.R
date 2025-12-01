######################################## RUN CARS ########################################

priv_CARS = function(X_train, Y_train, elbow_est_char){
  features = colnames(X_train)
  elbow_est = function(z){
    if (elbow_est_char == "MED"){
      n = length(z); p1 = c(1, z[1]); p2 = c(n, z[n])
      # Distance from each point to the line p1-p2
      distances = sapply(1:n, \(i){
        p = c(i, z[i])
        cross = abs((p2[1] - p1[1]) * (p1[2] - p[2]) - (p1[1] - p[1]) * (p2[2] - p1[2]))
        norm = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
        return(cross / norm)
      })
      return(which.max(distances))
    }else{
      get_cum_err = function(z, split){
        df1 = data.frame(y = z[1:split], x = 1:split)
        df2 = data.frame(y = z[(split+1):length(z)], x = (split+1):length(z))
        lm1 = lm(y ~ x, data = df1); lm2 = lm(y ~ x, data = df2)
        return(sum(c(lm1$residuals^6, lm2$residuals^6)))
      }
      cum_errs = c(Inf, sapply(2:(length(z)-2), \(i)get_cum_err(z, i)))
      return(which.min(cum_errs))
    }
  }
  nonconst_preds = apply(X_train, 2, \(col)length(unique(col))>1)
  X = X_train[,nonconst_preds]
  cars_scores_raw = carSurvScore(Y_train[,1], Y_train[,2], X)
  names(cars_scores_raw) = colnames(X)
  cars_scores = rep(0, ncol(X_train))
  names(cars_scores) = colnames(X_train)
  cars_scores[names(cars_scores_raw)] = cars_scores_raw
  cars_curve = sort(abs(cars_scores), decreasing = TRUE)
  threshold = cars_curve[elbow_est(cars_curve)]
  screen_cols = abs(cars_scores) >= threshold
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

CARS_MED = function(X_train, Y_train)priv_CARS(X_train, Y_train, elbow_est_char = "MED")
CARS_MSR = function(X_train, Y_train)priv_CARS(X_train, Y_train, elbow_est_char = "MSR")
