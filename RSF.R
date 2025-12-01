######################################### RUN RSF ########################################

priv_RSF = function(X_train, Y_train, screen){
  best_nodesize = 8
  best_mtry = 403
  features_char = colnames(X_train)
  screen_cols = if (screen){
    sapply(1:ncol(X_train), \(i){
      replaceMissing(summary(coxph(Y_train~X_train[,i]))$coefficients[5], 1)<0.1
    })
  } else rep(TRUE, ncol(X_train))
  train_df = data.frame(X_train[,screen_cols], time = Y_train[,1], status = Y_train[,2])
  colnames(train_df)[1:sum(screen_cols)] = colnames(X_train)[screen_cols]
  my_rsf = rfsrc(Surv(time, status) ~ ., train_df, mtry = best_mtry, nodesize = best_nodesize)
  # prepare outputs
  predict = function(X_test){
    p_obj = predict.rfsrc(my_rsf, newdata = X_test[,screen_cols])
    step_times = p_obj$time.interest
    sf_mat = p_obj$survival
    get_median = function(k)step_times[which(sf_mat[k,]<0.5)[1]]
    sapply(1:nrow(X_test), get_median)
  }
  inv_risk = predict
  bscore = function(X_test, Y_test, time_points){
    test_df = data.frame(X_test[,screen_cols], time = Y_test[,1], status = Y_test[,2])
    colnames(test_df)[1:sum(screen_cols)] = colnames(X_test)[screen_cols]
    my_pec = pec(my_rsf, Surv(time, status) ~ 1, data = test_df, times = time_points,
                 exact = TRUE, cens.model = "marginal")
    sapply(time_points, \(t)my_pec$AppErr$rfsrc[which.min(abs(my_pec$time-t))])
    
  }
  selected_char = var.select(object = my_rsf, verbose = FALSE)$topvars
  selected = features_char %in% selected_char
  list(predict = predict, inv_risk = inv_risk, bscore = bscore, selected = selected)
}

RSF = function(X_train, Y_train)priv_RSF(X_train, Y_train, screen = FALSE)
sRSF = function(X_train, Y_train)priv_RSF(X_train, Y_train, screen = TRUE)
