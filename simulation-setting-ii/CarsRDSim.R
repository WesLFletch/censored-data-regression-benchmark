######################################## CARS RDS ########################################

cars_rds = function(new_df, true_tte, true_sig_preds, residuals_elbow){
  
  library("carSurv")
  
  # split into train and test dataframes
  n_train = floor(9*nrow(new_df)/10)
  train_df = new_df[1:n_train,]
  test_df = new_df[-(1:n_train),]
  test_true_tte = true_tte[-(1:n_train)]
  
  # perform feature selection, estimate beta vector, measure computation time
  start_time = Sys.time()
  beta_hat = rep(0, ncol(train_df)-1)
  names(beta_hat) = colnames(train_df)[-1]
  cars_scores = carSurvScore(train_df$surv_time[,1], train_df$surv_time[,2],
                             as.matrix(train_df[,-1]))
  cars_curve = cars_curve = sort(abs(cars_scores), decreasing = TRUE)
  find_elbow = ifelse(residuals_elbow, # use either MED or MSR to find elbow
    function(z){ # the MSR algorithm
      get_cum_err = function(z, split){
        df1 = data.frame(y = z[1:split], x = 1:split)
        df2 = data.frame(y = z[(split+1):length(z)], x = (split+1):length(z))
        lm1 = lm(y ~ x, data = df1); lm2 = lm(y ~ x, data = df2)
        return(sum(c(lm1$residuals^6, lm2$residuals^6)))
      }
      cum_errs = c(Inf, sapply(2:(length(z)-2), \(i)get_cum_err(z, i)))
      return(which.min(cum_errs))
    },
    function(z){ # the MED algorithm
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
  sig_preds = colnames(train_df)[which(abs(cars_scores) >= threshold) + 1]
  ph_model = coxph(
    surv_time ~ .,
    data = train_df[,colnames(train_df) %in% c("surv_time", sig_preds)]
  )
  beta_hat[sig_preds] = coef(ph_model)
  comp_time = difftime(Sys.time(), start_time, units = "secs")
  
  # estimate median survival times on test dataframe
  get_median = function(i){
    surv_fit = survfit(ph_model, newdata = test_df[i, sig_preds])
    summary(surv_fit)$table["median"]
  }
  median_preds = unname(sapply(1:nrow(test_df), get_median))
  
  # get risk scores on test dataframe
  test_df$risk_scores = as.vector(exp(-(as.matrix(test_df[,-1]) %*% beta_hat)))
  
  # get metric values
  CI = concordance_survival(test_df, surv_time, risk_scores)$.estimate
  RMSE = sqrt(mean(log(median_preds/test_true_tte)^2, na.rm = TRUE))
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

CARS_med = function(new_df, true_tte, true_sig_preds){
  cars_rds(new_df, true_tte, true_sig_preds, residuals_elbow = FALSE)
}

CARS_msr = function(new_df, true_tte, true_sig_preds){
  cars_rds(new_df, true_tte, true_sig_preds, residuals_elbow = TRUE)
}
