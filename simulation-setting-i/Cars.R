####################################### CARS FILTER ######################################

run_Cars = function(data, n, p, s, residuals_elbow){
  
  library("carSurv")
  library("doBy")
  library("dplyr")
  library("survival")
  library("WGCNA")
  library("yardstick")
  
  # prep data
  train_bound = (2*n/3)
  train_data = data$df[1:train_bound,]
  test_data = data$df[(train_bound+1):n,]
  test_true_tte = data$true_tte[(train_bound+1):n]
  
  # create beta_hat (have to make it before removing predictors)
  beta_hat = rep(0, p)
  names(beta_hat) = names(train_data)[1:p]
  
  # start computation timer
  start_time = Sys.time()
  
  # get CARS scores and perform variable selection
  cars_scores = carSurvScore(train_data$surv_time[,1], train_data$surv_time[,2],
                             train_data[,1:p])
  cars_curve = cars_curve = sort(abs(cars_scores), decreasing = TRUE)
  find_elbow = ifelse(residuals_elbow, # branch on residuals_elbow
    function(z){ # this is estimation by MSR
      get_cum_err = function(z, split){
        df1 = data.frame(y = z[1:split], x = 1:split)
        df2 = data.frame(y = z[(split+1):length(z)], x = (split+1):length(z))
        lm1 = lm(y ~ x, data = df1); lm2 = lm(y ~ x, data = df2)
        return(sum(c(lm1$residuals^6, lm2$residuals^6)))
      }
      cum_errs = c(Inf, sapply(2:(length(z)-2), \(i)get_cum_err(z, i)))
      return(which.min(cum_errs))
    },
    function(z){ # this is estimation by MED
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
  sig_preds = which(abs(cars_scores) >= threshold)
  
  # trim down train_data
  train_data = train_data[,c(sig_preds,p+1)]
  
  # fit PH model to trimmed data
  coxph_obj = coxph(surv_time ~ ., data = train_data, x = TRUE)
  beta_hat[names(coef(coxph_obj))] = coef(coxph_obj)
  
  # stop computation timer
  computation_time = difftime(Sys.time(), start_time, units = "secs")
  
  # predict risk scores using parameter estimates
  test_data$preds = as.vector(
    exp(-(as.matrix(test_data[,1:p]) %*% as.matrix(beta_hat)))
  )
  
  # get survival times using median estimate from proportional hazards model
  get_median = function(k){
    surv_fit = survfit(coxph_obj, newdata = test_data[k, which(beta_hat != 0)])
    summary(surv_fit)$table["median"]
  }
  median_preds = unname(sapply(1:nrow(test_data), get_median))
  
  # get number of TP, TN, FP, and FN
  TP = sum(beta_hat[1:s] != 0)
  FN = s - TP
  FP = sum(beta_hat[(s+1):p] != 0)
  TN = p - s - FP
  
  # get metrics
  CI = concordance_survival(test_data, surv_time, preds)$.estimate
  precision = replaceMissing(TP / (TP + FP)) # higher is better
  recall = TP / (TP + FN) # = power, higher is better
  FDR = 1 - precision # lower is better
  f1_score = replaceMissing(2/((1/recall)+(1/precision)))
  rmse1 = sqrt(mean(log(median_preds/test_true_tte)^2, na.rm = TRUE))
  
  return(c(CI, f1_score, FDR, computation_time, rmse1))
  
}

####################################### VARIATIONS #######################################

CARS = function(data, n, p, s){ # synonymous with CARS with MED
  run_Cars(data, n, p, s, residuals_elbow = FALSE)
}

CARS_residuals = function(data, n, p, s){ # synonymous with CARS with MSR
  run_Cars(data, n, p, s, residuals_elbow = TRUE)
}
