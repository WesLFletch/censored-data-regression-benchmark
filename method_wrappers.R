########################################## LASSO #########################################

# fLASSO(X, Y, cens)
# Wrapper to perform LASSO for estimating log-hazards ratios.
#
# INPUTS:
#   X: `matrix` or `data.frame` object of data features with dimensions n by p.
#   Y: length n vector of times-to-event.
#   cens: length n binary vector, 0:censored, 1:uncensored.
#
# OUTPUT:
#   length p vector of estimated log-hazards ratios.
library("glmnet")
library("survival")
fLASSO = function(X, Y, cens){
  as.numeric(
    coef(
      cv.glmnet(as.matrix(X), Surv(Y, cens), family = "cox", type.measure = "C"),
      s = "lambda.min"
    )
  )
}

######################################### ALASSO #########################################

# fALASSO(X, Y, cens)
# Wrapper to perform ALASSO for estimating log-hazards ratios.
#
# INPUTS:
#   X: `matrix` or `data.frame` object of data features with dimensions n by p.
#   Y: length n vector of times-to-event.
#   cens: length n binary vector, 0:censored, 1:uncensored.
#
# OUTPUT:
#   length p vector of estimated log-hazards ratios.
library("Coxnet")
library("survival")
fALASSO = function(X, Y, cens){
  Coxnet(as.matrix(X), Surv(Y, cens), penalty = "Lasso", adaptive = TRUE, nfolds = 10)$Beta
}

######################################## CoxBoost ########################################

# fCoxBoost(X, Y, cens)
# Wrapper to perform CoxBoost for estimating log-hazards ratios.
#
# INPUTS:
#   X: `matrix` or `data.frame` object of data features with dimensions n by p.
#   Y: length n vector of times-to-event.
#   cens: length n binary vector, 0:censored, 1:uncensored.
#
# OUTPUT:
#   length p vector of estimated log-hazards ratios.
library("CoxBoost")
library("survival")
fCoxBoost = function(X, Y, cens){
  coef(CoxBoost(Y, cens, as.matrix(X)))
}

########################################### RSF ##########################################

# fRSF(X, Y, cens)
# Wrapper to perform RSF for creating a survival forest.
#
# INPUTS:
#   X: `matrix` or `data.frame` object of data features with dimensions n by p.
#   Y: length n vector of times-to-event.
#   cens: length n binary vector, 0:censored, 1:uncensored.
#
# OUTPUT:
#   `rfsrc` object fit using all p features.
library("randomForestSRC")
library("survival")
fRSF = function(X, Y, cens){
  full_df = data.frame(surv_time = Surv(Y, cens), X)
  rfsrc(surv_time ~ ., full_df, nodesize = 3)
}

########################################## sRSF ##########################################

# fsRSF(X, Y, cens)
# Wrapper to perform sRSF for 1) reducing dimensionality of X using the univariate Cox
# regression and 2) creating a survival forest using the remaining data features.
#
# INPUTS:
#   X: `matrix` or `data.frame` object of data features with dimensions n by p.
#   Y: length n vector of times-to-event.
#   cens: length n binary vector, 0:censored, 1:uncensored.
#
# OUTPUT:
#   `rfsrc` object fit using a subset of the original p features.
library("randomForestSRC")
library("survival")
fsRSF = function(X, Y, cens){
  pvals = apply(X, 2, \(col)summary(coxph(Surv(Y, cens) ~ col))$coefficients[5])
  X_red = X[,pvals < 0.1]
  full_df = data.frame(surv_time = Surv(Y, cens), X_red)
  rfsrc(surv_time ~ ., full_df, nodesize = 3)
}

###################################### BH procedure ######################################

# fBH_procedure(X, Y, cens)
# Wrapper to perform the BH procedure for 1) selecting features with adjusted p-values
# from the univariate Cox regression less than 0.05 and 2) returning estimated log-hazards
# ratios from the multivariate Cox regression using the selected features.
#
# INPUTS:
#   X: `matrix` or `data.frame` object of data features with dimensions n by p.
#   Y: length n vector of times-to-event.
#   cens: length n binary vector, 0:censored, 1:uncensored.
#
# OUTPUT:
#   length p vector of estimated log-hazards ratios.
library("mutoss")
library("survival")
fBH_procedure = function(X, Y, cens){
  X = as.data.frame(X)
  pvals = apply(X, 2, \(col)summary(coxph(Surv(Y, cens) ~ col))$coefficients[5])
  BH_rejected = BH(pvals, 0.05, TRUE)$rejected
  X_red = X[,BH_rejected]
  full_df = data.frame(surv_time = Surv(Y, cens), X_red)
  coefs = coef(coxph(surv_time ~ ., full_df))
  beta_hat = rep(0, ncol(X))
  names(beta_hat) = colnames(X)
  beta_hat[names(coefs)] = coefs
  return(beta_hat)
}

######################################### Q-value ########################################

# fQ_value(X, Y, cens)
# Wrapper to perform the Q-value procedure for 1) selecting features with q-values
# from the univariate Cox regression less than 0.05 and 2) returning estimated log-hazards
# ratios from the multivariate Cox regression using the selected features.
#
# INPUTS:
#   X: `matrix` or `data.frame` object of data features with dimensions n by p.
#   Y: length n vector of times-to-event.
#   cens: length n binary vector, 0:censored, 1:uncensored.
#
# OUTPUT:
#   length p vector of estimated log-hazards ratios.
library("survival")
library("WGCNA")
fQ_value = function(X, Y, cens){
  X = as.data.frame(X)
  pvals = apply(X, 2, \(col)summary(coxph(Surv(Y, cens) ~ col))$coefficients[5])
  QV_rejected = qvalue(pvals, fdr.level = 0.05)$significant
  X_red = X[,QV_rejected]
  full_df = data.frame(surv_time = Surv(Y, cens), X_red)
  coefs = coef(coxph(surv_time ~ ., full_df))
  beta_hat = rep(0, ncol(X))
  names(beta_hat) = colnames(X)
  beta_hat[names(coefs)] = coefs
  return(beta_hat)
}

####################################### CARS (MED) #######################################

# fCARS_MED(X, Y, cens)
# Wrapper to perform CARS with MED for 1) selecting features with large absolute CARS
# scores with threshold chosen by MED and 2) returning estimated log-hazards ratios from
# the multivariate Cox regression using the selected features.
#
# INPUTS:
#   X: `matrix` or `data.frame` object of data features with dimensions n by p.
#   Y: length n vector of times-to-event.
#   cens: length n binary vector, 0:censored, 1:uncensored.
#
# OUTPUT:
#   length p vector of estimated log-hazards ratios.
library("carSurv")
library("survival")
fCARS_MED = function(X, Y, cens){
  MED = function(z){
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
  X = as.data.frame(X)
  cars_scores = carSurvScore(Y, cens, as.matrix(X))
  sorted_abs_cars_scores = sort(abs(cars_scores), decreasing = TRUE)
  elbow_idx = MED(sorted_abs_cars_scores)
  threshold = sorted_abs_cars_scores[elbow_idx]
  selected = abs(cars_scores) >= threshold
  X_red = X[,selected]
  full_df = data.frame(surv_time = Surv(Y, cens), X_red)
  coefs = coef(coxph(surv_time ~ ., full_df))
  beta_hat = rep(0, ncol(X))
  names(beta_hat) = colnames(X)
  beta_hat[names(coefs)] = coefs
  return(beta_hat)
}

####################################### CARS (MSR) #######################################

# fCARS_MSR(X, Y, cens)
# Wrapper to perform CARS with MSR for 1) selecting features with large absolute CARS
# scores with threshold chosen by MSR and 2) returning estimated log-hazards ratios from
# the multivariate Cox regression using the selected features.
#
# INPUTS:
#   X: `matrix` or `data.frame` object of data features with dimensions n by p.
#   Y: length n vector of times-to-event.
#   cens: length n binary vector, 0:censored, 1:uncensored.
#
# OUTPUT:
#   length p vector of estimated log-hazards ratios.
library("carSurv")
library("survival")
fCARS_MSR = function(X, Y, cens){
  MSR = function(z){
    get_cum_err = function(z, split){
      df1 = data.frame(y = z[1:split], x = 1:split)
      df2 = data.frame(y = z[(split+1):length(z)], x = (split+1):length(z))
      lm1 = lm(y ~ x, data = df1); lm2 = lm(y ~ x, data = df2)
      return(sum(c(lm1$residuals^6, lm2$residuals^6)))
    }
    cum_errs = c(Inf, sapply(2:(length(z)-2), \(i)get_cum_err(z, i)))
    return(which.min(cum_errs))
  }
  X = as.data.frame(X)
  cars_scores = carSurvScore(Y, cens, as.matrix(X))
  sorted_abs_cars_scores = sort(abs(cars_scores), decreasing = TRUE)
  elbow_idx = MSR(sorted_abs_cars_scores)
  threshold = sorted_abs_cars_scores[elbow_idx]
  selected = abs(cars_scores) >= threshold
  X_red = X[,selected]
  full_df = data.frame(surv_time = Surv(Y, cens), X_red)
  coefs = coef(coxph(surv_time ~ ., full_df))
  beta_hat = rep(0, ncol(X))
  names(beta_hat) = colnames(X)
  beta_hat[names(coefs)] = coefs
  return(beta_hat)
}
