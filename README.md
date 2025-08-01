# Censored Data Regression Benchmark Materials

Supplementary code for performing the `LASSO`, `ALASSO`, `CoxBoost`, `RSF`, `sRSF`, `BH` and `Q-value` procedures, and `CARS` with `MED` and `MSR` methods as done in our benchmark study.

## User Guide

### Function Definitions

The file `method_wrappers.R` contains the R code required to perform all methods on a survival time dataset. All functions require the following arguments:

- `X`, a feature matrix with dimensions $n$ by $p$
- `Y`, a vector of times-to-event with length $n$
- `cens`, a binary $n$ dimensional vector, where a `0` and `1` denote a censored and an uncensored observation, respectively

Most methods are parametric and return log-hazard ratio estimates ($\widehat\beta$) for the proportional hazards model; The RSF/sRSF methods are nonparametric and return a survival forest (`rfsrc`) object.

### Usage Examples

```{R}
##### generate synthetic survival time dataset

library("MASS")
set.seed(1234)
# set data dimensions
n = 100
p = 300
# X follows N(0, Sigma)
Sigma = matrix(0.3, nrow = p, ncol = p)
diag(Sigma) = 1
X = mvrnorm(n, mu = rep(0, p), Sigma)
# first ten features are effect variables
beta = c(rep(1, 5), rep(-1, 5), rep(0, 290))
lin_pred = as.vector(X %*% beta)
# baseline hazard is exponential distribution
failure_time = rexp(n, exp(lin_pred))
# censorship is independent of features
cens_time = runif(n, min = quantile(failure_time, 0.2), max = quantile(failure_time, 0.9))
# get event time with censorship indicator
Y = pmin(failure_time, cens_time)
cens = as.numeric(failure_time <= cens_time) # 83% uncensored observations
# resulting dataset is (X, Y, cens)

##### perform methods on dataset

source("method_wrappers.R")

# get three beta estimates from parametric methods
beta_hat_LASSO = fLASSO(X, Y, cens)
beta_hat_Q_value = fQ_value(X, Y, cens)
beta_hat_CARS_MSR = fCARS_MSR(X, Y, cens)

# create survival forests
my_RSF = fRSF(X, Y, cens)
my_sRSF = fsRSF(X, Y, cens)

# determine features selected by all methods
selected_preds_LASSO = which(beta_hat_LASSO != 0)
selected_preds_Q_value = which(beta_hat_Q_value != 0)
selected_preds_CARS_MSR = which(beta_hat_CARS_MSR != 0)
selected_preds_RSF = sort(as.numeric(gsub("X", "", # cast to numeric
       var.select(my_RSF, verbose = FALSE)$topvars
)))
selected_preds_sRSF = sort(as.numeric(gsub("X", "", # cast to numeric
  var.select(my_sRSF, verbose = FALSE)$topvars
)))
```
