############ COMPARE STABILITY OF MED AND MSR EBLOW POINT ESTIMATOR ALGORITHMS ###########

library("carSurv")
library("survival")

make_data = function(seed, alpha, n, p, s, sig_str){
  # make data
  set.seed(seed = seed)
  Sigma = matrix(alpha, p, p)
  diag(Sigma) = 1
  A = t(chol(Sigma))
  beta = rep(0, p)
  beta[1:floor(s/2)] = sig_str
  beta[(floor(s/2)+1):s] = -sig_str
  z = matrix(rnorm(n*p), nrow = p)
  x = t(A %*% matrix(rnorm(n*p), nrow = p))
  ty = rexp(n, exp(-(x %*% as.matrix(beta))))
  cy = runif(n, quantile(ty, prob=0.5), quantile(ty, prob=0.90))
  df = as.data.frame(x)
  df$surv_time = Surv(apply(cbind(ty,cy), 1, min), as.numeric(ty<cy))
  return(list(x = x, y = df$surv_time))
}

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

# parameters match those in setting-i of our simulations.
# n=200, the number of observations in the training sets in setting-i.
# datasets may not be exactly the same as in setting-i, but the same data characteristics
# and a large number of replications are used so the distributions should be nearly
# identical.

num_trials = 200
n = 200
p = 1000
s_vals = c(20, 50, 100)
alpha_vals = c(0, 0.5)
sig_str_vals = c(0.5, 1, 2)

sig_str_list = lapply(sig_str_vals, \(i)rep(NA, num_trials))
names(sig_str_list) = as.character(sig_str_vals)
alpha_list = lapply(alpha_vals, \(i)sig_str_list)
names(alpha_list) = as.character(alpha_vals)
s_list = lapply(s_vals, \(i)alpha_list)
names(s_list) = as.character(s_vals)
ep_list = lapply(1:2, \(i)s_list)
names(ep_list) = c("MSR", "MED")

for (i in 1:length(s_vals)){
  for (j in 1:length(alpha_vals)){
    for (k in 1:length(sig_str_vals)){
      ep_vals = sapply(1:num_trials, \(seed){
        data = make_data(seed, alpha_vals[j], n, p, s_vals[i], sig_str_vals[k])
        x = data$x
        y = data$y
        cars_scores = carSurvScore(y[,1], y[,2], x)
        cars_curve = sort(abs(cars_scores), decreasing = TRUE)
        return(c(msr = MSR(cars_curve), med = MED(cars_curve)))
      })
      ep_list$MSR[[i]][[j]][[k]] = ep_vals["msr",]
      ep_list$MED[[i]][[j]][[k]] = ep_vals["med",]
    }
  }
}

saveRDS(ep_list, "elbows.RDS")

############################ PLOT THE ELBOW POINT ESTIMATIONS ############################

library("ggplot2")
ep_list = readRDS("elbows.RDS")
plot_df = data.frame(Method = character(0), sparsity = character(0), value = numeric(0))
s_vals = factor(c("s==0.02", "s==0.05", "s==0.1"), levels = c("s==0.02", "s==0.05", "s==0.1"))
alpha_vals = c(0, 0.5)
sig_str_vals = c(0.5, 1, 2)
for (i in 1:length(s_vals)){
  for (j in 1:length(alpha_vals)){
    for (k in 1:length(sig_str_vals)){
      plot_df = rbind(
        plot_df,
        data.frame(Method = "MSR", sparsity = s_vals[i], value = ep_list$MSR[[i]][[j]][[k]]),
        data.frame(Method = "MED", sparsity = s_vals[i], value = ep_list$MED[[i]][[j]][[k]])
      )
    }
  }
}
plot_obj = ggplot(plot_df) +
  geom_boxplot(aes(y = value, fill = Method), outlier.size = 0.4) +
  facet_wrap(vars(sparsity), labeller = label_parsed) +
  ylab("Selected features") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold")
  )
ggsave("elbow_stability.jpeg", plot_obj, width = 1300, height = 800, units = "px")
