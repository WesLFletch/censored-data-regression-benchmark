########################## PLOT ELBOW POINT ESTIMATION EXAMPLES ##########################

library("carSurv")
library("ggplot2")
library("latex2exp")
library("survival")

standardize = function(full_df, static_idxs){
  full_df[,-static_idxs] = scale(full_df[,-static_idxs])
  full_df = full_df[,colSums(is.na(full_df))==0]
  return(full_df)
}
full_df = standardize(readRDS("BLCA.RDS"), 1:4)
prelim_mrna = readRDS("CARS3000.RDS")
full_df = full_df[,colnames(full_df) %in% c(colnames(full_df)[1:4], prelim_mrna)]

cars_scores = carSurvScore(full_df$surv_time[,1], full_df$surv_time[,2], full_df[,-(1:4)])
cars_curve = sort(abs(cars_scores), decreasing = TRUE)

fe_med = function(z){
  n = length(z); p1 = c(1, z[1]); p2 = c(n, z[n])
  distances = sapply(1:n, \(i){
    p = c(i, z[i])
    cross = abs((p2[1] - p1[1]) * (p1[2] - p[2]) - (p1[1] - p[1]) * (p2[2] - p1[2]))
    norm = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
    return(cross / norm)
  })
  return(which.max(distances))
}
fe_msr = function(z){
  get_cum_err = function(z, split){
    df1 = data.frame(y = z[1:split], x = 1:split)
    df2 = data.frame(y = z[(split+1):length(z)], x = (split+1):length(z))
    lm1 = lm(y ~ x, data = df1); lm2 = lm(y ~ x, data = df2)
    return(sum(c(lm1$residuals^6, lm2$residuals^6)))
  }
  cum_errs = c(Inf, sapply(2:(length(z)-2), \(i)get_cum_err(z, i)))
  return(which.min(cum_errs))
}

elbow_med = fe_med(cars_curve)
elbow_msr = fe_msr(cars_curve)
threshold_med = cars_curve[elbow_med]
threshold_msr = cars_curve[elbow_msr]

ep_plot = ggplot(data.frame(x = 1:length(cars_curve), y = cars_curve), aes(x, y)) +
  geom_path(
    aes(x, y, linetype = "MED"),
    data.frame(x = rep(elbow_med, 2), y = c(threshold_med, -0.01)), linewidth = 0.8
  ) + # for elbow_med
  geom_path(
    aes(x, y, linetype = "MSR"),
    data.frame(x = rep(elbow_msr, 2), y = c(threshold_msr, -0.01)), linewidth = 0.8
  ) + # for elbow_msr
  geom_path(
    aes(x, y, linetype = "MED"),
    data.frame(x = c(elbow_med, -500), y = rep(threshold_med, 2)), linewidth = 0.8
  ) + # for threshold_med
  geom_path(
    aes(x, y, linetype = "MSR"),
    data.frame(x = c(elbow_msr, -500), y = rep(threshold_msr, 2)), linewidth = 0.8
  ) + # for threshold_msr
  geom_point(aes(x, y), shape = 16, size = 0.8) + # for drawing cars_curve
  geom_point(
    aes(x, y), data.frame(x = elbow_med, y = threshold_med), shape = 18, size = 3
  ) + # for drawing med point
  geom_point(
    aes(x, y), data.frame(x = elbow_msr, y = threshold_msr), shape = 18, size = 3
  ) + # for drawing msr point
  labs(
    x = TeX("$i$"),
    y = TeX("$\\hat{abs(\\theta)}_{(i)}$"),
    linetype = "Elbow point approximation"
  ) +
  coord_cartesian(
    xlim = c(0, length(cars_curve)),
    ylim = c(0, cars_curve[1])
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(family = "Arial", face = "bold", size = 13, color = "black"),
    axis.text = element_text(family = "Arial", face = "bold", size = 10, color = "black"),
    legend.title = element_text(family = "Arial", face = "bold", size = 10, color = "black"),
    legend.text = element_text(family = "Arial", face = "bold", size = 10, color = "black"),
    legend.position = "bottom"
  )

ggsave("elbow_point.jpeg", ep_plot, width = 2000, height = 1000, units = "px")

######################## PLOT APC HISTOGRAMS BEFORE AND AFTER PFS ########################

library("ggplot2")
library("patchwork")
library("scales")
library("survival")

standardize = function(full_df, static_idxs){
  full_df[,-static_idxs] = scale(full_df[,-static_idxs])
  full_df = full_df[,colSums(is.na(full_df))==0]
  return(full_df)
}
mrna_mat = as.matrix(standardize(readRDS("BLCA.RDS"), 1:4)[,-(1:4)])
cor_mat = cor(mrna_mat)
cor_vals_abs = data.frame(
  vals = abs(as.vector(cor_mat[upper.tri(cor_mat)]))
)
apc_plot = ggplot(cor_vals_abs, aes(x = vals)) +
  geom_histogram(color = "black", breaks = (0:20)/20) +
  stat_bin(
    aes(label = comma(after_stat(count))),
    breaks = (0:20)/20,
    geom = "text",
    hjust = -0.2,
    vjust = 0.3,
    angle = 90,
    size = 3
  ) +
  scale_x_continuous(breaks = (0:20)/20) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 8e+7)
  ) +
  labs(
    x = "Absolute Pairwise Correlations",
    y = "Frequency"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(face = "bold", size = 8),
  )
trimmed_mrna = readRDS("CARS3000.RDS")
cor_mat_red = cor_mat[trimmed_mrna,trimmed_mrna]
cor_vals_abs_red = data.frame(
  vals = abs(as.vector(cor_mat_red[upper.tri(cor_mat_red)]))
)
apc_plot_red = ggplot(cor_vals_abs_red, aes(x = vals)) +
  geom_histogram(color = "black", breaks = (0:20)/20) +
  stat_bin(
    aes(label = comma(after_stat(count))),
    breaks = (0:20)/20,
    geom = "text",
    hjust = -0.2,
    vjust = 0.3,
    angle = 90,
    size = 3
  ) +
  scale_x_continuous(breaks = (0:20)/20) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 2.1e+6)
  ) +
  labs(
    x = "Absolute Pairwise Correlations",
    y = "Frequency"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(face = "bold", size = 8),
  )
apc_plot_duo = apc_plot / apc_plot_red
ggsave("pairwise_cors_duo.jpeg", apc_plot_duo, width = 2000, height = 3000, units = "px")

############################### PLOT SURVIVAL CURVE BY SEX ###############################

library("ggplot2")
library("survival")
library("survminer")

standardize = function(full_df, static_idxs){
  full_df[,-static_idxs] = scale(full_df[,-static_idxs])
  full_df = full_df[,colSums(is.na(full_df))==0]
  return(full_df)
}
full_df = standardize(readRDS("BLCA.RDS"), 1:4)
full_df$is_male = ifelse(full_df$is_male == 1, "Male", "Female")
plot_obj = ggsurvplot(survfit(surv_time ~ as.factor(is_male), data = full_df),
                      pval = TRUE, conf.int = TRUE, linetype = c("solid", "dotted"),
                      palette = c("black", "black"),
                      surv.median.line = "v", ggtheme = theme_bw(),
                      legend.title = "Sex", legend.labs = c("Female", "Male"),
                      censor.shape = 124, censor.size = 7,
                      pval.size = 10, pval.coord = c(4300, 0.9))
plot_obj$plot = plot_obj$plot + ylab("") + xlab("Days past diagnosis")
png(file = "CSF.jpeg", width=1000, height=600)
plot_obj$plot + theme(
  axis.text = element_text(size = 16, family = "Arial", face = "bold"),
  axis.title.x = element_text(size = 20, family = "Arial", face = "bold"),
  axis.title.y = element_text(size = 20, family = "Arial", face = "bold"),
  legend.text = element_text(size = 16, family = "Arial", face = "bold"),
  legend.title = element_text(size = 20, family = "Arial", face = "bold")
) + guides(color = "none", fill = "none") + labs("is_male" = "Sex")
dev.off()
