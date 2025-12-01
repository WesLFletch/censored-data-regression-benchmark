#################### READ ALL RESULTS AND CONSTRUCT FIGURES AND TABLES ###################

setwd("/FILEPATH/") # MODIFY

##################################### MAKE S1 FIGURES ####################################

# three figures are made, one for the two FS metrics, one for the two predictive metrics,
# and the last for computation time, and also create accompanying table for summary
# statistics across all methods

# first make figures

library("ggplot2")
library("patchwork")

s1_results = readRDS("Setting1Results.RDS")
methods = factor(c(
  "LASSO", "ALASSO", "ENET", "CB", "RSF", "sRSF", "BH", "QV", "CARS (MED)", "CARS (MSR)"
), levels = c(
  "LASSO", "ALASSO", "ENET", "CB", "RSF", "sRSF", "BH", "QV", "CARS (MED)", "CARS (MSR)"
))
sparsities = factor(c("s == 0.02", "s == 0.05", "s == 0.1"),
                    levels = c("s == 0.02", "s == 0.05", "s == 0.1"))
alphas = factor(c("alpha == 0", "alpha == 0.5"), levels = c("alpha == 0", "alpha == 0.5"))
sig_strengths = factor(c("gamma == 0.5", "gamma == 1", "gamma == 2"),
                       levels = c("gamma == 0.5", "gamma == 1", "gamma == 2"))
metrics = factor(c("FDR", "F1-Score", "CI", "Brier score", "RMSE", "Computation time"),
                 levels = c("FDR", "F1-Score", "CI", "Brier score", "RMSE", "Computation time"))
{
  s1_raw_df = data.frame()
  for (i in 1:length(methods)){
    for (j in 1:length(sparsities)){
      for (k in 1:length(alphas)){
        for (l in 1:length(sig_strengths)){
          for (m in 1:length(metrics)){
            s1_raw_df = rbind(s1_raw_df, data.frame(
              Method = methods[i],
              sparsity = sparsities[j],
              alpha = alphas[k],
              gamma = sig_strengths[l],
              metric = metrics[m],
              value = s1_results[[i]][[j]][[k]][[l]][[m]]
            ))
  } } } } }
}
make_plot = function(metric){
  plot_df = s1_raw_df[s1_raw_df$metric==metric,]
  ggplot(plot_df) +
    geom_boxplot(
      aes(y = value, fill = Method),
      outlier.size = 0.01, whisker.linewidth = 0.1, staple.linewidth = 0.1,
      box.linewidth = 0.1, median.linewidth = 0.2
    ) +
    facet_grid(gamma ~ alpha + sparsity, scales = "free", labeller = label_parsed) +
    ylab("") +
    xlab("") +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(face = "bold", size = 5),
      axis.ticks.x = element_blank(),
      axis.title.x = element_text(face = "bold", size = 10),
      legend.position = "right",
      legend.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      panel.spacing = unit(0.05, "lines"),
      plot.title = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold")
    )
}

FDR_plot = make_plot(metrics[1]) + ylab("FDR") + theme(legend.position = "none") +
  ggtitle("Setting-I: Feature selection metrics and computation time")
F1_plot = make_plot(metrics[2]) + ylab("F1-score")
CI_plot = make_plot(metrics[3]) + ylab("CI") + theme(legend.position = "none") +
  ggtitle("Setting-I: Predictive metrics")
brier_plot = make_plot(metrics[4]) + ylab("Brier score")
RMSE_plot = make_plot(metrics[5]) + ylab("RMSE") + theme(legend.position = "none")
comptime_plot = make_plot(metrics[6]) + ylab("Computation time") +
  theme(legend.position = "none")

FSCT_plots = FDR_plot / F1_plot / comptime_plot
pred_plots = CI_plot / brier_plot / RMSE_plot
ggsave("S1FSCTPlots.jpeg", FSCT_plots, width = 2000, height = 2500, units = "px")
ggsave("S1PredPlots.jpeg", pred_plots, width = 2000, height = 2500, units = "px")

##################################### MAKE S1 TABLES #####################################

# next make table of summary statistics of the form `median(IQR)`

# row is methods inside of metrics
# col is signal strengths inside of sparsities inside of alphas

s1_results = readRDS("Setting1Results.RDS")
methods = c(
  "LASSO", "ALASSO", "ENET", "CB", "RSF", "sRSF", "BH", "QV", "CARS (MED)", "CARS (MSR)"
)
sparsities = c("0.02", "0.05", "0.1")
alphas = c("0", "0.5")
sig_strengths = c("0.5", "1", "2")
metrics = c("FDR", "F1-Score", "CI", "Brier score", "RMSE", "Computation time")
{
  s1_table_df = as.data.frame(matrix(
    "",
    nrow = length(methods)*length(metrics),
    ncol = length(alphas)*length(sparsities)*length(sig_strengths)
  ))
  for (i in 1:length(methods)){
    for (j in 1:length(sparsities)){
      for (k in 1:length(alphas)){
        for (l in 1:length(sig_strengths)){
          for (m in 1:length(metrics)){
            r = (m-1)*length(methods)+i
            c = (k-1)*length(sparsities)*length(sig_strengths)+(j-1)*length(sig_strengths)+l
            s1_table_df[r,c] = paste0(
              round(median(s1_results[[i]][[j]][[k]][[l]][[m]], na.rm = TRUE), 2),
              "(",
              round(IQR(s1_results[[i]][[j]][[k]][[l]][[m]], na.rm = TRUE), 2),
              ")"
            )
  } } } } }
}
write.csv(s1_table_df, file = "S1AggTable.csv")

##################################### MAKE S2 TABLES #####################################

{
  s2_results = readRDS("Setting2Results.RDS")
  methods = names(s2_results)
  metrics = names(s2_results[[1]])
  for (i in methods){
    for (j in metrics){
      s2_results[[i]][[j]] = c(
        median = median(s2_results[[i]][[j]], na.rm = TRUE),
        IQR = IQR(s2_results[[i]][[j]], na.rm = TRUE)
      )
  } }
}
{
  s2_agg_df = as.data.frame(matrix(nrow = length(methods), ncol = length(metrics)))
  rownames(s2_agg_df) = methods
  colnames(s2_agg_df) = metrics
  for (i in 1:length(methods)){
    for (j in 1:length(metrics)){
      s2_agg_df[i,j] = paste0(
        round(s2_results[[i]][[j]][1], 2),
        "(",
        round(s2_results[[i]][[j]][2], 2),
        ")"
      )
  } }
}
write.csv(s2_agg_df, file = "S2AggTable.csv")

##################################### MAKE RDA TABLES ####################################

rda_raw = readRDS("RDAResults.RDS")
rda_features = readRDS("RDASelectedFeatures.RDS")
methods = names(rda_raw)
metrics = names(rda_raw[[1]])
nfolds = length(rda_features[[1]])
get_dice_coef = function(a, b) 2*sum(a&b)/(sum(a)+sum(b))
dice_mats = lapply(1:length(methods), \(m){
  out_mat = matrix(NA, nfolds, nfolds)
  for (i in 1:(nfolds-1)){
    for (j in (i+1):nfolds){
      out_mat[i,j] = get_dice_coef(rda_features[[m]][[i]], rda_features[[m]][[j]])
    } }
  out_mat
})
dice_chars = sapply(dice_mats, \(mat)paste0(
  round(median(mat[upper.tri(mat)], na.rm = TRUE), digits = 2), "(",
  round(IQR(mat[upper.tri(mat)], na.rm = TRUE), digits = 2), ")"
))
{
  rda_agg_df = as.data.frame(matrix(nrow = length(methods), ncol = length(metrics)))
  rownames(rda_agg_df) = methods
  colnames(rda_agg_df) = metrics
  for (i in 1:length(methods)){
    for (j in 1:length(metrics)){
      rda_agg_df[i,j] = paste0(
        round(median(rda_raw[[i]][[j]], na.rm = TRUE), 2),
        "(",
        round(IQR(rda_raw[[i]][[j]], na.rm = TRUE), 2),
        ")"
      )
  } }
  rda_agg_df$"Dice coefficient" = dice_chars
  rda_agg_df = rda_agg_df[,c(1,2,7,3,4,5,6)]
}
write.csv(rda_agg_df, file = "RDAAggTable.csv")

########################## MAKE CUMULATIVE METHOD RANKINGS TABLE #########################

s1_results = readRDS("Setting1Results.RDS")
s2_results = readRDS("Setting2Results.RDS")
methods = c(
  "LASSO", "ALASSO", "ENET", "CB", "RSF", "sRSF", "BH", "QV", "CARS (MED)", "CARS (MSR)"
)
sparsities = c("0.02", "0.05", "0.1")
alphas = c("0", "0.5")
sig_strengths = c("0.5", "1", "2")
metrics = c("FDR", "F1-Score", "CI", "Brier score", "RMSE", "Computation time")
{
  s1_medians = array(dim = c(length(methods), length(sparsities), length(alphas),
                             length(sig_strengths), length(metrics)))
  for (i in 1:length(methods)){
    for (j in 1:length(sparsities)){
      for (k in 1:length(alphas)){
        for (l in 1:length(sig_strengths)){
          for (m in 1:length(metrics)){
            s1_medians[i,j,k,l,m] = median(s1_results[[i]][[j]][[k]][[l]][[m]], na.rm = TRUE)
  } } } } }
}
{
  s2_medians = array(dim = c(length(methods), length(metrics)))
  for (i in 1:length(methods)){
    for (j in 1:length(metrics)){
      s2_medians[i,j] = median(s2_results[[i]][[j]], na.rm = TRUE)
  } }
}
{
  s1_ranks_nonagg = array(dim = c(length(methods), length(sparsities), length(alphas),
                                  length(sig_strengths), length(metrics)))
  for (j in 1:length(sparsities)){
    for (k in 1:length(alphas)){
      for (l in 1:length(sig_strengths)){
        for (m in 1:length(metrics)){
          # some metrics have high is good, some low is good
          high_is_good = m %in% c(2, 3)
          s1_ranks_nonagg[,j,k,l,m] = rank(s1_medians[,j,k,l,m])
          if (high_is_good)
            s1_ranks_nonagg[,j,k,l,m] = length(methods)+1-s1_ranks_nonagg[,j,k,l,m]
  } } } }
}
{
  s2_ranks = array(dim = c(length(methods), length(metrics)))
  for (j in 1:length(metrics)){
    # some metrics have high is good, some low is good
    high_is_good = j %in% c(2, 3)
    s2_ranks[,j] = rank(s2_medians[,j])
    if (high_is_good) s2_ranks[,j] = length(methods)+1-s2_ranks[,j]
  }
  s2_ranks = as.data.frame(s2_ranks)
  rownames(s2_ranks) = methods
  colnames(s2_ranks) = metrics
}
s1_ranksums = apply(s1_ranks_nonagg, c(1,5), sum)
{
  s1_ranksum_ranks = as.data.frame(apply(s1_ranksums, 2, rank))
  rownames(s1_ranksum_ranks) = methods
  colnames(s1_ranksum_ranks) = metrics
}
sim_ranks_agg = rbind(s1_ranksum_ranks, s2_ranks)
write.csv(sim_ranks_agg, file = "SimRanksAgg.csv")
