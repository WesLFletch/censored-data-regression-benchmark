############################ PERFORM PFS SENSITIVITY ANALYSIS ############################

library("carSurv")
library("doBy")
library("dplyr")
library("survival")

# set parameters
cohort_file = "BLCA.RDS"
PFS_file = "CARS3000.RDS"
pos_controls_char = c(
  "IGF2BP1", "MMP19", "MMP20", "MMP27", "MMP8", "FOXA3", "FOXF1", "FOXF2", "FOXG1",
  "FOXI1", "FOXK2", "FOXL2", "FOXN4", "FOXP3", "FOXR1", "FGF14", "FGF22", "FGF5",
  "FGFRL1", "TOP1", "TOP3B", "PSMD4", "PSMD13", "SQLE", "PLA2G3", "PLA2G4A", "PLA2G4D",
  "ATG5", "IGF1R", "DLEU1"
)
nfolds = 10

# read and prepare unscaled full data
unscaled_df = readRDS(cohort_file)[,-1] %>%
  {\(df){
    df[,-(1:3)] %>%
      # drop features with 1-10 unique values (will likely be constant in 9/10 folds)
      {\(x)x[,!apply(x,2,\(col)length(unique(col))%in%1:10)]}() %>%
      data.frame(surv_time = df$surv_time, .)
  }}()

# determine the folds (matching those of the RDA, using the same seed)
set.seed(123)
n = nrow(unscaled_df)
p = ncol(unscaled_df)-1
row_nums = sample(1:n)
fold_idxs = lapply(1:nfolds, \(f){
  sort(row_nums[(floor((f-1)*n/nfolds)+1):floor(f*n/nfolds)])
})

# perform PFS on the ten folds, record the selected features
raw_list = lapply(fold_idxs, \(oob){
  X = scale(as.matrix(unscaled_df[-oob,-1]))
  Y = unscaled_df$surv_time[-oob]
  selected_features = carSurvScore(Y[,1], Y[,2], X) %>%
    abs() %>%
    which.maxn(3000) %>%
    {\(idx)colnames(X)[idx]}()
})

# compute pairwise Dice coefficient between the ten folds
pairwise_dice_coefs = matrix(NA, nfolds, nfolds)
for (i in 1:(nfolds-1)){
  for (j in (i+1):nfolds){
    pairwise_dice_coefs[i,j] = length(intersect(raw_list[[i]], raw_list[[j]]))/3000
} }
# compute Dice coefficients between the ten folds and the set selected by full data PFS
comparison_dice_coefs = sapply(raw_list, \(features){
  length(intersect(features, readRDS(PFS_file)))/3000
})
# compute number of true positives per fold
fold_tp = sapply(raw_list, \(features){
  length(intersect(features, pos_controls_char))
})
