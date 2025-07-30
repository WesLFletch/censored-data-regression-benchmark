############################### REAL DATA ANALYSIS CONTROL ###############################

source("LassoRDA.R")
source("AlassoRDA.R")
source("CoxBoostRDA.R")
source("RandomForestRDA.R")
source("BHProcedureRDA.R")
source("QvalueRDA.R")
source("CarsRDA.R")

run_RDA = function(methods, files, prelim){
  cluster = makeCluster(detectCores() - 1)
  registerDoParallel(cluster)
  methods_list = foreach(i = 1:length(methods), .packages = c(
    "doBy", "doParallel", "dplyr", "survival", "tidyr", "yardstick")) %dopar% {
      source("AnalysisControl.R")
      lapply(1:length(files), \(j)run_method_on_file(methods[i], files[j], prelim))
    }
  stopCluster(cluster)
  files_list = lapply(1:length(files), \(i)NA)
  names(files_list) = files
  for (j in 1:length(files_list)){
    files_list[[j]] = lapply(1:length(methods), \(i)methods_list[[i]][[j]])
    names(files_list[[j]]) = methods
  }
  return(files_list)
}

run_method_on_file = function(method_name, file, prelim){
  method = method_names_list[[method_name]]
  full_df = standardize(readRDS(file), 1:4)
  full_df = full_df[,c(rep(TRUE, 4), sapply(5:ncol(full_df), \(i)IQR(full_df[,i])!=0))]
  n = nrow(full_df)
  fold_idxs = lapply(1:10, \(i){
    (floor((i-1)*n/10)+1):floor(i*n/10)
  })
  return(method(full_df, fold_idxs, readRDS(prelim)))
}

standardize = function(full_df, static_idxs){
  full_df[,-static_idxs] = scale(full_df[,-static_idxs])
  full_df = full_df[,colSums(is.na(full_df))==0]
  return(full_df)
}

method_names_list = list(
  "LASSO" = LASSO, "ALASSO" = ALASSO, "CB" = CB, "RF_md" = RF_md,
  "RF_md_screen" = RF_md_screen, "BHP" = BHP, "QV" = QV, "CARS_med" = CARS,
  "CARS_msr" = CARS_residuals
)
