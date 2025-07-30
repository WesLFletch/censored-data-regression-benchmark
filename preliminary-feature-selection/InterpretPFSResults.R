####################### INTERPRET THE RESULTS OF THE PRELIMINARY FS ######################

setwd("YOUR-DIRECTORY-HERE") # FIXME: set to appropriate directory

prelim_fs = readRDS("PrelimFSResults.RDS")
prelim_fs_clean = prelim_fs[rowSums(is.na(prelim_fs)) == 0,]
mrna_to_keep = prelim_fs_clean$gene[prelim_fs_clean$adj_r_sq < 0.85]

saveRDS(mrna_to_keep, file = "UncorMRNA.RDS")
