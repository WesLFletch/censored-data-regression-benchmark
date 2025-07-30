########################## PART I: INSTALLATION INTO LOCAL FILES #########################

setwd("YOUR-DIRECTORY-HERE") # FIXME: set to appropriate directory

# required packages
library("UCSCXenaTools")
library("dplyr")
library("survival")

# show all TCGA cohorts
# as.vector(unique(XenaData[XenaData$XenaHostNames=="tcgaHub", 3]))
# show BLCA datasets
# as.vector(XenaData[XenaData$XenaCohorts=="TCGA Bladder Cancer (BLCA)", 4:5])

# cart BLCA dataset
blca_mrna_cart=XenaGenerate(subset=XenaHostNames=="tcgaHub") %>%
  XenaFilter(filterDatasets="BLCA") %>%
  XenaFilter(filterDatasets="HiSeqV2")
blca_surv_cart=XenaGenerate(subset=XenaHostNames=="tcgaHub") %>%
  XenaFilter(filterDatasets="BLCA") %>%
  XenaFilter(filterDatasets="survival")
blca_clin_cart=XenaGenerate(subset=XenaHostNames=="tcgaHub") %>%
  XenaFilter(filterDatasets="BLCA") %>%
  XenaFilter(filterDatasets="clinical")
# install BLCA dataset (MAY TAKE A WHILE)
blca_mrna_download=XenaQuery(blca_mrna_cart) %>%
  XenaDownload(.)
blca_surv_download=XenaQuery(blca_surv_cart) %>%
  XenaDownload(.)
blca_clin_download=XenaQuery(blca_clin_cart) %>%
  XenaDownload(.)
# get data into R environment
blca_mrna_prep=XenaPrepare(blca_mrna_download)
blca_surv_prep=XenaPrepare(blca_surv_download)
blca_clin_prep=XenaPrepare(blca_clin_download)
# write tables for later use
write.table(blca_mrna_prep[[3]], "blca_mrna.txt")
# write.table(blca_mrna_prep[[4]], "blca_mrna_exon.txt")
write.table(blca_surv_prep, "blca_surv.txt")
write.table(blca_clin_prep, "blca_clin.txt")

############################### PART II: DATA MANIPULATIONS ##############################

# read and reformat blca_mrna_df
blca_mrna_df=as.data.frame(t(read.table("blca_mrna.txt", header=TRUE)))
colnames(blca_mrna_df)=blca_mrna_df[1,]
blca_mrna_df=blca_mrna_df[-1,]
blca_mrna_df$sample=rownames(blca_mrna_df)
rownames(blca_mrna_df)=NULL
blca_mrna_df=blca_mrna_df[,c(ncol(blca_mrna_df),1:(ncol(blca_mrna_df)-1))]
blca_mrna_df[,-1]=lapply(blca_mrna_df[,-1], as.numeric)

# read and reformat blca_surv_df
blca_surv_df=read.table("blca_surv.txt", header=TRUE)
blca_surv_df$sample=gsub("-", ".", blca_surv_df$sample)
blca_surv_df=blca_surv_df[,c("sample", "OS", "OS.time")]
blca_surv_df[,-1]=lapply(blca_surv_df[,-1], as.numeric)

# read and reformat blca_clin_df
blca_clin_df=read.table("blca_clin.txt", header=TRUE)
blca_clin_df$sampleID=gsub("-", ".", blca_clin_df$sampleID)
blca_clin_df=blca_clin_df[,c("sampleID", "days_to_birth", "gender")]
blca_clin_df$days_to_birth = -blca_clin_df$days_to_birth
blca_clin_df$gender=as.numeric(blca_clin_df$gender=="MALE")
colnames(blca_clin_df)[1]="sample"
colnames(blca_clin_df)[2]="age_in_days"
colnames(blca_clin_df)[3]="is_male"

# make full data frame
# join data frames together by sample id and re-arrange feature order
blca_full_df=inner_join(blca_surv_df, blca_clin_df, by=c("sample"="sample"))
blca_full_df=inner_join(blca_full_df, blca_mrna_df, by=c("sample"="sample"))
blca_full_df$surv_time=Surv(blca_full_df$OS.time, blca_full_df$OS)
blca_full_df=blca_full_df[,c(1,ncol(blca_full_df),4:(ncol(blca_full_df)-1))]
blca_full_df=blca_full_df[!(rowSums(is.na(blca_full_df)) > 0),]
blca_full_df=blca_full_df[blca_full_df$surv_time[,1] > 0,]

# STRUCTURE OF blca_full_df:
# col 1: patient's sample id
# col 2: `Surv` object containing patient's time to event and censorship indicator
# col 3: patient's age in days at the time of diagnosis
# col 4: patient sex as binary feature, 0:female, 1:male
# cols 5+: mRNA features after log2 transformation (not yet normalized to mean 0 var 1)

saveRDS(blca_full_df, file="BLCA.RDS")
