
setwd("/home/jsjukara/")
library(data.table)
library(tidyverse)

##export_hr, data frame containing columns:
#"rsid": identifier of the genetic exposure
#"cause": the disease
#"loghr_m": the log hazard ratio (meta-analyzed for FinnGen and UK Biobank)
export_hr <- readRDS(file="/home/jsjukara/finngen_data/export_hr.rda")


##reg_df, data frame containing columns:
#"rsid": identifier of the genetic exposure
#"cause": the disease
#"shrunk": indicator of whether the genetic exposure-disease HR was shrunk in the shrinkage
reg_df <- readRDS(file="/home/jsjukara/finngen_data/export_hr.rda")

#convert to wide format
hr_wide <- export_hr %>%
    select(rsid, cause, loghr_m) %>%
    pivot_wider(names_from="cause", values_from="loghr_m")
rownames(hr_wide) <- hr_wide$rsid


#estimate correlation matrix of betas
corr_mat <- matrix(0, ncol=length(causes), nrow=length(causes))
colnames(corr_mat) <- causes
rownames(corr_mat) <- causes

for (i in 1:nrow(corr_mat)) {
    for (j in 1:ncol(corr_mat)) {
        causeA <- colnames(corr_mat)[i]
        causeB <- colnames(corr_mat)[j]
        
        tempA <- reg_df %>%
            filter(cause == causeA, shrunk==1)
        rsidA <- tempA$rsid
        
        tempB <- reg_df %>%
            filter(cause == causeB, shrunk==1)
        rsidB <- tempB$rsid
        
        rsid_not_shrunk <- intersect(rsidA, rsidB)
        corr_data <- hr_wide %>%
            filter(rsid %in% rsid_not_shrunk) %>%
            select(causeA, causeB)
        corr_mat[i,j] <- cor(corr_data[,causeA], corr_data[,causeB])
    }
}

saveRDS(corr_mat, "/home/jsjukara/finngen_data/beta_correlation_matrix.rds")
