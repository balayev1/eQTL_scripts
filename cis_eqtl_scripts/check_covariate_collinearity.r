##########################
##################################### Check covariate collineary for cis-EQTL run

### load R libraries
suppressPackageStartupMessages({
    library(data.table)})

cov.eur <- as.data.frame(fread("/home/venteicher_30050/balay011/ciseqtl_analysis/eur/eur_covariates_matrix.txt"))
cov.eas <- as.data.frame(fread("/home/venteicher_30050/balay011/ciseqtl_analysis/eas/eas_covariates_matrix.txt"))
colnames(cov.eur)[1] <- "id"
colnames(cov.eas)[1] <- "id"

rownames(cov.eur) <- cov.eur$id
cov.eur$id <- NULL
rownames(cov.eas) <- cov.eas$id
cov.eas$id <- NULL

cov.eur <- t(cov.eur) 
cov.eas <- t(cov.eas) 

cor_eur_mat <- cor(cov.eur, use = "complete.obs")
cor_eas_mat <- cor(as.data.frame(cov.eas), use = "complete.obs")