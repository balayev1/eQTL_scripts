suppressPackageStartupMessages({
  library(data.table)})

args <- commandArgs(trailingOnly = TRUE)
eur_pca_file <- args[1]
eas_pca_file <- args[2]
eur_tech_file <- args[3]
eas_tech_file <- args[4]
eur_expr_file <- args[5]
eas_expr_file <- args[6]
bai_anno_wgs <- args[7]
k <- args[8]

# set up output directories
eur_ciseqtl_dir <- dirname(eur_pca_file)
eas_ciseqtl_dir <- dirname(eas_pca_file)

# EUR
## load gwas covariates
covar.eur <- fread(eur_tech_file)

## rename sample IDs
covar.eur$FID <- NULL
covar.eur$IID <- sub("[-_][Bb]lood.*", "", covar.eur$IID)
covar.eur$IID <- gsub("chord", "Chord_", covar.eur$IID)

## load vst-normalized expression matrix
exprs.eur <- fread(eur_expr_file)

## subset expression samples
common_samples <- colnames(exprs.eur)[-1]
covar.common <- covar.eur[IID %in% common_samples]
covar.common <- covar.common[match(common_samples, IID)]
covar_mat <- t(as.matrix(covar.common[, !"IID"]))
colnames(covar_mat) <- common_samples
rownames(covar_mat)[1:20] <- paste0(rownames(covar_mat)[1:20], "geno")

## load expression PCs and merge with gwas covariates
k <- 20 # number of PCs
pca.eur <- fread(eur_pca_file)

## add rna-seq batch covariate
pca.eur$batch <- c(rep(1, 33), 2, rep(1, 3), rep(2, 20))

select.cols <- c("sample_id", paste0("PC", 1:20), "batch")
pca.eur.t <- as.data.frame(pca.eur[, select.cols, with = FALSE])
rownames(pca.eur.t) <- pca.eur.t$sample_id
pca.eur.t$sample_id <- NULL
pca.eur.t <- t(pca.eur.t)
rownames(pca.eur.t) <- paste0(rownames(pca.eur.t), "exprs") 
covar_mat_final <- rbind(covar_mat, pca.eur.t)

fwrite(covar_mat_final, file = file.path(eur_ciseqtl_dir, "eur_covariates_matrix.txt"), sep = "\t", quote = FALSE, row.names = TRUE)

# EAS
## load gwas covariates
covar.eas <- fread(eas_tech_file)

## rename sample IDs
covar.eas$FID <- NULL
covar.eas$IID <- sub("[-_][Bb]lood.*", "", covar.eas$IID)
covar.eas$IID <- gsub("chord", "Chord_", covar.eas$IID)
covar.eas$IID <- gsub("_sorted_MD_BQ", "", covar.eas$IID)

map_id2subj_wgs <- fread(bai_anno_wgs)
map_id2subj_wgs$V1 <- gsub(".bqsr.bam", "", basename(map_id2subj_wgs$V1))

for (i in seq_along(covar.eas$IID)) {
  if (covar.eas$IID[i] %in% map_id2subj_wgs$V1) {
    covar.eas$IID[i] <- paste0("T", map_id2subj_wgs$V3[which(map_id2subj_wgs$V1 == covar.eas$IID[i])])
  }
}

## load vst-normalized expression matrix
exprs.eas <- fread(eas_expr_file)

## subset expression samples
common_samples <- colnames(exprs.eas)[-1]
covar.common <- covar.eas[IID %in% common_samples]
covar.common <- covar.common[match(common_samples, IID)]
covar_mat <- t(as.matrix(covar.common[, !"IID"]))
colnames(covar_mat) <- common_samples
rownames(covar_mat)[1:20] <- paste0(rownames(covar_mat)[1:20], "geno")

## load expression PCs and merge with gwas covariates
k <- 20 # number of PCs
pca.eas <- fread(eas_pca_file)

## add rna-seq batch covariate
pca.eas$batch <- c(rep(1, 2), rep(2, 27))

select.cols <- c("sample_id", paste0("PC", 1:20), "batch")
pca.eas.t <- as.data.frame(pca.eas[, select.cols, with = FALSE])
rownames(pca.eas.t) <- pca.eas.t$sample_id
pca.eas.t$sample_id <- NULL
pca.eas.t <- t(pca.eas.t)
rownames(pca.eas.t) <- paste0(rownames(pca.eas.t), "exprs") 
covar_mat_final <- rbind(covar_mat, pca.eas.t)

fwrite(covar_mat_final, file = file.path(eas_ciseqtl_dir, "eas_covariates_matrix.txt"), sep = "\t", quote = FALSE, row.names = TRUE)
