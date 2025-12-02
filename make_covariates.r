suppressPackageStartupMessages({
  library(data.table)})

args <- commandArgs(trailingOnly = TRUE)
eur_geno_pca_file <- args[1]
eas_geno_pca_file <- args[2]
eur_exprs_pca_file <- args[3]
eas_exprs_pca_file <- args[4]
eur_tech_file <- args[5]
eas_tech_file <- args[6]
eur_expr_file <- args[7]
eas_expr_file <- args[8]
bai_anno_wgs <- args[9]
npc_geno <- args[10]
npc_exprs <- args[11]

# eur_geno_pca_file <- "/home/venteicher_30050/balay011/ciseqtl_analysis/eur/eur_pca_results.eigenvec"
# eas_geno_pca_file <- "/home/venteicher_30050/balay011/ciseqtl_analysis/eas/eas_pca_results.eigenvec"
# eur_exprs_pca_file <- "/home/venteicher_30050/balay011/ciseqtl_analysis/eur/eur_expr_pcs.txt"
# eas_exprs_pca_file <- "/home/venteicher_30050/balay011/ciseqtl_analysis/eas/eas_expr_pcs.txt"
# eur_tech_file <- "/home/venteicher_30050/balay011/gwas_essentials/gwas_covariates/gwas_covariates_eur.txt"
# eas_tech_file <- "/home/venteicher_30050/balay011/gwas_essentials/gwas_covariates/gwas_covariates_eas.txt"
# eur_expr_file <- "/home/venteicher_30050/balay011/ciseqtl_analysis/eur/eur_expression_matrix_vst.txt"
# eas_expr_file <- "/home/venteicher_30050/balay011/ciseqtl_analysis/eas/eas_expression_matrix_vst.txt"
# bai_anno_wgs <- "/projects/standard/aventeic/balay011/scripts/WGS_chordomaBai_anno.txt"
# npc_geno <- 10
# npc_exprs <- 10

npc_geno <- as.numeric(npc_geno)
npc_exprs <- as.numeric(npc_exprs)

# set up output directories
eur_ciseqtl_dir <- dirname(eur_exprs_pca_file)
eas_ciseqtl_dir <- dirname(eas_exprs_pca_file)

# EUR
## load gwas covariates
covar.eur <- fread(eur_tech_file)

## select only SEX and COVERAGE
covar.eur <- covar.eur[, c("IID", "SEX", "COVERAGE")]

## rename sample IDs
covar.eur$IID <- sub("[-_][Bb]lood.*", "", covar.eur$IID)
covar.eur$IID <- gsub("chord", "Chord_", covar.eur$IID)

## load vst-normalized expression matrix
exprs.eur <- fread(eur_expr_file)

## subset samples within expression matrix
common_samples <- colnames(exprs.eur)[-1]
covar.common <- covar.eur[IID %in% common_samples]
covar.common <- covar.common[match(common_samples, IID)]
covar_mat <- t(as.matrix(covar.common[, !"IID"]))
colnames(covar_mat) <- common_samples

## load genotype PCs and merge with gwas covariates
pca.geno.eur <- fread(eur_geno_pca_file)

## rename sample IDs
pca.geno.eur$IID <- sub("[-_][Bb]lood.*", "", pca.geno.eur$IID)
pca.geno.eur$IID <- gsub("chord", "Chord_", pca.geno.eur$IID)
pca.geno.eur$'#FID' <- NULL
colnames(pca.geno.eur)[2:(npc_geno+1)] <- paste0(colnames(pca.geno.eur)[2:(npc_geno+1)], "geno")
pca.geno.eur <- pca.geno.eur[, 1:(npc_geno+1)]

if (all(pca.geno.eur$IID %in% colnames(covar_mat))){
  covar_mat_t <- as.data.frame(t(covar_mat))
  covar_mat_t$IID <- rownames(covar_mat_t)
  covar_mat_t <- merge(covar_mat_t, pca.geno.eur, by = "IID")
}

## load expression PCs and merge with gwas covariates
pca.exprs.eur <- fread(eur_exprs_pca_file)

## add rna-seq batch covariate
pca.exprs.eur$batch <- c(rep(1, 33), 2, rep(1, 3), rep(2, 20))

select.cols <- c("sample_id", paste0("PC", 1:npc_exprs), "batch")
pca.exprs.eur <- as.data.frame(pca.exprs.eur[, select.cols, with = FALSE])
rownames(pca.exprs.eur) <- pca.exprs.eur$sample_id
pca.exprs.eur$sample_id <- NULL
colnames(pca.exprs.eur) <- paste0(colnames(pca.exprs.eur), "exprs") 
pca.exprs.eur$IID <- rownames(pca.exprs.eur)
covar_mat_t <- as.data.frame(merge(covar_mat_t, pca.exprs.eur, by = "IID"))
rownames(covar_mat_t) <- covar_mat_t$IID
covar_mat_t$IID <- NULL
covar_mat_t <- t(covar_mat_t)

fwrite(covar_mat_t, file = file.path(eur_ciseqtl_dir, "eur_covariates_matrix.txt"), sep = "\t", quote = FALSE, row.names = TRUE)

# EAS
## load gwas covariates
covar.eas <- fread(eas_tech_file)

## select only SEX and COVERAGE
covar.eas <- covar.eas[, c("IID", "SEX", "COVERAGE")]

## rename sample IDs
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

## load genotype PCs and merge with gwas covariates
pca.geno.eas <- fread(eas_geno_pca_file)

## rename sample IDs
pca.geno.eas$IID <- sub("[-_][Bb]lood.*", "", pca.geno.eas$IID)
pca.geno.eas$IID <- gsub("chord", "Chord_", pca.geno.eas$IID)
pca.geno.eas$IID <- gsub("_sorted_MD_BQ", "", pca.geno.eas$IID)
pca.geno.eas$'#FID' <- NULL

for (i in seq_along(pca.geno.eas$IID)) {
  if (pca.geno.eas$IID[i] %in% map_id2subj_wgs$V1) {
    pca.geno.eas$IID[i] <- paste0("T", map_id2subj_wgs$V3[which(map_id2subj_wgs$V1 == pca.geno.eas$IID[i])])
  }
}

colnames(pca.geno.eas)[2:(npc_geno+1)] <- paste0(colnames(pca.geno.eas)[2:(npc_geno+1)], "geno")
pca.geno.eas <- pca.geno.eas[, 1:(npc_geno+1)]

if (all(pca.geno.eas$IID %in% colnames(covar_mat))){
  covar_mat_t <- as.data.frame(t(covar_mat))
  covar_mat_t$IID <- rownames(covar_mat_t)
  covar_mat_t <- merge(covar_mat_t, pca.geno.eas, by = "IID")
}


## load expression PCs and merge with gwas covariates
pca.exprs.eas <- fread(eas_exprs_pca_file)

## add rna-seq batch covariate
pca.exprs.eas$batch <- c(rep(1, 2), rep(2, 27))

select.cols <- c("sample_id", paste0("PC", 1:npc_exprs), "batch")
pca.exprs.eas <- as.data.frame(pca.exprs.eas[, select.cols, with = FALSE])
rownames(pca.exprs.eas) <- pca.exprs.eas$sample_id
pca.exprs.eas$sample_id <- NULL
colnames(pca.exprs.eas) <- paste0(colnames(pca.exprs.eas), "exprs") 
pca.exprs.eas$IID <- rownames(pca.exprs.eas)
covar_mat_t <- as.data.frame(merge(covar_mat_t, pca.exprs.eas, by = "IID"))
rownames(covar_mat_t) <- covar_mat_t$IID
covar_mat_t$IID <- NULL
covar_mat_t <- t(covar_mat_t)

fwrite(covar_mat_t, file = file.path(eas_ciseqtl_dir, "eas_covariates_matrix.txt"), sep = "\t", quote = FALSE, row.names = TRUE)