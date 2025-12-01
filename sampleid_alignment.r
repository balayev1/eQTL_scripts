suppressPackageStartupMessages({
    library(data.table)})

args <- commandArgs(trailingOnly = TRUE)
eur_genotype_counts  <- args[1]
eas_genotype_counts  <- args[2]
asv_counts <- args[3]
bai_counts <- args[4]
ciseqtl_dir <- args[5]
rna_bai_anno <- args[6]
wgs_bai_anno <- args[7]

### working matrix--EQTL directory
eur_ciseqtl_dir <- file.path(ciseqtl_dir, "eur")
eas_ciseqtl_dir <- file.path(ciseqtl_dir, "eas")

if (!dir.exists(eur_ciseqtl_dir)) {
  dir.create(eur_ciseqtl_dir, recursive = TRUE)
}
if (!dir.exists(eas_ciseqtl_dir)) {
  dir.create(eas_ciseqtl_dir, recursive = TRUE)
}

# EUR
### load SNP allele count statistics
snps.eur <- fread(eur_genotype_counts, check.names = FALSE, header = TRUE)

sample_ids <- snps.eur$IID
sample_ids <- sub("[-_][Bb]lood.*", "", sample_ids)
sample_ids <- gsub("chord", "Chord_", sample_ids)

snps.eur <- as.matrix(snps.eur[, !c("FID","IID","PAT","MAT","SEX","PHENOTYPE"), with = FALSE])
mode(snps.eur) <- "numeric"
snps.eur.t <- t(snps.eur)
colnames(snps.eur.t) <- sample_ids

snps.eur.df <- data.frame(
  snpid = rownames(snps.eur.t),
  snps.eur.t,
  check.names = FALSE,
  row.names = NULL)

rm(snps.eur, snps.eur.t)

### load gene expression matrix
counts.chordasv <- fread(asv_counts)
sample_ids_counts <- colnames(counts.chordasv)[8:ncol(counts.chordasv)]
sample_ids_counts <- basename(sample_ids_counts)
sample_ids_counts <- sub("Aligned\\.sortedByCoord\\.out\\.bam$", "", sample_ids_counts)
sample_ids_counts <- sub(".*(Chord.*)", "\\1", sample_ids_counts)
sample_ids_counts <- gsub("-", "_", sample_ids_counts)
colnames(counts.chordasv)[8:ncol(counts.chordasv)] <- sample_ids_counts

#### remove low-quality and redundant samples (keep duplicate copy with higher read depth)
counts.chordasv.sub <- counts.chordasv[, c(intersect(sample_ids_counts, colnames(counts.chordasv))), with = FALSE]
# colSums(counts.chordasv.sub)
low.quality.samples <- c("Chord_97_S17", "Chord_66_S78", "Chord_96b_S16", "Chord_37_S47", "Chord_1_S1", "Chord_22_S30", "Chord_77_S87",
    "Chord_20a_S27", "Chord_67_S79", "Chord185_ffpe_S13", "Chord_35_S45", "Chord_40_S50", "Chord_62_S74", "Chord_78_S88", "Chord188_ffpe_S16",
    "Chord190b_S18", "Chord136_M1_S26")
duplicate.samples <- c("Chord_1_S1", "Chord_104_ffpe_S14", "Chord_14a_S17", "Chord_26_S35", "Chord_28_S37", "Chord_32_S42", "Chord_55b_ffpe_S9",
    "Chord10e_TP19_P530Tissue_S85", "Chord150_Sacral_S2", "Chord178_I_S3", "Chord178_M_S4", "Chord178_V_S6", "Chord195_M1_ffpe_S23",
    "Chord195_M2_ffpe_S1", "Chord195_M3_ffpe_S24", "Chord_99_ffpe_S21", "Chord136_M2_S27", "Chord136_M4_S29", "Chord136_M5_S30", "Chord136_M6_S31",
    "Chord136_M7_S32", "Chord136_M8_S33", "Chord136_M9_S34", "Chord141_M1_S39", "Chord141_M2_S40")

counts.chordasv.sub <- counts.chordasv.sub[, !c(low.quality.samples, duplicate.samples), with = FALSE]
clean_ids <- grep("^Chord", colnames(counts.chordasv.sub), value = TRUE)
counts.chordasv.sub <- counts.chordasv.sub[, c(clean_ids), with = FALSE]
setnames(counts.chordasv.sub, old = clean_ids, new = sub("^(Chord)_?([0-9]+[A-Za-z]*).*", "Chord_\\2", clean_ids))
counts.chordasv.sub <- data.table(geneid = counts.chordasv$Geneid, counts.chordasv.sub)

#### align sample IDs between genotype and expression matrices
expr_samples <- colnames(counts.chordasv.sub)[-1]
geno_samples <- colnames(snps.eur.df)[-1]
common_samples <- expr_samples[expr_samples %in% geno_samples]
expr.common <- counts.chordasv.sub[, c("geneid", common_samples), with = FALSE]
geno.common <- snps.eur.df[, c("snpid", common_samples)]

write.table(expr.common,
            file = file.path(eur_ciseqtl_dir, "eur_expression_matrix.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)
write.table(geno.common,
            file = file.path(eur_ciseqtl_dir, "eur_genotype_matrix.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

# EAS
### load SNP allele count statistics
snps.eas <- fread(eas_genotype_counts, check.names = FALSE, header = TRUE)
sample_ids <- snps.eas$IID
sample_ids <- sub("[-_][Bb]lood.*", "", sample_ids)
sample_ids <- gsub("chord", "Chord_", sample_ids)
sample_ids <- gsub("_sorted_MD_BQ", "", sample_ids)

snps.eas <- as.matrix(snps.eas[, !c("FID","IID","PAT","MAT","SEX","PHENOTYPE"), with = FALSE])
mode(snps.eas) <- "numeric"
snps.eas.t <- t(snps.eas)
colnames(snps.eas.t) <- sample_ids

snps.eas.df <- data.frame(
  snpid = rownames(snps.eas.t),
  snps.eas.t,
  check.names = FALSE,
  row.names = NULL)

rm(snps.eas, snps.eas.t)

### load gene expression matrix
counts.chordbai <- fread(bai_counts)
sample_ids_counts <- colnames(counts.chordbai)[8:ncol(counts.chordbai)]
sample_ids_counts <- basename(sample_ids_counts)
sample_ids_counts <- sub("Aligned\\.sortedByCoord\\.out\\.bam$", "", sample_ids_counts)

map_id2subj_rna <- fread(rna_bai_anno)
sample_ids_counts <- map_id2subj_rna$biospecimen_repository_sample_id[match(sample_ids_counts, map_id2subj_rna$Run)]
colnames(counts.chordbai)[8:ncol(counts.chordbai)] <- sample_ids_counts

#### align sample IDs between genotype and expression matrices
map_id2subj_wgs <- fread(wgs_bai_anno)
map_id2subj_wgs$V1 <- gsub(".bqsr.bam", "", basename(map_id2subj_wgs$V1))

for (i in seq_along(colnames(snps.eas.df)[-1])) {
  if (colnames(snps.eas.df)[-1][i] %in% map_id2subj_wgs$V1) {
    colnames(snps.eas.df)[-1][i] <- paste0("T", map_id2subj_wgs$V3[which(map_id2subj_wgs$V1 == colnames(snps.eas.df)[-1][i])])
  }
}
geno_samples <- colnames(snps.eas.df)[-1]

asv_samples <- colnames(counts.chordasv.sub)[-1][colnames(counts.chordasv.sub)[-1] %in% geno_samples]
bai_samples <- colnames(counts.chordbai)[8:ncol(counts.chordbai)]
expr_samples <- c(asv_samples, bai_samples)
expr.common <- cbind(counts.chordasv.sub[, c("geneid", asv_samples), with = FALSE],
                      counts.chordbai[, bai_samples, with = FALSE])

all(expr_samples %in% colnames(snps.eas.df))
# [1] TRUE
geno.common <- snps.eas.df[, c("snpid", expr_samples)]

write.table(expr.common,
            file = file.path(eas_ciseqtl_dir, "eas_expression_matrix.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)
write.table(geno.common,
            file = file.path(eas_ciseqtl_dir, "eas_genotype_matrix.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)