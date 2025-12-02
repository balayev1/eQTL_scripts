# How to identify cis-eQTLs for GWAS variants?

## Some steps of the pipeline are highly memory-consuming, so make sure to have enough supplied resources
srun --mem-per-cpu=128GB --time=12:00:00 --pty --cpus-per-task=1 --x11 -p interactive bash

CIS_EQTL_OUTDIR="/home/venteicher_30050/balay011/ciseqtl_analysis"

## We start by first generating allele dosage matrix from plink2 GWAS files.
### First, we do it for samples of european ancestry
EUR_GWAS_SUMMARY="/home/venteicher_30050/balay011/gwas_essentials/gwas_output/gwas_results_eur_autosomes.add.TUMOR.glm.logistic.hybrid"
EUR_GWAS_IN="/scratch.global/balay011/gwas_analysis/eur/eur_all_4gwas"
awk 'NR>1 {print $3"\t"$7}' \
  $EUR_GWAS_SUMMARY \
  > "$CIS_EQTL_OUTDIR/eur_A1_alleles_for_export.txt"

EUR_PHENO_TSV="/panfs/jay/groups/9/venteicher_30050/balay011/gwas_essentials/pheno_files/pheno_eur.txt"
conda run -n gwas_env \
  plink2 --bfile "$EUR_GWAS_IN" \
        --pheno $EUR_PHENO_TSV --pheno-name TUMOR \
        --keep-if TUMOR==2 \
         --export A \
         --export-allele $CIS_EQTL_OUTDIR/eur_A1_alleles_for_export.txt \
         --out $CIS_EQTL_OUTDIR/eur_allele_dosage

### Next, we do it for samples of eastern asian ancestry
EAS_GWAS_SUMMARY="/home/venteicher_30050/balay011/gwas_essentials/gwas_output/gwas_results_eas_autosomes.add.TUMOR.glm.logistic.hybrid"
EAS_GWAS_IN="/scratch.global/balay011/gwas_analysis/eas/eas_all_4gwas"
awk 'NR>1 {print $3"\t"$7}' \
  $EAS_GWAS_SUMMARY \
  > $CIS_EQTL_OUTDIR/eas_A1_alleles_for_export.txt

EAS_PHENO_TSV="/panfs/jay/groups/9/venteicher_30050/balay011/gwas_essentials/pheno_files/pheno_eas.txt"

conda run -n gwas_env \
  plink2 --bfile "$EAS_GWAS_IN" \
        --pheno $EAS_PHENO_TSV --pheno-name TUMOR \
        --keep-if TUMOR==2 \
         --export A \
         --export-allele $CIS_EQTL_OUTDIR/eas_A1_alleles_for_export.txt \
         --out $CIS_EQTL_OUTDIR/eas_allele_dosage

## Now, we have to align sample IDs in allele dosage matrix and expression matrix for both ancestries
BAI_ANNO_RNA="/projects/standard/aventeic/balay011/scripts/RNA_Seq_chordomaBai_anno.txt"
BAI_ANNO_WGS="/projects/standard/aventeic/balay011/scripts/WGS_chordomaBai_anno.txt"
conda run -n deseq_env \
  Rscript /scratch.global/balay011/ciseqtl_scripts/sampleid_alignment.r $CIS_EQTL_OUTDIR/eur_allele_dosage.raw $CIS_EQTL_OUTDIR/eas_allele_dosage.raw "/home/venteicher_30050/balay011/RNA_Seq_counts/CHORDASV/Counts/raw_counts.matrix.csv" "/home/venteicher_30050/balay011/RNA_Seq_counts/RNASeq_CHORDBAI/Counts/raw_counts.matrix.csv" $CIS_EQTL_OUTDIR  $BAI_ANNO_RNA $BAI_ANNO_WGS

## Next we use sample ID aligned expression matrices for normalization using vst from DESeq2
conda run -n deseq_env \
  Rscript /scratch.global/balay011/ciseqtl_scripts/apply_vst.r $CIS_EQTL_OUTDIR/eur/eur_expression_matrix.txt $CIS_EQTL_OUTDIR/eas/eas_expression_matrix.txt

## Next we extract PCs in the vst-normalized expression matrices
conda run -n deseq_env \
  Rscript /scratch.global/balay011/ciseqtl_scripts/run_pca.r $CIS_EQTL_OUTDIR/eur/eur_expression_matrix_vst.txt $CIS_EQTL_OUTDIR/eas/eas_expression_matrix_vst.txt


## Next we make table of covariates: 10 genotype PCs, 10 expression PCs, sex, whole-genome coverage and RNA-Seq batch
EUR_GENO_COVS="$CIS_EQTL_OUTDIR/eur/eur_pca_results.eigenvec"
EAS_GENO_COVS="$CIS_EQTL_OUTDIR/eas/eas_pca_results.eigenvec"
EUR_EXPRS_PCA_COVS="$CIS_EQTL_OUTDIR/eur/eur_expr_pcs.txt"
EAS_EXPRS_PCA_COVS="$CIS_EQTL_OUTDIR/eas/eas_expr_pcs.txt"
EUR_TECH_COVS="/home/venteicher_30050/balay011/gwas_essentials/gwas_covariates/gwas_covariates_eur.txt"
EAS_TECH_COVS="/home/venteicher_30050/balay011/gwas_essentials/gwas_covariates/gwas_covariates_eas.txt"
NPC_GENO=10
NPC_EXPRS=10

conda run -n deseq_env \
  Rscript /scratch.global/balay011/ciseqtl_scripts/make_covariates.r $EUR_GENO_COVS $EAS_GENO_COVS $EUR_EXPRS_PCA_COVS $EAS_EXPRS_PCA_COVS $EUR_TECH_COVS $EAS_TECH_COVS $CIS_EQTL_OUTDIR/eur/eur_expression_matrix_vst.txt $CIS_EQTL_OUTDIR/eas/eas_expression_matrix_vst.txt $BAI_ANNO_WGS $NPC_GENO $NPC_EXPRS

## Last we run MatrixEQTL to find cis-EQTLs
conda run -p $MSIPROJECT/balay011/.conda/envs/r433_env \
  Rscript /scratch.global/balay011/ciseqtl_scripts/run_matrixeqtl.r $CIS_EQTL_OUTDIR/eur/eur_genotype_matrix.txt $CIS_EQTL_OUTDIR/eas/eas_genotype_matrix.txt $CIS_EQTL_OUTDIR/eur/eur_expression_matrix_vst.txt $CIS_EQTL_OUTDIR/eas/eas_expression_matrix_vst.txt /projects/standard/aventeic/balay011/references/reference_genome_anno/gencode.v45.primary_assembly.annotation.gtf