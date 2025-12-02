#!/usr/bin/env bash

############# extract PCs for genotype data

### request some resources
srun --mem-per-cpu=16GB --time=12:00:00 --pty --cpus-per-task=8 --x11 -p interactive bash

# Threads
THREADS=8

### Input tumor SNP and INDEL VCF files
WORK_DIR="/scratch.global/balay011/gwas_analysis"

SNP_EUR_TUMOR="${WORK_DIR}/cohort.eur.cases.snps.vcf.gz"
SNP_EAS_TUMOR="${WORK_DIR}/cohort.eas.cases.snps.vcf.gz"
INDEL_EUR_TUMOR="${WORK_DIR}/cohort.eur.cases.indels.vcf.gz"
INDEL_EAS_TUMOR="${WORK_DIR}/cohort.eas.cases.indels.vcf.gz"

### Tools/modules
module load bcftools/1.16-gcc-8.2.0-5d4xg4y
module load conda
source activate gwas_env

BCFTOOLS="bcftools"
PLINK2="plink2"
PLINK="plink"

### Output directory
mkdir -p "${WORK_DIR}"/{qc,tmp,eur,eas}
cd "$WORK_DIR"

### Other files
PSAM_TUMOR_EUR="$WORK_DIR/qc/tumors.eur.sex.psam"         # sex info for EUR cases
PSAM_TUMOR_EAS="$WORK_DIR/qc/tumors.eas.sex.psam"       # sex info for EAS cases
EUR_PHENO_TSV="$WORK_DIR/eur/pheno_eur.txt"   # FID\tIID\tTUMOR(1=control,2=case)
EAS_PHENO_TSV="$WORK_DIR/eas/pheno_eas.txt"   # FID\tIID\tTUMOR(1=control,2=case)

######################## QC thresholds
MIND=0.1      # sample missingness
GENO_SNP=0.05   # SNP missingness
GENO_INDEL=0.05 # INDEL missingness
MAF_CUTOFF=0.01 # MAF cutoff
HWE_CUTOFF=1e-6  # HWE cutoff
DUPLICATE_KING_FILE=/scratch.global/balay011/gwas_qc_struct/logs/king_cohort_duplicates.king.cutoff.out.id # KING duplicate sample IDs
REF_FA=$MSIPROJECT/balay011/references/reference_genome/GRCh38.primary_assembly.genome.fa

######################## Convert tumor VCFs to PLINK2 PGEN format keeping only autosome variants
$PLINK2 --threads $THREADS \
  --psam $PSAM_TUMOR_EUR \
  --pheno $EUR_PHENO_TSV --pheno-name TUMOR --keep $EUR_PHENO_TSV \
  --vcf $SNP_EUR_TUMOR \
  --autosome \
  --make-pgen \
  --out $WORK_DIR/eur/cohort1_eur_cases_snps

$PLINK2 --threads $THREADS \
  --psam $PSAM_TUMOR_EUR \
  --pheno $EUR_PHENO_TSV --pheno-name TUMOR --keep $EUR_PHENO_TSV \
  --vcf $INDEL_EUR_TUMOR \
  --autosome \
  --make-pgen \
  --out $WORK_DIR/eur/cohort1_eur_cases_indels

$PLINK2 --threads $THREADS \
  --psam $PSAM_TUMOR_EAS \
  --pheno $EAS_PHENO_TSV --pheno-name TUMOR --keep $EAS_PHENO_TSV \
  --vcf $SNP_EAS_TUMOR \
  --autosome \
  --make-pgen \
  --out $WORK_DIR/eas/cohort1_eas_cases_snps

$PLINK2 --threads $THREADS \
  --psam $PSAM_TUMOR_EAS \
  --pheno $EAS_PHENO_TSV --pheno-name TUMOR --keep $EAS_PHENO_TSV \
  --vcf $INDEL_EAS_TUMOR \
  --autosome \
  --make-pgen \
  --out $WORK_DIR/eas/cohort1_eas_cases_indels

# Remove chondrosarcoma sample, nHets/HomAlt outlier, low-quality RNA-Seq and duplicate samples
echo "Chord_7b_Blood_S54_sorted_MD_BQ Chord_7b_Blood_S54_sorted_MD_BQ" >> $DUPLICATE_KING_FILE
echo "Chord_55b_Blood_S50_sorted_MD_BQ Chord_55b_Blood_S50_sorted_MD_BQ" >> $DUPLICATE_KING_FILE
echo "Chord_20a_Blood_S26_sorted_MD_BQ Chord_20a_Blood_S26_sorted_MD_BQ" >> $DUPLICATE_KING_FILE
echo "Chord_67_Blood_S51_sorted_MD_BQ Chord_67_Blood_S51_sorted_MD_BQ" >> $DUPLICATE_KING_FILE
echo "Chord_22_Blood_S31_sorted_MD_BQ Chord_22_Blood_S31_sorted_MD_BQ" >> $DUPLICATE_KING_FILE
echo "Chord_135_Blood_S15_sorted_MD_BQ Chord_135_Blood_S15_sorted_MD_BQ" >> $DUPLICATE_KING_FILE

# Include the list of SRR accessions that have both RNA-Seq and WGS
SRR_ACCESSIONS=(
    "SRR14097754" "SRR14097756" "SRR14097853" "SRR14097761" "SRR14097831"
    "SRR14097760" "SRR14097851" "SRR14097711" "SRR14097856" "SRR14097858"
    "SRR14097838" "SRR14097860" "SRR14097697" "SRR14097863" "SRR14097865"
    "SRR14097867" "SRR14097869" "SRR14097871" "SRR14097699" "SRR14097701"
    "SRR14097703" "SRR14097705" "SRR14097707" "SRR14097709" "SRR14097713"
    "SRR14097715" "SRR14097717" "SRR14097813" "SRR14097719" "SRR14097721"
    "SRR14097723" "SRR14097725" "SRR14097814" "SRR14097727" "SRR14097815"
    "SRR14097729" "SRR14097816" "SRR14097738" "SRR14097818" "SRR14097820"
    "SRR14097742" "SRR14097744" "SRR14097746" "SRR14097748" "SRR14097821"
    "SRR14097752" "SRR14097823" "SRR14097828" "SRR14097829" "SRR14097830"
    "SRR14097835" "SRR14097844" "SRR14097847"
)

# Loop through the accessions, add the suffix, duplicate the name, and append to the file
for id in "${SRR_ACCESSIONS[@]}"; do
    SUFFIXED_ID="${id}_sorted_MD_BQ"
    echo "${SUFFIXED_ID} ${SUFFIXED_ID}" >> "$DUPLICATE_KING_FILE"
done

$PLINK2 --threads $THREADS \
  --pfile $WORK_DIR/eur/cohort1_eur_cases_snps \
  --remove $DUPLICATE_KING_FILE \
  --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 800 \
  --rm-dup force-first \
  --output-chr 26 \
  --make-bed \
  --out $WORK_DIR/eur/cohort1_eur_cases_snps_unrel
# 57 cases and 0 controls remaining after main filters.

$PLINK2 --threads $THREADS \
  --pfile $WORK_DIR/eur/cohort1_eur_cases_indels \
  --remove $DUPLICATE_KING_FILE \
  --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 800 \
  --rm-dup force-first \
  --output-chr 26 \
  --make-bed \
  --out $WORK_DIR/eur/cohort1_eur_cases_indels_unrel
# 57 cases and 0 controls remaining after main filters.

echo ">>> Implementing filters"
######################## Missingness filters (variant then sample)
# EUR TUMOR - merge & filter
${PLINK} --bfile $WORK_DIR/eur/cohort1_eur_cases_snps_unrel \
  --bmerge $WORK_DIR/eur/cohort1_eur_cases_indels_unrel \
  --make-bed \
  --out $WORK_DIR/eur/cohort1_eur_cases_merged
# 27619945 variants and 57 people pass filters and QC.

${PLINK2} --threads $THREADS \
--bfile $WORK_DIR/eur/cohort1_eur_cases_merged \
--geno $GENO_SNP \
--maf $MAF_CUTOFF \
--make-bed \
--out $WORK_DIR/eur/cohort2_eur_cases_merged
# 7065284 variants remaining after main filters.

${PLINK2} --threads $THREADS \
--bfile $WORK_DIR/eur/cohort2_eur_cases_merged \
--mind $MIND \
--make-bed \
--out $WORK_DIR/eur/cohort3_eur_cases_merged
# 0 samples removed due to missing genotype data (--mind).

echo "Number of SNPs:"
awk 'length($5)==1 && length($6)==1' $WORK_DIR/eur/cohort3_eur_cases_merged.bim | wc -l
# 6504979

echo "Number of INDELs:"
awk 'length($5)!=1 || length($6)!=1' $WORK_DIR/eur/cohort3_eur_cases_merged.bim | wc -l
# 560305

echo ">>> Extracting biallelic SNPs"
$PLINK2 --threads $THREADS \
  --bfile $WORK_DIR/eur/cohort3_eur_cases_merged \
  --snps-only just-acgt \
  --max-alleles 2 \
  --ref-from-fa $REF_FA \
  --make-bed \
  --out $WORK_DIR/eur/eur_cases_snps_4pca
# 6503945 variants remaining after main filters.

echo ">>> Running PCA"
$PLINK2 --bfile $WORK_DIR/eur/eur_cases_snps_4pca \
  --pca approx 20 allele-wts \
  --out $WORK_DIR/eur/eur_pca_results

cp $WORK_DIR/eur/eur_pca_results* /home/venteicher_30050/balay011/ciseqtl_analysis/eur

# EAS TUMOR - merge & filter
${PLINK2} --threads ${THREADS} \
  --pfile ${WORK_DIR}/eas/cohort1_eas_cases_snps \
  --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 800 \
  --remove $DUPLICATE_KING_FILE \
  --rm-dup force-first \
  --output-chr 26 \
  --make-bed \
  --out ${WORK_DIR}/eas/cohort1_eas_cases_snps_unrel
# 29 cases and 0 controls remaining after main filters.

${PLINK2} --threads ${THREADS} \
  --pfile ${WORK_DIR}/eas/cohort1_eas_cases_indels \
  --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 800 \
  --remove $DUPLICATE_KING_FILE \
  --rm-dup force-first \
  --output-chr 26 \
  --make-bed \
  --out ${WORK_DIR}/eas/cohort1_eas_cases_indels_unrel
# 29 cases and 0 controls remaining after main filters.

${PLINK} --bfile $WORK_DIR/eas/cohort1_eas_cases_snps_unrel \
  --bmerge $WORK_DIR/eas/cohort1_eas_cases_indels_unrel \
  --make-bed \
  --out $WORK_DIR/eas/cohort1_eas_cases_merged
# Total genotyping rate is 0.635864.
# 27619945 variants and 29 people pass filters and QC.

${PLINK2} --threads ${THREADS} \
  --bfile ${WORK_DIR}/eas/cohort1_eas_cases_merged \
  --geno ${GENO_SNP} \
  --maf ${MAF_CUTOFF} \
  --make-bed \
  --out ${WORK_DIR}/eas/cohort2_eas_cases_merged
# 7091534 variants remaining after main filters.

${PLINK2} --threads ${THREADS} \
  --bfile ${WORK_DIR}/eas/cohort2_eas_cases_merged \
  --mind ${MIND} \
  --make-bed \
  --out ${WORK_DIR}/eas/cohort3_eas_cases_merged
# 0 samples removed due to missing genotype data (--mind).

echo "Number of SNPs:"
awk 'length($5)==1 && length($6)==1' ${WORK_DIR}/eas/cohort3_eas_cases_merged.bim | wc -l
# 6492922

echo "Number of INDELs:"
awk 'length($5)!=1 || length($6)!=1' ${WORK_DIR}/eas/cohort3_eas_cases_merged.bim | wc -l
# 598612

echo ">>> Extracting biallelic SNPs"
$PLINK2 --threads $THREADS \
  --bfile ${WORK_DIR}/eas/cohort3_eas_cases_merged \
  --snps-only just-acgt \
  --max-alleles 2 \
  --ref-from-fa $REF_FA \
  --make-bed \
  --out $WORK_DIR/eas/eas_cases_snps_4pca
# 6490685 variants remaining after main filters.

$PLINK2 --bfile $WORK_DIR/eas/eas_cases_snps_4pca \
    --freq \
    --out $WORK_DIR/eas/eas_cases_snps_4pca.frq

echo ">>> Running PCA"
$PLINK2 --bfile $WORK_DIR/eas/eas_cases_snps_4pca \
    --read-freq $WORK_DIR/eas/eas_cases_snps_4pca.frq.afreq \
    --pca approx 20 allele-wts \
    --out $WORK_DIR/eas/eas_pca_results

cp $WORK_DIR/eas/eas_pca_results* /home/venteicher_30050/balay011/ciseqtl_analysis/eas/