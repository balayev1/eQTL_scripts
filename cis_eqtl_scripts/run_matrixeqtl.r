### load R libraries
suppressPackageStartupMessages({
    library(data.table)
    library(MatrixEQTL)})

args <- commandArgs(trailingOnly = TRUE)
eur_in_geno <- args[1]
eas_in_geno <- args[2]
eur_in_vst <- args[3]
eas_in_vst <- args[4]
gtf_file <- args[5]

eur_in_geno <- "/projects/standard/venteicher_30050/balay011/ciseqtl_analysis/eur/eur_genotype_matrix.txt"
eas_in_geno <- "/projects/standard/venteicher_30050/balay011/ciseqtl_analysis/eas/eas_genotype_matrix.txt"
eur_in_vst <- "/projects/standard/venteicher_30050/balay011/ciseqtl_analysis/eur/eur_expression_matrix_vst.txt"
eas_in_vst <- "/projects/standard/venteicher_30050/balay011/ciseqtl_analysis/eas/eas_expression_matrix_vst.txt"
gtf_file <- "/projects/standard/aventeic/balay011/references/reference_genome_anno/gencode.v45.primary_assembly.annotation.gtf"

# EUR
geno.common <- fread(eur_in_geno, header = TRUE)
expr.common <- fread(eur_in_vst, header = TRUE)

### generate SNP and gene position files
snpinfo <- data.frame(snpid = geno.common$snpid)
split_fields <- do.call(rbind, strsplit(snpinfo$snpid, ":", fixed = TRUE))
snpinfo$chr <- paste0("chr", split_fields[,1])
snpinfo$pos <- split_fields[,2]
snpinfo <- snpinfo[, c("snpid", "chr", "pos")]

expr_genes <- expr.common$geneid
skip_n <- system(paste("grep -n -m 1 -v '^#'", gtf_file, "| cut -d: -f1", sep = " "), intern = TRUE)
skip_n <- as.integer(skip_n) - 1
gtf <- fread(
  gtf_file,
  sep = "\t",
  header = FALSE,
  skip = skip_n,
  col.names = c("chr","source","feature","start","end","score","strand","frame","attr"))
gtf_gene <- gtf[feature == "gene"]
gtf_gene[, gene_name := sub('.*gene_name "([^"]+)".*', '\\1', attr)]
gtf_gene[, gene_id   := sub('.*gene_id "([^"]+)".*',   '\\1', attr)]
gtf_gene[gene_name == attr | gene_name == "", gene_name := gene_id]
gene_map <- gtf_gene[gene_name %in% expr_genes,
                     .(geneid = gene_name,
                       chr,
                       left = start,
                       right = end)]
setkey(gene_map, geneid)
gene_map <- gene_map[.(expr_genes)]
gene_map_unique <- gene_map[!duplicated(geneid)]

### set MatrixEQTL parameters
base.dir = dirname(eur_in_vst)
useModel = modelLINEAR
SNP_file_name = eur_in_geno
expression_file_name = eur_in_vst
covariates_file_name = paste0(base.dir, "/eur_covariates_matrix.txt")
output_file_name = ""
output_file_name.cis = paste0(base.dir, "/eur_cis_eqtl_results.txt")
cisDist = 1e6
pvOutputThreshold <- 0
pvOutputThreshold.cis <- 1e-2
snpspos <- snpinfo
snpspos$pos <- as.integer(snpspos$pos)
genepos <- gene_map_unique
pvalue.hist <- "qqplot"

snps = SlicedData$new()
snps$fileDelimiter = "\t"      # the TAB character
snps$fileOmitCharacters = "NA" # denote missing values
snps$fileSkipRows = 1          # one row of column labels
snps$fileSkipColumns = 1       # one column of row labels
snps$fileSliceSize = 2000      # read file in pieces of 2,000 rows
snps$LoadFile(SNP_file_name)

gene = SlicedData$new()
gene$fileDelimiter = "\t"      # the TAB character
gene$fileOmitCharacters = "NA" # denote missing values
gene$fileSkipRows = 1          # one row of column labels
gene$fileSkipColumns = 1       # one column of row labels
gene$fileSliceSize = 2000      # read file in pieces of 2,000 rows
gene$LoadFile(expression_file_name)

cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"      # the TAB character
cvrt$fileOmitCharacters = "NA" # denote missing values
cvrt$fileSkipRows = 1          # one row of column labels
cvrt$fileSkipColumns = 1       # one column of row labels
cvrt$fileSliceSize = 2000      # read file in pieces of 2,000 rows
cvrt$LoadFile(covariates_file_name)

### the main Matrix eQTL function execution
eqtl_eur_res = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    output_file_name.cis = output_file_name.cis,
    pvOutputThreshold = pvOutputThreshold,
    pvOutputThreshold.cis = pvOutputThreshold.cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    useModel = useModel,
    errorCovariance = numeric(),
    verbose = TRUE,
    pvalue.hist = pvalue.hist,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)

### save outputs
saveRDS(eqtl_eur_res, file = file.path(base.dir, "eur_ciseqtl_results.rds"))


# EAS
geno.common <- fread(eas_in_geno)
expr.common <- fread(eas_in_vst)

### generate SNP and gene position files
snpinfo <- data.frame(snpid = geno.common$snpid)
split_fields <- do.call(rbind, strsplit(snpinfo$snpid, ":", fixed = TRUE))
snpinfo$chr <- paste0("chr", split_fields[,1])
snpinfo$pos <- split_fields[,2]
snpinfo <- snpinfo[, c("snpid", "chr", "pos")]

expr_genes <- expr.common$geneid
skip_n <- system(paste("grep -n -m 1 -v '^#'", gtf_file, "| cut -d: -f1", sep = " "), intern = TRUE)
skip_n <- as.integer(skip_n) - 1
gtf <- fread(
  gtf_file,
  sep = "\t",
  header = FALSE,
  skip = skip_n,
  col.names = c("chr","source","feature","start","end","score","strand","frame","attr"))
gtf_gene <- gtf[feature == "gene"]
gtf_gene[, gene_name := sub('.*gene_name "([^"]+)".*', '\\1', attr)]
gtf_gene[, gene_id   := sub('.*gene_id "([^"]+)".*',   '\\1', attr)]
gtf_gene[gene_name == attr | gene_name == "", gene_name := gene_id]
gene_map <- gtf_gene[gene_name %in% expr_genes,
                     .(geneid = gene_name,
                       chr,
                       left = start,
                       right = end)]
setkey(gene_map, geneid)
gene_map <- gene_map[.(expr_genes)]
gene_map_unique <- gene_map[!duplicated(geneid)]

### set MatrixEQTL parameters
base.dir = dirname(eas_in_vst)
useModel = modelLINEAR
SNP_file_name = eas_in_geno
expression_file_name = eas_in_vst
covariates_file_name = paste0(base.dir, "/eas_covariates_matrix.txt")
output_file_name = ""
output_file_name.cis = paste0(base.dir, "/eas_cis_eqtl_results.txt")
cisDist = 1e6
pvOutputThreshold <- 0
pvOutputThreshold.cis <- 1e-2
snpspos <- snpinfo
snpspos$pos <- as.integer(snpspos$pos)
genepos <- gene_map_unique
pvalue.hist <- "qqplot"

snps = SlicedData$new()
snps$fileDelimiter = "\t"      # the TAB character
snps$fileOmitCharacters = "NA" # denote missing values
snps$fileSkipRows = 1          # one row of column labels
snps$fileSkipColumns = 1       # one column of row labels
snps$fileSliceSize = 2000      # read file in pieces of 2,000 rows
snps$LoadFile(SNP_file_name)

gene = SlicedData$new()
gene$fileDelimiter = "\t"      # the TAB character
gene$fileOmitCharacters = "NA" # denote missing values
gene$fileSkipRows = 1          # one row of column labels
gene$fileSkipColumns = 1       # one column of row labels
gene$fileSliceSize = 2000      # read file in pieces of 2,000 rows
gene$LoadFile(expression_file_name)

cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"      # the TAB character
cvrt$fileOmitCharacters = "NA" # denote missing values
cvrt$fileSkipRows = 1          # one row of column labels
cvrt$fileSkipColumns = 1       # one column of row labels
cvrt$fileSliceSize = 2000      # read file in pieces of 2,000 rows
cvrt$LoadFile(covariates_file_name)

### the main Matrix eQTL function execution
eqtl_eas_res = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    output_file_name.cis = output_file_name.cis,
    pvOutputThreshold = pvOutputThreshold,
    pvOutputThreshold.cis = pvOutputThreshold.cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    useModel = useModel,
    errorCovariance = numeric(),
    verbose = TRUE,
    pvalue.hist = pvalue.hist,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)

### save outputs
saveRDS(eqtl_eas_res, file = file.path(base.dir, "eas_ciseqtl_results.rds"))