suppressPackageStartupMessages({
    library(data.table)
    library(DESeq2)})

args <- commandArgs(trailingOnly = TRUE)
eur_expr_counts  <- args[1]
eas_expr_counts  <- args[2]
gtf_file <- args[3]

eur_ciseqtl_dir <- dirname(eur_expr_counts)
eas_ciseqtl_dir <- dirname(eas_expr_counts)

filter_genes <- function(count_matrix, gtf_file){
  # read GTF file
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
  gtf_gene[, gene_type := sub('.*gene_type "([^"]+)".*', '\\1', attr)]
  gtf_gene[gene_name == attr | gene_name == "", gene_name := gene_id]

  count_matrix <- as.data.frame(count_matrix)
  rownames(count_matrix) <- count_matrix$geneid
  count_matrix$geneid <- NULL

  keep_types <- c("lncRNA", "protein_coding", "IG_C_gene", "TR_C_gene")  
  keep_names <- gtf_gene$gene_name[gtf_gene$gene_type %in% keep_types]

  ensg_ids <- gtf_gene$gene_name[grep("ENSG", gtf_gene$gene_name)]
  keep_names <- keep_names[!(keep_names %in% ensg_ids)]

  count_matrix_filtered <- count_matrix[rownames(count_matrix) %in% keep_names, ]

  return(count_matrix_filtered)
}

apply_vst <- function(count_matrix) {

  dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                colData = data.frame(row.names = colnames(count_matrix),
                                                     condition = rep("A", ncol(count_matrix))),
                                design = ~1)

  keep <- rowSums(counts(dds) >= 10) >= round(0.1 * ncol(counts(dds)))
  dds <- dds[keep, ]

  vsd <- vst(dds, blind = TRUE)
  vst_mat <- assay(vsd)

  vst_df <- data.frame(
    geneid = rownames(vst_mat),
    vst_mat,
    check.names = FALSE,
    row.names = NULL)

  return(vst_df)
}

eur.expr.common <- fread(eur_expr_counts)
eur.expr.common.filtered <- filter_genes(eur.expr.common, gtf_file)
eur.expr.common.vst <- apply_vst(eur.expr.common.filtered)
write.table(eur.expr.common.vst,
            file = file.path(eur_ciseqtl_dir, "eur_expression_matrix_vst.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

eas.expr.common <- fread(eas_expr_counts)
eas.expr.common.filtered <- filter_genes(eas.expr.common, gtf_file)
eas.expr.common.vst <- apply_vst(eas.expr.common.filtered)
write.table(eas.expr.common.vst,
            file = file.path(eas_ciseqtl_dir, "eas_expression_matrix_vst.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)