suppressPackageStartupMessages({
    library(data.table)
    library(DESeq2)})

args <- commandArgs(trailingOnly = TRUE)
eur_expr_counts  <- args[1]
eas_expr_counts  <- args[2]

eur_ciseqtl_dir <- dirname(eur_expr_counts)
eas_ciseqtl_dir <- dirname(eas_expr_counts)

apply_vst <- function(count_matrix) {

  count_matrix <- as.data.frame(count_matrix)
  rownames(count_matrix) <- count_matrix$geneid
  count_matrix$geneid <- NULL
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
eur.expr.common.vst <- apply_vst(eur.expr.common)
write.table(eur.expr.common.vst,
            file = file.path(eur_ciseqtl_dir, "eur_expression_matrix_vst.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

eas.expr.common <- fread(eas_expr_counts)
eas.expr.common.vst <- apply_vst(eas.expr.common)
write.table(eas.expr.common.vst,
            file = file.path(eas_ciseqtl_dir, "eas_expression_matrix_vst.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)