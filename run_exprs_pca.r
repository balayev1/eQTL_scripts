suppressPackageStartupMessages({
  library(data.table)})

args <- commandArgs(trailingOnly = TRUE)
eur_in_vst <- args[1]
eas_in_vst <- args[2]

# set up output directories
eur_ciseqtl_dir <- dirname(eur_in_vst)
eas_ciseqtl_dir <- dirname(eas_in_vst)

apply_pca <- function(in_vst) {
    expr_dt <- fread(in_vst)
    expr <- t(expr_dt[,-1])

    prcompResult<-prcomp(expr,center=TRUE,scale.=TRUE)
    PCs<-prcompResult$x

    importanceTable <- summary(prcompResult)$importance
    PVEs <- importanceTable[2,]
    plot_filename <- sub(".txt$", "_pve_plot.pdf", basename(in_vst))

    # Save the plot
    plot_path <- file.path(dirname(in_vst), plot_filename)
    pdf(plot_path)
    plot(PVEs, 
         xlab = "PC Index", 
         ylab = "Proportion of Variance Explained (PVE)",
         main = paste("PCA PVE for", basename(in_vst)))
    dev.off()
    
    message(paste("PVE plot saved to:", plot_path))

    pca_dt <- data.table(
        sample_id = rownames(PCs),
        PCs,
        check.names = FALSE)
        
    return(pca_dt)
}

# EUR
eur_pca_dt <- apply_pca(eur_in_vst)
eur_filename <- "eur_expr_pcs.txt"
fwrite(eur_pca_dt, file = file.path(eur_ciseqtl_dir, eur_filename), sep = "\t", quote = FALSE)

# EAS
eas_pca_dt <- apply_pca(eas_in_vst)
eas_filename <- "eas_expr_pcs.txt"
fwrite(eas_pca_dt, file = file.path(eas_ciseqtl_dir, eas_filename), sep = "\t", quote = FALSE)