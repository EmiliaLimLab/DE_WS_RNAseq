# QC: plot and investigate PCA output from nf-core/rnaseq

out_dir <- "output/PCA/"

#-------------------------------------------------------------------------------
# Load libraries
#-------------------------------------------------------------------------------
library(tidyverse)
library(DESeq2)

#-------------------------------------------------------------------------------
# Read files
#-------------------------------------------------------------------------------
# Bias corrected count data: length scaled gene counts
bias_corrected_counts_mergeFA <- read.delim("data/nfcore-rnaseq/GSC-2581_mergeFA_ref-hg38-no-alt_trimgalore_salmon-gene-counts-length-scaled.tsv", 
                                            check.names=FALSE) %>% 
  column_to_rownames(var="gene_id")

#-------------------------------------------------------------------------------
# PCA: recompute principal components on length scaled gene counts
#-------------------------------------------------------------------------------
# nf-core/rnaseq uses length scaled gene counts and DESeq2 to perform PCA
# https://github.com/nf-core/rnaseq/blob/4053b2ec173fc0cfdd13ff2b51cfaceb59f783cb/workflows/rnaseq/main.nf#L681
#
# PCA code extracted from https://github.com/nf-core/rnaseq/blob/master/bin/deseq2_qc.r 

# Create DESEQDataSet (DDS) object
# * Note: Length-scaled gene counts are not integers, hence need to round
# * `design=~1` creates intercept-only model, equivalent to setting `blind=TRUE` for transformation.
samples.vec     <- colnames(bias_corrected_counts_mergeFA %>% 
                              select(-gene_name))
coldata         <- data.frame(samples.vec, sample=samples.vec, row.names=1) 
dds_pca <- DESeqDataSetFromMatrix(countData=round(bias_corrected_counts_mergeFA %>% 
                                                    select(-gene_name)), 
                                  colData=coldata, design=~1)
dds_pca <- estimateSizeFactors(dds_pca)
vst_name <- "rlog"
rld <- rlog(dds_pca)
assay(dds_pca, vst_name) <- assay(rld)

# Pre-process 
plotPCA_vst <- function (object,  ntop = 500, assay=length(assays(object))) {
  rv         <- rowVars(assay(object, assay))
  select     <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca        <- prcomp(t(assay(object, assay)[select, ]), center=TRUE, scale=FALSE)
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  df         <- cbind( as.data.frame(colData(object)), pca$x)
  #Order points so extreme samples are more likely to get label
  ord        <- order(abs(rank(df$PC1)-median(df$PC1)), abs(rank(df$PC2)-median(df$PC2)))
  df         <- df[ord,]
  attr(df, "percentVar") <- data.frame(PC=seq(along=percentVar), percentVar=100*percentVar)
  return(df)
}

# Perform PCA 
ntop <- c(500, Inf)
for (n_top_var in ntop) {
  pca.data      <- plotPCA_vst(dds_pca, assay=vst_name, ntop=n_top_var)
  percentVar    <- round(attr(pca.data, "percentVar")$percentVar)
} # at end of loop, we'll be using the user-defined ntop if any, else all genes

# Extract PC1 vs PC2 values to file
pca.vals <- pca.data[,c("PC1","PC2")]
colnames(pca.vals) <- paste0(colnames(pca.vals), ": ", percentVar[1:2], '% variance')
pca.vals <- rownames_to_column(pca.vals, var="sample")

#-------------------------------------------------------------------------------
# Write PC1 and PC2 to file
#-------------------------------------------------------------------------------
# Create ouptut directory if it doesn't exist
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

# Write to file
separate_wider_delim(pca.vals, cols="sample", delim="-", 
                     names=c("Participant", "Exposure")) %>% 
  write.table(paste0(out_dir, "GSC-2581_mergeFA_PCA.tsv"), 
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

