# Differential expression analysis using DESeq2
# * Compare DE vs. FA (DE) and WS vs FA (WS)

out_dir <- "output/DEA/"

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

# CIBERSORTX output
cibersort_basal_cell_counts <- read.delim("data/CIBERSORTx/GSC-2581_mergeFA_cibersortx.tsv") %>% 
  select(Mixture, Basal)

#-------------------------------------------------------------------------------
# Data wrangling
#-------------------------------------------------------------------------------
# Generate metadata from cell counts matrix and add basal cell proportions, 
# as determined by CIBERSORTx
# * ensure same order as bias_corrected_counts_mergeFA
# * separate sample into participant, condition
metadata_mergeFA_withBasalCounts <- data.frame("sample"=select(bias_corrected_counts_mergeFA, -gene_name) %>% 
                                                 colnames()) %>% 
  left_join(cibersort_basal_cell_counts, by=join_by("sample" == "Mixture")) %>% 
  separate_wider_delim(cols=sample, delim="-", names=c("subject", "condition"), cols_remove=FALSE) %>% 
  column_to_rownames(var="sample")

# Set factor levels
metadata_mergeFA_withBasalCounts <- mutate(metadata_mergeFA_withBasalCounts,
                                           subject=factor(subject, 
                                                          levels = c("6", "1", "2", "3", "4", "5")))
  
#-------------------------------------------------------------------------------
# Differential expression analysis
#-------------------------------------------------------------------------------
# Function definition: create DESEQDataSet (DDS) object
# * Note: Length-scaled gene counts are not integers, hence need to round
create_dds <- function(counts_matrix, metadata) {
  dds <- DESeqDataSetFromMatrix(countData=round(counts_matrix), 
                                colData=metadata, 
                                ~ subject + Basal + condition)
  
  # Pre-filtering low count genes
  # * Recommended starting values from developer: X=10, Y=smallest_group_sample_size
  # * May want to adjust based on data
  X <- 10
  Y <- metadata %>% 
    group_by(condition) %>% 
    summarize(n=n()) %>% 
    pull(n) %>% 
    min()
  keep <- rowSums(counts(dds) >= X) >= Y
  dds_filtered_lc <- dds[keep, ]
  
  return(dds_filtered_lc)
}

# Run DESeq: 
# * DE and WS vs. FA
# * FA vs. IC
deseq_results <- DESeq(create_dds(select(bias_corrected_counts_mergeFA, -gene_name), metadata_mergeFA_withBasalCounts))
deseq_results_mergeFA_WS <- results(deseq_results, contrast=c("condition", "WS", "FA")) %>% 
  data.frame() %>% 
  arrange(padj)
deseq_results_mergeFA_DE <- results(deseq_results, contrast=c("condition", "DE", "FA")) %>% 
  data.frame() %>% 
  arrange(padj)
deseq_results_IC_mergeFA <- results(deseq_results, contrast=c("condition", "FA", "IC")) %>% 
  data.frame() %>% 
  arrange(padj)

#-------------------------------------------------------------------------------
# Write DESeq2 results to file
#-------------------------------------------------------------------------------
# Create ouptut directory if it doesn't exist
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

# Function definition: write DESeq2 results tables to file
write_deseq_results <- function(deseq_results, out_fn) {
  deseq_results %>% 
    rownames_to_column(var="gene_id") %>% 
    left_join(rownames_to_column(bias_corrected_counts_mergeFA, var="gene_id") %>% 
                select(gene_id, gene_name), 
              by="gene_id") %>% 
    relocate(gene_name, .after=gene_id) %>% 
    write.table(out_fn, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}

# Write tables
write_deseq_results(deseq_results_mergeFA_WS, paste0(out_dir, "GSC-2581_mergeFA_deseq-inclBasal-WSvsFA.tsv"))
write_deseq_results(deseq_results_mergeFA_DE, paste0(out_dir, "GSC-2581_mergeFA_deseq-inclBasal-DEvsFA.tsv"))
write_deseq_results(deseq_results_IC_mergeFA, paste0(out_dir, "GSC-2581_mergeFA_deseq-inclBasal-FAvsIC.tsv"))

