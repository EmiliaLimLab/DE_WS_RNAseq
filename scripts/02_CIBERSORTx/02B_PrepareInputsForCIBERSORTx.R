# CIBERSORTx: Prepare inputs for CIBERSORTx - CPM and scRNA-seq reference matrices

out_dir <- "output/CIBERSORTx/"

#-------------------------------------------------------------------------------
# Load libraries
#-------------------------------------------------------------------------------
library(biomaRt)
library(Seurat)
library(tidyverse)
library(edgeR)

#-------------------------------------------------------------------------------
# Read files
#-------------------------------------------------------------------------------
# Gene counts
gene_counts_mergeFA <- read.delim("https://www.bcgsc.ca/downloads/lim_lab/DE_WS_RNAseq/GSC-2581_mergeFA_ref-hg38-no-alt_trimgalore_salmon-gene-counts.tsv", 
                                  check.names=FALSE)

# Seurat object: subsampled HCLA dataset
hcla_seurat_epithelial_subset <- readRDS(gzcon(url("https://www.bcgsc.ca/downloads/lim_lab/DE_WS_RNAseq/20241118_HLCA_subset.rds")))

#-------------------------------------------------------------------------------
# Prepare CPM matrix
# * transform gene count matrix to CPM matrix
# * keep genes of interest only
#-------------------------------------------------------------------------------
# Transform gene counts to CPM
cpm_matrix_mergeFA <- filter(gene_counts_mergeFA, !is.na(gene_name)) %>% 
  column_to_rownames(var="gene_id") %>% 
  dplyr::select(-gene_name) %>% 
  as.matrix() %>% 
  cpm() %>% 
  as.data.frame()

# Drop gene if CPM is missing in 80% or greater of the participants
X <- 0  # minimum CPM
Y <- ncol(cpm_matrix_mergeFA)  # number of participants
keep <- rowSums(cpm_matrix_mergeFA > X) >= Y * 0.2
cpm_matrix_mergeFA_gt0.2subj <- cpm_matrix_mergeFA[keep, ]

# Ensure ENSEMBL gene IDs are properly formatted
cpm_matrix_mergeFA_gt0.2subj <- rownames_to_column(cpm_matrix_mergeFA_gt0.2subj, var="ensembl_gene_id") %>%
  mutate(ensembl_gene_id=sub('\\.[0-9]*$', '', ensembl_gene_id))

# Get corresponding gene IDs from BioMart
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=109)
BM_geneIDs <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "chromosome_name"), 
                    filters="ensembl_gene_id", 
                    values=pull(cpm_matrix_mergeFA_gt0.2subj, ensembl_gene_id), 
                    mart=ensembl)

# Identify genes of interest - exclude any:
# * mitochondrial genes and ENSEMBL gene IDs without corresponding HGNC symbol
# * Duplicated HGNC symbols
BM_geneIDs_filtered <- filter(BM_geneIDs, 
                              ! chromosome_name %in% c("MT", "GL000194.1") 
                              & hgnc_symbol != "") %>% 
  distinct(hgnc_symbol, .keep_all = TRUE)

# Keep genes of interest only in CPM matrix
cpm_matrix_mergeFA_gt0.2subj_filtered <- filter(cpm_matrix_mergeFA_gt0.2subj, 
                                               ensembl_gene_id %in% pull(BM_geneIDs_filtered, ensembl_gene_id))

#-------------------------------------------------------------------------------
# Prepare scRNA-seq data (used to create signature matrix)
# * only include genes that are detected in the cultured epithelium in our 
#   experiment
#-------------------------------------------------------------------------------
# Extract genes present in our experiment
genes_of_interest <- pull(cpm_matrix_mergeFA_gt0.2subj_filtered, ensembl_gene_id)

# Extract normalized data from HCLA scRNA-seq subset Seurat object and filter 
# to only include genes of interest
Idents(hcla_seurat_epithelial_subset) <- hcla_seurat_epithelial_subset$ann_level_3 # make sure that level 3 is set as default
hcla_epithelial_subset_norm <- as.data.frame(as.matrix(hcla_seurat_epithelial_subset@assays$RNA@data))
hcla_epithelial_subset_norm_filtered <- filter(hcla_epithelial_subset_norm, 
       rownames(hcla_epithelial_subset_norm) %in% genes_of_interest)

#-------------------------------------------------------------------------------
# Write CIBERSORTx inputs to file
#-------------------------------------------------------------------------------
# Create ouptut directory if it doesn't exist
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

# Function definition: write CIBERSORTx scRNA-seq reference matrix to file
write_cibersort_scRef <- function(seurat_obj, norm_data, out_fn) {
  # Write column names first
  write.table(t(c('ensembl_gene_id', as.character(seurat_obj@active.ident))), 
	      out_fn, 
	      row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

  # Write normalized values
  write.table(norm_data, out_fn, 
	      row.names=TRUE, col.names=FALSE, sep="\t", quote=FALSE, append=TRUE)

}

# Write to file
write.table(cpm_matrix_mergeFA_gt0.2subj_filtered, 
            paste0(out_dir, "GSC-2581_mergeFA_CIBERSORTx_geneMatrix.txt"), 
              row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

write_cibersort_scRef(hcla_seurat_epithelial_subset, hcla_epithelial_subset_norm_filtered, 
		      paste0(out_dir, "CIBERSORTx_scRef_HLCA_Epi.txt"))


