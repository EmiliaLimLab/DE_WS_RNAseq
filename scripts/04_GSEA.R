# GSEA using clusterProfiler
# * identify enriched KEGG pathways

out_dir <- "output/GSEA/"

#-------------------------------------------------------------------------------
# Load libraries
#-------------------------------------------------------------------------------
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyverse)
library(clusterProfiler)

#-------------------------------------------------------------------------------
# Read files
#-------------------------------------------------------------------------------
# DEA results
deseq_results_mergeFA_DE <- read.delim("data/DEA/GSC-2581_mergeFA_deseq-inclBasal-DEvsFA.tsv")
deseq_results_mergeFA_WS <- read.delim("data/DEA/GSC-2581_mergeFA_deseq-inclBasal-WSvsFA.tsv")

#-------------------------------------------------------------------------------
# Gene set enrichment analysis (GSEA) using gene lists ranked by DESeq2 log2FC
# * performed on WS and DE vs. FA only
#-------------------------------------------------------------------------------
# Create gene lists inputs and sort in decreasing order
deseq_results_mergeFA_WS_gene_list <- pull(deseq_results_mergeFA_WS, log2FoldChange)
names(deseq_results_mergeFA_WS_gene_list) <- pull(deseq_results_mergeFA_WS, gene_id)
deseq_results_mergeFA_WS_gene_list <- na.omit(deseq_results_mergeFA_WS_gene_list)
deseq_results_mergeFA_WS_gene_list_sorted <- sort(deseq_results_mergeFA_WS_gene_list, decreasing = TRUE)

deseq_results_mergeFA_DE_gene_list <- pull(deseq_results_mergeFA_DE, log2FoldChange)
names(deseq_results_mergeFA_DE_gene_list) <- pull(deseq_results_mergeFA_DE, gene_id)
deseq_results_mergeFA_DE_gene_list <-na.omit(deseq_results_mergeFA_DE_gene_list)
deseq_results_mergeFA_DE_gene_list_sorted <- sort(deseq_results_mergeFA_DE_gene_list, decreasing = TRUE)

# Map Ensembl to Entrez IDs
# * note that there are many Ensembl IDs without a corresponding Entrez ID - NA entries removed
deseq_results_mergeFA_WS_gene_list_entrez<- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                                                  keys = names(deseq_results_mergeFA_WS_gene_list_sorted),
                                                                  column = "ENTREZID",
                                                                  keytype = "ENSEMBL",
                                                                  multiVals = "first") %>% 
  enframe() %>% 
  dplyr::filter(!is.na(value)) %>% 
  left_join(enframe(deseq_results_mergeFA_WS_gene_list_sorted), 
            by="name") %>% 
  dplyr::select(-name) %>% 
  deframe()
deseq_results_mergeFA_WS_gene_list_entrez_sorted <- sort(deseq_results_mergeFA_WS_gene_list_entrez, decreasing = TRUE)

deseq_results_mergeFA_DE_gene_list_entrez<- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                                                  keys = names(deseq_results_mergeFA_DE_gene_list_sorted),
                                                                  column = "ENTREZID",
                                                                  keytype = "ENSEMBL",
                                                                  multiVals = "first") %>% 
  enframe() %>% 
  dplyr::filter(!is.na(value)) %>% 
  left_join(enframe(deseq_results_mergeFA_DE_gene_list_sorted), 
            by="name") %>% 
  dplyr::select(-name) %>% 
  deframe()
deseq_results_mergeFA_DE_gene_list_entrez_sorted <- sort(deseq_results_mergeFA_DE_gene_list_entrez, decreasing = TRUE)

# Make the KEGG enrichment object (gene ontology)
# * calculate gene ratio (https://github.com/YuLab-SMU/clusterProfiler/issues/364)
mergeFA_WS_gseKEGG <- gseKEGG(geneList=deseq_results_mergeFA_WS_gene_list_entrez_sorted, 
                              organism = "hsa", 
                              keyType = "kegg", 
                              pvalueCutoff = 1, 
                              minGSSize = 3, 
                              maxGSSize = 800, 
                              verbose = TRUE, 
                              pAdjustMethod = "BH",
                              by = "fgsea") %>% 
  as.data.frame() %>% 
  mutate(gene_ratio = (str_count("/") + 1) / setSize, 
         direction = ifelse(NES < 0, "suppressed", "activated"))

mergeFA_DE_gseKEGG <- gseKEGG(geneList=deseq_results_mergeFA_DE_gene_list_entrez_sorted, 
                              organism = "hsa", 
                              keyType = "kegg", 
                              pvalueCutoff = 1, 
                              minGSSize = 3, 
                              maxGSSize = 800, 
                              verbose = TRUE, 
                              pAdjustMethod = "BH",
                              by = "fgsea") %>% 
  as.data.frame() %>% 
  mutate(gene_ratio = (str_count("/") + 1) / setSize, 
         direction = ifelse(NES < 0, "suppressed", "activated"))

#-------------------------------------------------------------------------------
# Write GSEA KEGG results to file
#-------------------------------------------------------------------------------
# Create ouptut directory if it doesn't exist
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

# Write to file
dplyr::select(mergeFA_WS_gseKEGG, -c(gene_ratio, direction)) %>% 
  write.table(paste0(out_dir, "GSC-2581_mergeFA_deseq-inclBasal-WSvsFA_GSEA-KEGG.tsv"), 
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
dplyr::select(mergeFA_DE_gseKEGG, -c(gene_ratio, direction)) %>% 
  write.table(paste0(out_dir, "GSC-2581_mergeFA_deseq-inclBasal-DEvsFA_GSEA-KEGG.tsv"), 
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

