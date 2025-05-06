# Generate plots in figure 2 for DE-WS RNA-seq paper

out_dir <- "output/"

#-------------------------------------------------------------------------------
# Load libraries
#-------------------------------------------------------------------------------
library(tidyverse)
library(colorspace)
library(ggrepel)
library(ggpp)
library(ggtext)
library(eulerr)
library(patchwork)

#-------------------------------------------------------------------------------
# Read files
#-------------------------------------------------------------------------------
# Sample PCA
deseq_pca_with_metadata <- read.delim("data/QC/GSC-2581_mergeFA_PCA.tsv", 
                                      check.names=FALSE) %>% 
  mutate(Participant=as.character(Participant))

# DEA results
deseq_results_mergeFA_DE <- read.delim("data/DEA/GSC-2581_mergeFA_deseq-inclBasal-DEvsFA.tsv")
deseq_results_mergeFA_WS <- read_delim("data/DEA/GSC-2581_mergeFA_deseq-inclBasal-WSvsFA.tsv")

# Select genes involved in oxidative stress response (PMID31984210)
# * specifically, those in Fig 1F, Nrf2 induced only
OS_genes <- c("GPX2", "GCLC", "GCLM", "NQO1", "HMOX1", 
              "G6PD", "GSTA1", "GSTM1", "TXN", "NFE2L2", 
              "GSTP1", "CYP1B1")

# Select genes related to antiviral response
# * identified based on search of significant DEGs against anti-viral DEG lists 
#   (PMID22398282)
antiviral_genes <- c("HSH2D", "IL1RL1", "ISG15", "IFI27", "IFI6", 
                     "RSAD2", "DDX60", "IFIT1", "OASL", "ST3GAL3", 
                     "PLOD2", "CMPK2", "DPYSL3", "BCYRN1", "CYP4F11", 
                     "FCGBP", "VPS13D", "BHLHE40", "OAS1", "IFI35", 
                     "RPSA", "FXYD3")

#-------------------------------------------------------------------------------
# Plot A: PCA of length-scaled gene counts across all samples
#-------------------------------------------------------------------------------
sample_pca_plot <- deseq_pca_with_metadata %>% 
  ggplot(aes(x=`PC1: 37% variance`, y=`PC2: 13% variance`, 
             shape=Exposure, 
             color=Participant, fill=Participant)) + 
  geom_point(size=4) + 
  scale_shape_manual(values=c("IC" = 10, "FA" = 12, 
                              "DE" = 23, "WS" = 24), 
                     limits=c("IC", "FA", "DE", "WS")) + 
  scale_fill_discrete_qualitative(palette="dark3") + 
  scale_color_discrete_qualitative(palette="dark3") + 
  labs(shape="Exposure", color="Participant") + 
  theme_bw() + 
  theme(text=element_text(size=14)) + 
  guides(fill="none")

#-------------------------------------------------------------------------------
# Plot B: Volcano plot of DE & WS vs. FA DEGs
#-------------------------------------------------------------------------------
# Function definition: retain genes padj < 0.05 and label for significant OS genes
# + any outliers (thresholds identified by eye - padj < 0.0005, abs(log2FC) > 2.5)
add_plot_labels <- function(deseq_results, padj_threshold) {
  significant_genes <- deseq_results %>% 
    filter(padj < padj_threshold) %>% 
    pull(gene_id)
  
  deseq_results %>% 
    mutate(significant_gene_label=ifelse(gene_id %in% significant_genes & gene_name %in% OS_genes |
                                           gene_id %in% significant_genes & gene_name %in% antiviral_genes,
                                         gene_name, ""),  
           fill_label=ifelse(gene_id %in% significant_genes, 
                             ifelse(log2FoldChange < 0, 
                                    "Downregulated", "Upregulated"), 
                             NA))
}

# Add labels and combine into one df
deseq_results_mergeFA_WS_labelled <- add_plot_labels(deseq_results_mergeFA_WS, 0.05) %>% 
  mutate(facet_label="WS vs. FA")
deseq_results_mergeFA_DE_labelled <- add_plot_labels(deseq_results_mergeFA_DE, 0.05) %>% 
  mutate(facet_label="DE vs. FA")
deseq_results_mergeFA_volcano_data <- rbind(deseq_results_mergeFA_DE_labelled, 
                                            deseq_results_mergeFA_WS_labelled) %>% 
  mutate(log10_pval = -log10(pvalue), 
         log10_padj = -log10(padj))


# Create volcano plots
# PDF (landscape): 12in x 6in
deseq_results_mergeFA_volcano <- ggplot(deseq_results_mergeFA_volcano_data, 
                                        aes(x=log2FoldChange, y=log10_padj, 
                                            fill=fill_label, label=significant_gene_label)) + 
  geom_point(shape=21, size=4, alpha=0.7, colour="black") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed") + 
  geom_vline(xintercept=c(-1,1), linetype="dashed") + 
  labs(x="log<sub>2</sub>FC", y="-log<sub>10</sub>*P*<sub>adj</sub>") + 
  scale_fill_manual(values=c("Upregulated" = "#C7828C",
                             "Downregulated" = "#3D5F75", 
                             "NA" = "#55525B"), 
                    limits=c("Downregulated", "Upregulated")) + 
  xlim(-6, 6) + 
  geom_label_repel(position = position_nudge_center(x = 0.4, y = 0.7,
                                                    center_x = 0, center_y = 0), 
                   label.size=NA, label.padding = 0.1, 
                   box.padding=0.5, max.overlaps=Inf, min.segment.length=0, 
                   fill=rgb(red = 1, green = 1, blue = 1, alpha = 0.75)) + 
  facet_grid(. ~ facet_label) + 
  theme_bw() + 
  theme(text=element_text(size=14), 
        axis.title=element_markdown(), 
        legend.title=element_blank(), 
        legend.position="top", 
        legend.text=element_markdown())

#-------------------------------------------------------------------------------
# Plot C: Venn diagram of DEGs betwen DE & WS
#-------------------------------------------------------------------------------
# Data wrangling: pull lists of DEGs
deseq_results_mergeFA_DEGs <- list(DE=filter(deseq_results_mergeFA_DE_labelled, !is.na(fill_label)) %>% 
                                     pull(gene_id),
                                   WS=filter(deseq_results_mergeFA_WS_labelled, !is.na(fill_label)) %>% 
                                     pull(gene_id))

# Create Venn diagram
deseq_results_mergeFA_venn <- plot(euler(deseq_results_mergeFA_DEGs), 
                                   quantities=TRUE, 
                                   as_ggplot=TRUE)

#-------------------------------------------------------------------------------
# Generate panel figure
#-------------------------------------------------------------------------------
# PDF (landscape): 7in x 13in
layout <- "
AABBBB
AABBBB
CCBBBB
"

wrap_plots(A=sample_pca_plot, B=deseq_results_mergeFA_volcano, 
           C=deseq_results_mergeFA_venn, design=layout) + 
  plot_annotation(tag_level = "A") & 
  theme(plot.tag = element_text(size = 16,
                                hjust = 0,
                                vjust = 1))

# Create ouptut directory if it doesn't exist
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

# Write to file
ggsave(paste0(out_dir, "Figure2.pdf"), plot=panel_figure, 
       width=13, height=7, units="in", device="pdf")

