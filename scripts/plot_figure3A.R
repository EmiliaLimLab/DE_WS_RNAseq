# Generate plot in figure 3A for DE-WS RNA-seq paper

out_dir <- "output/"

#-------------------------------------------------------------------------------
# Read files
#-------------------------------------------------------------------------------
# GSEA results
deseq_results_mergeFA_DE_gseKEGG <- read.delim("data/GSEA/GSC-2581_mergeFA_deseq-inclBasal-DEvsFA_GSEA-KEGG.tsv")
deseq_results_mergeFA_WS_gseKEGG <- read.delim("data/GSEA/GSC-2581_mergeFA_deseq-inclBasal-WSvsFA_GSEA-KEGG.tsv")

#-------------------------------------------------------------------------------
# Plot A: KEGG pathway enrichment results
#-------------------------------------------------------------------------------
# Data wrangling: calculate gene ratio
deseq_results_mergeFA_DE_gseKEGG <- mutate(deseq_results_mergeFA_DE_gseKEGG, 
                                           gene_ratio = (str_count("/") + 1) / setSize, 
                                           direction = ifelse(NES < 0, "suppressed", "activated"), 
                                           condition="DE vs. FA")

deseq_results_mergeFA_WS_gseKEGG <- mutate(deseq_results_mergeFA_WS_gseKEGG, 
                                           gene_ratio = (str_count("/") + 1) / setSize, 
                                           direction = ifelse(NES < 0, "suppressed", "activated"), 
                                           condition="WS vs. FA")


# Generate plot
gseKEGG_dotplot <- rbind(deseq_results_mergeFA_DE_gseKEGG, 
                         deseq_results_mergeFA_WS_gseKEGG) %>% 
  group_by(condition) %>% 
  slice_head(n=10) %>% 
  ggplot(aes(x=gene_ratio, y=str_wrap(Description, 30), fill=p.adjust, size=setSize)) + 
  geom_point(pch=21, color='black') + 
  facet_grid(condition ~ direction, scales="free_y", space="free_y") + 
  scale_x_continuous(limits=c(0, 0.15), breaks=c(0, 0.05, 0.10, 0.15)) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0.05) + 
  labs(x="Gene ratio", y="Description", fill="*P*<sub>adj</sub>", size="Set size") + 
  theme_bw() + 
  theme(text=element_text(size=14), 
        axis.text.y=element_text(lineheight = 0.6), 
        legend.title=element_markdown(), 
        panel.spacing.x = unit(0.4, "cm", data = NULL))

#-------------------------------------------------------------------------------
# Generate panel figure
#-------------------------------------------------------------------------------
# Generate plots
layout <- "
AABBB
CCCCC
CCCCC
"

panel_figure <- wrap_plots(A=gseKEGG_dotplot, B=plot_spacer(), 
           C=plot_spacer(), design=layout) + 
  plot_annotation(tag_level = "A")

# Create ouptut directory if it doesn't exist
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

# Write to file
ggsave(paste0(out_dir, "Figure3.pdf"), plot=panel_figure, 
       width=11, height=11, units="in", device="pdf")

