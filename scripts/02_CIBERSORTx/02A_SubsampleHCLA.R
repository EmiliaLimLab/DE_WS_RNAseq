# CIBERSORTx: Downsample HCLA data
# * scRNA-seq data used downloaded from https://github.com/LungCellAtlas/HLCA

out_dir <- "output/CIBERSORTx/"

#-------------------------------------------------------------------------------
# Load libraries
#-------------------------------------------------------------------------------
library(Seurat)

#-------------------------------------------------------------------------------
# Read files
#-------------------------------------------------------------------------------
hcla_seurat <- readRDS("/path/to/HCLA_core_dataset/b351804c-293e-4aeb-9c4c-043db67f4540.rds")

#-------------------------------------------------------------------------------
# Filter out non-epithelial cells, nasal cells, and data from active smokers
#-------------------------------------------------------------------------------
hcla_seurat_epithelial <- subset(hcla_data, 
                               ann_level_1 == "Epithelial" & 
                                 smoking_status != "active" & 
                                 tissue_level_2 != "inferior turbinate")

#-------------------------------------------------------------------------------
# Randomly subsample 1000 cells per cell type at ann_level_3
#-------------------------------------------------------------------------------
Idents(hcla_seurat_epithelial) <- hcla_data_epithelial$ann_level_3
hcla_seurat_epithelial_subset <- subset(hcla_data_epithelial, downsample=1000)

#-------------------------------------------------------------------------------
# Write HCLA scRNA-seq subset Seurat object to file
#-------------------------------------------------------------------------------
# Create ouptut directory if it doesn't exist
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

# Write to file
currentDate <- format(Sys.Date(), "%Y%m%d")
saveRDS(hcla_seurat_epithelial_subset, 
        paste0(out_dir, currentDate, "_HLCA_subset.rds"))


