# Featureplot
library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)
library(tidyverse)
library(GenomicRanges)
library(future)
library(DoubletFinder)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")


setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/v10/v10.2.1/data")
gexatac_merge <- readRDS("4_gexatac_merge_links.rds")


# Featureplot
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/v10/v10.2.1/Featureplot")
DefaultAssay(gexatac_merge) <- "RNA"

goi <- "Pitx2"
for(i in 1:length(goi)){
  g1 <- FeaturePlot(gexatac_merge, 
                    features = goi[i], reduction = "umap.rpca")
  g1
  ggsave(file = paste(goi[i], ".pdf", sep = ""), plot = g1, width = 4, height = 4)  
}



DefaultAssay(gexatac_merge) <- "ATAC"  
