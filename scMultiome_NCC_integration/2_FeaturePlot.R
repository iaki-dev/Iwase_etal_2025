

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


setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered")
NCC_merge <- readRDS("NCC_merge2_RPCA_velo_region.rds")
DimPlot(NCC_merge, reduction = "umap.rpca")
DefaultAssay(NCC_merge)
dim(NCC_merge)
# [1] 33282 67150

# Featureplot
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/FeaturePlot")
DefaultAssay(NCC_merge) <- "RNA"

goi <- c("Wnt1","Sox2","Dlx1","Dlx3","Dlx5","Ebf2","Shox2","Sox9","Scx","Osr1","Myh11","Acta2","Plp1","Phox2a","Phox2b","Dct")
goi <- "Rgs5"
for(i in 1:length(goi)){
  g1 <- FeaturePlot(NCC_merge, 
                    features = goi[i], reduction = "umap.rpca", pt.size = 0.001, order = T)
  g1
  ggsave(file = paste(goi[i], ".pdf", sep = ""), plot = g1, width = 4, height = 4)  
}


