
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
# options(future.globals.maxSize = 1e9)
# future.globals.maxSize を 2 GiB に増加（必要に応じてサイズを調整）
options(future.globals.maxSize = 5 * 1024^3)



# C1
setwd("/Volumes/Pegasus32R8/workspace/omics_data/scRNA-seq/20200403_Seuratv3_based2/Seurat/dim8_k4/subset_cluster4")
query<- readRDS("C1.rds")
query@meta.data$sub.cluster <- as.factor(query@meta.data$sub.cluster)
query@meta.data$sub.cluster <- factor(query@meta.data$sub.cluster, 
                                         levels = c("0","1","2_0","2_1","3","4_0","4_1","5", "6","7","8","9","10","11","12"))
Idents(query) <- query@meta.data$sub.cluster 
DefaultAssay(query)
# [1] "SCT"

query@meta.data$integration <- paste("C1", query@meta.data$sub.cluster, sep = "_")
query@meta.data$method <- "C1" 
query@meta.data$sub.cluster.Fluidigm <- query@meta.data$sub.cluster 


setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered")
NCC_merge <- readRDS("NCC_merge2_RPCA_velo_region.rds")
Idents(NCC_merge) <- NCC_merge@meta.data$region2
NCC_merge <- RunUMAP(NCC_merge, reduction = "integrated.rpca", dims = 1:50, reduction.name = "umap.rpca", return.model = T)
DimPlot(NCC_merge, reduction = "umap.rpca", group.by = "region2")
NCC_merge@meta.data$method <- "10x"

DefaultAssay(query) <- "RNA"
query <- NormalizeData(query)
query <- FindVariableFeatures(query, selection.method = "vst")
query <- ScaleData(query)
query <- RunPCA(query, features = VariableFeatures(query, layer = "counts"))

query.anchors <- FindTransferAnchors(
  reference = NCC_merge,
  query = query,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  reduction = "rpca" 
)

query.mapped <- MapQuery(
  anchorset = query.anchors,
  query = query,
  reference = NCC_merge,
  refdata = list(celltype = "integrated",
                 region = "region2"), 
  reference.reduction = "integrated.rpca",          
  reduction.model = "umap.rpca"               
)

DimPlot(query.mapped, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE)
DimPlot(query.mapped, reduction = "ref.umap", group.by = "sub.cluster.Fluidigm", label = TRUE)
query.mapped@meta.data$predicted.region <- as.factor(query.mapped@meta.data$predicted.region )
query.mapped@meta.data$predicted.region  <- factor(query.mapped@meta.data$predicted.region ,
                                                   levels = c("Pharyngeal","Transitional", "Cardiac_Cushion","Subvalvular",
                                                              "AP_septum","SMC_GA","SMC_DA", "SMC_CA", "Neural","Glial",
                                                              "Melanocyte","Neural_tube", "Cardiomyocyte"))
DimPlot(query.mapped, reduction = "ref.umap", group.by = "predicted.region", label = F)
ggsave("umap_C1_region.pdf", width = 7, height = 5)


setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/with_Fluidigm/Projection")
saveRDS(query.mapped, "query.mapped.rds")








