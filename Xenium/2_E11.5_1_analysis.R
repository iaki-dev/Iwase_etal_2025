# Xenium

library(Seurat)
library(future)
plan("multisession", workers = 10)
library(ggplot2)
library(tidyverse)
library(spacexr)
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
library(Nebulosa)
library(BPCells)
library(data.table)
options(future.globals.maxSize = 2e9)
options(Seurat.object.assay.version = "v5")


setwd("/Users/iwaseakiyasu/workspace/data/Xenium_data/analysis/Wnt1_Cre_E115_1")
path <- "/Users/iwaseakiyasu/workspace/data/Xenium_data/original_data/output-XETG00053__0016856__Wnt1_Cre_E115_1__20240124__080118"
# Load the Xenium data
xenium.obj <- LoadXenium(path, fov = "fov")
# remove cells with 0 counts
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)
dim(xenium.obj)
# [1]   334 77263
VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ggsave("0_QC.pdf", width = 5, height = 4)

xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")
xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))
xenium.obj <- RunUMAP(xenium.obj, dims = 1:30)
xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:30)
xenium.obj <- FindClusters(xenium.obj, resolution = 0.8) 
DimPlot(xenium.obj, label = T)
ggsave("2_umap_clusters.pdf", width = 5, height = 5)

FeaturePlot(xenium.obj, features = c("Pecam1", "Tcf24", "Csf1r"))

# segmaentation
DefaultBoundary(xenium.obj[["fov"]]) <- "segmentation"
ImageDimPlot(xenium.obj, fov = "fov", axes = F, border.color = "white", border.size = 0.01, cols = "polychrome",
             coord.fixed = FALSE) +
  scale_y_reverse()
# molecules = c("Sox9", "Myh11", "Pecam1", "Csf1r"), nmols = 10000,  mols.size = 0.01)
ggsave("1_ImageDim_clusters.pdf", width = 7, height = 5)

# save
# saveRDS(xenium.obj, "Wnt1_E11.5_1.rds")

setwd("/Volumes/Pegasus32R8/workspace/omics_data/Xenium/analysis/Wnt1_Cre_E115_1")
xenium.obj <- readRDS("Wnt1_E11.5_1.rds")
xenium.obj@meta.data$cell_id <- rownames(xenium.obj@meta.data)

meta <- read_csv("meta_celltype.full.csv")
colnames(meta)[2] <- "celltype_full"


setwd("/Volumes/Pegasus32R8/workspace/omics_data/Xenium/analysis/RCTD_visualization/")
df <- read.table("df.txt", header = T)
tmp.df <- df %>% dplyr::filter(orig.ident %in%  c("E11.5_NCC", "E11.5_nonNCC")) 
tmp.df$x %>% min() %>% round()
# Caution x and y
cropped.coords <- Crop(xenium.obj[["fov"]], y = c(tmp.df$x %>% min() %>% round(),
                                                  tmp.df$x %>% max() %>% round()),
                                                  x = c(tmp.df$y %>% min() %>% round(),
                                                        tmp.df$y %>% max() %>% round()), coords = "plot")

DefaultAssay(xenium.obj) <- "Xenium"
xenium.obj[["zoom"]] <- cropped.coords
DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentation"


xenium.obj@meta.data <- left_join(xenium.obj@meta.data , meta, by = "cell_id")

Idents(xenium.obj) <- "celltype_full"

setwd("/Volumes/Pegasus32R8/workspace/omics_data/Xenium/analysis/Wnt1_Cre_E115_1/ImageDimPlot")

# pharynegal
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.001, cols = "polychrome",
             coord.fixed = T, molecules = c("Dkk2", "Efemp1", "Foxd1", "Osr1", "Ebf2"), nmols = 20000, mols.size = 0.1, mols.alpha = 1)
ggsave("E11.5_pharyngeal.pdf", height = 5, width = 7)
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.001, cols = "polychrome",
             coord.fixed = T, molecules = c("Dkk2", "Efemp1", "Foxd1", "Osr1", "Ebf2"), nmols = 20000, mols.size = 0.1, mols.alpha = 1) +
theme_void() +
  theme(
    panel.background = element_rect(fill = "black"),  
    plot.background = element_rect(fill = "black")
  ) +
  NoLegend()
ggsave("E11.5_pharyngeal_NoLegend.pdf", height = 5, width = 7)


# OFT
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.001, cols = "polychrome",
             coord.fixed = T, molecules = c("Tbx20","Igfbp7", "Sox9", "Sall3"), nmols = 20000, mols.size = 0.001, mols.alpha = 1) 
ggsave("E11.5_OFT.pdf", height = 5, width = 7)
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.001, cols = "polychrome",
             coord.fixed = T, molecules = c("Tbx20","Igfbp7", "Sox9", "Sall3"), nmols = 20000, mols.size = 0.001, mols.alpha = 1) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black"),  
    plot.background = element_rect(fill = "black")
  ) +
  NoLegend()
ggsave("E11.5_OFT_NoLegend.pdf", height = 5, width = 7)


# SMC
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.001, cols = "polychrome",
             coord.fixed = T, molecules = c("Rgs5","Acta2","Eln"), nmols = 20000, mols.size = 0.001, mols.alpha = 1) 
ggsave("E11.5_SMC.pdf", height = 5, width = 7)
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.001, cols = "polychrome",
             coord.fixed = T, molecules = c("Rgs5","Acta2","Eln"), nmols = 20000, mols.size = 0.001, mols.alpha = 1) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black"),  
    plot.background = element_rect(fill = "black")
  ) +
  NoLegend()
ggsave("E11.5_SMC_NoLegend.pdf", height = 5, width = 7)


# Tcf24
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.001, cols = "polychrome",
             coord.fixed = T, molecules = c("Tcf24"), nmols = 20000, mols.size = 0.001, mols.alpha = 1) 
ggsave("E11.5_Tcf24.pdf", height = 5, width = 7)
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.001, cols = "polychrome",
             coord.fixed = T, molecules = c("Tcf24"), nmols = 20000, mols.size = 0.001, mols.alpha = 1) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black"),  
    plot.background = element_rect(fill = "black")
  ) +
  NoLegend()
ggsave("E11.5_Tcf24_NoLegend.pdf", height = 5, width = 7)
xenium.obj@meta.data %>% head


# Focus
cropped.coords <- Crop(xenium.obj[["fov"]], y = c(1100,1600),
                       x = c(1800,2400), coords = "plot")

DefaultAssay(xenium.obj) <- "Xenium"
xenium.obj[["zoom"]] <- cropped.coords
DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentation"

setwd("/Volumes/Pegasus32R8/workspace/omics_data/Xenium/analysis/Wnt1_Cre_E115_1/ImageDimPlot")
# pharynegal
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = T, molecules = c("Dkk2", "Efemp1", "Foxd1", "Osr1", "Ebf2"), nmols = 20000, mols.size = 0.1, mols.alpha = 1)
ggsave("E11.5_pharyngeal_zoom.pdf", height = 5, width = 7)
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = T, molecules = c("Dkk2", "Efemp1", "Foxd1", "Osr1", "Ebf2"), nmols = 20000, mols.size = 0.1, mols.alpha = 1) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black"), 
    plot.background = element_rect(fill = "black")
  ) +
  NoLegend()
ggsave("E11.5_pharyngeal_zoom_NoLegend.pdf", height = 5, width = 7)


# OFT
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = T, molecules = c("Tbx20","Igfbp7", "Sox9", "Sall3"), nmols = 20000, mols.size = 0.001, mols.alpha = 1) 
ggsave("E11.5_OFT_zoom.pdf", height = 5, width = 7)
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = T, molecules = c("Tbx20","Igfbp7", "Sox9", "Sall3"), nmols = 20000, mols.size = 0.001, mols.alpha = 1) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black"),  
    plot.background = element_rect(fill = "black")
  ) +
  NoLegend()
ggsave("E11.5_OFT_zoom_NoLegend.pdf", height = 5, width = 7)


# SMC
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = T, molecules = c("Rgs5","Acta2","Eln"), nmols = 20000, mols.size = 0.001, mols.alpha = 1) 
ggsave("E11.5_zoom_SMC.pdf", height = 5, width = 7)
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = T, molecules = c("Rgs5","Acta2","Eln"), nmols = 20000, mols.size = 0.001, mols.alpha = 1) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black"),  
    plot.background = element_rect(fill = "black")
  ) +
  NoLegend()
ggsave("E11.5_SMC_zoom_NoLegend.pdf", height = 5, width = 7)

# Tcf24
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = T, molecules = c("Tcf24"), nmols = 20000, mols.size = 0.001, mols.alpha = 1) 
ggsave("E11.5_Tcf24_zoom.pdf", height = 5, width = 7)
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = T, molecules = c("Tcf24"), nmols = 20000, mols.size = 0.001, mols.alpha = 1) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black"),  
    plot.background = element_rect(fill = "black")
  ) +
  NoLegend()
ggsave("E11.5_Tcf24_zoom_NoLegend.pdf", height = 5, width = 7)



