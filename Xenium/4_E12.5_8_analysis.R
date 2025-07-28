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
options(future.globals.maxSize = 2e9)
options(Seurat.object.assay.version = "v5")

setwd("/Users/iwaseakiyasu/workspace/data/Xenium_data/analysis/Wnt1_Cre_E125_8")
path <- "/Users/iwaseakiyasu/workspace/data/Xenium_data/original_data/output-XETG00053__0016856__Wnt1_Cre_E125_8__20240124__080118"
# Load the Xenium data
xenium.obj <- LoadXenium(path, fov = "fov")
# remove cells with 0 counts
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)
dim(xenium.obj)
# [1]   334 103509
VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ggsave("0_QC.pdf", width = 5, height = 4)

xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")
xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))
xenium.obj <- RunUMAP(xenium.obj, dims = 1:30)
xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:30)
xenium.obj <- FindClusters(xenium.obj, resolution = 0.8) 
DimPlot(xenium.obj, label = T)
ggsave("2_umap_clusters.pdf", width = 5, height = 5)

FeaturePlot(xenium.obj, features = c("Pecam1", "Tcf24"))

# segmaentation
DefaultBoundary(xenium.obj[["fov"]]) <- "segmentation"
ImageDimPlot(xenium.obj, fov = "fov", axes = F, border.color = "white", border.size = 0.01, cols = "polychrome",
             coord.fixed = FALSE) +
  scale_y_reverse()
ggsave("1_ImageDim_clusters.pdf", width = 7, height = 5)

# save
# saveRDS(xenium.obj, "Wnt1_E12.5_8.rds")
setwd("/Users/iwaseakiyasu/workspace/data/Xenium_data/analysis/Wnt1_Cre_E125_8")
xenium.obj <- readRDS("Wnt1_E12.5_8.rds")


setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/v11/data")
gexatac_merge <- readRDS("gexatac_merge_RPCA.rds")
DimPlot(gexatac_merge, reduction = "umap.rpca", group.by = "seurat_clusters", label = T)


setwd("/Users/iwaseakiyasu/workspace/data/Xenium_data/analysis/merge")
xenium_merge <- readRDS("xenium_merge.full.rds")

coi <- xenium_merge@meta.data %>% dplyr::filter(section == "E12.5_8")
xenium_merge.sub <- xenium_merge[, rownames(coi)]
xenium_merge.sub@meta.data %>% rownames() %>% head
# [1] "aaaaeako-1_E12.5_8" "aaabhimh-1_E12.5_8" "aaacbago-1_E12.5_8" "aaadcjee-1_E12.5_8" "aaaelkoh-1_E12.5_8" "aaaemdkm-1_E12.5_8"

xenium.obj@meta.data %>% rownames() %>% head
# [1] "aaaaeako-1" "aaabhimh-1" "aaacbago-1" "aaadcjee-1" "aaaelkoh-1" "aaaemdkm-1"

tmp.name <- xenium_merge.sub@meta.data %>% rownames() %>% str_split(pattern = "_", simplify = T)

tmp.name[,1] %>% head

xenium.obj@meta.data$celltype.full <- xenium_merge.sub@meta.data$celltype.full
xenium.obj@meta.data$celltype.full <- factor(xenium.obj@meta.data$celltype.full, 
                                             levels = c(0:38))
xenium.obj@meta.data$celltype.full
DimPlot(xenium.obj, group.by = "celltype.full", label = T)
ggsave("3_umap_celltype.full.pdf", width = 5, height = 5)


ImageDimPlot(xenium.obj, fov = "fov", axes = F, border.color = "white", border.size = 0.01, cols = "polychrome",
             coord.fixed = FALSE, group.by = "celltype.full") +
  scale_y_reverse()
ggsave("3_ImageDim_celltype.full.pdf", width = 7, height = 5)


meta <- xenium.obj@meta.data
write.table(meta, "meta.txt", quote = F, sep = "\t", col.names = NA)



# RTCD cell type decomposition
query.counts <- GetAssayData(xenium.obj, assay = "Xenium", layer = "counts")[, Cells(xenium.obj[["fov"]])]
coords <- GetTissueCoordinates(xenium.obj[["fov"]], which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))

# E12.5(scMultiome) NCC
coi <- gexatac_merge@meta.data %>% dplyr::filter(sample %in% c("E12.5_NCC"))
E12.5 <- gexatac_merge[, rownames(coi)]
E12.5$seurat_clusters <- E12.5$seurat_clusters %>% as.character()

table.cluster_lineage <- E12.5@meta.data$seurat_clusters %>% table %>% as.data.frame()
table.cluster_lineage <- table.cluster_lineage %>% dplyr::filter(Freq > 5)
using.cluster_lineages <- table.cluster_lineage$. %>% as.character()
using.cluster_lineages

coi <- E12.5@meta.data %>% dplyr::filter(seurat_clusters %in% using.cluster_lineages)
counts <- GetAssayData(E12.5[, rownames(coi)], layer = "RNA", slot = "counts")

counts@x <- round(counts@x)

cluster <- as.factor(E12.5[, rownames(coi)]$seurat_clusters)
names(cluster) <- colnames(E12.5[, rownames(coi)])
nUMI <- E12.5[, rownames(coi)]$nCount_RNA
names(nUMI) <- colnames(E12.5[, rownames(coi)])
nUMI <- colSums(counts)
reference <- Reference(counts, cluster, nUMI)

# run RCTD with many cores
RCTD <- create.RCTD(query, reference, max_cores = 8,
                    CELL_MIN_INSTANCE = 1)

# doublet mode, which assigns 1-2 cell types per spot and is recommended for technologies with
# high spatial resolution such as Slide-seq and MERFISH
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
setwd("/Users/iwaseakiyasu/workspace/data/Xenium_data/analysis/Wnt1_Cre_E125_8/RCTD/NCC")
# saveRDS(RCTD, "RCTD_E12.5_celltype_double_NCC.rds")

RCTD_NCC <- readRDS("RCTD_E12.5_celltype_double_NCC.rds")

annotations.df <- RCTD_NCC@results$results_df
annotations.df %>% ggplot(aes(x=singlet_score)) +
  geom_histogram()
ggsave("hist_siglet_score.pdf", height = 4, width = 5)
annotations.df$singlet_score %>% summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.376  25.401  37.707  45.547  67.880 118.723 

df_weights <- RCTD_NCC@results$weights %>% as.data.frame()
df_weights <- cbind(rownames(df_weights) ,df_weights)
colnames(df_weights)[1] <- "cell_id"
df_weights %>% dplyr::filter(cell_id == "khhhkaob-1") %>% t() %>% as.data.frame() %>% View

# first_typeだとrejectされたものも含まれるので注意して変更する
annotations <- annotations.df$first_type
names(annotations) <- rownames(annotations.df)
xenium.obj$predicted.celltype <- annotations
keep.cells <- Cells(xenium.obj)[!is.na(xenium.obj$predicted.celltype)]
xenium.obj <- subset(xenium.obj, cells = keep.cells)

ImageDimPlot(xenium.obj, fov =  "fov", group.by = "predicted.celltype", size = 1.5, border.size = 0.01, cols = "polychrome", dark.background = T) +
  ggtitle("Cell type") +
  scale_y_reverse()
ggsave("4_decomp_E12.5_NCC.pdf", width = 7, height = 5)
write.table(xenium.obj@meta.data, "4_meta.predicted.celltype.txt", quote = F, sep = "\t", col.names = NA)




setwd("/Volumes/Pegasus32R8/workspace/omics_data/Xenium/analysis/Wnt1_Cre_E125_8")
xenium.obj <- readRDS("Wnt1_E12.5_8.rds")
xenium.obj@meta.data$cell_id <- rownames(xenium.obj@meta.data)

meta <- read_csv("meta_celltype.full.csv")
colnames(meta)[2] <- "celltype_full"


setwd("/Volumes/Pegasus32R8/workspace/omics_data/Xenium/analysis/RCTD_visualization/")
df <- read.table("df.txt", header = T)
tmp.df <- df %>% dplyr::filter(orig.ident %in%  c("E12.5_NCC", "E12.5_nonNCC")) 
tmp.df$x %>% min() %>% round()
# Seuratでx, yが逆になるので注意
cropped.coords <- Crop(xenium.obj[["fov"]], y = c(tmp.df$x %>% min() %>% round(),
                                                  tmp.df$x %>% max() %>% round()),
                       x = c(tmp.df$y %>% min() %>% round(),
                             tmp.df$y %>% max() %>% round()), coords = "plot")

DefaultAssay(xenium.obj) <- "Xenium"
xenium.obj[["zoom"]] <- cropped.coords
DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentation"


xenium.obj@meta.data <- left_join(xenium.obj@meta.data , meta, by = "cell_id")

Idents(xenium.obj) <- "celltype_full"

setwd("/Volumes/Pegasus32R8/workspace/omics_data/Xenium/analysis/Wnt1_Cre_E125_8/ImageDimPlot")

# pharynegal
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = T, molecules = c("Dkk2", "Efemp1", "Foxd1", "Osr1", "Ebf2"), nmols = 20000, mols.size = 0.1, mols.alpha = 1)
ggsave("E12.5_pharyngeal.pdf", height = 5, width = 7)
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = T, molecules = c("Dkk2", "Efemp1", "Foxd1", "Osr1", "Ebf2"), nmols = 20000, mols.size = 0.1, mols.alpha = 1) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black"),  
    plot.background = element_rect(fill = "black")
  ) +
  NoLegend()
ggsave("E12.5_pharyngeal_NoLegend.pdf", height = 5, width = 7)


# OFT
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = T, molecules = c("Tbx20","Igfbp7", "Sox9", "Sall3"), nmols = 20000, mols.size = 0.001, mols.alpha = 1) 
ggsave("E12.5_OFT.pdf", height = 5, width = 7)
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = T, molecules = c("Tbx20","Igfbp7", "Sox9", "Sall3"), nmols = 20000, mols.size = 0.001, mols.alpha = 1) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black"),  
    plot.background = element_rect(fill = "black")
  ) +
  NoLegend()
ggsave("E12.5_OFT_NoLegend.pdf", height = 5, width = 7)


# SMC
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = T, molecules = c("Rgs5","Acta2", "Eln"), nmols = 20000, mols.size = 0.001, mols.alpha = 1) 
ggsave("E12.5_SMC.pdf", height = 5, width = 7)
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = T, molecules = c("Rgs5","Acta2", "Eln"), nmols = 20000, mols.size = 0.001, mols.alpha = 1) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black"),  
    plot.background = element_rect(fill = "black")
  ) +
  NoLegend()
ggsave("E12.5_SMC_NoLegend.pdf", height = 5, width = 7)


ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
            coord.fixed = T, molecules = c("Tcf24"), nmols = 20000, mols.size = 0.001, mols.alpha = 1) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black"),  
    plot.background = element_rect(fill = "black")
  ) +
  NoLegend()
ggsave("E12.5_Tcf24_NoLegend.pdf", height = 5, width = 7)



# Focus
cropped.coords <- Crop(xenium.obj[["fov"]], y = c(2500,3500),
                       x = c(3000,3800), coords = "plot")

DefaultAssay(xenium.obj) <- "Xenium"
xenium.obj[["zoom"]] <- cropped.coords
DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentation"

setwd("/Volumes/Pegasus32R8/workspace/omics_data/Xenium/analysis/Wnt1_Cre_E125_8/ImageDimPlot")
# pharynegal
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = T, molecules = c("Dkk2", "Efemp1", "Foxd1", "Osr1", "Ebf2"), nmols = 20000, mols.size = 0.1, mols.alpha = 1)
ggsave("E12.5_pharyngeal_zoom.pdf", height = 5, width = 7)
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = T, molecules = c("Dkk2", "Efemp1", "Foxd1", "Osr1", "Ebf2"), nmols = 20000, mols.size = 0.1, mols.alpha = 1) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black"), 
    plot.background = element_rect(fill = "black")
  ) +
  NoLegend()
ggsave("E12.5_pharyngeal_zoom_NoLegend.pdf", height = 5, width = 7)


# OFT
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = T, molecules = c("Tbx20","Igfbp7", "Sox9", "Sall3"), nmols = 20000, mols.size = 0.001, mols.alpha = 1) 
ggsave("E12.5_OFT_zoom.pdf", height = 5, width = 7)
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = T, molecules = c("Tbx20","Igfbp7", "Sox9", "Sall3"), nmols = 20000, mols.size = 0.001, mols.alpha = 1) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black"), 
    plot.background = element_rect(fill = "black")
  ) +
  NoLegend()
ggsave("E12.5_OFT_zoom_NoLegend.pdf", height = 5, width = 7)


# SMC
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = T, molecules = c("Rgs5","Acta2","Eln"), nmols = 20000, mols.size = 0.001, mols.alpha = 1) 
ggsave("E12.5_zoom_SMC.pdf", height = 5, width = 7)
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = T, molecules = c("Rgs5","Acta2","Eln"), nmols = 20000, mols.size = 0.001, mols.alpha = 1) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black"),  
    plot.background = element_rect(fill = "black")
  ) +
  NoLegend()
ggsave("E12.5_SMC_zoom_NoLegend.pdf", height = 5, width = 7)

# Tcf24
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = T, molecules = c("Tcf24"), nmols = 20000, mols.size = 0.001, mols.alpha = 1) 
ggsave("E12.5_Tcf24_zoom.pdf", height = 5, width = 7)
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = T, molecules = c("Tcf24"), nmols = 20000, mols.size = 0.001, mols.alpha = 1) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black"),  
    plot.background = element_rect(fill = "black")
  ) +
  NoLegend()
ggsave("E12.5_Tcf24_zoom_NoLegend.pdf", height = 5, width = 7)




# cell typeごとの領域を保存する
Idents(xenium.obj) <- "predicted.celltype"
for(i in 0:36){
  ImageDimPlot(xenium.obj, fov = "fov", cells = WhichCells(xenium.obj, idents = i), border.size = 0.01, cols = "yellow") +
    ggtitle(i) +
    scale_y_reverse()
  ggsave(paste("4_decomp_E12.5_", i, ".pdf", sep = ""),  width = 7, height = 5)
}



