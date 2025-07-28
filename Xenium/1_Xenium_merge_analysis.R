# Xenium
# sketch integration

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
options(future.globals.maxSize = 3e9)
options(Seurat.object.assay.version = "v5")

setwd("/Users/iwaseakiyasu/workspace/data/Xenium_data/analysis/merge/")
path <- "/Users/iwaseakiyasu/workspace/data/Xenium_data/original_data/output-XETG00053__0016856__Wnt1_Cre_E115_1__20240124__080118"
# Load the Xenium data
E11.5_1 <- LoadXenium(path, fov = "fov")
# remove cells with 0 counts
E11.5_1 <- subset(E11.5_1, subset = nCount_Xenium > 0)
dim(E11.5_1)
# [1]   334 77263
VlnPlot(E11.5_1, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)

path <- "/Users/iwaseakiyasu/workspace/data/Xenium_data/original_data/output-XETG00053__0016856__Wnt1_Cre_E115_3__20240124__080118"
# Load the Xenium data
E11.5_3<- LoadXenium(path, fov = "fov")
# remove cells with 0 counts
E11.5_3 <- subset(E11.5_3, subset = nCount_Xenium > 0)
dim(E11.5_3)
# [1]   334 77263
VlnPlot(E11.5_3, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)

path <- "/Users/iwaseakiyasu/workspace/data/Xenium_data/original_data/output-XETG00053__0016856__Wnt1_Cre_E125_5__20240124__080118"
# Load the Xenium data
E12.5_5 <- LoadXenium(path, fov = "fov")
# remove cells with 0 counts
E12.5_5 <- subset(E12.5_5, subset = nCount_Xenium > 0)
dim(E12.5_5)
# [1]   334 103509
VlnPlot(E12.5_5, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)

path <- "/Users/iwaseakiyasu/workspace/data/Xenium_data/original_data/output-XETG00053__0016856__Wnt1_Cre_E125_7__20240124__080118"
# Load the Xenium data
E12.5_7 <- LoadXenium(path, fov = "fov")
# remove cells with 0 counts
E12.5_7 <- subset(E12.5_7, subset = nCount_Xenium > 0)
dim(E12.5_7)
# [1]   334 103509
VlnPlot(E12.5_7, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)


path <- "/Users/iwaseakiyasu/workspace/data/Xenium_data/original_data/output-XETG00053__0016856__Wnt1_Cre_E125_8__20240124__080118"
# Load the Xenium data
E12.5_8 <- LoadXenium(path, fov = "fov")
# remove cells with 0 counts
E12.5_8 <- subset(E12.5_8, subset = nCount_Xenium > 0)
dim(E12.5_8)
# [1]   334 103509
VlnPlot(E12.5_8, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)


path <- "/Users/iwaseakiyasu/workspace/data/Xenium_data/original_data/output-XETG00053__0017041__Wnt1_Cre_E125_9__20240124__080118"
# Load the Xenium data
E12.5_9 <- LoadXenium(path, fov = "fov")
# remove cells with 0 counts
E12.5_9 <- subset(E12.5_9, subset = nCount_Xenium > 0)
dim(E12.5_9)
# [1]   334 103509
VlnPlot(E12.5_9, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)



path <- "/Users/iwaseakiyasu/workspace/data/Xenium_data/original_data/output-XETG00053__0017041__Csf1r_Cre_E115_10__20240124__080118"
# Load the Xenium data
E11.5_10 <- LoadXenium(path, fov = "fov")
# remove cells with 0 counts
E11.5_10 <- subset(E11.5_10, subset = nCount_Xenium > 0)
dim(E11.5_10)
# [1]   334 103509
VlnPlot(E11.5_10, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)


path <- "/Users/iwaseakiyasu/workspace/data/Xenium_data/original_data/output-XETG00053__0017041__Csf1r_Cre_E125_13__20240124__080119"
# Load the Xenium data
E12.5_13 <- LoadXenium(path, fov = "fov")
# remove cells with 0 counts
E12.5_13 <- subset(E12.5_13, subset = nCount_Xenium > 0)
dim(E12.5_13)
# [1]   334 103509
VlnPlot(E12.5_13, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)


path <- "/Users/iwaseakiyasu/workspace/data/Xenium_data/original_data/output-XETG00053__0017041__Csf1r_Cre_E125_14__20240124__080119"
# Load the Xenium data
E12.5_14 <- LoadXenium(path, fov = "fov")
# remove cells with 0 counts
E12.5_14 <- subset(E12.5_14, subset = nCount_Xenium > 0)
dim(E12.5_14)
# [1]   334 103509
VlnPlot(E12.5_14, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)


path <- "/Users/iwaseakiyasu/workspace/data/Xenium_data/original_data/output-XETG00053__0017041__Csf1r_Cre_E125_15__20240124__080119"
# Load the Xenium data
E12.5_15 <- LoadXenium(path, fov = "fov")
# remove cells with 0 counts
E12.5_15 <- subset(E12.5_15, subset = nCount_Xenium > 0)
dim(E12.5_15)
# [1]   334 103509
VlnPlot(E12.5_15, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)


path <- "/Users/iwaseakiyasu/workspace/data/Xenium_data/original_data/output-XETG00053__0017041__Csf1r_Cre_E125_16__20240124__080119"
# Load the Xenium data
E12.5_16 <- LoadXenium(path, fov = "fov")
# remove cells with 0 counts
E12.5_16 <- subset(E12.5_16, subset = nCount_Xenium > 0)
dim(E12.5_16)
# [1]   334 103509
VlnPlot(E12.5_16, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)


library(data.table)
df1 <- E11.5_1@assays$Xenium$counts %>% as.data.table()
colnames(df1) <- paste(colnames(df1), "E11.5_1", sep = "_")
df1$feature <- rownames(E11.5_1@assays$Xenium$counts)

df2 <- E11.5_3@assays$Xenium$counts %>% as.data.table()
colnames(df2) <- paste(colnames(df2), "E11.5_3", sep = "_")
df2$feature <- rownames(E11.5_3@assays$Xenium$counts)

df3 <- E11.5_10@assays$Xenium$counts %>% as.data.table()
colnames(df3) <- paste(colnames(df3), "E11.5_10", sep = "_")
df3$feature <- rownames(E11.5_10@assays$Xenium$counts)

df4 <- E12.5_5@assays$Xenium$counts %>% as.data.table()
colnames(df4) <- paste(colnames(df4), "E12.5_5", sep = "_")
df4$feature <- rownames(E12.5_5@assays$Xenium$counts)

df5 <- E12.5_7@assays$Xenium$counts %>% as.data.table()
colnames(df5) <- paste(colnames(df5), "E12.5_7", sep = "_")
df5$feature <- rownames(E12.5_7@assays$Xenium$counts)

df6 <- E12.5_8@assays$Xenium$counts %>% as.data.table()
colnames(df6) <- paste(colnames(df6), "E12.5_8", sep = "_")
df6$feature <- rownames(E12.5_8@assays$Xenium$counts)

df7 <- E12.5_9@assays$Xenium$counts %>% as.data.table()
colnames(df7) <- paste(colnames(df7), "E12.5_9", sep = "_")
df7$feature <- rownames(E12.5_9@assays$Xenium$counts)

df8 <- E12.5_13@assays$Xenium$counts %>% as.data.table()
colnames(df8) <- paste(colnames(df8), "E12.5_13", sep = "_")
df8$feature <- rownames(E12.5_13@assays$Xenium$counts)

df9 <- E12.5_14@assays$Xenium$counts %>% as.data.table()
colnames(df9) <- paste(colnames(df9), "E12.5_14", sep = "_")
df9$feature <- rownames(E12.5_14@assays$Xenium$counts)

df10 <- E12.5_15@assays$Xenium$counts %>% as.data.table()
colnames(df10) <- paste(colnames(df10), "E12.5_15", sep = "_")
df10$feature <- rownames(E12.5_15@assays$Xenium$counts)

df11 <- E12.5_16@assays$Xenium$counts %>% as.data.table()
colnames(df11) <- paste(colnames(df11), "E12.5_16", sep = "_")
df11$feature <- rownames(E12.5_16@assays$Xenium$counts)

df_merge <- merge.data.table(df1, df2, by = "feature", all = T) %>% 
  merge.data.table(df3, by = "feature", all = T) %>% 
  merge.data.table(df4, by = "feature", all = T) %>% 
  merge.data.table(df5, by = "feature", all = T) %>% 
  merge.data.table(df6, by = "feature", all = T) %>% 
  merge.data.table(df7, by = "feature", all = T) %>% 
  merge.data.table(df8, by = "feature", all = T) %>% 
  merge.data.table(df9, by = "feature", all = T) %>% 
  merge.data.table(df10, by = "feature", all = T) %>% 
  merge.data.table(df11, by = "feature", all = T) 
dim(df_merge)
# [1]     334 1003928
# save
# saveRDS(df_merge, "df_merge.matrix.rds")
df_merge <- readRDS("df_merge.matrix.rds")

meta <- colnames(df_merge) %>% as.data.table()
meta <- meta[-1,]
colnames(meta) <- "sample"
section <- meta$sample %>% str_split(pattern = "_", simplify = T) %>% as.data.frame()
section$section <- paste(section$V2, section$V3, sep = "_")
colnames(section) <- c("id", "stage", "number", "section")
rownames(section) <- meta$sample
  

df_merge <- as.data.frame(df_merge)
rownames(df_merge) <- df_merge$feature
df_merge <- df_merge[, -1]

xenium_merge <- CreateSeuratObject(counts = df_merge, meta.data = section,
                                   assay = "RNA")

# saveRDS(xenium_merge, "xenium_merge.rds")

xenium_merge <- readRDS("xenium_merge.rds")

xenium_merge <- NormalizeData(xenium_merge)

# split assay 
xenium_merge[["RNA"]] <- split(xenium_merge[["RNA"]], f = xenium_merge$section)
xenium_merge <- FindVariableFeatures(xenium_merge, verbose = FALSE)

# Sample representative cells from each datase
xenium_merge <- SketchData(object = xenium_merge, ncells = 5000, method = "LeverageScore", sketched.assay = "sketch")
xenium_merge

# Perform integration on the sketched cells across sample
DefaultAssay(xenium_merge) <- "sketch"
xenium_merge <- FindVariableFeatures(xenium_merge, verbose = F)
xenium_merge <- ScaleData(xenium_merge, verbose = F)
xenium_merge <- RunPCA(xenium_merge, verbose = F)
# integrate the datasets
xenium_merge <- IntegrateLayers(xenium_merge, method = RPCAIntegration, orig = "pca", new.reduction = "integrated.rpca",
                          dims = 1:30, k.anchor = 20, 
                          reference = which(Layers(xenium_merge, search = "data") %in% c("data.E12.5_8")),
                          verbose = F)
# cluster the integrated data
xenium_merge <- FindNeighbors(xenium_merge, reduction = "integrated.rpca", dims = 1:30)
xenium_merge <- FindClusters(xenium_merge, resolution = 1.5)

xenium_merge <- RunUMAP(xenium_merge, reduction = "integrated.rpca", dims = 1:30, return.model = T, verbose = F)
DimPlot(xenium_merge, group.by = "seurat_clusters", reduction = "umap", label = T)
ggsave("umap.xenium_merge.pdf", width = 5, height = 5)
xenium_merge@meta.data$section <- as.factor(xenium_merge@meta.data$section)
xenium_merge@meta.data$section <- factor(xenium_merge@meta.data$section,
                                         levels = c("E11.5_1", "E11.5_3", "E11.5_10", 
                                                    "E12.5_5", "E12.5_7", "E12.5_8", "E12.5_9",
                                                    "E12.5_13", "E12.5_14", "E12.5_15", "E12.5_16"
                                                    ))
DimPlot(xenium_merge, group.by = "section", reduction = "umap", label = F)
ggsave("umap.xenium_merge_section.pdf", width = 5, height = 5)


# you can now rejoin the layers in the sketched assay this is required to perform differential
# expression
xenium_merge[["sketch"]] <- JoinLayers(xenium_merge[["sketch"]])


all.markers <- FindAllMarkers(object = xenium_merge, max.cells.per.ident = 500, only.pos = TRUE)
write.table(all.markers, "all.markers.txt", sep = "\t", quote = F, row.names = F)


FeaturePlot(xenium_merge, features = c("Tcf24", "Acta2", "Myh11", "Ebf2", "Osr1", "Tbx20",
                                       "Cyp26b1", "Dlx5", "Runx2", 
                                       "Csf1r", "Ccr2", "Mrc1", "Itgb2","Wt1", "Aldh1a2", 
                                       "Cdh5", "Prox1", "Lyve1", "Aplnr", "Apold1","Tek", "Flt4", "Gata3",
                                       "Col2a1", "Sox9", "Twist1",
                                       "Myog", "Tnnt2", "Nppa",
                                       "Nrxn1", "Sox10", "Gpnmb","Elavl3", "Nrxn1", "Sall3","Sox2","Cdh1", "Tbx3",
                                       "Foxo3",
                                       "Igfbp7", "Abcc9", "Kdr","Sall3","Pdgfa"))
ggsave("Featureplot_markers.pdf", height = 20, width = 12)


# saveRDS(xenium_merge, "xenium_merge_sketch_subset.rds")

# library(loupeR)
# setup()
# xenium_merge.tmp <- xenium_merge
# xenium_merge.tmp[["sketch"]] <- JoinLayers(xenium_merge.tmp[["sketch"]])
# create_loupe_from_seurat(xenium_merge.tmp)




# Integrate the full datasets
# resplit the sketched cell assay into layers this is required to project the integration onto
# all cells
xenium_merge[["sketch"]] <- split(xenium_merge[["sketch"]], f = xenium_merge$section)

xenium_merge <- ProjectIntegration(object = xenium_merge, sketched.assay = "sketch", assay = "RNA", reduction = "integrated.rpca")


xenium_merge <- ProjectData(object = xenium_merge, sketched.assay = "sketch", assay = "RNA", 
                            sketched.reduction = "integrated.rpca.full",
                      full.reduction = "integrated.rpca.full", dims = 1:30, refdata = list(celltype.full = "seurat_clusters"))

xenium_merge <- RunUMAP(xenium_merge, reduction = "integrated.rpca.full", dims = 1:30, reduction.name = "umap.full",
                  reduction.key = "UMAP_full_")

DimPlot(xenium_merge, reduction = "umap.full", group.by = "section", alpha = 0.1)
ggsave("umap.full.xenium_merge_section.pdf", width = 5, height = 5)


order_factor <- xenium_merge@meta.data$celltype.full %>% unique() %>% as.numeric() %>% sort()
xenium_merge@meta.data$celltype.full <- as.factor(xenium_merge@meta.data$celltype.full)
xenium_merge@meta.data$celltype.full <- factor(xenium_merge@meta.data$celltype.full,
                                               levels = order_factor)
DimPlot(xenium_merge, reduction = "umap.full", group.by = "celltype.full", alpha = 0.1, label = T)
ggsave("umap.full.xenium_merge_celltype.pdf", width = 5, height = 5)

# saveRDS(xenium_merge, "xenium_merge.full.rds")

# xenium_mergeのアノテーション
new.cluster.ids <- c("0_EBf2+Mes",
                     "1_Aplnr+EC",
                     "2_Twist1+Mes",
                     "3_Myocyte",
                     "4_Macrophage",
                     "5_Mes",
                     "6_Osteoblast",
                     "7_ambiguous",
                     "8_Foxd1+Mes" ,
                     "9_Spinal_mantle",
                     "10_Chondroblast1",
                     "11_Epithelium",
                     "12_Thoracic_wall",
                     "13_Ccl2+Mes",
                     "14_Myh7+Myocyte",
                     "15_SMC",
                     "16_Mesothelium",
                     "17_Schwann",
                     "18_Blood",
                     "19_Ventricle_CM",
                     "20_Marginal_layer",
                     "21_Chondroblast2",
                     "22_Monocyte",
                     "23_Trachea_Esophagus_Mes",
                     "24_Spinal_ependymal",
                     "25_Cushion_Mes",
                     "26_Thoracic_Mes",
                     "27_Efemp1+Mes",
                     "28_Pericytes",
                     "29_Endocarium",
                     "30_Arterial_EC",
                     "31_Atrium_Venticle_CM",
                     "32_Ganglion",
                     "33_LEC",
                     "34_ambious",
                     "35_Sox10+Sox9+Chondroblast",
                     "36_Trachea_Esophagus_epithelium",
                     "37_Melanocyte",
                     "38_Neuron"
)
names(new.cluster.ids) <- levels(xenium_merge)
xenium_merge <- RenameIdents(xenium_merge, new.cluster.ids)
xenium_merge@meta.data$celltype_label <- Idents(xenium_merge)
DimPlot(xenium_merge, reduction = "umap.full", group.by = "celltype_label", alpha = 0.1, label = F)
ggsave("umap.full.xenium_merge_celltype_laabel.pdf", width = 8, height = 5)








coi <- gexatac_merge@meta.data %>% dplyr::filter(sample == "E12.5_nonNCC")
tmp <- gexatac_merge[, rownames(coi)]
FeaturePlot(tmp, features = c("Tbx20", "Sox10", "Tubb3"), reduction = "umap.rpca", order = T)

# Pseudobulk analysis
setwd("/Users/iwaseakiyasu/workspace/data/Xenium_data/analysis/merge")
setwd("/Volumes/Pegasus32 R8/workspace/omics_data/Xenium/analysis/merge")
xenium_merge <- readRDS("xenium_merge.full.rds")
all.markers <- read.table("all.markers.txt", header = T)
all.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10


bulk <- AggregateExpression(xenium_merge, return.seurat = T, slot = "data", assays = "sketch",
                            group.by = c("celltype.full"))
bulk@meta.data$celltype.full <- bulk@meta.data$celltype.full %>% str_replace_all(pattern = "g", replacement = "")
order_factor <- bulk@meta.data$celltype.full %>% unique() %>% as.numeric() %>% sort()
bulk@meta.data$celltype.full <- as.factor(bulk@meta.data$celltype.full)
bulk@meta.data$celltype.full <- factor(bulk@meta.data$celltype.full,
                                       levels = order_factor)
Idents(bulk) <- bulk@meta.data$celltype.full
tail(Cells(bulk))

DoHeatmap(bulk, features = top10$gene, draw.lines = FALSE)  + NoLegend()
ggsave("heatmap_markeres_log2FC_top10_pseudobulk.png", height = 25, width = 15)














markers_rna %>% dplyr::filter(avg_log2FC > 0.25, pct.1 > 0.3) %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
top10$cluster <- as.numeric(top10$cluster)
top10 <- top10 %>% arrange(cluster)
DefaultAssay(gexatac_merge) <- "RNA"
DoHeatmap(gexatac_merge, features = top10$gene) + NoLegend()
ggsave("heatmap.top10_log2FC.png", height = 25, width = 20)
DoHeatmap(cluster.averages, features = top10$gene, draw.lines = FALSE)  + NoLegend()
ggsave("heatmap_markeres_log2FC_top10_pseudobulk.png", height = 30, width = 8)
DoHeatmap(cluster.averages, features = top10$gene, draw.lines = FALSE)  + NoLegend()  + RotatedAxis() + coord_flip()
ggsave("heatmap_markeres_log2FC_top10_pseudobulk_rotated.png", height = 4, width = 35)
