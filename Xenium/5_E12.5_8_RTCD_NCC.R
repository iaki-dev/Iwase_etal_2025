# Xenium
# RCTD_allv1ではNCCの領域性を絞ってマッピングしてみる
# さらに縦隔領域をを抜き出すこともしてみる

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
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

# v10.2.1
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/v10/v10.2.1/data/")
gexatac_merge <- readRDS("gexatac_merge_RPCA.rds")
DimPlot(gexatac_merge, reduction = "umap.rpca", group.by = "sub.cluster", label = T)
DefaultAssay(gexatac_merge) <- "RNA"

# setwd("/Users/iwaseakiyasu/workspace/data/Xenium_data/analysis/merge")
setwd("/Volumes/Pegasus32R8/workspace/omics_data/Xenium/analysis/merge")
xenium_merge <- readRDS("xenium_merge.full.rds")

# setwd("/Users/iwaseakiyasu/workspace/data/Xenium_data/analysis/Wnt1_Cre_E125_8")
setwd("/Volumes/Pegasus32R8/workspace/omics_data/Xenium/analysis/Wnt1_Cre_E125_8")
xenium.obj <- readRDS("Wnt1_E12.5_8.rds")
setwd("/Volumes/Pegasus32R8/workspace/omics_data/Xenium/analysis/Wnt1_Cre_E125_8/RCTD")
# setwd("/Users/iwaseakiyasu/workspace/data/Xenium_data/analysis/Wnt1_Cre_E125_8/RCTD")
coi <- xenium_merge@meta.data %>% dplyr::filter(section == "E12.5_8")
xenium_merge.sub <- xenium_merge[, rownames(coi)]
xenium_merge.sub@meta.data %>% rownames() %>% head
# [1] "aaaaeako-1_E12.5_8" "aaabhimh-1_E12.5_8" "aaacbago-1_E12.5_8" "aaadcjee-1_E12.5_8" "aaaelkoh-1_E12.5_8" "aaaemdkm-1_E12.5_8"

xenium.obj@meta.data %>% rownames() %>% head
# [1] "aaaaeako-1" "aaabhimh-1" "aaacbago-1" "aaadcjee-1" "aaaelkoh-1" "aaaemdkm-1"

# xenium_merge.subから名前の不要部分を削除する
tmp.name <- xenium_merge.sub@meta.data %>% rownames() %>% str_split(pattern = "_", simplify = T)

tmp.name[,1] %>% head
# 順番が揃っている

xenium.obj@meta.data$celltype.full <- xenium_merge.sub@meta.data$celltype.full
xenium.obj@meta.data$celltype.full <- factor(xenium.obj@meta.data$celltype.full, 
                                             levels = c(0:38))
xenium.obj@meta.data$celltype.full
DimPlot(xenium.obj, group.by = "celltype.full", label = T)
ggsave("3_umap_celltype.full.pdf", width = 5, height = 5)


# さらにNCCが分布するサンプリングした縦隔領域を抽出する
# つまり、celltype.fullの0, 15, 17 , 25, 27, 32, 38
coi <- xenium.obj@meta.data %>% dplyr::filter(celltype.full %in% c("0", "15", "17" , "25", "27", "32", "38"))
xenium.obj <- xenium.obj[, rownames(coi)]

# さらに、物理的にサンプリングした縦隔領域もフィルタリングに加える
com <- read_csv("Mediastinum_cells_stats.csv", skip = 2)
com <- com$`Cell ID`
xenium.obj <- xenium.obj[, com]

ImageDimPlot(xenium.obj, fov = "fov", axes = F, border.color = "white", border.size = 0.01, cols = "polychrome",
             coord.fixed = FALSE, group.by = "celltype.full") +
  scale_y_reverse()
setwd("/Users/iwaseakiyasu/workspace/data/Xenium_data/analysis/Wnt1_Cre_E125_8/RCTD/NCC")
ggsave("3_ImageDim_celltype.full_NCC_Mediastinum.pdf", width = 7, height = 5)


meta <- xenium.obj@meta.data
write.table(meta, "meta.txt", quote = F, sep = "\t", col.names = NA)



# RTCD cell type decomposition
query.counts <- GetAssayData(xenium.obj, assay = "Xenium", layer = "counts")[, Cells(xenium.obj[["fov"]])]
coords <- GetTissueCoordinates(xenium.obj[["fov"]], which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))


# NCCの属するクラスタを定義する
meta.tmp <- gexatac_merge@meta.data$sample %>%
  str_replace_all(pattern = "E11.5_NCC", replacement = "NCC") %>%
  str_replace_all(pattern = "E12.5_NCC", replacement = "NCC") %>%
  str_replace_all(pattern = "E14.5_NCC", replacement = "NCC") %>%
  str_replace_all(pattern = "E17.5_NCC", replacement = "NCC") %>%
  str_replace_all(pattern = "E11.5_nonNCC", replacement = "nonNCC") %>%
  str_replace_all(pattern = "E12.5_nonNCC", replacement = "nonNCC") %>%
  str_replace_all(pattern = "E14.5_nonNCC", replacement = "nonNCC") %>%
  str_replace_all(pattern = "E17.5_nonNCC", replacement = "nonNCC")
gexatac_merge <- AddMetaData(gexatac_merge, metadata = meta.tmp, col.name = "lineage")
DimPlot(gexatac_merge, group.by = "lineage", reduction = "umap.rpca")
# ggsave("umap_lineage.pdf", height = 6, width = 6)

# sampleごとに分けて表示できるようにする
table_lineage <- gexatac_merge@meta.data %>% dplyr::select(sub.cluster, sample) %>% table %>% as.data.frame()
table_lineage <- table_lineage %>% pivot_wider(names_from = sample, values_from = Freq)
table_lineage <- table_lineage %>% mutate(
  total_cells = table_lineage %>% dplyr::select(-sub.cluster) %>% apply(1, sum)
)
table_lineage %>% head()

# NCCに着目する
using.cluster_lineages <- table_lineage %>% dplyr::select(sub.cluster, E11.5_NCC, E12.5_NCC, total_cells)
sum_col <- using.cluster_lineages %>% dplyr::select(E11.5_NCC, E12.5_NCC) %>% apply(1, sum)
using.cluster_lineages$sum_col <- sum_col
using.cluster_lineages <- using.cluster_lineages %>% mutate(ratio_NCC =  (sum_col / total_cells * 100))
# 5%以上で切る
using.cluster_lineages <- using.cluster_lineages %>% dplyr::filter(ratio_NCC > 5)
using.cluster_lineages <- using.cluster_lineages$sub.cluster %>% as.character()
using.cluster_lineages

# E11.5とE12.5 (scMultiome) のうち5%以上のNCCを含むデータを解析対象とする
coi <- gexatac_merge@meta.data %>% dplyr::filter(sample %in% c("E11.5_NCC", "E12.5_NCC"),
                                                 seurat_clusters %in% using.cluster_lineages)
NCC <- gexatac_merge[, rownames(coi)]
NCC$sub.cluster <- NCC$sub.cluster %>% as.character()
order_factor <- NCC$sub.cluster %>% as.numeric() %>% unique() %>% sort()
NCC$sub.cluster <- as.factor(NCC$sub.cluster)
NCC$sub.cluster <- factor(NCC$sub.cluster,
                            levels = order_factor)
Idents(NCC) <- NCC@meta.data$sub.cluster
DimPlot(NCC)
counts <- GetAssayData(NCC[, rownames(coi)], layer = "RNA", slot = "counts")
# counts@xが整数でないとreferenceがworkしない
# NCC_NCC, NCC_NCCがworkしなかったので、これはsoupXのせい
# roundで整数値する
counts@x <- round(counts@x)

# cluster <- as.factor(NCC$cluster_lineage)
cluster <- as.factor(NCC[, rownames(coi)]$sub.cluster)
names(cluster) <- colnames(NCC[, rownames(coi)])
nUMI <- NCC[, rownames(coi)]$nCount_RNA
names(nUMI) <- colnames(NCC[, rownames(coi)])
nUMI <- colSums(counts)
reference <- Reference(counts, cluster, nUMI)

# run RCTD with many cores
RCTD <- create.RCTD(query, reference, max_cores = 8,
                    CELL_MIN_INSTANCE = 1)

# doublet mode, which assigns 1-2 cell types per spot and is recommended for technologies with
# high spatial resolution such as Slide-seq and MERFISH
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
# saveRDS(RCTD, "RCTD_NCC_celltype_doublet_NCC_Mediastinum.rds")

RCTD <- readRDS("RCTD_NCC_celltype_doublet_NCC_Mediastinum.rds")

annotations.df <- RCTD@results$results_df
annotations.df %>% ggplot(aes(x=singlet_score)) +
  geom_histogram()
ggsave("hist_siglet_score.pdf", height = 4, width = 5)
annotations.df$singlet_score %>% summary()

df_weights <- RCTD@results$weights %>% as.data.frame()
df_weights <- cbind(rownames(df_weights) ,df_weights)
colnames(df_weights)[1] <- "cell_id"
df_weights %>% dplyr::filter(cell_id == "khhhkaob-1") %>% t() %>% as.data.frame() %>% View

# first_typeだとrejectされたものも含まれるので注意して変更する
# annotations <- annotations.df$first_type
# names(annotations) <- rownames(annotations.df)
# xenium.obj$predicted.celltype <- annotations
# keep.cells <- Cells(xenium.obj)[!is.na(xenium.obj$predicted.celltype)]
# xenium.obj <- subset(xenium.obj, cells = keep.cells)
# rejectとコールされたものだけ除く
# annotations.df <- annotations.df %>% dplyr::filter(spot_class == "singlet")
# xenium.obj <- xenium.obj[, rownames(annotations)]

annotations <- annotations.df$first_type
names(annotations) <- rownames(annotations.df)
xenium.obj$predicted.celltype <- annotations
keep.cells <- Cells(xenium.obj)[!is.na(xenium.obj$predicted.celltype)]
xenium.obj <- subset(xenium.obj, cells = keep.cells)

xenium.obj@meta.data$predicted.celltype %>% class
arrange_factor <- xenium.obj@meta.data$predicted.celltype %>% as.character() %>% unique() %>% as.numeric() %>% sort()
xenium.obj@meta.data$predicted.celltype <- factor(xenium.obj@meta.data$predicted.celltype,
                                                  levels = arrange_factor)
ImageDimPlot(xenium.obj, fov =  "fov", group.by = "predicted.celltype", size = 1.5, border.size = 0.01, cols = "polychrome", dark.background = T) +
  ggtitle("Cell type")  +
  scale_y_reverse()
ggsave("4_decomp_NCC_NCC.pdf", width = 7, height = 7)
write.table(xenium.obj@meta.data, "4_meta.predicted.celltype.txt", quote = F, sep = "\t", col.names = NA)

# save
saveRDS(xenium.obj, "xenium.obj_NCC.rds")


# # multi mode, an extension of doublet mode that can discover more than two cell types per spot
# # as an alternative option to full mode
# RCTD_NCC_multi <- run.RCTD(RCTD, doublet_mode = "multi") workしなかった
# saveRDS(RCTD, "RCTD_E11.5_celltype_multi.rds")
# df_weights <- RCTD@results$weights %>% as.data.frame()
# df_weights <- cbind(rownames(df_weights) ,df_weights)
# colnames(df_weights)[1] <- "cell_id"
# df_weights %>% dplyr::filter(cell_id == "kchcmkka-1") %>% t() %>% as.data.frame() %>% View


ここまで






# cell typeごとの領域を保存する
Idents(xenium.obj) <- "predicted.celltype"
for(i in 0:36){
  ImageDimPlot(xenium.obj, fov = "fov", cells = WhichCells(xenium.obj, idents = i), border.size = 0.01, cols = "yellow") +
    ggtitle(i) +
    scale_y_reverse()
  ggsave(paste("4_decomp_E11.5_", i, ".pdf", sep = ""),  width = 7, height = 5)
}







# Niche
xenium.obj <- BuildNicheAssay(object = xenium.obj, fov = "fov", group.by = "predicted.celltype",
                              niches.k = 5, neighbors.k = 10)

niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches") +
  scale_fill_manual(values = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B"))
celltype.plot | niche.plot

table(xenium.obj$predicted.celltype, xenium.obj$niches)


meta <- xenium.obj@meta.data
write.table()



# # 拡大する
# cropped.coords <- Crop(xenium.obj[["fov"]], x = c(1200, 2900), y = c(3750, 4550), coords = "plot")
# xenium.obj[["zoom"]] <- cropped.coords
# # visualize cropped area with cell segmentations & selected molecules
# DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentation"
# ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
#              coord.fixed = FALSE, molecules = c("Gad1", "Sst", "Npy2r", "Pvalb", "Nrn1"), nmols = 10000)


