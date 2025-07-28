

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
options(future.globals.maxSize = 5 * 1024^3)


s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

s.genes <- paste(s.genes %>% str_sub(start = 1, end = 1), s.genes %>% str_sub(start = 2) %>% str_to_lower(), sep = "")
g2m.genes <- paste(g2m.genes %>% str_sub(start = 1, end = 1), g2m.genes %>% str_sub(start = 2) %>% str_to_lower(), sep = "")


setwd("/Volumes/Pegasus32 R8/workspace/public_data/NCC_NatureComm2023/analysis")
E8.5 <- readRDS("E8.5_doublet_velo.rds")
E9.5 <- readRDS("E9.5_doublet_velo.rds")
E10.5_1 <- readRDS("E10.5_1_doublet_velo.rds")
E10.5_2 <- readRDS("E10.5_2_doublet_velo.rds")
setwd("/Volumes/Pegasus32 R8/workspace/public_data/NCC_EMBO2021/analysis/")
E10.5_3 <- readRDS("E10.5_3_doublet_velo.rds")
E11.5 <- readRDS("E11.5_doublet_velo.rds")
E12.5 <- readRDS("E12.5_doublet_velo.rds")
E13.5 <- readRDS("E13.5_doublet_velo.rds")
E14.5 <- readRDS("E14.5_doublet_velo.rds")
E17.5 <- readRDS("E17.5_doublet_velo.rds")
P1 <- readRDS("P1_doublet_velo.rds")
P7 <- readRDS("P7_doublet_velo.rds")


coi <- E8.5@meta.data %>% dplyr::filter(DF == "Singlet")
E8.5 <- E8.5[, rownames(coi)]
E8.5 <- NormalizeData(E8.5)
E8.5 <- FindVariableFeatures(E8.5, selection.method = "vst")
E8.5 <- ScaleData(E8.5, features = rownames(E8.5))
E8.5 <- RunPCA(E8.5, features = VariableFeatures(E8.5, layer = "counts"))
E8.5 <- CellCycleScoring(E8.5, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

coi <- E9.5@meta.data %>% dplyr::filter(DF == "Singlet")
E9.5 <- E9.5[, rownames(coi)]
E9.5 <- NormalizeData(E9.5)
E9.5 <- FindVariableFeatures(E9.5, selection.method = "vst")
E9.5 <- ScaleData(E9.5, features = rownames(E9.5))
E9.5 <- RunPCA(E9.5, features = VariableFeatures(E9.5, layer = "counts"))
E9.5 <- CellCycleScoring(E9.5, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

coi <- E10.5_1@meta.data %>% dplyr::filter(DF == "Singlet")
E10.5_1 <- E10.5_1[, rownames(coi)]
E10.5_1 <- NormalizeData(E10.5_1)
E10.5_1 <- FindVariableFeatures(E10.5_1, selection.method = "vst")
E10.5_1 <- ScaleData(E10.5_1, features = rownames(E10.5_1))
E10.5_1 <- RunPCA(E10.5_1, features = VariableFeatures(E10.5_1, layer = "counts"))
E10.5_1 <- CellCycleScoring(E10.5_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


coi <- E10.5_2@meta.data %>% dplyr::filter(DF == "Singlet")
E10.5_2 <- E10.5_2[, rownames(coi)]
E10.5_2 <- NormalizeData(E10.5_2)
E10.5_2 <- FindVariableFeatures(E10.5_2, selection.method = "vst")
E10.5_2 <- ScaleData(E10.5_2, features = rownames(E10.5_2))
E10.5_2 <- RunPCA(E10.5_2, features = VariableFeatures(E10.5_2, layer = "counts"))
E10.5_2 <- CellCycleScoring(E10.5_2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

coi <- E10.5_3@meta.data %>% dplyr::filter(DF == "Singlet")
E10.5_3 <- E10.5_3[, rownames(coi)]
E10.5_3 <- NormalizeData(E10.5_3)
E10.5_3 <- FindVariableFeatures(E10.5_3, selection.method = "vst")
E10.5_3 <- ScaleData(E10.5_3, features = rownames(E10.5_3))
E10.5_3 <- RunPCA(E10.5_3, features = VariableFeatures(E10.5_3, layer = "counts"))
E10.5_3 <- CellCycleScoring(E10.5_3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

coi <- E11.5@meta.data %>% dplyr::filter(DF == "Singlet")
E11.5 <- E11.5[, rownames(coi)]
E11.5 <- NormalizeData(E11.5)
E11.5 <- FindVariableFeatures(E11.5, selection.method = "vst")
E11.5 <- ScaleData(E11.5, features = rownames(E11.5))
E11.5 <- RunPCA(E11.5, features = VariableFeatures(E11.5, layer = "counts"))
E11.5 <- CellCycleScoring(E11.5, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

coi <- E12.5@meta.data %>% dplyr::filter(DF == "Singlet")
E12.5 <- E12.5[, rownames(coi)]
E12.5 <- NormalizeData(E12.5)
E12.5 <- FindVariableFeatures(E12.5, selection.method = "vst")
E12.5 <- ScaleData(E12.5, features = rownames(E12.5))
E12.5 <- RunPCA(E12.5, features = VariableFeatures(E12.5, layer = "counts"))
E12.5 <- CellCycleScoring(E12.5, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

coi <- E13.5@meta.data %>% dplyr::filter(DF == "Singlet")
E13.5 <- E13.5[, rownames(coi)]
E13.5 <- NormalizeData(E13.5)
E13.5 <- FindVariableFeatures(E13.5, selection.method = "vst")
E13.5 <- ScaleData(E13.5, features = rownames(E13.5))
E13.5 <- RunPCA(E13.5, features = VariableFeatures(E13.5, layer = "counts"))
E13.5 <- CellCycleScoring(E13.5, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

coi <- E14.5@meta.data %>% dplyr::filter(DF == "Singlet")
E14.5 <- E14.5[, rownames(coi)]
E14.5 <- NormalizeData(E14.5)
E14.5 <- FindVariableFeatures(E14.5, selection.method = "vst")
E14.5 <- ScaleData(E14.5, features = rownames(E14.5))
E14.5 <- RunPCA(E14.5, features = VariableFeatures(E14.5, layer = "counts"))
E14.5 <- CellCycleScoring(E14.5, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

coi <- E17.5@meta.data %>% dplyr::filter(DF == "Singlet")
E17.5 <- E17.5[, rownames(coi)]
E17.5 <- NormalizeData(E17.5)
E17.5 <- FindVariableFeatures(E17.5, selection.method = "vst")
E17.5 <- ScaleData(E17.5, features = rownames(E17.5))
E17.5 <- RunPCA(E17.5, features = VariableFeatures(E17.5, layer = "counts"))
E17.5 <- CellCycleScoring(E17.5, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

coi <- P1@meta.data %>% dplyr::filter(DF == "Singlet")
P1 <- P1[, rownames(coi)]
P1 <- NormalizeData(P1)
P1 <- FindVariableFeatures(P1, selection.method = "vst")
P1 <- ScaleData(P1, features = rownames(P1))
P1 <- RunPCA(P1, features = VariableFeatures(P1, layer = "counts"))
P1 <- CellCycleScoring(P1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

coi <- P7@meta.data %>% dplyr::filter(DF == "Singlet")
P7 <- P7[, rownames(coi)]
P7 <- NormalizeData(P7)
P7 <- FindVariableFeatures(P7, selection.method = "vst")
P7 <- ScaleData(P7, features = rownames(P7))
P7 <- RunPCA(P7, features = VariableFeatures(P7, layer = "counts"))
P7 <- CellCycleScoring(P7, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


coi <- E8.5@meta.data %>% dplyr::filter(DF == "Singlet")
E8.5 <- E8.5[, rownames(coi)]

coi <- E9.5@meta.data %>% dplyr::filter(DF == "Singlet")
E9.5 <- E9.5[, rownames(coi)]

coi <- E10.5_1@meta.data %>% dplyr::filter(DF == "Singlet")
E10.5_1 <- E10.5_1[, rownames(coi)]

coi <- E10.5_2@meta.data %>% dplyr::filter(DF == "Singlet")
E10.5_2 <- E10.5_2[, rownames(coi)]

coi <- E10.5_3@meta.data %>% dplyr::filter(DF == "Singlet")
E10.5_3 <- E10.5_3[, rownames(coi)]

coi <- E11.5@meta.data %>% dplyr::filter(DF == "Singlet")
E11.5 <- E11.5[, rownames(coi)]

coi <- E12.5@meta.data %>% dplyr::filter(DF == "Singlet")
E12.5 <- E12.5[, rownames(coi)]

coi <- E13.5@meta.data %>% dplyr::filter(DF == "Singlet")
E13.5 <- E13.5[, rownames(coi)]

coi <- E14.5@meta.data %>% dplyr::filter(DF == "Singlet")
E14.5 <- E14.5[, rownames(coi)]

coi <- E17.5@meta.data %>% dplyr::filter(DF == "Singlet")
E17.5 <- E17.5[, rownames(coi)]

coi <- P1@meta.data %>% dplyr::filter(DF == "Singlet")
P1 <- P1[, rownames(coi)]

coi <- P7@meta.data %>% dplyr::filter(DF == "Singlet")
P7 <- P7[, rownames(coi)]


setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/v10/v10.2.1/data")
gexatac_merge <- readRDS("4_gexatac_merge_links.rds")
gexatac_merge[["RNA"]] <- split(gexatac_merge[["RNA"]], f = gexatac_merge$sample)

setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/")
NCC_merge <- merge(E8.5, y = c(E9.5, E10.5_1, E10.5_2, E10.5_3, E11.5, E12.5, E13.5, E14.5, E17.5, P1, P7,
                               gexatac_merge), 
                   project = "NCC_merge")

NCC_merge@meta.data$orig.ident %>% table


coi <- NCC_merge@meta.data %>% dplyr::filter(!(orig.ident %in% c("E11.5_non_EYFP", "E12.5_non_EYFP", "E14.5_non_EYFP", "E17.5_non_EYFP")))
NCC_merge <- NCC_merge[, rownames(coi)]
NCC_merge@meta.data$orig.ident %>% table

meta <- NCC_merge@meta.data$orig.ident 
NCC_merge@meta.data <- NCC_merge@meta.data %>% 
  mutate(
    sample = case_when(
      orig.ident == "E10.5" ~ "E10.5_3",  
      TRUE ~ orig.ident                 
    )
  )

NCC_merge@meta.data$sample %>% table

NCC_merge@meta.data$sample <- as.factor(NCC_merge@meta.data$sample)
NCC_merge@meta.data$sample <- factor(NCC_merge@meta.data$sample,
                                         levels = c("E8.5","E9.5","E10.5_1","E10.5_2","E10.5_3","E11.5","E11.5_EYFP","E12.5","E12.5_EYFP",
                                                    "E13.5","E14.5","E14.5_EYFP","E17.5", "E17.5_EYFP","P1","P7"))

rm(E8.5)
rm(E9.5)
rm(E10.5_1)
rm(E10.5_2)
rm(E10.5_3)
rm(E11.5)
rm(E12.5)
rm(E13.5)
rm(E14.5)
rm(E17.5)
rm(P1)
rm(P7)
rm(gexatac_merge)

NCC_merge[["RNA"]] <- JoinLayers(NCC_merge[["RNA"]])
NCC_merge[["RNA"]] <- split(NCC_merge[["RNA"]], f = NCC_merge$sample)
NCC_merge

DefaultAssay(NCC_merge)
# [1] "RNA"


setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/")
# saveRDS(NCC_merge, "NCC_merge_unintegrated.rds")
# NCC_merge <- readRDS("NCC_merge_unintegrated.rds")


NCC_merge <- ScaleData(NCC_merge, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(NCC_merge))
NCC_merge <- RunPCA(NCC_merge)
NCC_merge <- FindNeighbors(NCC_merge, dims = 1:30, reduction = "pca")
NCC_merge <- FindClusters(NCC_merge, resolution = 2, cluster.name = "unintegrated_clusters")
NCC_merge <- RunUMAP(NCC_merge, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(NCC_merge, reduction = "umap.unintegrated", group.by = c("sample"))

# RPCA
DefaultAssay(NCC_merge)
NCC_merge <- RunPCA(NCC_merge)
NCC_merge <- IntegrateLayers(
  object = NCC_merge, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca", dims = 1:50,
  verbose = FALSE
)

setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll")
NCC_merge <- RunUMAP(NCC_merge, reduction = "integrated.rpca", dims = 1:50, reduction.name = "umap.rpca")
DimPlot(NCC_merge, reduction = "umap.rpca", group.by = c("sample"), combine = FALSE)
ggsave("1_umap_RPCA_sample.pdf", height = 6, width = 6)

NCC_merge <- FindNeighbors(NCC_merge, reduction = "integrated.rpca", dims = 1:50,  k.param = 20)
NCC_merge <- FindClusters(NCC_merge, resolution = 0.8, cluster.name = "rpca_clusters")
NCC_merge@meta.data$rpca_clusters <- as.character(NCC_merge@meta.data$rpca_clusters) %>% as.numeric()
Idents(NCC_merge)
DimPlot(NCC_merge, reduction = "umap.rpca", group.by = c("rpca_clusters"), combine = F, label = T)
# NCC_merge@meta.data %>% dplyr::select(seurat_clusters, sample) %>% table
ggsave("1_umap_RPCA_clusters.pdf", height = 6, width = 6)

FeaturePlot(NCC_merge, features = c("Dlx1"), reduction = "umap.rpca")


NCC_merge <- JoinLayers(NCC_merge)
NCC_merge


# saveRDS(NCC_merge, "NCC_merge_RPCA.rds")

library(loupeR)
clusters <- select_clusters(NCC_merge)
projections <- select_projections(NCC_merge)
create_loupe(count_mat = NCC_merge@assays$RNA$counts, clusters = clusters,
             projections = projections, output_name = "NCC_merge_RNA_counts")

clusters <- select_clusters(NCC_merge)
projections <- select_projections(NCC_merge)
create_loupe(count_mat = NCC_merge@assays$RNA$data, clusters = clusters,
             projections = projections, output_name = "NCC_merge_RNA_data")



setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/")
NCC_merge <- readRDS("NCC_merge2_RPCA_velo_region.rds")
Idents(NCC_merge) <- "integrated"
DimPlot(NCC_merge, reduction = "umap.rpca", label = T)

markers <- FindAllMarkers(NCC_merge, only.pos = TRUE)
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/")
write.table(markers, "markers.txt", quote = F, sep = "\t", row.names = F)
markers <- read.table("markers.txt", header = T)

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

cluster.averages <- AverageExpression(NCC_merge, return.seurat = TRUE, assays = "RNA")
cluster.averages

DoHeatmap(cluster.averages, features = top10$gene, draw.lines = FALSE)  + NoLegend()
ggsave("heatmap_markeres_log2FC_top10_pseudobulk.png", height = 25, width = 8)



coi <- NCC_merge@meta.data %>% dplyr::filter(!seurat_clusters %in% c("27", "28"))
NCC_merge2 <- NCC_merge[, rownames(coi)]

NCC_merge2[["RNA"]] <- split(NCC_merge2[["RNA"]], f = NCC_merge2$sample)
NCC_merge2

NCC_merge2 <- IntegrateLayers(
  object = NCC_merge2, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca", dims = 1:50,
  verbose = FALSE
)


NCC_merge2 <- RunUMAP(NCC_merge2, reduction = "integrated.rpca", dims = 1:50, reduction.name = "umap.rpca")
DimPlot(NCC_merge2, reduction = "umap.rpca", group.by = c("sample"), combine = FALSE)
ggsave("1_umap_RPCA_sample.pdf", height = 6, width = 6)

NCC_merge2 <- FindNeighbors(NCC_merge2, reduction = "integrated.rpca", dims = 1:50,  k.param = 20)
NCC_merge2 <- FindClusters(NCC_merge2, resolution = 0.8, cluster.name = "rpca_clusters")
NCC_merge2@meta.data$rpca_clusters <- as.character(NCC_merge2@meta.data$rpca_clusters) %>% as.numeric()
Idents(NCC_merge2)
DimPlot(NCC_merge2, reduction = "umap.rpca", group.by = c("rpca_clusters"), combine = F, label = T)
# NCC_merge2@meta.data %>% dplyr::select(seurat_clusters, sample) %>% table
ggsave("1_umap_RPCA_clusters.pdf", height = 6, width = 6)

FeaturePlot(NCC_merge2, features = c("Cdh5"), reduction = "umap.rpca")

NCC_merge2 <- JoinLayers(NCC_merge2)
NCC_merge2

# saveRDS(NCC_merge2, "NCC_merge2_RPCA_velo.rds")
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered")
NCC_merge <- readRDS("NCC_merge2_velo.rds")
DimPlot(NCC_merge, reduction = "umap.rpca")
DimPlot(NCC_merge, reduction = "umap.rpca", group.by = "Phase")
ggsave("1_umap_RPCA_Phase.pdf", height = 6, width = 6)


Idents(NCC_merge)
NCC_merge <- FindSubCluster(NCC_merge, cluster = 0, graph.name = "RNA_snn", algorithm = 1, resolution = 0.15,
                            subcluster.name = "integrated")
DimPlot(NCC_merge, reduction = "umap.rpca", group.by = c("integrated"), combine = F, label = T)

Idents(NCC_merge) <- "integrated"
NCC_merge@meta.data$integrated <- NCC_merge@meta.data$integrated %>%
  str_replace_all(pattern = "0_0", replacement = "0") %>%
  str_replace_all(pattern = "0_1", replacement = "27") 

NCC_merge@meta.data$integrated <- as.factor(NCC_merge@meta.data$integrated)
NCC_merge@meta.data$integrated <- factor(NCC_merge@meta.data$integrated,
                                              levels = 0:27)
DimPlot(NCC_merge, reduction = "umap.rpca", group.by = c("integrated"), combine = F, label = T)
ggsave("1_umap_RPCA_integrated.pdf", height = 6, width = 6)


NCC_merge@meta.data$region <- ifelse(NCC_merge@meta.data$integrated %in% c("11","15", "5","8","13", "6"), "Pharyngeal",
                               ifelse(NCC_merge@meta.data$integrated %in% c("27"), "SMC_GA", 
                                      ifelse(NCC_merge@meta.data$integrated %in% c("4"), "SMC_DA", 
                                      ifelse(NCC_merge@meta.data$integrated %in% c("23"), "SMC_CA", 
                                             ifelse(NCC_merge@meta.data$integrated %in% c("0","7","14"), "Transitional", 
                                                    ifelse(NCC_merge@meta.data$integrated %in% c("16"), "AP septum", 
                                                           ifelse(NCC_merge@meta.data$integrated %in% c("9","12","25"), "Glial", 
                                                                  ifelse(NCC_merge@meta.data$integrated %in% c("17","20"), "Neural", 
                                                                         ifelse(NCC_merge@meta.data$integrated %in% c("24"), "Melanocyte", 
                                                                                ifelse(NCC_merge@meta.data$integrated %in% c("26"), "Cardiomyocyte", 
                                                                                       ifelse(NCC_merge@meta.data$integrated %in% c("1","3","19"), "Neural tube",      
                                                           ifelse(NCC_merge@meta.data$integrated %in% c("2","21","22","10","18"), "Cardiac Cushion", NA))))))))))))

NCC_merge@meta.data$region <- as.factor(NCC_merge@meta.data$region)
NCC_merge@meta.data$region<- factor(NCC_merge@meta.data$region,
                              levels =  c("Pharyngeal","Transitional","Cardiac Cushion", "AP septum", "SMC_GA","SMC_DA","SMC_CA",
                                          "Neural","Glial","Melanocyte", "Neural tube", "Cardiomyocyte"
                                          ))
DimPlot(NCC_merge, reduction = "umap.rpca", group.by = "region",label = F) 

ggsave("1_umap.NCC_merge_region.pdf", height = 6, width = 6)

saveRDS(NCC_merge, "NCC_merge2_RPCA_velo_region.rds")

setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/")
NCC_merge <- readRDS("NCC_merge2_RPCA_velo_region.rds")

NCC_merge@meta.data <- NCC_merge@meta.data %>%
  mutate(stage = case_when(
    sample %in% c("E10.5_1", "E10.5_2", "E10.5_3") ~ "E10.5",
    sample %in% c("E11.5", "E11.5_EYFP")           ~ "E11.5",
    sample %in% c("E12.5", "E12.5_EYFP")           ~ "E12.5",
    sample %in% c("E14.5", "E14.5_EYFP")           ~ "E14.5",
    sample %in% c("E17.5", "E17.5_EYFP")           ~ "E17.5",
    TRUE                                           ~ sample  
  ))
NCC_merge@meta.data$stage <- as.factor(NCC_merge@meta.data$stage) 
NCC_merge@meta.data$stage <- factor(NCC_merge@meta.data$stage,
                                     levels = c("E8.5","E9.5","E10.5","E11.5",
                                                "E12.5", "E13.5","E14.5","E17.5","P1", "P7")) 
DimPlot(NCC_merge, reduction = "umap.rpca", group.by = "region2",label = F,
        split.by = "stage", ncol = 5)
ggsave("1_umap_RPCA_region2_facet_stage.pdf", height = 6, width = 12)

stages <- unique(NCC_merge$stage)
colors <- rev(rainbow(length(stages)))

DimPlot(NCC_merge, reduction = "umap.rpca", group.by = "stage", label = FALSE) +
  scale_color_manual(values = setNames(colors, stages))
ggsave("1_umap_RPCA_stage.pdf", height = 6, width = 6)




df <- FetchData(NCC_merge, vars = c("stage", "region2", "Phase"))

colnames(df) <- c("stage", "region", "cell_cycle")

df_summary <- df %>%
  group_by(stage, region) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(stage) %>%
  mutate(prop = n / sum(n))  # 割合に変換
ggplot(df_summary, aes(x = stage, y = prop, fill = region)) +
  geom_bar(stat = "identity", position = "stack") +
  ylab("Proportion") +
  theme_minimal() +
  theme(
    legend.position = "right",         
    legend.box = "vertical",       
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  guides(fill = guide_legend(ncol = 1)) 
ggsave("bar_propotion_stage_region2.pdf", width = 6, height = 4)


df <- FetchData(NCC_merge, vars = c("stage", "integrated", "Phase"))

colnames(df) <- c("stage", "cluster", "cell_cycle")

df_summary <- df %>%
  group_by(cluster, cell_cycle) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n))  

ggplot(df_summary, aes(x = cluster, y = prop, fill = cell_cycle)) +
  geom_bar(stat = "identity", position = "stack") +
  ylab("Proportion") +
  theme_minimal() +
  theme(
    legend.position = "right",        
    legend.box = "vertical",          
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  guides(fill = guide_legend(ncol = 1)) 
ggsave("bar_propotion_stage_cellcycle.pdf", width = 6, height = 4)






coi <- NCC_merge@meta.data %>% dplyr::filter(integrated %in% c("21","2","22","10","18","16",
                                                               "23","27","4","0","7","14",
                                                               "6","13","11","15","5","8","13"))
NCC_merge <- NCC_merge[, rownames(coi)]
DimPlot(NCC_merge, reduction = "umap.rpca", group.by = "integrated",label = F)
ggsave("1_umap_RPCA_subset.pdf", height = 6, width = 6)


 




library(loupeR)
clusters <- select_clusters(NCC_merge2)
projections <- select_projections(NCC_merge2)
create_loupe(count_mat = NCC_merge2@assays$RNA$counts, clusters = clusters,
             projections = projections, output_name = "NCC_merge2_RNA_counts")

clusters <- select_clusters(NCC_merge)
projections <- select_projections(NCC_merge)
create_loupe(count_mat = NCC_merge@assays$RNA$data, clusters = clusters,
             projections = projections, output_name = "NCC_merge2_RNA_data")







