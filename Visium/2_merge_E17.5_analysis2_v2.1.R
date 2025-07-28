#merge_E17.5_analysis

library(hdf5r)
library(Seurat)

library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
library(SeuratObject)

setwd("/Volumes/Neon-HD/Visium/200818_original_data/SpaceRanger_EYFP/3_E17_5_1_EYFP")
E17.5_1 <- Load10X_Spatial(data.dir = "./outs",
                           filename = "filtered_feature_bc_matrix.h5",
                           slice = "slice3")



E17.5_1@meta.data$orig.ident <- "E17.5_1"
E17.5_1@meta.data$orig.ident <- as.factor(E17.5_1@meta.data$orig.ident)
E17.5_1@meta.data$slice <- c("E17.5_1")


img <- Seurat::Read10X_Image(image.dir = "./outs/spatial")

Seurat::DefaultAssay(object = img) <- 'Spatial'

E17.5_1

SpatialDimPlot(E17.5_1,
               cells.highlight = WhichCells(E17.5_1, 
                                            expression = 
                                              slice3_imagecol < 400))

dim(E17.5_1)
# [1] 24422   836
E17.5_1 <- subset(E17.5_1, slice3_imagecol < 400,invert = F)
dim(E17.5_1)
# [1] 24422   835


setwd("/Volumes/Neon-HD/Visium/200818_original_data/SpaceRanger_EYFP/4_E17_5_2_EYFP")
E17.5_2 <- Load10X_Spatial(data.dir = "./outs",
                           filename = "filtered_feature_bc_matrix.h5",
                           slice = "slice4")


E17.5_2@meta.data$orig.ident <- "E17.5_2"
E17.5_2@meta.data$orig.ident <- as.factor(E17.5_2@meta.data$orig.ident)
E17.5_2@meta.data$slice <- c("E17.5_2")

img <- Seurat::Read10X_Image(image.dir = "./outs/spatial")

Seurat::DefaultAssay(object = img) <- 'Spatial'

E17.5_2
# An object of class Seurat 
# 24422 features across 595 samples within 1 assay 
# Active assay: Spatial (24422 features, 0 variable features))
#これで除外すべき細胞を判定できる
SpatialDimPlot(E17.5_2,
               cells.highlight = WhichCells(E17.5_2, 
                                            expression =  
                                              slice4_imagecol < 400 & 
                                              slice4_imagecol > 170))
dim(E17.5_2)
# [1] 24422   595
E17.5_2 <- subset(E17.5_2, slice4_imagecol < 400 & slice4_imagecol > 170,invert = F)
dim(E17.5_2)
# [1] 24422   593


setwd("~/Downloads/Visium/analysis/merge_E17.5")

E17.5_1 <- SCTransform(E17.5_1, assay = "Spatial", verbose = FALSE)
E17.5_2 <- SCTransform(E17.5_2, assay = "Spatial", verbose = FALSE)

#merge
E17.5.merge <- merge(E17.5_1, E17.5_2)

########################################
#保存
#saveRDS(E17.5.merge, "E17.5.merge.rds")
########################################



DefaultAssay(E17.5.merge) <- "SCT"

い

VariableFeatures(E17.5.merge) <- c(VariableFeatures(E17.5_1), VariableFeatures(E17.5_2))
VariableFeatures(E17.5.merge) %>% length()


E17.5.merge <- RunPCA(E17.5.merge, verbose = FALSE)
E17.5.merge <- FindNeighbors(E17.5.merge, dims = 1:30, k.param = 4)

E17.5.merge <- FindClusters(E17.5.merge, verbose = FALSE, resolution = 1)
E17.5.merge <- RunUMAP(E17.5.merge, dims = 1:30)
#spatial cluster
SpatialDimPlot(E17.5.merge, crop=T, pt.size.factor = 3, label = T, label.box = F,label.size = 3)
ggsave("1_spatial_cluster_label.pdf", height = 6, width = 8)
SpatialDimPlot(E17.5.merge, crop=T, pt.size.factor = 3, label = F)
ggsave("1_spatial_cluster.pdf", height = 6, width = 8)
#umap cluster
DimPlot(E17.5.merge, reduction = "umap", group.by = c("ident", "orig.ident"))
ggsave("1_umap_cluster.pdf", height = 6, width = 16)
DimPlot(E17.5.merge, reduction = "umap", group.by = c("ident", "orig.ident"), 
        label =T)
ggsave("1_umap_cluster_label.pdf", height = 6, width = 16)

SpatialDimPlot(E17.5.merge, crop=T, pt.size.factor = 3, label = T, 
               label.box = F,label.size = 3, image = "slice3")
ggsave("1_spatial_cluster_slice3_label.pdf", height = 6, width = 8)
SpatialDimPlot(E17.5.merge, crop=T, pt.size.factor = 3, label = T, 
               label.box = F,label.size = 3, image = "slice4")
ggsave("1_spatial_cluster_slice4_label.pdf", height = 6, width = 8)

DimPlot(E17.5.merge, reduction = "umap", group.by = c("orig.ident"), 
        label =F)
ggsave("1_dimplot_batch.pdf", height = 6, width = 8)

SpatialDimPlot(E17.5.merge, 
               cells.highlight = 
                 CellsByIdentities(object = E17.5.merge, 
                                   idents = c(5,8,9)), 
               image = "slice4",
               facet.highlight = TRUE, ncol = 3)

########################################
#保存
setwd("~/Downloads/Visium/analysis/merge_E17.5/201117_scRNA-seq_prediciton/")
#saveRDS(E17.5.merge, "E17.5.merge2.rds")
########################################

setwd("~/Downloads/Visium/analysis/merge_E17.5/201117_scRNA-seq_prediciton/")
VlnPlot(E17.5.merge, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
ggsave("2_nCount_violin.pdf", height = 6, width = 8)

SpatialFeaturePlot(E17.5.merge, features = "nCount_Spatial",pt.size.factor=3)
ggsave("2_nCount_spatial.pdf", height = 6, width = 8)  


SpatialFeaturePlot(E17.5.merge, features = "nCount_Spatial",pt.size.factor=3,
                   images = "slice3") +theme(legend.position = "right")
ggsave("2_nCount_spatial_slice3.pdf", height = 6, width = 8)  
SpatialFeaturePlot(E17.5.merge, features = "nCount_Spatial",pt.size.factor=3,
                   images = "slice4") +theme(legend.position = "right")
ggsave("2_nCount_spatial_slice4.pdf", height = 6, width = 8)  



setwd("~/Downloads/Visium/analysis/merge_E17.5/201117_scRNA-seq_prediciton")
E17.5.merge <- readRDS("E17.5.merge2.rds")

df <- E17.5.merge@assays$SCT@counts %>% as.data.frame()
write.table(df, "E17.5.merge_SCT.txt", sep = "\t", quote =F, col.names = NA)
meta <- E17.5.merge@meta.data %>% as.data.frame()
write.table(meta, "E17.5.merge_meta.txt", sep = "\t", quote =F, col.names = NA)



setwd("~/Downloads/Visium/analysis/merge_E17.5/201117_scRNA-seq_prediciton")
E17.5.merge <- readRDS("E17.5.merge2.rds")

DimPlot(E17.5.merge, reduction = "umap", 
        label = T) + NoLegend()

levels(E17.5.merge)
# [1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15"
# [17] "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30" "31"
SpatialDimPlot(E17.5.merge, images  = "slice3",
               cells.highlight = CellsByIdentities(object = E17.5.merge, 
                                                   idents = c(0:31)),
               facet.highlight = TRUE, ncol = 6)
SpatialDimPlot(E17.5.merge, images  = "slice4",
               cells.highlight = CellsByIdentities(object = E17.5.merge, 
                                                   idents = c(0:31)),
               facet.highlight = TRUE, ncol = 6)
# 10x10で保存
new.cluster.ids <- 
  c("0.Left ventricle 1", 
    "1.Interventricular septum 1",
    "2.Atrium/Epicardium", 
    "3.Right ventricle 1",
    "4.Right atrium 1", 
    "5.Aorta/Pulmonary artery 1",
    "6.Right atrium 2", 
    "7.Right ventricle 2",
    "8.Mediastinal tissue", 
    "9.Periaortic tissue 1",
    "10.Trachea1", 
    "11.Interventricular septum 2",
    "12.Aorta/Pulmonary artery 2",
    "13.Left ventricle 2",
    "14.Right venticle 3", 
    "15.Left atrium", 
    "16.Periaortic tissue 2",
    "17.Ventricle/Epicardium", 
    "18.Ventricle/Endocardial cushion",
    "19.Left ventricle 3",
    "20.Left ventricle 4", 
    "21.Left ventricle/Coronary artery", 
    "22.Ganglion", 
    "23.Interventricular septum 3",
    "24.Left ventricle 5",
    "25.Trachea 2", 
    "26.Right atrium/Superior vena cava", 
    "27.Trachea 3",
    "28.Left ventricle 6",
    "29.Skeletal muscle", 
    "30.Right atrium/Sinoatrial node", 
    "31.Pulmonary valve")
names(new.cluster.ids) <- levels(E17.5.merge)
E17.5.merge <- RenameIdents(E17.5.merge, new.cluster.ids)
DimPlot(E17.5.merge, reduction = "umap", label = TRUE, pt.size = 0.5)



ggsave("2_umap_annotation.pdf", height = 6, width = 8)  
SpatialDimPlot(E17.5.merge, label = TRUE, label.size = 3)




setwd("~/Downloads/Visium/analysis/merge_E17.5/201117_scRNA-seq_prediciton/goi")
goi <- c("Cd68", "Csf1r", "Adgre1","Csf1", "Ccl2", "Ccr2",
         "Sox9", "Scx", "Col2a1",
         "Acta2", "Tagln", "Myh11", "Cnn1",
         "Klf2", "Klf4", "Foxs1", "Osr1", "Sost",
         "Abcc9", "Pdgfrb", "Cspg4", "Dlk1", "Myh10",
         "Tnnt2", "Myh7", "Myl7",
         "Cdh5", "Pecam1", "Prox1", "Lyve1",
         "Fstl1", "Postn", "Bgn", "Dcn", "Mgp",
         "Sox10", "S100b", "Ngfr",
         "Notch1", "Notch2", "Notch3", "Dll4", "Jag1",
         "EYFP")

for(i in 1 : length(goi)){
  SpatialFeaturePlot(E17.5.merge, features = goi[i], pt.size.factor=3)
  ggsave(paste(goi[i], ".pdf", sep =""), height = 6, width = 8)
}



################################################
#Identification of Spatially Variable Features
################################################
setwd("~/Downloads/Visium/analysis/merge_E17.5/201117_scRNA-seq_prediciton")
E17.5.merge <- readRDS("E17.5.merge2.rds")

all.markers <- FindAllMarkers(E17.5.merge, only.pos = TRUE, min.pct = 0.25, 
                              logfc.threshold = 0.25)
setwd("~/Downloads/Visium/analysis/merge_E17.5/201117_scRNA-seq_prediciton/markers")
write.table(all.markers, "all.markers.txt", quote = F, sep = "\t", col.names = NA)

setwd("~/Downloads/Visium/analysis/merge_E17.5/201117_scRNA-seq_prediciton/markers")
all.markers <- read.table("all.markers.txt", header = T)
all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

all.genes <- rownames(E17.5.merge)
E17.5.merge <- ScaleData(E17.5.merge, features = all.genes)
top <-  all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(E17.5.merge, features = top$gene, label = T, size = 3) + NoLegend()
ggsave("top10.png", height = 30, width = 15)

SpatialFeaturePlot(object = E17.5.merge, 
                   features = rownames(all.markers)[1:3], alpha = c(0.1, 1), ncol = 3)



setwd("/Volumes/iwaseakiyasu/Downloads/scRNA-seq/20200403_Seuratv3_based2/Seurat/dim8_k4/subset_cluster4")
Seuratcds <- readRDS("C1.rds")
DimPlot(Seuratcds, reduction = "umap", label = T) 

Seuratcds@meta.data$celltype <- Seuratcds@meta.data$sub.cluster %>% as.factor()
Seuratcds@meta.data %>% head


Seuratcds <- subset(Seuratcds, stage %in% c("E17.5"))
table(Seuratcds@meta.data$stage)



########################################
#VisiumからEYFP > 0
########################################
setwd("/Volumes/iwaseakiyasu/Downloads/Visium/analysis/merge_E17.5/201117_scRNA-seq_prediciton")
E17.5.merge <- readRDS("E17.5.merge2.rds")
SpatialDimPlot(E17.5.merge, crop=T, pt.size.factor = 3, label = T, label.box = F,label.size = 3)

counts_sct <- E17.5.merge$SCT@counts %>% as.data.frame() %>% t() %>% as.data.frame()
dim(counts_sct)
# [1]   1428 15729
subset(counts_sct, EYFP >  0) %>% dim
# [1]   536 15729
eyfp_spot <- subset(counts_sct, EYFP >  0) %>% rownames()

E17.5.merge_EYFP <- E17.5.merge[, eyfp_spot]
dim(E17.5.merge_EYFP)
# [1] 15729   536
SpatialDimPlot(E17.5.merge_EYFP, crop=T, pt.size.factor = 3, label = T, label.box = F,label.size = 3)
ggsave("1_spatial_cluster_EYFP_label.pdf", height = 6, width = 8)

E17.5.merge_EYFP <- SCTransform(E17.5.merge_EYFP, 
                                assay = "Spatial",
                                variable.features.n = rownames(E17.5.merge_EYFP) %>% length) %>%
  RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30)
DimPlot(E17.5.merge_EYFP, reduction = "umap", group.by = c("ident", "orig.ident"))


Seuratcds2 <- Seuratcds 
Seuratcds2@meta.data$celltype
Seuratcds2 <- SCTransform(Seuratcds2, 
                          assay = "RNA",
                          variable.features.n = rownames(Seuratcds) %>% length) %>%
  RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:4) %>% 
  DimPlot(Seuratcds2, reduction = "umap", group.by = c("ident", "stage"))
Seuratcds2@meta.data$celltype

anchors <- FindTransferAnchors(reference = Seuratcds2, 
                               query = E17.5.merge_EYFP, normalization.method = "SCT")


predictions.assay <- TransferData(anchorset = anchors, 
                                  refdata = Seuratcds2$celltype, prediction.assay = TRUE, 
                                  weight.reduction = E17.5.merge_EYFP[["pca"]])

library(stringr)
rownames(predictions.assay@data) <- rownames(predictions.assay@data) %>% str_replace(pattern = "-", replacement = "_")
E17.5.merge_EYFP[["predictions"]] <- predictions.assay
DefaultAssay(E17.5.merge_EYFP) <- "predictions"


##############
#可視化と自動保存
##############
setwd("/Volumes/iwaseakiyasu/Downloads/Visium/analysis/merge_E17.5/201117_scRNA-seq_prediciton/EYFP_prediction2/")

celltypelist <- Seuratcds2@meta.data %>% distinct(celltype)
celltypelist$celltype %>% as.character() 
celltypelist<- celltypelist$celltype %>% as.character() %>% as.list()



library(RColorBrewer)
for(i in 1:(celltypelist %>% length())){
  celltypelist[i]
  SpatialFeaturePlot(E17.5.merge_EYFP, 
                     features = celltypelist[i] %>% as.character(), 
                     pt.size.factor = 3, ncol = 2, crop = TRUE,
                     image = "slice3") + 
    scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")), 
                         limits = c(0, 
                                    E17.5.merge_EYFP@assays$predictions@data[i,] %>%
                                      max()))
  ggsave(paste("spatial_celltype_estimation_slice3_", celltypelist[i] %>% as.character(), ".pdf", sep=""),
         height = 6, width = 8)
  SpatialFeaturePlot(E17.5.merge_EYFP, 
                     features = celltypelist[i] %>% as.character(), 
                     pt.size.factor = 3, ncol = 2, crop = TRUE,
                     image = "slice4") + 
    scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")), 
                         limits = c(0, 
                                    E17.5.merge_EYFP@assays$predictions@data[i,] %>%
                                      max()))
  ggsave(paste("spatial_celltype_estimation_slice4_", celltypelist[i] %>% as.character(), ".pdf", sep=""),
         height = 6, width = 8)
  
}




########################################
#保存
setwd("~/Downloads/Visium/analysis/merge_E17.5/201117_scRNA-seq_prediciton/")
#saveRDS(E17.5.merge_EYFP, "E17.5.merge_EYFP.rds")
########################################


########################################
#VisiumからEYFP > 0
########################################
setwd("~/Downloads/Visium/analysis/merge_E17.5/201117_scRNA-seq_prediciton")
E17.5.merge <- readRDS("E17.5.merge2.rds")
SpatialDimPlot(E17.5.merge, crop=T, pt.size.factor = 3, label = T, label.box = F,label.size = 3)

counts_sct <- E17.5.merge$SCT@counts %>% as.data.frame() %>% t() %>% as.data.frame()
dim(counts_sct)
# [1]   1428 15729
subset(counts_sct, EYFP >  0) %>% dim
# [1]   536 15729
eyfp_spot <- subset(counts_sct, EYFP >  0) %>% rownames()

E17.5.merge_EYFP <- E17.5.merge[, eyfp_spot]
dim(E17.5.merge_EYFP)
# [1] 15729   536
setwd("~/Downloads/Visium/analysis/merge_E17.5/201117_scRNA-seq_prediciton/goi_EYFP/")
# setwd("./goi_EYFP/")
SpatialDimPlot(E17.5.merge_EYFP, crop=T, pt.size.factor = 3, label = T, 
               label.box = F,label.size = 3)
goi <- c("Cd68", "Csf1r", "Adgre1","Csf1", "Ccl2", "Ccr2",
         "Sox9", "Scx", "Col2a1",
         "Acta2", "Tagln", "Myh11", "Cnn1",
         "Klf2", "Klf4", "Foxs1", "Osr1", "Sost",
         "Abcc9", "Pdgfrb", "Cspg4", "Dlk1", "Myh10",
         "Tnnt2", "Myh7", "Myl7",
         "Cdh5", "Pecam1", "Prox1", "Lyve1",
         "Fstl1", "Postn", "Bgn", "Dcn", "Mgp",
         "Sox10", "S100b", "Ngfr",
         "Notch1", "Notch2", "Notch3", "Dll4", "Jag1",
         "EYFP")
goi <- c("Twist1")
goi <- c("Acan", "Fmod")
goi <- c("Hmga2", "Mki67", "Top2a", "Ccnb2")
goi <- c("Hapln1", "Vcan", "Tbx20", "Postn", "Dcn", "Lum", "Bgn", "Fmod",
         "Fbln5", "Myom1", "Lmod1", "Jag1" , "Notch3")
goi <- c("Plp1")
goi <- c("L1cam", "Gjc3")
goi <- c("Meox1", "Col3a1", "Dcn")


for(i in 1 : length(goi)){
  SpatialFeaturePlot(E17.5.merge_EYFP, features = goi[i], pt.size.factor=3) +
  ggsave(paste(goi[i], ".pdf", sep =""), height = 6, width = 8)
}
###







