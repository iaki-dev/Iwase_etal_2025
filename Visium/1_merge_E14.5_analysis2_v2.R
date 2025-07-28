r#merge_E14.5_analysis2_v2

library(hdf5r)
library(Seurat)

library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)



setwd("/Volumes/Neon-HD/Visium/200818_original_data/SpaceRanger_EYFP/1_E14_5_1_EYFP")
E14.5_1 <- Load10X_Spatial(data.dir = "./outs",
                           filename = "filtered_feature_bc_matrix.h5",
                           slice = "slice1")

E14.5_1@meta.data$orig.ident <- "E14.5_1"
E14.5_1@meta.data$orig.ident <- as.factor(E14.5_1@meta.data$orig.ident)
E14.5_1@meta.data$slice <- c("E14.5_1")


img <- Seurat::Read10X_Image(image.dir = "./outs/spatial")

Seurat::DefaultAssay(object = img) <- 'Spatial'

E14.5_1

SpatialDimPlot(E14.5_1,
               cells.highlight = WhichCells(E14.5_1, 
                                            expression = slice1_imagerow > 10 | 
                                              slice1_imagecol < 170))



setwd("/Volumes/Neon-HD/Visium/200818_original_data/SpaceRanger_EYFP/2_E14_5_2_EYFP")
E14.5_2 <- Load10X_Spatial(data.dir = "./outs",
                           filename = "filtered_feature_bc_matrix.h5",
                           slice = "slice2")

E14.5_2@meta.data$orig.ident <- "E14.5_2"
E14.5_2@meta.data$orig.ident <- as.factor(E14.5_2@meta.data$orig.ident)
E14.5_2@meta.data$slice <- c("E14.5_2")


img <- Seurat::Read10X_Image(image.dir = "./outs/spatial")

Seurat::DefaultAssay(object = img) <- 'Spatial'

E14.5_2

SpatialDimPlot(E14.5_2,
               cells.highlight = WhichCells(E14.5_2, 
                                            expression =  
                                              slice2_imagecol < 400))
dim(E14.5_2)
# [1] 24422   291
E14.5_2 <- subset(E14.5_2, slice2_imagecol < 400,invert = F)
dim(E14.5_2)
# [1] 24422   29



setwd("~/Downloads/Visium/analysis/merge_E14.5")

E14.5_1 <- SCTransform(E14.5_1, assay = "Spatial", verbose = FALSE)
E14.5_2 <- SCTransform(E14.5_2, assay = "Spatial", verbose = FALSE)

#merge
E14.5.merge <- merge(E14.5_1, E14.5_2)

########################################
#保存
setwd("~/Downloads/Visium/analysis/merge_E14.5/201117_scRNA-seq_prediction")
#saveRDS(E14.5.merge, "E14.5.merge.rds")
########################################

setwd("~/Downloads/Visium/analysis/merge_E14.5")
E14.5.merge <- readRDS("E14.5.merge.rds")

DefaultAssay(E14.5.merge) <- "SCT"

VariableFeatures(E14.5.merge) <- c(VariableFeatures(E14.5_1), VariableFeatures(E14.5_2))
VariableFeatures(E14.5.merge) %>% length()

E14.5.merge <- RunPCA(E14.5.merge, verbose = FALSE)
E14.5.merge <- FindNeighbors(E14.5.merge, dims = 1:30, k.param = 4)

E14.5.merge <- FindClusters(E14.5.merge, verbose = FALSE, resolution = 1)
E14.5.merge <- RunUMAP(E14.5.merge, dims = 1:30)
#spatial cluster
SpatialDimPlot(E14.5.merge, crop=T, pt.size.factor = 4, label = T, label.box = F,label.size = 3)
ggsave("1_spatial_cluster_label.pdf", height = 6, width = 8)
SpatialDimPlot(E14.5.merge, crop=T, pt.size.factor = 4, label = F)
ggsave("1_spatial_cluster.pdf", height = 6, width = 8)
#umap cluster
DimPlot(E14.5.merge, reduction = "umap", group.by = c("ident", "orig.ident"))
ggsave("1_umap_cluster.pdf", height = 6, width = 16)
DimPlot(E14.5.merge, reduction = "umap", group.by = c("ident", "orig.ident"), 
        label =T)
ggsave("1_umap_cluster_label.pdf", height = 6, width = 16)

SpatialDimPlot(E14.5.merge, crop=T, pt.size.factor = 4, label = T, 
               label.box = F,label.size = 3, image = "slice1")
ggsave("1_spatial_cluster_slice1_label.pdf", height = 6, width = 8)
SpatialDimPlot(E14.5.merge, crop=T, pt.size.factor = 4, label = T, 
               label.box = F,label.size = 3, image = "slice2")
ggsave("1_spatial_cluster_slice2_label.pdf", height = 6, width = 8)

DimPlot(E14.5.merge, reduction = "umap", group.by = c("orig.ident"), 
        label =F)
ggsave("1_dimplot_batch.pdf", height = 6, width = 8)


########################################
#保存
setwd("~/Downloads/Visium/analysis/merge_E14.5/201117_scRNA-seq_prediction/")
#saveRDS(E14.5.merge, "E14.5.merge2.rds")
########################################
setwd("~/Downloads/Visium/analysis/merge_E14.5/201117_scRNA-seq_prediction/")
E14.5.merge <- readRDS("E14.5.merge2.rds")

DimPlot(E14.5.merge, reduction = "umap",
        label = T) + NoLegend()

VlnPlot(E14.5.merge, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
ggsave("2_nCount_violin.pdf", height = 6, width = 8)

SpatialFeaturePlot(E14.5.merge, features = "nCount_Spatial",pt.size.factor=4,
                   images = "slice1") +theme(legend.position = "right")
ggsave("2_nCount_spatial_slice1.pdf", height = 6, width = 8)  
SpatialFeaturePlot(E14.5.merge, features = "nCount_Spatial",pt.size.factor=4,
                   images = "slice2") +theme(legend.position = "right")
ggsave("2_nCount_spatial_slice.pdf", height = 6, width = 8)  


setwd("~/Downloads/Visium/analysis/merge_E14.5/201117_scRNA-seq_prediction")
df <- E14.5.merge@assays$SCT@counts %>% as.data.frame()
write.table(df, "E14.5.merge_SCT.txt", sep = "\t", quote =F, col.names = NA)
meta <- E14.5.merge@meta.data %>% as.data.frame()
write.table(meta, "E14.5.merge_meta.txt", sep = "\t", quote =F, col.names = NA)



levels(E14.5.merge)
# [1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15"
# [17] "16" "17" "18" "19" "20" "21" "22"
SpatialDimPlot(E14.5.merge, images  = "slice1",
               cells.highlight = CellsByIdentities(object = E14.5.merge, 
                                                   idents = c(0:22)),
               facet.highlight = TRUE, ncol = 6)
SpatialDimPlot(E14.5.merge, images  = "slice2",
               cells.highlight = CellsByIdentities(object = E14.5.merge, 
                                                   idents = c(0:22)),
               facet.highlight = TRUE, ncol = 6)
# 10x10で保存


new.cluster.ids <- 
  c("0.Ventricle/Epicardium",
    "1.Right venticle 1", 
    "2.Left venticle 1", 
    "3.Aorta/Pulmonary artery",
    "4.Semilunar valve/Endocardial cushion ",
    "5.Atrium", 
    "6.Interventricular septum", 
    "7.Mediastinal tissue 1", 
    "8.Trachea 1", 
    "9.Aorticopulmonary septum", 
    "10.Ventricle/Endocardial cushion",
    "11.Left venticle 2", 
    "12.Right atrium", 
    "13.Right ventricle 2", 
    "14.Left ventricle 2",
    "15.Mediastinal tissue 2",
    "16.Left atrium", 
    "17.Right ventricle 3", 
    "18.Right ventricle 4",
    "19.Esophagus", 
    "20.Right ventricle 5", 
    "21.Mediastinal tissue 3",
    "22.Trachea 2")
names(new.cluster.ids) <- levels(E14.5.merge)
E14.5.merge <- RenameIdents(E14.5.merge, new.cluster.ids)
DimPlot(E14.5.merge, reduction = "umap", label = TRUE, pt.size = 0.5,
        repel = T) + NoLegend()
ggsave("2_umap_annotation.pdf", height = 6, width = 8)  

SpatialDimPlot(E14.5.merge, crop=T, pt.size.factor = 4, label = T, 
               images = "slice1",
               label.box = F,label.size = 3) 
SpatialDimPlot(E14.5.merge, crop=T, pt.size.factor = 4, label = T, 
               images = "slice2",
               label.box = F,label.size = 3) 
SpatialDimPlot(E14.5.merge, crop=T, pt.size.factor = 4, label = F, 
               label.box = F,label.size = 7,
               ncol = 2) + NoLegend()
# 4x6で保存
              






setwd("~/Downloads/Visium/analysis/merge_E14.5/201117_scRNA-seq_prediction/goi")
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
goi <- c("Lmod1", "Tnfsf15", "Aim1", "Nell1", "Mab21l2", "Sptssb", "Hoxd3")
goi <- c("Shisa2", "Ociad2", "Lmo7", "Myl2", "Ldha", "Myh7", "Nupr1", "Cd200", "Kcnj3")
setwd("~/Downloads/Visium/analysis/merge_E14.5/201117_scRNA-seq_prediction/goi_tmp")
goi <- c("Apc","Axin1","Ruvbl2","Snai2","Stk3","Cthrc1","Nkd2","Bicc1","Wwtr1","Scyl2","Actg1","Fzd9","Lef1","Smad7","Dlg2","Sulf2")
goi <- c("Scx","Meox1","Sall3", "Snai2","Sox6", "Sox9")
for(i in 1 : length(goi)){
  SpatialFeaturePlot(E14.5.merge, features = goi[i], pt.size.factor=4, images = "slice2")
  ggsave(paste(goi[i], ".pdf", sep =""), height = 6, width = 8)
}




load("/Users/iwaseakiyasu/Downloads/scRNA-seq/20200403_Seuratv3_based2/Seurat/dim8_k4/R3.RData")
DimPlot(Seuratcds, reduction = "umap", label = T) 

Seuratcds@meta.data$celltype <- Seuratcds@meta.data$seurat_clusters %>% as.factor()
Seuratcds@meta.data %>% head

Seuratcds <- subset(Seuratcds, stage %in% c("E14.5"))
table(Seuratcds@meta.data$stage)


########################################
#VisiumからEYFP > 0
########################################
setwd("~/Downloads/Visium/analysis/merge_E14.5/201117_scRNA-seq_prediction/")
E14.5.merge <- readRDS("E14.5.merge2.rds")
SpatialDimPlot(E14.5.merge, crop=T, pt.size.factor = 4, label = T, label.box = F,label.size = 3)

counts_sct <- E14.5.merge$SCT@counts %>% as.data.frame() %>% t() %>% as.data.frame()
dim(counts_sct)
# [1]   702 15076
subset(counts_sct, EYFP >  0) %>% dim
# [1]   313 15076
eyfp_spot <- subset(counts_sct, EYFP >  0) %>% rownames()

E14.5.merge_EYFP <- E14.5.merge[, eyfp_spot]
dim(E14.5.merge_EYFP)
# [1] 15076   313
SpatialDimPlot(E14.5.merge_EYFP, crop=T, pt.size.factor = 4, label = T, label.box = F,label.size = 3)
ggsave("1_spatial_cluster_EYFP_label.pdf", height = 6, width = 8)



E14.5.merge_EYFP <- SCTransform(E14.5.merge_EYFP, 
                                assay = "Spatial",
                                variable.features.n = rownames(E14.5.merge_EYFP) %>% length) %>%
  RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30)
DimPlot(E14.5.merge_EYFP, reduction = "umap", group.by = c("ident", "orig.ident"))

Seuratcds2 <- Seuratcds 
Seuratcds2 <- SCTransform(Seuratcds2, 
                          assay = "RNA",
                          variable.features.n = rownames(Seuratcds) %>% length) %>%
  RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:4) %>% 
  DimPlot(Seuratcds2, reduction = "umap", group.by = c("ident", "stage"))

anchors <- FindTransferAnchors(reference = Seuratcds2, 
                               query = E14.5.merge_EYFP, normalization.method = "SCT")

# dimを指定したversion
predictions.assay <- TransferData(anchorset = anchors, 
                                  refdata = Seuratcds2$celltype, prediction.assay = TRUE, 
                                  weight.reduction = E14.5.merge_EYFP[["pca"]], dims = 1:8)

E14.5.merge_EYFP[["predictions"]] <- predictions.assay
DefaultAssay(E14.5.merge_EYFP) <- "predictions"


##############
# visualization
##############
setwd("~/Downloads/Visium/analysis/merge_E14.5/201117_scRNA-seq_prediction")

celltypelist <- Seuratcds2@meta.data %>% distinct(celltype)
celltypelist$celltype %>% as.character() 
celltypelist<- celltypelist$celltype %>% as.character() %>% as.list()


for(i in 1:(celltypelist %>% length())){
  celltypelist[i]
  SpatialFeaturePlot(E14.5.merge_EYFP, 
                   features = celltypelist[i] %>% as.character(), 
                   pt.size.factor = 3, ncol = 2, crop = TRUE,
                   image = "slice1") + 
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")), 
                       limits = c(0, 
                                  E14.5.merge_EYFP@assays$predictions@data[i,] %>%
                                    max()))
  ggsave(paste("spatial_celltype_estimation_slice1_", celltypelist[i] %>% as.character(), ".pdf", sep=""),
         height = 6, width = 8)
  SpatialFeaturePlot(E14.5.merge_EYFP, 
                     features = celltypelist[i] %>% as.character(), 
                     pt.size.factor = 3, ncol = 2, crop = TRUE,
                     image = "slice2") + 
    scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")), 
                         limits = c(0, 
                                    E14.5.merge_EYFP@assays$predictions@data[i,] %>%
                                      max()))
  ggsave(paste("spatial_celltype_estimation_slice2_", celltypelist[i] %>% as.character(), ".pdf", sep=""),
         height = 6, width = 8)
  
}







##goi
SpatialFeaturePlot(E14.5.merge_EYFP, 
                   features = "EYFP", 
                   pt.size.factor = 3, ncol = 2, crop = TRUE)

########################################
#保存
setwd("~/Downloads/Visium/analysis/merge_E14.5/201117_scRNA-seq_prediction")
#saveRDS(E14.5.merge_EYFP, "E14.5.merge_EYFP.rds")
########################################
E14.5.merge_EYFP <- readRDS("E14.5.merge_EYFP.rds")


########################################
#VisiumからEYFP > 0
########################################
setwd("~/Downloads/Visium/analysis/merge_E14.5/201117_scRNA-seq_prediction/")
E14.5.merge <- readRDS("E14.5.merge2.rds")
SpatialDimPlot(E14.5.merge, crop=T, pt.size.factor = 4, label = T, label.box = F,label.size = 3)

counts_sct <- E14.5.merge$SCT@counts %>% as.data.frame() %>% t() %>% as.data.frame()
dim(counts_sct)
# [1]   702 15076
subset(counts_sct, EYFP >  0) %>% dim
# [1]   313 15076
eyfp_spot <- subset(counts_sct, EYFP >  0) %>% rownames()

E14.5.merge_EYFP <- E14.5.merge[, eyfp_spot]
dim(E14.5.merge_EYFP)
# [1] 15076   313
SpatialDimPlot(E14.5.merge_EYFP, crop=T, pt.size.factor = 4, label = T, label.box = F,label.size = 3)
setwd("~/Downloads/Visium/analysis/merge_E14.5/201117_scRNA-seq_prediction/goi_EYFP")
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
goi <- c("Hmga2", "Mki67", "Top2a", "Ccnb2")
goi <- c("Hapln1", "Vcan", "Tbx20", "Postn", "Dcn", "Lum", "Bgn", "Fmod",
         "Fbln5", "Myom1", "Lmod1", "Jag1" , "Notch3")
goi <- c("Plp1")
for(i in 1 : length(goi)){
  SpatialFeaturePlot(E14.5.merge_EYFP, features = goi[i], pt.size.factor=4) +
    ggsave(paste(goi[i], ".pdf", sep =""), height = 6, width = 8)
}









################################################
#Identification of Spatially Variable Features
################################################
all.markers <- FindAllMarkers(E14.5.merge, only.pos = TRUE, min.pct = 0.25, 
                              logfc.threshold = 0.25)
setwd("~/Downloads/Visium/analysis/merge_E14.5/201117_scRNA-seq_prediction/markers")
write.table(all.markers, "all.markers.txt", quote = F, sep = "\t", col.names = NA)

setwd("~/Downloads/Visium/analysis/merge_E14.5/201117_scRNA-seq_prediction/markers")
all.markers <- read.table("all.markers.txt", header = T)

all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

all.genes <- rownames(E14.5.merge)
E14.5.merge <- ScaleData(E14.5.merge, features = all.genes)
top <-  all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(E14.5.merge, features = top$gene, label = T,
          size = 5) + NoLegend()
ggsave("top10.png", height = 30, width = 15, limitsize = FALSE)
# 1500x1500

SpatialFeaturePlot(object = E14.5.merge, 
                   features = rownames(all.markers)[1:3], alpha = c(0.1, 1), ncol = 3)

all.markers <- read.table("all.markers.txt", header = T)
subset(all.markers, cluster == 19) %>% View
submarkers <- subset(all.markers, cluster == 19)



setwd("~/Downloads/Visium/analysis/merge_E14.5/scRNAseq_integration")

celltypelist <- Seuratcds2@meta.data %>% distinct(celltype)
celltypelist$celltype %>% as.character() 
celltypelist<- celltypelist$celltype %>% as.character() %>% as.list()

for(i in 1:(celltypelist %>% length())){
  SpatialFeaturePlot(E14.5.merge2, 
                     features = celltypelist[i] %>% as.character(), 
                     pt.size.factor = 3, ncol = 2, crop = TRUE)
  ggsave(paste("spatial_celltype_estimation_", celltypelist[i] %>% as.character(), ".pdf", sep=""),
         height = 6, width = 8)
}







