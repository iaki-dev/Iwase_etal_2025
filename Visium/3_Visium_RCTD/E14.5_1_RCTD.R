
library(hdf5r)
library(Seurat)
#library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
library(SeuratObject)


setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered")
NCC_merge <- readRDS("NCC_merge_RPCA_velo_integrated.rds")
Idents(NCC_merge) <- NCC_merge@meta.data$integrated
DimPlot(NCC_merge, reduction = "umap.rpca", group.by = "rpca_clusters", label = T)


coi <- NCC_merge@meta.data %>% 
  dplyr::filter(orig.ident %in% c("E14.5_EYFP", "E14.5"),
                integrated %in% c("21","22","10","2","16","18","23","27","0","7","4","6","14")) %>% rownames()
NCC_merge <- NCC_merge[,coi]
table(NCC_merge@meta.data$orig.ident)
DimPlot(NCC_merge, reduction = "umap.rpca", label = T) 

########################################
#EYFP > 0
########################################
setwd("/Volumes/Pegasus32R8/workspace/omics_data/Visium/Visium1_embryonic_heart/analysis/merge_E14.5/201117_scRNA-seq_prediciton")
setwd("/Users/iwaseakiyasu/Downloads/Visium_RCTD/250428_RCTD_region/")
E14.5.merge <- readRDS("E14.5.merge2.rds")
SpatialDimPlot(E14.5.merge, crop=T, pt.size.factor = 3, label = T, label.box = F,label.size = 3)

coi <- E14.5.merge@meta.data %>% dplyr::filter(orig.ident == "E14.5_1")
E14.5.merge <- E14.5.merge[, rownames(coi)]
SpatialDimPlot(E14.5.merge, crop=T, pt.size.factor = 3, label = T, label.box = F,label.size = 3, images = "slice1")

# interactive
SpatialDimPlot(E14.5.merge, interactive = TRUE,  images = "slice1")

coi <- E14.5.merge@meta.data %>% dplyr::filter(!seurat_clusters %in% c("8", "19"))
E14.5.merge <- E14.5.merge[, rownames(coi)]

SpatialDimPlot(E14.5.merge, crop=T, pt.size.factor = 3, label = T, label.box = F,label.size = 3,  images = "slice1")


setwd("/Users/iwaseakiyasu/Downloads/Visium_RCTD/250428_RCTD_region/E14.5_1")


counts_sct <- E14.5.merge$SCT@counts %>% as.data.frame() %>% t() %>% as.data.frame()
dim(counts_sct)
# [1]   384 15076
subset(counts_sct, EYFP >  0) %>% dim
# [1]    142 15729
eyfp_spot <- subset(counts_sct, EYFP >  0) %>% rownames()

E14.5.merge_EYFP <- E14.5.merge[, eyfp_spot]
dim(E14.5.merge_EYFP)
# [1] 15729   142
SpatialDimPlot(E14.5.merge_EYFP, crop=T, pt.size.factor = 3, label = T, label.box = F,label.size = 3, images = "slice1")

E14.5.merge_EYFP <- SCTransform(E14.5.merge_EYFP,
                                assay = "Spatial") %>%
  RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30)
DimPlot(E14.5.merge_EYFP, reduction = "umap", group.by = c("ident", "orig.ident"))

# saveRDS(E14.5.merge_EYFP, "E14.5_1_EYFP_cut.rds")
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/Visium_RCTD/250428_RCTD_region/E14.5_1")
E14.5.merge_EYFP <- readRDS("E14.5_1_EYFP_cut.rds")


NCC_merge2 <- NCC_merge 
NCC_merge2@meta.data$rpca_clusters %>% table


library(spacexr)
ref <- NCC_merge2
ref <- UpdateSeuratObject(ref)
ref@meta.data$integrated <- as.character(ref@meta.data$integrated) %>% as.numeric()
ref@meta.data$integrated <- as.factor(ref@meta.data$integrated)
Idents(ref) <- "integrated"
# extract information to pass to the RCTD Reference function
counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$integrated)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster, nUMI)

# set up query with the RCTD function SpatialRNA
coi <- E14.5.merge_EYFP@meta.data %>% dplyr::filter(orig.ident == "E14.5_1")
E14.5.merge_EYFP <- E14.5.merge_EYFP[, rownames(coi)]
counts <- E14.5.merge_EYFP[["Spatial"]]$counts
coords <- GetTissueCoordinates(E14.5.merge_EYFP, image = "slice1")
colnames(coords) <- c("x", "y")
coords[is.na(colnames(coords))] <- NULL
query <- SpatialRNA(coords, counts, colSums(counts))

RCTD <- create.RCTD(query, reference, max_cores = 8, CELL_MIN_INSTANCE = 3)
RCTD <- run.RCTD(RCTD, doublet_mode = "full")
RCTD_res <- RCTD@results$weights %>% as.data.frame()
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/Visium_RCTD/250428_RCTD_region/E14.5_1")
saveRDS(RCTD, "RCTD_E14.5_1.rds")

E14.5.merge_EYFP <- AddMetaData(E14.5.merge_EYFP, metadata = RCTD_res)

saveRDS(E14.5.merge_EYFP, "E14.5_1_EYFP_RCTD.rds")


setwd("/Users/iwaseakiyasu/Downloads/Visium_RCTD/250428_RCTD_region/E14.5_1")
E14.5_1_EYFP <- readRDS("E14.5_1_EYFP_RCTD.rds")

# aspect ratio
coord <- GetTissueCoordinates(object = E14.5_1_EYFP@images$slice1)
colnames(coord) <- c("x", "y")
# calculate the aspect ratio of rows to columns
ratio_E14.5_1<- (max(coord$x) - min(coord$x)) / (max(coord$y) - min(coord$y))
# force the image into the right aspect ratio
SpatialDimPlot(E14.5_1_EYFP, crop = TRUE, images = "slice1") + theme(aspect.ratio = ratio_E14.5_1)

SpatialFeaturePlot(E14.5_1_EYFP, features = "18",  pt.size.factor = 3, ncol = 2, crop = TRUE, images = "slice1") +
  theme(aspect.ratio = ratio_E14.5_1)






