
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
  dplyr::filter(orig.ident %in% c("E17.5_EYFP", "E17.5"),
                integrated %in% c("21","22","10","2","16","18","23","27","0","7","4","6","14")) %>% rownames()
NCC_merge <- NCC_merge[,coi]
table(NCC_merge@meta.data$orig.ident)
DimPlot(NCC_merge, reduction = "umap.rpca", label = T) 

########################################
#EYFP > 0る
########################################
setwd("/Volumes/Pegasus32R8/workspace/omics_data/Visium/Visium1_embryonic_heart/analysis/merge_E17.5/201117_scRNA-seq_prediciton")
setwd("/Users/iwaseakiyasu/Downloads/Visium_RCTD/250428_RCTD_region/")

E17.5.merge <- readRDS("E17.5.merge2.rds")
SpatialDimPlot(E17.5.merge, crop=T, pt.size.factor = 3, label = T, label.box = F,label.size = 3)

coi <- E17.5.merge@meta.data %>% dplyr::filter(orig.ident == "E17.5_1")
E17.5.merge <- E17.5.merge[, rownames(coi)]
SpatialDimPlot(E17.5.merge, crop=T, pt.size.factor = 3, label = T, label.box = F,label.size = 3, images = "slice3")

# interactive
SpatialDimPlot(E17.5.merge, interactive = TRUE,  images = "slice3")

coi <- E17.5.merge@meta.data %>% dplyr::filter(!seurat_clusters %in% c("10", "27", "25"))
E17.5.merge <- E17.5.merge[, rownames(coi)]


SpatialDimPlot(E17.5.merge, crop=T, pt.size.factor = 3, label = T, label.box = F,label.size = 3,  images = "slice3")

setwd("/Users/iwaseakiyasu/Downloads/Visium_RCTD/250428_RCTD_region/E17.5_1")
# saveRDS(E17.5.merge, "E17.5_1_cut.rds")

counts_sct <- E17.5.merge$SCT@counts %>% as.data.frame() %>% t() %>% as.data.frame()
dim(counts_sct)
# [1]   762 15729
subset(counts_sct, EYFP >  0) %>% dim
# [1]    298 15729
eyfp_spot <- subset(counts_sct, EYFP >  0) %>% rownames()

E17.5.merge_EYFP <- E17.5.merge[, eyfp_spot]
dim(E17.5.merge_EYFP)
# [1] 15729   298
SpatialDimPlot(E17.5.merge_EYFP, crop=T, pt.size.factor = 3, label = T, label.box = F,label.size = 3, images = "slice3")

E17.5.merge_EYFP <- SCTransform(E17.5.merge_EYFP,
                                assay = "Spatial") %>%
  RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30)
DimPlot(E17.5.merge_EYFP, reduction = "umap", group.by = c("ident", "orig.ident"))

# saveRDS(E17.5.merge_EYFP, "E17.5_1_EYFP_cut.rds")
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/Visium_RCTD/250428_RCTD_region/E17.5_1")
E17.5.merge_EYFP <- readRDS("E17.5_1_EYFP_cut.rds")


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
coi <- E17.5.merge_EYFP@meta.data %>% dplyr::filter(orig.ident == "E17.5_1")
E17.5.merge_EYFP <- E17.5.merge_EYFP[, rownames(coi)]
counts <- E17.5.merge_EYFP[["Spatial"]]$counts
coords <- GetTissueCoordinates(E17.5.merge_EYFP, image = "slice3")
colnames(coords) <- c("x", "y")
coords[is.na(colnames(coords))] <- NULL
query <- SpatialRNA(coords, counts, colSums(counts))

RCTD <- create.RCTD(query, reference, max_cores = 8, CELL_MIN_INSTANCE = 3)
RCTD <- run.RCTD(RCTD, doublet_mode = "full")
RCTD_res <- RCTD@results$weights %>% as.data.frame()
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/Visium_RCTD/250428_RCTD_region/E17.5_1")
saveRDS(RCTD, "RCTD_E17.5_1.rds")

E17.5.merge_EYFP <- AddMetaData(E17.5.merge_EYFP, metadata = RCTD_res)

saveRDS(E17.5.merge_EYFP, "E17.5_1_EYFP_RCTD.rds")


setwd("/Users/iwaseakiyasu/Downloads/Visium_RCTD/250428_RCTD_region/E17.5_1")
E17.5_1_EYFP <- readRDS("E17.5_1_EYFP_RCTD.rds")

# aspect ratio
coord <- GetTissueCoordinates(object = E17.5_1_EYFP@images$slice4)
colnames(coord) <- c("x", "y")
# calculate the aspect ratio of rows to columns
ratio_E17.5_1<- (max(coord$x) - min(coord$x)) / (max(coord$y) - min(coord$y))
# force the image into the right aspect ratio
SpatialDimPlot(E17.5_1_EYFP, crop = TRUE, images = "slice4") + theme(aspect.ratio = ratio_E17.5_1)

SpatialFeaturePlot(E17.5_1_EYFP, features = "6",  pt.size.factor = 3, ncol = 2, crop = TRUE, images = "slice4") +
  theme(aspect.ratio = ratio_E17.5_1)






