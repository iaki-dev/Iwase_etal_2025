

library(hdf5r)
library(Seurat)
#library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
library(SeuratObject)
library(tidyverse)



setwd("/Users/iwaseakiyasu/Downloads/Visium_RCTD/250428_RCTD_region/E17.5_1")
E17.5_1_EYFP <- readRDS("E17.5_1_EYFP_RCTD.rds")
# aspect ratio
coord <- GetTissueCoordinates(object = E17.5_1_EYFP@images$slice3)
colnames(coord) <- c("x", "y")
# calculate the aspect ratio of rows to columns
ratio_E17.5_1 <- (max(coord$x) - min(coord$x)) / (max(coord$y) - min(coord$y))
# force the image into the right aspect ratio
SpatialDimPlot(E17.5_1_EYFP, crop = TRUE, images = "slice3") + theme(aspect.ratio = ratio_E17.5_1)

SpatialFeaturePlot(E17.5_1_EYFP, features = "21", images = "slice3") + theme(aspect.ratio = ratio_E17.5_1)

setwd("/Users/iwaseakiyasu/Downloads/Visium_RCTD/250428_RCTD_region/E17.5_2")
E17.5_2_EYFP <- readRDS("E17.5_2_EYFP_RCTD.rds")
# aspect ratio
coord <- GetTissueCoordinates(object = E17.5_2_EYFP@images$slice4)
colnames(coord) <- c("x", "y")
# calculate the aspect ratio of rows to columns
ratio_E17.5_2<- (max(coord$x) - min(coord$x)) / (max(coord$y) - min(coord$y))
# force the image into the right aspect ratio
SpatialDimPlot(E17.5_2_EYFP, crop = TRUE, images = "slice4") + theme(aspect.ratio = ratio_E17.5_2)

SpatialFeaturePlot(E17.5_2_EYFP, features = "21", images = "slice4") + theme(aspect.ratio = ratio_E17.5_2)


meta1 <- E17.5_1_EYFP@meta.data
meta2 <- E17.5_2_EYFP@meta.data
meta_merge <- rbind(meta1, meta2)
meta_merge <- meta_merge[,9:ncol(meta_merge)]
colnames(meta_merge) <- paste("cluster", colnames(meta_merge), sep = "")


setwd("/Users/iwaseakiyasu/Downloads/Visium_RCTD/250428_RCTD_region")
library(RColorBrewer)

celltypelist <- colnames(meta_merge)
for(i in 1:(celltypelist %>% length())){
  celltypelist[i]
  tmp.name <- str_replace(celltypelist, pattern = "cluster", replacement = "")[i]
  SpatialFeaturePlot(E17.5_1_EYFP, features = tmp.name, pt.size.factor = 3, images = "slice3") + 
    theme(aspect.ratio = ratio_E17.5_1) +
    scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")), 
                         limits = c(0, 
                                    meta_merge[i] %>% 
                                      max()))
  ggsave(paste("spatial_celltype_estimation_slice3_", celltypelist[i] %>% as.character(), ".pdf", sep=""),
         height = 6, width = 8)
  
  SpatialFeaturePlot(E17.5_2_EYFP, features = tmp.name, pt.size.factor = 3, images = "slice4") + 
    theme(aspect.ratio = ratio_E17.5_2) +
    scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")), 
                         limits = c(0, 
                                    meta_merge[i] %>% 
                                      max()))
  ggsave(paste("spatial_celltype_estimation_slice4_", celltypelist[i] %>% as.character(), ".pdf", sep=""),
         height = 6, width = 8)
  
}



