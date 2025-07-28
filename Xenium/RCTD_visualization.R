
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


setwd("/Volumes/Pegasus32R8/workspace/omics_data/Xenium/analysis/Wnt1_Cre_E115_1/RCTD/NCC")
xenium_E11.5_NCC <- readRDS("xenium.obj_NCC.rds")
df_E11.5_NCC <- xenium_E11.5_NCC@images$fov$centroids@coords %>% as.data.frame()
df_E11.5_NCC <- cbind(df_E11.5_NCC, xenium_E11.5_NCC@meta.data)
df_E11.5_NCC$orig.ident <- "E11.5_NCC"

setwd("/Volumes/Pegasus32R8/workspace/omics_data/Xenium/analysis/Wnt1_Cre_E115_1/RCTD/nonNCC")
xenium_E11.5_nonNCC <- readRDS("xenium.obj_nonNCC.rds")
df_E11.5_nonNCC <- xenium_E11.5_nonNCC@images$fov$centroids@coords %>% as.data.frame()
df_E11.5_nonNCC <- cbind(df_E11.5_nonNCC, xenium_E11.5_nonNCC@meta.data)
df_E11.5_nonNCC$orig.ident <- "E11.5_nonNCC"

setwd("/Volumes/Pegasus32R8/workspace/omics_data/Xenium/analysis/Wnt1_Cre_E125_8/RCTD/NCC")
xenium_E12.5_NCC <- readRDS("xenium.obj_NCC.rds")
df_E12.5_NCC <- xenium_E12.5_NCC@images$fov$centroids@coords %>% as.data.frame()
df_E12.5_NCC <- cbind(df_E12.5_NCC, xenium_E12.5_NCC@meta.data)
df_E12.5_NCC$orig.ident <- "E12.5_NCC"

setwd("/Volumes/Pegasus32R8/workspace/omics_data/Xenium/analysis/Wnt1_Cre_E125_8/RCTD/nonNCC")
xenium_E12.5_nonNCC <- readRDS("xenium.obj_nonNCC.rds")
df_E12.5_nonNCC <- xenium_E12.5_nonNCC@images$fov$centroids@coords %>% as.data.frame()
df_E12.5_nonNCC <- cbind(df_E12.5_nonNCC, xenium_E12.5_nonNCC@meta.data)
df_E12.5_nonNCC$orig.ident <- "E12.5_nonNCC"

df <- rbind(df_E11.5_NCC[c("x","y","orig.ident","predicted.celltype")], 
            df_E12.5_NCC[c("x","y","orig.ident","predicted.celltype")], 
            df_E11.5_nonNCC[c("x","y","orig.ident","predicted.celltype")], 
            df_E12.5_nonNCC[c("x","y","orig.ident","predicted.celltype")])

ggplot(df, aes(x=x, y=y, color = predicted.celltype)) +
  geom_point() +
  facet_wrap(~orig.ident)



df_transformed <- df %>%
  group_by(orig.ident) %>%
  mutate(
    x_center = mean(x),
    y_center = mean(y),
    x = x + (2500 - x_center),
    y = y + (3000 - y_center)
  )

ggplot(df_transformed, aes(x=x, y=y, color = predicted.celltype)) +
  geom_point() +
  facet_wrap(~orig.ident)


setwd("/Volumes/Pegasus32R8/workspace/omics_data/Xenium/analysis/RCTD_visualization/")
write.table(df, "df.txt", quote = F, sep = "\t", col.names = NA)
write.table(df_transformed, "df_transformed.txt", quote = F, sep = "\t", col.names = NA)


df_transformed$predicted.celltype <- as.character(df_transformed$predicted.celltype) %>% as.numeric()
df_transformed$predicted.celltype <- as.factor(df_transformed$predicted.celltype)
df_transformed$predicted.celltype <- factor(df_transformed$predicted.celltype ,
                                            levels = 0:20)
df_transformed %>% dplyr::filter(orig.ident == "E11.5_NCC") %>% 
ggplot(aes(x=x, y=y, color = predicted.celltype)) +
  geom_point(size = 0.01) + NoLegend() +
  xlim(1500, 4000) +
  ylim(1500, 4000) +
  theme(
    panel.background = element_rect(fill = "black"),  
    plot.background = element_rect(fill = "black"),   
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    axis.text = element_blank(),
    axis.title =element_blank()
  ) +
  scale_y_reverse() 
ggsave("decomp_E11.5_NCC_adjust.pdf", width = 7, height = 5)

df_transformed %>% dplyr::filter(orig.ident == "E12.5_NCC") %>% 
  ggplot(aes(x=x, y=y, color = predicted.celltype)) +
  geom_point(size = 0.01) + NoLegend() +
  xlim(1500, 4000) +
  ylim(1500, 4000) +
  theme(
    panel.background = element_rect(fill = "black"),  
    plot.background = element_rect(fill = "black"),   
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    axis.text = element_blank(),
    axis.title =element_blank()
  ) +
  scale_x_reverse() +
  coord_flip()
ggsave("decomp_E12.5_NCC_adjust.pdf", width = 7, height = 5)

df_transformed %>% dplyr::filter(orig.ident == "E11.5_nonNCC") %>% 
  ggplot(aes(x=x, y=y, color = predicted.celltype)) +
  geom_point(size = 0.01) + NoLegend() +
  xlim(1500, 4000) +
  ylim(1500, 4000) +
  theme(
    panel.background = element_rect(fill = "black"),  
    plot.background = element_rect(fill = "black"), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.text = element_blank(),
    axis.title =element_blank()
  ) +
  scale_y_reverse() 
ggsave("decomp_E11.5_nonNCC_adjust.pdf", width = 7, height = 5)

df_transformed %>% dplyr::filter(orig.ident == "E12.5_nonNCC") %>% 
  ggplot(aes(x=x, y=y, color = predicted.celltype)) +
  geom_point(size = 0.01) + NoLegend() +
  xlim(1500, 4000) +
  ylim(1500, 4000) +
  theme(
    panel.background = element_rect(fill = "black"),  
    plot.background = element_rect(fill = "black"),   # 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.text = element_blank(),
    axis.title =element_blank()
  ) +
  scale_x_reverse() +
  coord_flip()
ggsave("decomp_E12.5_nonNCC_adjust.pdf", width = 7, height = 5)


df_transformed %>% dplyr::filter(orig.ident == "E12.5_nonNCC") %>% 
  ggplot(aes(x=x, y=y, color = predicted.celltype)) +
  geom_point(size = 0.01) + 
  xlim(1500, 4000) +
  ylim(1500, 4000) +
  theme(
    panel.background = element_rect(fill = "black"),  
    plot.background = element_rect(fill = "black"),   
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.background = element_rect(fill = "black"),  
    legend.key = element_rect(fill = "black"),
    legend.text = element_text(color = "white"),  
    legend.title = element_text(color = "white")  
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))
ggsave("decomp_E12.5_nonNCC_adjust_legend.pdf", width = 7, height = 5)





