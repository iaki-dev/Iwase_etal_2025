

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
library(pheatmap)

setwd("/Volumes/Pegasus32R8/workspace/omics_data/scRNA-seq/20200403_Seuratv3_based2/SiGN-BN/200514_test/4thv2/modularity_overlapward2")
community <- read.table("community.txt", header = T)
community$node <- community$node %>% str_replace(pattern = ":", replacement = "")
colnames(community)[2] <- "community"
write.table(community, "community.txt", sep = "\t", quote = F,row.names = F)

# SiGN-BN
signbn <- read.table("../result.txt", header = T)

setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered")
NCC_merge <- readRDS("NCC_merge_RPCA_velo_integrated.rds")
Idents(NCC_merge) <- NCC_merge@meta.data$integrated
DimPlot(NCC_merge, reduction = "umap.rpca", group.by = "integrated", label = T)
coi <- NCC_merge@meta.data %>% 
  dplyr::filter(integrated %in% c("21","22","10","2","16","18","23","27","0","7","4","6","14")) %>% rownames()
NCC_merge <- NCC_merge[,coi]
table(NCC_merge@meta.data$orig.ident)
DimPlot(NCC_merge, reduction = "umap.rpca", label = T) 

exp_mat <- NCC_merge@assays$RNA@layers$data %>% as.data.frame()
dim(exp_mat)
# [1] 33282 32982
rownames(exp_mat) <- rownames(NCC_merge)
colnames(exp_mat) <- colnames(NCC_merge)

meta3 <- NCC_merge@meta.data$integrated %>% as.character() %>% as.numeric()



comm_mix<- matrix(nrow = 28, ncol = 1)
comm_mix[,1] <- seq(from = 0, to = 27)
colnames(comm_mix) <- c("cluster")
comm_mix <- as.data.frame(comm_mix)
comm_mix <- comm_mix %>% dplyr::filter(cluster %in% c("21","22","10","2","16","18","23","27","0","7","4","6","14"))

setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/OC")

nonzero_genes <- rownames(exp_mat)[rowSums(exp_mat) != 0]

total_cluster = 109
for(i in 1: total_cluster ){
  community_number = i
  inodes <- community %>% subset(community == (community_number)) 
  inodes <- inodes[,1]
  inodes <- as.character(inodes)
  inodes
  
  child <- signbn %>% dplyr::filter(Parent %in% inodes)
  inodes <- c(inodes, child$Child) %>% unique()
 
  inodes <- intersect(inodes, nonzero_genes)
  
  sub.exp_mat <- exp_mat[inodes, ]
  sub.exp_mat <- t(sub.exp_mat) %>% scale() %>% t() %>% as.data.frame()
 
  tidydata <- cbind(meta3, as.data.frame(t(sub.exp_mat))) %>% as.data.frame()
  comm <- tidydata %>% 
    group_by(.data$meta3) %>% summarise_all(list(~ mean(.))) 

  comm <- comm[, -1]
  comm.sum <- apply(comm,1, mean) %>% as.data.frame()
  colnames(comm.sum) <- paste("OC", (i), sep = "")
  comm_mix <- cbind(comm_mix, as.data.frame(comm.sum))
}
rownames(comm_mix) <- comm_mix$cluster
write.table(comm_mix, "OC_mean_by_cluster.txt", sep ="\t",row.names = F, quote = F)


# heatmap
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/OC")
comm_mix <- read.table("OC_mean_by_cluster.txt", header =T, row.names = 1)
rownames(comm_mix) <- paste("C", rownames(comm_mix), sep = "")
comm_mix %>% t() %>% scale()  %>% as.data.frame() -> comm_mix2
scale_max = 3
scale_min = -3
comm_mix2[comm_mix2 > scale_max] = scale_max
comm_mix2[comm_mix2 < scale_min] = scale_min
pheatmap(comm_mix2, fontsize_row =8, angle_col = 0)
write.table(comm_mix2, "comm_mix2.txt", sep = "\t", quote = F, col.names = NA)


comm_mix2 <- comm_mix2[c("OC12","OC18","OC24","OC26","OC45","OC25","OC86","OC91","OC8","OC23","OC34",
                         "OC4","OC33","OC14","OC58","OC60","OC66","OC44","OC71"),]
pheatmap(comm_mix2, fontsize_row =8, angle_col = 0)


dim(NCC_merge)
# [1] 33282 32982
comm_cell_mix <- matrix(nrow = 32982, ncol =1) %>% as.data.frame()
comm_cell_mix$V1 <- NCC_merge@meta.data$integrated
colnames(comm_cell_mix) <- c("cluster")

total_cluster = 109
for(i in 1: total_cluster ){
  community_number = i
  inodes <- community %>% subset(community == (community_number)) 
  inodes <- inodes[,1]
  inodes <- as.character(inodes)
  inodes

  child <- signbn %>% dplyr::filter(Parent %in% inodes)
  inodes <- c(inodes, child$Child) %>% unique()

  inodes <- intersect(inodes, nonzero_genes)

  sub.exp_mat <- exp_mat[inodes, ]
  sub.exp_mat <- t(sub.exp_mat) %>% scale() %>% as.data.frame()
  

  comm_cell <- apply(sub.exp_mat, 1, mean) %>% as.data.frame()
  colnames(comm_cell) <- paste("OC", (i), sep = "")
  comm_cell_mix <- cbind(comm_cell_mix, as.data.frame(comm_cell))
}
write.table(comm_cell_mix, "comm_cell_mix.txt", sep ="\t",col.names = NA, quote = F)



comm_cell_mix <- read.table("comm_cell_mix.txt", header = T)
sub_hm <- comm_cell_mix[, -c(1)] %>% t() %>%scale()  %>% as.data.frame()
write.table(sub_hm, "comm_cell_mix_scaled.txt", sep = "\t", quote = F, col.names = NA)


setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/OC")
gf <- read_tsv("comm_cell_mix_scaled.txt")
dim(gf)
#[1] 109 32983
tmp.names <- paste("OC", 1:nrow(gf), sep = "")
gf <- as.data.frame(gf)
rownames(gf) <- tmp.names
gf <- gf[-c(1)]
gf <- gf %>% as.data.frame()

OC <- CreateAssayObject(counts = gf)

# Add Assay 
NCC_merge[["OC"]] <- OC

DefaultAssay(NCC_merge) <- "OC"
NCC_merge@meta.data$integrated
DimPlot(NCC_merge, reduction = "umap.rpca")



setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/OC/goi")
goi <- rownames(NCC_merge)
for(i in 1:length(goi)){
  g1 <- FeaturePlot(NCC_merge, 
                    features = goi[i], reduction = "umap.rpca", slot = "counts", min.cutoff = -1,  max.cutoff = 1)　+ # OCはZ-scored済みなのでcounts slot
    scale_color_gradient2(low = 'grey90', mid = 'grey90', high = "darkgreen", midpoint = 0)
  g1
  ggsave(file = paste(goi[i], ".pdf", sep = ""), plot = g1, width = 4, height = 4)  
}



library(scales)
library(patchwork)
show_col(hue_pal()(13))
usingcols <- hue_pal()(13)

# Pharygeal and Transitional
g1 <- FeaturePlot(NCC_merge, features = "OC100", reduction = "umap.rpca", slot = "counts", min.cutoff = -1,  max.cutoff = 1,
                  cols = c("grey90", usingcols[1])) +
  scale_color_gradient2(low = 'grey90', mid = 'grey90', high = usingcols[1],midpoint = 0)
g2 <- FeaturePlot(NCC_merge, features = "OC71", reduction = "umap.rpca", slot = "counts", min.cutoff = -1,  max.cutoff = 1,
                  cols = c("grey90", usingcols[1])) +
  scale_color_gradient2(low = 'grey90', mid = 'grey90', high = usingcols[1],midpoint = 0)
g3 <- FeaturePlot(NCC_merge, features = "OC45", reduction = "umap.rpca", slot = "counts", min.cutoff = -1,  max.cutoff = 1,
                  cols = c("grey90", usingcols[1])) +
  scale_color_gradient2(low = 'grey90', mid = 'grey90', high = usingcols[1],midpoint = 0)
g4 <- FeaturePlot(NCC_merge, features = "OC18", reduction = "umap.rpca", slot = "counts", min.cutoff = -1,  max.cutoff = 1,
                  cols = c("grey90", usingcols[1])) +
  scale_color_gradient2(low = 'grey90', mid = 'grey90', high = usingcols[1],midpoint = 0)

# cushion
g5 <- FeaturePlot(NCC_merge, features = "OC86", reduction = "umap.rpca", slot = "counts", min.cutoff = -1,  max.cutoff = 1,
                  cols = c("grey90", usingcols[3])) +
  scale_color_gradient2(low = 'grey90', mid = 'grey90', high = usingcols[3],midpoint = 0)
g6 <- FeaturePlot(NCC_merge, features = "OC58", reduction = "umap.rpca", slot = "counts", min.cutoff = -1, max.cutoff = 1,
                  cols = c("grey90", usingcols[3])) +
  scale_color_gradient2(low = 'grey90', mid = 'grey90', high = usingcols[3],midpoint = 0)
g7 <- FeaturePlot(NCC_merge, features = "OC25", reduction = "umap.rpca", slot = "counts", min.cutoff = -1, max.cutoff = 1,
                  cols = c("grey90", usingcols[3])) +
  scale_color_gradient2(low = 'grey90', mid = 'grey90', high = usingcols[3], midpoint = 0)

# AP septum
g8 <- FeaturePlot(NCC_merge, features = "OC60", reduction = "umap.rpca", slot = "counts", min.cutoff = -1, max.cutoff = 1,
                  cols = c("grey90", usingcols[5])) +
  scale_color_gradient2(low = 'grey90', mid = 'grey90', high = usingcols[5], midpoint = 0)

# subvalvular
g9 <- FeaturePlot(NCC_merge, features = "OC23", reduction = "umap.rpca", slot = "counts", min.cutoff = -1, max.cutoff = 1,
                   cols = c("grey90", usingcols[4])) +
  scale_color_gradient2(low = 'grey90', mid = 'grey90', high = usingcols[4], midpoint = 0)

# SMC
g10 <- FeaturePlot(NCC_merge, features = "OC33", reduction = "umap.rpca", slot = "counts", min.cutoff = -1,  max.cutoff = 1,
                  cols = c("grey90", usingcols[7])) +
  scale_color_gradient2(low = 'grey90', mid = 'grey90', high = usingcols[7],midpoint = 0)
g11 <- FeaturePlot(NCC_merge, features = "OC4", reduction = "umap.rpca", slot = "counts",  min.cutoff = -1, max.cutoff = 1,
                  cols = c("grey90", usingcols[7])) +
  scale_color_gradient2(low = 'grey90', mid = 'grey90', high = usingcols[7],midpoint = 0)

# CA
g12 <- FeaturePlot(NCC_merge, features = "OC8", reduction = "umap.rpca", slot = "counts", min.cutoff = -1, max.cutoff = 1,
                  cols = c("grey90", usingcols[8])) +
  scale_color_gradient2(low = 'grey90', mid = 'grey90', high = usingcols[8],midpoint = 0)

(g1 | g2 | g3 | g5 |g5 |g6) / ( g7 | g8 |g9| g10 | g11 | g12 )

ggsave("Featureplot_OC.png", height = 6, width = 16)








