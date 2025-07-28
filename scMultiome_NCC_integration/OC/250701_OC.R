# 統合データでcommunityを計算する

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


#任意のコミュニティに属する遺伝子を抽出する
#まずseurat_clustersラベルした行列を作る
#ただし、クラスタは0から12までで、11は除く
comm_mix<- matrix(nrow = 28, ncol = 1)
comm_mix[,1] <- seq(from = 0, to = 27)
colnames(comm_mix) <- c("cluster")
comm_mix <- as.data.frame(comm_mix)
# Fluidigmでサンプリングした領域に相当するクラスタだけにフィルタリングする
comm_mix <- comm_mix %>% dplyr::filter(cluster %in% c("21","22","10","2","16","18","23","27","0","7","4","6","14"))

setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/OC")

# exp_mat_nonzero <- exp_mat[rowSums(exp_mat) != 0, ]
nonzero_genes <- rownames(exp_mat)[rowSums(exp_mat) != 0]
# 発現量が0の遺伝子があるとエラーになるので除外する

total_cluster = 109
for(i in 1: total_cluster ){
  community_number = i
  inodes <- community %>% subset(community == (community_number)) 
  inodes <- inodes[,1]
  inodes <- as.character(inodes)
  inodes
  # OCなのでsignbnの情報を使う
  child <- signbn %>% dplyr::filter(Parent %in% inodes)
  inodes <- c(inodes, child$Child) %>% unique()
  # 発現がある遺伝子のみにする
  inodes <- intersect(inodes, nonzero_genes)
  #サブデータフレーム作成と遺伝子ごとにz-scale化
  sub.exp_mat <- exp_mat[inodes, ]
  sub.exp_mat <- t(sub.exp_mat) %>% scale() %>% t() %>% as.data.frame()
  #細胞を合わせる（元から合っているけれども）
  #sub.pd <- stpseu[colnames(exp_mat),]
  #tidyな形へ合わせる
  tidydata <- cbind(meta3, as.data.frame(t(sub.exp_mat))) %>% as.data.frame()
  comm <- tidydata %>% 
    group_by(.data$meta3) %>% summarise_all(list(~ mean(.))) #クラスタごとにZ-score化した値の平均値算出

  comm <- comm[, -1]
  comm.sum <- apply(comm,1, mean) %>% as.data.frame()#行ごとに平均する
  colnames(comm.sum) <- paste("OC", (i), sep = "")
  comm_mix <- cbind(comm_mix, as.data.frame(comm.sum))
}
rownames(comm_mix) <- comm_mix$cluster
write.table(comm_mix, "OC_mean_by_cluster.txt", sep ="\t",row.names = F, quote = F)
#clusterごとのcommunity平均値算出完了

#生のheatmap
# comm_mix <- read.table("community_mean_by_cluster.txt", header =T, row.names = 1)
# comm_mix <- t(comm_mix) %>% as.data.frame()
# pheatmap(comm_mix, fontsize_row = 4, angle_col = 45)

#Z-score化したheat map
#どのcommunityが高いかscaleをかける（クラスタ内でどのコミュニティが高いかscale）
#comm_mixは列がクラスタ、行がコミュニティ
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/OC")
comm_mix <- read.table("OC_mean_by_cluster.txt", header =T, row.names = 1)
rownames(comm_mix) <- paste("C", rownames(comm_mix), sep = "")
comm_mix %>% t() %>% scale()  %>% as.data.frame() -> comm_mix2
scale_max = 3
scale_min = -3
comm_mix2[comm_mix2 > scale_max] = scale_max
comm_mix2[comm_mix2 < scale_min] = scale_min
pheatmap(comm_mix2, fontsize_row =8, angle_col = 0)
# 14x6で保存
write.table(comm_mix2, "comm_mix2.txt", sep = "\t", quote = F, col.names = NA)

# 任意の部分だけのOC heatmap
comm_mix2 <- comm_mix2[c("OC12","OC18","OC24","OC26","OC45","OC25","OC86","OC91","OC8","OC23","OC34",
                         "OC4","OC33","OC14","OC58","OC60","OC66","OC44","OC71"),]
pheatmap(comm_mix2, fontsize_row =8, angle_col = 0)
# 4z6で保存


# 次に細胞ごとのOCの値を算出する
#細胞を合わせる（元から合っているけれども）
dim(NCC_merge)
# [1] 33282 32982
comm_cell_mix <- matrix(nrow = 32982, ncol =1) %>% as.data.frame()
comm_cell_mix$V1 <- NCC_merge@meta.data$integrated
colnames(comm_cell_mix) <- c("cluster")

total_cluster = 109
#これはlinkcomm_at$nodeclusters %>% tailで事前に調べておく
for(i in 1: total_cluster ){
  community_number = i
  inodes <- community %>% subset(community == (community_number)) 
  inodes <- inodes[,1]
  inodes <- as.character(inodes)
  inodes
  # OCなのでsignbnの情報を使う
  child <- signbn %>% dplyr::filter(Parent %in% inodes)
  inodes <- c(inodes, child$Child) %>% unique()
  # 発現がある遺伝子のみにする
  inodes <- intersect(inodes, nonzero_genes)
  #サブデータフレーム作成と遺伝子ごとにz-scale化
  sub.exp_mat <- exp_mat[inodes, ]
  sub.exp_mat <- t(sub.exp_mat) %>% scale() %>% as.data.frame()
  
  #tidyな形へ合わせる
  #comm_cell <- sub.tpm[,-c(1,2)]
  comm_cell <- apply(sub.exp_mat, 1, mean) %>% as.data.frame()
  colnames(comm_cell) <- paste("OC", (i), sep = "")
  comm_cell_mix <- cbind(comm_cell_mix, as.data.frame(comm_cell))
}
#cellごとのcommunity平均値算出完了
write.table(comm_cell_mix, "comm_cell_mix.txt", sep ="\t",col.names = NA, quote = F)


#heatmap, 縦に（cellごとにCommunityを）z-score化
comm_cell_mix <- read.table("comm_cell_mix.txt", header = T)
sub_hm <- comm_cell_mix[, -c(1)] %>% t() %>%scale()  %>% as.data.frame()
# scale_max = 3
# scale_min = -3
# sub_hm[sub_hm > scale_max] = scale_max
# sub_hm[sub_hm < scale_min] = scale_min
# pheatmap(sub_hm, fontsize_row = 4, angle_col = 45, show_colnames = F,
#          annotation_col = anocol)
write.table(sub_hm, "comm_cell_mix_scaled.txt", sep = "\t", quote = F, col.names = NA)


# #heatmap, 横に（Communityごとにcellを）z-score化
# anocol <- comm_cell_mix[, c(1)]  %>% as.data.frame()
# rownames(anocol) <- rownames(comm_cell_mix)
# colnames(anocol) <- c("cluster")
# anocol$cluster <- as.factor(anocol$cluster)
# 
# sub_hm <- comm_cell_mix[, -c(1)] %>%  scale() %>% t() %>% as.data.frame()
# scale_max = 3
# scale_min = -3
# sub_hm[sub_hm > scale_max] = scale_max
# sub_hm[sub_hm < scale_min] = scale_min
# pheatmap(sub_hm, fontsize_row = 4, angle_col = 45, show_colnames = F,
#          annotation_col = anocol)
# 
# 
# 



# OCを読み込む
#これは、細胞内でどのコミュニティをscale化したもの
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

# assayを追加
NCC_merge[["OC"]] <- OC

DefaultAssay(NCC_merge) <- "OC"
NCC_merge@meta.data$integrated
DimPlot(NCC_merge, reduction = "umap.rpca")

# OCのAUCを算出する
# markers_OC <- presto:::wilcoxauc.Seurat(X = NCC_merge, group_by = 'seurat_clusters',
#                                         assay = 'counts', seurat_assay = 'OC')
# write.table(markers_OC, "markers_OC_intC.txt", quote = F, sep ="\t", row.names = F)
# 
# markers_OC <- presto:::wilcoxauc.Seurat(X = NCC_merge, group_by = 'sub.cluster',
#                                         assay = 'counts', seurat_assay = 'OC')
# write.table(markers_OC, "markers_OC_sub.cluster.txt", quote = F, sep ="\t", row.names = F)


# 109 OCを可視化
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/OC/goi")
goi <- rownames(NCC_merge)
for(i in 1:length(goi)){
  g1 <- FeaturePlot(NCC_merge, 
                    features = goi[i], reduction = "umap.rpca", slot = "counts", min.cutoff = -1,  max.cutoff = 1)　+ # OCはZ-scored済みなのでcounts slot
    scale_color_gradient2(low = 'grey90', mid = 'grey90', high = "darkgreen", midpoint = 0)
  g1
  ggsave(file = paste(goi[i], ".pdf", sep = ""), plot = g1, width = 4, height = 4)  
}


# マーカーのクラスタごとの色
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
# 6x12で保存
ggsave("Featureplot_OC.png", height = 6, width = 16)








