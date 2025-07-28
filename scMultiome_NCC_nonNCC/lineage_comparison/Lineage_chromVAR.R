
library(chromVAR)
library(JASPAR2022)
library(TFBSTools)
library(motifmatchr)
library(biovizBase)
library(GenomeInfoDb)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(Seurat)
library(cowplot)
library(tidyverse)
library(Signac)


# lineage_chromVAR
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/v10/v10.2.1/data")
gexatac_merge <- readRDS("4_gexatac_merge_links.rds")
DefaultAssay(gexatac_merge) <- "ATAC"



markers.volcano <- c()
table_df2 <- gexatac_merge@meta.data %>% dplyr::select(sub.cluster, sample) %>% table() %>% as.data.frame()
table_df2 <- table_df2 %>% dplyr::filter(sample %in% c("E11.5_NCC", "E11.5_nonNCC"))
table_df2 <- table_df2 %>% pivot_wider(names_from = sample, values_from = Freq)
colnames(table_df2) <-c("sub.cluster", "NCC", "non_NCC")
table_df2$threshold_NCC <- table_df2$NCC / sum(table_df2$NCC) * 100
table_df2$threshold_non_NCC <- table_df2$non_NCC / sum(table_df2$non_NCC) * 100

cluster_of_interest <- table_df2 %>% dplyr::filter(threshold_NCC > 0.5, threshold_non_NCC > 0.5) %>% pull(sub.cluster) %>% as.character()
cluster_of_interest
for(i in 1:length(cluster_of_interest)){
  
  cluster_of_interest.tmp <- cluster_of_interest[i]
  coi <- gexatac_merge@meta.data %>% dplyr::filter(sub.cluster %in% cluster_of_interest.tmp)
  coi <- coi %>% dplyr::filter(sample %in% c("E11.5_NCC", "E11.5_nonNCC"))
  tmp <- gexatac_merge[, rownames(coi)]
  
  Idents(tmp) <- "lineage"
  
  DefaultAssay(tmp) <- "chromvar"

  markers.volcano.tmp <- presto:::wilcoxauc.Seurat(X = tmp, group_by = 'lineage', assay = 'data', seurat_assay = 'chromvar')
  
  markers.volcano.tmp <- markers.volcano.tmp %>% mutate(sub.cluster = cluster_of_interest[i]) %>%
    mutate(stage = "E11.5")
  
  markers.volcano <- rbind(markers.volcano, markers.volcano.tmp)
}

table_df2 <- gexatac_merge@meta.data %>% dplyr::select(sub.cluster, sample) %>% table() %>% as.data.frame()
table_df2 <- table_df2 %>% dplyr::filter(sample %in% c("E12.5_NCC", "E12.5_nonNCC"))
table_df2 <- table_df2 %>% pivot_wider(names_from = sample, values_from = Freq)
colnames(table_df2) <-c("sub.cluster", "NCC", "non_NCC")
table_df2$threshold_NCC <- table_df2$NCC / sum(table_df2$NCC) * 100
table_df2$threshold_non_NCC <- table_df2$non_NCC / sum(table_df2$non_NCC) * 100

cluster_of_interest <- table_df2 %>% dplyr::filter(threshold_NCC > 0.5, threshold_non_NCC > 0.5) %>% pull(sub.cluster) %>% as.character()
cluster_of_interest
for(i in 1:length(cluster_of_interest)){
  
  cluster_of_interest.tmp <- cluster_of_interest[i]
  coi <- gexatac_merge@meta.data %>% dplyr::filter(sub.cluster %in% cluster_of_interest.tmp)
  coi <- coi %>% dplyr::filter(sample %in% c("E12.5_NCC", "E12.5_nonNCC"))
  tmp <- gexatac_merge[, rownames(coi)]
  
  Idents(tmp) <- "lineage"
  
  DefaultAssay(tmp) <- "chromvar"
  
  markers.volcano.tmp <- presto:::wilcoxauc.Seurat(X = tmp, group_by = 'lineage', assay = 'data', seurat_assay = 'chromvar')
  
  markers.volcano.tmp <- markers.volcano.tmp %>% mutate(sub.cluster = cluster_of_interest[i]) %>%

    mutate(stage = "E12.5")
  
  markers.volcano <- rbind(markers.volcano, markers.volcano.tmp)
}

table_df2 <- gexatac_merge@meta.data %>% dplyr::select(sub.cluster, sample) %>% table() %>% as.data.frame()
table_df2 <- table_df2 %>% dplyr::filter(sample %in% c("E14.5_NCC", "E14.5_nonNCC"))
table_df2 <- table_df2 %>% pivot_wider(names_from = sample, values_from = Freq)
colnames(table_df2) <-c("sub.cluster", "NCC", "non_NCC")
table_df2$threshold_NCC <- table_df2$NCC / sum(table_df2$NCC) * 100
table_df2$threshold_non_NCC <- table_df2$non_NCC / sum(table_df2$non_NCC) * 100

cluster_of_interest <- table_df2 %>% dplyr::filter(threshold_NCC > 0.5, threshold_non_NCC > 0.5) %>% pull(sub.cluster) %>% as.character()
cluster_of_interest
for(i in 1:length(cluster_of_interest)){
  
  cluster_of_interest.tmp <- cluster_of_interest[i]
  coi <- gexatac_merge@meta.data %>% dplyr::filter(sub.cluster %in% cluster_of_interest.tmp)
  coi <- coi %>% dplyr::filter(sample %in% c("E14.5_NCC", "E14.5_nonNCC"))
  tmp <- gexatac_merge[, rownames(coi)]
  
  Idents(tmp) <- "lineage"
  
  DefaultAssay(tmp) <- "chromvar"
 
  markers.volcano.tmp <- presto:::wilcoxauc.Seurat(X = tmp, group_by = 'lineage', assay = 'data', seurat_assay = 'chromvar')
  
  markers.volcano.tmp <- markers.volcano.tmp %>% mutate(sub.cluster = cluster_of_interest[i]) %>%

    mutate(stage = "E14.5")
  
  markers.volcano <- rbind(markers.volcano, markers.volcano.tmp)
}

table_df2 <- gexatac_merge@meta.data %>% dplyr::select(sub.cluster, sample) %>% table() %>% as.data.frame()
table_df2 <- table_df2 %>% dplyr::filter(sample %in% c("E17.5_NCC", "E17.5_nonNCC"))
table_df2 <- table_df2 %>% pivot_wider(names_from = sample, values_from = Freq)
colnames(table_df2) <-c("sub.cluster", "NCC", "non_NCC")
table_df2$threshold_NCC <- table_df2$NCC / sum(table_df2$NCC) * 100
table_df2$threshold_non_NCC <- table_df2$non_NCC / sum(table_df2$non_NCC) * 100

cluster_of_interest <- table_df2 %>% dplyr::filter(threshold_NCC > 0.5, threshold_non_NCC > 0.5) %>% pull(sub.cluster) %>% as.character()
cluster_of_interest
for(i in 1:length(cluster_of_interest)){
  
  cluster_of_interest.tmp <- cluster_of_interest[i]
  coi <- gexatac_merge@meta.data %>% dplyr::filter(sub.cluster %in% cluster_of_interest.tmp)
  coi <- coi %>% dplyr::filter(sample %in% c("E17.5_NCC", "E17.5_nonNCC"))
  tmp <- gexatac_merge[, rownames(coi)]
  
  Idents(tmp) <- "lineage"
  
  DefaultAssay(tmp) <- "chromvar"
  
  markers.volcano.tmp <- presto:::wilcoxauc.Seurat(X = tmp, group_by = 'lineage', assay = 'data', seurat_assay = 'chromvar')
  
  markers.volcano.tmp <- markers.volcano.tmp %>% mutate(sub.cluster = cluster_of_interest[i]) %>%

    mutate(stage = "E17.5")
  
  markers.volcano <- rbind(markers.volcano, markers.volcano.tmp)
}
# pct_diff
markers.volcano$pct_diff <- markers.volcano$pct_in - markers.volcano$pct_out
# pct_FC
markers.volcano$pct_FC <- markers.volcano$pct_in / markers.volcano$pct_out
markers.volcano$gene <- markers.volcano$feature
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/v10/v10.2.1/lineage/chromVAR")
write.table(markers.volcano, "markers_volcano_NCC_nonNCC_subcluster_stage_chromvar.txt", quote = F, sep = "\t", row.names = F)






# load 
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/v10/v10.2.1/lineage/DAP")
markers.DAP <- read.table("markers_volcano_NCC_nonNCC_DAP_filt.txt", header = T)

setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/v10/v10.2.1/lineage/DEG")
markers_rna <- read.table("markers_volcano_NCC_nonNCC_DEG_filt.txt", header = T)
markers_rna <- subset(markers_rna, !grepl('^Hbb', gene))
markers_rna <- subset(markers_rna, !grepl('^Hba', gene))
write.table(markers_rna, "markers_volcano_NCC_nonNCC_DEG_filt.txt", quote = F, sep = "\t", row.names = F)


setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/v10/v10.2.1/lineage/chromVAR")
markers_motifs <- read.table("markers_volcano_NCC_nonNCC_subcluster_stage_chromvar.txt", header = T)

motif.names <- markers_motifs$feature
colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))
colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))
markers_rna$gene <- markers_rna$RNA.feature
markers_motifs$gene <- ConvertMotifID(gexatac_merge, id = motif.names)

markers_motifs$gene  <- paste(markers_motifs$gene  %>% str_sub(start = 1, end = 1),
                              markers_motifs$gene  %>% str_sub(start = 2) %>% str_to_lower(), sep = "")
# In Signac v1.10.0, when performing topTFs with inner join for motif uniqueness, some results are missing, so debug here.
markers_motifs$gene_unique <-markers_motifs$gene
markers_motifs$gene <- markers_motifs$gene %>% str_replace_all(pattern = "\\.\\d", replacement = "")


topTFs_lineage <- function(celltype, lineage, padj.cutoff = 1e-2, stageofinterest) { 
  ctmarkers_rna <- dplyr::filter(
    markers_rna, RNA.stage == stageofinterest, RNA.group == lineage, RNA.pval < padj.cutoff, RNA.logFC > 0, RNA.sub.cluster == celltype) %>%
    arrange(-RNA.auc)
  ctmarkers_motif <- dplyr::filter(
    markers_motifs, motif.stage == stageofinterest, motif.group == lineage, motif.pval < padj.cutoff, motif.logFC > 0, motif.sub.cluster == celltype) %>%
    arrange(-motif.auc)
  
  top_tfs <- inner_join(
    x = ctmarkers_rna[, c(17,2, 11, 6, 7)],
    y = ctmarkers_motif[, c(15, 2, 1, 11, 6, 7)], by = "gene"
  )
  top_tfs$avg_auc <- (top_tfs$RNA.auc + top_tfs$motif.auc) / 2
  top_tfs <- arrange(top_tfs, -avg_auc)
  return(top_tfs)
}




opts <- list()
opts[["tax_group"]] <-  "vertebrates"
opts[["all_versions"]] <- FALSE


PFMatrixList <- getMatrixSet(JASPAR2022, opts)
PFMatrixList

saveRDS(PFMatrixList, "PFMatrixList_JASPAR2022.rds")

setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/v10/v10.2.1/peaks/chromvar")
tfs <- readRDS("tfs.rds")
tfs2 <- readRDS("tfs2.rds")

celltype.names <- levels(gexatac_merge)
for(i in 1:length(celltype.names)){
  setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/v10/v10.2.1/lineage/chromVAR/")
  dir.create(celltype.names[i])
  setwd(celltype.names[i])
  
  tmp.genelist <- topTFs_lineage(celltype = celltype.names[i], 
                                 lineage = "nonNCC", stageofinterest = "E12.5",
                                 padj.cutoff = 5e-2)
  
  if(nrow(tmp.genelist) > 0){
    for(j in 1:nrow(tmp.genelist)){
      goi <- tmp.genelist$gene[j]
      DefaultAssay(gexatac_merge) <- "ATAC"
      
      motif.name <- tfs2 %>% dplyr::filter(name == str_to_upper(goi))
      motif.name <- motif.name$ID
      
     
      for(k in 1:length(motif.name)){
        DefaultAssay(gexatac_merge) <- "RNA"
        gene_plot <- FeaturePlot(gexatac_merge, features = goi, reduction = 'umap.rpca',order = T)
        DefaultAssay(gexatac_merge) <- "chromvar"
        motif_plot <- FeaturePlot(gexatac_merge, features = motif.name[k], min.cutoff = 0,  max.cutoff = "q80",
                                  cols = c("lightgrey", "darkred"), reduction = 'umap.rpca', order = T)
        gene_plot | motif_plot
        ggsave(paste(goi, "_",k, "_RNA_chormvar.pdf", sep = ""), height = 4, width = 8)
      } 
    } 
  } 
} 



topTFs_lineage(celltype = "0", lineage = "NCC", stageofinterest = "E11.5", padj.cutoff = 1e-3)



# heatmap
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/v10/v10.2.1/lineage//chromVAR/")
# NCC
df_toptf_NCC_E11.5 <- c()
celltype.names <- markers_rna$RNA.sub.cluster %>% unique() %>% sort()
for(i in 1:length(celltype.names)){
  tmp.df_toptfs <- head(topTFs_lineage(celltype = celltype.names[i], lineage = "NCC", stageofinterest = "E11.5", padj.cutoff = 1e-3),
                        n = 5)
    df_toptf_NCC_E11.5 <- rbind(df_toptf_NCC_E11.5, tmp.df_toptfs)
}
df_toptf_NCC_E12.5 <- c()
celltype.names <- markers_rna$RNA.sub.cluster %>% unique() %>% sort()
for(i in 1:length(celltype.names)){
  tmp.df_toptfs <- head(topTFs_lineage(celltype = celltype.names[i], lineage = "NCC", stageofinterest = "E12.5", padj.cutoff = 1e-3),
                        n = 5)
  df_toptf_NCC_E12.5 <- rbind(df_toptf_NCC_E12.5, tmp.df_toptfs)
}
df_toptf_NCC_E14.5 <- c()
celltype.names <- markers_rna$RNA.sub.cluster %>% unique() %>% sort()
for(i in 1:length(celltype.names)){
  tmp.df_toptfs <- head(topTFs_lineage(celltype = celltype.names[i], lineage = "NCC", stageofinterest = "E14.5", padj.cutoff = 1e-3),
                        n = 5)
  df_toptf_NCC_E14.5 <- rbind(df_toptf_NCC_E14.5, tmp.df_toptfs)
}
df_toptf_NCC_E17.5 <- c()
celltype.names <- markers_rna$RNA.sub.cluster %>% unique() %>% sort()
for(i in 1:length(celltype.names)){
  tmp.df_toptfs <- head(topTFs_lineage(celltype = celltype.names[i], lineage = "NCC", stageofinterest = "E17.5", padj.cutoff = 1e-3),
                        n = 5)
  df_toptf_NCC_E17.5 <- rbind(df_toptf_NCC_E17.5, tmp.df_toptfs)
}
# nonNCC
df_toptf_nonNCC_E11.5 <- c()
celltype.names <- markers_rna$RNA.sub.cluster %>% unique() %>% sort()
for(i in 1:length(celltype.names)){
  tmp.df_toptfs <- head(topTFs_lineage(celltype = celltype.names[i], lineage = "nonNCC", stageofinterest = "E11.5", padj.cutoff = 1e-3),
                        n = 5)
  df_toptf_nonNCC_E11.5 <- rbind(df_toptf_nonNCC_E11.5, tmp.df_toptfs)
}
df_toptf_nonNCC_E12.5 <- c()
celltype.names <- markers_rna$RNA.sub.cluster %>% unique() %>% sort()
for(i in 1:length(celltype.names)){
  tmp.df_toptfs <- head(topTFs_lineage(celltype = celltype.names[i], lineage = "nonNCC", stageofinterest = "E12.5", padj.cutoff = 1e-3),
                        n = 5)
  df_toptf_nonNCC_E12.5 <- rbind(df_toptf_nonNCC_E12.5, tmp.df_toptfs)
}
df_toptf_nonNCC_E14.5 <- c()
celltype.names <- markers_rna$RNA.sub.cluster %>% unique() %>% sort()
for(i in 1:length(celltype.names)){
  tmp.df_toptfs <- head(topTFs_lineage(celltype = celltype.names[i], lineage = "nonNCC", stageofinterest = "E14.5", padj.cutoff = 1e-3),
                        n = 5)
  df_toptf_nonNCC_E14.5 <- rbind(df_toptf_nonNCC_E14.5, tmp.df_toptfs)
}
df_toptf_nonNCC_E17.5 <- c()
celltype.names <- markers_rna$RNA.sub.cluster %>% unique() %>% sort()
for(i in 1:length(celltype.names)){
  tmp.df_toptfs <- head(topTFs_lineage(celltype = celltype.names[i], lineage = "nonNCC", stageofinterest = "E17.5", padj.cutoff = 1e-3),
                        n = 5)
  df_toptf_nonNCC_E17.5 <- rbind(df_toptf_nonNCC_E17.5, tmp.df_toptfs)
}


df_toptf <- rbind(df_toptf_NCC_E11.5, df_toptf_NCC_E12.5, df_toptf_NCC_E14.5, df_toptf_NCC_E17.5,
                  df_toptf_nonNCC_E11.5, df_toptf_nonNCC_E12.5, df_toptf_nonNCC_E14.5, df_toptf_nonNCC_E17.5)

write.table(df_toptf, "df_toptf_lineage.txt", quote = F,sep = "\t", row.names = F)

df_toptf　<- read.table("df_toptf_lineage.txt", header = T)

df_toptf$order <- rownames(df_toptf)
df_toptf$order <- as.numeric(df_toptf$order)
gene_order <- df_toptf %>% arrange(order) %>% pull(gene) %>% unique()
gene_order
gene_order_df <- as.data.frame(gene_order)
gene_order_df$order <- 1:nrow(gene_order_df)

setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/v10/v10.2.1/lineage/DEG/")
markers_rna <- read.table("markers_volcano_NCC_nonNCC_subcluster_stage_AUC.txt", header = T)
df_markers.rna <- dplyr::filter(markers_rna, gene %in% df_toptf$gene)

df_markers.rna$stage_lineage_cluster <- paste(df_markers.rna$stage,　df_markers.rna$group, df_markers.rna$sub.cluster, sep = "_")



setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/v10/v10.2.1/lineage/chromVAR")
markers_motifs <- read.table("markers_volcano_NCC_nonNCC_subcluster_stage_chromvar.txt", header = T)
motif.names <- markers_motifs$feature
markers_motifs$gene <- ConvertMotifID(gexatac_merge, id = motif.names)

markers_motifs$gene  <- paste(markers_motifs$gene  %>% str_sub(start = 1, end = 1),
                              markers_motifs$gene  %>% str_sub(start = 2) %>% str_to_lower(), sep = "")

markers_motifs$gene_unique <-markers_motifs$gene
markers_motifs$gene <- markers_motifs$gene %>% str_replace_all(pattern = "\\.\\d", replacement = "")

df_markers.motifs <- dplyr::filter(markers_motifs, gene %in% df_toptf$gene)
df_markers.motifs$motif.gene <- paste(df_markers.motifs$feature, df_markers.motifs$gene, sep ="_")

df_markers.motifs$stage_lineage_cluster <- paste(df_markers.motifs$stage,　df_markers.motifs$group, df_markers.motifs$sub.cluster, sep = "_")


df_markers.rna$gene <- as.factor(df_markers.rna$gene)
df_markers.rna$gene <- factor(df_markers.rna$gene,
                              levels = c(gene_order))
df_markers.rna$sub.cluster <- as.factor(df_markers.rna$sub.cluster)
df_markers.rna$sub.cluster  <- factor(df_markers.rna$sub.cluster,
                                    levels = c(celltype.names))
df_markers.motifs$sub.cluster <- as.factor(df_markers.motifs$sub.cluster)
df_markers.motifs$sub.cluster  <- factor(df_markers.motifs$sub.cluster,
                                         levels = c(celltype.names))
df_markers.motifs %>% head
tmp.df_toptf <- df_toptf %>% dplyr::select(gene, order)
tmp.df_toptf %>% head
df_markers.motifs$order <- NA

for(i in 1:nrow(df_markers.motifs)){
  for(j in 1:nrow(gene_order_df)){
    if(df_markers.motifs$gene[i] == gene_order_df$gene_order[j]){
      df_markers.motifs$order[i] <- gene_order_df$order[j]
    }
  }
}

df_markers.motifs <- df_markers.motifs %>% arrange(order)
df_markers.motifs %>% head
motif_order <- df_markers.motifs %>% arrange(order) %>% pull(motif.gene) %>% unique()
motif_order
df_markers.motifs$motif.gene <- as.factor(df_markers.motifs$motif.gene)
df_markers.motifs$motif.gene  <- factor(df_markers.motifs$motif.gene,
                                        levels = c(motif_order))


# logFC
p1 <- ggplot(df_markers.rna, aes(x=stage_lineage_cluster, y=gene, fill = logFC)) +
  geom_tile() +
  # scale_fill_gradient(low = "white", high = "red")
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  
  ggsave("hetamp_RNA.top5.pdf", height = 6, width = 12)


p2 <- ggplot(df_markers.motifs, aes(x=stage_lineage_cluster, y=motif.gene, fill = logFC)) +
  geom_tile() +
  # scale_fill_gradient(low = "white", high = "red")
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "motif.logFC", 
                       limits = c(-3, 3), oob = scales::squish  
                       ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  
ggsave("hetamp_motif.top5.pdf", height = 6, width = 12)


p1 / p2
ggsave("hetamp_merge.pdf", height = 12, width = 12)


write.table(df_markers.motifs, "df_markers.motifs_ggplot.txt", quote = F, sep = "\t", row.names = F)
write.table(df_markers.rna, "df_markers.rna_ggplot.txt", quote = F, sep = "\t", row.names = F)

ggplot(df_markers.motifs, aes(x=stage_lineage_cluster, y = avgExpr)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave("scatter_motif.pdf", height = 4, width = 12)

ggplot(df_markers.rna, aes(x=stage_lineage_cluster, y = avgExpr)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave("scatter_rna.pdf", height = 4, width = 12)


df_markers.motifs %>% dplyr::filter(avgExpr > 5) %>% View
df_markers.motifs %>%View



# avgExpr
p1 <- ggplot(df_markers.motifs, aes(x=stage_lineage_cluster, y=motif.gene, fill = avgExpr)) +
  geom_tile() +
  # scale_fill_gradient(low = "white", high = "red")
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "motif.avgExpr", 
                       limits = c(-3, 3), oob = scales::squish  
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  
ggsave("hetamp_motif.top5.pdf", height = 6, width = 12)

p2 <- ggplot(df_markers.rna, aes(x=stage_lineage_cluster, y=gene, fill = avgExpr)) +
  geom_tile() +
  # scale_fill_gradient(low = "white", high = "red")
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  
ggsave("hetamp_RNA.top5.pdf", height = 6, width = 12)

p1 / p2
ggsave("hetamp_merge.pdf", height = 12, width = 12)





setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v10/v10.2.1/lineage/chromVAR")
df_markers.rna <- read_tsv("df_markers.rna_ggplot.txt")
df_markers.motifs <- read_tsv("df_markers.motifs_ggplot.txt")

setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v10/v10.2.1/lineage/chromVAR/250611_lineage_chromVAR")
# 順番の調整
df_markers.rna$gene <- as.factor(df_markers.rna$gene)
df_markers.rna$gene <- factor(df_markers.rna$gene,
                              levels = c(gene_order))
df_markers.motifs <- df_markers.motifs %>% arrange(order)
df_markers.motifs %>% head
motif_order <- df_markers.motifs %>% arrange(order) %>% pull(motif.gene) %>% unique()
motif_order
df_markers.motifs$motif.gene <- as.factor(df_markers.motifs$motif.gene)
df_markers.motifs$motif.gene  <- factor(df_markers.motifs$motif.gene,
                                        levels = c(motif_order))

df_markers.rna_NCC <- df_markers.rna %>% dplyr::filter(stage %in% c("E11.5", "E12.5"), group == "NCC")
df_markers.motifs_NCC <- df_markers.motifs %>% dplyr::filter(stage %in% c("E11.5", "E12.5"),  group == "NCC")
df_markers.rna_nonNCC <- df_markers.rna %>% dplyr::filter(stage %in% c("E11.5", "E12.5"), group == "nonNCC")
df_markers.motifs_nonNCC <- df_markers.motifs %>% dplyr::filter(stage %in% c("E11.5", "E12.5"),  group == "nonNCC")

df_markers.rna_NCC$stage_lineage_cluster  <- as.factor(df_markers.rna_NCC$stage_lineage_cluster)
df_markers.rna_NCC$stage_lineage_cluster <- factor(df_markers.rna_NCC$stage_lineage_cluster,
                                     levels = c("E11.5_NCC_0","E11.5_NCC_1", "E11.5_NCC_2", "E11.5_NCC_3", "E11.5_NCC_4",
                                                "E11.5_NCC_6","E11.5_NCC_15", "E11.5_NCC_20",
                                                "E12.5_NCC_0","E12.5_NCC_1", "E12.5_NCC_2", "E12.5_NCC_3", "E12.5_NCC_4",
                                                "E12.5_NCC_6","E12.5_NCC_15", "E12.5_NCC_20"
                                                ))
df_markers.rna_nonNCC$stage_lineage_cluster  <- as.factor(df_markers.rna_nonNCC$stage_lineage_cluster)
df_markers.rna_nonNCC$stage_lineage_cluster <- factor(df_markers.rna_nonNCC$stage_lineage_cluster,
                                                         levels = c("E11.5_nonNCC_0","E11.5_nonNCC_1", "E11.5_nonNCC_2", "E11.5_nonNCC_3", "E11.5_nonNCC_4",
                                                "E11.5_nonNCC_6","E11.5_nonNCC_15", "E11.5_nonNCC_20",
                                                "E12.5_nonNCC_0","E12.5_nonNCC_1", "E12.5_nonNCC_2", "E12.5_nonNCC_3", "E12.5_nonNCC_4",
                                                "E12.5_nonNCC_6","E12.5_nonNCC_15", "E12.5_nonNCC_20"
                                                         ))

df_markers.motifs_NCC$stage_lineage_cluster  <- as.factor(df_markers.motifs_NCC$stage_lineage_cluster)
df_markers.motifs_NCC$stage_lineage_cluster <- factor(df_markers.motifs_NCC$stage_lineage_cluster,
                                                   levels = c("E11.5_NCC_0","E11.5_NCC_1", "E11.5_NCC_2", "E11.5_NCC_3", "E11.5_NCC_4",
                                                              "E11.5_NCC_6","E11.5_NCC_15", "E11.5_NCC_20",
                                                              "E12.5_NCC_0","E12.5_NCC_1", "E12.5_NCC_2", "E12.5_NCC_3", "E12.5_NCC_4",
                                                              "E12.5_NCC_6","E12.5_NCC_15", "E12.5_NCC_20"
                                                   ))
df_markers.motifs_nonNCC$stage_lineage_cluster  <- as.factor(df_markers.motifs_nonNCC$stage_lineage_cluster)
df_markers.motifs_nonNCC$stage_lineage_cluster <- factor(df_markers.motifs_nonNCC$stage_lineage_cluster,
                                                      levels = c("E11.5_nonNCC_0","E11.5_nonNCC_1", "E11.5_nonNCC_2", "E11.5_nonNCC_3", "E11.5_nonNCC_4",
                                                                 "E11.5_nonNCC_6","E11.5_nonNCC_15", "E11.5_nonNCC_20",
                                                                 "E12.5_nonNCC_0","E12.5_nonNCC_1", "E12.5_nonNCC_2", "E12.5_nonNCC_3", "E12.5_nonNCC_4",
                                                                 "E12.5_nonNCC_6","E12.5_nonNCC_15", "E12.5_nonNCC_20"
                                                      ))


p1 <- ggplot(df_markers.rna_NCC, aes(x=stage_lineage_cluster, y=gene, fill = avgExpr)) +
  geom_tile() +
  # scale_fill_gradient(low = "white", high = "red")
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),                 
        axis.ticks.x = element_blank())  +
  labs(y=NULL) +
  labs(x=NULL)

p2 <- ggplot(df_markers.motifs_NCC, aes(x=stage_lineage_cluster, y=motif.gene, fill = avgExpr)) +
  geom_tile() +
  # scale_fill_gradient(low = "white", high = "red")
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "motif.avgExpr", 
                       limits = c(-3, 3), oob = scales::squish  
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none") +
  labs(x=NULL)

p3 <- ggplot(df_markers.rna_nonNCC, aes(x=stage_lineage_cluster, y=gene, fill = avgExpr)) +
  geom_tile() +
  # scale_fill_gradient(low = "white", high = "red")
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme(axis.text.y = element_blank(),                  
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),                 
        axis.ticks.x = element_blank()) +
  labs(y=NULL) +
  labs(x=NULL)

p4 <- ggplot(df_markers.motifs_nonNCC, aes(x=stage_lineage_cluster, y=motif.gene, fill = avgExpr)) +
  geom_tile() +
  # scale_fill_gradient(low = "white", high = "red")
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "motif.avgExpr", 
                       limits = c(-3, 3), oob = scales::squish  
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
           axis.text.y = element_blank(),                 
           axis.ticks.y = element_blank()) +
  labs(y=NULL)+
  labs(x=NULL)

(p1+p3) / (p2+p4)

ggsave("hetamp_E11.5_12.5_merge.pdf", height = 8, width = 8)






setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v10/v10.2.1/lineage/chromVAR")
df_markers.rna <- read_tsv("df_markers.rna_ggplot.txt")
df_markers.motifs <- read_tsv("df_markers.motifs_ggplot.txt")

setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v10/v10.2.1/lineage/chromVAR/250611_lineage_chromVAR")

df_markers.rna$gene <- as.factor(df_markers.rna$gene)
df_markers.rna$gene <- factor(df_markers.rna$gene,
                              levels = c(gene_order))
df_markers.motifs <- df_markers.motifs %>% arrange(order)
df_markers.motifs %>% head
motif_order <- df_markers.motifs %>% arrange(order) %>% pull(motif.gene) %>% unique()
motif_order
df_markers.motifs$motif.gene <- as.factor(df_markers.motifs$motif.gene)
df_markers.motifs$motif.gene  <- factor(df_markers.motifs$motif.gene,
                                        levels = c(motif_order))

df_markers.rna_NCC <- df_markers.rna %>% dplyr::filter(stage %in% c("E14.5", "E17.5"), group == "NCC")
df_markers.motifs_NCC <- df_markers.motifs %>% dplyr::filter(stage %in% c("E14.5", "E17.5"),  group == "NCC")
df_markers.rna_nonNCC <- df_markers.rna %>% dplyr::filter(stage %in% c("E14.5", "E17.5"), group == "nonNCC")
df_markers.motifs_nonNCC <- df_markers.motifs %>% dplyr::filter(stage %in% c("E14.5", "E17.5"),  group == "nonNCC")

df_markers.rna_NCC$stage_lineage_cluster  <- as.factor(df_markers.rna_NCC$stage_lineage_cluster)
df_markers.rna_NCC$stage_lineage_cluster <- factor(df_markers.rna_NCC$stage_lineage_cluster,
                                                   levels = c("E14.5_NCC_0","E14.5_NCC_1", "E14.5_NCC_2", "E14.5_NCC_3", "E14.5_NCC_4",
                                                              "E14.5_NCC_6","E14.5_NCC_15", "E14.5_NCC_20",
                                                              "E17.5_NCC_0","E17.5_NCC_1", "E17.5_NCC_2", "E17.5_NCC_3", "E17.5_NCC_4",
                                                              "E17.5_NCC_6","E17.5_NCC_15", "E17.5_NCC_20"
                                                   ))
df_markers.rna_nonNCC$stage_lineage_cluster  <- as.factor(df_markers.rna_nonNCC$stage_lineage_cluster)
df_markers.rna_nonNCC$stage_lineage_cluster <- factor(df_markers.rna_nonNCC$stage_lineage_cluster,
                                                      levels = c("E14.5_nonNCC_0","E14.5_nonNCC_1", "E14.5_nonNCC_2", "E14.5_nonNCC_3", "E14.5_nonNCC_4",
                                                                 "E14.5_nonNCC_6","E14.5_nonNCC_15", "E14.5_nonNCC_20",
                                                                 "E17.5_nonNCC_0","E17.5_nonNCC_1", "E17.5_nonNCC_2", "E17.5_nonNCC_3", "E17.5_nonNCC_4",
                                                                 "E17.5_nonNCC_6","E17.5_nonNCC_15", "E17.5_nonNCC_20"
                                                      ))

df_markers.motifs_NCC$stage_lineage_cluster  <- as.factor(df_markers.motifs_NCC$stage_lineage_cluster)
df_markers.motifs_NCC$stage_lineage_cluster <- factor(df_markers.motifs_NCC$stage_lineage_cluster,
                                                      levels = c("E14.5_NCC_0","E14.5_NCC_1", "E14.5_NCC_2", "E14.5_NCC_3", "E14.5_NCC_4",
                                                                 "E14.5_NCC_6","E14.5_NCC_15", "E14.5_NCC_20",
                                                                 "E17.5_NCC_0","E17.5_NCC_1", "E17.5_NCC_2", "E17.5_NCC_3", "E17.5_NCC_4",
                                                                 "E17.5_NCC_6","E17.5_NCC_15", "E17.5_NCC_20"
                                                      ))
df_markers.motifs_nonNCC$stage_lineage_cluster  <- as.factor(df_markers.motifs_nonNCC$stage_lineage_cluster)
df_markers.motifs_nonNCC$stage_lineage_cluster <- factor(df_markers.motifs_nonNCC$stage_lineage_cluster,
                                                         levels = c("E14.5_nonNCC_0","E14.5_nonNCC_1", "E14.5_nonNCC_2", "E14.5_nonNCC_3", "E14.5_nonNCC_4",
                                                                    "E14.5_nonNCC_6","E14.5_nonNCC_15", "E14.5_nonNCC_20",
                                                                    "E17.5_nonNCC_0","E17.5_nonNCC_1", "E17.5_nonNCC_2", "E17.5_nonNCC_3", "E17.5_nonNCC_4",
                                                                    "E17.5_nonNCC_6","E17.5_nonNCC_15", "E17.5_nonNCC_20"
                                                         ))


p1 <- ggplot(df_markers.rna_NCC, aes(x=stage_lineage_cluster, y=gene, fill = avgExpr)) +
  geom_tile() +
  # scale_fill_gradient(low = "white", high = "red")
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),                 
        axis.ticks.x = element_blank())  +
  labs(x=NULL)

p2 <- ggplot(df_markers.motifs_NCC, aes(x=stage_lineage_cluster, y=motif.gene, fill = avgExpr)) +
  geom_tile() +
  # scale_fill_gradient(low = "white", high = "red")
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "motif.avgExpr", 
                       limits = c(-3, 3), oob = scales::squish  
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none") +
  labs(x=NULL)

p3 <- ggplot(df_markers.rna_nonNCC, aes(x=stage_lineage_cluster, y=gene, fill = avgExpr)) +
  geom_tile() +
  # scale_fill_gradient(low = "white", high = "red")
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme(axis.text.y = element_blank(),                 
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),                  
        axis.ticks.x = element_blank()) +
  labs(y=NULL) +
  labs(x=NULL)

p4 <- ggplot(df_markers.motifs_nonNCC, aes(x=stage_lineage_cluster, y=motif.gene, fill = avgExpr)) +
  geom_tile() +
  # scale_fill_gradient(low = "white", high = "red")
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "motif.avgExpr", 
                       limits = c(-3, 3), oob = scales::squish  
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_blank(),                  
        axis.ticks.y = element_blank()) +
  labs(y=NULL)+
  labs(x=NULL)

(p1+p3) / (p2+p4)

ggsave("hetamp_E14.5_17.5_merge.pdf", height = 8, width = 8)
