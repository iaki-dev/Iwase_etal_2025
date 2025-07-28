
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


setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered")
NCC_merge <- readRDS("NCC_merge2_RPCA_velo_region.rds")
DimPlot(NCC_merge, reduction = "umap.rpca", group.by = "region", label = TRUE)



coi <- NCC_merge@meta.data %>% dplyr::filter(!seurat_clusters %in% c("24","17","20","12","9","25","26","1","3","19"))
mes <- NCC_merge[, rownames(coi)]
DimPlot(mes, reduction = "umap.rpca", label = T)
DimPlot(mes, reduction = "umap.rpca", group.by = "region2",label = F, pt.size = 0.5) 

setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/mes/region_specific")
ggsave("umap.mes_region_all.pdf", height = 4, width = 5)

coi <- mes@meta.data %>% dplyr::filter(sample %in% c("E11.5_EYFP","E12.5_EYFP","E14.5_EYFP","E17.5_EYFP"))
mes <- mes[, rownames(coi)]
Idents(mes) <- mes@meta.data$region
DimPlot(mes, reduction = "umap.rpca", group.by =  "region", label = F)
ggsave("umap.mes_region.pdf", height = 4, width = 5)


setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/mes/region_specific/chromVAR")
markers_motifs <- presto:::wilcoxauc.Seurat(X = mes, group_by = 'region2', assay = 'data', seurat_assay = 'chromvar')
# write.table(markers_motifs, "markers_NCC_motifs.txt", sep = "\t", quote = F, row.names = F)
markers_rna <- presto:::wilcoxauc.Seurat(X = mes, group_by = 'region2', assay = 'data', seurat_assay = 'RNA')
# write.table(markers_rna, "markers_NCC_rna.txt", sep = "\t", quote = F, row.names = F)
# markers_atac <- presto:::wilcoxauc.Seurat(X = NCC_merge, group_by = 'sub.cluster', assay = 'data', seurat_assay = 'ATAC')
# write.table(markers_atac, "markers_atac.txt", sep = "\t", quote = F, row.names = F)


# ロード
markers_rna <- read_tsv("markers_NCC_rna.txt")
markers_motifs <- read_tsv("markers_NCC_motifs.txt")

DefaultAssay(mes) <- "ATAC"
motif.names <- markers_motifs$feature
colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))
colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))
markers_rna$gene <- markers_rna$RNA.feature
markers_motifs$gene <- ConvertMotifID(mes, id = motif.names)

markers_motifs$gene  <- paste(markers_motifs$gene  %>% str_sub(start = 1, end = 1),
                              markers_motifs$gene  %>% str_sub(start = 2) %>% str_to_lower(), sep = "")
markers_motifs$gene_unique <-markers_motifs$gene
markers_motifs$gene <- markers_motifs$gene %>% str_replace_all(pattern = "\\.\\d", replacement = "")




topTFs <- function(celltype, padj.cutoff = 1e-2) { 
  ctmarkers_rna <- dplyr::filter(
    markers_rna, RNA.group == celltype, RNA.pval < padj.cutoff, RNA.logFC > 0) %>%
    arrange(-RNA.auc)
  ctmarkers_motif <- dplyr::filter(
    markers_motifs, motif.group == celltype, motif.pval < padj.cutoff, motif.logFC > 0) %>%
    arrange(-motif.auc)

  top_tfs <- inner_join(
    x = ctmarkers_rna[, c(2, 11, 6, 7)],
    y = ctmarkers_motif[, c(2, 1, 11, 6, 7)], by = "gene"
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

# saveRDS(PFMatrixList, "PFMatrixList_JASPAR2022.rds")


setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v10/v10.2.1/peaks/CARTA_v0.7.1")
tfs <- readRDS("tfs.rds")
tfs2 <- readRDS("tfs2.rds")


Idents(mes) <- mes@meta.data$region2
celltype.names <- levels(mes)
for(i in 1:length(celltype.names)){
  setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/mes/region_specific/chromVAR")
  dir.create(celltype.names[i])
  setwd(celltype.names[i])

  tmp.genelist <- topTFs(celltype = celltype.names[i], padj.cutoff = 5e-2)

  if(nrow(tmp.genelist) > 0){
    for(j in 1:nrow(tmp.genelist)){
      goi <- tmp.genelist$gene[j]
      DefaultAssay(mes) <- "ATAC"
     
      motif.name <- tfs2 %>% dplyr::filter(name == str_to_upper(goi))
      motif.name <- motif.name$ID

     
      for(k in 1:length(motif.name)){
        DefaultAssay(mes) <- "RNA"
        gene_plot <- FeaturePlot(mes, features = goi, reduction = 'umap.rpca')
        DefaultAssay(mes) <- "chromvar"
        motif_plot <- FeaturePlot(mes, features = motif.name[k], min.cutoff = 0,  max.cutoff = "q80",
                                  cols = c("lightgrey", "darkred"), reduction = 'umap.rpca')
        gene_plot | motif_plot
        ggsave(paste(goi, "_",k, "_RNA_chormvar.pdf", sep = ""), height = 4, width = 8)
      } 
    } 
  } 
} 



topTFs(celltype = "Pharyngeal", padj.cutoff = 5e-1)
topTFs(celltype = "4", padj.cutoff = 1e-3)
topTFs(celltype = "9_0", padj.cutoff = 1e-1)


head(topTFs("0"), 5)
head(topTFs("13"), 3)



setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/mes/region_specific/chromVAR")
df_toptf <- c()
celltype.names <- levels(mes)
for(i in 1:length(celltype.names)){
  tmp.df_toptfs <- head(topTFs(celltype.names[i]), 20)
  df_toptf <- rbind(df_toptf, tmp.df_toptfs)
}
write.table(df_toptf, "df_toptf.txt", quote = F,sep = "\t", row.names = F)


df_toptf <- df_toptf %>% dplyr::filter(gene %in% c("Meis1","Meis2","Hoxa2", "Hoxb3", "Tbx20","Gata4", "Mef2c"))

df_toptf$order <- rownames(df_toptf)
df_toptf$order <- as.numeric(df_toptf$order)
gene_order <- c("Mef2c","Tbx20","Gata4", "Meis2", "Meis1","Hoxb3", "Hoxa2")

df_markers.rna <- dplyr::filter(markers_rna, gene %in% df_toptf$gene)
df_markers.motifs <- dplyr::filter(markers_motifs, gene %in% df_toptf$gene)
df_markers.motifs$motif.gene <- paste(df_markers.motifs$motif.feature, df_markers.motifs$gene, sep ="_")

# 順番の調整
df_markers.rna$gene <- as.factor(df_markers.rna$gene)
df_markers.rna$gene <- factor(df_markers.rna$gene,
                              levels = c(gene_order))
df_markers.rna$RNA.group <- as.factor(df_markers.rna$RNA.group)
df_markers.rna$RNA.group  <- factor(df_markers.rna$RNA.group,
                                    levels = c(celltype.names))
df_markers.motifs$motif.group <- as.factor(df_markers.motifs$motif.group)
df_markers.motifs$motif.group  <- factor(df_markers.motifs$motif.group,
                                    levels = c(celltype.names))
df_markers.motifs %>% head
tmp.df_toptf <- df_toptf %>% select(gene, order)
tmp.df_toptf %>% head
df_markers.motifs <- inner_join(df_markers.motifs, tmp.df_toptf, by = "gene")
df_markers.motifs <- df_markers.motifs %>% arrange(order)
df_markers.motifs %>% head
# motif_order
motif_order <- c("MA0497.1_Mef2c" ,"MA0689.1_Tbx20","MA0482.2_Gata4",
                 "MA0774.1_Meis2", "MA1640.1_Meis2" ,"MA0498.2_Meis1", "MA1639.1_Meis1",  "MA0903.1_Hoxb3","MA0900.2_Hoxa2"
                 )
df_markers.motifs$motif.gene <- as.factor(df_markers.motifs$motif.gene)
df_markers.motifs$motif.gene  <- factor(df_markers.motifs$motif.gene,
                                         levels = c(motif_order))

df_markers.rna$gene
ggplot(df_markers.rna, aes(x=RNA.group, y=gene, fill = RNA.avgExpr)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 1), oob = scales::squish) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())  
ggsave("hetamp_RNA.interest.pdf", height = 4, width = 4)

df_markers.motifs <- df_markers.motifs %>% mutate(limit_motif.logFC = case_when(
  motif.logFC > 3 ~ 3,
  motif.logFC < -3 ~ -3,
))
for(i in 1:nrow(df_markers.motifs)){
  if(is.na(df_markers.motifs$limit_motif.logFC[i]) == T){
    df_markers.motifs$limit_motif.logFC[i] <- df_markers.motifs$motif.logFC[i]
  }
}


ggplot(df_markers.motifs, aes(x=motif.group, y=motif.gene, fill = motif.avgExpr)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,  
                       limits = c(-3, 3),  oob = scales::squish) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) 
ggsave("hetamp_motif.interest.pdf", height = 4, width = 5)






