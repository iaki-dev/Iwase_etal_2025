
library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)
library(tidyverse)
library(GenomicRanges)
library(future)
library(ggplot2)
library(patchwork)
library(BSgenome.Mmusculus.UCSC.mm10)
library(motifmatchr)
set.seed(1234)
library(data.table)
library(progress)

"%not.in%" <- Negate("%in%")

setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered")
NCC_merge <- readRDS("gexatac_merge_NCC_mes.rds")
DimPlot(gexatac_merge, reduction = "umap.rpca")

links.df <- data.frame(gexatac_merge[["ATAC"]]@links)
links.df <- links.df %>% dplyr::filter(score > 0)

Idents(gexatac_merge) <- gexatac_merge@meta.data$region2
markers_atac <- FindAllMarkers(gexatac_merge, assay = "ATAC")
write_tsv(markers_atac, "markers_atac.txt")


# Pharyngeal
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/v6/")
CARTA_Net <- read_tsv("CARTA_Net_conservation_score.txt")  
dim(CARTA_Net)
# [1] 950035     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(avg_conservation_score > 0.8)
dim(CARTA_Net)
# [1] 14612     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(peak %in% links.df$peak)
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/")
markers <- read_tsv("markers_region2.txt")
markers <- markers %>% dplyr::filter(avg_log2FC > 0.5,
                                     p_val < 1e-4)
markers_sub <- markers %>% dplyr::filter(cluster == "Pharyngeal")
df_avg_atac <- AverageExpression(object = gexatac_merge,  group.by = "region2")$ATAC %>% as.data.frame()
df_avg_atac <- df_avg_atac %>% dplyr::filter(`Cardiac-Cushion` > 0.5)
markers_atac <- read_tsv("markers_atac.txt")
markers_atac <- markers_atac %>% dplyr::filter(cluster == "Pharyngeal", avg_log2FC > 0)
df_avg_atac <- df_avg_atac[markers_atac$gene,]

df_avg_rna <- AverageExpression(object = gexatac_merge,  group.by = "region2")$RNA %>% as.data.frame()
df_avg_rna <- df_avg_rna %>% dplyr::filter(`Cardiac-Cushion` > 0.3)
CARTA_Net_Pharyngeal <- CARTA_Net %>% dplyr::filter(target %in% markers_sub$gene, 
                                                         TF %in% rownames(df_avg_rna),
                                                         peak %in% rownames(df_avg_atac))
CARTA_Net_Pharyngeal$name_ID <- paste(CARTA_Net_Pharyngeal$TF, CARTA_Net_Pharyngeal$MotifID, sep = "_")
CARTA_Net_Pharyngeal$region <- "Pharyngeal"


# Transitional
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/v6/")
CARTA_Net <- read_tsv("CARTA_Net_conservation_score.txt")  
dim(CARTA_Net)
# [1] 950035     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(avg_conservation_score > 0.8)
dim(CARTA_Net)
# [1] 141712     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(peak %in% links.df$peak)
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/")
markers <- read_tsv("markers_region2.txt")
markers <- markers %>% dplyr::filter(avg_log2FC > 0.5,
                                     p_val < 1e-4)
markers_sub <- markers %>% dplyr::filter(cluster == "Transitional")
df_avg_atac <- AverageExpression(object = gexatac_merge,  group.by = "region2")$ATAC %>% as.data.frame()
df_avg_atac <- df_avg_atac %>% dplyr::filter(`Cardiac-Cushion` > 0.5)
markers_atac <- read_tsv("markers_atac.txt")
markers_atac <- markers_atac %>% dplyr::filter(cluster == "Transitional", avg_log2FC > 0)
df_avg_atac <- df_avg_atac[markers_atac$gene,]

df_avg_rna <- AverageExpression(object = gexatac_merge,  group.by = "region2")$RNA %>% as.data.frame()
df_avg_rna <- df_avg_rna %>% dplyr::filter(`Cardiac-Cushion` > 0.3)

CARTA_Net_Transitional <- CARTA_Net %>% dplyr::filter(target %in% markers_sub$gene, 
                                                    TF %in% rownames(df_avg_rna),
                                                    peak %in% rownames(df_avg_atac))
CARTA_Net_Transitional$name_ID <- paste(CARTA_Net_Transitional$TF, CARTA_Net_Transitional$MotifID, sep = "_")
CARTA_Net_Transitional$region <- "Transitional"


# Cardiac_Cushion
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/v6/")
CARTA_Net <- read_tsv("CARTA_Net_conservation_score.txt")  
dim(CARTA_Net)
# [1] 950035     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(avg_conservation_score > 0.8)
dim(CARTA_Net)
# [1] 14612     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(peak %in% links.df$peak)
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/")
markers <- read_tsv("markers_region2.txt")
markers <- markers %>% dplyr::filter(avg_log2FC > 0.5,
                                     p_val < 1e-4)
markers_sub <- markers %>% dplyr::filter(cluster == "Cardiac_Cushion")
df_avg_atac <- AverageExpression(object = gexatac_merge,  group.by = "region2")$ATAC %>% as.data.frame()
df_avg_atac <- df_avg_atac %>% dplyr::filter(`Cardiac-Cushion` > 0.5)
markers_atac <- read_tsv("markers_atac.txt")
markers_atac <- markers_atac %>% dplyr::filter(cluster == "Cardiac_Cushion", avg_log2FC > 0)
df_avg_atac <- df_avg_atac[markers_atac$gene,]

df_avg_rna <- AverageExpression(object = gexatac_merge,  group.by = "region2")$RNA %>% as.data.frame()
df_avg_rna <- df_avg_rna %>% dplyr::filter(`Cardiac-Cushion` > 0.3)

CARTA_Net_Cardiac_Cushion <- CARTA_Net %>% dplyr::filter(target %in% markers_sub$gene, 
                                                    TF %in% rownames(df_avg_rna),
                                                    peak %in% rownames(df_avg_atac))
CARTA_Net_Cardiac_Cushion$name_ID <- paste(CARTA_Net_Cardiac_Cushion$TF, CARTA_Net_Cardiac_Cushion$MotifID, sep = "_")
CARTA_Net_Cardiac_Cushion$region <- "Cardiac_Cushion"


# Subvalvular
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/v6/")
CARTA_Net <- read_tsv("CARTA_Net_conservation_score.txt")  
dim(CARTA_Net)
# [1] 950035     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(avg_conservation_score > 0.8)
dim(CARTA_Net)
# [1] 14612     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(peak %in% links.df$peak)
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/")
markers <- read_tsv("markers_region2.txt")
markers <- markers %>% dplyr::filter(avg_log2FC > 0.5,
                                     p_val < 1e-4)
markers_sub <- markers %>% dplyr::filter(cluster == "Subvalvular")

df_avg_atac <- AverageExpression(object = gexatac_merge,  group.by = "region2")$ATAC %>% as.data.frame()
df_avg_atac <- df_avg_atac %>% dplyr::filter(`Cardiac-Cushion` > 0.5)
markers_atac <- read_tsv("markers_atac.txt")
markers_atac <- markers_atac %>% dplyr::filter(cluster == "Subvalvular", avg_log2FC > 0)
df_avg_atac <- df_avg_atac[markers_atac$gene,]
                           
df_avg_rna <- AverageExpression(object = gexatac_merge,  group.by = "region2")$RNA %>% as.data.frame()
df_avg_rna <- df_avg_rna %>% dplyr::filter(`Cardiac-Cushion` > 0.3)

CARTA_Net_Subvalvular <- CARTA_Net %>% dplyr::filter(target %in% markers_sub$gene, 
                                                    TF %in% rownames(df_avg_rna),
                                                    peak %in% rownames(df_avg_atac))
CARTA_Net_Subvalvular$name_ID <- paste(CARTA_Net_Subvalvular$TF, CARTA_Net_Subvalvular$MotifID, sep = "_")
CARTA_Net_Subvalvular$region <- "Subvalvular"


# AP_septum
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/v6/")
CARTA_Net <- read_tsv("CARTA_Net_conservation_score.txt")  
dim(CARTA_Net)
# [1] 950035     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(avg_conservation_score > 0.8)
dim(CARTA_Net)
# [1] 14612     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(peak %in% links.df$peak)
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/")
markers <- read_tsv("markers_region2.txt")
markers <- markers %>% dplyr::filter(avg_log2FC > 0.5,
                                     p_val < 1e-4)
markers_sub <- markers %>% dplyr::filter(cluster == "AP_septum")
df_avg_atac <- AverageExpression(object = gexatac_merge,  group.by = "region2")$ATAC %>% as.data.frame()
df_avg_atac <- df_avg_atac %>% dplyr::filter(`Cardiac-Cushion` > 0.5)
markers_atac <- read_tsv("markers_atac.txt")
markers_atac <- markers_atac %>% dplyr::filter(cluster == "AP_septum", avg_log2FC > 0)
df_avg_atac <- df_avg_atac[markers_atac$gene,]

df_avg_rna <- AverageExpression(object = gexatac_merge,  group.by = "region2")$RNA %>% as.data.frame()
df_avg_rna <- df_avg_rna %>% dplyr::filter(`Cardiac-Cushion` > 0.3)

CARTA_Net_AP_septum <- CARTA_Net %>% dplyr::filter(target %in% markers_sub$gene, 
                                                    TF %in% rownames(df_avg_rna),
                                                    peak %in% rownames(df_avg_atac))
CARTA_Net_AP_septum$name_ID <- paste(CARTA_Net_AP_septum$TF, CARTA_Net_AP_septum$MotifID, sep = "_")
CARTA_Net_AP_septum$region <- "AP_septum"


# SMC_GA
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/v6/")
CARTA_Net <- read_tsv("CARTA_Net_conservation_score.txt")  
dim(CARTA_Net)
# [1] 950035     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(avg_conservation_score > 0.8)
dim(CARTA_Net)
# [1] 14612     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(peak %in% links.df$peak)
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/")
markers <- read_tsv("markers_region2.txt")
markers <- markers %>% dplyr::filter(avg_log2FC > 0.5,
                                     p_val < 1e-4)
markers_sub <- markers %>% dplyr::filter(cluster == "SMC_GA")
df_avg_atac <- AverageExpression(object = gexatac_merge,  group.by = "region2")$ATAC %>% as.data.frame()
df_avg_atac <- df_avg_atac %>% dplyr::filter(`Cardiac-Cushion` > 0.5)
markers_atac <- read_tsv("markers_atac.txt")
markers_atac <- markers_atac %>% dplyr::filter(cluster == "SMC_GA", avg_log2FC > 0)
df_avg_atac <- df_avg_atac[markers_atac$gene,]

df_avg_rna <- AverageExpression(object = gexatac_merge,  group.by = "region2")$RNA %>% as.data.frame()
df_avg_rna <- df_avg_rna %>% dplyr::filter(`Cardiac-Cushion` > 0.3)

CARTA_Net_SMC_GA <- CARTA_Net %>% dplyr::filter(target %in% markers_sub$gene, 
                                                    TF %in% rownames(df_avg_rna),
                                                    peak %in% rownames(df_avg_atac))
CARTA_Net_SMC_GA$name_ID <- paste(CARTA_Net_SMC_GA$TF, CARTA_Net_SMC_GA$MotifID, sep = "_")
CARTA_Net_SMC_GA$region <- "SMC_GA"


# SMC_DA
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/v6/")
CARTA_Net <- read_tsv("CARTA_Net_conservation_score.txt")  
dim(CARTA_Net)
# [1] 950035     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(avg_conservation_score > 0.8)
dim(CARTA_Net)
# [1] 14612     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(peak %in% links.df$peak)
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/")
markers <- read_tsv("markers_region2.txt")
markers <- markers %>% dplyr::filter(avg_log2FC > 0.5,
                                     p_val < 1e-4)
markers_sub <- markers %>% dplyr::filter(cluster == "SMC_DA")

df_avg_atac <- AverageExpression(object = gexatac_merge,  group.by = "region2")$ATAC %>% as.data.frame()
df_avg_atac <- df_avg_atac %>% dplyr::filter(`Cardiac-Cushion` > 0.5)
markers_atac <- read_tsv("markers_atac.txt")
markers_atac <- markers_atac %>% dplyr::filter(cluster == "SMC_DA", avg_log2FC > 0)
df_avg_atac <- df_avg_atac[markers_atac$gene,]

df_avg_rna <- AverageExpression(object = gexatac_merge,  group.by = "region2")$RNA %>% as.data.frame()
df_avg_rna <- df_avg_rna %>% dplyr::filter(`Cardiac-Cushion` > 0.3)

CARTA_Net_SMC_DA <- CARTA_Net %>% dplyr::filter(target %in% markers_sub$gene, 
                                                    TF %in% rownames(df_avg_rna),
                                                    peak %in% rownames(df_avg_atac))
CARTA_Net_SMC_DA$name_ID <- paste(CARTA_Net_SMC_DA$TF, CARTA_Net_SMC_DA$MotifID, sep = "_")
CARTA_Net_SMC_DA$region <- "SMC_DA"


# SMC_CA
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/v6/")
CARTA_Net <- read_tsv("CARTA_Net_conservation_score.txt")  
dim(CARTA_Net)
# [1] 950035     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(avg_conservation_score > 0.8)
dim(CARTA_Net)
# [1] 14612     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(peak %in% links.df$peak)
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/")
markers <- read_tsv("markers_region2.txt")
markers <- markers %>% dplyr::filter(avg_log2FC > 0.5,
                                     p_val < 1e-4)
markers_sub <- markers %>% dplyr::filter(cluster == "SMC_CA")
df_avg_atac <- AverageExpression(object = gexatac_merge,  group.by = "region2")$ATAC %>% as.data.frame()
df_avg_atac <- df_avg_atac %>% dplyr::filter(`Cardiac-Cushion` > 0.5)
markers_atac <- read_tsv("markers_atac.txt")
markers_atac <- markers_atac %>% dplyr::filter(cluster == "SMC_CA", avg_log2FC > 0)
df_avg_atac <- df_avg_atac[markers_atac$gene,]

df_avg_rna <- AverageExpression(object = gexatac_merge,  group.by = "region2")$RNA %>% as.data.frame()
df_avg_rna <- df_avg_rna %>% dplyr::filter(`Cardiac-Cushion` > 0.3)

CARTA_Net_SMC_CA <- CARTA_Net %>% dplyr::filter(target %in% markers_sub$gene, 
                                                    TF %in% rownames(df_avg_rna),
                                                    peak %in% rownames(df_avg_atac))
CARTA_Net_SMC_CA$name_ID <- paste(CARTA_Net_SMC_CA$TF, CARTA_Net_SMC_CA$MotifID, sep = "_")
CARTA_Net_SMC_CA$region <- "SMC_CA"


setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/v6/")
CARTA_Net_summarize <- bind_rows(CARTA_Net_Pharyngeal,CARTA_Net_Transitional,CARTA_Net_Cardiac_Cushion,
                                 CARTA_Net_Subvalvular,CARTA_Net_AP_septum,
                                 CARTA_Net_SMC_GA, CARTA_Net_SMC_DA, CARTA_Net_SMC_CA)
write_tsv(CARTA_Net_summarize, "CARTA_Net_summarize_filt_conservation.txt")




