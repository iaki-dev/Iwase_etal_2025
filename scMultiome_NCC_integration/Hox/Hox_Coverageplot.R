
library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)
library(tidyverse)
library(GenomicRanges)
library(future)
# library(DoubletFinder)
# library(SeuratData)
# library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(BSgenome.Mmusculus.UCSC.mm10)
library(motifmatchr)
# library(gtools) # 順列計算
set.seed(1234)
library(data.table)
library(progress)

"%not.in%" <- Negate("%in%")

setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered")
gexatac_merge <- readRDS("gexatac_merge_NCC_mes.rds")
DimPlot(gexatac_merge, reduction = "umap.rpca")

setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/Hox")
DefaultAssay(gexatac_merge) <- "ATAC"
da_peaks <- FindMarkers(
  object = gexatac_merge,
  ident.1 = c("Pharyngeal", "Transitional"),
  ident.2 = c("Cardiac_Cushion"),
  only.pos = TRUE
)

# get top differentially accessible peaks
top.da.peak <- da_peaks
top.da.peak$peak <- rownames(top.da.peak)
top.da.peak <- top.da.peak %>%
  separate(peak, into = c("seqnames", "start", "end"), sep = "-")
top.da.peak$start <- as.numeric(top.da.peak$start)
top.da.peak$end <- as.numeric(top.da.peak$end)

pharyngeal_peak <- top.da.peak %>% dplyr::filter(avg_log2FC > 0)
cushion_peak <- top.da.peak %>% dplyr::filter(avg_log2FC < 0)

# Hoxa2 - Hoxa5 (Hoxa3)
HoxA_peak <- pharyngeal_peak %>% dplyr::filter(seqnames == "chr6",
                                               start > 52162500,
                                               end < 52220000)
write.table(HoxA_peak, "HoxA_peak.txt", quote = F, sep = "\t", col.names = NA)
# Hoxb2 - Hoxa5 (Hoxa3)
HoxB_peak <- pharyngeal_peak %>% dplyr::filter(seqnames == "chr11",
                                               start > 96302000,
                                               end < 96355000)
write.table(HoxB_peak, "HoxB_peak.txt", quote = F, sep = "\t", col.names = NA)

enriched.motifs <- FindMotifs(
  object = gexatac_merge,
  features = c(rownames(HoxA_peak), rownames(HoxB_peak))
)
write.table(enriched.motifs, "enriched.motif_Hox_region_pharyngeal.txt", quote = F, sep = "\t", row.names = F)


# RNA in Pharyngeal
df_avg_rna <- AverageExpression(object = gexatac_merge,  group.by = "region2")$RNA %>% as.data.frame()
df_avg_rna <- df_avg_rna %>% dplyr::filter(Pharyngeal > 0.3)
exp.genes <- toupper(rownames(df_avg_rna))

# enriched.motifsのmortif.nameが重複阻止されているので修正する
enriched.motifs <- enriched.motifs %>% mutate(motif.name_correct =　str_replace(motif.name, "\\..*$", ""))
enriched.motifs <- enriched.motifs %>% dplyr::filter(motif.name_correct %in% exp.genes)
write.table(enriched.motifs, "enriched.motif_Hox_region_pharyngeal_expressed.txt", quote = F, sep = "\t", col.names = NA)

MotifPlot(
  object = gexatac_merge,
  motifs = head(rownames(enriched.motifs))
)
ggsave("Motif_Hox_region_pharyngeal_expressed.pdf", width = 8, height = 6)



# Meis1/2の結合サイトがHoxの領域にあるか？
tmp.df <- df_all %>% dplyr::filter(motif_id %in% c("MA1640.1", "MA1639.1"),
                                   seqnames == "chr6",
                                   start > 52162500,
                                   end < 52220000) 
tmp.peaks <- paste(tmp.df$seqnames, tmp.df$start, tmp.df$end, sep = "-")
ranges.show <- StringToGRanges(tmp.peaks)
ranges.show$color <- "red"

ranges.show2 <- StringToGRanges(rownames(HoxA_peak))
ranges.show2$color <- "green"

tmp.df <- df_all %>% dplyr::filter(motif_id %in% c("MA0774.1", "MA0498.2"),
                                   seqnames == "chr6",
                                   start > 52162500,
                                   end < 52220000) 
tmp.peaks <- paste(tmp.df$seqnames, tmp.df$start, tmp.df$end, sep = "-")
ranges.show3 <- StringToGRanges(tmp.peaks)
ranges.show3$color <- "blue"

ranges.all <- c(ranges.show2, ranges.show, ranges.show3)
CoveragePlot(gexatac_merge, region = "Hoxa3", extend.upstream = 20000, extend.downstream = 20000,
             features = c("Hoxa2","Hoxa3","Hoxa4", "Hoxa5"),
             region.highlight = ranges.all)
ggsave("Coverage_HoxA_Meis1and2_both_sites.pdf", width = 20, height = 8)


tmp.df <- df_all %>% dplyr::filter(motif_id %in% c("MA1640.1", "MA1639.1"),
                                  seqnames == "chr11",
                                  start > 96302000,
                                  end < 96355000) 
tmp.peaks <- paste(tmp.df$seqnames, tmp.df$start, tmp.df$end, sep = "-")
ranges.show <- StringToGRanges(tmp.peaks)
ranges.show$color <- "red"


ranges.show2 <- StringToGRanges(rownames(HoxB_peak))
ranges.show2$color <- "green"

tmp.df <- df_all %>% dplyr::filter(motif_id %in% c("MA0774.1", "MA0498.2"),
                                   seqnames == "chr11",
                                   start > 96302000,
                                   end < 96355000) 
tmp.peaks <- paste(tmp.df$seqnames, tmp.df$start, tmp.df$end, sep = "-")
ranges.show3 <- StringToGRanges(tmp.peaks)
ranges.show3$color <- "blue"

ranges.all <- c(ranges.show2, ranges.show, ranges.show3)
CoveragePlot(gexatac_merge, region = "Hoxb3", extend.upstream = 30000, extend.downstream = 20000,
             features = c("Hoxb2","Hoxb3","Hoxb4", "Hoxb5"),
             region.highlight = ranges.all)
ggsave("Coverage_HoxB_Meis1and2_both_sites.pdf", width = 20, height = 8)




# HoxA
enriched.motifs <- FindMotifs(
  object = gexatac_merge,
  features = c(rownames(HoxA_peak))
)
# RNAでPharyngealで発現しているもの
df_avg_rna <- AverageExpression(object = gexatac_merge,  group.by = "region2")$RNA %>% as.data.frame()
df_avg_rna <- df_avg_rna %>% dplyr::filter(Pharyngeal > 0.3)
exp.genes <- toupper(rownames(df_avg_rna))
# enriched.motifsのmortif.nameが重複阻止されているので修正する
enriched.motifs <- enriched.motifs %>% mutate(motif.name_correct =　str_replace(motif.name, "\\..*$", ""))
enriched.motifs <- enriched.motifs %>% dplyr::filter(motif.name_correct %in% exp.genes)
MotifPlot(
  object = gexatac_merge,
  motifs = head(rownames(enriched.motifs))
)
ggsave("Motif_HoxA_region_pharyngeal_expressed.pdf", width = 8, height = 6)

# HoxB
enriched.motifs <- FindMotifs(
  object = gexatac_merge,
  features = c(rownames(HoxB_peak))
)
# RNAでPharyngealで発現しているもの
df_avg_rna <- AverageExpression(object = gexatac_merge,  group.by = "region2")$RNA %>% as.data.frame()
df_avg_rna <- df_avg_rna %>% dplyr::filter(Pharyngeal > 0.3)
exp.genes <- toupper(rownames(df_avg_rna))
# enriched.motifsのmortif.nameが重複阻止されているので修正する
enriched.motifs <- enriched.motifs %>% mutate(motif.name_correct =　str_replace(motif.name, "\\..*$", ""))
enriched.motifs <- enriched.motifs %>% dplyr::filter(motif.name_correct %in% exp.genes)
MotifPlot(
  object = gexatac_merge,
  motifs = head(rownames(enriched.motifs))
)
ggsave("Motif_HoxB_region_pharyngeal_expressed.pdf", width = 8, height = 6)





