
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
library(ggbeeswarm)

setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered")
NCC_merge <- readRDS("NCC_merge2_RPCA_velo_region.rds")
Idents(NCC_merge) <- "integrated"
DimPlot(NCC_merge, reduction = "umap.rpca")
DefaultAssay(NCC_merge)
dim(NCC_merge)


coi <- NCC_merge@meta.data %>% dplyr::filter(sample %in% c("E11.5_EYFP","E12.5_EYFP","E14.5_EYFP","E17.5_EYFP"))
gexatac_merge <- NCC_merge[, rownames(coi)]
DimPlot(gexatac_merge, reduction = "umap.rpca", label = T)

gexatac_merge

coi <- gexatac_merge@meta.data %>% dplyr::filter(!seurat_clusters %in% c("24","17","20","12","9","25","26","1","3","19"))
mes <- gexatac_merge[, rownames(coi)]
DimPlot(mes, reduction = "umap.rpca", label = T)

DimPlot(mes, reduction = "umap.rpca", group.by = "region",label = F, pt.size = 0.5) 
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/mes/region_specific")
ggsave("umap.mes_region.pdf", height = 4, width = 5)




library(tidyplots)
df_rna <- NCC_merge[["RNA"]]$data[c("Hoxa2","Hoxb2",
                              "Hoxa3","Hoxb3","Hoxd3",
                              "Hoxa4","Hoxb4", "Hoxc4","Hoxd4",
                              "Hoxa5","Hoxb5", "Hoxc5"
                              ),] %>% as.data.frame() %>% t() %>% as.data.frame

df_rna <- df_rna %>%
  mutate(Hox_status = case_when(
    (Hoxa2 > 0 | Hoxb2 > 0) & 
      (Hoxa3 == 0 & Hoxb3 == 0 & Hoxd3 == 0) &
      (Hoxa4 == 0 & Hoxb4 == 0 & Hoxc4 == 0 & Hoxd4 == 0) &
      (Hoxa5 == 0 & Hoxb5 == 0 & Hoxc5 == 0) ~ "Hox2",
    
    (Hoxa3 > 0 | Hoxb3 > 0 | Hoxd3 > 0) & 
      (Hoxa4 == 0 & Hoxb4 == 0 & Hoxc4 == 0 & Hoxd4 == 0) &
      (Hoxa5 == 0 & Hoxb5 == 0 & Hoxc5 == 0) ~ "Hox3",
    
    (Hoxa4 > 0 | Hoxb4 > 0 | Hoxc4 > 0 | Hoxd4 > 0) & 
      (Hoxa5 == 0 & Hoxb5 == 0 & Hoxc5 == 0) ~ "Hox4",
    
    (Hoxa5 > 0 | Hoxb5 > 0 | Hoxc5 > 0) ~ "Hox5",
    
    TRUE ~ "Hox_Neg"
  ))

df <- cbind(df_rna, NCC_merge@meta.data[c(71,73)])

df <- df %>% dplyr::filter(region2 %in%  c("Pharyngeal","Transitional","Cardiac_Cushion", "Subvalvular","AP septum", "SMC_GA","SMC_DA","SMC_CA"))
df$region2 <- as.factor(df$region2)
df$region2 <- factor(df$region2,
                    levels =  c("Pharyngeal","Transitional","Cardiac_Cushion", "Subvalvular","AP septum", "SMC_GA","SMC_DA","SMC_CA"))
df %>%
  ggplot(aes(x = region2, fill = Hox_status)) +
  geom_bar(position = "fill") +
  labs(y = "Proportion", x = "region", fill = "Hox_status") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())
ggsave("barplot_region_Hox.pdf", width = 4, height = 4)


df %>% head

df <- df %>%
  mutate(Hox2_mean = (Hoxa2 + Hoxb2) / 2,
         Hox3_mean = (Hoxa3 + Hoxb3 + Hoxd3) / 3,
         Hox4_mean = (Hoxa4 + Hoxb4 + Hoxc4 + Hoxd4) / 4,
         Hox5_mean = (Hoxa5 + Hoxb5 + Hoxc5) / 3,
         ) 

goi <- "Hox4_mean"
df %>% 
  tidyplot(x = region, y = !!sym(goi), color = region) %>% 
  add_mean_bar(alpha = 0.3) |>  
  add_sd_errorbar() |>  
  add_test_asterisks(ref.group = 1, hide_info = TRUE, method = "wilcox_test", p.adjust.method = "bonferroni") +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())  
ggsave(paste(goi, ".pdf", sep = ""), width = 4, height = 3, bg = "white")  

df_long <- df %>%
  group_by(region2) %>%
  summarise(Hox2_mean = mean(Hox2_mean, na.rm = TRUE),
            Hox3_mean = mean(Hox3_mean, na.rm = TRUE),
            Hox4_mean = mean(Hox4_mean, na.rm = TRUE),
            Hox5_mean = mean(Hox5_mean, na.rm = TRUE)) %>%
  pivot_longer(cols = c(Hox2_mean, Hox3_mean, Hox4_mean, Hox5_mean), 
               names_to = "Gene", 
               values_to = "Expression")
df_long <- df_long %>%
  mutate(Gene = factor(Gene, levels = c("Hox5_mean", "Hox4_mean", "Hox3_mean", "Hox2_mean")))

# heatmap
ggplot(df_long, aes(x = region2, y = Gene, fill = Expression)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  theme(axis.title.x  = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/plots/region_Hox")
ggsave("heatmap_Hox.pdf", width = 4, height = 4, bg = "white")
