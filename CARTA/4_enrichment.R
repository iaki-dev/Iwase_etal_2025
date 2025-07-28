
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
set.seed(1234)
library(data.table)
library(progress)

"%not.in%" <- Negate("%in%")


setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/v6")
CARTA_Net_summarize <- read_tsv("CARTA_Net_summarize_filt_conservation.txt")
dim(CARTA_Net_summarize)
# [1] 80869     26


enrichment_all <- CARTA_Net_summarize %>%
  group_by(region, name_ID) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(region) %>%
  mutate(
    n_target = n_distinct(CARTA_Net_summarize$target[CARTA_Net_summarize$region == unique(region)]),
    enrichment = count / n_target
  ) %>%
  ungroup() %>%
  arrange(region, desc(enrichment))
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/v6/enrichment/")
write_tsv(enrichment_all, "enrichment_all.txt")


df <- enrichment_all %>% dplyr::filter(name_ID %in% c("Meis2_MA1640.1","Meis2_MA0774.1"))

df_wide <- df %>% dplyr::select(region, name_ID, enrichment) %>% 
  pivot_wider(names_from = name_ID, values_from = c(enrichment))
df_ratio <- df_wide %>%
  mutate(
    enrichment_ratio = Meis2_MA0774.1 / Meis2_MA1640.1
  )
df_ratio %>% dplyr::filter(region %in% c("Cardiac_Cushion", "Pharyngeal")) %>% 
  ggplot(aes(x = region, y = enrichment_ratio)) +
  geom_col(aes(fill = region)) +
  labs(y = "Meis2 Motif Ratio (Hexamer / Octamer)") +
  theme_classic() +
  theme(
    legend.position = "bottom", 
    legend.text = element_text(size = 10),
    axis.text.x = element_blank() ,
    axis.title.x = element_blank()
  ) +
  scale_fill_manual(values = c(
    "Pharyngeal" = "#F8766D",       
    "Cardiac_Cushion" = "#8CAB00" 
  ))
ggsave("Meis2_Hexameric_Octamerc_ratio.pdf", width = 4, height = 6)

