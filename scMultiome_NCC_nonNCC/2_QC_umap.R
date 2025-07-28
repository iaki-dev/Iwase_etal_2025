

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
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")


setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/data")
E11.5_NCC <- readRDS("E11.5_EYFP_raw.v1.rds")
E12.5_NCC <- readRDS("E12.5_EYFP_raw.v1.rds")
E14.5_NCC <- readRDS("E14.5_EYFP_raw.v1.rds")
E17.5_NCC <- readRDS("E17.5_EYFP_raw.v1.rds")
E11.5_nonNCC <- readRDS("E11.5_non_EYFP_raw.v1.rds")
E12.5_nonNCC <- readRDS("E12.5_non_EYFP_raw.v1.rds")
E14.5_nonNCC <- readRDS("E14.5_non_EYFP_raw.v1.rds")
E17.5_nonNCC <- readRDS("E17.5_non_EYFP_raw.v1.rds")

# SoupXデータ  E12.5_nonNCC
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/v8/E12.5_nonNCC")
rna_counts <- Read10X(data.dir = "strainedCounts")
E12.5_nonNCC_soup <- CreateSeuratObject(counts = rna_counts, project = "E12.5_nonNCC")
E12.5_nonNCC_soup[["percent.mt"]] <- PercentageFeatureSet(E12.5_nonNCC_soup, pattern = "^mt-")
dim(E12.5_nonNCC_soup)
dim(E12.5_nonNCC[["RNA"]])

# SoupXデータ  E11.5_nonNCC
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/v8/E11.5_nonNCC")
rna_counts <- Read10X(data.dir = "strainedCounts")
E11.5_nonNCC_soup <- CreateSeuratObject(counts = rna_counts, project = "E11.5_nonNCC")
E11.5_nonNCC_soup[["percent.mt"]] <- PercentageFeatureSet(E11.5_nonNCC_soup, pattern = "^mt-")
dim(E11.5_nonNCC_soup)
dim(E11.5_nonNCC[["RNA"]])


E11.5_nonNCC[["RNA"]] <-E11.5_nonNCC_soup[["RNA"]]
E12.5_nonNCC[["RNA"]] <-E12.5_nonNCC_soup[["RNA"]]


setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/v10")

###############
#  Doublet finder
###############
DefaultAssay(E11.5_NCC) <- "RNA"
DefaultAssay(E12.5_NCC) <- "RNA"
DefaultAssay(E14.5_NCC) <- "RNA"
DefaultAssay(E17.5_NCC) <- "RNA"
DefaultAssay(E11.5_nonNCC) <- "RNA"
DefaultAssay(E12.5_nonNCC) <- "RNA"
DefaultAssay(E14.5_nonNCC) <- "RNA"
DefaultAssay(E17.5_nonNCC) <- "RNA"

# doubletfinder
# remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", upgrade = F)
suppressMessages(require(DoubletFinder))
# E11.5_NCC
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_E11.5_NCC <- paramSweep_v3(E11.5_NCC, PCs = 1:20, sct = T)
sweep.stats_E11.5_NCC <- summarizeSweep(sweep.res.list_E11.5_NCC, GT = FALSE)
#  ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results
bcmvn_E11.5_NCC <- find.pK(sweep.stats_E11.5_NCC)
# barplot(bcmvn_E11.5_NCC$BCmetric, names.arg = bcmvn_E11.5_NCC$pK, las=2)
ggplot(bcmvn_E11.5_NCC, aes(x=pK, y=BCmetric)) +
  geom_point()
# 3x10で保存
ggsave("bcmvn_E11.5_NCC.pdf", height = 3, width = 10)

dim(E11.5_NCC)
# [1] 32286    876
annotations <- E11.5_NCC@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.007*nrow(E11.5_NCC@meta.data))  ## Assuming 0.7% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# pk automation
pk <- bcmvn_E11.5_NCC %>% dplyr::filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
pk <- as.numeric(as.character(pk[[1]]))


E11.5_NCC <- doubletFinder_v3(E11.5_NCC, PCs = 1:20, pN = 0.25, pK = pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
DF.name = colnames(E11.5_NCC@meta.data)[grepl("DF.classification", colnames(E11.5_NCC@meta.data))]
DF.name

lastnum <- E11.5_NCC@meta.data %>% ncol()
E11.5_NCC@meta.data$DF <- E11.5_NCC@meta.data[,lastnum]
E11.5_NCC@meta.data %>% head


# E11.5_nonNCC
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_E11.5_nonNCC <- paramSweep_v3(E11.5_nonNCC, PCs = 1:20, sct = T)
sweep.stats_E11.5_nonNCC <- summarizeSweep(sweep.res.list_E11.5_nonNCC, GT = FALSE)
#  ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results
bcmvn_E11.5_nonNCC <- find.pK(sweep.stats_E11.5_nonNCC)
# barplot(bcmvn_E11.5_nonNCC$BCmetric, names.arg = bcmvn_E11.5_nonNCC$pK, las=2)
ggplot(bcmvn_E11.5_nonNCC, aes(x=pK, y=BCmetric)) +
  geom_point()
# 3x10で保存
ggsave("bcmvn_E11.5_nonNCC.pdf", height = 3, width = 10)

dim(E11.5_nonNCC)
# [1] 32286    924
annotations <- E11.5_nonNCC@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.008*nrow(E11.5_nonNCC@meta.data))  ## Assuming 0.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# pk automation
pk <- bcmvn_E11.5_nonNCC %>% dplyr::filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
pk <- as.numeric(as.character(pk[[1]]))


E11.5_nonNCC <- doubletFinder_v3(E11.5_nonNCC, PCs = 1:20, pN = 0.25, pK = pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
DF.name = colnames(E11.5_nonNCC@meta.data)[grepl("DF.classification", colnames(E11.5_nonNCC@meta.data))]
DF.name

lastnum <- E11.5_nonNCC@meta.data %>% ncol()
E11.5_nonNCC@meta.data$DF <- E11.5_nonNCC@meta.data[,lastnum]
E11.5_nonNCC@meta.data %>% head



# E12.5_NCC
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_E12.5_NCC <- paramSweep_v3(E12.5_NCC, PCs = 1:20, sct = T)
sweep.stats_E12.5_NCC <- summarizeSweep(sweep.res.list_E12.5_NCC, GT = FALSE)
#  ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results
bcmvn_E12.5_NCC <- find.pK(sweep.stats_E12.5_NCC)
# barplot(bcmvn_E12.5_NCC$BCmetric, names.arg = bcmvn_E12.5_NCC$pK, las=2)
ggplot(bcmvn_E12.5_NCC, aes(x=pK, y=BCmetric)) +
  geom_point()
# 3x10で保存
ggsave("bcmvn_E12.5_NCC.pdf", height = 3, width = 10)

dim(E12.5_NCC)
# [1] 32286 12648
annotations <- E12.5_NCC@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.1*nrow(E12.5_NCC@meta.data))  ## Assuming 10% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# pk automation
pk <- bcmvn_E12.5_NCC %>% dplyr::filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
pk <- as.numeric(as.character(pk[[1]]))


E12.5_NCC <- doubletFinder_v3(E12.5_NCC, PCs = 1:20, pN = 0.25, pK = pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
DF.name = colnames(E12.5_NCC@meta.data)[grepl("DF.classification", colnames(E12.5_NCC@meta.data))]
DF.name

lastnum <- E12.5_NCC@meta.data %>% ncol()
E12.5_NCC@meta.data$DF <- E12.5_NCC@meta.data[,lastnum]
E12.5_NCC@meta.data %>% head


# E12.5_nonNCC
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_E12.5_nonNCC <- paramSweep_v3(E12.5_nonNCC, PCs = 1:20, sct = T)
sweep.stats_E12.5_nonNCC <- summarizeSweep(sweep.res.list_E12.5_nonNCC, GT = FALSE)
#  ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results
bcmvn_E12.5_nonNCC <- find.pK(sweep.stats_E12.5_nonNCC)
# barplot(bcmvn_E12.5_nonNCC$BCmetric, names.arg = bcmvn_E12.5_nonNCC$pK, las=2)
ggplot(bcmvn_E12.5_nonNCC, aes(x=pK, y=BCmetric)) +
  geom_point()
# 3x10で保存
ggsave("bcmvn_E12.5_nonNCC.pdf", height = 3, width = 10)

dim(E12.5_nonNCC)
# [1]  32286 3289 
annotations <- E12.5_nonNCC@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.025*nrow(E12.5_nonNCC@meta.data))  ## Assuming 2.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# pk automation
pk <- bcmvn_E12.5_nonNCC %>% dplyr::filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
pk <- as.numeric(as.character(pk[[1]]))


E12.5_nonNCC <- doubletFinder_v3(E12.5_nonNCC, PCs = 1:20, pN = 0.25, pK = pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
DF.name = colnames(E12.5_nonNCC@meta.data)[grepl("DF.classification", colnames(E12.5_nonNCC@meta.data))]
DF.name

lastnum <- E12.5_nonNCC@meta.data %>% ncol()
E12.5_nonNCC@meta.data$DF <- E12.5_nonNCC@meta.data[,lastnum]
E12.5_nonNCC@meta.data %>% head


# E14.5_NCC
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_E14.5_NCC <- paramSweep_v3(E14.5_NCC, PCs = 1:20, sct = T)
sweep.stats_E14.5_NCC <- summarizeSweep(sweep.res.list_E14.5_NCC, GT = FALSE)
#  ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results
bcmvn_E14.5_NCC <- find.pK(sweep.stats_E14.5_NCC)
# barplot(bcmvn_E14.5_NCC$BCmetric, names.arg = bcmvn_E14.5_NCC$pK, las=2)
ggplot(bcmvn_E14.5_NCC, aes(x=pK, y=BCmetric)) +
  geom_point()
# 3x10で保存
ggsave("bcmvn_E14.5_NCC.pdf", height = 3, width = 10)

dim(E14.5_NCC)
# [1] 32286    593
annotations <- E14.5_NCC@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.005*nrow(E14.5_NCC@meta.data))  ## Assuming 0.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# pk automation
pk <- bcmvn_E14.5_NCC %>% dplyr::filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
pk <- as.numeric(as.character(pk[[1]]))


E14.5_NCC <- doubletFinder_v3(E14.5_NCC, PCs = 1:20, pN = 0.25, pK = pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
DF.name = colnames(E14.5_NCC@meta.data)[grepl("DF.classification", colnames(E14.5_NCC@meta.data))]
DF.name

lastnum <- E14.5_NCC@meta.data %>% ncol()
E14.5_NCC@meta.data$DF <- E14.5_NCC@meta.data[,lastnum]
E14.5_NCC@meta.data %>% head


# E14.5_nonNCC
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_E14.5_nonNCC <- paramSweep_v3(E14.5_nonNCC, PCs = 1:20, sct = T)
sweep.stats_E14.5_nonNCC <- summarizeSweep(sweep.res.list_E14.5_nonNCC, GT = FALSE)
#  ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results
bcmvn_E14.5_nonNCC <- find.pK(sweep.stats_E14.5_nonNCC)
# barplot(bcmvn_E14.5_nonNCC$BCmetric, names.arg = bcmvn_E14.5_nonNCC$pK, las=2)
ggplot(bcmvn_E14.5_nonNCC, aes(x=pK, y=BCmetric)) +
  geom_point()
# 3x10で保存
ggsave("bcmvn_E14.5_nonNCC.pdf", height = 3, width = 10)

dim(E14.5_nonNCC)
# [1] 32286    1145
annotations <- E14.5_nonNCC@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.012*nrow(E14.5_nonNCC@meta.data))  ## Assuming 1.2% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# pk automation
pk <- bcmvn_E14.5_nonNCC %>% dplyr::filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
pk <- as.numeric(as.character(pk[[1]]))


E14.5_nonNCC <- doubletFinder_v3(E14.5_nonNCC, PCs = 1:20, pN = 0.25, pK = pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
DF.name = colnames(E14.5_nonNCC@meta.data)[grepl("DF.classification", colnames(E14.5_nonNCC@meta.data))]
DF.name

lastnum <- E14.5_nonNCC@meta.data %>% ncol()
E14.5_nonNCC@meta.data$DF <- E14.5_nonNCC@meta.data[,lastnum]
E14.5_nonNCC@meta.data %>% head


# E17.5_NCC
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_E17.5_NCC <- paramSweep_v3(E17.5_NCC, PCs = 1:20, sct = T)
sweep.stats_E17.5_NCC <- summarizeSweep(sweep.res.list_E17.5_NCC, GT = FALSE)
#  ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results
bcmvn_E17.5_NCC <- find.pK(sweep.stats_E17.5_NCC)
# barplot(bcmvn_E17.5_NCC$BCmetric, names.arg = bcmvn_E17.5_NCC$pK, las=2)
ggplot(bcmvn_E17.5_NCC, aes(x=pK, y=BCmetric)) +
  geom_point()
# 3x10で保存
ggsave("bcmvn_E17.5_NCC.pdf", height = 3, width = 10)

dim(E17.5_NCC)
# [1] 32286    344
annotations <- E17.5_NCC@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.004*nrow(E17.5_NCC@meta.data))  ## Assuming 0.4% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# pk automation
pk <- bcmvn_E17.5_NCC %>% dplyr::filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
pk <- as.numeric(as.character(pk[[1]]))


E17.5_NCC <- doubletFinder_v3(E17.5_NCC, PCs = 1:20, pN = 0.25, pK = pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
DF.name = colnames(E17.5_NCC@meta.data)[grepl("DF.classification", colnames(E17.5_NCC@meta.data))]
DF.name

lastnum <- E17.5_NCC@meta.data %>% ncol()
E17.5_NCC@meta.data$DF <- E17.5_NCC@meta.data[,lastnum]
E17.5_NCC@meta.data %>% head


# E17.5_nonNCC
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_E17.5_nonNCC <- paramSweep_v3(E17.5_nonNCC, PCs = 1:20, sct = T)
sweep.stats_E17.5_nonNCC <- summarizeSweep(sweep.res.list_E17.5_nonNCC, GT = FALSE)
#  ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results
bcmvn_E17.5_nonNCC <- find.pK(sweep.stats_E17.5_nonNCC)
# barplot(bcmvn_E17.5_nonNCC$BCmetric, names.arg = bcmvn_E17.5_nonNCC$pK, las=2)
ggplot(bcmvn_E17.5_nonNCC, aes(x=pK, y=BCmetric)) +
  geom_point()
# 3x10で保存
ggsave("bcmvn_E17.5_nonNCC.pdf", height = 3, width = 10)

dim(E17.5_nonNCC)
# [1] 32286   1771
annotations <- E17.5_nonNCC@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.014*nrow(E17.5_nonNCC@meta.data))  ## Assuming 1.4% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# pk automation
pk <- bcmvn_E17.5_nonNCC %>% dplyr::filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
pk <- as.numeric(as.character(pk[[1]]))


E17.5_nonNCC <- doubletFinder_v3(E17.5_nonNCC, PCs = 1:20, pN = 0.25, pK = pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
DF.name = colnames(E17.5_nonNCC@meta.data)[grepl("DF.classification", colnames(E17.5_nonNCC@meta.data))]
DF.name

lastnum <- E17.5_nonNCC@meta.data %>% ncol()
E17.5_nonNCC@meta.data$DF <- E17.5_nonNCC@meta.data[,lastnum]
E17.5_nonNCC@meta.data %>% head

# 保存
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/v10/data")
# saveRDS(E11.5_NCC, "E11.5_NCC_raw.v1.rds")
# saveRDS(E11.5_nonNCC, "E11.5_nonNCC_raw.v1.rds")
# saveRDS(E12.5_NCC, "E12.5_NCC_raw.v1.rds")
# saveRDS(E12.5_nonNCC, "E12.5_nonNCC_raw.v1.rds")
# saveRDS(E17.5_NCC, "E17.5_NCC_raw.v1.rds")
# saveRDS(E14.5_NCC, "E14.5_NCC_raw.v1.rds")
# saveRDS(E17.5_nonNCC, "E17.5_nonNCC_raw.v1.rds")
# saveRDS(E14.5_nonNCC, "E14.5_nonNCC_raw.v1.rds")
# ロード
E11.5_NCC <- readRDS("E11.5_NCC_raw.v1.rds")
E12.5_NCC <- readRDS("E12.5_NCC_raw.v1.rds")
E14.5_NCC <- readRDS("E14.5_NCC_raw.v1.rds")
E17.5_NCC <- readRDS("E17.5_NCC_raw.v1.rds")
E11.5_nonNCC <- readRDS("E11.5_nonNCC_raw.v1.rds")
E12.5_nonNCC <- readRDS("E12.5_nonNCC_raw.v1.rds")
E14.5_nonNCC <- readRDS("E14.5_nonNCC_raw.v1.rds")
E17.5_nonNCC <- readRDS("E17.5_nonNCC_raw.v1.rds")


# Convert for Assay5
E11.5_NCC[["RNA"]] <- as(object = E11.5_NCC[["RNA"]], Class = "Assay5")
E12.5_NCC[["RNA"]] <- as(object = E12.5_NCC[["RNA"]], Class = "Assay5")
E14.5_NCC[["RNA"]] <- as(object = E14.5_NCC[["RNA"]], Class = "Assay5")
E17.5_NCC[["RNA"]] <- as(object = E17.5_NCC[["RNA"]], Class = "Assay5")
E11.5_nonNCC[["RNA"]] <- as(object = E11.5_nonNCC[["RNA"]], Class = "Assay5")
E12.5_nonNCC[["RNA"]] <- as(object = E12.5_nonNCC[["RNA"]], Class = "Assay5")
E14.5_nonNCC[["RNA"]] <- as(object = E14.5_nonNCC[["RNA"]], Class = "Assay5")
E17.5_nonNCC[["RNA"]] <- as(object = E17.5_nonNCC[["RNA"]], Class = "Assay5")


E11.5_NCC <- subset(
  x = E11.5_NCC,   subset =
    nCount_RNA > 3000 & nCount_RNA < 5e4&
    nCount_ATAC > 3000 & nCount_ATAC < 5e5 &
    TSS.enrichment > 2 & nucleosome_signal < 4 &
    percent.mt < 25
)
E12.5_NCC <- subset(
  x = E12.5_NCC,   subset =
    nCount_RNA > 3000 & nCount_RNA < 5e4&
    nCount_ATAC > 3000 & nCount_ATAC < 5e5 &
    TSS.enrichment > 2 & nucleosome_signal < 4 &
    percent.mt < 25
)
E14.5_NCC <- subset(
  x = E14.5_NCC,   subset =
    nCount_RNA > 3000 & nCount_RNA < 5e4&
    nCount_ATAC > 3000 & nCount_ATAC < 5e5 &
    TSS.enrichment > 2 & nucleosome_signal < 4 &
    percent.mt < 25
)
E17.5_NCC <- subset(
  x = E17.5_NCC,   subset =
    nCount_RNA > 3000 & nCount_RNA < 5e4&
    nCount_ATAC > 3000 & nCount_ATAC < 5e5 &
    TSS.enrichment > 2 & nucleosome_signal < 4 &
    percent.mt < 25
)
E11.5_nonNCC <- subset(
  x = E11.5_nonNCC,   subset =
    nCount_RNA > 3000 & nCount_RNA < 5e4&
    nCount_ATAC > 3000 & nCount_ATAC < 5e5 &
    TSS.enrichment > 2 & nucleosome_signal < 4 &
    percent.mt < 25
)
E12.5_nonNCC <- subset(
  x = E12.5_nonNCC,   subset =
    nCount_RNA > 3000 & nCount_RNA < 5e4&
    nCount_ATAC > 3000 & nCount_ATAC < 5e5 &
    TSS.enrichment > 2 & nucleosome_signal < 4 &
    percent.mt < 25
)
E14.5_nonNCC <- subset(
  x = E14.5_nonNCC,   subset =
    nCount_RNA > 3000 & nCount_RNA < 5e4&
    nCount_ATAC > 3000 & nCount_ATAC < 5e5 &
    TSS.enrichment > 2 & nucleosome_signal < 4 &
    percent.mt < 25
)
E17.5_nonNCC <- subset(
  x = E17.5_nonNCC,   subset =
    nCount_RNA > 3000 & nCount_RNA < 5e4&
    nCount_ATAC > 3000 & nCount_ATAC < 5e5 &
    TSS.enrichment > 2 & nucleosome_signal < 4 &
    percent.mt < 25
)

# Extract Singlet
coi <- E11.5_NCC@meta.data %>% dplyr::filter(DF == "Singlet")
E11.5_NCC <- E11.5_NCC[, rownames(coi)]
coi <- E12.5_NCC@meta.data %>% dplyr::filter(DF == "Singlet")
E12.5_NCC <- E12.5_NCC[, rownames(coi)]
coi <- E14.5_NCC@meta.data %>% dplyr::filter(DF == "Singlet")
E14.5_NCC <- E14.5_NCC[, rownames(coi)]
coi <- E17.5_NCC@meta.data %>% dplyr::filter(DF == "Singlet")
E17.5_NCC <- E17.5_NCC[, rownames(coi)]
coi <- E11.5_nonNCC@meta.data %>% dplyr::filter(DF == "Singlet")
E11.5_nonNCC <- E11.5_nonNCC[, rownames(coi)]
coi <- E12.5_nonNCC@meta.data %>% dplyr::filter(DF == "Singlet")
E12.5_nonNCC <- E12.5_nonNCC[, rownames(coi)]
coi <- E14.5_nonNCC@meta.data %>% dplyr::filter(DF == "Singlet")
E14.5_nonNCC <- E14.5_nonNCC[, rownames(coi)]
coi <- E17.5_nonNCC@meta.data %>% dplyr::filter(DF == "Singlet")
E17.5_nonNCC <- E17.5_nonNCC[, rownames(coi)]


# ribosomal_genes ratio
E11.5_NCC <- PercentageFeatureSet(E11.5_NCC, "^Rp[sl]", col.name = "percent_ribo")
E12.5_NCC <- PercentageFeatureSet(E12.5_NCC, "^Rp[sl]", col.name = "percent_ribo")
E14.5_NCC <- PercentageFeatureSet(E14.5_NCC, "^Rp[sl]", col.name = "percent_ribo")
E17.5_NCC <- PercentageFeatureSet(E17.5_NCC, "^Rp[sl]", col.name = "percent_ribo")
E11.5_nonNCC <- PercentageFeatureSet(E11.5_nonNCC, "^Rp[sl]", col.name = "percent_ribo")
E12.5_nonNCC <- PercentageFeatureSet(E12.5_nonNCC, "^Rp[sl]", col.name = "percent_ribo")
E14.5_nonNCC <- PercentageFeatureSet(E14.5_nonNCC, "^Rp[sl]", col.name = "percent_ribo")
E17.5_nonNCC <- PercentageFeatureSet(E17.5_nonNCC, "^Rp[sl]", col.name = "percent_ribo")



tmp <- E11.5_NCC[["RNA"]]$counts["EYFP",]
tmp <- tmp[tmp > 0] %>% as.data.frame()
E11.5_NCC <- E11.5_NCC[, rownames(tmp)]

tmp <- E12.5_NCC[["RNA"]]$counts["EYFP",]
tmp <- tmp[tmp > 0] %>% as.data.frame()
E12.5_NCC <- E12.5_NCC[, rownames(tmp)]

tmp <- E14.5_NCC[["RNA"]]$counts["EYFP",]
tmp <- tmp[tmp > 0] %>% as.data.frame()
E14.5_NCC <- E14.5_NCC[, rownames(tmp)]

tmp <- E17.5_NCC[["RNA"]]$counts["EYFP",]
tmp <- tmp[tmp > 0] %>% as.data.frame()
E17.5_NCC <- E17.5_NCC[, rownames(tmp)]

tmp <- E11.5_nonNCC[["RNA"]]$counts["EYFP",]
tmp <- tmp[tmp == 0] %>% as.data.frame()
E11.5_nonNCC <- E11.5_nonNCC[, rownames(tmp)]

tmp <- E12.5_nonNCC[["RNA"]]$counts["EYFP",]
tmp <- tmp[tmp == 0] %>% as.data.frame()
E12.5_nonNCC <- E12.5_nonNCC[, rownames(tmp)]

tmp <- E14.5_nonNCC[["RNA"]]$counts["EYFP",]
tmp <- tmp[tmp == 0] %>% as.data.frame()
E14.5_nonNCC <- E14.5_nonNCC[, rownames(tmp)]

tmp <- E17.5_nonNCC[["RNA"]]$counts["EYFP",]
tmp <- tmp[tmp ==  0] %>% as.data.frame()
E17.5_nonNCC <- E17.5_nonNCC[, rownames(tmp)]



s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

s.genes <- paste(s.genes %>% str_sub(start = 1, end = 1), s.genes %>% str_sub(start = 2) %>% str_to_lower(), sep = "")
g2m.genes <- paste(g2m.genes %>% str_sub(start = 1, end = 1), g2m.genes %>% str_sub(start = 2) %>% str_to_lower(), sep = "")

E11.5_NCC <- NormalizeData(E11.5_NCC)
E11.5_NCC <- FindVariableFeatures(E11.5_NCC, selection.method = "vst")
E11.5_NCC <- ScaleData(E11.5_NCC, features = rownames(E11.5_NCC))
E11.5_NCC <- RunPCA(E11.5_NCC, features = VariableFeatures(E11.5_NCC, layer = "counts"))
E11.5_NCC <- CellCycleScoring(E11.5_NCC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
E11.5_NCC <- ScaleData(E11.5_NCC, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(E11.5_NCC))
E11.5_NCC <- RunPCA(E11.5_NCC, features = VariableFeatures(E11.5_NCC, layer = "counts"))

E12.5_NCC <- NormalizeData(E12.5_NCC)
E12.5_NCC <- FindVariableFeatures(E12.5_NCC, selection.method = "vst")
E12.5_NCC <- ScaleData(E12.5_NCC, features = rownames(E12.5_NCC))
E12.5_NCC <- RunPCA(E12.5_NCC, features = VariableFeatures(E12.5_NCC, layer = "counts"))
E12.5_NCC <- CellCycleScoring(E12.5_NCC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
E12.5_NCC <- ScaleData(E12.5_NCC, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(E12.5_NCC))
E12.5_NCC <- RunPCA(E12.5_NCC, features = VariableFeatures(E12.5_NCC, layer = "counts"))

E14.5_NCC <- NormalizeData(E14.5_NCC)
E14.5_NCC <- FindVariableFeatures(E14.5_NCC, selection.method = "vst")
E14.5_NCC <- ScaleData(E14.5_NCC, features = rownames(E14.5_NCC))
E14.5_NCC <- RunPCA(E14.5_NCC, features = VariableFeatures(E14.5_NCC, layer = "counts"))
E14.5_NCC <- CellCycleScoring(E14.5_NCC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
E14.5_NCC <- ScaleData(E14.5_NCC, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(E14.5_NCC))
E14.5_NCC <- RunPCA(E14.5_NCC, features = VariableFeatures(E14.5_NCC, layer = "counts"))

E17.5_NCC <- NormalizeData(E17.5_NCC)
E17.5_NCC <- FindVariableFeatures(E17.5_NCC, selection.method = "vst")
E17.5_NCC <- ScaleData(E17.5_NCC, features = rownames(E17.5_NCC))
E17.5_NCC <- RunPCA(E17.5_NCC, features = VariableFeatures(E17.5_NCC, layer = "counts"))
E17.5_NCC <- CellCycleScoring(E17.5_NCC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
E17.5_NCC <- ScaleData(E17.5_NCC, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(E17.5_NCC))
E17.5_NCC <- RunPCA(E17.5_NCC, features = VariableFeatures(E17.5_NCC, layer = "counts"))

E11.5_nonNCC <- NormalizeData(E11.5_nonNCC)
E11.5_nonNCC <- FindVariableFeatures(E11.5_nonNCC, selection.method = "vst")
E11.5_nonNCC <- ScaleData(E11.5_nonNCC, features = rownames(E11.5_nonNCC))
E11.5_nonNCC <- RunPCA(E11.5_nonNCC, features = VariableFeatures(E11.5_nonNCC, layer = "counts"))
E11.5_nonNCC <- CellCycleScoring(E11.5_nonNCC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
E11.5_nonNCC <- ScaleData(E11.5_nonNCC, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(E11.5_nonNCC))
E11.5_nonNCC <- RunPCA(E11.5_nonNCC, features = VariableFeatures(E11.5_nonNCC, layer = "counts"))

E12.5_nonNCC <- NormalizeData(E12.5_nonNCC)
E12.5_nonNCC <- FindVariableFeatures(E12.5_nonNCC, selection.method = "vst")
E12.5_nonNCC <- ScaleData(E12.5_nonNCC, features = rownames(E12.5_nonNCC))
E12.5_nonNCC <- RunPCA(E12.5_nonNCC, features = VariableFeatures(E12.5_nonNCC, layer = "counts"))
E12.5_nonNCC <- CellCycleScoring(E12.5_nonNCC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
E12.5_nonNCC <- ScaleData(E12.5_nonNCC, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(E12.5_nonNCC))
E12.5_nonNCC <- RunPCA(E12.5_nonNCC, features = VariableFeatures(E12.5_nonNCC, layer = "counts"))

E14.5_nonNCC <- NormalizeData(E14.5_nonNCC)
E14.5_nonNCC <- FindVariableFeatures(E14.5_nonNCC, selection.method = "vst")
E14.5_nonNCC <- ScaleData(E14.5_nonNCC, features = rownames(E14.5_nonNCC))
E14.5_nonNCC <- RunPCA(E14.5_nonNCC, features = VariableFeatures(E14.5_nonNCC, layer = "counts"))
E14.5_nonNCC <- CellCycleScoring(E14.5_nonNCC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
E14.5_nonNCC <- ScaleData(E14.5_nonNCC, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(E14.5_nonNCC))
E14.5_nonNCC <- RunPCA(E14.5_nonNCC, features = VariableFeatures(E14.5_nonNCC, layer = "counts"))

E17.5_nonNCC <- NormalizeData(E17.5_nonNCC)
E17.5_nonNCC <- FindVariableFeatures(E17.5_nonNCC, selection.method = "vst")
E17.5_nonNCC <- ScaleData(E17.5_nonNCC, features = rownames(E17.5_nonNCC))
E17.5_nonNCC <- RunPCA(E17.5_nonNCC, features = VariableFeatures(E17.5_nonNCC, layer = "counts"))
E17.5_nonNCC <- CellCycleScoring(E17.5_nonNCC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
E17.5_nonNCC <- ScaleData(E17.5_nonNCC, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(E17.5_nonNCC))
E17.5_nonNCC <- RunPCA(E17.5_nonNCC, features = VariableFeatures(E17.5_nonNCC, layer = "counts"))



setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/v10/v10.2/data")
# saveRDS(E11.5_NCC, "E11.5_NCC_regress_cc.v1.rds")
# saveRDS(E11.5_nonNCC, "E11.5_nonNCC_regress_cc.v1.rds")
# saveRDS(E12.5_NCC, "E12.5_NCC_regress_cc.v1.rds")
# saveRDS(E12.5_nonNCC, "E12.5_nonNCC_regress_cc.v1.rds")
# saveRDS(E17.5_NCC, "E17.5_NCC_regress_cc.v1.rds")
# saveRDS(E14.5_NCC, "E14.5_NCC_regress_cc.v1.rds")
# saveRDS(E17.5_nonNCC, "E17.5_nonNCC_regress_cc.v1.rds")
# saveRDS(E14.5_nonNCC, "E14.5_nonNCC_regress_cc.v1.rds")


gexatac_merge <- merge(E11.5_NCC, y = c(E12.5_NCC, E14.5_NCC, E17.5_NCC,
                                        E11.5_nonNCC ,E12.5_nonNCC, E14.5_nonNCC, E17.5_nonNCC),
                       add.cell.ids = c("E11.5_NCC", "E12.5_NCC", "E14.5_NCC", "E17.5_NCC",
                                        "E11.5_nonNCC", "E12.5_nonNCC", "E14.5_nonNCC", "E17.5_nonNCC"), project = "gexatac_merge")
sample <- gexatac_merge@meta.data$orig.ident %>%
  str_replace_all(pattern = "non_EYFP", replacement = "nonNCC") %>%
  str_replace_all(pattern = "EYFP", replacement = "NCC")
gexatac_merge <- AddMetaData(
  object = gexatac_merge,
  metadata = sample,
  col.name = 'sample'
)
gexatac_merge@meta.data$sample %>% table

gexatac_merge[["RNA"]] <- JoinLayers(gexatac_merge[["RNA"]])
gexatac_merge[["RNA"]] <- split(gexatac_merge[["RNA"]], f = gexatac_merge$sample)
gexatac_merge

gexatac_merge@meta.data$sample <- as.factor(gexatac_merge@meta.data$sample)
gexatac_merge@meta.data$sample <- factor(gexatac_merge@meta.data$sample,
                                         levels = c("E11.5_NCC", "E12.5_NCC", "E14.5_NCC", "E17.5_NCC",
                                                    "E11.5_nonNCC", "E12.5_nonNCC", "E14.5_nonNCC", "E17.5_nonNCC"))
# unnintegrated
gexatac_merge <- NormalizeData(gexatac_merge)
gexatac_merge <- FindVariableFeatures(gexatac_merge)
intersect(VariableFeatures(gexatac_merge), "EYFP")

gexatac_merge <- ScaleData(gexatac_merge, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(gexatac_merge))
gexatac_merge <- RunPCA(gexatac_merge)
gexatac_merge <- FindNeighbors(gexatac_merge, dims = 1:30, reduction = "pca")
gexatac_merge <- FindClusters(gexatac_merge, resolution = 2, cluster.name = "unintegrated_clusters")
gexatac_merge <- RunUMAP(gexatac_merge, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(gexatac_merge, reduction = "umap.unintegrated", group.by = c("sample"))
ggsave("1_umap_unintegrated_sample.pdf", height = 6, width = 6)

setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/v10/v10.2/data")
gexatac_merge <- readRDS("gexatac_merge.unintegrated.rds")

gexatac_merge <- IntegrateLayers(
  object = gexatac_merge, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca", dims = 1:50,
  verbose = FALSE
)


setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/v10/v10.2.1")

# RPCA
gexatac_merge <- RunUMAP(gexatac_merge, reduction = "integrated.rpca", dims = 1:50, reduction.name = "umap.rpca")
DimPlot(gexatac_merge, reduction = "umap.rpca", group.by = c("sample"), combine = FALSE)
ggsave("1_umap_RPCA_sample.pdf", height = 6, width = 6)
gsave("1_umap_RPCA_sample_resize.pdf", height = 4, width = 6)
DimPlot(gexatac_merge, reduction = "umap.rpca", group.by = c("sample"), combine = FALSE, split.by = "sample", ncol=4)
ggsave("1_umap_RPCA_sample_facet.pdf", height = 8, width = 12)

gexatac_merge <- FindNeighbors(gexatac_merge, reduction = "integrated.rpca", dims = 1:50, k.param = 25)
gexatac_merge <- FindClusters(gexatac_merge, resolution = 0.5, cluster.name = "rpca_clusters", algorithm = 1)
DimPlot(gexatac_merge, reduction = "umap.rpca", group.by = c("rpca_clusters"), combine = F, label = T)
gexatac_merge@meta.data %>% dplyr::select(seurat_clusters, sample) %>% table
ggsave("1_umap_RPCA_clusters.pdf", height = 6, width = 6)

gexatac_merge <- FindSubCluster(gexatac_merge, cluster = 9, graph.name = "RNA_snn", algorithm = 1, resolution = 0.2)
DimPlot(gexatac_merge, reduction = "umap.rpca", group.by = c("sub.cluster"), combine = F, label = T)
Idents(gexatac_merge) <- "sub.cluster"
gexatac_merge <- FindSubCluster(gexatac_merge, cluster = 4, graph.name = "RNA_snn", algorithm = 1, resolution = 0.3)
DimPlot(gexatac_merge, reduction = "umap.rpca", group.by = c("sub.cluster"), combine = F, label = T)
Idents(gexatac_merge) <- "sub.cluster"
gexatac_merge@meta.data$sub.cluster <- gexatac_merge@meta.data$sub.cluster %>%
  str_replace_all(pattern = "9_0", replacement = "9") %>%
  str_replace_all(pattern = "9_1", replacement = "19") %>% 
  str_replace_all(pattern = "4_1", replacement = "4") %>%
  str_replace_all(pattern = "4_0", replacement = "20") 
  
gexatac_merge@meta.data$sub.cluster <- as.factor(gexatac_merge@meta.data$sub.cluster)
gexatac_merge@meta.data$sub.cluster <- factor(gexatac_merge@meta.data$sub.cluster,
                                              levels = 0:20)
DimPlot(gexatac_merge, reduction = "umap.rpca", group.by = c("sub.cluster"), combine = F, label = T)
ggsave("1_umap_RPCA_subclusters.pdf", height = 6, width = 6)


FeaturePlot(gexatac_merge, features = c("Sox9", "Scx", "Twist1", "Foxp2", "Sox5", "Dlx5","Col2a1", "Runx2", "Lef1", "Pitx2", "Tbx4", "Tcf24", "Tcf21",
                                               "Acta2","Myh11","Dlk1", "Abcc9","Osr1",
                                               "Tnnt2", "Myh3", "Myod1",
                                               "Nr2f2","Cdh5", "Jag1", "Cdh11", "Aplnr","Prox1","Csf1r", "Runx1", "Hbb-y", "Wt1", "Sox10", "Tubb3",
                                               "Mki67",
                                               "Epcam" ,"EYFP"), reduction = "umap.rpca", ncol = 6, pt.size = 0.01)
ggsave("Featureplot_markers.pdf", height = 16, width = 16)


gexatac_merge <- JoinLayers(gexatac_merge)
gexatac_merge

gexatac_merge$sub.cluster <- gexatac_merge$sub.cluster %>% as.character()
order_factor <- gexatac_merge$sub.cluster %>% as.numeric() %>% unique() %>% sort()
gexatac_merge$sub.cluster <- as.factor(gexatac_merge$sub.cluster)
gexatac_merge$sub.cluster <- factor(gexatac_merge$sub.cluster,
                             levels = order_factor)
Idents(gexatac_merge) <- gexatac_merge@meta.data$sub.cluster

# 保存
# saveRDS(gexatac_merge, "./data/gexatac_merge_RPCA.rds")
gexatac_merge <- readRDS("gexatac_merge_RPCA.rds")
DimPlot(gexatac_merge, reduction = "umap.rpca", group.by = "sub.cluster", label = T)

DimPlot(gexatac_merge, reduction = "umap.rpca", group.by = c("sub.cluster"), combine = FALSE,label = T,
        split.by = "lineage", ncol=2)
ggsave("1_umap_RPCA_sub.cluster_facet_lineage.pdf", height = 4, width = 8)


# Table
table1 <- gexatac_merge@meta.data %>% select(sub.cluster, sample) %>% table() %>% as.data.frame()
table1 %>%
  ggplot(aes(x = sub.cluster, y = Freq, fill = sample)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_classic() +
  scale_fill_brewer(palette = "Set1")
ggsave("table_clusters_sample.pdf", height = 4, width = 8)


meta <- gexatac_merge@meta.data 
meta %>% dplyr::filter(nFeature_RNA < 3000) %>% select(seurat_clusters) %>% table



# All markers
library(presto)
markers_rna <- presto:::wilcoxauc.Seurat(X = gexatac_merge, group_by = 'sub.cluster', assay = 'data', seurat_assay = 'RNA')
write.table(markers_rna, "markers_rna.txt", quote = F, sep = "\t", row.names = F)


Idents(gexatac_merge)
markers <- FindAllMarkers(gexatac_merge, only.pos = TRUE)
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/v10/v10.2.1/markers")
write.table(markers, "markers.txt", quote = F, sep = "\t", row.names = F)

setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/v10/v10.2.1/markers")
markers <- read.table("markers.txt", header = T)
markers.filt <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
write.table(markers.filt, "markers.filt.txt", quote = F, sep = "\t", row.names = F)


markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10


cluster.averages <- AverageExpression(gexatac_merge, return.seurat = TRUE, assays = "RNA")
cluster.averages

DoHeatmap(cluster.averages, features = top10$gene, draw.lines = FALSE)  + NoLegend()
ggsave("heatmap_markeres_log2FC_top10_pseudobulk.png", height = 25, width = 8)
DoHeatmap(cluster.averages, features = top10$gene, draw.lines = FALSE)  + NoLegend()  + RotatedAxis() + coord_flip()
ggsave("heatmap_markeres_log2FC_top10_pseudobulk_rotated.png", height = 4, width = 35)

cluster14.markers <- FindMarkers(gexatac_merge, ident.1 = 14, ident.2 = 11)

cluster14.markers %>% dplyr::filter(avg_log2FC > 1) %>%  slice_head(n = 10) 
head(cluster14.markers, n = 10)


new.cluster.ids <- c("0_EBf2+_Mes",
                     "1_Barx1_Mes",
                     "2_Dkk2+_Mes",
                     "3_VIC",
                     "4_Mature_SMC",
                     "5_Epicardium",
                     "6_Foxf1+_Mes",
                     "7_Ednrb+_EC",
                     "8_SC",
                     "9_Aplnr+_EC",
                     "10_Myocyte",
                     "11_Macrpphage",
                     "12_Neuron",
                     "13_Cardiomyocyte",
                     "14_Lymph_progenitor",
                     "15_Irx1+_Mes",
                     "16_RBC",
                     "17_Arterial_EC",
                     "18_Epithelium",
                     "19_LEC",
                     "20_Immature_SMC"
                     )
names(new.cluster.ids) <- levels(gexatac_merge)
gexatac_merge <- RenameIdents(gexatac_merge, new.cluster.ids)
gexatac_merge@meta.data$celltype <- Idents(gexatac_merge)
DimPlot(gexatac_merge, reduction = "umap.rpca", group.by = "celltype",label = TRUE, pt.size = 0.5) 
  # guides(color=guide_legend(ncol=1)
ggsave("1_umap_RPCA_celltype.pdf", height = 6, width = 8)






