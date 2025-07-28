


library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)
library(tidyverse)
library(GenomicRanges)
library(future)



########################
# E11.5_EYFP
# the 10x hdf5 file contains both data types.
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E11.5/NCC/original_data")
inputdata.10x <- Read10X_h5("filtered_feature_bc_matrix.h5")
# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks
# Create Seurat object
E11.5_EYFP <- CreateSeuratObject(counts = rna_counts)
E11.5_EYFP[["percent.mt"]] <- PercentageFeatureSet(E11.5_EYFP, pattern = "^mt-")

########################
# E11.5_non_EYFP
# the 10x hdf5 file contains both data types.
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E11.5/non_NCC/original_data")
inputdata.10x <- Read10X_h5("filtered_feature_bc_matrix.h5")
# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks
# Create Seurat object
E11.5_non_EYFP <- CreateSeuratObject(counts = rna_counts)
E11.5_non_EYFP[["percent.mt"]] <- PercentageFeatureSet(E11.5_non_EYFP, pattern = "^mt-")

########################
# E12.5_EYFP
# the 10x hdf5 file contains both data types.
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E12.5/NCC/original_data")
inputdata.10x <- Read10X_h5("filtered_feature_bc_matrix.h5")
# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks
# Create Seurat object
E12.5_EYFP <- CreateSeuratObject(counts = rna_counts)
E12.5_EYFP[["percent.mt"]] <- PercentageFeatureSet(E12.5_EYFP, pattern = "^mt-")

########################
# E12.5_non_EYFP
# the 10x hdf5 file contains both data types.
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E12.5/non_NCC/original_data")
inputdata.10x <- Read10X_h5("filtered_feature_bc_matrix.h5")
# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks
# Create Seurat object
E12.5_non_EYFP <- CreateSeuratObject(counts = rna_counts)
E12.5_non_EYFP[["percent.mt"]] <- PercentageFeatureSet(E12.5_non_EYFP, pattern = "^mt-")


########################
# E17.5_EYFP
# the 10x hdf5 file contains both data types.
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E17.5/NCC/original_data")
inputdata.10x <- Read10X_h5("filtered_feature_bc_matrix.h5")
# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks
# Create Seurat object
E17.5_EYFP <- CreateSeuratObject(counts = rna_counts)
E17.5_EYFP[["percent.mt"]] <- PercentageFeatureSet(E17.5_EYFP, pattern = "^mt-")

#########################
#E14.5_EYFP
# the 10x hdf5 file contains both data types.
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E14.5/NCC/original_data")
inputdata.10x <- Read10X_h5("filtered_feature_bc_matrix.h5")
# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks
# Create Seurat object
E14.5_EYFP <- CreateSeuratObject(counts = rna_counts)
E14.5_EYFP[["percent.mt"]] <- PercentageFeatureSet(E14.5_EYFP, pattern = "^mt-")

#########################
# E17.5_non_EYFP
# the 10x hdf5 file contains both data types.
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E17.5/non_NCC/original_data")
inputdata.10x <- Read10X_h5("filtered_feature_bc_matrix.h5")
# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks
# Create Seurat object
E17.5_non_EYFP <- CreateSeuratObject(counts = rna_counts)
E17.5_non_EYFP[["percent.mt"]] <- PercentageFeatureSet(E17.5_non_EYFP, pattern = "^mt-")

#########################
#E14.5_non_EYFP
# the 10x hdf5 file contains both data types.
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E14.5/non_NCC/original_data")
inputdata.10x <- Read10X_h5("filtered_feature_bc_matrix.h5")
# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks
# Create Seurat object
E14.5_non_EYFP <- CreateSeuratObject(counts = rna_counts)
E14.5_non_EYFP[["percent.mt"]] <- PercentageFeatureSet(E14.5_non_EYFP, pattern = "^mt-")


rm(atac_counts)
rm(rna_counts)
rm(inputdata.10x)
# 
# Creating a common peak set
# 
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

# read in peak sets
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E11.5/NCC/original_data")
peaks.E11.5_EYFP <- read.table(
  file = "atac_peaks.bed",
  col.names = c("chr", "start", "end")
)
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E11.5/non_NCC/original_data")
peaks.E11.5_non_EYFP <- read.table(
  file = "atac_peaks.bed",
  col.names = c("chr", "start", "end")
)
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E12.5/NCC/original_data")
peaks.E12.5_EYFP <- read.table(
  file = "atac_peaks.bed",
  col.names = c("chr", "start", "end")
)
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E12.5/non_NCC/original_data")
peaks.E12.5_non_EYFP <- read.table(
  file = "atac_peaks.bed",
  col.names = c("chr", "start", "end")
)
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E14.5/NCC/original_data")
peaks.E14.5_EYFP <- read.table(
  file = "atac_peaks.bed",
  col.names = c("chr", "start", "end")
)
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E14.5/non_NCC/original_data")
peaks.E14.5_non_EYFP <- read.table(
  file = "atac_peaks.bed",
  col.names = c("chr", "start", "end")
)
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E17.5/NCC/original_data")
peaks.E17.5_EYFP <- read.table(
  file = "atac_peaks.bed",
  col.names = c("chr", "start", "end")
)
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E17.5/non_NCC/original_data")
peaks.E17.5_non_EYFP <- read.table(
  file = "atac_peaks.bed",
  col.names = c("chr", "start", "end")
)


# convert to genomic ranges
gr.E11.5_EYFP <- makeGRangesFromDataFrame(peaks.E11.5_EYFP)
gr.E11.5_non_EYFP <- makeGRangesFromDataFrame(peaks.E11.5_non_EYFP)
gr.E12.5_EYFP <- makeGRangesFromDataFrame(peaks.E12.5_EYFP)
gr.E12.5_non_EYFP <- makeGRangesFromDataFrame(peaks.E12.5_non_EYFP)
gr.E17.5_EYFP <- makeGRangesFromDataFrame(peaks.E17.5_EYFP)
gr.E14.5_EYFP <- makeGRangesFromDataFrame(peaks.E14.5_EYFP)
gr.E17.5_non_EYFP <- makeGRangesFromDataFrame(peaks.E17.5_non_EYFP)
gr.E14.5_non_EYFP <- makeGRangesFromDataFrame(peaks.E14.5_non_EYFP)


gr.E11.5_EYFP <- gr.E11.5_EYFP[seqnames(gr.E11.5_EYFP) %in% standardChromosomes(gr.E11.5_EYFP)]
gr.E11.5_non_EYFP <- gr.E11.5_non_EYFP[seqnames(gr.E11.5_non_EYFP) %in% standardChromosomes(gr.E11.5_non_EYFP)]
gr.E12.5_EYFP <- gr.E12.5_EYFP[seqnames(gr.E12.5_EYFP) %in% standardChromosomes(gr.E12.5_EYFP)]
gr.E12.5_non_EYFP <- gr.E12.5_non_EYFP[seqnames(gr.E12.5_non_EYFP) %in% standardChromosomes(gr.E12.5_non_EYFP)]
gr.E17.5_EYFP <- gr.E17.5_EYFP[seqnames(gr.E17.5_EYFP) %in% standardChromosomes(gr.E17.5_EYFP)]
gr.E14.5_EYFP <- gr.E14.5_EYFP[seqnames(gr.E14.5_EYFP) %in% standardChromosomes(gr.E14.5_EYFP)]
gr.E17.5_non_EYFP <- gr.E17.5_non_EYFP[seqnames(gr.E17.5_non_EYFP) %in% standardChromosomes(gr.E17.5_non_EYFP)]
gr.E14.5_non_EYFP <- gr.E14.5_non_EYFP[seqnames(gr.E14.5_non_EYFP) %in% standardChromosomes(gr.E14.5_non_EYFP)]

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- GenomicRanges::reduce(x = c(gr.E11.5_EYFP, gr.E11.5_non_EYFP, gr.E12.5_EYFP, gr.E12.5_non_EYFP,
                                              gr.E14.5_EYFP, gr.E14.5_non_EYFP, gr.E17.5_EYFP, gr.E17.5_non_EYFP))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
peakwidths %>% summary
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/data")
# 保存
saveRDS(combined.peaks, "combined.peaks.rds")

# 
# Create Fragment objects
# 
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E11.5/NCC/original_data")
frags.E11.5_EYFP <- CreateFragmentObject(
  path = "atac_fragments.tsv.gz"
)
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E11.5/non_NCC/original_data")
frags.E11.5_non_EYFP <- CreateFragmentObject(
  path = "atac_fragments.tsv.gz"
)

setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E12.5/NCC/original_data")
frags.E12.5_EYFP <- CreateFragmentObject(
  path = "atac_fragments.tsv.gz"
)
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E12.5/non_NCC/original_data")
frags.E12.5_non_EYFP <- CreateFragmentObject(
  path = "atac_fragments.tsv.gz"
)

setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E14.5/NCC/original_data")
frags.E14.5_EYFP <- CreateFragmentObject(
  path = "atac_fragments.tsv.gz"
)
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E14.5/non_NCC/original_data")
frags.E14.5_non_EYFP <- CreateFragmentObject(
  path = "atac_fragments.tsv.gz"
)

setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E17.5/NCC/original_data")
frags.E17.5_EYFP <- CreateFragmentObject(
  path = "atac_fragments.tsv.gz"
)
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E17.5/non_NCC/original_data")
frags.E17.5_non_EYFP <- CreateFragmentObject(
  path = "atac_fragments.tsv.gz"
)


# 
# Quantify peaks in each dataset
# 
atac.counts_E11.5_EYFP <- FeatureMatrix(
  fragments = frags.E11.5_EYFP,
  features = combined.peaks
)
atac.counts_E11.5_non_EYFP <- FeatureMatrix(
  fragments = frags.E11.5_non_EYFP,
  features = combined.peaks
)
atac.counts_E12.5_EYFP <- FeatureMatrix(
  fragments = frags.E12.5_EYFP,
  features = combined.peaks
)
atac.counts_E12.5_non_EYFP <- FeatureMatrix(
  fragments = frags.E12.5_non_EYFP,
  features = combined.peaks
)
atac.counts_E17.5_EYFP <- FeatureMatrix(
  fragments = frags.E17.5_EYFP,
  features = combined.peaks
)
atac.counts_E14.5_EYFP <- FeatureMatrix(
  fragments = frags.E14.5_EYFP,
  features = combined.peaks
)
atac.counts_E17.5_non_EYFP <- FeatureMatrix(
  fragments = frags.E17.5_non_EYFP,
  features = combined.peaks
)
atac.counts_E14.5_non_EYFP <- FeatureMatrix(
  fragments = frags.E14.5_non_EYFP,
  features = combined.peaks
)



# カウント値の細胞数をフィルタリグする　
atac.counts_E11.5_EYFP <- atac.counts_E11.5_EYFP[, colnames(E11.5_EYFP)]
atac.counts_E11.5_non_EYFP <- atac.counts_E11.5_non_EYFP[, colnames(E11.5_non_EYFP)]
atac.counts_E12.5_EYFP <- atac.counts_E12.5_EYFP[, colnames(E12.5_EYFP)]
atac.counts_E12.5_non_EYFP <- atac.counts_E12.5_non_EYFP[, colnames(E12.5_non_EYFP)]
atac.counts_E17.5_EYFP <- atac.counts_E17.5_EYFP[, colnames(E17.5_EYFP)]
atac.counts_E14.5_EYFP <- atac.counts_E14.5_EYFP[, colnames(E14.5_EYFP)]
atac.counts_E17.5_non_EYFP <- atac.counts_E17.5_non_EYFP[, colnames(E17.5_non_EYFP)]
atac.counts_E14.5_non_EYFP <- atac.counts_E14.5_non_EYFP[, colnames(E14.5_non_EYFP)]

# dim(atac.counts_E11.5_EYFP)
# [1] 292070    876
# dim(atac.counts_E11.5_non_EYFP)
# [1] 292070    924
# dim(atac.counts_E12.5_EYFP)
# [1] 292070    12648
# dim(atac.counts_E12.5_non_EYFP)
# [1] 292070    3289
# dim(atac.counts_E17.5_EYFP)
# # [1] 126035    344
# dim(atac.counts_E14.5_EYFP)
# # [1] 126035    593
# dim(atac.counts_E17.5_non_EYFP)
# # [1] 156465   1771
# dim(atac.counts_E14.5_non_EYFP)
# # [1] 175400   1145


# 
# Create the objects
# 
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79) #mouse
seqlevelsStyle(annotations) <- "UCSC" 
genome(annotations) <- "mm10"

setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E11.5/NCC/original_data")
frag.file <- "atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac.counts_E11.5_EYFP,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = frag.file,
  min.cells = 5,
  annotation = annotations
)
E11.5_EYFP[["ATAC"]] <- chrom_assay

setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E11.5/non_NCC/original_data")
frag.file <- "atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac.counts_E11.5_non_EYFP,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = frag.file,
  min.cells = 5,
  annotation = annotations
)
E11.5_non_EYFP[["ATAC"]] <- chrom_assay

setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E12.5/NCC/original_data")
frag.file <- "atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac.counts_E12.5_EYFP,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = frag.file,
  min.cells = 5,
  annotation = annotations
)
E12.5_EYFP[["ATAC"]] <- chrom_assay

setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E12.5/non_NCC/original_data")
frag.file <- "atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac.counts_E12.5_non_EYFP,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = frag.file,
  min.cells = 5,
  annotation = annotations
)
E12.5_non_EYFP[["ATAC"]] <- chrom_assay

setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E14.5/NCC/original_data")
frag.file <- "atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac.counts_E14.5_EYFP,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = frag.file,
  min.cells = 5,
  annotation = annotations
)
E14.5_EYFP[["ATAC"]] <- chrom_assay

setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E14.5/non_NCC/original_data")
frag.file <- "atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac.counts_E14.5_non_EYFP,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = frag.file,
  min.cells = 5,
  annotation = annotations
)
E14.5_non_EYFP[["ATAC"]] <- chrom_assay

setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E17.5/NCC/original_data")
frag.file <- "atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac.counts_E17.5_EYFP,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = frag.file,
  min.cells = 5,
  annotation = annotations
)
E17.5_EYFP[["ATAC"]] <- chrom_assay

setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E17.5/non_NCC/original_data")
frag.file <- "atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac.counts_E17.5_non_EYFP,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = frag.file,
  min.cells = 5,
  annotation = annotations
)
E17.5_non_EYFP[["ATAC"]] <- chrom_assay


E11.5_EYFP@meta.data$orig.ident <- "E11.5_EYFP"
Idents(E11.5_EYFP) <- E11.5_EYFP@meta.data$orig.ident 
E11.5_non_EYFP@meta.data$orig.ident <- "E11.5_non_EYFP"
Idents(E11.5_non_EYFP) <- E11.5_non_EYFP@meta.data$orig.ident 
E12.5_EYFP@meta.data$orig.ident <- "E12.5_EYFP"
Idents(E12.5_EYFP) <- E12.5_EYFP@meta.data$orig.ident 
E12.5_non_EYFP@meta.data$orig.ident <- "E12.5_non_EYFP"
Idents(E12.5_non_EYFP) <- E12.5_non_EYFP@meta.data$orig.ident 
E17.5_EYFP@meta.data$orig.ident <- "E17.5_EYFP"
Idents(E17.5_EYFP) <- E17.5_EYFP@meta.data$orig.ident
E14.5_EYFP@meta.data$orig.ident <- "E14.5_EYFP"
Idents(E14.5_EYFP) <-E14.5_EYFP@meta.data$orig.ident
E17.5_non_EYFP@meta.data$orig.ident <- "E17.5_non_EYFP"
Idents(E17.5_non_EYFP) <- E17.5_non_EYFP@meta.data$orig.ident
E14.5_non_EYFP@meta.data$orig.ident <- "E14.5_non_EYFP"
Idents(E14.5_non_EYFP) <-E14.5_non_EYFP@meta.data$orig.ident


# save
saveRDS(E11.5_EYFP, "E11.5_EYFP_raw.v1.rds")
saveRDS(E11.5_non_EYFP, "E11.5_non_EYFP_raw.v1.rds")
saveRDS(E12.5_EYFP, "E12.5_EYFP_raw.v1.rds")
saveRDS(E12.5_non_EYFP, "E12.5_non_EYFP_raw.v1.rds")
saveRDS(E17.5_EYFP, "E17.5_EYFP_raw.v1.rds")
saveRDS(E14.5_EYFP, "E14.5_EYFP_raw.v1.rds")
saveRDS(E17.5_non_EYFP, "E17.5_non_EYFP_raw.v1.rds")
saveRDS(E14.5_non_EYFP, "E14.5_non_EYFP_raw.v1.rds")


# read velocyo data
library(velocyto.R)
library(Seurat)
library(SeuratWrappers)
# E11.5
setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E11.5/NCC/original_data/velocyto")
ldat <- ReadVelocity(file = "gex_possorted_bam_MDLLT.loom")
ldat <- lapply(ldat,function(x) {
  colnames(x) <-  gsub("gex_possorted_bam_MDLLT:","",gsub(".*:","",colnames(x)))
  colnames(x) <-  gsub("x","-1",gsub(".*:","",colnames(x)))
  x
})
E11.5_EYFP_velo<- as.Seurat(x = ldat)
dim(E11.5_EYFP_velo)
# [1] 32286   876
# spliced
spliced <- E11.5_EYFP_velo@assays$spliced@counts %>% t() %>% as.data.frame()
spliced <- cbind(rownames(spliced), spliced)
colnames(spliced)[1] <- "cell"
spliced <- spliced %>% arrange(cell) %>% t() %>% as.data.frame()
spliced <- spliced[-1,]
spliced <- CreateAssayObject(spliced)
spliced
E11.5_EYFP[["spliced"]] <- spliced
# unspliced
unspliced <- E11.5_EYFP_velo@assays$unspliced@counts %>% t() %>% as.data.frame()
unspliced <- cbind(rownames(unspliced), unspliced)
colnames(unspliced)[1] <- "cell"
unspliced <- unspliced %>% arrange(cell) %>% t() %>% as.data.frame()
unspliced <- unspliced[-1,]
unspliced <- CreateAssayObject(unspliced)
unspliced
E11.5_EYFP[["unspliced"]] <- unspliced
# ambiguous
ambiguous <- E11.5_EYFP_velo@assays$ambiguous@counts %>% t() %>% as.data.frame()
ambiguous <- cbind(rownames(ambiguous), ambiguous)
colnames(ambiguous)[1] <- "cell"
ambiguous <- ambiguous %>% arrange(cell) %>% t() %>% as.data.frame()
ambiguous <- ambiguous[-1,]
ambiguous <- CreateAssayObject(ambiguous)
ambiguous
E11.5_EYFP[["ambiguous"]] <- ambiguous

E11.5_EYFP
# これでvelocityデータが搭載された

setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E11.5/non_NCC/original_data/velocyto")
ldat <- ReadVelocity(file = "gex_possorted_bam_FK6VC.loom")
ldat <- lapply(ldat,function(x) {
  colnames(x) <-  gsub("gex_possorted_bam_FK6VC:","",gsub(".*:","",colnames(x)))
  colnames(x) <-  gsub("x","-1",gsub(".*:","",colnames(x)))
  x
})
E11.5_non_EYFP_velo<- as.Seurat(x = ldat)
dim(E11.5_non_EYFP_velo)
# [1]  32286  924
# spliced
spliced <- E11.5_non_EYFP_velo@assays$spliced@counts %>% t() %>% as.data.frame()
spliced <- cbind(rownames(spliced), spliced)
colnames(spliced)[1] <- "cell"
spliced <- spliced %>% arrange(cell) %>% t() %>% as.data.frame()
spliced <- spliced[-1,]
spliced <- CreateAssayObject(spliced)
spliced
E11.5_non_EYFP[["spliced"]] <- spliced
# unspliced
unspliced <- E11.5_non_EYFP_velo@assays$unspliced@counts %>% t() %>% as.data.frame()
unspliced <- cbind(rownames(unspliced), unspliced)
colnames(unspliced)[1] <- "cell"
unspliced <- unspliced %>% arrange(cell) %>% t() %>% as.data.frame()
unspliced <- unspliced[-1,]
unspliced <- CreateAssayObject(unspliced)
unspliced
E11.5_non_EYFP[["unspliced"]] <- unspliced
# ambiguous
ambiguous <- E11.5_non_EYFP_velo@assays$ambiguous@counts %>% t() %>% as.data.frame()
ambiguous <- cbind(rownames(ambiguous), ambiguous)
colnames(ambiguous)[1] <- "cell"
ambiguous <- ambiguous %>% arrange(cell) %>% t() %>% as.data.frame()
ambiguous <- ambiguous[-1,]
ambiguous <- CreateAssayObject(ambiguous)
ambiguous
E11.5_non_EYFP[["ambiguous"]] <- ambiguous

E11.5_non_EYFP

setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E12.5/NCC/original_data/velocyto")
ldat <- ReadVelocity(file = "gex_possorted_bam_25MPX.loom")
ldat <- lapply(ldat,function(x) {
  colnames(x) <-  gsub("gex_possorted_bam_25MPX:","",gsub(".*:","",colnames(x)))
  colnames(x) <-  gsub("x","-1",gsub(".*:","",colnames(x)))
  x
})
E12.5_EYFP_velo<- as.Seurat(x = ldat)
dim(E12.5_EYFP_velo)
# [1] 32286 12648
# spliced
spliced <- E12.5_EYFP_velo@assays$spliced@counts %>% t() %>% as.data.frame()
spliced <- cbind(rownames(spliced), spliced)
colnames(spliced)[1] <- "cell"
spliced <- spliced %>% arrange(cell) %>% t() %>% as.data.frame()
spliced <- spliced[-1,]
spliced <- CreateAssayObject(spliced)
spliced
E12.5_EYFP[["spliced"]] <- spliced
# unspliced
unspliced <- E12.5_EYFP_velo@assays$unspliced@counts %>% t() %>% as.data.frame()
unspliced <- cbind(rownames(unspliced), unspliced)
colnames(unspliced)[1] <- "cell"
unspliced <- unspliced %>% arrange(cell) %>% t() %>% as.data.frame()
unspliced <- unspliced[-1,]
unspliced <- CreateAssayObject(unspliced)
unspliced
E12.5_EYFP[["unspliced"]] <- unspliced
# ambiguous
ambiguous <- E12.5_EYFP_velo@assays$ambiguous@counts %>% t() %>% as.data.frame()
ambiguous <- cbind(rownames(ambiguous), ambiguous)
colnames(ambiguous)[1] <- "cell"
ambiguous <- ambiguous %>% arrange(cell) %>% t() %>% as.data.frame()
ambiguous <- ambiguous[-1,]
ambiguous <- CreateAssayObject(ambiguous)
ambiguous
E12.5_EYFP[["ambiguous"]] <- ambiguous

E12.5_EYFP


setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E12.5/non_NCC/original_data/velocyto")
ldat <- ReadVelocity(file = "gex_possorted_bam_M1RC8.loom")
ldat <- lapply(ldat,function(x) {
  colnames(x) <-  gsub("gex_possorted_bam_M1RC8:","",gsub(".*:","",colnames(x)))
  colnames(x) <-  gsub("x","-1",gsub(".*:","",colnames(x)))
  x
})
E12.5_non_EYFP_velo<- as.Seurat(x = ldat)
dim(E12.5_non_EYFP_velo)
# [1]  32286  3289
# spliced
spliced <- E12.5_non_EYFP_velo@assays$spliced@counts %>% t() %>% as.data.frame()
spliced <- cbind(rownames(spliced), spliced)
colnames(spliced)[1] <- "cell"
spliced <- spliced %>% arrange(cell) %>% t() %>% as.data.frame()
spliced <- spliced[-1,]
spliced <- CreateAssayObject(spliced)
spliced
E12.5_non_EYFP[["spliced"]] <- spliced
# unspliced
unspliced <- E12.5_non_EYFP_velo@assays$unspliced@counts %>% t() %>% as.data.frame()
unspliced <- cbind(rownames(unspliced), unspliced)
colnames(unspliced)[1] <- "cell"
unspliced <- unspliced %>% arrange(cell) %>% t() %>% as.data.frame()
unspliced <- unspliced[-1,]
unspliced <- CreateAssayObject(unspliced)
unspliced
E12.5_non_EYFP[["unspliced"]] <- unspliced
# ambiguous
ambiguous <- E12.5_non_EYFP_velo@assays$ambiguous@counts %>% t() %>% as.data.frame()
ambiguous <- cbind(rownames(ambiguous), ambiguous)
colnames(ambiguous)[1] <- "cell"
ambiguous <- ambiguous %>% arrange(cell) %>% t() %>% as.data.frame()
ambiguous <- ambiguous[-1,]
ambiguous <- CreateAssayObject(ambiguous)
ambiguous
E12.5_non_EYFP[["ambiguous"]] <- ambiguous

E12.5_non_EYFP



setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E14.5/NCC/original_data/velocyto")    
ldat <- ReadVelocity(file = "gex_possorted_bam_FVUJ5.loom")
ldat <- lapply(ldat,function(x) {
  colnames(x) <-  gsub("gex_possorted_bam_FVUJ5:","",gsub(".*:","",colnames(x)))
  colnames(x) <-  gsub("x","-1",gsub(".*:","",colnames(x)))
  x
})
E14.5_EYFP_velo<- as.Seurat(x = ldat)
dim(E14.5_EYFP_velo)
# [1] 32285   593
# spliced
spliced <- E14.5_EYFP_velo@assays$spliced@counts %>% t() %>% as.data.frame()
spliced <- cbind(rownames(spliced), spliced)
colnames(spliced)[1] <- "cell"
spliced <- spliced %>% arrange(cell) %>% t() %>% as.data.frame()
spliced <- spliced[-1,]
spliced <- CreateAssayObject(spliced)
spliced
E14.5_EYFP[["spliced"]] <- spliced
# unspliced
unspliced <- E14.5_EYFP_velo@assays$unspliced@counts %>% t() %>% as.data.frame()
unspliced <- cbind(rownames(unspliced), unspliced)
colnames(unspliced)[1] <- "cell"
unspliced <- unspliced %>% arrange(cell) %>% t() %>% as.data.frame()
unspliced <- unspliced[-1,]
unspliced <- CreateAssayObject(unspliced)
unspliced
E14.5_EYFP[["unspliced"]] <- unspliced
# ambiguous
ambiguous <- E14.5_EYFP_velo@assays$ambiguous@counts %>% t() %>% as.data.frame()
ambiguous <- cbind(rownames(ambiguous), ambiguous)
colnames(ambiguous)[1] <- "cell"
ambiguous <- ambiguous %>% arrange(cell) %>% t() %>% as.data.frame()
ambiguous <- ambiguous[-1,]
ambiguous <- CreateAssayObject(ambiguous)
ambiguous
E14.5_EYFP[["ambiguous"]] <- ambiguous

E14.5_EYFP


setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E17.5/NCC/original_data/velocyto")
ldat <- ReadVelocity(file = "gex_possorted_bam_OYKGS.loom")
ldat <- lapply(ldat,function(x) {
  colnames(x) <-  gsub("gex_possorted_bam_OYKGS:","",gsub(".*:","",colnames(x)))
  colnames(x) <-  gsub("x","-1",gsub(".*:","",colnames(x)))
  x
})
E17.5_EYFP_velo<- as.Seurat(x = ldat)
dim(E17.5_EYFP_velo)
# [1] 32285   344
# spliced
spliced <- E17.5_EYFP_velo@assays$spliced@counts %>% t() %>% as.data.frame()
spliced <- cbind(rownames(spliced), spliced)
colnames(spliced)[1] <- "cell"
spliced <- spliced %>% arrange(cell) %>% t() %>% as.data.frame()
spliced <- spliced[-1,]
spliced <- CreateAssayObject(spliced)
spliced
intersect(colnames(spliced) %>% head,  colnames(E17.5_EYFP) %>% head)

E17.5_EYFP[["spliced"]] <- spliced
# unspliced
unspliced <- E17.5_EYFP_velo@assays$unspliced@counts %>% t() %>% as.data.frame()
unspliced <- cbind(rownames(unspliced), unspliced)
colnames(unspliced)[1] <- "cell"
unspliced <- unspliced %>% arrange(cell) %>% t() %>% as.data.frame()
unspliced <- unspliced[-1,]
unspliced <- CreateAssayObject(unspliced)
unspliced
E17.5_EYFP[["unspliced"]] <- unspliced
# ambiguous
ambiguous <- E17.5_EYFP_velo@assays$ambiguous@counts %>% t() %>% as.data.frame()
ambiguous <- cbind(rownames(ambiguous), ambiguous)
colnames(ambiguous)[1] <- "cell"
ambiguous <- ambiguous %>% arrange(cell) %>% t() %>% as.data.frame()
ambiguous <- ambiguous[-1,]
ambiguous <- CreateAssayObject(ambiguous)
E17.5_EYFP[["ambiguous"]] <- ambiguous

E17.5_EYFP


setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E14.5/non_NCC/original_data/velocyto")
ldat <- ReadVelocity(file = "gex_possorted_bam_R7GUY.loom")
ldat <- lapply(ldat,function(x) {
  colnames(x) <-  gsub("gex_possorted_bam_R7GUY:","",gsub(".*:","",colnames(x)))
  colnames(x) <-  gsub("x","-1",gsub(".*:","",colnames(x)))
  x
})
E14.5_non_EYFP_velo<- as.Seurat(x = ldat)
dim(E14.5_non_EYFP_velo)
# [1] 32285   1145
# spliced
spliced <- E14.5_non_EYFP_velo@assays$spliced@counts %>% t() %>% as.data.frame()
spliced <- cbind(rownames(spliced), spliced)
colnames(spliced)[1] <- "cell"
spliced <- spliced %>% arrange(cell) %>% t() %>% as.data.frame()
spliced <- spliced[-1,]
spliced <- CreateAssayObject(spliced)
spliced
E14.5_non_EYFP[["spliced"]] <- spliced
# unspliced
unspliced <- E14.5_non_EYFP_velo@assays$unspliced@counts %>% t() %>% as.data.frame()
unspliced <- cbind(rownames(unspliced), unspliced)
colnames(unspliced)[1] <- "cell"
unspliced <- unspliced %>% arrange(cell) %>% t() %>% as.data.frame()
unspliced <- unspliced[-1,]
unspliced <- CreateAssayObject(unspliced)
unspliced
E14.5_non_EYFP[["unspliced"]] <- unspliced
# ambiguous
ambiguous <- E14.5_non_EYFP_velo@assays$ambiguous@counts %>% t() %>% as.data.frame()
ambiguous <- cbind(rownames(ambiguous), ambiguous)
colnames(ambiguous)[1] <- "cell"
ambiguous <- ambiguous %>% arrange(cell) %>% t() %>% as.data.frame()
ambiguous <- ambiguous[-1,]
ambiguous <- CreateAssayObject(ambiguous)
ambiguous
E14.5_non_EYFP[["ambiguous"]] <- ambiguous

E14.5_non_EYFP


setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Multiome_Wnt1_E17.5/non_NCC/original_data/velocyto")
ldat <- ReadVelocity(file = "gex_possorted_bam_PTX3M.loom")
ldat <- lapply(ldat,function(x) {
  colnames(x) <-  gsub("gex_possorted_bam_PTX3M:","",gsub(".*:","",colnames(x)))
  colnames(x) <-  gsub("x","-1",gsub(".*:","",colnames(x)))
  x
})
E17.5_non_EYFP_velo<- as.Seurat(x = ldat)
dim(E17.5_non_EYFP_velo)
# [1] 32285   1771
# spliced
spliced <- E17.5_non_EYFP_velo@assays$spliced@counts %>% t() %>% as.data.frame()
spliced <- cbind(rownames(spliced), spliced)
colnames(spliced)[1] <- "cell"
spliced <- spliced %>% arrange(cell) %>% t() %>% as.data.frame()
spliced <- spliced[-1,]
spliced <- CreateAssayObject(spliced)
spliced
intersect(colnames(spliced) %>% head,  colnames(E17.5_non_EYFP) %>% head)

E17.5_non_EYFP[["spliced"]] <- spliced
# unspliced
unspliced <- E17.5_non_EYFP_velo@assays$unspliced@counts %>% t() %>% as.data.frame()
unspliced <- cbind(rownames(unspliced), unspliced)
colnames(unspliced)[1] <- "cell"
unspliced <- unspliced %>% arrange(cell) %>% t() %>% as.data.frame()
unspliced <- unspliced[-1,]
unspliced <- CreateAssayObject(unspliced)
unspliced
E17.5_non_EYFP[["unspliced"]] <- unspliced
# ambiguous
ambiguous <- E17.5_non_EYFP_velo@assays$ambiguous@counts %>% t() %>% as.data.frame()
ambiguous <- cbind(rownames(ambiguous), ambiguous)
colnames(ambiguous)[1] <- "cell"
ambiguous <- ambiguous %>% arrange(cell) %>% t() %>% as.data.frame()
ambiguous <- ambiguous[-1,]
ambiguous <- CreateAssayObject(ambiguous)
E17.5_non_EYFP[["ambiguous"]] <- ambiguous

E17.5_non_EYFP


setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/data")
# TSSSなど
# compute nucleosome signal score per cell
DefaultAssay(E11.5_EYFP) <- "ATAC"
E11.5_EYFP <- NucleosomeSignal(object = E11.5_EYFP)
# compute TSS enrichment score per cell
E11.5_EYFP <- TSSEnrichment(object = E11.5_EYFP, fast = FALSE)
VlnPlot(E11.5_EYFP, features = c("TSS.enrichment", "nucleosome_signal"))
ggsave("0_QC_violin_TSS_E11.5_EYFP.pdf", height = 4, width = 5)

DefaultAssay(E11.5_non_EYFP) <- "ATAC"
E11.5_non_EYFP <- NucleosomeSignal(object = E11.5_non_EYFP)
# compute TSS enrichment score per cell
E11.5_non_EYFP <- TSSEnrichment(object = E11.5_non_EYFP, fast = FALSE)
VlnPlot(E11.5_non_EYFP, features = c("TSS.enrichment", "nucleosome_signal"))
ggsave("0_QC_violin_TSS_E11.5_non_EYFP.pdf", height = 4, width = 5)


DefaultAssay(E12.5_EYFP) <- "ATAC"
E12.5_EYFP <- NucleosomeSignal(object = E12.5_EYFP)
# compute TSS enrichment score per cell
E12.5_EYFP <- TSSEnrichment(object = E12.5_EYFP, fast = FALSE)
VlnPlot(E12.5_EYFP, features = c("TSS.enrichment", "nucleosome_signal"))
ggsave("0_QC_violin_TSS_E12.5_EYFP.pdf", height = 4, width = 5)

DefaultAssay(E12.5_non_EYFP) <- "ATAC"
E12.5_non_EYFP <- NucleosomeSignal(object = E12.5_non_EYFP)
# compute TSS enrichment score per cell
E12.5_non_EYFP <- TSSEnrichment(object = E12.5_non_EYFP, fast = FALSE)
VlnPlot(E12.5_non_EYFP, features = c("TSS.enrichment", "nucleosome_signal"))
ggsave("0_QC_violin_TSS_E12.5_non_EYFP.pdf", height = 4, width = 5)


DefaultAssay(E14.5_EYFP) <- "ATAC"
E14.5_EYFP <- NucleosomeSignal(object = E14.5_EYFP)
# compute TSS enrichment score per cell
E14.5_EYFP <- TSSEnrichment(object = E14.5_EYFP, fast = FALSE)
VlnPlot(E14.5_EYFP, features = c("TSS.enrichment", "nucleosome_signal"))
ggsave("0_QC_violin_TSS_E14.5_EYFP.pdf", height = 4, width = 5)

DefaultAssay(E14.5_non_EYFP) <- "ATAC"
E14.5_non_EYFP <- NucleosomeSignal(object = E14.5_non_EYFP)
# compute TSS enrichment score per cell
E14.5_non_EYFP <- TSSEnrichment(object = E14.5_non_EYFP, fast = FALSE)
VlnPlot(E14.5_non_EYFP, features = c("TSS.enrichment", "nucleosome_signal"))
ggsave("0_QC_violin_TSS_E14.5_non_EYFP.pdf", height = 4, width = 5)


DefaultAssay(E17.5_EYFP) <- "ATAC"
E17.5_EYFP <- NucleosomeSignal(object = E17.5_EYFP)
# compute TSS enrichment score per cell
E17.5_EYFP <- TSSEnrichment(object = E17.5_EYFP, fast = FALSE)
VlnPlot(E17.5_EYFP, features = c("TSS.enrichment", "nucleosome_signal"))
ggsave("0_QC_violin_TSS_E17.5_EYFP.pdf", height = 4, width = 5)

DefaultAssay(E17.5_non_EYFP) <- "ATAC"
E17.5_non_EYFP <- NucleosomeSignal(object = E17.5_non_EYFP)
# compute TSS enrichment score per cell
E17.5_non_EYFP <- TSSEnrichment(object = E17.5_non_EYFP, fast = FALSE)
VlnPlot(E17.5_non_EYFP, features = c("TSS.enrichment", "nucleosome_signal"))
ggsave("0_QC_violin_TSS_E17.5_non_EYFP.pdf", height = 4, width = 5)


setwd("/Users/iwaseakiyasu/workspace/data/Multiome_Wnt1/Merge_Wnt1/data")
saveRDS(E11.5_EYFP, "E11.5_EYFP_raw.v1.rds")
saveRDS(E11.5_non_EYFP, "E11.5_non_EYFP_raw.v1.rds")
saveRDS(E12.5_EYFP, "E12.5_EYFP_raw.v1.rds")
saveRDS(E12.5_non_EYFP, "E12.5_non_EYFP_raw.v1.rds")
saveRDS(E17.5_EYFP, "E17.5_EYFP_raw.v1.rds")
saveRDS(E14.5_EYFP, "E14.5_EYFP_raw.v1.rds")
saveRDS(E17.5_non_EYFP, "E17.5_non_EYFP_raw.v1.rds")
saveRDS(E14.5_non_EYFP, "E14.5_non_EYFP_raw.v1.rds")






