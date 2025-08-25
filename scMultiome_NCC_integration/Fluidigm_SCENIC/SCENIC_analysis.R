# SCENIC

library(SCENIC)
library(SCopeLoomR)
library(dplyr)
library(ggplot2)
library(tidyverse)

library(BiocParallel)
library(AUCell)
library(RcisTarget)

######################
#Expression matrix
######################
setwd("~/Downloads/scRNA-seq/20200403_Seuratv3_based2/TPM")
data <- read.table("tpm_wo_cl11.txt", header =T, row.names =1)

datat <- t(data)
cell_name <- as.data.frame(rownames(datat))
colnames(cell_name) <- c("cell")
datat <- cbind(cell_name, datat)
colnames(datat) %>% head
data2 <- datat %>% arrange(cell) 
rownames(data2) <- as.character(data2[,1])
data2 <- data2 %>% t() %>% as.data.frame()
exprMar <- data2[-1,]
dim(exprMar)
#[1] 15124  284
cellInfo <- read.table("../Seurat/Seuratcds.metadata.txt",
                       header = T, row.names= 1, sep = "\t")
COI <- rownames(cellInfo)
class(COI)
GOI <- rownames(exprMar) %>% as.character()
class(GOI)
# Filtering
exprMar <- exprMar[GOI, COI]
dim(exprMar)
#[1] 15124   284



setwd("~/Downloads/scRNA-seq/20200403_Seuratv3_based2/SCENIC5")
#保存
saveRDS(cellInfo, file="int/cellInfo.Rds")



#https://qiita.com/hoxo_b/items/c569da6dbf568032e04a
ggColorHue <- function(n, l=65) {
  hues <- seq(15, 375, length=n+1)
  hcl(h=hues, l=l, c=100)[1:n]
}
ggColorHue(11)
#[1] "#F8766D" "#DB8E00" "#AEA200" "#64B200" "#00BD5C" "#00C1A7" "#00BADE" "#00A6FF"
#[9] "#B385FF" "#EF67EB" "#FF63B6"

colVars <- list(cluster=c("0" = "#F8766D",
                          "1" = "#DB8E00",
                          "2" = "#AEA200",
                          "3" = "#64B200",
                          "4" = "#00BD5C",
                          "5" = "#00C1A7",
                          "6" = "#00BADE",
                          "7" = "#00A6FF",
                          "8" = "#B385FF",
                          "9" = "#EF67EB",
                          "10" = "#FF63B6"))
#colVars$cluster <- colVars$cluster[intersect(names(colVars$cluster), cellInfo$seurat_clusters)]
plot.new(); legend(0,1, fill=colVars$cluster, legend=names(colVars$cluster))
#dir.create("int")
saveRDS(colVars, file="int/colVars.Rds")



###############################
#Initialize SCENIC settings
#https://rawcdn.githack.com/aertslab/SCENIC/0a4c96ed8d930edd8868f07428090f9dae264705/inst/doc/SCENIC_Running.html
###############################
org="mgi" # or hgnc, or dmel
#featheフォルダを~/Downloads/scRNA-seq/20190819_iws_TPM/TF解析/SCENIC/feather
dbDir="./feather" # RcisTarget databases location
myDatasetTitle="SCENIC example on Trajectory" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, 
                                  datasetTitle=myDatasetTitle, nCores=6) 
# Motif databases selected: 
#   mm9-500bp-upstream-7species.mc9nr.feather 
# mm9-tss-centered-10kb-7species.mc9nr.feather
# Modify if needed
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
# Databases:
# scenicOptions@settings$dbs <- c("mm9-5kb-mc8nr"="mm9-tss-centered-5kb-10species.mc8nr.feather")
# scenicOptions@settings$db_mcVersion <- "v8"

# Save to use at a later time...
saveRDS(scenicOptions, file="int/scenicOptions.Rds")


###############################
#Co-expression network
###############################


###############################
#Gene filter/selection
# (Adjust minimum values according to your dataset)
###############################
exprMat <- read.table("../TPM/tpm_wo_cl11.txt", header=T, row.names = 1)
exprMat <- as.matrix(exprMat)
dim(exprMat) 
#[1]  15124   284
exprMat[,1] %>% class
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
# Maximum value in the expression matrix: 158707.078368769
# Ratio of detected vs non-detected: 0.46
# Number of counts (in the dataset units) per gene:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0     830    3598   18778   11607 4492410 
# Number of cells in which each gene is detected:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0    20.0    72.0    89.5   141.0   284.0 
# 
# Number of genes left after applying the following filters (sequential):
#   14894	genes with counts per gene > 8.52
# 14801	genes detected in more than 2.84 cells
# 13268	genes available in RcisTarget database
# Gene list saved in int/1.1_genesKept.Rds
#前回は12095遺伝子がkeptだったが、13268遺伝子に変化

interestingGenes <- c("Sox9", "Scx", "Osr1", "Klf2", "Myh11", "Sox10")
interestingGenes[which(!interestingGenes %in% genesKept)]

exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)
#[1] 13268    284
rm(exprMat)

######################
#Correlation
######################
runCorrelation(exprMat_filtered, scenicOptions)


######################
## If launched in a new session, you will need to reload...
# setwd("...")
# loomPath <- "..."
# loom <- open_loom(loomPath, mode="r")
# exprMat <- get_dgem(loom)
# close_loom(loom)
# genesKept <- loadInt(scenicOptions, "genesKept")
# exprMat_filtered <- exprMat[genesKept,]
# library(SCENIC)
# scenicOptions <- readRDS("int/scenicOptions.Rds")
######################
# Optional: add log (if it is not logged/normalized already)
exprMat_filtered <- log2(exprMat_filtered+1) 

# Run GENIE3
runGenie3(exprMat_filtered, scenicOptions)
#14:26-
###int1.RDataとして保存
#save.image("~/Downloads/scRNA-seq/20200403_Seuratv3_based2/SCENIC2/int1.RData")


######################
#Build and score the GRN (runSCENIC_…)
#SCENICのワークフローは以下
#Build the gene regulatory network: 
#1. Get co-expression modules 
#2. Get regulons (with RcisTarget): TF motif analysis)
#Identify cell states: 
#3. Score GRN (regulons) in the cells (with AUCell) 
#4. Cluster cells according to the GRN activity
#必要であればリロードせよとのこと
#loom <- open_loom(loomPath, mode="r")
#exprMat <- get_dgem(loom)
#close_loom(loom)
## Optional: log expression (for TF expression plot, it does not affect any other calculation)
#logMat <- log2(exprMat+1)
#dim(exprMat)
######################
dim(exprMat_filtered)
#[1] 8231   284
library(SCENIC)
setwd("~/Downloads/scRNA-seq/20200403_Seuratv3_based2/SCENIC5")
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 6
scenicOptions@settings$seed <- 123

.openDev <- function(fileName, devType, ...)
{
  if(devType=="pdf")
    pdf(paste0(fileName, ".pdf"), ...)
  
  if(devType=="png")
    png(paste0(fileName, ".png"), ...)
  
  if(devType=="cairo_pfd") # similar to Cairo::CairoPDF?
    grDevices::cairo_pdf(paste0(fileName, ".pdf"), ...)
}

.openDevHeatmap <- function(fileName, devType)
{
  if(devType!="pdf") 
  {
    if(devType=="png") .openDev(fileName=fileName, devType=devType, width=1200,height=1200)
    if(devType!="png") .openDev(fileName=fileName, devType=devType)
    fileName <- NA
  }else{
    fileName <- paste0(fileName,".pdf")
  }
  return(fileName)
}

.closeDevHeatmap <- function(devType)
{
  if(devType!="pdf") 
  {
    dev.off()
  }
}


linkList <- loadInt(scenicOptions, "genie3ll")
a <- subset(linkList , TF =="Sox9")。

if(!all(colnames(linkList) == c("TF", "Target", "weight"))) stop('The link list colnames should be "TF", "Target", "weight"')

uniquePairs <- nrow(unique(linkList[,c("TF", "Target")]))
if(uniquePairs < nrow(linkList)) 
  stop("There are duplicated regulator-target (gene id/name) pairs in the input link list.")

msg <- paste0(format(Sys.time(), "%H:%M"), "\tCreating TF modules")
if(getSettings(scenicOptions, "verbose")) message(msg)

print(quantile(linkList$weight, probs=c(0.75, 0.90)))
# 75%         90% 
#   0.002053706 0.003024414 
.openDev(fileName=getIntName(scenicOptions, "genie3weighPlot"), 
         devType=getSettings(scenicOptions, "devType"))
plot(linkList$weight[1:1000000], type="l", ylim=c(0, max(linkList$weight)), main="Weight of the links",
     ylab="Weight", xlab="Links sorted decreasingly")
abline(h=0.001, col="blue") # Threshold
#sum(linkList$weight>0.001)/nrow(linkList)
dev.off()


# Keep only genes with weight > threshold
linkList_001 <- linkList[which(linkList[,"weight"]>getSettings(scenicOptions, "modules/weightThreshold")),]
if(getSettings(scenicOptions, "verbose")) message("Number of links between TFs and targets: ", nrow(linkList_001))

#### Create the gene-sets & save:
tfModules <- list()

linkList_001$TF <- as.character(linkList_001$TF)
linkList_001$Target <- as.character(linkList_001$Target)

### Create TF-modules:
# 1: Weight > 0.001 (filtered in previous step)
tfModules[["w001"]] <- split(linkList_001$Target, factor(linkList_001$TF))

# 2: Weight > 0.005
llminW <- linkList_001[which(linkList_001[,"weight"]>0.005),]
tfModules[["w005"]] <- split(llminW$Target, factor(llminW$TF))

# 3: Top 50 targets for each TF
# ("w001" should be ordered decreasingly by weight)
tfModules[["top50"]] <- lapply(tfModules[["w001"]], function(x) x[1:(min(length(x), 50))])

# 4-6: Top regulators per target
# (linkList_001 should be ordered by weight!)
linkList_001_byTarget <- split(linkList_001, factor(linkList_001$Target))

nTopTfs <- c(5, 10, 50)
nTopTfs <- setNames(nTopTfs, paste("top", nTopTfs, "perTarget", sep=""))

#library(reshape2); library(data.table)
topTFsperTarget <- lapply(linkList_001_byTarget, function(llt) {
  nTFs <- nTopTfs[which(nTopTfs <= nrow(llt))]
  reshape2::melt(lapply(nTFs, function(x) llt[1:x,"TF"]))
})


topTFsperTarget <- topTFsperTarget[which(!sapply(sapply(topTFsperTarget, nrow), is.null))]
topTFsperTarget.asDf <-  data.frame(data.table::rbindlist(topTFsperTarget, idcol=TRUE))


# topTFsperTarget.asDf <- apply(topTFsperTarget.asDf, 2, as.character)
colnames(topTFsperTarget.asDf) <- c("Target", "TF", "method")


# Merge the all the gene-sets:
tfModules.melted <- reshape2::melt(tfModules)
colnames(tfModules.melted) <- c("Target", "TF", "method")
tfModules <- rbind(tfModules.melted, topTFsperTarget.asDf)
rm(tfModules.melted); rm(topTFsperTarget.asDf)
tfModules$TF <- as.character(tfModules$TF)
tfModules$Target <- as.character(tfModules$Target)

# Basic counts:  #TODO add comment
if(getSettings(scenicOptions, "verbose")) 
  print(
    rbind(nTFs=length(unique(tfModules$TF)),
          nTargets=length(unique(tfModules$Target)),
          nGeneSets=nrow(unique(tfModules[,c("TF","method")])),
          nLinks=nrow(tfModules))
  )
# [,1]
# nTFs          789
# nTargets     8231
# nGeneSets    4734
# nLinks    3928451

### Add correlation to split into positive- and negative-correlated targets
corrMat <- loadInt(scenicOptions, "corrMat")
# Keep only correlation between TFs and potential targets
tfs <- unique(tfModules$TF)
missingTFs <- tfs[which(!tfs %in% rownames(corrMat))]
if(length(missingTFs) >0 ) 
{ 
  warning("The following TFs are missing from the correlation matrix: ", paste(missingTFs, collapse=", "))
  
  tfs <- tfs[which(tfs %in% rownames(corrMat))]
  corrMat <- corrMat[tfs,]
}

# Add correlation to the table
# "corr" column: 1 if the correlation between the TF and the target is > 0.03, -1 if the correlation is < -0.03 and 0 otherwise.(デフォルト)
#[TF, target]閾値の条件検討をする
corrMat[c("Sox9", "Foxp2", "Tbx20", "Twist1", "Scx", "Klf2","Klf4", "Osr1","Foxs1", "Plagl1"), 
        c("Col2a1", "Scx", "Fmod", "Myh11", "Dlk1", "Fgfr2", "Sox9")]
# Col2a1          Scx        Fmod      Myh11        Dlk1       Fgfr2        Sox9
# Sox9    0.2422023  0.212116953  0.01040008 -0.2009097 -0.10442442  0.28333575  1.00000000
# Foxp2   0.2573957  0.073208810 -0.15226220 -0.1325460 -0.13889556  0.27977901  0.24455341
# Tbx20   0.2053332  0.199206636  0.15153416 -0.1539279 -0.06314064  0.28040246  0.31374585
# Twist1  0.1373538  0.004240417 -0.21382936 -0.2505399 -0.15932055  0.25194801  0.14412630
# Scx     0.1099113  1.000000000  0.39662474  0.1368374  0.01775373  0.25117710  0.21211695
# Klf2   -0.2827657  0.161213581  0.45708977  0.5736602  0.46556598 -0.12957617 -0.20456438
# Klf4   -0.1177469  0.127171903  0.29712285  0.3017404  0.31363351  0.01972795 -0.04398116
# Osr1   -0.2378249  0.051769749  0.18497216  0.3785965  0.27167023 -0.14309025 -0.18120791
# Foxs1  -0.2069784  0.158972383  0.41702754  0.5413131  0.20974639 -0.05870651 -0.19920460
# Plagl1 -0.1028338 -0.018014087  0.14000513  0.2515350  0.31892630 -0.13477703 -0.04989005
#######

tfModules_byTF <- split(tfModules, as.factor(tfModules$TF))
tfModules_withCorr_byTF <- lapply(tfModules_byTF[tfs[1 : length(tfs)]], function(tfGeneSets)
{
  tf <- as.character(unique(tfGeneSets$TF))
  targets <- as.character(tfGeneSets$Target)
  cbind(tfGeneSets, corr=c(as.numeric(corrMat[tf,targets] > 0.03) - as.numeric(corrMat[tf,targets] < -0.03)))
})
tfModules_withCorr <- data.frame(data.table::rbindlist(tfModules_withCorr_byTF))
if(length(missingTFs) >0 ) 
{ 
  tfModules_withCorr <- rbind(tfModules_withCorr, data.frame(tfModules[tfModules$TF %in% missingTFs,], corr=NA))
}

tfModules_withCorr %>% dim
saveRDS(tfModules_withCorr, file=getIntName(scenicOptions, "tfModules_asDF"))

scenicOptions@settings$nCores <- 6
minGenes=20
coexMethod=NULL

#Load co-expression modules and databases:
nCores <- getSettings(scenicOptions, "nCores")
tfModules_asDF <- loadInt(scenicOptions, "tfModules_asDF")

#if(!is.null(coexMethod)) tfModules_asDF <-tfModules_asDF[which(tfModules_asDF$method %in% coexMethod),]
if(nrow(tfModules_asDF)==0) stop("The co-expression modules are empty.")

# Set cores for RcisTarget::addMotifAnnotation(). The other functions use foreach package.
if("BiocParallel" %in% installed.packages()) library(BiocParallel); register(MulticoreParam(nCores), default=TRUE) 

msg <- paste0(format(Sys.time(), "%H:%M"), "\tStep 2. Identifying regulons")
if(getSettings(scenicOptions, "verbose")) message(msg)

### Check org and load DBs
if(is.na(getDatasetInfo(scenicOptions, "org"))) stop('Please provide an organism (scenicOptions@inputDatasetInfo$org).')
library(AUCell)
library(RcisTarget)
motifAnnot <- getDbAnnotations(scenicOptions)


if(is.null(names(getSettings(scenicOptions, "dbs")))) 
{
  names(scenicOptions@settings$"dbs") <- scenicOptions@settings$"dbs"
  tmp <- sapply(strsplit(getSettings(scenicOptions, "dbs"),"-", fixed=T), function(x) x[grep("bp|kb",x)])
  if(all(lengths(tmp)>0)) names(scenicOptions@settings$"dbs") <- tmp
}

loadAttempt <- sapply(getDatabases(scenicOptions), dbLoadingAttempt)
if(any(!loadAttempt)) stop("It is not possible to load the following databses: \n",
                           paste(dbs[which(!loadAttempt)], collapse="\n"))

genesInDb <- unique(unlist(lapply(getDatabases(scenicOptions), function(x)
  names(feather::feather_metadata(x)[["types"]]))))

# Remove genes missing from RcisTarget databases
#  (In case the input matrix wasn't already filtered)
tfModules_asDF$TF <- as.character(tfModules_asDF$TF)
tfModules_asDF$Target <- as.character(tfModules_asDF$Target)
allTFs <- getDbTfs(scenicOptions)

tfModules_asDF <- tfModules_asDF[which(tfModules_asDF$TF %in% allTFs),]
geneInDb <- tfModules_asDF$Target %in% genesInDb
missingGene <- sort(unique(tfModules_asDF[which(!geneInDb),"Target"]))
if(length(missingGene)>0) 
  warning(paste0("Genes in co-expression modules not available in RcisTargetDatabases: ", 
                 paste(missingGene, collapse=", ")))
tfModules_asDF <- tfModules_asDF[which(geneInDb),]

# Targets with positive correlation
tfModules_Selected <- tfModules_asDF[which(tfModules_asDF$corr==1),]

# Add a column with the geneSet name (TF_method)
tfModules_Selected <- cbind(tfModules_Selected, geneSetName=paste(tfModules_Selected$TF, tfModules_Selected$method, sep="_"))
tfModules_Selected$geneSetName <- factor(as.character(tfModules_Selected$geneSetName))
# head(tfModules_Selected)
allGenes <- unique(tfModules_Selected$Target)
as.matrix(allGenes) %>% View

# Split into tfModules (TF-modules, with several methods)
tfModules <- split(tfModules_Selected$Target, tfModules_Selected$geneSetName)

# Add TF to the gene set (used in the following steps, careful if editing)
tfModules <- setNames(lapply(names(tfModules), function(gsn) {
  tf <- strsplit(gsn, "_")[[1]][1]
  unique(c(tf, tfModules[[gsn]]))
}), names(tfModules))

# Keep gene sets with at least 'minGenes' genes
tfModules <- tfModules[which(lengths(tfModules)>=minGenes)]
saveRDS(tfModules, file=getIntName(scenicOptions, "tfModules_forEnrichment")) 
#TODO as geneset? & previous step?

if(getSettings(scenicOptions, "verbose")) {
  tfModulesSummary <- t(sapply(strsplit(names(tfModules), "_"), function(x) x[1:2]))
  message("tfModulesSummary:")
  print(sort(table(tfModulesSummary[,2])))
}

# top5perTarget top10perTarget          top50           w005 top50perTarget           w001 
# 645            781            783            783            789            789 

####Motif enrichment analysis & identifying direct targets
#1. Calculate motif enrichment for each TF-module (Run RcisTarget)
#1.1 Calculate enrichment
msg <- paste0(format(Sys.time(), "%H:%M"), "\tRcisTarget: Calculating AUC")
if(getSettings(scenicOptions, "verbose")) message(msg)

motifs_AUC <- lapply(getDatabases(scenicOptions), function(rnkName) {
  ranking <- importRankings(rnkName, columns=allGenes)
  message("Scoring database: ", ranking@description)
  RcisTarget::calcAUC(tfModules, ranking, aucMaxRank=0.03*getNumColsInDB(ranking), nCores=nCores, verbose=FALSE)})

saveRDS(motifs_AUC, file=getIntName(scenicOptions, "motifs_AUC"))

#1.2 Convert to table, filter by NES & add the TFs to which the motif is annotated
# motifs_AUC <- loadInt(scenicOptions, "motifs_AUC") # to start from here
# (For each database...)

msg <- paste0(format(Sys.time(), "%H:%M"), "\tRcisTarget: Adding motif annotation")
message(msg)
motifEnrichment <- lapply(motifs_AUC, function(aucOutput)
{
  # Extract the TF of the gene-set name (i.e. MITF_w001):
  tf <- sapply(setNames(strsplit(rownames(aucOutput), "_"), rownames(aucOutput)), function(x) x[[1]])
  
  # Calculate NES and add motif annotation (provide tf in 'highlightTFs'):
  addMotifAnnotation(aucOutput, 
                     nesThreshold=3, digits=3, 
                     motifAnnot=motifAnnot,
                     motifAnnot_highConfCat=c("directAnnotation", "inferredBy_Orthology"),
                     motifAnnot_lowConfCat=c("inferredBy_MotifSimilarity",
                                             "inferredBy_MotifSimilarity_n_Orthology"), 
                     highlightTFs=tf)
})
# Merge both tables, adding a column that contains the 'motifDb'
motifEnrichment <- do.call(rbind, lapply(names(motifEnrichment), function(dbName){
  cbind(motifDb=dbName, motifEnrichment[[dbName]])
}))
saveRDS(motifEnrichment, file=getIntName(scenicOptions, "motifEnrichment_full"))
msg <- paste0("Number of motifs in the initial enrichment: ", nrow(motifEnrichment))
if(getSettings(scenicOptions, "verbose")) message(msg)
#motifEnrichmentの中にはScxデータが
#subset(motifEnrichment, motifEnrichment$highlightedTFs %in% "Scx") %>% headとしてある。

#1.3 Keep only the motifs annotated to the initial TF
# motifEnrichment <- loadInt(scenicOptions, "motifEnrichment_full")
motifEnrichment_selfMotifs <- motifEnrichment[which(motifEnrichment$TFinDB != ""),, drop=FALSE]
dim(motifEnrichment)
#[1] 2119013       9
dim(motifEnrichment_selfMotifs)
#[1] 17496     9

msg <- paste0("Number of motifs annotated to the corresponding TF: ", nrow(motifEnrichment))
if(getSettings(scenicOptions, "verbose")) message(msg)
#rm(motifEnrichment)
#Number of motifs annotated to the corresponding TF: 2119013

if(nrow(motifEnrichment)==0) 
  stop("None of the co-expression modules present enrichment of the TF motif: There are no regulons.")

#2. Prune targets
msg <- paste0(format(Sys.time(), "%H:%M"), "\tRcisTarget: Prunning targets")
if(getSettings(scenicOptions, "verbose")) message(msg)

dbNames <- getDatabases(scenicOptions)
motifEnrichment_wGenes <- lapply(names(dbNames), function(motifDbName){
  ranking <- importRankings(dbNames[motifDbName], columns=allGenes)
  addSignificantGenes(resultsTable=motifEnrichment[motifEnrichment$motifDb==motifDbName,],
                      geneSets=tfModules,
                      rankings=ranking,
                      maxRank=5000, method="aprox", nCores=nCores)
})
#17:34- 35GB
suppressPackageStartupMessages(library(data.table))
motifEnrichment_wGenes <- rbindlist(motifEnrichment_wGenes)
#保存する名前がmotifEnrichment_wGenesからmotifEnrichment_selfMotifs_wGenesに変更したので注意
saveRDS(motifEnrichment_wGenes, file=getIntName(scenicOptions, "motifEnrichment_selfMotifs_wGenes"))
#iMac O/N




if(getSettings(scenicOptions, "verbose")) 
{
  # TODO messages/print
  message(format(Sys.time(), "%H:%M"), "\tNumber of motifs that support the regulons: ", nrow(motifEnrichment_wGenes))
  motifEnrichment_wGenes[order(motifEnrichment_wGenes$NES,decreasing=TRUE),][1:5,(1:ncol(motifEnrichment_wGenes)-1), with=F] 
}
#Save motif enrichment results as text and HTML (optional):
# motifEnrichment_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes") # to start from here

# Text:
if(!file.exists("output")) dir.create("output") 
write.table(motifEnrichment_wGenes, file=getOutName(scenicOptions, "s2_motifEnrichment"),
            sep="\t", quote=FALSE, row.names=FALSE)

# HTML
if("DT" %in% installed.packages() && nrow(motifEnrichment_wGenes)>0)
{
  nvm <- tryCatch({
    colsToShow <- c("motifDb", "logo", "NES", "geneSet", "TF_highConf", "TF_lowConf")
    motifEnrichment_2html <- viewMotifs(motifEnrichment_wGenes, colsToShow=colsToShow, options=list(pageLength=100))
    
    fileName <- getOutName(scenicOptions, "s2_motifEnrichmentHtml")
    
    dirName <- dirname(fileName)
    fileName <- basename(fileName)
    suppressWarnings(DT::saveWidget(motifEnrichment_2html, fileName))
    file.rename(fileName, file.path(dirName, fileName))
    if(getSettings(scenicOptions, "verbose")) message("Preview of motif enrichment saved as: ", file.path(dirName, fileName))
  }, error = function(e) print(e$message))
}



#Format regulons & save
motifEnrichment.asIncidList <- apply(motifEnrichment_wGenes, 1, function(oneMotifRow) {
  genes <- strsplit(oneMotifRow["enrichedGenes"], ";")[[1]]
  oneMotifRow <- data.frame(rbind(oneMotifRow), stringsAsFactors=FALSE)
  data.frame(oneMotifRow[rep(1, length(genes)),c("NES", "motif", "highlightedTFs", "TFinDB")], genes, stringsAsFactors = FALSE)
})
motifEnrichment.asIncidList <- data.table::rbindlist(motifEnrichment.asIncidList)
colnames(motifEnrichment.asIncidList) <- c("NES", "motif", "TF", "annot", "gene")
motifEnrichment.asIncidList <- data.frame(motifEnrichment.asIncidList, stringsAsFactors = FALSE)
#subset(motifEnrichment.asIncidList, gene %in% "Scx") %>% head
#subset(motifEnrichment.asIncidList, TF %in% "Scx") %>% View
#Scxが出てくる

# Get targets for each TF, but keep info about best motif/enrichment
# (directly annotated motifs are considered better)
regulonTargetsInfo <- lapply(split(motifEnrichment.asIncidList, motifEnrichment.asIncidList$TF), function(tfTargets){
  # print(unique(tfTargets$TF))
  tfTable <- as.data.frame(do.call(rbind, lapply(split(tfTargets, tfTargets$gene), function(enrOneGene){
    highConfAnnot <- "**" %in% enrOneGene$annot
    enrOneGeneByAnnot <- enrOneGene
    if(highConfAnnot) enrOneGeneByAnnot <- enrOneGeneByAnnot[which(enrOneGene$annot == "**"),]
    bestMotif <- which.max(enrOneGeneByAnnot$NES)
    
    cbind(TF=unique(enrOneGene$TF), gene=unique(enrOneGene$gene), nMotifs=nrow(enrOneGene),
          bestMotif=as.character(enrOneGeneByAnnot[bestMotif,"motif"]), NES=as.numeric(enrOneGeneByAnnot[bestMotif,"NES"]),
          highConfAnnot=highConfAnnot)
  })), stringsAsFactors=FALSE)
  tfTable[order(tfTable$NES, decreasing = TRUE),]
})
rm(motifEnrichment.asIncidList)
regulonTargetsInfo <- data.table::rbindlist(regulonTargetsInfo)
colnames(regulonTargetsInfo) <- c("TF", "gene", "nMotifs", "bestMotif", "NES", "highConfAnnot")
#subset(regulonTargetsInfo, TF  %in% "Scx")
#Scxあり

#Optional: Add GENIE3 score to export
linkList <- loadInt(scenicOptions, "genie3ll", ifNotExists="null")
if(!is.null(linkList) & ("weight" %in% colnames(linkList)))
{
  if(is.data.table(linkList)) linkList <- as.data.frame(linkList)
  
  uniquePairs <- nrow(unique(linkList[,c("TF", "Target")]))
  if(uniquePairs == nrow(linkList)) {
    linkList <- linkList[which(linkList$weight>=getSettings(scenicOptions, "modules/weightThreshold")),]  # TODO: Will not work with GRNBOOST!
    rownames(linkList) <- paste(linkList$TF, linkList$Target,sep="__")
    regulonTargetsInfo <- cbind(regulonTargetsInfo, Genie3Weight=linkList[paste(regulonTargetsInfo$TF, regulonTargetsInfo$gene,sep="__"),"weight"])
  }else {
    warning("There are duplicated regulator-target (gene id/name) pairs in the co-expression link list.",
            "\nThe co-expression weight was not added to the regulonTargetsInfo table.")
  }
}else warning("It was not possible to add the weight to the regulonTargetsInfo table.")

saveRDS(regulonTargetsInfo, file=getIntName(scenicOptions, "regulonTargetsInfo"))

write.table(regulonTargetsInfo, file=getOutName(scenicOptions, "s2_regulonTargetsInfo"),
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
rm(linkList)

#Split into regulons according to the motif annotation
regulonTargetsInfo_splitByAnnot <- split(regulonTargetsInfo, regulonTargetsInfo$highConfAnnot)
regulons <- NULL
if(!is.null(regulonTargetsInfo_splitByAnnot[["TRUE"]]))
{
  regulons <- lapply(split(regulonTargetsInfo_splitByAnnot[["TRUE"]], regulonTargetsInfo_splitByAnnot[["TRUE"]][,"TF"]), function(x) sort(as.character(unlist(x[,"gene"]))))
}
regulons_extended <- NULL
if(!is.null(regulonTargetsInfo_splitByAnnot[["FALSE"]]))
{
  regulons_extended <- lapply(split(regulonTargetsInfo_splitByAnnot[["FALSE"]],regulonTargetsInfo_splitByAnnot[["FALSE"]][,"TF"]), function(x) unname(unlist(x[,"gene"])))
  regulons_extended <- setNames(lapply(names(regulons_extended), function(tf) sort(unique(c(regulons[[tf]], unlist(regulons_extended[[tf]]))))), names(regulons_extended))
  names(regulons_extended) <- paste(names(regulons_extended), "_extended", sep="")
}
regulons <- c(regulons, regulons_extended)
saveRDS(regulons, file=getIntName(scenicOptions, "regulons"))




