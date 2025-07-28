

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
library(caret)
library(tidyverse)
library(doParallel)
detectCores()
#[1] 10
library(randomForest)

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



setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/OC")
gf <- read_tsv("comm_cell_mix_scaled.txt")
dim(gf)
#[1] 109 32983
tmp.names <- paste("OC", 1:nrow(gf), sep = "")
gf <- as.data.frame(gf)
rownames(gf) <- tmp.names
gf <- gf[-c(1)]
gf <- t(gf) %>% as.data.frame()

gf <- cbind(NCC_merge@meta.data["integrated"], gf)
dim(gf)
#[1] 32982   110
colnames(gf)[1] <- "cluster"
gf$cluster <- as.character(gf$cluster)
gf$cluster <- as.factor(gf$cluster)

# Divide into training and test
set.seed(423)
inTrain <- createDataPartition(y = gf$cluster, p = 0.7, list = F)
gf.train <- gf[inTrain,]
gf.test <- gf[-inTrain,]
cat("train=", nrow(gf.train), " test=", nrow(gf.test))
#train= 23093  test= 9889

gf.train$cluster %>% table
gf.test$cluster %>% table



set.seed(0)
features <- gf.train[2:ncol(gf.train)]
labels <- gf.train[1]
labels <- as.factor( labels[[1]] )

modelDrive <- randomForest( x=features, y=labels, importance=T, proximity=T )

print( importance( modelDrive ))

varImpPlot( modelDrive )

plot( modelDrive )

set.seed(423)
rfTuning <- tuneRF( x = features,y = labels,   ntreeTry=500,
                    stepFactor = 2, improve = 0.05, trace = TRUE, plot = TRUE, doBest = T )
modelDrive <- randomForest( x = features, y = labels, mtry = 10,
                            ntree = 500, importance = T )

prdct <- predict( modelDrive, newdata = gf.test )

correctAns <- 0
for ( i in 1:nrow( table( prdct,  gf.test$cluster )))
  correctAns <- correctAns + table(prdct, gf.test$cluster)[i,i]
correctAns / nrow( gf.test )


model <- randomForest(cluster ~ ., data=gf.train, 
                      importance=TRUE, ntree=500, mtry = 10, do.trace=100)
model
reprtree:::plot.getTree(model)

saveRDS(model, "model.rds")


gfRF <- train(
  cluster ~ (.), 
  data = gf.train, 
  method = "rf", 
  tuneLength = 4, 
  preProcess = c('center', 'scale'),
  trControl = trainControl(method = "cv")
)

predgfRF <- predict(gfRF, gf.test)

varImp(gfRF,scale = T) %>% plot()


confusionMatrix(data = predgfRF, gf.test$cluster)

saveRDS(gfRF, "gfRF.rds")
saveRDS(predgfRF, "predgfRF.rds")


#これをggplotでヒートマップにする
confusionMatrix(data = predgfRF, gf.test$cluster)$table %>% 
  as.data.frame() -> table1
ggplot(table1, aes(Prediction, Reference, group=Reference))	+
  geom_tile(aes(fill = Freq)) +
  geom_text(aes(fill = table1$Freq, label = table1$Freq))+
  scale_fill_gradient(low = "white", high = "red") 


