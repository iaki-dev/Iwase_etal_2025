###ランダムフォレストによる分類
# データの分割方法を変える
# caretではなく、ランダムフォレストする
# 統合データでcommunityを計算する

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
#並列化演算による高速化パッケージインストール
# install.packages("doParallel")
library(doParallel)
#コア数を調べおく
detectCores()
#[1] 10
# cl <- makePSOCKcluster(16)
# registerDoParallel(cl)
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





#データの読み込み
#これは、細胞内でどのコミュニティをscale化したもの
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/OC")
gf <- read_tsv("comm_cell_mix_scaled.txt")
dim(gf)
#[1] 109 32983
tmp.names <- paste("OC", 1:nrow(gf), sep = "")
gf <- as.data.frame(gf)
rownames(gf) <- tmp.names
gf <- gf[-c(1)]
gf <- t(gf) %>% as.data.frame()
#clusterラベル用データ
gf <- cbind(NCC_merge@meta.data["integrated"], gf)
dim(gf)
#[1] 32982   110
colnames(gf)[1] <- "cluster"
gf$cluster <- as.character(gf$cluster)
gf$cluster <- as.factor(gf$cluster)

#訓練データとテストデータに分割
set.seed(423)
inTrain <- createDataPartition(y = gf$cluster, p = 0.7, list = F)
gf.train <- gf[inTrain,]
gf.test <- gf[-inTrain,]
cat("train=", nrow(gf.train), " test=", nrow(gf.test))
#train= 23093  test= 9889

gf.train$cluster %>% table
gf.test$cluster %>% table


# ランダムフォレスト（RandomForestパッケージ）
set.seed(0)
#データセットを特徴量とラベルに分割する
features <- gf.train[2:ncol(gf.train)]
labels <- gf.train[1]
labels <- as.factor( labels[[1]] )
#モデルを作成する
modelDrive <- randomForest( x=features, y=labels, importance=T, proximity=T )
#説明変数の重要度を表示する
print( importance( modelDrive ))
#説明変数(分類に寄与した変数)の重要度をプロット
varImpPlot( modelDrive )
# 学習の収束状況をプロットする
plot( modelDrive )
# グリッドサーチを行う
set.seed(423)
rfTuning <- tuneRF( x = features,y = labels,   ntreeTry=500,
                    stepFactor = 2, improve = 0.05, trace = TRUE, plot = TRUE, doBest = T )
modelDrive <- randomForest( x = features, y = labels, mtry = 10,
                            ntree = 500, importance = T )
# テストデータにモデルを適用する
prdct <- predict( modelDrive, newdata = gf.test )
# 分類の正確さの確認
correctAns <- 0
for ( i in 1:nrow( table( prdct,  gf.test$cluster )))
  correctAns <- correctAns + table(prdct, gf.test$cluster)[i,i]
correctAns / nrow( gf.test )

# これでもいける
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
#ランダムフォレスト予測
predgfRF <- predict(gfRF, gf.test)
#重要度
varImp(gfRF,scale = T) %>% plot()
# 16x6で保存
#評価
confusionMatrix(data = predgfRF, gf.test$cluster)

library(rattle)
fancyRpartPlot(gfRF$finalModel)


model <- randomForest(Species ~ ., data=iris, 
                      importance=TRUE, ntree=500, mtry = 2, do.trace=100)
reprtree:::plot.getTree(model)


model <- randomForest(cluster ~ ., data=gf.train, 
                      importance=TRUE, ntree=500, mtry = 2, do.trace=100)



reprtree:::plot.getTree(gfRF$finalModel)

saveRDS(gfRF, "gfRF.rds")
saveRDS(predgfRF, "predgfRF.rds")





#クラスタを調べる
# linkcomm_at$nodeclusters %>% subset(cluster == 91)
#ここまで実施

#confusinmatrixはconfusionMatrix(data = predgfRF, gf.test$cluster)$tableにある
#これをggplotでヒートマップにする
confusionMatrix(data = predgfRF, gf.test$cluster)$table %>% 
  as.data.frame() -> table1
ggplot(table1, aes(Prediction, Reference, group=Reference))	+
  geom_tile(aes(fill = Freq)) +
  geom_text(aes(fill = table1$Freq, label = table1$Freq))+
  scale_fill_gradient(low = "white", high = "red") 


setwd("~/Downloads/scRNA-seq/20200403_Seuratv3_based2/SiGN-BN/200514_test/4thv2/modularity_overlapward2/RFs_v8")
importance <- read.table("importance.txt")

res <- matrix(nrow = 8, ncol = 1)
for(i in 0:nrow(importance)){
  coi <- importance$V1[1:i] %>% as.character()
  coi <- c("cluster", coi)
  tmp.train <- gf.train[,coi]
  tmp.test <- gf.test[,coi]
  
  # ランダムフォレスト（RandomForestパッケージ）
  set.seed(0)
  rm(tmpRF)
  tmpRF <- train(
    cluster ~ (.), 
    data = tmp.train, 
    method = "rf", 
    tuneLength = 4, 
    preProcess = c('center', 'scale'),
    trControl = trainControl(method = "cv")
  )
  
  #ランダムフォレスト予測
  rm(pred.tmp)
  pred.tmp <- predict(tmpRF, tmp.test)
  x <- confusionMatrix(data = pred.tmp, tmp.test$cluster)
  tmp <- x$overall %>% as.data.frame()
  top <- i %>% data.frame()
  rownames(top) <- "top"
  tmp <- rbind(tmp, top)
  
  res <- cbind(res,  tmp)
  # ここにAccuracyのデータが存在している（再ランしても51% Accuracy）
  
}
# 09:40-10:01

res2 <- res %>% t() %>% as.data.frame()
res2 <- res2[-1,]
res2 <- res2[-1,]
res2 %>% arrange(Accuracy) %>% View
# これだとmaxで55.5%
# 保存
setwd("~/Downloads/scRNA-seq/20200403_Seuratv3_based2/SiGN-BN/200514_test/4thv2/modularity_overlapward2/RFs_v8")
write.table(res2, "res2.txt", quote = F, sep = "\t", row.names = F)

res2 <- read.table("res2.txt", header = T)
res2 %>% arrange(Accuracy) %>% View


#confusinmatrixはconfusionMatrix(data = predgfRF, gf.test$cluster)$tableにある
#これをggplotでヒートマップにする
confusionMatrix(data = pred.tmp, tmp.test$cluster)$table %>% 
  as.data.frame() -> table1
ggplot(table1, aes(Prediction, Reference, group=Reference))	+
  geom_tile(aes(fill = Freq)) +
  geom_text(aes(fill = table1$Freq, label = table1$Freq))+
  scale_fill_gradient(low = "white", high = "red") 
# 4x6で保存

gf.train$cluster %>% table
gf.test$cluster %>% table


