

library(tidyverse)
library(igraph)
library(linkcomm)
library(pheatmap)


setwd("~/Downloads/scRNA-seq/20200403_Seuratv3_based2/SiGN-BN/200514_test/4thv2")
net1 <- read.table("result.txt", header =T, row.names=NULL)

regTar <- read.table("regTar16v2.txt", header =F, row.names =NULL)
regTar %>% distinct(V1, .keep_all = F) %>% dim

regTar %>% distinct(V1, .keep_all = F) -> TFs
TFs$V1 %>% as.character() -> TFs
subset(net1, Child %in% TFs) %>% dim

subset(net1, Child %in% TFs) -> TFsnet
TFsnet <- TFsnet[c(1,2)]

a <- as.character(TFsnet$Parent) %>% as.data.frame()
b <- as.character(TFsnet$Child) %>% as.data.frame()
c <- rbind(a, b)
colnames(c) <- "gene"
c %>% distinct(gene) %>% dim





g <- simplify(graph.data.frame(TFsnet,directed=T),
              remove.multiple=T,remove.loops=T)

V(g)$name <- paste(V(g)$name,V(g)$Faction,sep=":")

el <- get.edgelist(g)

el <- gsub(" ","_",el)



linkcomm <- getLinkCommunities(el, hcmethod="ward.D2", use.all.edges=T,
                               directed = T)

plot(linkcomm, type="graph")
text(0.95*par("usr")[1],0.95*par("usr")[4],label="LinkComm",pos=4)

plot(linkcomm, type="members")
text(0.95*par("usr")[1],0.95*par("usr")[4],label="LinkComm",pos=4)

linkcomm$pdmax 

linkcomm_at <- newLinkCommsAt(linkcomm, cutat=linkcomm$pdmax)
linkcomm_at

plot(linkcomm_at, type="graph")
text(0.95*par("usr")[1],0.95*par("usr")[4],label="LinkComm",pos=4)


setwd("~/Downloads/scRNA-seq/20200403_Seuratv3_based2/SiGN-BN/200514_test/4thv2/modularity_overlapward2")
colour_code <- read.delim("colour_code.txt", header =T, row.names = NULL)
colnames(colour_code) <- c("code", "colour")
colour_code$code[1:109]

plotLinkCommGraph(linkcomm_at, pal = colour_code$code[1:109], 
                  vlabel.cex = "none")
plotLinkCommSumm(linkcomm_at, pal = colour_code$code[1:109])


linkcomm_at$nodeclusters %>% View
linkcomm_at$nodeclusters %>% subset(cluster == 78)
linkcomm_at$nodeclusters %>% subset(cluster == 99) -> comm

########################################################################################
#save.image("~/Downloads/scRNA-seq/20200403_Seuratv3_based2/SiGN-BN/200514_test/4thv2/modularity_overlapward2/G1.RData")

########################################################################################





