setwd("/Users/leizhang/Desktop/Bioinformatics/cDC_processing")

library(Seurat)
library(dplyr)
library(monocle)
options(stringsAsFactors=FALSE)
library(reticulate)
library(SingleR)
library(SingleCellExperiment)
library(scran)

cDC <- readRDS("seurat.obj.md.RDS")

cluster <- as.vector(unique(cDC@meta.data[["cell_label"]]))
train_clusters = c(cluster[grep("DC-C1-CD1C",cluster)], cluster[grep("DC-C2-CLEC9A",cluster)]) ##c('M03_cDC1_CLEC9A','M04_cDC2_CD1C')
test_cluster = cluster[grep("DC-C3-LAMP3",cluster)]  ##c('M05_cDC3_LAMP3')


trainMatrix <- cDC@assays[["RNA"]][,rownames(cDC@meta.data[cDC@meta.data[["cell_label"]] %in% train_clusters,])]
testMatrix <- cDC@assays[["RNA"]][,rownames(cDC@meta.data[cDC@meta.data[["cell_label"]] %in% test_cluster,])]

train_label = as.vector(cDC@meta.data[cDC@meta.data[["cell_label"]] %in% train_clusters,]["cell_label"])
train_labels <-  t(train_label)
out <- pairwiseTTests(trainMatrix, train_labels, direction = "up")
markers <- getTopMarkers(out$statistics, out$pairs, n=100)
trained <- trainSingleR(trainMatrix, labels=train_labels, genes=markers)
pred2b <- classifySingleR(testMatrix, trained)
pred_table <- data.frame(SingleR_Label = pred2b$labels)
rownames(pred_table) <- pred2b@rownames

source_test <- as.vector(cDC@meta.data[cDC@meta.data[["cell_label"]] %in% test_cluster,]["sample_type"])
pred_table <- cbind(pred_table,source_test)
Her2statu <- as.vector(cDC@meta.data[cDC@meta.data[["cell_label"]] %in% test_cluster,]["Her2_statu"])
ERstatu <- as.vector(cDC@meta.data[cDC@meta.data[["cell_label"]] %in% test_cluster,]["ER_statu"])
pred_table <- cbind(pred_table,ERstatu,Her2statu)

table(pred_table$source)
table(pred_table$SingleR_Label)
write.csv(pred_table, "pred_table.csv")

pre_Tumor <- pred_table[pred_table$sample_type ==  "Tumor",]
pre_Lymph <- pred_table[pred_table$sample_type ==  "Lymph Node",]
pre_T <- data.frame(table(pre_Tumor$SingleR_Label)) 
pre_T$source <- c("T","T")
x <- as.numeric(sum(pre_T$Freq)) 
pre_T$Freqs <- pre_T$Freq/x
pre_L <- data.frame(table(pre_Lymph$SingleR_Label)) 
pre_L$source <- c("L","L")
x <- as.numeric(sum(pre_L$Freq)) 
pre_L$Freqs <- pre_L$Freq/x
pre_test <- rbind(pre_T,pre_L)

preHerposiv <- pred_table[pred_table$Her2_statu == "positive",]
preHerneg <- pred_table[pred_table$Her2_statu == "Negative",]
preHerP <- data.frame(table(preHerposiv$SingleR_Label))
preHerP$Her2Statu <- c("Pos","Pos")
x <- as.numeric(sum(preHerP$Freq)) 
preHerP$Freqs <- preHerP$Freq/x
preHerN <- data.frame(table(preHerneg$SingleR_Label))
preHerN$Her2Statu <- c("Neg","Neg")
x <- as.numeric(sum(preHerN$Freq)) 
preHerN$Freqs <- preHerN$Freq/x
pre_HEr2 <- rbind(preHerP,preHerN)

preHerposiv <- pred_table[pred_table$Her2_statu == "positive",]
preHpT <- preHerposiv[preHerposiv$sample_type == "Tumor",]
preHpL <- preHerposiv[preHerposiv$sample_type == "Lymph Node",]
preHT <- data.frame(table(preHpT$SingleR_Label))
preHT$source <- c("T","T")
x <- as.numeric(sum(preHT$Freq)) 
preHT$Freqs <- preHT$Freq/x
preHL <- data.frame(table(preHerposiv$SingleR_Label))
preHL$source <- c("L","L")
x <- as.numeric(sum(preHL$Freq)) 
preHL$Freqs <- preHL$Freq/x
preHp <- rbind(preHT,preHL)


p1 <- ggplot(pre_test, aes(x = source ,y = Freqs, fill = Var1))+
  ####position="stack"¶Ñµþ×´
  geom_bar(stat ="identity",width = 0.6,position ="stack")+
  labs( title = "tumorVslymph")

p2 <- ggplot(pre_HEr2, aes(x = Her2Statu ,y = Freqs, fill = Var1))+
  ####position="stack"¶Ñµþ×´
  geom_bar(stat ="identity",width = 0.6,position ="stack")+
  labs( title = "Her2_statu")  

p3 <- ggplot(preHp, aes(x = source ,y = Freqs, fill = Var1))+
  ####position="stack"¶Ñµþ×´
  geom_bar(stat ="identity",width = 0.6,position ="stack")+
  labs( title = "Her2+_tumorVslymph")
pdf("cDCLAMP3.pdf",9,9)
print(p1+p2+p3)
dev.off()

