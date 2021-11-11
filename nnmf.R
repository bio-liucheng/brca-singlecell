#my nnmf
library(Seurat)
library(NNLM)
library(dplyr)
library(pheatmap)
options(stringsAsFactors = F)
meta <-read.csv("meta.csv")
patient_id = "patient_2"
meta_patient <- meta[meta$patient_id == patient_id,]
seurat.obj.md <- seurat.obj[,meta_patient$X]
seurat.obj.md <- NormalizeData(seurat.obj.md)
seurat.obj.md <- FindVariableFeatures(seurat.obj.md, nfeatures = 3000)
seurat.obj.md <- ScaleData(seurat.obj.md, features = rownames(seurat.obj.md))

data <- as(object = seurat.obj.md[["RNA"]]@data, Class = "dgTMatrix")
#data <- seurat.obj.md[["RNA"]]@scale.data
ave.data <-  Matrix::rowSums(data) / Matrix::rowSums(data > 0)

data@x <- data@x - ave.data[data@i + 1]


data@x[data@x <0] <- 0

patient_nnmf <- nnmf(as.matrix(data), 60, verbose = 0)

W <- patient_nnmf$W
colnames(W) <- paste0("p", 1:dim(W)[2])
H <- patient_nnmf$H
rownames(H) <- paste0("p", 1:dim(H)[1])

all.genes <- rownames(W)
sel.W <- (W > quantile(W, 1 - 50/dim(W)[1]))
for(pi in 1:dim(W)[2]){
  if(sum(sel.W[, pi]) > 10){
    tmp <- data.frame(program = colnames(W)[pi], gene = all.genes[sel.W[, pi]], value = W[sel.W[, pi], pi])
    tmp <- tmp[order(tmp$value, decreasing = T), ]
  }else{
    tmp <- data.frame(program = colnames(W)[pi], gene = all.genes, value = W[, pi])
    tmp <- tmp[order(tmp$value, decreasing = T), ]
    tmp <- tmp[1:10, ]
  }
  if(pi == 1){
    program.gene.value <- tmp
  }else{
    program.gene.value <- rbind(program.gene.value, tmp)
  }
}


up.bound <- quantile(as.matrix(H), 0.995)
H <- limitData(H, max = up.bound)

savepath = "patient_2"
meta_patient <- read.csv(file.path(savepath, "meta_patient.csv"))

meta_patient <- meta_patient %>% arrange(sample_type, seurat_clusters)
cluster <- data.frame(cluster =paste0("cluster_", meta_patient$seurat_clusters), sample_type = meta_patient$sample_type)
rownames(cluster) <- meta_patient$X
pheatmap(H[,meta_patient$X], show_colnames = F, cluster_cols = T,cluster_rows = T, annotation_col = cluster)
s <- apply(H, 1, sd)
s[order(s)]
programs = names(s)[s >0.01]


program = c("p2", "p4")


write.csv(W[, colnames(W) %in% program], file.path(savepath, "W-gene-program-manual.csv"))
write.csv(H[rownames(H) %in% program, ], file.path(savepath, "H-gene-program-manual.csv"))
write.csv(program.gene.value[program.gene.value$program %in% program,], file.path(savepath, "program.gene.value.manual.csv"))
write.csv(meta_patient,  file.path(savepath, "meta_patient.csv"))

n= 10

top_20 <- program.gene.value %>% filter(program %in% programs)  #%>% 
top_20 <- top_20 %>% mutate(rank = 1:nrow(top_20)) %>% group_by(program) %>% top_n(-n, rank)

gene = unique(top_20$gene)
pheatmap(data[gene,meta_patient$X], show_colnames = F,cluster_cols = T,cluster_rows = T, annotation_col = cluster)

gene = program.gene$gene[program.gene$program %in% program]



gene = gene[gene %in% rownames(matrix)]

pheatmap(matrix[gene,meta_patient$X], show_colnames = F,cluster_cols = T, cluster_rows = T, annotation_col = cluster)
