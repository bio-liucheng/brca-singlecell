library(pheatmap)
library(dplyr)
library(NNLM)

limitData <- function(data, min = NULL, max = NULL){
  data2 <- data
  if(!is.null(min)){
    data2[data2 < min] <- min
  }
  if(!is.null(max)){
    data2[data2 > max] <- max
  }
  return(data2)
}
patient_list = paste("patient", c(1,2,4,5,6,7,8), sep = '_')

for(i in 1:length(patient_list)){
  
seurat_object_edgeR <- subset(seurat.obj.md,  patient_id %in% patient_list[i])

seurat_object_edgeR <- NormalizeData(seurat_object_edgeR)
seurat_object_edgeR <- FindVariableFeatures(seurat_object_edgeR, nfeatures = 5000)
seurat_object_edgeR <- ScaleData(seurat_object_edgeR)

seurat_object_edgeR <- RunPCA(seurat_object_edgeR)
seurat_object_edgeR <- FindNeighbors(seurat_object_edgeR, reduction = "pca", dims = 1:20)
seurat_object_edgeR <- FindClusters(seurat_object_edgeR, resolution = 0.8)



data <- as(object = seurat_object_edgeR[["RNA"]]@data, Class = "dgTMatrix")
#data <- seurat.obj.md[["RNA"]]@scale.data
ave.data <-  Matrix::rowSums(data) / Matrix::rowSums(data > 0)

data@x <- data@x - ave.data[data@i + 1]


data@x[data@x <0] <- 0

patient_nnmf <- nnmf(as.matrix(data), 100, verbose = 0)

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


savepath = patient_list[i]
savepath = "all_100"
if(!dir.exists(savepath)){
  dir.create(savepath)
}

meta_s <- seurat_object_edgeR[[]]
write.csv(W, file.path(savepath, "W-gene-program.csv"))
write.csv(H, file.path(savepath, "H-gene-program.csv"))
write.csv(program.gene.value, file.path(savepath, "program.gene.value.csv"))
write.csv(meta_s,  file.path(savepath, "meta_patient.csv"))

t <- cor(W[,programe_sd>0.01])
t[t>0.8] <- 0.8
t[t<0.1] <- 0
p <- pheatmap(t)

pdf(file.path(savepath, paste0("programe", savepath, ".pdf")), width = 8, height = 8)
grid::grid.newpage()
grid::grid.draw(p$gtable)
dev.off()

}

`up.bound <- quantile(as.matrix(s$H), 0.995)
H <- limitData(s$H, max = up.bound)
meta_s <- seurat_object_edgeR[[]]
meta_s$cell_name <- rownames(meta_s)
meta_s <- meta_s %>% arrange(seurat_clusters)
cluster <- data.frame(cluster = meta_s$seurat_clusters, 
                      sample_type = meta_s$sample_type,
                      patient_id = meta_s$patient_id)
rownames(cluster) <- meta_s$cell_name

t_all <- c()
for(i in 1:length(t)){
  all <- sum(t[1:i])
  t_all <- c(t_all, all)
}
pheatmap(H[,meta_s$cell_name], show_colnames = F,cluster_cols = F, annotation_col = cluster, gas_col = t_all[-length(t_all)])
