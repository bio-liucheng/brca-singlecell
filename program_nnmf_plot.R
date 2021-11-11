library(pheatmap)
library(dplyr)
library(ComplexHeatmap)
options(stringsAsFactors = F)
savepath = "."
write.csv(patient_nnmf$W, file.path(savepath, "W-gene-program.csv"))
write.csv(patient_nnmf$H, file.path(savepath, "H-gene-program.csv"))
write.csv(program.gene.value, file.path(savepath, "program.gene.value.csv"))
write.csv(meta,  file.path(savepath, "meta_patient.csv"))

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
savepath = "all_100"

H <- read.csv(file.path(savepath, "H-gene-program.csv"))

rownames(H) <- H$X
H <- H[,-1]

W <- read.csv(file.path(savepath, "W-gene-program.csv"))
rownames(W) <- W$X
W <- W[,-1]

meta <- read.csv(file.path(savepath, "meta_patient.csv"))
program.gene <- read.csv(file.path(savepath, "program.gene.value.csv"))
up.bound <- quantile(as.matrix(H), 0.995)
H <- limitData(H, max = up.bound)

meta <- meta %>% arrange(seurat_clusters)
cluster <- data.frame(cluster = paste0("cluster_",meta$seurat_clusters), sample_type = meta$sample_type)
rownames(cluster) <- meta$X
t <- cor(t(H))
pheatmap(t)


programe_sd <- apply(H, 1, sd)

t <- cor(W[,programe_sd>0.01])
pheatmap(t)


pheatmap(H[,meta$X], show_colnames = F,cluster_cols = F,cluster_rows = F, annotation_col = cluster)
#Heatmap(H[,meta$cell_name], show_column_names = F, show_row_names = F)

average_cluster <- function(x, cluster){
  s <- tapply(x, cluster, mean)
}
t <- apply(H[,meta$X], 1, tapply, INDEX = meta$seurat_clusters, FUN = mean)
s <- apply(H[,meta$X], 1, average_cluster, cluster = meta$seurat_clusters)

s <- apply(s, 2, function(x) x/sum(x) )
#entropy
e <- apply(s, 2, function(x) -sum(x/log(x)))
