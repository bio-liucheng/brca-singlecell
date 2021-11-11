
t <- table(meta$leiden, meta$new_label)
t <- as.data.frame(t)
t <- melt(t)
t <- dcast(t, Var1 ~ Var2)
t <- t[,-1]
s  <- apply(t, 1, function(x) {which(x == max(x))})
new_name <- colnames(t)[s]
meta$new_label <- plyr::mapvalues(meta$leiden, 
                                  1:nrow(t), as.vector(new_name))


library(reshape2)

Idents(seurat.obj.filter) <- meta2$leiden

seurat.obj$cell_label <- " "
seurat.obj$cell_label <- plyr::mapvalues(colnames(seurat.obj), as.vector(meta2$X), as.vector(meta2$cell_label))
seurat.obj$cell_label[grepl("_",seurat.obj$cell_label)] <- "A"
cell_label <- unique(meta$cell_label)
for(i in 1:length(cell_label)){
  meta$cell_label[colnames(seurat.obj.md)[seurat.obj.md$cell_label == cell_label[i]]] <-cell_label[i]
  
}
meta$cell_label[meta$cell_label == 0] <- "CDT-mono"


meta <- read.csv("meta.csv")

s <- table(meta$leiden.res.2,meta$cell_label)
library(reshape2)
s <- as.data.frame(s)
s <- dcast(s, Var1~Var2)
rownames(s) <- s$Var1
s <- s[,-1]
t <- apply(s, 1, function(x){which(x == max(x))})
new_label <- colnames(s)[t]
meta$cell_label <- plyr::mapvalues(as.vector(meta$leiden.res.2), 
                                   1:(nrow(s)),
                                   as.vector(new_label))
meta$cell_label[meta$leiden %in% c(12, 14, 19)] <- "CD4-C4"
meta$cell_label[meta$cell_label == "CAFs"] <- "CAF-un"
meta$new_label[meta$new_label == "Plasmacytoid dendritic cell"] <- "Myeloid cells"
write.csv(meta, "meta.csv", row.names = F)

Idents(seurat.obj.md) <- meta$leiden
meta$cell_label <- plyr::mapvalues(meta$X, metao$X, as.character(metao$cell_label))
meta$cell_label <- as.character(meta$cell_label)
meta$cell_label[grepl("_", meta$cell_label)] <- "A"

s <- table(seurat.obj$cell_label, seurat.obj$orig.ident)
library(reshape2)
s <- as.data.frame(s)
s <- dcast(s, Var1~Var2)
rownames(s) <- s$Var1
s <- s[,-1]

s <- table(meta2$cell_label, meta2$orig.ident)
library(reshape2)
s <- as.data.frame(s)
s <- dcast(s, Var1~Var2)
rownames(s) <- s$Var1
s <- s[,-1]
