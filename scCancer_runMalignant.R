
seurat.obj.md <- subset(seurat.obj, cell_label == "Epithelial Cells")
seurat.obj.md <- seurat.obj[,meta$X]
seurat.obj.md <- NormalizeData(seurat.obj.md)
seurat.obj.md <- FindVariableFeatures(seurat.obj.md)
seurat.obj.md <- ScaleData(seurat.obj.md)
seurat.obj.md <- RunPCA(seurat.obj.md)
seurat.obj.md <- FindNeighbors(seurat.obj.md, reduction = "pca", dims = 1:20)


cell.annotation$cellName <- cell.annotation$barcodes
rownames(cell.annotation) <- cell.annotation$cellName


cell.annotation <- data.frame(barcodes = colnames(seurat.obj.md), Cell.Type = "Observation")


umap <- read.csv("scanpy_umap.csv")

colnames(umap) <- c('num', 'UMAP_1', 'UMAP_2')

umap <- umap[,c(2,3)]

cell.annotation <- data.frame(cell.annotation, tSNE_1 = umap$UMAP_1, tSNE_2 = umap$UMAP_2)

cell.annotation$Cluster <- paste("cluster_", meta$leiden)

colnames(list) <- c("EnsemblID", "ENTREZID", "SYMBOL")
list <- list[!duplicated(list$SYMBOL),]
rownames(list) <- list$SYMBOL

s <- runMalignancy(seurat.obj.md, list, cell.annotation = cell.annotation, savePath = "direct")
ggsave(filename = file.path(savePath, "figures/cellCycle-point.png"))
