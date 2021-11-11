library(harmony)

seurat.obj <- CreateSeuratObject(seurat.obj, project = "adjust", min.cells = 0, min.features = 0)

seurat.obj <- NormalizeData(seurat.obj)
seurat.obj <- FindVariableFeatures(seurat.obj, nfeatures = 3000)

seurat.obj[["percent.ribo"]] <- PercentageFeatureSet(seurat.obj, 
                                                     pattern = "^RPL|^RPS")




#seurat.obj <- ScaleData(seurat.obj, vars.to.regress = c("percent.mt"))
seurat.obj <- ScaleData(seurat.obj)
seurat.obj <- RunPCA(seurat.obj)
seurat.obj <- RunHarmony(seurat.obj, "orig.ident", dims.use =1:50,theta = 2, max.iter.harmony = 10)

seurat.obj <- FindNeighbors(seurat.obj, reduction = "pca", dims = 1:30)
seurat.obj <- FindNeighbors(seurat.obj, reduction = "harmony", dims = 1:30)
library(leiden)
library(igraph)
graph_object <- graph_from_adjacency_matrix(seurat.obj@graphs$RNA_nn, mode = "directed")
adjacency_matrix <- igraph::as_adjacency_matrix(graph_object)
partition1 <- leiden(adjacency_matrix, resolution_parameter = 1)
partition1.2 <- leiden(adjacency_matrix, resolution_parameter = 1.2)
partition1.5 <- leiden(adjacency_matrix, resolution_parameter = 1.5)
partition2 <- leiden(adjacency_matrix, resolution_parameter = 2)
seurat.obj$leiden.res.1 <- partition1
seurat.obj$leiden.res.1.2 <- partition1.2
seurat.obj$leiden.res.1.5 <- partition1.5
seurat.obj$leiden.res.2 <- partition2

seurat.obj <- RunUMAP(seurat.obj,reduction = "harmony", dims = 1:20, min.dist = 0.3)

seurat.obj <- RunUMAP(seurat.obj,graph = "RNA_snn")
seurat.obj <- RunTSNE(seurat.obj,reduction = "pca", dims = 1:20,tsne.method = "FIt-SNE" )
seurat.obj <- RunUMAP(seurat.obj,reduction = "pca", dims = 1:30, min.dist = 0.3)


seurat.obj <- FindClusters(seurat.obj, resolution = 1)

Idents(seurat.obj) <- seurat.obj$leiden.res.1.2
seurat.obj.maker_1 <- FindAllMarkers(seurat.obj, only.pos = TRUE, min.pct = 0.25, 
                                     logfc.threshold = 0.25, max.cells.per.ident = 2000, min.diff.pct = 0.25)

Idents(seurat.obj) <- seurat.obj$leiden.res.2
seurat.obj.maker_2 <- FindAllMarkers(seurat.obj.md, only.pos = TRUE, min.pct = 0.25, 
                                     logfc.threshold = 0.25, max.cells.per.ident = 3000, min.diff.pct = 0.25)

seurat.obj.maker_2 <- FindMarkers(seurat.obj,ident.1 = "4.5-IL6ST", only.pos = TRUE, min.pct = 0.25, 
                                     logfc.threshold = 0.25,min.diff.pct = 0.25)


seurat.obj$main_label[seurat.obj$leiden.res.1.2 == 13] <- "Fibro-DCN"
seurat.obj$main_label[seurat.obj$leiden.res.1.2 == 12] <- "Fibro-SOD3"

DimPlot(seurat.obj, reduction = "umap", group.by = "big_label")
FeaturePlot(seurat.obj, features =c("DCN"))
VlnPlot(seurat.obj.tm, features = c("ACTA2", "CD4"))

seurat.obj$main_label[colnames(seurat.obj.im)[seurat.obj.im$leiden.res.1 == 16]] <- "Mast cells"
saveRDS(seurat.obj, "seurat.obj.filter.RDS")
seurat.obj$big_label[colnames(seurat.obj)[seurat.obj$big_label == "KI67 high cells"]] <- "KI67+ cells"

seurat.obj$cell_label <- plyr::mapvalues(as.vector(seurat.obj.im$leiden.res.1.5), 
                                            as.vector(top_20$cluster), 
                                            as.vector(top_20$label))

seurat.obj.im$cell_label_tmp <- plyr::mapvalues(as.vector(colnames(seurat.obj.im)), 
                                         as.vector(colnames(seurat.obj)), 
                                         as.vector(seurat.obj$leiden.res.1.8))
