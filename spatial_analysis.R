# spatial analysis
# load packages

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

space.data <- Load10X_Spatial(data.dir = "./A_1938345_11/")

space.data <- Load10X_Spatial(data.dir = "./B_1938529_9/")

space.data <- Load10X_Spatial(data.dir = "./C_2000752_23/")

space.data <- Load10X_Spatial(data.dir = "./D_2000910_33/")

brain <- SCTransform(space.data, assay = "Spatial", verbose = FALSE)

DefaultAssay(brain) = "Spatial"
SpatialFeaturePlot(brain, features = c("CTLA4", "EPCAM","FOXP3","CD8A"))


SpatialFeaturePlot(space.data_d, features = c("CD79A", "EPCAM","CD3E", "CCR7"))

brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE, resolution = 0.5)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)
p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1+p2


brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000],
                                       selection.method = "markvariogram")

top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "markvariogram"), 6)
SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1))

# reference data
DimPlot(Reference, group.by = "cell_label")

library(dplyr)
Reference <- SCTransform(Reference, ncells = 11370, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)

anchors <- FindTransferAnchors(reference = Reference, query = brain, normalization.method = "SCT")

predictions.assay <- TransferData(anchorset = anchors, refdata = Reference$cell_label, prediction.assay = TRUE,
                                  weight.reduction = brain[["pca"]], dims = 1:30)
brain[["predictions"]] <- predictions.assay
DefaultAssay(brain) <- "predictions"
SpatialFeaturePlot(brain, features = c("CD4-C5-FOXP3"), pt.size.factor = 1.6, crop = TRUE)
ggsave("test.pdf", width = 4.5, height = 4.5)


for(i in unique(Reference$cell_label)){
  SpatialFeaturePlot(brain, features =i, pt.size.factor = 1.6, crop = TRUE)
  ggsave(paste0("patient_8_", i, ".png"), width = 4, height = 4)
}


cortex <- FindSpatiallyVariableFeatures(cortex, assay = "predictions", selection.method = "markvariogram",
                                        features = rownames(cortex), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(cortex), 4)
SpatialPlot(object = cortex, features = top.clusters, ncol = 2)
SpatialFeaturePlot(cortex, features = c("Astro", "L2/3 IT", "L4", "L5 PT", "L5 IT", "L6 CT", "L6 IT",
                                        "L6b", "Oligo"), pt.size.factor = 1, ncol = 2, crop = FALSE, alpha = c(0.1, 1))