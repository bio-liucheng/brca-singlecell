seurat.obj <- readRDS("D:/single_cell/analysis_all/seurat.obj.filter.patient_id.RDS")

metadata <- read.csv("resolution2.csv")
Idents(seurat.obj) <- metadata$leiden

umap <- read.csv("all_scanpy_umap.csv")
umap <- umap[,-1]
umap <- as.matrix(umap)
rownames(umap) <- colnames(seurat.obj)
colnames(umap) <- paste0("UMAP_", 1:2)
seurat.obj[["umap"]] <- CreateDimReducObject(embeddings = umap, 
                                      key = "UMAP_", assay = DefaultAssay(seurat.obj))

DimPlot(seurat.obj, reduction = "umap")


meta_all$cell_label <- " "
cell_label <- unique(meta$cell_label)
for(i in 1:length(cell_label)){
  meta_all[rownames(meta)[meta$cell_label == cell_label[i]],]$cell_label <- cell_label[i]
}

write.csv(meta_all, "meta.csv", row.names = F)
