
#subsample
cell_pool <- unique(seurat.obj$cell_label)
cell_type <- grepl("CD8-", cell_pool)

cell_select <- WhichCells(seurat.obj.immue, expression = cell_label %in% cell_pool[cell_type] )

seurat.obj <- seurat.obj.immue[,cell_select]
seurat.obj <- process_seurat(seurat.obj)
#seurat to scanpy
library(Matrix)
writeMM(seurat.obj[["RNA"]]@counts, "cell_counts_matrix.mtx")
write.csv(seurat.obj[[]], "seurat_obj_immu.csv")
write.csv(Embeddings(seurat.obj, "harmony"), "immu_harmony.csv")
write.csv(Embeddings(seurat.obj, "umap"), "immu_umap.csv")
write.csv(colnames(seurat.obj), "immu_cellnames.csv")
write.csv(rownames(seurat.obj), "immu_genenames.csv")

saveRDS(seurat.obj, "seurat.other.programe.RDS")
