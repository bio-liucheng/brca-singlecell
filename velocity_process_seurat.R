seurat.obj <- subset(seurat.obj, cells = myelod_cell)
seurat.obj <- subset(seurat.obj, cell_label_modify3 != "T")
seurat.obj <- subset(seurat.obj, leiden.res.1.8 == 23)

library(harmony)

process_seurat <- function(seurat.obj){
seurat.obj <- NormalizeData(seurat.obj)
seurat.obj <- FindVariableFeatures(seurat.obj, dispersion.cutoff = c(0.5, Inf))
seurat.obj <- ScaleData(seurat.obj)

seurat.obj <- RunPCA(seurat.obj)
seurat.obj <- RunHarmony(seurat.obj, "patient_id",  theta = 2, dims.use =1:50)

seurat.obj <- FindNeighbors(seurat.obj, reduction = "harmony", dims = 1:20)

#seurat.obj <- FindClusters(seurat.obj, resolution = 0.8, algorithm = 4)
#seurat.obj <- FindClusters(seurat.obj, resolution = 1.5)
seurat.obj <- RunUMAP(seurat.obj,reduction = "harmony", 
                          dims = 1:20)
return(seurat.obj)
}

process_velocity <- function(seurat.obj.md){
seurat.obj.md <- SCTransform(seurat.obj.md, assay = "sf", new.assay.name = "spliced")
seurat.obj.md <- SCTransform(seurat.obj.md, assay = "uf", new.assay.name = "unspliced")
DefaultAssay(seurat.obj.md) <- "spliced"
seurat.obj.md <- RunPCA(seurat.obj.md)
seurat.obj.md <-FindNeighbors(seurat.obj.md)
seurat.obj.md <- FindClusters(seurat.obj.md)
seurat.obj.md <- RunUMAP(seurat.obj.md, dims = 1:20)
seurat.obj.md <- RunVelocity(seurat.obj.md, ncores = 1, reduction = "pca", verbose = FALSE)
}





spliced_minmax <- 0.2
spliced <- filter.genes.by.cluster.expression(
  velocity.obj[["sf"]]@counts, seurat.obj.md$cell_label,
  min.max.cluster.average = spliced_minmax
)

unspliced_minmax <- 0.05
unspliced <- filter.genes.by.cluster.expression(
  velocity.obj[["uf"]]@counts, seurat.obj.md$cell_label,
  min.max.cluster.average = unspliced_minmax
)

velocity <- gene.relative.velocity.estimates(
  spliced, unspliced,
  deltaT       = 1,
  kCells       = 30,
  fit.quantile = 0.02,
  n.cores      = 1
)

umap_embedding <- show.velocity.on.embedding.cor(
  Embeddings(seurat.obj.md, "umap"), velocity,
  n.cores = 1, show.grid.flow = TRUE, return.details = TRUE
)


process_velocity <- function(seurat.obj.md){
  seurat.obj.md <- SCTransform(seurat.obj.md, assay = "sf", new.assay.name = "spliced")
  seurat.obj.md <- SCTransform(seurat.obj.md, assay = "uf", new.assay.name = "unspliced")
  DefaultAssay(seurat.obj.md) <- "spliced"
  seurat.obj.md <- RunPCA(seurat.obj.md)
  seurat.obj.md <-FindNeighbors(seurat.obj.md)
  seurat.obj.md <- FindClusters(seurat.obj.md)
  seurat.obj.md <- RunUMAP(seurat.obj.md, dims = 1:20)
  seurat.obj.md <- RunVelocity(seurat.obj.md, ncores = 1, reduction = "pca", verbose = FALSE)
}

umap_arrows <- umap_embedding$garrows %>%
  as.data.frame() %>%
  mutate(x2 = x0 + (x1 - x0) * 2,
         y2 = y0 + (y1 - y0) * 2)



show_seurat_marker <- function(seurat.obj){
seurat.obj.maker <- FindAllMarkers(seurat.obj, only.pos = TRUE, 
                                        min.pct = 0.25, 
                                        ogfc.threshold = 0.25, 
                                        max.cells.per.ident = 1000)

library(dplyr)
n <- 10
top_20 <- seurat.obj.maker %>% arrange(cluster, desc(avg_logFC))  %>% group_by(cluster) %>% top_n(n, avg_logFC)
DoHeatmap(seurat.obj, features = top_20$gene)

}
