seurat.obj2[["tsne"]] <- CreateDimReducObject(embeddings = Embeddings(seurat.obj, "tsne"), 
                                                key = "TSNE_", 
                                                assay = DefaultAssay(seurat.obj2))

x <- HVFInfo(seurat.obj)

ggplot(x, aes(mean, variance.standardized)) + geom_point() 
