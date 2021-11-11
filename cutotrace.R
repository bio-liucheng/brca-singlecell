meta <- meta[grepl("_1|_6|_8|_7", meta$patient_id),]

library(ggplot2)
ggplot(meta, aes(sample_type, cytoTrance)) + geom_boxplot() + facet_wrap(~patient_id) + geom_jitter()


library(CytoTRACE)
meta <- read.csv("meta.csv")
seurat.obj.md <- seurat.obj[,meta$X]
seurat.obj.md <- NormalizeData(seurat.obj.md)
mat <- seurat.obj.md[["RNA"]]@data
mat <- as.matrix(mat)
s <- CytoTRACE(mat, enableFast = F, ncores = 1)

meta$CytoTrace <- s$CytoTRACE
ggplot(meta, aes(cell_label, CytoTrace)) + geom_boxplot() + geom_jitter()


