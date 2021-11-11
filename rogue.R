library(ggplot2)
library(tidyverse)
library(ROGUE)

expr <- as.matrix(seurat.obj.md[["RNA"]]@counts)
ent.res <- SE_fun(expr)
SEplot(ent.res)
varible.gene <- ent.res$Gene[ent.res$p.value < 0.1]
rogue.res <- rogue(expr, labels = seurat.obj.im$leiden.res.1, samples = seurat.obj.im$patient_id, platform = "UMI", span = 0.6)


seurat.obj.md <- NormalizeData(seurat.obj.md)
seurat.obj.md <- FindVariableFeatures(seurat.obj.md)
seurat.obj.md <- ScaleData(seurat.obj.md, features = varible.gene)

seurat.obj.md <- RunPCA(seurat.obj.md)
seurat.obj.md <- RunHarmony(seurat.obj.md, "orig.ident", dims.use =1:20)

seurat.obj.md <- FindNeighbors(seurat.obj.md, reduction = "harmony", dims = 1:20)
seurat.obj.md <- FindClusters(seurat.obj.md, resolution = 1)
seurat.obj.md <- RunUMAP(seurat.obj.md,reduction = "harmony", dims = 1:20)

seurat.obj.md.maker_1 <- FindAllMarkers(seurat.obj.md, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
