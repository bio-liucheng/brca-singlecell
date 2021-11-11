library(plyr)

cluster<- read.csv("cluster.csv")
seurat.obj$main_label <- mapvalues(as.vector(seurat.obj$RNA_snn_res.1.2), 
                                   as.vector(cluster$cluster), as.vector(cluster$main_label))
seurat.obj$maker_gene <- mapvalues(as.vector(seurat.obj$RNA_snn_res.1.2), 
                                   as.vector(cluster$cluster), as.vector(cluster$maker))
seurat.obj$putative_label <- mapvalues(as.vector(seurat.obj$RNA_snn_res.1.2), 
                                   as.vector(cluster$cluster), as.vector(cluster$putative_label))
options(stringsAsFactors = F)
label <- read.csv("tmp.csv")
seurat.obj$cell_label_modify2 <- plyr::mapvalues(as.vector(seurat.obj$cell_label_modify), 
                                                 as.vector(label$x), 
                                                 as.vector(label$x.1))
seurat.obj$median_label <- plyr::mapvalues(as.vector(seurat.obj$cell_label_modify), 
                                           as.vector(label$x), 
                                           as.vector(label$median_label))
seurat.obj$figure_label <- plyr::mapvalues(as.vector(seurat.obj$cell_label_modify), 
                                           as.vector(label$x), 
                                           as.vector(label$figure_label))
seurat.obj.md$immue_label <- plyr::mapvalues(as.vector(seurat.obj.md$leiden.res.1.2), 
                                             as.vector(label$X), 
                                             as.vector(label$name))
