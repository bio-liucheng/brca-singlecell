library(Seurat)
library(harmony)
unique(seurat.obj$main_label)

seurat.obj.md <- subset(seurat.obj, new_label == "Malignant cells" )

seurat.obj.md <- subset(seurat.obj, main_label == "NK cells" |
                          main_label == "CD8 T" |
                           main_label == "CD4 T" |
                           main_label == "Immue cells"
                        )

seurat.obj.md <- subset(seurat.obj, main_label == "Myeloid cells")

seurat.obj.md <- subset(seurat.obj, new_label  == "CAFs" |
                          new_label  == "NK cells" |
                          new_label  == "CD8 T" | 
                          new_label  == "CD4 T"
                          
                        )
meta <- read.csv("meta_with_umap.csv")

seurat.obj.md <- seurat.obj[,meta$X]
seurat.obj.md <- subset(seurat.obj, cells = myelod_cell)
seurat.obj.md <- subset(seurat.obj.md, cell_label_modify3 != "T")
seurat.obj.md <- subset(seurat.obj, leiden.res.1.8 == 23)
seurat.obj.md[["percent.ribo"]] <- PercentageFeatureSet(seurat.obj.md, 
                                                        pattern = "^RPL|^RPS")

seurat.obj.md <- NormalizeData(seurat.obj.md)


seurat.obj.md <- FindVariableFeatures(seurat.obj.md)
#seurat.obj.md <- ScaleData(seurat.obj.md,  
#                           vars.to.regress = c("nUMI","percent.ribo", "percent.mt"))
#seurat.obj.md <- ScaleData(seurat.obj.md, vars.to.regress = "percent.ribo")

seurat.obj.md <- ScaleData(seurat.obj.md)


seurat.obj.md <- RunPCA(seurat.obj.md)
seurat.obj.md <- RunHarmony(seurat.obj.md, theta = 4,"group_patient", dims.use =1:50, max.iter.harmony = 20)

seurat.obj.md <- FindNeighbors(seurat.obj.md, reduction = "harmony", dims = 1:20)

seurat.obj.md <- FindClusters(seurat.obj.md, resolution = 0.5)
#seurat.obj.md <- FindClusters(seurat.obj.md, resolution = 1.5)
seurat.obj.md <- RunUMAP(seurat.obj.md, reduction = "pca", dims = 1:20, min.dist = 0.3)



library(leiden)
library(igraph)
graph_object <- graph_from_adjacency_matrix(seurat.obj.md@graphs$RNA_nn, mode = "directed")
adjacency_matrix <- igraph::as_adjacency_matrix(graph_object)
partition <- leiden(adjacency_matrix, resolution_parameter = 0.3)
seurat.obj.md$leiden.res.1 <- partition
Idents(seurat.obj.md) <- seurat.obj.md$leiden.res.1
partition <- leiden(adjacency_matrix, resolution_parameter = 1.5)

seurat.obj.md$leiden.res.1.5 <- partition
Idents(seurat.obj.md) <- seurat.obj.md$leiden.res.1.5

partition <- leiden(adjacency_matrix, resolution_parameter = 2)
seurat.obj.md$leiden.res.2 <- partition
Idents(seurat.obj.md) <- seurat.obj.md$leiden.res.2


seurat.obj.md.maker_1 <- FindAllMarkers(seurat.obj, only.pos = TRUE, 
                                       min.pct = 0.25, logfc.threshold = 0.25,max.cells.per.ident = 2000
                                      )

seurat.obj.md.maker_1 <- FindAllMarkers(seurat.obj.md, only.pos = TRUE, 
                                        min.pct = 0.25, logfc.threshold = 0.25,max.cells.per.ident = 2000
)

seurat.obj.md <- seurat.obj[,meta$X]
Idents(seurat.obj.md) <- meta$cell_label
Idents(seurat.obj.md) <- factor(meta$cell_label, levels = color_map_s$cell_label)
seurat.obj.md.maker_ng <- FindAllMarkers(monocle_prsudotime, only.pos = TRUE, min.diff.pct = 0.25,
                                      min.pct = 0.25, logfc.threshold = 0.25, 
                                      max.cells.per.ident = 2000
)

maker1 <- FindMarkers(seurat.obj.md, ident.1 ="CAF-C1-POSTN", ident.2 = "CAF-C3-PLA2G2A",
                                       min.pct = 0.25, logfc.threshold = 0.25,
                                       max.cells.per.ident = 20000)
                      
library(dplyr)
n <- 4

top_20 <- seurat.obj.md.maker_ng  %>% mutate(rank = 1:nrow(seurat.obj.md.maker_ng)) %>% group_by(cluster) %>% top_n(-n, rank)
top_20 <- seurat.obj.md.maker_ng   %>% group_by(cluster) %>% top_n(-n, p_val_adj)
top_20 <- seurat.obj.md.maker_ng  %>% group_by(cluster) %>% top_n(n, avg_logFC)
top_20 <- deg  %>% group_by(cluster) %>% top_n(n, avg_logFC)
num <- rep(1:n,length(unique(top_20$cluster)))
top_20$rank <- num
library(reshape2)
df <- melt(top_20, id.vars = c("cluster", "rank"), measure.vars = "gene")
df <- dcast(df, rank ~ cluster)
write.csv(top_20, "top_20_fc.csv")
write.csv(df, "top_20_im.csv")

#比较新cluster与之前的cluster
t <- table(seurat.obj.md$leiden.res.1, seurat.obj.md$cell_label_modify2)
t <- as.data.frame(t)
t <- dcast(t, Var1 ~ Var2)
t <- t[,-1]
s  <- apply(t, 1, function(x) {which(x == max(x))})
write.csv(data.frame(cluster = unique(seurat.obj.md$cell_label_modify2)), "tmp.csv")


#annotation 
top_20 <- read.csv("top_20_im.csv")
top_20 <- top_20[,c(-1,-2)]
top_20 <- t(top_20)
top_20 <- as.data.frame(top_20)
top_20$cluster <- unique(meta$leiden)
#top_20$label <- paste(top_20$V1, top_20$V2, sep = " ")

top_20$label <- top_20$V1
meta$cell_label <- plyr::mapvalues(as.vector(meta$leiden), 
                                            as.vector(top_20$cluster), 
                                            as.vector(top_20$label))


write.csv(meta, "meta.csv", row.names = F)

new_label <- read.csv("tmp.csv")
seurat.obj.md$cell_label_modify3 <- plyr::mapvalues(as.vector(seurat.obj.md$cell_label_modify2), 
                                            as.vector(new_label$cluster), 
                                            as.vector(new_label$Now_type))

seurat.obj.md$cell_label_modify3 <- "Mast cells"
seurat.obj$cell_label <- " "
cell_label <- unique(seurat.obj.md$cell_label)
for(i in 1:length(cell_label)){
  seurat.obj.md$cell_label_modify3[colnames(seurat.obj.md)[seurat.obj.md$cell_label == cell_label[i]]] <-cell_label[i]

}

for(i in 1:length(cell_label)){
  seurat.obj$cell_label[colnames(seurat.obj.md)[seurat.obj.md$cell_label == cell_label[i]]] <-cell_label[i]
  
}

seurat.obj.md <- RunUMAP(seurat.obj.md,reduction = "harmony", dims = 1:20, min.dist = 0.5)
seurat.obj.md <- RunTSNE(seurat.obj.md,reduction = "harmony", dims = 1:20)

DimPlot(seurat.obj.md, reduction = "umap", group.by = "sample_type")
DimPlot(seurat.obj.md, reduction = "umap", 
        group.by = "leiden.res.1", label = T, pt.size = 1.2)

FeaturePlot(seurat.obj.md, features =c("CD8A"))
VlnPlot(seurat.obj.md, features =c("CXCL11"))
DimPlot(seurat.obj.md, reduction = "umap", 
        group.by = c("sample_type","leiden.res.1"), label = T)

unique(seurat.obj$RNA_snn_res.1)
seurat.obj$figure_label[as.vector(seurat.obj$cell_label_modify2) == "Mast cells"] <- "Mast"
seurat.obj$main_label[colnames(seurat.obj.md)[seurat.obj.md$leiden.res.1.5 == 15]] <-"Imuue_cancer" 

saveRDS(seurat.obj.md, "seurat.obj.md.patient.RDS")
seurat.obj.md$myeloid_label[as.vector(seurat.obj.md$myeloid_label) == "Macro-CD1C"] <- "DC-CD1C"



colnames(seurat.obj.md.ma[[]])


flag <- grepl("-un|IFI44L|-KI67|-SCGB2A2", meta_all$cell_label2)

meta_all <- meta_all[!flag,]
