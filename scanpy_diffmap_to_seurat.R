
#add diffusionmap

diff <- read.csv("scanpy_diffmap.csv")
meta <- read.csv("meta.csv")

mds <- diff[,c(3,4)]
mds <- as.matrix(mds)
colnames(mds) <- paste0("DIFF_", 1:2)
rownames(mds) <- meta$X
seurat.obj.md[["diff"]] <- CreateDimReducObject(embeddings = mds, 
                                                key = "DIFF_", 
                                                assay = DefaultAssay(seurat.obj.md))

DimPlot(seurat.obj.md, group.by = "sample_type", 
        reduction = "diff", split.by = "patient_id", pt.size = 1.5, ncol = 4)

#add diffusionmap

diff <- read.csv("scanpy_pca.csv")
meta <- read.csv("meta.csv")

mds <- diff[,c(2,3,4)]
mds <- as.matrix(mds)
colnames(mds) <- paste0("pca_", 1:3)
rownames(mds) <- meta$X
meta <- data.frame(meta, mds)

gg_scatter <- ggplot(meta, aes(pca_1, pca_2, color = factor(sample_type, levels = color_map_s$cell_label))) + 
  geom_point(size = 0.8) + theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0.5, 0.5), "cm")) +
  theme(legend.position="none") +
  scale_color_manual(values = color_map_s$cell_color)

gg_1 <- ggplot(meta, aes(pca_1, fill = factor(sample_type, levels = color_map_s$cell_label))) + 
  geom_density(alpha=.5)+ 
  theme_classic()+
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
          axis.text=element_blank(),
          axis.line=element_blank(),
          axis.ticks=element_blank(),
          axis.title.y = element_blank()) + 
  theme(plot.margin = unit(c(0.5, 0, 0, 0.7), "cm")) +
  scale_fill_manual(values = color_map_s$cell_color)

gg_2 <- ggplot(meta, aes(pca_2, fill = factor(sample_type, levels = color_map_s$cell_label))) + 
  geom_density(alpha=.5)+ 
  theme_classic() +
  theme(legend.position="none") + 
  theme(axis.title.x=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(), 
        axis.title.y = element_blank()) + coord_flip() +
  theme(plot.margin = unit(c(0, 0.5, 0.5, 0), "cm"))+
  scale_fill_manual(values = color_map_s$cell_color)


first_col = plot_grid(gg_1, gg_scatter, ncol = 1, rel_heights = c(1, 4))
second_col = plot_grid(NULL, gg_2, ncol = 1, rel_heights = c(1, 4))
perfect = plot_grid(first_col, second_col, ncol = 2, rel_widths = c(4, 1))
perfect

diff$X <- meta$sample_type

cell_label <- unique(meta$cell_label)

for(i in 1:length(cell_label)){
s <- hotelling.test(.~X, data = diff[meta$cell_label %in% cell_label[i], ])
print(s)
print(cell_label[i])
print("-------------")
}


for(i in 1:length(cell_label)){
  s <- hotelling.test(.~X, data = diff[meta$cell_label %in% cell_label[i], ])
  print(s)
  print(cell_label[i])
  print("-------------")
}


knn.sample_type <- function(index, index_name){
  index_name <- index_name[index]
  return(tumor = sum(index_name == "Tumor"))
}

for(i in 1:length(cell_label)){

matrix.cell_label <- diff[meta$cell_label %in% cell_label[i], ]
k <- sample(1:length(matrix.cell_label$X), 1000, replace = T)
k_n = 20
#knn.out <- queryKNN(matrix.cell_label[,-1], subset(matrix.cell_label, X != "Tumor")[,-1], k = k_n)
knn.out <- queryKNN(matrix.cell_label[,-1], matrix.cell_label[k,-1], k = k_n)
knn.index <- knn.out[["index"]]

t <- apply(knn.index, 1, knn.sample_type, index_name = matrix.cell_label$X)

print(cell_label[i])
print(paste0("knn_index: ", mean(t)))

print(sum(matrix.cell_label$X == "Tumor")/length(matrix.cell_label$X) * (k_n-1))

print("---------------")
}
