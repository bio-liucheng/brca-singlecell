library(ggplot2)
library(ggsci)
library(scales)
cl <- c(pal_nejm("default", alpha = 0.5)(8), pal_npg("nrc", 0.5)(10), pal_jama("default", 0.5)(7), pal_simpsons(palette = "springfield", 0.5)(16))
cl2 <- c(pal_nejm("default", alpha = 0.8)(8), pal_npg("nrc", 0.8)(10), pal_jama("default", 1)(7), pal_simpsons(palette = "springfield", 0.8)(16))


#Ö×ÁöÁÜ°ÍÑÕÉ«¶¨Òå
tumor_color <- pal_nejm("default", 0.8)(2)[1]
lynode_color <- pal_nejm("default", 0.8)(2)[2]
tu_vs_ly_color <- c(tumor_color, lynode_color)
#fig 1
DimPlot(seurat.obj, reduction = "umap", group.by = "big_label", label = T) + 
  theme_void() + scale_color_manual(values = cl[1:11])
#fig 2
DimPlot(seurat.obj, reduction = "umap", group.by = "sample_type") + 
  theme_void() + scale_color_manual(values = c(tumor_color, lynode_color)) + theme()

seurat.obj$sample_type[colnames(seurat.obj)[seurat.obj$sample_type == "CA"]] <- "Tumor"

immue_color <- cl[11:25]
immue_label <- as.factor(unique(seurat.obj$immue_label))
immue_color_to_label <- data.frame(immue_label = levels(immue_label), immue_color = immue_color)
#fig3
DimPlot(seurat.obj.md, reduction = "umap", group.by = "immue_label", cols = immue_color) + 
  theme_void() 
#fig4
DimPlot(seurat.obj.md, reduction = "umap", group.by = "sample_type", cols = immue_color) + 
  theme_void() 

#fig5

myeloid_color <- cl2[27:36]
myeloid_label <- as.factor(unique(seurat.obj.md2$myeloid_label))
myeloid_color_to_label <- data.frame(myeloid_label = levels(myeloid_label), myeloid_color = myeloid_color)

DimPlot(seurat.obj.md2, reduction = "umap", 
        group.by = "cell_label_modify3",
        cols = myeloid_color, label = T, repel = T
        ) +theme_void() 

#fig5
DimPlot(seurat.obj.md2, reduction = "umap", group.by = "sample_type",cols = tu_vs_ly_color ) +theme_void() 
library(Matrix)


