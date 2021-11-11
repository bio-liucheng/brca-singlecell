options(stringsAsFactors = F)
meta <- read.csv("meta.csv")
meta <- read.csv("meta_tcr.csv")
library(ggplot2)
library(ggpubr)

#subset
library(Seurat)
seurat.obj.md <- seurat.obj[,meta$X]
seurat.obj.md <- NormalizeData(seurat.obj.md)
seurat.obj.md <- ScaleData(seurat.obj.md, features = rownames(seurat.obj.md))
seurat.obj.md$cell_label <- meta$cell_label

seurat.obj.md <- subset(seurat.obj.md, cell_label %in% c("Macro-C4", "Macro-C3"))
#get signature
setwd("G:/scRNA-seq/LC/TCGA/BRCA")
signature <- read.csv("signature.csv")
colnames(signature)
sig_select = "exhausted"

sig_select <- signature[,sig_select]

sig_select <- sig_select[sig_select!= ""]

#sig_select <- c("UGDH", "RGN", "ALDH2")

#sig_select <- ascor$SYMBOL

#normal_type <- unique(normal$cluster)
#sig_select <- normal$gene[normal$cluster == normal_type[1]]
  
  
expr_scale <- seurat.obj.md[["RNA"]]@data
expr_scale <- expr_scale[intersect(sig_select, rownames(expr_scale)),]
flag <- apply(expr_scale, 1, function(x) sum(x>0.1))


sig_cell <- apply(expr_scale[flag>100,], 2, mean)
seurat.obj.md$signature <- sig_cell
seurat.obj.md$cell_label <- meta$cell_label
#FeaturePlot(seurat.obj.md, features = "signature", reduction = "ddr")
VlnPlot(seurat.obj.md, features = "signature", group.by = "cell_label", split.by = "sample_type")

sig_frame <- data.frame(sig_score = sig_cell, cell_names = colnames(seurat.obj.md),
                        cell_label = meta$cell_label, 
                        sample_type = seurat.obj.md$sample_type)
sig_frame <- data.frame(M1= m1, M2 = m2, cell_names = colnames(seurat.obj.md),
                        cell_label = meta$cell_label, 
                        sample_type = meta$sample_type)
ggplot(sig_frame, aes(M1, M2, color = cell_label)) + geom_point() + stat_cor(label.x.npc = 0.7, label.y.npc = 0.3) +  theme_classic() + 
  scale_color_manual(values = color_map_s$cell_color) + labs(x = "M1 signature", y = "M2 signature") +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank())

ggplot(sig_frame, aes(sample_type, sig_score, color = sample_type, group = sample_type)) + 
  geom_violin(width = 0.7, scale = "width", kernel = "gaussian") + theme_classic() +
  theme(legend.position = NULL) 
#  stat_compare_means(method = "wilcox.test",label = "p.format")
#beautiful plot
library(ggpubr)

library(ggbeeswarm)
library(ggpubr)

sig_frame$sample_type <- factor(sig_frame$sample_type, levels= c("Tumor", "Lymph Node"))
select_cluster <- c("CD4-C5-FOXP3", "CD4-C6-CXCL13", "CD8-C3-GZMK", "CD8-C5-CXCL13")
flag <- sig_frame$cell_label %in% select_cluster
flag_outline <- 1:nrow(sig_frame) %in% e$iRight
ggplot(sig_frame, aes(cell_label, sig_score, color = cell_label, group = cell_label)) + 
  geom_boxplot(width = 0.6, outlier.size = -1) + geom_quasirandom(width = 0.2, cex= 0.3, alpha = 0.5) + theme_classic() +
  theme(legend.position = "none", 
        legend.title= element_blank(), 
        legend.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 30, hjust = 1), 
        axis.text.y = element_text(size = 12),
        strip.background = element_blank()
        ) + scale_color_manual(values = color_map_s$cell_color) + 
#  scale_color_manual(values = color_map_s$cell_color) +
  labs(y = "Anti inflammatory score")
ggsave("CD8_Anti_inflammatory_score.pdf", width = 5, height = 3, dpi = 300)


ggplot(sig_frame, aes(sample_type, sig_score, color = sample_type, group = sample_type)) + 
  geom_boxplot(width = 0.6, outlier.size = -1) + geom_quasirandom(width = 0.2, cex= 0.3, alpha = 0.5) + theme_classic() +
  theme(legend.position = "none", 
        legend.title= element_blank(), 
        legend.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 30, hjust = 1), 
        axis.text.y = element_text(size = 12),
        strip.background = element_blank()
  )  +
  scale_color_manual(values = color_map_s2$cell_color) +
  labs(y = "CD8 T activation signature")
ggsave("sample_type_CD8 T activation signature_score_no_pvalue.png", width = 3, height = 3, dpi = 300)
