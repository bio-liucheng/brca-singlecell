library(reshape2)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
seurat.obj.md.md <- seurat.obj[,meta[meta$clonotype == "6292",]$X]
table(seurat.obj.md.md$sample_type, seurat.obj.md.md$cell_label)
expr_scale <- seurat.obj.md.md[["RNA"]]@data
expr_scale <- expr_scale[intersect(sig_select, rownames(expr_scale)),]
sig_cell <- apply(expr_scale, 2, mean)
seurat.obj.md.md$signature <- sig_cell
VlnPlot(seurat.obj.md.md, features = "IFNG", group.by = "sample_type")


seurat.obj.md.md <- seurat.obj[,meta[meta$clonotype == "6273",]$X]

meta_cd4_cxcl13 <- meta[meta$cell_label == "CD4-C5-CXCL13",]
meta_cd4_cxcl13 <- meta_cd4_cxcl13[order(meta_cd4_cxcl13$ln_colonal + meta_cd4_cxcl13$tumor_colonal, decreasing = T),]


s <- table(meta_ifng$clonotype, meta_ifng$sample_type)
s <- as.data.frame(s)
s <- dcast(s, Var1 ~ Var2)

meta_ifng$tumor_colonal <- plyr::mapvalues(meta_ifng$clonotype, s$Var1, s$Tumor)
meta_ifng$ln_colonal <- plyr::mapvalues(meta_ifng$clonotype, s$Var1, s$`Lymph Node`)


gene <- "GPX4"
ifng <- seurat.obj.md[["RNA"]]@data[gene,]
meta_ifng <- data.frame(meta_migration, ifng = ifng)
meta_ifng <- meta_ifng[grepl("CD8-C4-HSPA1A", meta_ifng$cell_label),]
meta_ifng <- meta_ifng[meta_ifng$tumor_colonal >1 & meta_ifng$ln_colonal >1,]
meta_ifng_s <- select(meta_ifng, clonotype, sample_type, cell_label, tumor_colonal, ln_colonal, ifng)
ifng_dcast <- dcast(meta_ifng_s, clonotype~sample_type, mean)
ifng_result <- melt(ifng_dcast)
ifng_result <- melt(ifng_dcast, id.vars = "clonotype")
ggplot(ifng_result, aes(variable, value)) + geom_boxplot() + geom_line(aes(group = clonotype))  + stat_compare_means(comparisons = list(c("Lymph Node", "Tumor")))

ggplot(ifng_result, aes(variable, value)) + geom_boxplot() + geom_quasirandom() +
 stat_compare_means(comparisons = list(c("Lymph Node", "Tumor")))



table(ifng_dcast$Tumor- ifng_dcast$`Lymph Node` >0)

meta_migration <- meta[meta$migration,]
meta_migration <- meta_migration[grepl("CD4", meta_migration$cell_label),]
gene <- "CTLA4"
gene_expr <- seurat.obj[,meta_migration$X][["RNA"]]@data[gene,]
ifng_result <- select(meta_migration, clonotype, sample_type, cell_label)
ifng_result <- data.frame(meta_migration, gene_expr = gene_expr)

ggplot(ifng_result[grepl("CD8", ifng_result$cell_label),], aes(factor(sample_type, levels = color_map_s$cell_label), gene_expr, color = factor(sample_type, levels = color_map_s$cell_label)))  + 
  geom_boxplot(outlier.size = 0) +geom_quasirandom(cex = 0.3)+ 
  scale_colour_manual(values = color_map_s$cell_color) +
  facet_wrap(~cell_label, nrow = 1) +
  stat_compare_means(comparisons = list(c("Lymph Node", "Tumor")),label.y = 3.8) + theme_classic() +
  theme(legend.title = element_blank(), strip.background = element_blank(), 
        axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 10)
        ) + labs(y = gene) + ylim(c(0,4))
ggsave("CD8_cell_TGFB1.jpeg", width = 7, height = 3)

ifng_result[grepl("CD8", ifng_result$cell_label),]

ggplot(subset(ifng_result, cell_label == "CD4-C5-FOXP3"), aes(factor(sample_type, levels = color_map_s$cell_label), gene_expr, color = factor(sample_type, levels = color_map_s$cell_label)))  + 
  geom_boxplot(outlier.size = 0, width = 0.4) +geom_quasirandom(cex = 0.3, width = 0.4)+ 
  scale_colour_manual(values = color_map_s$cell_color) +
  stat_compare_means(comparisons = list(c("Lymph Node", "Tumor")),label.y =3.8) + theme_classic() +
  theme(legend.title = element_blank(), strip.background = element_blank(), 
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
        axis.title.x = element_blank(),
        legend.position = "none",
  ) + labs(y = gene, title = "CD4-C5-FOXP3")
ggsave( paste0("CD4-C5-FOXP3_", gene, ".png"), width = 2.5, height = 3, dpi = 300)
