#Dotplot
library(ggplot2)
DotPlot(seurat.obj, features = unique(top_20$gene), cols = tu_vs_ly_color, split.by = "sample_type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
DotPlot(seurat.obj.md, features = unique(top_20$gene), cols = c('white',"dark red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#counts
analysis_label = "cell_label"

metadata <- seurat.obj[[]]
counts <- table(meta$analysis_label, meta$sample_type)
counts <- table(meta[,analysis_label], meta[,sample_type])
counts <- as.data.frame(counts)
#¶ÑµþÖùÐÎÍ¼
ggplot(meta, aes(factor(sample_type, levels = c("Tumor", "Lymph Node")))) + 
  geom_bar(aes(fill = cell_label), width = 0.5, position = "fill")  + theme_classic() +
  scale_fill_manual(values = cl[1:11]) +
  theme(axis.ticks = element_blank(),
        axis.title.x =element_blank(),
        legend.title =element_blank(),
        axis.text.x = element_text(face="bold",size=10),
        axis.text.y = element_text(face="bold",size=10),
        axis.title.y =element_blank(),
        strip.text.x = element_text(face="bold",size=10),
        
        legend.text = element_text(size=10)
  )



#dot percent
library(reshape2)

main_type <- unique(meta_bk$cell_label)
meta <- meta_bk[meta_bk$cell_label %in% "Myeloid cells", ]

counts <- table(meta$main_label, meta$sample_type)
counts <- as.data.frame(counts)
counts <- dcast(counts, Var1~Var2)
rownames(counts) <- counts$Var1
counts <- counts[,-1]
counts_normal <- apply(counts, 2, function(x) {x/sum(x)*100})
counts_normal <- as.data.frame(counts_normal)
colnames(counts_percent) <- paste0(colnames(counts_percent), "_M")

counts_percent <- apply(counts_normal, 1, function(x) {x/sum(x)*100})
counts_percent <- t(counts_percent)
counts_percent <- as.data.frame(counts_percent)
colnames(counts_percent) <- paste0(colnames(counts_percent), "_X")

counts_all <- cbind(counts, counts_percent)
counts_all$cell_label <- rownames(counts_all)

level = counts_all$cell_label[order(counts_all$`Tumor_X`)]

ggplot(counts_all,aes(factor(cell_label, levels = level), Tumor_X)) +geom_point(colour = "grey", aes(size = Tumor)) +
  geom_point(data = counts_all, aes(factor(cell_label, levels = level), `Lymph Node_X`, size = `Lymph Node`),colour = "grey")+
  theme_classic() +
  scale_size_continuous(breaks = c(100, 1000,10000), trans = "sqrt")+
  theme(axis.ticks = element_blank(),
        axis.title.x =element_blank(),
        legend.title =element_blank(),
        axis.text.x = element_text(face="bold",size=10,angle = 40, hjust = 1),
        axis.text.y = element_text(face="bold",size=10),
        axis.title.y =element_blank(),
        strip.text.x = element_text(face="bold",size=10),
        legend.text = element_text(size=10)
        #    axis.text.x = element_text(angle = 30, hjust = 1))
        
  )


#For continuous scales, the name of a transformation
#object or the object itself. Built-in transformations 
#include "asn", "atanh", "boxcox", "date", "exp", "hms", "identity", 
#"log", "log10", "log1p", "log2", "logit", "modulus", "probability", 
#"probit", "pseudo_log", "reciprocal", "reverse", "sqrt" and "time".


#important dot plot
ggplot(counts_all,aes(factor(cell_label, levels  = color_map_s$cell_label), Tumor_X)) +geom_point(colour = "#F0A19D", aes(size = Tumor))+ 
  geom_point(data = counts_all, aes(factor(cell_label, levels = color_map_s$cell_label), `Lymph Node_X`, size = `Lymph Node`),colour = "#80C9AE")+
  theme_classic() +
  scale_size_continuous(breaks = c(100, 1000,10000), trans = "sqrt")+
  theme(axis.ticks = element_blank(),
        axis.title.x =element_blank(),
        legend.title =element_blank(),
        axis.text.x = element_text(face="bold",size=10,angle = 90, hjust = 1),
#        axis.text.x = element_blank(),
        axis.text.y = element_text(face="bold",size=10),
        axis.title.y =element_blank(),
        strip.text.x = element_text(face="bold",size=10),
        legend.position = "None",
        legend.text = element_text(size=10),
        #    axis.text.x = element_text(angle = 30, hjust = 1))
        axis.line = element_line(colour='black')
        
  ) +geom_point(data = color_map_s, aes(x = factor(cell_label, levels  = color_map_s$cell_label), y = -10, color = factor(cell_label, levels  = color_map_s$cell_label)), size = 5) +
  scale_color_manual(values = color_map_s$cell_color) +
  set_font(p, family="Arial")

ggsave("dot_plot_all_cell_label.pdf", width = 10, height = 5)


DotPlot(seurat.obj, features = unique(top_20$gene))

level = 
p <- ggplot(meta, aes(factor(sample_type, levels = c("Tumor", "Lymph Node")))) + 
  geom_bar(aes(fill = cell_label), width = 0.6, position = "fill") + 
  scale_fill_manual(values = cl[1:11]) + theme_classic() +
  scale_y_continuous(labels = percent) + facet_wrap(~patient_id, scales = "free_x", nrow = 2) +
  theme(axis.ticks = element_blank(),
        axis.title.x =element_blank(),
        legend.title =element_blank(),
        axis.text.x = element_text(face="bold",size=10,angle = 30, hjust = 1),
        axis.text.y = element_text(face="bold",size=10),
        axis.title.y =element_blank(),
        strip.text.x = element_text(face="bold",size=10),
        legend.text = element_text(size=10)
#    axis.text.x = element_text(angle = 30, hjust = 1))

)

tu_vs_ly_color = c("#FF7F0E","#1F77B4")
ggplot(meta, aes(cell_label)) + 
  geom_bar(aes(fill =factor(sample_type, levels = c("Tumor", "Lymph Node") )), 
               position = "fill") + 
  scale_fill_manual(values = tu_vs_ly_color) + theme_classic() +
  scale_y_continuous(labels = percent) +
  theme(axis.ticks = element_blank(),
        axis.title.x =element_blank(),
        legend.title =element_blank(),
        axis.text.x = element_text(face="bold",size=10,angle = 40, hjust = 1),
        axis.text.y = element_text(face="bold",size=10),
        axis.title.y =element_blank(),
        strip.text.x = element_text(face="bold",size=10),
        legend.text = element_text(size=10)
        #    axis.text.x = element_text(angle = 40, hjust = 1))
        
  )

tu_vs_ly_color = c("#1F77B4","#FF7F0E")
library(dplyr)
meta %>% group_by(cell_label) %>% summarize(frac = n())


df <- meta %>% group_by(sample_type, cell_label) %>% 
  summarize(frac = n()) 
df2 <- meta %>% group_by(sample_type) %>% 
  summarize(frac = n()) 
df$frac2 <- plyr::mapvalues(df$sample_type, df2$sample_type, df2$frac)
df$frac2 <- as.integer(df$frac2)
df <- mutate(df, radio = frac/frac2)

ggplot(df, aes(factor(cell_label, levels  = rev(color_map_s$cell_label)), radio, fill = factor(sample_type, levels = c("Tumor", "Lymph Node") ))) +
  geom_bar(stat = 'identity', position = "fill")  + theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.title.x =element_blank(),
        legend.title =element_blank(),
        axis.text.x = element_text(face="bold",size=10, hjust = 1),
        axis.text.y = element_text(face="bold",size=10),
        axis.title.y =element_blank(),
        strip.text.x = element_text(face="bold",size=10),
        legend.text = element_text(face="bold",size=10)
        #    axis.text.x = element_text(angle = 30, hjust = 1))
        
  ) + coord_flip()+ scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(values = color_map_s2$cell_color)

ggsave("sample_type_cell_label_scale.pdf", width = 6.5, height = 4)
#fig 7b

percent <- table(meta$cell_label, meta$sample_type)
percent <- as.data.frame(percent)

ggplot(meta, aes(factor(cell_label, levels  = rev(color_map_s$cell_label)))) + 
  geom_bar(aes(fill = factor(Her2_statu, levels = c("positive", "Negative"))), width = 0.85, position = "fill")  + theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.title.x =element_blank(),
        legend.title =element_blank(),
        axis.title.y =element_blank(),
        axis.text= element_text(size = 12),
        legend.text = element_text(size = 12)
  ) + coord_flip() + 
  scale_y_continuous(labels=scales::percent)

ggsave("HER2_cell_label.pdf", width = 6.5, height = 4)


#fig S2b
meta_s <- meta %>% filter(sample_type == "Lymph Node")

percent <- table(meta_s$main_label, meta_s$patient_id)
percent <- as.data.frame(percent)

ggplot(meta_s, aes(factor(patient_id, levels = rev(unique(meta_s$patient_id))))) + 
  geom_bar(aes(fill = factor(main_label, levels = color_map_s$cell_label)), width = 0.85, position = "fill")  + theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.title.x =element_blank(),
        legend.title =element_blank(),
        axis.title.y =element_blank(),
        axis.text= element_text(size = 12),
        legend.text = element_text(size = 12)
  ) + coord_flip() + scale_fill_manual(values = color_map_s$cell_color) + scale_y_continuous(labels=scales::percent)

ggsave("ln_main_label.pdf", width = 6, height = 3.5)


#fig 7b

percent <- table(meta$leiden, meta$sample_type)
percent <- as.data.frame(percent)

ggplot(meta, aes(factor(seurat_clusters, levels = c(2,1,0,3)))) + 
  geom_bar(aes(fill = factor(sample_type, levels = color_map_s$cell_label)), width = 0.85, position = "fill")  + theme_classic() +
  theme(axis.ticks = element_blank(),
        legend.title =element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        axis.text= element_text(size = 12),
        legend.text = element_text(size = 12)
  ) + coord_flip() + scale_fill_manual(values = color_map_s$cell_color)+ scale_y_continuous(labels=scales::percent)

ggsave("percent_cxcl13_leiden.pdf", width = 6, height = 4)


percent <- table(meta$main_label, meta$sample_type)
percent <- as.data.frame(percent)
ggplot(meta, aes(factor(sample_type, levels = rev(c("Tumor", "Lymph Node"))))) + 
  geom_bar(aes(fill = leiden), width = 0.85, position = "fill")  + theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.title.x =element_blank(),
        legend.title =element_blank(),
        axis.title.y =element_blank(),
        legend.position = "bottom",
        axis.text= element_text(size = 13),
        legend.text = element_text(size = 13)
  ) + scale_y_continuous(labels=scales::percent)+ coord_flip() 

ggsave("cluster_transfer_sample_type.pdf", width = 6, height =2)

