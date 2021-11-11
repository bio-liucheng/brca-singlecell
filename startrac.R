options(stringsAsFactors = F)
meta <- read.csv("meta_tcr.csv")
library(Startrac)
library(ggplot2)
library("ggpubr")
library("ComplexHeatmap")
library("circlize")

flag_t <- grepl("CD4|CD8", meta$cell_label)
meta <- meta[flag_t,]

flag_tcr <- meta$has_tcr == "True"

meta <- meta[flag_tcr,]


tcr <- data.frame(Cell_Name = as.character(meta$X),
                  clone.id = paste(meta$patient_id, meta$clonotype, sep = ":"),
                  patient = as.character(meta$patient_id),
                  majorCluster = as.character(meta$cell_label),
                  loc = as.character(meta$sample_type))

out <- Startrac.run(tcr, proj="brca",verbose=F)
plot(out,index.type="cluster.all",byPatient=F)
plot(out,index.type="cluster.all",byPatient=T)
out_tumor <- Startrac.run(subset(tcr, loc == "Tumor"), proj="brca",verbose=F)
out_ln <- Startrac.run(subset(tcr, loc == "Lymph Node"), proj="brca",verbose=F)

cluster_tumor <- out_tumor@cluster.data
cluster_ln <- out_ln@cluster.data

cluster_tumor <- cluster_tumor[grepl("patient", cluster_tumor$aid),]
cluster_tumor$sample_type = "Tumor"

cluster_ln <- cluster_ln[grepl("patient", cluster_ln$aid),]
cluster_ln$sample_type = "Lymph Node"

cluster_all <- rbind(cluster_tumor, cluster_ln)

ggplot(cluster_all, aes(sample_type, expa,colour = aid, group = aid)) +geom_point() + geom_line()+
  facet_wrap(~majorCluster, nrow = 3, scales = "free")


cluster_all$sample_type <- factor(cluster_all$sample_type, levels = color_map_s$cell_label)


ggboxplot(cluster_all, x = "majorCluster", 
          y = "expr", color = "majorCluster", 
          add = "jitter") + 
#  scale_color_manual(values = color_map_s$cell_color)  +
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank()) +
  labs(y = "STARTRAC-transition")
ggsave("CD8-transition.pdf", width = 5, height = 3)
# all cell label expansition and transition
ggboxplot(cluster_merge, x = "majorCluster", y = "migr", color = "majorCluster", add = "jitter")  + scale_color_manual(values = color_map_s$cell_color) +
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10), 
        axis.text = element_text(size = 12)) + 
  labs(y = "STARTRAC-migration")
ggsave("cell_label_STARTRAC-migration.pdf", width = 7, height = 4)

ggplot(cluster_merge, aes(majorCluster, expa,colour = majorCluster)) +
  geom_boxplot() + geom_jitter() + ylim(c(0, 0.25)) + scale_color_manual(values = color_map_s$cell_color) +
  theme_classic() +
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank()) + 
  labs(y = "STARTRAC-expansition")

# CD4 or CD8 boxplot Fig 4A
ggboxplot(cluster_all[grepl("CD4", cluster_all$majorCluster),], x = "majorCluster", 
          y = "expa", color = "sample_type", 
          add = "jitter") + 
  stat_compare_means(aes(group = sample_type), 
                     label =  "p.format", method = "wilcox.test", label.y = 0.08) +
  scale_color_manual(values = color_map_s2$cell_color)  +
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank()) + 
  labs(y = "STARTRAC-expansition")
ggsave("CD4-expansition.pdf", width = 7, height = 3)


ggboxplot(cluster_tumor[grepl("CD8", cluster_tumor$majorCluster),], x = "sample_type", 
          y = "expa", color = "sample_type", 
          add = "jitter") + facet_grid(~majorCluster)+
  stat_compare_means(aes(group = sample_type), 
                     label =  "p.format", method = "t.test") +
  scale_color_manual(values = color_map_s$cell_color) +
  theme(axis.title.x = element_blank()
#        axis.text.x = element_blank()
        ) + labs(y = "STARTRAC???trasition")


migration <- out@cluster.data
migration <- migration[grepl("patient", migration$aid) & !grepl("_2|_4", migration$aid),]
migration$er <- plyr::mapvalues(aes)

add_col <- colnames(add_meta)[-1]
for(i in 1:length(add_col)){
  cluster_tumor[,add_col[i]] <- plyr::mapvalues(cluster_tumor$aid, 
                                       as.vector(add_meta$patient_id),
                                       as.vector(add_meta[,add_col[i]]))
}


ggboxplot(migration, x = "majorCluster", 
          y = "migr", add = "jitter", color = "Grade") +
  scale_color_manual(values = color_map_s$cell_color)
