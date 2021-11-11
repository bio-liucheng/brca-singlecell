
meta$Grade <- paste0("Grade_", meta$Grade)
library(dplyr)
library(Rmisc)
library(ggpubr)

meta <- read.csv("meta2.csv")

meta <- meta[!grepl("2|4", meta$orig.ident),]
meta <- meta[!grepl("CD-un", meta$cell_label),]
meta <- meta[grepl("CD4|CD8", meta$cell_label),]

meta$sample_type <- factor(meta$sample_type, levels = c("Tumor", "Lymph Node"))

patient_ncells <- meta %>% group_by(orig.ident)%>% dplyr::summarise(n = n())

patient_nlabel <- meta %>% group_by(patient_id,orig.ident, sample_type, cell_label) %>% dplyr::summarise(n = n())

patient_nlabel$ncell <- plyr::mapvalues(patient_nlabel$orig.ident, 
                                        as.vector(patient_ncells$orig.ident),
                                        as.vector(patient_ncells$n))

patient_nlabel$ncell <- as.vector(patient_nlabel$ncell)
patient_nlabel$ncell <- as.integer(patient_nlabel$ncell)
patient_nlabel <- mutate(patient_nlabel, percent  = n/ncell)
patient_nlabel <- patient_nlabel[!grepl("CD-un", patient_nlabel$cell_label),]


ggplot(patient_nlabel, aes(Grade, percent)) + 
  geom_boxplot() + facet_grid(~cell_label,  scales = "free")



sp <- patient_nlabel %>% filter(sample_type == "Lymph Node") %>% summarySE( measurevar = "percent", groupvars = c("ER_statu", "cell_label"))
patient_nlabel_s  <- patient_nlabel %>% filter(sample_type == "Lymph Node")

sp <- patient_nlabel %>% filter(!(orig.ident %in% c("CA2", "CA4") & cell_label %in% c("CD4+ T", "CD8+ T"))) %>% summarySE( measurevar = "percent", groupvars = c("sample_type", "cell_label"))
patient_nlabel_s  <- patient_nlabel %>%  filter(!(orig.ident %in% c("CA2", "CA4") & cell_label %in% c("CD4+ T", "CD8+ T")))

patient_nlabel<- patient_nlabel[grepl("B|CD4|Epithelial Cells|Mast cells", patient_nlabel$main_label),]
patient_nlabel<- patient_nlabel[!grepl("plasma", patient_nlabel$main_label),]

sp <- patient_nlabel  %>% summarySE(measurevar = "percent", groupvars = c("sample_type",  "cell_label"))
patient_nlabel_s  <- patient_nlabel 



ggplot(sp, aes(x=sample_type, y=percent)) + 
  geom_bar(position=position_dodge(), stat="identity", width = 0.5, fill = "grey") +
  geom_errorbar(aes(ymin=percent-se, ymax=percent+se),
                width=.1,                    # Width of the error bars
                position=position_dodge(.9)) +
  facet_wrap(~cell_label,nrow = 2,  scales = "free") +
  geom_point(data = patient_nlabel_s, 
             aes(x = sample_type, y = percent, color = factor(patient_id, levels = color_map_s$cell_label)))  +
  theme_classic() + scale_color_manual(values = color_map_s$cell_color) +
  theme(legend.title = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 11),
        legend.text = element_text(size = 11),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.x = element_blank()) + labs(y = "% of all cells")+
  stat_compare_means(data = patient_nlabel_s, na.rm = T,
                     aes(group = sample_type),
                     label = "p.format", method = "t.test") +
  guides(colour = guide_legend(override.aes = list(size = 3)))

ggsave("cell_label_patient percent.pdf", width = 9, height = 5.5)

color_map_s <- color_map[color_map$cell_label %in% unique(patient_nlabel_s$patient_id), ]


ggplot(sp, aes(x=factor(sample_type, levels = c("Tumor", "Lymph Node")), y=percent)) + 
  geom_bar(position=position_dodge(), stat="identity", width = 0.3, fill = "grey") +
  geom_errorbar(aes(ymin=percent-se, ymax=percent+se),
                width=.1,                    # Width of the error bars
                position=position_dodge(.9)) +
  facet_wrap(~cell_label,nrow = 2,  scales = "free") +
  geom_point(data = patient_nlabel_s, 
             aes(x = factor(sample_type, levels = c("Tumor", "Lymph Node")), y = percent,color = factor(patient_id, levels = color_map_s$cell_label)))  +
  theme_classic() + scale_color_manual(values = color_map_s$cell_color) +
  theme(legend.title = element_blank(), strip.background = element_blank(), 
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.x = element_blank()) + labs(y = "% of myeloid cells")+geom_line(data = patient_nlabel_s,aes(group = patient_id))+
  stat_compare_means(data = patient_nlabel_s, 
                     aes(group = factor(sample_type, levels = c("Tumor", "Lymph Node"))),
                     method = "wilcox.test", label = "p.format", paired = T) +
  guides(colour = guide_legend(override.aes = list(size = 3)))
