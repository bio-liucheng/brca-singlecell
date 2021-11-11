meta <- meta[!grepl("_2|_4", meta$patient_id),]

meta$Grade <- paste0("Grade_", meta$Grade)
library(dplyr)
library(ggplot2)
library(ggpubr)
patient_ncells <- meta  %>% group_by(orig.ident)%>% dplyr::summarise(n = n())

patient_nlabel <- meta  %>% group_by(orig.ident, cell_label, patient_id, sample_type) %>% dplyr::summarise(n = n())

patient_nlabel$ncell <- plyr::mapvalues(patient_nlabel$orig.ident, 
                                        as.vector(patient_ncells$orig.ident),
                                        patient_ncells$n)

patient_nlabel$ncell <- as.integer(patient_nlabel$ncell)
patient_nlabel <- mutate(patient_nlabel, percent  = n/ncell)

ggplot(sp, aes(x=Grade, y=percent)) + 
  geom_bar(position=position_dodge(), stat="identity", width = 0.3, fill = "grey") +
  geom_errorbar(aes(ymin=percent-se, ymax=percent+se),
                width=.1,                    # Width of the error bars
                position=position_dodge(.9)) +
  facet_wrap(~cell_label, ncol = 4, scales = "free") +
  geom_point(data = patient_nlabel, aes(x = Grade, y = percent,color = patient_id))  +
  theme_classic()+
  stat_compare_means(data = patient_nlabel,  method = "t.test")


ggplot(data = patient_nlabel, aes(x = factor(sample_type, levels = c("Tumor", "Lymph Node")), y = percent,color = Grade, group = patient_id))+
geom_point()  + geom_line() +
  facet_wrap(~cell_label, ncol = 4, scales = "free")+
  theme_classic()+
  stat_compare_means(method = "t.test")

ggplot(data = cxcl13, aes(x = factor(sample_type, levels = c("Tumor", "Lymph Node")), y = percent, group = patient_id))+
  geom_point()  + geom_line() +
  theme_classic()+
  stat_compare_means(paired = T, comparisons = )

ggplot(data = cd8b, aes(x = factor(sample_type, levels = c("Tumor", "Lymph Node")), y = percent, group = patient_id))+
  geom_point(aes(color = patient_id), size = 3)  + geom_line() + 
  scale_colour_manual(values = color_map_s$cell_color) + theme_classic() +
  stat_compare_means(paired = T, comparisons = my_comparision) + 
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        ) + labs(x = "", y = "% percent in T cells")
