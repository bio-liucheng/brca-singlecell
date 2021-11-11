tcr_files <- list.files()
for(i in 1:length(tcr_files)){
  tcr <- read.csv(file.path(tcr_files[i], "filtered_contig_annotations.csv"))
  tcr$colon_id <- paste(tcr_files[i], tcr$clonotype_id, sep = "_")
  tcr$sample_type <- substr(tcr_files[i], 1,2)
  tcr$sample_type <- ifelse(tcr$sample_type == "CA", "Tumor", "Lymph Node")
  tcr$patient_id <- paste("patient", substr(tcr_files[i], 3,3), sep = "_")
  if(i ==1)
    tcr_all <- tcr
  else
    tcr_all <- rbind(tcr_all, tcr)
}

tcr_all$colon_id <- paste(tcr_all$colon_id,tcr_all$raw_clonotype_id, sep = '')
tcr_all$contig_id <- as.character(tcr_all$contig_id)
tig_id <- strsplit(tcr_all$contig_id, "_contig_")
tig_id <- unlist(lapply(tig_id, "[", 2))
tig_id <- as.integer(tig_id)
tcr_all$contig_id_every_cell <- tig_id
s_all$colon_num <- plyr::mapvalues(s_all$colon_id, 
                                     as.vector(trc_colon$colon_id), trc_colon$frequency)
library(scales)
#fig tcr.1
tcr_all$percent_tc <- ifelse(tcr_all$frequency ==1, "n=1", "n¡Ý2")

level_percent_tc <- c("n¡Ý2","n=1")
ggplot(tcr_all, aes(patient_id)) + 
  geom_bar(aes(fill = factor(percent_tc, level_percent_tc)), position = "fill") + 
  facet_wrap(~factor(sample_type, levels = level_sample_type))+ theme_classic() +
  scale_fill_manual(values = pal_npg("nrc", 0.8)(10)) +
  theme(axis.ticks = element_blank(),
        axis.title.x =element_blank(),
        legend.title =element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1)
  )+
  scale_y_continuous(labels = percent)


#fig tcr.2
ggplot(subset(s_all, colon_num >1), aes(x = colon_num)) + 
  geom_histogram(binwidth = 10, aes(fill = TR_type)) + 
  facet_grid(~sample_type) + theme_classic()

library(ggsci)
#fig tcr.2.2
level_tcr <- c("single TRA", "single TRB", "Double TRB", "Double TRA",
            "Two TRA One TRB", "One TRA Two TRB","TRA TRB")
level_sample_type <- c("Tumor", "Lymph Node")
ggplot(subset(s_all, colon_num >1), aes(x = colon_num)) + 
  geom_histogram(binwidth = 10, aes(fill = factor(TR_type, levels = level_tcr))) + 
  facet_wrap(~factor(sample_type, level_sample_type)) + theme_classic() + 
  scale_fill_manual(values = pal_npg("nrc", 0.8)(10)) +
  theme(axis.ticks = element_blank(),
      axis.title.x =element_blank(),
      legend.title =element_blank(),
)
#fig tcr.2.3

level_tcr <- c("single TRA", "single TRB", "Double TRA", "Double TRB",
               "Two TRA One TRB", "One TRA Two TRB","TRA TRB")
level_sample_type <- c("Tumor", "Lymph Node")
ggplot(subset(s_all, colon_num ==1), aes(x = factor(sample_type, level_sample_type))) + 
  geom_bar(aes(fill = factor(TR_type, levels = level_tcr)), position = "fill", width = 0.3,) + 
  theme_classic() + 
  scale_fill_manual(values = pal_npg("nrc", 0.8)(10)) +
  theme(axis.ticks = element_blank(),
        axis.title.x =element_blank(),
        legend.title =element_blank(),
  )+scale_y_continuous(labels = percent)
