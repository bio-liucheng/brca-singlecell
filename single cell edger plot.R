#»ðÉ½Í¼
library(ggplot2)
library(ggrepel)


#EDGR ¼ÆËã
i = 1
for(i in c(1,4,5)){


tt <- tt_list[[i]]
df <- tt$table
df$patient <- tt$patient
df$change <- ifelse(abs(df$logFC) > 0.5 & df$FDR < 0.05, 
                    "Significant", 
                    "Not Significant")
df$up_in <- ifelse(df$logFC > 0.5 & df$FDR < 0.05, 
                   "Up", 
                   ifelse(df$logFC < -0.5 & df$FDR < 0.05, "Down", "Nochange"))

df$gene_name <- rownames(df)
ggplot(df,aes(x =-logFC,y= -log10(FDR)))+
  geom_point(aes(colour=as.factor(up_in)), size = 0.6)+
  scale_color_manual(values =c("#7cae00", "grey", "#f8766d"))+
  xlab('log2(Fold Change)')+
  ylab('-log10(pvalue)')+
  theme_classic() +  
  geom_text_repel(data = df[1:30,], 
       aes(x =-logFC,y= -log10(FDR),
       label = rownames(df[1:30,])), size = 3.5) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size =13),
        axis.title = element_text(size =13),
        axis.text = element_text(size =10)
        )+
  guides(colour = guide_legend(override.aes = list(size = 3)))

ggsave("CD4-CXCL13_cluster2 vs 3_diffgene.pdf", width = 6, height = 3.8)


if(i == 1)
  df_all <- df
else
  df_all <- rbind(df_all, df)
ggsave(paste0(tt$patient, ".pdf"))
}


library(dplyr)
library(reshape2)
library(UpSetR)
results <- df_all %>% group_by(up_in, gene_name) %>% summarise(a = n())

up_ln <- filter(df_all, up_in == "Up") %>% select(gene_name, patient, change)
res <- melt(up_ln, id.vars = c("gene_name", "patient"))
res <- dcast(res, gene_name~patient, margins = T)
rownames(res) <- res$gene_name
res <- res[-nrow(res),-ncol(res)]
upset(res)

ggsave("patient_2_4_two_distingted_low.pdf")

up_tumor<- filter(df_all, up_in == "Down") %>% select(gene_name, patient, change)
res <- melt(up_tumor, id.vars = c("gene_name", "patient"))
res <- dcast(res, gene_name~patient, margins = T)
rownames(res) <- res$gene_name
res <- res[-nrow(res),-ncol(res)]
upset(res)

flag <- res$patient_7 +res$patient_8 == 2
