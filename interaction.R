file_list = c("tec_caf", "T_cell", "myeloid")
for(i in 1:length(file_list)){
  meta <- read.csv(paste0(file_list[i], "/meta.csv"))
  meta_s <- data.frame(cell = meta$X, cell_label = meta$cell_label)
  if(i == 1)
    meta_all = meta_s
  else
    meta_all = rbind(meta_all, meta_s)
}

meta_all <- read.csv("meta.csv")
meta_all <- meta_all[grepl("CD8-|CD4-|NK-|Macro|DC-|CAF-C1|CAF-C2|CAF-C3|TEC|B", meta_all$cell_label),]
meta_all <- meta_all[!grepl("un", meta_all$cell_label),]


meta_all$cell_label[grepl("CD4", meta_all$cell_label)] <- "CD4"
meta_all$cell_label[grepl("CD8", meta_all$cell_label) & !grepl("CD8-C5-CXCL13", meta_all$cell_label)] <- "CD8"
meta_all$cell_label[grepl("NK", meta_all$cell_label)] <- "NK"
meta_all$cell_label[grepl("CD8-C5-CXCL13", meta_all$cell_label)] <- "CD8-EX"
meta_all$cell_label[grepl("Macro", meta_all$cell_label)& !grepl("Macro-C4-CXCL11", meta_all$cell_label)] <- "Macro"
meta_all$cell_label[grepl("Macro-C4-CXCL1", meta_all$cell_label)] <- "Macro-M1"
meta_all$cell_label[grepl("DC", meta_all$cell_label)] <- "DC"
meta_all$cell_label[grepl("TEC", meta_all$cell_label)] <- "TEC"

meta_all <- meta_all[!grepl("plasma", meta_all$cell_label),]
meta_all$cell_label <- gsub("-", "_", meta_all$cell_label)
cell_label = unique(meta_all$cell_label)
index_all = c()
for(i in 1:length(cell_label)){
  index = which(meta_all$cell_label %in% cell_label[i])
  if(length(index) <1000){
    index_all = c(index_all, index)
  }
  else
  {
    index_all <- c(index_all, sample(index, 1000, replace = F))
  }
}


meta_all_s <- meta_all[sort(index_all),]
seurat.obj.md <- seurat.obj[,meta_all_s$X]
seurat.obj.md$cell_label <- meta_all_s$cell_label

write.table(data.frame(Cell = Cells(seurat.obj.md), cell_type = seurat.obj.md$cell_label), "three_interation_annoation.txt", sep = '\t', quote = F)
write.table(as.matrix(seurat.obj.md[["RNA"]]@counts), "three_interation_matrix.txt", sep = '\t', quote = F)
#calculate

p <- read.csv("pvalues.csv")
mean <- read.csv("significant_means.csv")
p <- p[grepl("simple", p$partner_a),]
mean <- mean[grepl("simple", mean$partner_a),]

p <- p[!duplicated(p$interacting_pair),]
mean <- mean[!duplicated(mean$interacting_pair),]

rownames(p) <- p$interacting_pair
rownames(mean) <-mean$interacting_pair

inter_name <- intersect(colnames(p), colnames(mean))
inter_pair <- intersect(rownames(p), rownames(mean))

p <- p[inter_pair,]
mean <- mean[inter_pair,]

p <- p[,inter_name]
mean <- mean[, inter_name]

p <- p[, -c(1:11)]
mean <- mean[, -c(1:11)]

mean <- as.matrix(mean)
p <- as.matrix(p)

mean[is.na(mean)] <- 0
p[is.na(p)] <- 1




library(reshape2)
mean_melt <- melt(mean)
p_melt <- melt(p)

mean_melt$value <- as.vector(mean_melt$value)
p_melt$value <- as.vector(p_melt$value)

mean_melt$value <- as.numeric(mean_melt$value)
p_melt$value <- as.numeric(p_melt$value)

table(mean_melt$value >0.5 & p_melt$value <0.1)

flag <- mean_melt$value >0.5 & p_melt$value <0.1

select <- mean_melt[flag,]

cluster <- as.vector(select$Var2)
cluster <- strsplit(cluster, '.', fixed = T)
cluster_1 <- unlist(lapply(cluster, "[", 1))
cluster_2 <- unlist(lapply(cluster, "[", 2))

interactor <- as.vector(select$Var1)
interactor <- strsplit(interactor, '_', fixed = T)
interactor_1 <- unlist(lapply(interactor, "[", 1))
interactor_2 <- unlist(lapply(interactor, "[", 2))


flag_deg <- interactor_1 %in% deg$gene | interactor_2 %in% deg$gene
select_deg <- select[flag_deg,]
flag_select <- grepl("CAF", select_deg$Var2) & grepl("Macro|DC|CD", select_deg$Var2)
select_deg <- select_deg[flag_select,]
select_deg$value <- log2(select_deg$value+0.1)
select_deg$value <- select_deg$value + abs(select_deg$value)




#select interaction and plot
select_caf <- grepl("FN1_a4b1|OGN|PLA2G2A_a4b1|VSIR|DPP4", select_deg$Var1)
select_caf <- select_deg[select_caf,]
cluster <- as.vector(select_caf$Var2)
cluster <- strsplit(cluster, '.', fixed = T)
cluster_1 <- unlist(lapply(cluster, "[", 1))
cluster_2 <- unlist(lapply(cluster, "[", 2))

interactor <- as.vector(select_caf$Var1)
interactor <- strsplit(interactor, '_', fixed = T)
interactor_1 <- unlist(lapply(interactor, "[", 1))
interactor_2 <- unlist(lapply(interactor, "[", 2))

flag <- grepl("CAF", cluster_1)
cluster_tmp <- cluster_1
cluster_1[!flag] <- cluster_2[!flag]
cluster_2[!flag] <- cluster_tmp[!flag]

interactor_tmp <- interactor_1
interactor_1[!flag] <- interactor_2[!flag]
interactor_2[!flag] <- interactor_tmp[!flag]

select_caf$Var2 <- paste(cluster_1, cluster_2, sep = ":")
select_caf$Var1 <- paste(interactor_1, interactor_2, sep = ":")

library(reshape2)

dcast_caf <- dcast(select_caf,Var1~Var2)
rownames(dcast_caf) <- dcast_caf$Var1
dcast_caf <- dcast_caf[,-1]
dcast_caf <- dcast_caf[,grepl("Macro|DC|CD", colnames(dcast_caf))]
dcast_caf[is.na(dcast_caf)] <- 0

library(pheatmap)

pheatmap(dcast_caf, cluster_rows = F, cluster_cols = F)
library(ggplot2)
caf_immune <- melt(as.matrix(dcast_caf))

select_caf$value <- log2(select_caf$value +0.1)
select_caf$value <- select_caf$value + abs(min(select_caf$value))

cluster <- strsplit(as.vector(caf_immune$Var2),":")
cluster <- unlist(lapply(cluster, "[", 2))
caf_immune$Var2 <- cluster
ggplot(caf_immune,aes(x=Var1,y=Var2))+
  geom_vline(aes(xintercept = Var1), colour = "grey", linetype = "dashed")+
  geom_hline(aes(yintercept = Var2), colour = "grey", linetype = "dashed")+
  geom_point(aes(color=value, size= value))+
  coord_flip()+
  scale_color_gradient2(low = muted("blue"),mid = "white", high = muted("red"), midpoint = 6) +
  facet_grid(~cluster, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
  )

ggsave("caf_immune.pdf", width = 6.5, height = 4)

col_name <- colnames(dcast_caf)
col_name <- strsplit(col_name, '.', fixed = T)
col_name_1 <- unlist(lapply(col_name, "[", 1))
col_name_2 <- unlist(lapply(col_name, "[", 2))


