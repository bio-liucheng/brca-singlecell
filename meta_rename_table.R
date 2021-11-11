
metadata <- metadata[inter_name,]
s_all <- s_all[inter_name,]
metadata_tcr <- cbind(metadata, s_all)



meta_data_tcr_type <- subset(metadata_tcr, sample_type != "Tumor")

s <- table(meta$leiden2, meta$new_label)
library(reshape2)
s <- as.data.frame(s)
s <- dcast(s, Var1~Var2)
rownames(s) <- s$Var1
s <- s[,-1]
t <- apply(s, 1, function(x){which(x == max(x))})
new_label <- colnames(s)[t]
meta$new_label <- plyr::mapvalues(as.vector(meta$leiden2), 
                              0:41,
                              as.vector(new_label)
                                        )
meta$new_label[meta$leiden2 == 39] <- "Myeloid cells"
meta$new_label[meta$new_label == "CD T"] <- "B"

s <- s[grepl("CD", rownames(s)),]
sum <- apply(s, 1, function(x) {x/sum(x)})

pheatmap(sum, cluster_rows = F, cluster_cols = F, scale = "row")

['CD4-C1', 'CD8-C2', 'CD4-C2', 'CD4-C3', 'CD4-C2', 'CD8-C3', 
  'CD4-C5', 'CD4-C6', 'CD8-C4', 'CD4-C4', 'CD4-C4', 'CD4-C5', 
  'CD8-C2', 'CD4-C2', 'CD4-C2', 'NK-C1', 'CD8-C1', 'CD8-C4', 'NK-C3', 
  'CD4-C4', 'CD8-C5', 'NK-C2']