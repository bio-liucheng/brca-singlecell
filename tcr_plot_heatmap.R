meta <- read.csv("meta_tcr_count.csv")
rownames(meta) <- meta$X

matrix <- meta[,grepl('vdj_v_|vdj_j_|vdj_constant_', colnames(meta))]
matrix <- cbind(matrix, meta[,c("cell_label", "sample_type", "orig.ident")])

library(dplyr)
library(reshape2)

matrix_melt <- melt(matrix, id = c("cell_label","sample_type","orig.ident"))

matrix_dcast = dcast(matrix_melt,orig.ident+cell_label+sample_type~variable,sum)

flag <- !matrix_dcast$orig.ident %in% c("CA2", "CA4") & 
  !grepl("NK-",matrix_dcast$cell_label)

matrix_dcast <- matrix_dcast[flag,]
rownames(matrix_dcast) <- paste(matrix_dcast$orig.ident, matrix_dcast$cell_label, sep = '_')

matrix_dcast <- matrix_dcast[,-c(1,2)]


matrix_dcast <- matrix_dcast[,apply(matrix_dcast, 2, function(x) sum(x>0)) >5]

matrix_dcast_norm <- apply(matrix_dcast, 1, function(x){x/sum(x)*10})

matrix_dcast_norm <- apply(matrix_dcast_norm, 1, scale)

rownames(matrix_dcast_norm) <- rownames(matrix_dcast)



matrix_dcast_norm <- cbind(matrix_dcast_norm, matrix_dcast[,c(1,2,3)])


matrix_dcast_norm_melt <- melt(matrix_dcast_norm, id = c("cell_label","sample_type","orig.ident"))

gene <- grepl("TRA", matrix_dcast_norm_melt$variable)
ggplot(matrix_dcast_norm_melt[gene,], aes(x = variable, y = value, fill = cell_label))+ geom_boxplot()


matrix_dcast_norm <- scale(matrix_dcast_norm)

mydata.pca<-prcomp(t(matrix_dcast_norm),scale. = F) 


ggbiplot(mydata.pca, obs.scale = 1, var.scale = 1, var.axes = F, labels =colnames(matrix_dcast_norm))


#presentation
gene = "TRBV"
sample_type = ""

matrix_dcast = dcast(matrix_melt,sample_type~variable,sum)

flag <-!grepl("NK-",matrix_dcast$orig.ident)

matrix_dcast <- matrix_dcast[flag,]
rownames(matrix_dcast) <- matrix_dcast$sample_type

matrix_dcast <- matrix_dcast[,-c(1)]

colnames <- grepl(gene, colnames(matrix_dcast))

matrix_select <- matrix_dcast[,colnames]

matrix_select <- matrix_select[,sort(colnames(matrix_select))]
matrix_select <- apply(matrix_select, 1, function(x) x/sum(x))

pheatmap(t(matrix_select), cluster_rows = F, cluster_cols = F,
         color =  viridis(100),cellwidth = 20,cellheight = 20, 
         border_color = NA, filename = paste0("cell_label_", gene, ".pdf"))
pheatmap(t(matrix_select), cluster_rows = F, cluster_cols = F,
         color =  viridis(100),cellwidth = 20,cellheight = 20,  filename = "test.pdf")

