
#presentation
genes = c("TRBV", "TRBJ", "TRBC", "TRAV", "TRAJ")
names_sample = sort(unique(matrix_melt$orig.ident))

for(i in 1:length(genes)){
  gene = genes[i]
  for(j in 1:length(names_sample)){
sample_type_name = names_sample[j]

matrix_melt_f <- matrix_melt %>% filter(orig.ident == sample_type_name)
matrix_dcast = dcast(matrix_melt_f,cell_label~variable,sum)


matrix_dcast <- matrix_dcast[flag,]
rownames(matrix_dcast) <- matrix_dcast$cell_label

matrix_dcast <- matrix_dcast[,-c(1)]

colnames <- grepl(gene, colnames(matrix_dcast))

matrix_select <- matrix_dcast[,colnames]

matrix_select <- matrix_select[,sort(colnames(matrix_select))]
matrix_select <- apply(matrix_select, 1, function(x) x/sum(x))

pheatmap(t(matrix_select), cluster_rows = F, cluster_cols = F,
         color =  viridis(100),cellwidth = 20,cellheight = 20,  
         border_color = NA, filename = paste0("cell_label_", gene,"_",sample_type_name, ".pdf"))
}
}



#presentation
genes = c("TRBV", "TRBJ", "TRBC", "TRAV", "TRAJ")
names_sample = sort(unique(matrix_melt$orig.ident))

for(i in 1:length(genes)){
  gene = genes[i]
    matrix_dcast = dcast(matrix_melt,orig.ident~variable,sum)
    
    
    rownames(matrix_dcast) <- matrix_dcast[,1]
    
    matrix_dcast <- matrix_dcast[,-c(1)]
    
    colnames <- grepl(gene, colnames(matrix_dcast))
    
    matrix_select <- matrix_dcast[,colnames]
    matrix_select <- matrix_select[,sort(colnames(matrix_select))]
    matrix_select <- apply(matrix_select, 1, function(x) x/sum(x))
    colnames(matrix_select) <- rownames(matrix_dcast)
    matrix_select <- as.matrix(matrix_select)
    matrix_select <- matrix_select[, !colnames(matrix_select) %in% c("CA2", "CA4")]
    pheatmap(t(matrix_select), cluster_rows = F, cluster_cols = F,
             color =  viridis(100),cellwidth = 20,cellheight = 20,  
             border_color = NA, filename = paste0("sample_", gene,"_",sample_type_name, ".pdf"))
}
