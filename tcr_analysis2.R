data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df2 <- data_summary(matrix_dcast_norm_melt, varname="value", 
                    groupnames=c("cell_label", "variable", "sample_type"))

ggplot(df2, aes(x=variable, y=value, group=cell_label, color=cell_label)) + 
  geom_pointrange(aes(ymin=value-sd, ymax=value+sd))

clone_count <- meta$vdj_clone_count



sample = "vdj_clonotype_n_5"
tumor_ln <- "sample_type"
sample_type = sort(unique(meta$vdj_clonotype_n_5))
sample_type = sample_type[!grepl("NK-", sample_type)]


for(j in 1:length(genes)){

for(i in 1:length(sample_type)){
  

matrix <- meta[!duplicated(meta$vdj_clonotype),]
matrix <- matrix[matrix[,sample] %in% sample_type[i] ,]
#matrix <- matrix[matrix$vdj_clone_count <=5 ,]
#matrix <- matrix[matrix[,tumor_ln] %in% "Lymph Node", ]
matrix <- matrix[, grepl(genes[j], colnames(matrix))]
t2 <- apply(matrix, 2, function(x) sum(x>0))
t2 <- t2/sum(t2)
if(i == 1)
  t_all <- t2
else
  t_all <- rbind(t_all, t2)

  }
rownames(t_all) <- sample_type
pheatmap(t_all, cluster_rows = F, cluster_cols = F,
           color =  viridis(100),cellwidth = 20,cellheight = 20, border_color = NA, 
           filename = paste0("colonal_unique_clone_all",sample_type[i],genes[j],".pdf"))
  
}

genes = c("TRBV", "TRBJ", "TRBC", "TRAV", "TRAJ")
for(i in 1:length(genes)){
pheatmap(t_all[,grepl(genes[i], colnames(t_all))], cluster_rows = F, cluster_cols = F,
         color =  viridis(100),cellwidth = 20,cellheight = 20, border_color = NA, 
         filename = paste0("cell_label_unique_LN_",genes[i], ".pdf"))
}
