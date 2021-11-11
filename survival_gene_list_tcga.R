library(survival)
library(survminer)
library(survMisc)


gene_list = c("GPR157")
setwd("G:/scRNA-seq/LC/TCGA/BRCA")
options(stringsAsFactors = F)
clin <- read.delim("BRCA_clinicalMatrix")
expr <- read.delim("HiSeqV2.gz")
surv <- read.delim("BRCA_survival.txt.gz")
rownames(expr) <- expr[,1]
expr <- expr[,-1]
colnames(expr) <- gsub(".", '-', colnames(expr), fixed = T)
rownames(clin) <- clin$sampleID
#¼ÆËãtumor vs normal ²îÒì

clin <- clin[colnames(expr),]
tumor_flag <- clin$sample_type == "Primary Tumor"

normol_flag <- clin$sample_type == "Solid Tissue Normal"




ER_positive <- clin$ER_Status_nature2012 == "Positive"

#caculate survival HR
clin <- clin[tumor_flag,]
expr <- expr[,tumor_flag]
rownames(surv) <- surv$sample
surv <- surv[,-1]
surv <- surv[rownames(clin),]
clin <- cbind(clin, surv)




#filter clin
clin2 <- clin[!grepl("Her2|Lum", clin$PAM50Call_RNAseq),]
expr2 <- expr[,rownames(clin2)]

data <- as.matrix(expr)

surv_genelist <- function(data, clin_s, gene_list){

  inter_gene <- intersect(rownames(data), gene_list[1:20])
  expr_gene <- data[inter_gene,]
  expr_gene <- apply(expr_gene, 1, scale)
  rownames(expr_gene) <- colnames(data)
  expr_gene <- t(expr_gene)
  expr_gene <- apply(expr_gene, 2, sum)
  flag <- ifelse(expr_gene >median(expr_gene, na.rm = T), "high", "low")
  flag2 <- expr_gene < quantile(expr_gene, probs = 0.4)
  flag3 <- expr_gene > quantile(expr_gene, probs = 0.6)

  
  fit <- coxph(Surv(OS.time, OS) ~ flag ,subset = flag2 | flag3,  data = clin_s)
  sum_fit <- summary(fit)
  fit2 <- survfit(Surv(OS.time, OS) ~ flag, subset = flag2 | flag3, data = clin_s)
  ggsurvplot(fit2, data = clin_s, pval = TRUE, )
  sum_fit
}
gene_list


#calculate survival with gene_list average expression gene
surv_genelist <- function(data, clin, clin_type, sub_type, gene_list){
  if(is.null(clin_type)){
    sub_expr <- data
    sub_clin <- clin
  }else{
    flag_sub <- clin[,clin_type] %in% sub_type
    sub_clin <- clin[flag_sub,]
    sub_expr <- data[,flag_sub]
  }
  if(length(gene_list) >1){
    inter_gene <- intersect(rownames(sub_expr), gene_list)
    expr_gene <- sub_expr[inter_gene,]
    expr_gene <- apply(expr_gene, 2, mean)
  }else{
    expr_gene <- sub_expr[gene_list,]
    expr_gene <- as.numeric(expr_gene)
  }
  
  flag <- ifelse(expr_gene >median(expr_gene, na.rm = T), "high", "low")
  flag <- factor(flag, levels = c("high", "low"))
  flag_low <- expr_gene < quantile(expr_gene, probs = 0.45, na.rm = T)
  flag_high <- expr_gene > quantile(expr_gene, probs = 0.55, na.rm = T)
  sub_clin$flag <- flag
  sub_clin$flag_low <- flag_low
  sub_clin$flag_high <- flag_high
#  fit <- coxph(Surv(time = OS.time, event = OS) ~ flag, subset = flag_low |flag_high, data = sub_clin)
  fit2 <- survfit(Surv(time = OS.time, event = OS) ~ flag, subset = flag_low |flag_high, data = sub_clin)
  ggsurvplot(fit2, data = sub_clin, pval = TRUE)
#  sum_fit <- summary(fit)
}


gene_list_a <- gene_list$gene
gene_list_b <- c("LAMP3", "C1DC", "CLEC9A")

gene_a <- apply(data[gene_list_a,], 1, scale)
gene_a <- t(gene_a)
gene_a <- apply(gene_a, 2, mean)

gene_a <- apply(data[intersect(gene_list_a, rownames(data)),], 2, mean)
gene_b <- apply(data[gene_list_b,], 2, mean)
  
gene_a <- as.numeric(data["LAMP3",])
gene_b <- as.numeric(data["LAMP3",])

relative_gene <- scale(gene_a) / gene_b

relative_gene <-  gene_b / gene_a

relative_gene <- gene_a

flag <- ifelse(relative_gene >median(relative_gene, na.rm = T), "high", "low")
flag <- factor(flag, levels = c("high", "low"))
flag_low <- relative_gene < quantile(relative_gene, probs = 0.45, na.rm = T)
flag_high <- relative_gene > quantile(relative_gene, probs = 0.55, na.rm = T)
clin$flag <- flag
clin$flag_low <- flag_low
clin$flag_high <- flag_high
fit2 <- survfit(Surv(time = OS.time, event = OS) ~ flag, subset = flag_low |flag_high, data = clin)
ggsurvplot(fit2, data = clin, pval = TRUE)


if(T){
  #get a gene list 
  gene_cluster <- "CAF_C3_PLA2G2A" 
  gene_list <- cluster_significant_markers[cluster_significant_markers$cluster == gene_cluster,]$gene
  if(!is.null(intersect(gene_list, rownames(expr))))
    fit <- surv_genelist(expr, clin, "ER_Status_nature2012", "Positive", gene_list)
  
}

type <- unique(clin$ER_Status_nature2012)
type <- type[-c(1,4)]
library(dplyr)
#calculate survival
for(j in 1:length(type)){
  type_s <- type[j]
  cluster <- unique(final_marker_s$cluster)
  for(i in 1:length(cluster)){
    gene_list <- final_marker_s$gene[final_marker_s$cluster == cluster[i]]
    f <- surv_genelist(expr, clin, NULL, type_s, gene_list)
    f_table <-data.frame(coef = f[["conf.int"]][1], 
                         low_95 = f[["conf.int"]][3], 
                         high_95 = f[["conf.int"]][4],
                         pvalue = f[["coefficients"]][5])
    if(i == 1)
      f_table_final <- f_table
    else
      f_table_final <- rbind(f_table_final, f_table)
    
  }
  f_table_final$cluster <- cluster
  f_table_final$patient_type <- type_s
  
  if(j ==1)
    f_table_final_final <- f_table_final
  else
    f_table_final_final <- rbind(f_table_final_final, f_table_final)
}




#calculate survival with two gene (high low) expression gene
surv_TwoGene <- function(data, clin_s, two_gene){
  
  inter_gene <- intersect(rownames(data), two_gene)
  expr_gene <- data[inter_gene,]
  expr_gene <- apply(expr_gene, 1, scale)
  rownames(expr_gene) <- colnames(data)
  gene1_hgih <- ifelse(expr_gene[,1] > median(expr_gene[,1]),paste0(two_gene[1], "_high"),paste0(two_gene[1], "low"))
  gene2_high <- ifelse(expr_gene[,2] > median(expr_gene[,2]),paste0(two_gene[2], "_high"),paste0(two_gene[2], "low"))
  
  s <- paste(gene1_hgih, gene2_high, sep = '_')
  s <- factor(s, levels = unique(s)[c(2,1,3,4)])
  
  clin_s$s <- s
  
  
  fit <- coxph(Surv(OS.time, OS) ~ s, data = clin_s)
  sum_fit <- summary(fit)
  fit2 <- survfit(Surv(OS.time, OS) ~  s, data = clin_s)
  ggsurvplot(fit2, data = clin_s, pval = TRUE)
}


surv_genelist(expr, clin, gene_list)
unique(clin$PAM50Call_RNAseq)
