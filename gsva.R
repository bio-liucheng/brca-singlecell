library(GSVA)
library(Seurat)
library(GSVAdata)
library(GSEABase)
library(clusterProfiler)
data(c2BroadSets)

seurat.obj.md <- seurat.obj[,grepl("Macro", seurat.obj$cell_label)]

meta <- meta[grepl("CAF", meta$cell_label),]

seurat.obj.md <- seurat.obj[,meta$X]

seurat.obj.md$cell_label <- meta$cell_label

Idents(seurat.obj.md) <- "sample_type"

avera_cluster <- AverageExpression(seurat_object_edgeR)
avera_cluster <- avera_cluster$RNA
nGS <- 100    ## number of gene sets
min.sz <- 10  ## minimum gene set size
max.sz <- 100 ## maximum gene set size

#gs <- lapply(gs, function(n, p) sample(1:p, size=n, replace=FALSE), p) ## sample gene sets




#id transform
k=keys(org.Hs.eg.db,keytype = "ENSEMBL")

list=select(org.Hs.eg.db,keys=k,columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
ID_list=list[match(rownames(avera_cluster),list[,"SYMBOL"]),]
ID_list = na.omit(ID_list)

avera_cluster <- avera_cluster[ID_list$SYMBOL,]
rownames(avera_cluster) <- ID_list$ENTREZID


canonicalC2BroadSets <- c2BroadSets[c(grep("PATH", names(c2BroadSets)))]

#filter gene
flag <- apply(avera_cluster, 1, function(x) sum(x >0) >0.1 *ncol(avera_cluster))

#all enrich
gsva_es <- gsva(as.matrix(avera_cluster[flag,]), c2BroadSets, mx.diff=1)

# only pathway
gsva_es <- gsva(as.matrix(avera_cluster[flag,]), canonicalC2BroadSets, mx.diff=1)

mad_o <- apply(gsva_es, 1, mad)
library(pheatmap)
pheatmap(gsva_es[order(mad_o, decreasing = T)[1:15],], 
         show_rownames = T, cluster_rows = T, cluster_cols = F, cellwidth = 25, cellheight = 15,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))



#vlnplot intresest pathway

gene_set = "NIKOLSKY_BREAST_CANCER_1Q21_AMPLICON"
ascor <- c2BroadSets[c(grep(gene_set, names(c2BroadSets)))]
ascor <- ascor[[1]]
ascor <- ascor@geneIds
ascor <- bitr(ascor, fromType = "ENTREZID", toType=c("SYMBOL","ENSEMBL"), OrgDb="org.Hs.eg.db")

sig_select <- ascor$SYMBOL
expr_scale <- seurat.obj.md[["RNA"]]@data
expr_scale <- expr_scale[intersect(sig_select, rownames(expr_scale)),]
sig_cell <- apply(expr_scale, 2, mean)

sig_frame <- data.frame(sig_score = sig_cell, cell_names = colnames(seurat.obj.md),
  cell_label = seurat.obj.md$cell_label, sample_type = seurat.obj.md$sample_type)

library(ggplot2)
library(ggpubr)
ggplot(sig_frame, aes(sample_type, sig_score,  group = sample_type)) +
  geom_boxplot(width = 0.7) + geom_jitter(width = 0.2, size = 0.3) + theme_classic() +
  theme(legend.position = NULL) +
  stat_compare_means(method = "wilcox.test", ref.group = "Macro-C3", label = "p.format")


seurat.obj.md.er <- subset(seurat.obj.md, patient_id %in% c("patient_8"))

##pathway analysis use single cell
library(clusterProfiler)
seurat.obj.md <- NormalizeData(seurat.obj.md)
cell_matrix <- seurat.obj.md[["RNA"]]@data
gene_list <- bitr(rownames(cell_matrix), fromType = "SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
gene_list = na.omit(gene_list)

library(Matrix)

cell_matrix <- as.matrix(cell_matrix)
cell_matrix <- cell_matrix[gene_list$SYMBOL,]
rownames(cell_matrix) <- gene_list$ENTREZID

#filter gene
flag <- apply(cell_matrix, 1, function(x) sum(x >0 ) >0.1 *ncol(cell_matrix))

gsva_es <- gsva(as.matrix(cell_matrix[flag,]), canonicalC2BroadSets, mx.diff=1)


flag <- seurat.obj.md$Malign.type == "malignant" &(seurat.obj.md$monocle_state %in% c(1))
gsva_es_mag <- gsva_es[,flag]
state <- seurat.obj.md$sample_type[flag]
library(limma)
adjPvalueCutoff <- 0.001
logFCcutoff <- log2(2)
design <- model.matrix(~ factor(state))
design <- model.matrix(~ factor(ifelse(meta$cell_label == "CAF-C3-PLA2G2A", "NEG", "PO")))
colnames(design) <- c("ALL", "C4")
fit <- lmFit(gsva_es, design)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef="C4", number=Inf)
DEgeneSets <- topTable(fit, coef="C4", number=Inf,
                       p.value=adjPvalueCutoff, adjust="BH")
res <- decideTests(fit, p.value=adjPvalueCutoff)
summary(res)

n = 20
library(dplyr)
DEgeneSets$row_name <- rownames(DEgeneSets)
DEgeneSets <- DEgeneSets %>% arrange(adj.P.Val)
annotation_col <- data.frame(
                            cell_label = seurat.obj.md$cell_label, 
                            sample_type = seurat.obj.md$sample_type)
rownames(annotation_col) <- colnames(gsva_es)
pheatmap(gsva_es[DEgeneSets$row_name[1:30],], show_rownames = T, show_colnames = F, 
                annotation_col = annotation_col, cluster_cols = T)
