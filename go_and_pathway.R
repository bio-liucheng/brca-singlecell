library(ReactomePA)
library(GSVAdata)
data(c2BroadSets)
library(dplyr)
library(clusterProfiler)

gene_list <- seurat.obj.md.maker_1 %>% filter(cluster == "P", avg_logFC >0.5)
gene_list <- gene_list$gene

#tansfer gene id use bitr

gene_list <- bitr(gene_list, fromType = "SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")

file_name = "group_P"
if(!dir.exists(file_name)){
  dir.create(file_name)
}

ego <- enrichGO(gene         = gene_list$SYMBOL,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                ont = "ALL",
                 pvalueCutoff  = 1,
                 qvalueCutoff  = 0.05)
head(summary(ego))
write.csv(ego@result, file.path(file_name, "GO_term.csv"))
dotplot(ego, showCategory=30)
ggsave(file.path(file_name, "GO_term_dotplot.pdf"), width = 10, height = 8)
emapplot(ego)
ggsave(file.path(file_name, "GO_term_emapplot.pdf"), width = 8, height = 8)



x <- enrichPathway(gene=gene_list$ENTREZID,pvalueCutoff=0.1, readable=T)
head(as.data.frame(x))
write.csv(x@result, file.path(file_name, "pathway.csv"))

dotplot(x, showCategory=15)
ggsave(file.path(file_name, "pathway_dotplot.pdf"), width = 10, height = 6)


emapplot(x)
ggsave(file.path(file_name, "pathway_emapplot.pdf"), width = 9, height = 8)



cnetplot(x, categorySize="pvalue", foldChange=geneList)

library(ggplot2)
theme(legend.title = element_blank(), legend.text = element_text(size = 5))


#

