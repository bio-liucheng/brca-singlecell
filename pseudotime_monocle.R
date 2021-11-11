library(Seurat)
library(monocle)
options(stringsAsFactors = F)
meta <- read.csv("immue_cd4_meta.csv")
sub_name <- meta$X

seurat.obj.md <-seurat.obj[,meta$X[meta$cell_label %in% c("Macro-C3-VCAN", "Macro-C5-SPP1")]]
seurat.obj.md <-seurat.obj[,meta$X[meta$cell_label %in% c("Macro-C4-CXCL11","Macro-C3-VCAN", "Macro-C5-SPP1")]]

seurat.obj <- seurat.obj.immue[,sub_name]

cells <- Cells(seurat.obj)

cells_s <- sample(cells, size = 4000, replace = F)

seurat.obj.md <- seurat.obj[,cells_s]

raw_data <- seurat.obj.md[["RNA"]]@counts
print("Computing var genes by cell type...")
raw_data <- breast_epith.md[["RNA"]]@counts



cds = newCellDataSet(cellData = as.matrix(raw_data), phenoData=NULL, 
      featureData=NULL,lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())
print("printing cds made using newCellDataSet function")
print(cds)
cell.labels = seurat.obj.md$seurat_clusters
pData(cds)$Cluster    = cell.labels

print("printing cds after adding cluster to pdata")
print(cds)

pData(cds)$patient_id = seurat.obj.md$patient_id
print("running estimatesizefactors for cds")
cds                   = estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.2)

#disp_table <- dispersionTable(cds)

#unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1, qval)


expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 500))

table(unsup_clustering_genes$gene_id %in% row.names(cds))



clustering_DEG_genes <- differentialGeneTest(cds[expressed_genes,], 
                                             fullModelFormulaStr = '~Cluster')
ordering_genes <- row.names (subset(clustering_DEG_genes, qval < 0.01))



cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds, max_components = 2,
                            method = 'DDRTree')

cds <- orderCells(cds)
cds <- orderCells(cds, root_state = 3)
plot_cell_trajectory(cds, color_by = "Pseudotime")
plot_cell_trajectory(cds, color_by = "State")
plot_cell_trajectory(cds, color_by = "Cluster")


BEAM_res <- BEAM(cds, branch_point = 1, cores = 1)

BEAM_res <- BEAM_res[order(BEAM_res$qval),]


BEAM_res <- BEAM_res[,c( "pval", "qval")]

flag <- !grepl("^RP|^MT", rownames(BEAM_res))

flag <- rownames(BEAM_res) %in% tf$Symbol
surface <- surface[surface$UniProt.Cell.surface == "yes",]
flag <- rownames(BEAM_res) %in% surface$ENTREZ.gene.symbol

setwd("G:/scRNA-seq/source/TF")
tf <- read.delim("TF.txt")
flag_expressed <- rownames(BEAM_res) %in% expressed_genes

table(BEAM_res$qval <1e-100 & rownames(BEAM_res) %in% tf$Symbol)

tf_select <- intersect(row.names(subset(BEAM_res, qval < 1e-32)), tf$Symbol)

plot_genes_branched_heatmap(cds[tf_select,],branch_labels = c("Luminal", "HER2"),
                            branch_point = 1,
                            num_clusters = 5,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)

genes <- c( "GATA3", "ESR1", "ATF3","CREM","PHB", "NME2")

plot_genes_branched_pseudotime(cds[genes,],
                               branch_point = 1,
                               color_by = "State",
                               ncol = 1)

var.genes.total = c()
print('Computing variable genes ... ')
for (j in 1:length(cell.labels)){
  print(sprintf("Choice %s out of %s ... ", as.character(j), as.character(length(cell.labels))))
  choices = meta$new_label == cell.labels[j]
  var.genes = differentialGeneTest(cds[, choices], fullModelFormulaStr = "~sm.ns(Pseudotime)")
  var.genes = cbind(var.genes, data.frame(gene_id = rownames(var.genes)))
  var.genes.ch = var.genes %>% arrange(qval)
  var.genes.ch = as.vector(var.genes.ch$gene_id[1:100])
  var.genes.total = union(var.genes.total, var.genes.ch)
}

print("Computing var genes globally...")
var.genes = differentialGeneTest(cds, fullModelFormulaStr = "~sm.ns(Pseudotime)")
var.genes = cbind(var.genes, data.frame(gene_id = rownames(var.genes)))
var.genes.ch = var.genes %>% arrange(qval)
var.genes.ch = as.vector(var.genes.ch$gene_id[1:100])
var.genes.total = union(var.genes.total, var.genes.ch)
MT_genes = var.genes.total[grep("^MT-", x=var.genes.total, ignore.case=T)]
var.genes.total = setdiff(var.genes.total, MT_genes)

saveRDS(seurat.obj, "breast_epith.RDS")

heatmap <- data.frame(gene = p$ph$tree_row$labels, order = p$ph$tree_row$order)
heatmap <- heatmap[order(heatmap$order),]
heatmap$survival <- plyr::mapvalues(heatmap$gene, er_survival$genename, er_survival$expef)
