
options(stringsAsFactors = F)

setwd("F:/LC/TCGA/BRCA")
clin <- read.delim("BRCA_clinicalMatrix")

setwd("F:/zhang zemin")
expr <- read.delim("TCGA-BRCA/TCGA-BRCA.htseq_fpkm.tsv.gz")
col <- colnames(expr)
col <- gsub('.', '-', col, fixed = T)
col <- gsub('A$', '', col)
colnames(expr) <- col


inter_name <- intersect(col, clin$sampleID)

expr <- expr[,inter_name]
rownames(clin) <- clin$sampleID
clin <- clin[inter_name,]

ensembl <- expr$Ensembl_ID

expr <- expr[!duplicated(gene_name),]
rownames(expr) <- gene_name[!duplicated(gene_name)]

flag_tumor <- clin$sample_type == "Primary Tumor"
clin <- clin[flag_tumor,]
expr <- expr[,flag_tumor]

flag_pam <- clin$PAM50Call_RNAseq != ""
clin_pam <- clin[flag_pam,]
expr_pam <- expr[,flag_pam]
flag <- rowSums(expr_pam > 0.1) > 80
expr_pam <- expr_pam[flag,]

pd <- new("AnnotatedDataFrame", data = clin_pam)

cds = newCellDataSet(cellData = as.matrix(expr_pam), phenoData=NULL, 
                     featureData=NULL,lowerDetectionLimit = 0.5, 
                     expressionFamily = tobit(Lower = 0.1))

rpc_matrix <- relative2abs(cds, method = "num_genes")

cds = newCellDataSet(cellData = as.matrix(rpc_matrix), phenoData=pd, 
                     featureData=NULL,lowerDetectionLimit = 0.5, 
                     expressionFamily = negbinomial.size())


cds                   = estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- setOrderingFilter(cds, unique(top_20$gene))
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')

cds <- orderCells(cds)
cds <- orderCells(cds, root_state = 2)


plot_cell_trajectory(cds, color_by = "PAM50Call_RNAseq")
plot_cell_trajectory(cds, color_by = "ER_Status_nature2012")
plot_cell_trajectory(cds, color_by = "State")
plot_cell_trajectory(cds, color_by = "pathologic_M")


pData(cds)$Cluster    = cell.labels

ddr <- cds@reducedDimS
ddr <- t(ddr)
colnames(ddr) <- c("ddr_1", "ddr_2")
ddr <- as.data.frame(ddr)

gene = "TYMS"
value <- rpc_matrix[gene,]
value <- log(value +0.1)
value[value <0] <- 0


ggplot(ddr, aes(ddr_1, ddr_2, color = value)) + geom_point() 


