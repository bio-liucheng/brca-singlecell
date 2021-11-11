## subset Seurat object for the cluster
options(stringsAsFactors = F)
meta <- read.csv("meta.csv")
library(Seurat)
seurat.obj.md <- seurat.obj[,meta$X]
patient_list = paste("patient", c(1,2,3,4,5,6,7,8), sep = '_')

tt_list = list()
for(i in 1:length(patient_list)){

seurat_object_edgeR <- subset(seurat.obj.md,  patient_id %in% c("patient_3", "patient_5", "patient_6"))
seurat_object_edgeR <- subset(seurat.obj.md,  monocle_state == 2 & Malign.type == "malignant")
seurat_object_edgeR <- NormalizeData(seurat_object_edgeR)
seurat_object_edgeR <- FindVariableFeatures(seurat_object_edgeR)
seurat_object_edgeR <- ScaleData(seurat_object_edgeR)

seurat_object_edgeR <- RunPCA(seurat_object_edgeR)
seurat_object_edgeR <- FindNeighbors(seurat_object_edgeR, reduction = "pca", dims = 1:20)
seurat_object_edgeR <- FindClusters(seurat_object_edgeR, resolution = 0.8)
seurat_object_edgeR <- RunTSNE(seurat_object_edgeR, tsne.method = "FIt-SNE")
#seurat.obj.md <- seurat.obj[,meta$X]
#Idents(seurat.obj.md) <- seurat.obj.md$ER_statu
meta_sub <- meta[meta$patient_id == "patient_5",]
library(edgeR)
# count matrix 

seurat.obj.md$cell_label <- meta_sub$cell_label
cell_label  <- unique(seurat.obj.md$cell_label)
i = 2
seurat_object_edgeR <- subset(seurat.obj.md,  cell_label == "DC-C3-LAMP3")
seurat_object_edgeR <- subset(seurat.obj.md,  cell_label == "CD4-C6-CXCL13")
seurat_object_edgeR <- subset(seurat.obj.md,  cell_label %in% c("CAF-C1-POSTN","CAF-C3-PLA2G2A"))
seurat_object_edgeR <- subset(seurat.obj.md,  leiden %in% c(2,3))
counts <- as.matrix(seurat_object_edgeR@assays$RNA@counts)
counts <- counts[Matrix::rowSums(counts >= 1) >= ncol(counts)*0.1, ]
# subset the meta data fro filtered gene/cells
metadata <- seurat_object_edgeR@meta.data
metadata$temp <- ifelse(metadata$sample_type == "Tumor", "Tumor", "LN")

metadata$temp <- ifelse(metadata$cell_label == "CAF-C1-POSTN", "t", "m")
#metadata$temp <- ifelse(metadata$cell_label %in% c("Macro-C5-SPP1"), "SPP1", "OTHER")
metadata <- metadata[,c("nCount_RNA", "temp")]
metadata <- metadata[colnames(counts),]

dge <- DGEList(counts, group = metadata$temp)

dge <- calcNormFactors(dge)
design <- model.matrix(~metadata$temp)
dge <- estimateDisp(dge, design = design)
fit <- glmQLFit(dge, design = design)

qlf <- glmQLFTest(fit)
tt <- topTags(qlf, n = Inf)

tt_list[[i]] <- list(tt =tt, patient = patient_list[i])
}

col_annotation <- data.frame(sample_type = seurat_object_edgeR$sample_type)
rownames(col_annotation) <- colnames(seurat_object_edgeR)

gene <- VariableFeatures(seurat_object_edgeR)[1:100]
gene <- gene[!grepl("HLA|IG",gene)]
data <- seurat_object_edgeR[["RNA"]]@data[gene,]
data[data >4] <- 4
pheatmap(data, show_rownames = T, show_colnames = F, annotation_col = col_annotation)
