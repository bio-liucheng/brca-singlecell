library(Seurat)
library(Matrix)
library(infercnv)


seurat.obj <- Read10X("filtered_gene_bc_matrices/hg19/")
seurat.obj <- CreateSeuratObject(counts = seurat.obj, 
                                 project = 'pbmc', 
                                 min.cells = 3, 
                                 min.features = 200)

# save file two files
## Counts file
## Cell label files
saveRDS(as.matrix(seurat.obj[["RNA"]]@counts), "sc.10x.counts.matrix")
write.table(seurat.obj$Malign.type, "seurat.obj.label.txt", sep = "\t", quote = F, col.names = F)

# load file
matrix_counts <- readRDS("sc.10x.counts.matrix")
ref <- read.delim("seurat.obj.label.txt", header = F)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = matrix_counts,
                                    annotations_file="seurat.obj.label.txt",
                                    delim="\t",
                                    gene_order_file="gencode_v19_gene_pos.txt",
                                    ref_group_names="nonMalignant")

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="output_2", 
                             cluster_by_groups=TRUE, 
                             no_prelim_plot = T,
                             denoise=TRUE,
                             HMM=F)

