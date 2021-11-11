library(Seurat)
library(DoubletFinder)
samplesA <- paste("CA", 7:8, sep = "")
samplesB <-  paste("LN", 7:8, sep = "")


dl_all_metadata <- data.frame(cell.id=character() , dl_classification = character())

for(i in 1:2){
sample <- samplesA[i]
seurat.obj <- Read10X(paste0(sample, "_trans/filtered_feature_bc_matrix/"))
cell_names <- paste(sample, colnames(seurat.obj), sep = "_")
colnames(seurat.obj) <- cell_names 
seurat.obj <- CreateSeuratObject(counts = out, 
                           project = sample, 
                           min.cells = 0, 
                           min.features = 0)


if(i == 1){
  seurat.obj.all <- seurat.obj
} else{
  seurat.obj.all <- merge(seurat.obj.all, seurat.obj)
}


seurat.obj = NormalizeData(object = seurat.obj)
seurat.obj = ScaleData(object = seurat.obj)
seurat.obj = FindVariableFeatures(object = seurat.obj)
seurat.obj <- RunPCA(seurat.obj)
seurat.obj <- FindNeighbors(seurat.obj)
seurat.obj <- FindClusters(seurat.obj)
annotation <- Idents(seurat.obj)

## pK Identification (no ground-truth)
sweep.res.list <- paramSweep_v3(seurat.obj, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

#optimize pk
max_num <- order(bcmvn$BCmetric, decreasing = T)[1]
optimise_pk <- as.numeric(as.character(bcmvn$pK))[max_num]


## Homotypic Doublet Proportion Estimate

nExp_poi <- round(0.04*length(seurat.obj$orig.ident))
homotypic.prop <- modelHomotypic(annotation)
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

reuse.pANN = paste("pANN", 0.25, optimise_pk, nExp_poi, sep = '_', collapse = "")
## Run DoubletFinder with varying classification stringencies
seurat.obj <- doubletFinder_v3(seurat.obj, PCs = 1:10, 
                               pN = 0.25, pK = optimise_pk,
                               nExp = nExp_poi, reuse.pANN = FALSE,
                               sct = FALSE)

seurat.obj <- doubletFinder_v3(seurat.obj, PCs = 1:10, 
                               pN = 0.25, pK = optimise_pk,
                               nExp = nExp_poi.adj, reuse.pANN = reuse.pANN, 
                               sct = FALSE)

dl_metadata <- data.frame(cell.id = rownames(seurat.obj@meta.data),
                          dl_classification = seurat.obj@meta.data[,8])
dl_all_metadata <- rbind(dl_all_metadata, dl_metadata)
}
saveRDS(dl_all_metadata, "CA_dl_meta_sup.RDS")
saveRDS(seurat.obj.all, "CA_sup_seurat.obj.RDS")
rm(dl_all_metadata)
rm(seurat.obj.all)
gc()


dl_all_metadata <- data.frame(cell.id=character() , dl_classification = character())

for(i in 2){
  sample <- samplesB[i]
  seurat.obj <- Read10X(paste0(sample, "_trans/filtered_feature_bc_matrix/"))
  cell_names <- paste(sample, colnames(seurat.obj), sep = "_")
  colnames(seurat.obj) <- cell_names 
  seurat.obj <- CreateSeuratObject(counts = seurat.obj, 
                                   project = sample, 
                                   min.cells = 3, 
                                   min.features = 200)
  
  if(i == 1){
    seurat.obj.all <- seurat.obj
  } else{
    seurat.obj.all <- merge(seurat.obj.all, seurat.obj)
  }
  
  
  seurat.obj = NormalizeData(object = seurat.obj)
  seurat.obj = ScaleData(object = seurat.obj)
  seurat.obj = FindVariableFeatures(object = seurat.obj)
  seurat.obj <- RunPCA(seurat.obj)
  seurat.obj <- FindNeighbors(seurat.obj)
  seurat.obj <- FindClusters(seurat.obj)
  annotation <- Idents(seurat.obj)
  
  ## pK Identification (no ground-truth)
  sweep.res.list <- paramSweep_v3(seurat.obj, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  #optimize pk
  max_num <- order(bcmvn$BCmetric, decreasing = T)[1]
  optimise_pk <- as.numeric(as.character(bcmvn$pK))[max_num]
  
  
  ## Homotypic Doublet Proportion Estimate
  nExp_poi <- round(0.04*length(seurat.obj$orig.ident))
  homotypic.prop <- modelHomotypic(annotation)
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  reuse.pANN = paste("pANN", 0.25, optimise_pk, nExp_poi, sep = '_', collapse = "")
  
  ## Run DoubletFinder with varying classification stringencies
  seurat.obj <- doubletFinder_v3(seurat.obj, PCs = 1:10, 
                                 pN = 0.25, pK = optimise_pk,
                                 nExp = nExp_poi, reuse.pANN = FALSE, 
                                 sct = FALSE)
  seurat.obj <- doubletFinder_v3(seurat.obj, PCs = 1:10, 
                                 pN = 0.25, pK = optimise_pk,
                                 nExp = nExp_poi.adj, reuse.pANN = reuse.pANN, 
                                 sct = FALSE)
  dl_metadata <- data.frame(cell.id = rownames(seurat.obj@meta.data),
                            dl_classification = seurat.obj@meta.data[,8])
  dl_all_metadata <- rbind(dl_all_metadata, dl_metadata)
}
saveRDS(dl_all_metadata, "LN_dl_meta_sup.RDS")
saveRDS(seurat.obj.all, "LN_sup_seurat.obj.RDS")
rm(dl_all_metadata)
rm(seurat.obj.all)
gc()

seurat.obj.ca <- readRDS("CA_seurat.obj.RDS")
seurat.obj.ln <- readRDS("LN_seurat.obj.RDS")
seurat.obj <- merge(seurat.obj.ca, seurat.obj.ln)
saveRDS(seurat.obj, "seurat.obj.all.0.4.RDS")

dl_ca <- readRDS("CA_dl_meta.RDS")
dl_ln <- readRDS("LN_dl_meta.RDS")
df <- rbind(dl_ca, dl_ln)
seurat.obj$doublet <- df$dl_classification

seurat.obj$sample_type <- substr(seurat.obj$orig.ident, 1, 2)