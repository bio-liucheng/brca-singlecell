library(Seurat)
library(SoupX)

bg.spec.genes <- list(
  igGenes = c('IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4','IGHD','IGHE','IGHM', 'IGLC1','IGLC2','IGLC3','IGLC4','IGLC5','IGLC6','IGLC7', 'IGKC', 'IGLL5', 'IGLL1'),
  HLAGenes = c('HLA-DRA', 'HLA-DRB5', 'HLA-DRB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DQB1', 'HLA-DQA2', 'HLA-DQB2', 'HLA-DPA1', 'HLA-DPB1'),
  HBGenes = c("HBB","HBD","HBG1","HBG2", "HBE1","HBZ","HBM","HBA2", "HBA1","HBQ1"),
  SGGenes = c("SCGB2A2", "KRT19")
  
)

process_seurat <- function(seurat.obj){
  seurat.obj <- NormalizeData(seurat.obj)
  seurat.obj <- FindVariableFeatures(seurat.obj)
  seurat.obj <- ScaleData(seurat.obj)
  seurat.obj <- RunPCA(seurat.obj)
  seurat.obj <- FindNeighbors(seurat.obj, reduction = "pca", dims = 1:40)
  seurat.obj <- FindClusters(seurat.obj, resolution = 0.8)
  return(seurat.obj)
}

data_dir <- "trans_all"
save_dir <- "soup_results"
samples <- list.files(data_dir)

for(i in 8:(length(samples)-1)){
  sample <- samples[i]
  dataPath <- file.path(data_dir, sample)
  # A path used to save the results files
  savePath <- file.path(save_dir, sample)
  message("[", Sys.time(), "] -----: load data in scopx")
  sc <- load10X(dataPath)
  seurat.obj <- Read10X(paste0(dataPath, "/filtered_feature_bc_matrix/"))
  seurat.obj <- CreateSeuratObject(counts = seurat.obj, 
                                 project = sample, 
                                 min.cells = 0, 
                                 min.features = 0)
  message("[", Sys.time(), "] -----: cluster finder")
  seurat.obj <- process_seurat(seurat.obj)
  
  cluster <- seurat.obj$seurat_clusters
  names(cluster) <- colnames(seurat.obj)
  sc$metaData$clusters <- cluster
  useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = bg.spec.genes, 
                                        clusters = cluster)

  sc = calculateContaminationFraction(sc, bg.spec.genes, useToEst = useToEst)
  
  message("[", Sys.time(), "] -----: adjustCounts counting:")
  out = adjustCounts(sc, nCores = 1, roundToInt = T)
  cell_names <- paste(sample, colnames(out), sep = "_")
  colnames(out) <- cell_names 
  saveRDS(out, file.path(savePath, "SoupX_adjust_mtx.RDS"))
  saveRDS(sc$toc, file.path(savePath, "SoupX_pre_adjust_mtx.RDS"))
  write.csv(sc$metaData, file.path(savePath, "SoupX_meta.csv"))
}
