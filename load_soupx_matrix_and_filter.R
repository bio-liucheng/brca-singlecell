#load data from soupx adjust matrix

library(Seurat)

data_dir <- "trans_all"
save_dir <- "soup_results"
samples <- list.files(data_dir)

for(i in 1:(length(samples))){
  sample <- samples[i]
  savePath <- file.path(save_dir, sample)
  seurat.obj <- readRDS(file.path(savePath, "SoupX_adjust_mtx.RDS"))
  seurat.obj <- CreateSeuratObject(counts = seurat.obj, 
                                   project = sample, 
                                   min.cells = 3,
                                   min.features = 300)
  seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "^MT-")
  
  
  #filter features
  median.nf <- median(log10(seurat.obj$nFeature_RNA))
  mad.nf <- mad(log10(seurat.obj$nFeature_RNA))
  trd.nf <- median.nf-3*mad.nf
  flag.nf <- log10(seurat.obj$nFeature_RNA) > trd.nf
  
  #filter percet.mt
  median.pmt <- median(seurat.obj$percent.mt)
  mad.pmt <- mad(seurat.obj$percent.mt)
  trd.pmt <- median.pmt+3*mad.pmt
  flag.pmt <- seurat.obj$percent.mt < trd.pmt
  
  qc <- flag.nf & flag.pmt
  seurat.obj$qc <- qc
  
  if(i == 1){
    seurat.obj.all <- seurat.obj
  } else{
    seurat.obj.all <- merge(seurat.obj.all, seurat.obj)
  }
}

setwd("F:/scCancer_analysis")
setwd("F:/raw")

predict_doublet <- plyr::mapvalues(colnames(seurat.obj.all), as.vector(doublet$barcode), as.vector(doublet$predicted_doublets))
