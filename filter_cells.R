#add metadata to seurat.obj

library(plyr)
identical(colnames(seurat.obj), as.character(dl$cell.id))
seurat.obj$doublet <- dl$dl_classification

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

qc <- flag.nf & flag.pmt & seurat.obj$doublet == "Singlet"
unique(seurat.obj$doublet)

qc <- flag.nf & flag.pmt & seurat.obj$doublet == "Singlet"
seurat.obj$qc <- qc
seurat.obj$sample_type <- substr(seurat.obj$orig.ident, 1, 2)
seurat.obj$sample_type <- ifelse(seurat.obj$sample_type == "CA", "Tumor", "Lymph Node")
seurat.obj$patient_id <- substr(seurat.obj$orig.ident, 3, 3)
seurat.obj$patient_id <-paste0("patient_", seurat.obj$patient_id)

DimPlot(seurat.obj, reduction = "umap", group.by = "orig.ident")
FeaturePlot(seurat.obj, features =c("nFeature_RNA", "percent.mt"))

unique(seurat.obj.md$Malign.type)
