setwd("D:/ZLei/ZONE/Result_X101SC21062649-Z01-J001")
options(stringsAsFactors = F)
library(dplyr)
library(Seurat)
library(patchwork)
library(SingleR)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(celldex)
library(monocle)
library(Rtsne)
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(stringr)
library(limma)
library(VennDiagram)
library(rvest)
library(cowplot)
library(hdf5r)
library(clustree)

#########????????????
i <- 1
while (i<2) {
  if(i==1){
    names <- "A_1938345_11"
    Q <- 'A'
  }
  if(i==2){
    names <- "B_1938529_9"
    Q <- 'B'
  }
  if(i==3){
    names <- "C_2000752_23"
    Q <- 'C'
  }
  if(i==4){
    names <- "D_2000910_33"
    Q <- 'D'
  }
  {
    expr <- paste("1.Count/","/filtered_feature_bc_matrix.h5", sep = names)
    expr.mydata <- Seurat::Read10X_h5(filename =  expr )
    mydata <- Seurat::CreateSeuratObject(counts = expr.mydata, project = Q, assay = 'Spatial')
    mydata$slice <- 1
    mydata$region <- Q #
    imgpath <- paste("1.Count/","/spatial", sep = names)
    img <- Seurat::Read10X_Image(image.dir = imgpath)
    Seurat::DefaultAssay(object = img) <- 'Spatial'
    img <- img[colnames(x = mydata)]
    mydata[['image']] <- img
    mydata <- SCTransform(mydata, assay = "Spatial", verbose = FALSE)
    VFea <- VariableFeatures(mydata)
  }
  if(i==1){ALLdata <- mydata}
  if(i==1){fea <- VFea }
  #if(i!=1){ALLdata <- merge(ALLdata,mydata)}
  if(i!=1){fea <- c(fea, VFea)}
  i<- i+1
}
ALLdata <- mydata
setwd("../")
dir.create("C")
setwd("C")
#########spot??????????
plot1 <- VlnPlot(ALLdata, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(ALLdata, features = "nCount_Spatial") + theme(legend.position = "right")
pdf("Count.pdf", 10,10)
plot_grid(plot1, plot2)
dev.off()

plot1 <- VlnPlot(ALLdata, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(ALLdata, features = "nFeature_Spatial") + theme(legend.position = "right")
pdf("nFeature.pdf", 10,10)
plot_grid(plot1, plot2)
dev.off()

ALLdata[["percent.mt"]] <- PercentageFeatureSet(ALLdata, pattern = "^MT-")
plot1 <- VlnPlot(ALLdata, features = "percent.mt", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(ALLdata, features = "percent.mt") + theme(legend.position = "right")
plot_grid(plot1, plot2)

gene<- data.frame(ALLdata@assays[["Spatial"]]@data@Dimnames[[1]]) 

ALLdata2 <- subset(ALLdata, subset = nFeature_Spatial > 2000 & nFeature_Spatial <10000 & nCount_Spatial > 1000 & nCount_Spatial < 60000)
plot1 <- VlnPlot(ALLdata2, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(ALLdata2, features = "nFeature_Spatial") + theme(legend.position = "right")
plot_grid(plot1, plot2)

ALLdata2
ALLdata


SpatialFeaturePlot(ALLdata, features = c("PLA2G2A", "CD8A"), alpha =c(0.5,1))
DefaultAssay(ALLdata) <- "SCT"

ALLdata <- FindVariableFeatures(ALLdata, selection.method = "vst", nfeature=2000)
pdf("variFeature.pdf",5,5)
VariableFeaturePlot(ALLdata)
dev.off()
ALLdata <- RunPCA(ALLdata, assay = "SCT", verbose = FALSE)
ALLdata <- FindNeighbors(ALLdata, reduction = "pca", dims = 1:30)
pbmc3k.final <- FindClusters(
  object = ALLdata,
  resolution = c(seq(0,1.6,.2))
)
pdf("cluster.pdf",5,5)
clustree(pbmc3k.final@meta.data, prefix = "SCT_snn_res.")
dev.off()
rm(pbmc3k.final)

ElbowPlot(ALLdata)

ALLdata <- FindClusters(ALLdata, verbose = FALSE, resolution = 0.6)
ALLdata <- RunUMAP(ALLdata, reduction = "pca", dims = 1:30)
ALLdata <- RunTSNE(ALLdata, reduction = "pca" , dims = 1:30)

pdf("dimplot.pdf",5,5)
DimPlot(ALLdata, reduction = "tsne", label = TRUE)
dev.off()
pdf("Spatialdimplot.pdf",5,5)
SpatialDimPlot(ALLdata, label = TRUE, label.size = 3, pt.size.factor = 1)
dev.off()
p1 <- DimPlot(ALLdata, reduction = "tsne", label = TRUE)
p2 <- SpatialDimPlot(ALLdata, label = TRUE, label.size = 3, pt.size.factor = 1)
plot_grid(p1, p2)

#SpatialDimPlot(ALLdata, cells.highlight = CellsByIdentities(object = ALLdata, idents = c(1:9)), facet.highlight = TRUE, ncol = 3)
DefaultAssay(ALLdata) <- "Spatial"
exp.mat <- read.table(file = "../cell_cycle_vignette_files/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE, 
                      as.is = TRUE, row.names = 1)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ALLdata <- CellCycleScoring(ALLdata, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

p1 <- DimPlot(ALLdata, group.by = "Phase", reduction = "tsne", label = TRUE)
p2 <- SpatialDimPlot(ALLdata,group.by = "Phase", label.size = 3)
pdf("cycle.pdf",10,10)
plot_grid(p1, p2)
dev.off()
ALLdata@active.ident <- factor(ALLdata@meta.data$seurat_clusters)

ALLdata.markers <- FindAllMarkers(ALLdata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ALLdata.markers$sig <- ALLdata.markers$pct.1 - ALLdata.markers$pct.2

ALLdata.marker <- ALLdata.markers %>% group_by(cluster) %>% top_n(n = 4, wt = sig)
ALLdata.markerss <- ALLdata.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)

hpca.se <- HumanPrimaryCellAtlasData()
Blue.se <- BlueprintEncodeData() 

library(scater)
ALLdatas <- as.SingleCellExperiment(ALLdata)
clusters <-ALLdatas[['seurat_clusters']]

#?????????????????????SingleR???????????????????????????????????????main???fine??????fine????????????????????????????????????????????????????????????main????????????????????????

anno.cluster.main<- SingleR(test = ALLdatas, ref = hpca.se, labels = hpca.se$label.main, method= "cluster", clusters = clusters)
celltype <-data.frame(ClusterID=anno.cluster.main@rownames,celltype=anno.cluster.main@listData[["first.labels"]],stringsAsFactors = F)
ALLdata[['celltype']]<- celltype$celltype[match(Idents(ALLdata), celltype$ClusterID)]

pdf("celltypespatial.pdf",5,5)
SpatialDimPlot(ALLdata, group.by = "celltype", alpha = 0.5)
dev.off()
pred.ALLdata <- SingleR(ALLdatas, ref = hpca.se, assay.type.test=1, labels = hpca.se$label.main)
ALLdata@meta.data$labels<- pred.ALLdata@listData[["first.labels"]]

pred.ALLdatas <- SingleR(ALLdatas, ref = Blue.se, assay.type.test=1, labels = Blue.se$label.main)

ALLdata@meta.data$BLU_labels<- pred.ALLdatas@listData[["first.labels"]]





p1 <- DimPlot(ALLdata, group.by = "celltype", reduction = "tsne", label = TRUE)
p2 <- DimPlot(ALLdata, group.by = "seurat_clusters", reduction = "tsne", label = TRUE)
pdf("dimclustercelltype.pdf",10,5)
plot_grid(p1, p2)
dev.off()

stem <- c("CD24", "CD34", "CD44", "FLOT2", "NANOG", "NES", "POU5F1", "SOX2","ALDH1A1")
stems <- c("CD24", "CD34", "CD44", "FLOT2", "NES")
VlnPlot(ALLdata, features = stems, group.by = "seurat_clusters")
VlnPlot(ALLdata, features = stems, group.by = "celltype")
SpatialFeaturePlot(ALLdata, features = stems,  alpha = c(0.1, 1))
SpatialDimPlot(ALLdata,group.by = "seurat_clusters", label.size = 3)


p1 <- DimPlot(ALLdata, group.by = "celltype", reduction = "tsne", label = TRUE)
p2 <- SpatialDimPlot(ALLdata,group.by = "celltype", alpha = 0.5)
pdf("celltypedimspi.pdf",10,5)
plot_grid(p1, p2)
dev.off()
table(ALLdata@meta.data$labels)
x <- factor(ALLdata@meta.data$labels, levels = unique(ALLdata@meta.data$labels))


tissue_stem <- subset(ALLdata, subset = labels == "Tissue_stem_cells")

p1 <- FeaturePlot(ALLdata, features = stem, reduction = "umap")
p2 <- FeaturePlot(tissue_stem, features = stem, reduction = "tsne")
plot_grid(p1, p2)

DefaultAssay(ALLdata) <- "Spatial"

setwd("../")
saveRDS(ALLdata, "C.rds")
saveRDS(ALLdata.markers, "C_DEGingroup.rds")


Alldata <- readRDS("C.rds")
ALLdata.markers <- readRDS("C_DEGingroup.rds")


HCM<- read.delim2("../Human_cell_markers.txt")
Alldata <- ALLdata
head(Alldata@meta.data)
rm(ALLdata)
HCMX<- HCM[str_detect(HCM$geneSymbol,ALLdata.marker[4,]$gene),]
HCMy <- HCMX[!is.na(HCMX$speciesType),]
genex <- c("CD3D","CD79A","CD4","CD8A","NKG7","GNLY","ACTA2","DCN","PLVAP")
FeaturePlot(Alldata, features = unique(ALLdata.marker$gene), reduction = "tsne")
FeaturePlot(Alldata, features = "CD79A", reduction = "umap")
DotPlot(Alldata,features= unique(ALLdata.marker$gene), cluster.idents=T)
DimPlot(Alldata, group.by ="labels" , split.by = "seurat_clusters", reduction = "umap" )+ facet_wrap(seurat_clusters~.)

p <- SpatialFeaturePlot(Alldata, features =stems, alpha = c(0.1,1))
stems <- c("CD24", "CD34", "CD44", "FLOT2", "NES")
p
setwd("../")
dir.create("results")
setwd("./results")
pdf("SP_CSCmarker.pdf",10,5)
print(p)
dev.off()

stemgene <-  HCM[HCM$cellName%in%c("Cancer stem cell")&HCM$cancerType=="Breast Cancer",]$geneSymbol
stemgene <- strsplit(stemgene, split = ",")
i<-1
while (i<length(stemgene)+1) {
  ifelse(i==1, stemx <- stemgene[[i]],stemx <- c(stemx, stemgene[[i]])) 
  i<-i+1
}
stemx <- unique(stemx)
genex  <- ALLdata.markers[ALLdata.markers$gene%in%stemx,]
geney <- unique(genex$gene)
pdf("stemgenespatial.pdf",9,9)
SpatialFeaturePlot(Alldata, features =geney, alpha = c(0.1,1))
dev.off()

p1 <- DimPlot(Alldata, group.by = "labels", reduction = "tsne", label = TRUE)
p2 <- SpatialDimPlot(Alldata,group.by = "labels", label.size = 1)
plot_grid(p1, p2)

SpatialFeaturePlot(Alldata, features =c("CD151","ADIRF"), alpha = c(0.1,1))
genex <- read.csv("stemgene.csv")
geney <- intersect(genex$Symbol,ALLdata.markers$gene)
genemarker <- ALLdata.markers[ALLdata.markers$gene%in%geney,]%>% top_n(n=18,wt=sig)
pdf("genemarkerstem.pdf",15,10)
SpatialFeaturePlot(Alldata, features =genemarker$gene, alpha = c(0.1,1))
dev.off()

##################gsva


library(Seurat)
library(GSVA)
library(pheatmap)
library(patchwork)


expr <- as.data.frame(Alldata@assays[["Spatial"]]@counts) #????????????
meta <- Alldata@meta.data[,c(9,13)]#??????
m_df = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") #??????????????????
m_df <- m_df[str_detect(m_df$gs_name,"STEM_CELL")&!str_detect(m_df$gs_name,"HEMATOPOIETIC"),]
msigdbr_list = split(x = m_df$gene_symbol, f = m_df$gs_name)
expr=as.matrix(expr) 
kegg <- gsva(expr, msigdbr_list, kcdf="Poisson",method = "gsva",parallel.sz=40) #gsva
saveRDS(kegg, "BP.rds")
p<- pheatmap(kegg, show_rownames=T, show_colnames=F, annotation_col=meta, width=15, height=12,scale = "row",color =colorRampPalette(c("blue", "white","red"))(100))#????????????
pdf("kegggsva.pdf",20,10)
print(p)
dev.off()


es <- data.frame(t(kegg),stringsAsFactors=F)  #?????????????????????????????????????????????????????????umap?????????????????????????????????????????????????????????
scRNA <- AddMetaData(Alldata, es)
saveRDS(scRNA,"GOBPscRNA.rds")
p1 <- FeaturePlot(scRNA, features = "GOBP_SOMATIC_STEM_CELL_POPULATION_MAINTENANCE", reduction = 'tsne')
p2 <- SpatialFeaturePlot(scRNA, features ="GOBP_SOMATIC_STEM_CELL_POPULATION_MAINTENANCE", alpha = c(0.1,1))
pdf("BPsomaticstemmGSVA.pdf", 14,10)
plot_grid(p1, p2)
dev.off()
pdf("immunemarker.pdf", 14,10)
SpatialFeaturePlot(scRNA, features =c("PTPRC","ITGAM","CD247","CD8A","PDCD1","CD274"), alpha = c(0.1,1))
dev.off()


##???????????????????????????????????????
meta <- meta %>%arrange(meta$seurat_clusters)
data <- kegg[,rownames(meta)]
group <- factor(meta[,"seurat_clusters"],ordered = F)
data1 <-NULL
for(i in 0:(length(unique(group))-1)){
  ind <-which(group==i)
  dat <- apply(data[,ind], 1, mean)
  data1 <-cbind(data1,dat)
}
colnames(data1) <-c("C0","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11")
result<- t(scale(t(data1)))
result[result>2]=2
result[result<-2]=-2
library(pheatmap)
p <- pheatmap(result[1:20,],
              cluster_rows = T,
              cluster_cols = F,
              show_rownames = T,
              show_colnames = T,
              color =colorRampPalette(c("blue", "white","red"))(100),
              cellwidth = 10, cellheight = 15,
              fontsize = 10)
pdf(("gsva_cluster.pdf"),width = 10,height = 7)
print(p)
dev.off()



