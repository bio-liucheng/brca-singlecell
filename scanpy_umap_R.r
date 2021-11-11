options(stringsAsFactors = F)
library(ggplot2)
umap <- read.csv("scanpy_umap.csv")
colnames(umap) <- c('num', 'UMAP_1', 'UMAP_2')
meta <- read.csv("meta.csv")
umap <- umap[,c(2,3)]
rownames(umap) <- colnames(seurat.obj.md)
umap <- as.matrix(umap)
seurat.obj.md[["umap"]] <- CreateDimReducObject(embeddings = umap, 
                                                key = "UMAP_", 
                                                assay = DefaultAssay(seurat.obj.md))


seurat.obj.md <- seurat.obj[,meta$X]
meta <- cbind(meta, umap)
ggplot(meta, aes(UMAP_1, UMAP_2, color = new_label)) +geom_point(size =1, alpha = 0.8) + 
  theme_void() + scale_color_manual(values = cl[1:11]) + 
  theme(legend.title = element_blank(),
        legend.text = element_text(face="bold",size=10))

label_data <- meta %>%
  group_by(cell_label) %>%
  summarise(UMAP1 = mean(umap1),
            UMAP2 = mean(umap2))

  meta$label_color <- ifelse(meta$cell_label == targert, "red", "gray")
  
  
  require(scales) 
  n <- length(unique(meta$cell_label)) # number of colors 
  cols <- hue_pal(h = c(0, 360) + 15, 
                  c = 100, l = 65, 
                  h.start = 0, direction = 1)(n)[order(sample(1:n, n))] # color palette in random order 
  ggplot(data,aes(x,y,colour=category))+stat_smooth(alpha=0, size=2) + 
    scale_color_manual(values = cols) 
  
  
color <- sample()
 ggplot(meta) +
  geom_point(aes(x = umap1, y = umap2, colour = factor(new_label)), alpha = 0.5, size = 0.6)  + 
  geom_point(data = label_data, aes(x = UMAP1, y = UMAP2, colour = cell_label),
             shape = 21, size = 4, stroke = 1, fill = "white") +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_text(data = label_data,
            aes(x = UMAP1, y = UMAP2,label = cell_label)) +
  ggtitle("Cell velocity field")

 
 
plot_meta <- function(meta, target){
 
 meta$label_color <- ifelse(meta$cell_label == target, "red", "gray")
 ggplot(meta) +
   geom_point(aes(x = umap1, y = umap2, colour = label_color), alpha = 0.3, size = 0.6, color =  hue_pal()(4)[2]) +
   theme_minimal() +
   theme(legend.position = "none") +
   ggtitle("Cell velocity field") 
}

