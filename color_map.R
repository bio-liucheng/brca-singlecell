
setwd("G:/scRNA-seq/scanpy/ALL2/configure")

color_map <- read.csv("cell_label_color_map3.csv")

color_map_s <- color_map[color_map$cell_label %in% unique(meta$sample_type), ]
color_map_s <- color_map[color_map$cell_label %in% unique(meta$cell_label), ]
color_map_s <- color_map[color_map$cell_label %in% unique(meta$patient_id), ]
color_map_s <- color_map[color_map$cell_label %in% unique(meta$main_label), ]

color_map_s2 <- color_map[color_map$cell_label %in% unique(meta$sample_type),]

color_map_s2 <- color_map[color_map$cell_label %in% unique(meta_sort$sample_type), ]
color_map_s2 <- color_map[color_map$cell_label %in% unique(meta$patient_id), ]

color_map_s <- color_map[color_map$cell_label %in% unique(seurat.obj$cell_label), ]

seurat.obj$cell_label <- factor(seurat.obj$cell_label, levels = color_map_s$cell_label)



#randam color generate

library(randomcoloR)
library(scales)
r <- randomcoloR::randomColor(13, luminosity = "light")
show_col(r)

