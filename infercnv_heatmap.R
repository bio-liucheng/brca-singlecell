library(ComplexHeatmap)
library(circlize)

infer_cnv <- readRDS("run.final.infercnv_obj")
gene_pos <- read.delim("gencode_v19_gene_pos.txt", header = F)
gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr),]

new_cluster <- unique(gene_pos$V2)

top_color <- HeatmapAnnotation(
  cluster = anno_block( # ÉèÖÃÌî³äÉ«
                       labels = gsub("chr", "", new_cluster), 
                       gp = gpar(col = "white"),
                       labels_gp = gpar(cex = 1, col = "black"),height = unit(2,"mm"))) 


text_list = as.list(new_cluster)

p0 <- Heatmap(meta_sort$Malign.type[!flag], name = "Malign.type", col = c("grey"), width = unit(5, "mm")) 
p1 <- Heatmap(meta_sort$patient_id[!flag], name = "patient", col = color_map_s$cell_color, width = unit(5, "mm")) 
p2 <- Heatmap(meta_sort$sample_type[!flag], name = "sample type", col = rev(color_map_s2$cell_color), width = unit(5, "mm"))
p4 <- Heatmap(meta_sort$ER_statu[!flag], name = "ER", width = unit(5, "mm"), col = r)

p0 <- Heatmap(meta_sort$Malign.type[flag], name = "Malign.type", col = c("light blue"), width = unit(5, "mm")) 
p1 <- Heatmap(meta_sort$patient_id[flag], name = "patient", col = color_map_s$cell_color, width = unit(5, "mm")) 
p2 <- Heatmap(meta_sort$sample_type[flag], name = "sample Type", col = rev(color_map_s2$cell_color), width = unit(5, "mm"))
p4 <- Heatmap(meta_sort$ER_statu[flag], name = "ER", width = unit(5, "mm"), col = r)


flag <- meta$Malign.type == "malignant"

s <- which(meta$Malign.type == "malignant")
flag <- meta_sort$Malign.type == "malignant"
col_fun = colorRamp2(c(-0.4, 0, 0.4), c("navy", "white", "firebrick3"))

col_fun_mag = colorRamp2(c(1, 4), c("white", "red"))

pdf("chr11.pdf", width = 18, height = 8)
p3 <- Heatmap(as.matrix(t(log2(expr_s[,flag]))),
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = F,
        heatmap_legend_param = list(
          title = "Expression",
          title_position = "leftcenter-rot"
        ), 
        col = col_fun,
        column_split = factor(gene_pos$V2, new_cluster),
        bottom_annotation = top_color,
        width = unit(25, "cm"),
        height = unit(5, "cm"),
        column_title = NULL)
p0 + p1+p2+p4+p3
p3
dev.off()



pdf("chr_patient_7_cluster.pdf", width = 10, height = 5)
p3 <- Heatmap(as.matrix(t(log2(expr_s))),
              cluster_rows = T,
              cluster_columns = F,
              show_column_names = F,
              show_row_names = F,
              heatmap_legend_param = list(
                title = "Expression",
                title_position = "leftcenter-rot"
              ), 
              column_split = factor(gene_pos$V2, new_cluster),
              bottom_annotation = top_color,
              col = col_fun,
              column_title = NULL)
p2+p3
dev.off()






library(dplyr)
meta <-read.csv("meta.csv")
meta_sort  <- meta %>% arrange(ER_statu, patient_id, desc(sample_type))

#meta_sort  <- meta %>% arrange(ER_statu, patient_id, desc(sample_type))
gene_pos <- read.delim("gencode_v19_gene_pos.txt", header = F)
gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr), ]
setwd("G:/scRNA-seq/scCancer_analysis/soup_seurat/epi_cell/infercnv/output_2")
infercnv <- readRDS("run.final.infercnv_obj")
expr <- infercnv@expr.data

flag <- meta$patient_id == "patient_1"
#flag <- meta_sort$patient_id == "patient_8"& meta_sort$Malign.type == "malignant"
expr_s <- expr[,meta$X]


#flag_gene_pos <- grepl("chr19", gene_pos$V2)

expr_s <- expr_s[,flag]
s <- dist(t(log2(expr_s)))
s <- hclust(s)
#order <- s$order



st<- meta_sort[flag,]$sample_type
ms <- meta_sort[flag,]$Malign.score
#p0 <- Heatmap(meta_sort$Malign.type[flag], name = "Malign.type", width = unit(5, "mm")) 
#p1 <- Heatmap(meta_sort$patient_id[flag], name = "patient", width = unit(5, "mm")) 

p2 <- Heatmap(st[order], name = "Sample Type", col = rev(color_map_s2$cell_color), width = unit(5, "mm"))

pm <- Heatmap(ms[order], name = "Malign.score", cluster_rows = F, width = unit(5, "mm"))
pdf("chr_patient_8_cluster2.pdf", width = 8, height = 5)
p3 <- Heatmap(as.matrix(t(log2(expr_s))),
              cluster_rows = T,
              clustering_distance_rows = function(m) dist(m),
              cluster_columns = F,
              show_column_names = F,
              show_row_names = F,
              right_annotation = rowAnnotation(sample_type = meta_sort[flag,]$sample_type ,
                                               Malig_score = meta_sort[flag,]$Malign.score,
                                               col =  list(sample_type =c("Tumor" = "#F0A19D", 
                                                                          "Lymph Node" = "#80C9AE"),
                                                           Malig_score = col_fun_mag
                                                           )
                                               ),
              
              heatmap_legend_param = list(
                title = "Expression",
                title_position = "leftcenter-rot"
              ), 
              row_km = 5,
              column_split = factor(gene_pos$V2, new_cluster),
              col = col_fun,
              column_title = NULL)
p2+pm+p3
p3+p2
p3
dev.off()

p2 <- Heatmap(st, name = "Sample Type", col = rev(color_map_s2$cell_color), width = unit(5, "mm"))
p3 <- Heatmap(as.matrix(t(log2(expr_s[,order]))),
              cluster_rows = T,
              clustering_distance_rows = function(m) dist(m),
              cluster_columns = F,
              show_column_names = F,
              show_row_names = F,
              heatmap_legend_param = list(
                title = "Expression",
                title_position = "leftcenter-rot"
              ), 
              column_split = factor(gene_pos$V2, new_cluster),
              bottom_annotation = top_color,
              col = col_fun,
              column_title = NULL)

