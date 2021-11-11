
meta <- read.csv("meta_tcr.csv")
flag_t <- grepl("CD4|CD8", meta$cell_label)
meta <- meta[flag_t,]

flag_tcr <- meta$has_tcr == "True"

meta <- meta[flag_tcr,]


flag_t <- !grepl("_2|_4", meta$patient_id)
meta <- meta[flag_t,]

library(dplyr)

t <- meta %>% group_by(patient_id, sample_type, clonotype) %>% summarise(freq = n())
t$clonotype_id <- paste(t$patient_id, t$clonotype, sep = '_')

library(reshape2)

t_dcast <- dcast(t, clonotype_id~sample_type, value.var = "freq")


meta$clonotype_id <- paste(meta$patient_id, meta$clonotype, sep = '_')

meta$clonotype_tumor_count <- plyr::mapvalues(meta$clonotype_id, 
                                              as.vector(t_dcast$clonotype_id),
                                              as.vector(t_dcast$Tumor)
                                              )
meta$clonotype_ln_count <- plyr::mapvalues(meta$clonotype_id, 
                                              as.vector(t_dcast$clonotype_id),
                                              as.vector(t_dcast$`Lymph Node`)
)
meta$clonotype_tumor_count <- as.numeric(meta$clonotype_tumor_count)
meta$clonotype_ln_count <- as.numeric(meta$clonotype_ln_count)


x = meta$clonotype_tumor_count
y = meta$clonotype_ln_count

z = ifelse(x ==1 & is.na(y), "Tumor_singleton", 
       ifelse(x >1 & is.na(y), "Tumor_mutiplets", 
            
          ifelse(is.na(x) & y == 1, "LN_singleton", 
              ifelse(is.na(x) & y >1, "LN_mutiplets", "Dual_expandr")
              )))

meta$clonotype_type <- z

Idents(seurat.obj.md) <- factor(group, levels = sort(unique(group)))
Idents(seurat.obj.md) <- z
meta$clonotype_type <- z
group <- paste(meta$clonotype_type, meta$patient_id, sep = "_")


aver <- AverageExpression(seurat.obj.md, features = unique(top_20$gene), slot = "scale.data")

pheatmap(aver$RNA, cluster_rows = F, cluster_cols = F)
