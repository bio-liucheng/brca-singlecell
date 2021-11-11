#add meta data

setwd("F:/patient_info")
add_meta <- read.csv("patiend_add_meta.csv")

meta <- read.csv("meta.csv")

add_col <- colnames(add_meta)[-1]
for(i in 1:length(add_col)){
  meta[,add_col[i]] <- plyr::mapvalues(meta$patient_id, 
                                       as.vector(add_meta$patient_id),
                                       as.vector(add_meta[,add_col[i]]))
}


oi <- strsplit(seurat.obj$orig.ident, "_")
oi <- unlist(lapply(oi, "[", 1))

seurat.obj$patient_id <- paste0("patient_",substr(oi, 3,3))


meta <- meta[!grepl("CA2|CA4", meta$orig.ident),]
