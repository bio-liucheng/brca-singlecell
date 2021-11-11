files <- list.dirs()

files <- files[3:5]

for(i in 1:length(files)){
  meta <- read.csv(file.path(files[i], "meta.csv"))
  if(i == 1)
    meta_all = meta
  else
    meta_all = rbind(meta_all, meta)
}

meta <- read.csv("meta.csv")

meta$cell_label <- meta$leiden

meta$cell_label <- plyr::mapvalues(meta$X, meta_all$X, meta_all$cell_label)

meta$cell_label[grepl('_', meta$cell_label)] <- "A"

meta$cell_label[meta$leiden == 9] <- "plasma B"

meta$cell_label[meta$leiden == 0] <- "B"

meta$cell_label[meta$leiden == 13] <- "KI67+ cells"

write.csv(meta, "meta.csv", row.names = F)
