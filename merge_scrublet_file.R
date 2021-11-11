list.dirs()
dir <- list.files()
dir <- dir[grepl("matrix", dir)]

sample <- substr(dir, 1,3)
for(i in 1:length(dir)){
  file <- read.csv(paste0(dir[i], "/filtered_feature_bc_matrix/doublet.txt"))
  barcode <- file$barcode
  barcode <- substr(barcode, 1, 16)
  file$barcode <- paste(sample[i], barcode, sep = "_")
  if(i ==1)
    file_all = file
  else
    file_all = rbind(file_all, file)
}
