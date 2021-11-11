library(scCancer)

data_dir <- "trans_all"
save_dir <- "results"
samples <- list.files(data_dir)

# A path containing the cell ranger processed data
for(i in length(samples)){
  sample <- samples[i]
dataPath <- file.path(data_dir, sample)
# A path used to save the results files
savePath <- file.path(save_dir, sample)
# The sample name
sampleName <- sample
# The author name or a string used to mark the report.
authorName <- "LC@bj"

# Run scStatistics

stat.results <- runScStatistics(
  dataPath = dataPath,
  savePath = savePath,
  sampleName = sampleName,
  authorName = authorName,
  bool.runSoupx = F
)

anno.results <- runScAnnotation(
  dataPath = dataPath,
  statPath = savePath,
  savePath = savePath,
  authorName = authorName,
  sampleName = sampleName,
  genome = "hg38",
  geneSet.method = "GSVA"   # or "GSVA"
)

}



