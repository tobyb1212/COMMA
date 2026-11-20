library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(stringr)

setwd("/COMMA/results/scatac/arrow_files")

addArchRGenome("hg38")
addArchRThreads(threads = 8) 

ArrowFiles <- createArrowFiles(
  inputFiles = input_files,
  sampleNames = names(input_files),
  minTSS = 3, 
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = TRUE
)

arrow_names <- paste0(names(input_files), ".arrow")

doubletScores <- addDoubletScores(
    input = arrow_names,
    k = 10, 
    knnMethod = "UMAP",
    LSIMethod = 1
)

proj <- ArchRProject(
  ArrowFiles = arrow_names, 
  outputDirectory = "/COMMA/results/scatac/ArchR/initial",
  copyArrows = FALSE, 
  threads = 8
)

saveArchRProject(ArchRProj = proj, outputDirectory = "/COMMA/results/scatac/ArchR/initial", load = FALSE)