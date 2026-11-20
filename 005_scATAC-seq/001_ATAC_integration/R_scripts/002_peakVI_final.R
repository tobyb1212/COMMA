library(ArchR)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(reticulate)
library(anndata)
library(Seurat)

proj <- loadArchRProject(path = '/COMMA/results/scatac/ArchR/filtered', showLogo = FALSE)
print("Loaded ArchR Project")

# Subset ArchR Object for Only Multiome Cells
idxSample <- BiocGenerics::which(!is.na(proj$multiome_rna_celltype)) 
proj_multiome <- subsetArchRProject(ArchRProj = proj,cells = proj$cellNames[idxSample],outputDirectory = "/COMMA/results/scatac/ArchR/multiome_only",dropCells = TRUE,force = TRUE)

# Add Group Coverages
proj_multiome <- addGroupCoverages(ArchRProj = proj_multiome, groupBy = "multiome_rna_celltype", sampleLabels = "Sample", minReplicates = 3, maxReplicates = 3, minCells = 40, maxCells = 2000, force = TRUE)

proj_multiome <- addReproduciblePeakSet(
  ArchRProj = proj_multiome, 
  groupBy = "multiome_rna_celltype", 
  threads = 10,
  reproducibility = "2",
  pathToMacs2 = 'conda/skylake/envs/MACS2/bin/macs2')

# Transfer PeakSet & Generate Peak Matrix
multiome_peak_set <- getPeakSet(proj_multiome)
proj <- addPeakSet(ArchRProj = proj,peakSet = multiome_peak_set, force = TRUE)
proj = addPeakMatrix(proj)

proj@cellColData$multiome_rna_celltype <- NULL

arrow_files <- list.files(path = "/COMMA/results/scatac/ArchR/filtered/ArrowFiles/",pattern = "\\.arrow$",full.names = TRUE)
for (arrow in arrow_files) {
# Convert ArchR Peak Matrix to Seperate Components for AnnData
peak_matrix <- getMatrixFromArrow(ArrowFile=arrow, ArchRProj = proj, useMatrix = "PeakMatrix", binarize=TRUE)
counts <- t(assay(peak_matrix))
cell_metadata <- as.data.frame(colData(peak_matrix))
peak_metadata <- as.data.frame(rowData(peak_matrix))

# Prepare AnnData for PeakVI
adata <- AnnData(
  X = counts,
  obs = cell_metadata,
  var = peak_metadata
)
base <- gsub("\\.arrow$", "", basename(arrow))
outfile <- file.path("/COMMA/results/scatac/matrices/peaks", paste0(base, "_PeakMatrix.h5ad"))

# Save to .h5ad File
adata$write_h5ad(outfile)
cat("Saved: ", outfile)}