library(ArchR)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(reticulate)
library(anndata)
library(Seurat)

proj <- loadArchRProject(path = '/COMMA/results/scatac/ArchR/clustered', showLogo = FALSE)

# Filter for Vascular Cells
idxSample <- BiocGenerics::which(proj$broad_cell_type %in% c("Endothelial Cells", "Mural Cells"))
proj_vascular <- subsetArchRProject(ArchRProj = proj,cells = proj$cellNames[idxSample],outputDirectory = "/COMMA/results/scatac/ArchR/vascular",dropCells = TRUE,force = TRUE)

# Add Group Coverages
proj_vascular <- addGroupCoverages(ArchRProj = proj_vascular, groupBy = "broad_cell_type", sampleLabels = "Sample", minReplicates = 3, maxReplicates = 3, minCells = 40, maxCells = 2000, force = TRUE)

# Recall Peaks 
proj_vascular <- addReproduciblePeakSet(
  ArchRProj = proj_vascular, 
  groupBy = "broad_cell_type", 
  threads = 10,
  reproducibility = "2",
  pathToMacs2 = 'conda/skylake/envs/MACS2/bin/macs2')

# Transfer PeakSet & Generate Peak Matrix
vascular_peak_set <- getPeakSet(proj_vascular)
proj_vascular <- addPeakSet(ArchRProj = proj_vascular,peakSet = vascular_peak_set, force = TRUE)
proj_vascular <- addPeakMatrix(proj_vascular)

proj_vascular@cellColData$multiome_rna_celltype <- NULL

# Save ArchR Project
saveArchRProject(proj_vascular, '/COMMA/results/scatac/ArchR/vascular', load = FALSE)

arrow_files <- list.files(path = "/COMMA/results/scatac/ArchR/vascular/ArrowFiles/",pattern = "\\.arrow$",full.names = TRUE)
for (arrow in arrow_files) {
# Convert ArchR Peak Matrix to Seperate Components for AnnData
peak_matrix <- getMatrixFromArrow(ArrowFile=arrow, ArchRProj = proj_vascular, useMatrix = "PeakMatrix", binarize=TRUE)
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
outfile <- file.path("/COMMA/results/scatac/matrices/vascular/peaks", paste0(base, "_PeakMatrix.h5ad"))

# Save to .h5ad File
adata$write_h5ad(outfile)
cat("Saved: ", outfile)}
