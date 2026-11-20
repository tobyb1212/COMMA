library(ArchR)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(anndata)
library(Seurat, lib.loc = "R/4.3/skylake/seurat_v5.3.0")
Sys.setenv(RETICULATE_PYTHON = "/conda/skylake/envs/r-reticulate/bin/python")
library(reticulate)
library(optparse)
library(schard)

proj <- loadArchRProject(path = '/COMMA/results/scatac/ArchR/vsmc_clustered', showLogo = FALSE)
print("Loaded ArchR Project")

arrow_files <- list.files(path = "/COMMA/results/scatac/ArchR/vsmc_clustered/ArrowFiles/",pattern = "\\.arrow$",full.names = TRUE)
for (arrow in arrow_files) {
# Convert ArchR Peak Matrix to Seperate Components for AnnData
peak_matrix <- getMatrixFromArrow(ArrowFile=arrow, ArchRProj = proj, useMatrix = "PeakMatrix", binarize=TRUE)
counts <- t(assay(peak_matrix))
cell_metadata <- as.data.frame(colData(peak_matrix))
peaks <- rowRanges(peak_matrix)
peak_metadata <- data.frame(idx=seq_along(peaks),chr=seqnames(peaks), start=start(peaks), end=end(peaks))

# Prepare AnnData for PeakVI
adata <- AnnData(
  X = counts,
  obs = cell_metadata,
  var = peak_metadata
)
base <- gsub("\\.arrow$", "", basename(arrow))
outfile <- file.path("/COMMA/results/scatac/multiomic_integration/vsmc/peak_matrices/", paste0(base, "_PeakMatrix.h5ad"))

# Save to .h5ad File
adata$write_h5ad(outfile)
cat("Saved: ", outfile)}
