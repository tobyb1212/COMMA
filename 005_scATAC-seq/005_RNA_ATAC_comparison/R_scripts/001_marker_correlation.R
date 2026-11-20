library(ArchR)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(reticulate)
library(anndata)

proj <- loadArchRProject(path = '/COMMA/results/scatac/ArchR/vascular_clustered', showLogo = FALSE)
print("Loaded ArchR Project")

setwd("/COMMA/results/scatac/differential_analysis/genescores/vascular")

# Convert ArchR genescore Matrix to Seperate Components for AnnData
genescore_matrix <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix", binarize=FALSE)
counts <- t(assay(genescore_matrix))
cell_metadata <- as.data.frame(colData(genescore_matrix))
gene_names <- rowData(genescore_matrix)$name
colnames(counts) <- gene_names

# Prepare AnnData for GeneScore 
adata <- AnnData(
  X = counts,
  obs = cell_metadata,
  var = data.frame(index = gene_names, row.names = gene_names)
)

# Save to .h5ad File
adata$write_h5ad("vascular_genescore_matrix.h5ad")
