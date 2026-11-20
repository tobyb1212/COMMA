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
library(patchwork)

use_condaenv("r-reticulate")

setwd("/COMMA/results/scatac/matrices")

# Load ArchR
proj <- loadArchRProject(path = '/COMMA/results/scatac/ArchR/filtered', showLogo = FALSE)
print("Loaded ArchR Project")

# Load PeakVI .h5ad as single cell experiment
sce_filtered <- schard::h5ad2sce("/COMMA/results/scatac/matrices/Peak_matrix_final_integrated_filtered.h5ad")

# Visualise embeddings in sce object and export UMAP and PeakVI Coordinates
reducedDimNames(sce_filtered)
umap <- reducedDim(sce_filtered, "X_umap")
peakVI <- reducedDim(sce_filtered, "X_PeakVI")

# Export Leiden Clustering Metadata
leiden <- colData(sce_filtered)$leiden_0.8
names(leiden) <- colnames(sce_filtered)

# Find common cells between SCE and ArchR
common_cells <- intersect(rownames(umap), proj$cellNames)

# Subset embeddings and metadata to shared cells
umap_sub <- umap[common_cells, , drop = FALSE]
colnames(umap_sub) <- c("#UMAP1", "#UMAP2")
peakVI_sub <- peakVI[common_cells, , drop = FALSE]

leiden_sub <- leiden[common_cells]

# Prepare Leiden cluster data as a data.frame
leiden_df <- data.frame(leiden_PeakVI = leiden_sub)
rownames(leiden_df) <- common_cells

# Filter ArchR File
idxSample <- BiocGenerics::which(proj$cellNames %in% common_cells)
proj@reducedDims <- SimpleList()
proj_filtered <- subsetArchRProject(ArchRProj = proj, cells = proj$cellNames[idxSample],
  outputDirectory = "/COMMA/results/scatac/ArchR/filtered_leblanc", dropCells = TRUE, force = TRUE)
proj_filtered <-  loadArchRProject(path = '/COMMA/results/scatac/ArchR/filtered_leblanc', showLogo = FALSE)

# Add embeddings into ArchR project
proj_filtered@embeddings$UMAP <- SimpleList(df = umap_sub,  params = list())
rd <- SimpleList(matRD = peakVI_sub, params = list(corCutOff = NA), scaleDims = FALSE, corToDepth = NA,
  cellColData = proj_filtered@cellColData)
proj_filtered@reducedDims[["PeakVI"]] <- rd

# Add Leiden cluster assignment as metadata
proj_filtered <- addCellColData(
  ArchRProj = proj_filtered,
  data = leiden_df$leiden_PeakVI,
  name = "leiden_PeakVI",
  cells = rownames(leiden_df),
  force = TRUE)

# Subset ArchR Porject to Only RNA Corresponding Cells
valid_cells <- proj_filtered$cellNames[!is.na(proj_filtered$multiome_rna_celltype)]
proj_subset <- subsetArchRProject(ArchRProj = proj_filtered, cells = valid_cells,
  outputDirectory = "/COMMA/results/scatac/ArchR/rna_paired", dropCells = TRUE, force = TRUE)

# Plot UMAPs 
p1 <- plotEmbedding(ArchRProj = proj_subset, colorBy = "cellColData", name = c("multiome_rna_celltype"), embedding = "UMAP", size = 0.05)
ggsave(plot = p1, file="subset_rna_type.png")
p2 <- plotEmbedding(ArchRProj = proj_filtered, colorBy = "cellColData", name = c("leiden_PeakVI"), embedding = "UMAP", size = 0.05)
ggsave(plot = p2, file="subset_rna_cluster.png")

# Assign Clusters with 40% of Cells have a High Doublet Enrichment Score
doublet_score <- getCellColData(proj_filtered, "DoubletEnrichment")[, 1]
leiden <- getCellColData(proj_filtered, "leiden_PeakVI")[, 1]
top10_cutoff <- quantile(doublet_score, probs = 0.9, na.rm = TRUE)
cut_by_cluster <- split(doublet_score, leiden)
high_cut_frac <- sapply(cut_by_cluster,
  function(x) mean(x >= top10_cutoff, na.rm = TRUE))
doublet_cluster_ids <- names(high_cut_frac[high_cut_frac >= 0.4])
doublet_flag <- leiden %in% doublet_cluster_ids
names(doublet_flag) <- proj_filtered$cellNames
proj_filtered <- addCellColData(ArchRProj = proj_filtered, data = doublet_flag, name = "doublet_cluster", cells = names(doublet_flag), force = TRUE)

# Rename Map Creation
newClusterNames <- c(
  "Remove",
  "Fibroblasts", 
  "Myeloid Cells", 
  "Endothelial Cells", 
  "Cardiomyocytes", 
  "Cardiomyocytes", 
  "Mural Cells", 
  "Fibroblasts", 
  "Endocardial Endothelial Cells", 
  "Lymphoid Cells", 
  "Adipocytes", 
  "Cardiomyocytes", 
  "Cardiomyocytes", 
  "Glial Cells", 
  "Remove")

oldClusterNames <- as.character(sort(as.numeric(unique(proj_filtered$leiden_PeakVI))))
proj_filtered$broad_cell_type <- mapLabels(proj_filtered$leiden_PeakVI,oldLabels = oldClusterNames,newLabels = newClusterNames)

# Remove Cells where snATAC Cell Type != sc/snRNA Cell Type or Annotate as Remove/Doublet Cluster
keep_cells <- getCellNames(proj_filtered)[
  proj_filtered$broad_cell_type != "Remove" &
  (is.na(proj_filtered$broad_cell_type) |
    is.na(proj_filtered$multiome_rna_celltype) |
    as.character(proj_filtered$broad_cell_type) ==
    as.character(proj_filtered$multiome_rna_celltype)) & !as.logical(proj_filtered$doublet_cluster)]
proj_clustered <- subsetArchRProject(ArchRProj = proj_filtered, cells = keep_cells,
  outputDirectory = "/COMMA/results/scatac/ArchR/clustered", dropCells = TRUE, force = TRUE)

# Recompute UMAP after Cell Removal
proj_clustered <- addUMAP(ArchRProj = proj_clustered, reducedDims = "PeakVI", name = "UMAP", nNeighbors = 30, minDist = 0.05, force = TRUE)
p1 <- plotEmbedding(ArchRProj = proj_clustered, colorBy = "cellColData", name = "Study", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj_clustered, colorBy = "cellColData", name = "broad_cell_type", embedding = "UMAP")
p_combined <- p1 + p2
ggsave(filename = "UMAP_Study_and_BroadCellType_0.02.png", plot = p_combined,width = 12,height = 6,dpi = 300)

# Save ArchR Project
saveArchRProject(proj_clustered, '/COMMA/results/scatac/ArchR/clustered', load = FALSE)

# Export to H5AD
obs <- getCellColData(proj_clustered)
obs <- as.data.frame(obs)
umap <- as.data.frame(getEmbedding(proj_clustered, embedding = "UMAP"))
umap <- umap[rownames(obs), , drop = FALSE]
var = as.data.frame("1")
X <- matrix(0L, nrow = nrow(obs), ncol = length(var))

adata <- AnnData(
  X = X,
  obs = obs,
  var = var)
adata$obsm$umap <- umap

outfile = "/COMMA/results/scatac/whole_clustered.h5ad"
adata$write_h5ad(outfile)