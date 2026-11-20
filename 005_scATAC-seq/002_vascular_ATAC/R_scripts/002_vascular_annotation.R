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

setwd("/COMMA/results/scatac/matrices/vascular")

addArchRGenome("hg38")
addArchRThreads(threads = 10) 

proj <- loadArchRProject(path = '/COMMA/results/scatac/ArchR/vascular', showLogo = FALSE)

# Add Embeddings from PeakVI
sce <- schard::h5ad2sce("/COMMA/results/scatac/matrices/vascular/filtered_Peak_matrix_integrated.h5ad")
saveRDS(sce, "/COMMA/results/scatac/matrices/vascular/filtered_Peak_matrix_integrated.rds")
umap <- reducedDim(sce, "X_umap")
peakVI <- reducedDim(sce, "X_PeakVI")
leiden <- colData(sce)$leiden_0.8
names(leiden) <- colnames(sce)
common_cells <- intersect(rownames(umap), proj$cellNames)

# Filter Selewa 2023 Dataset
idxSample <- BiocGenerics::which(proj$cellNames %in% common_cells)
proj@reducedDims <- SimpleList()
proj_filtered <- subsetArchRProject(ArchRProj = proj, cells = proj$cellNames[idxSample],
  outputDirectory = "/COMMA/results/scatac/ArchR/vascular_filtered_selewa", dropCells = TRUE, force = TRUE)
proj <-  loadArchRProject(path = '/COMMA/results/scatac/ArchR/vascular_filtered_selewa', showLogo = FALSE)

# Add UMAP Embeddings
umap_sub <- umap[common_cells, , drop = FALSE]
colnames(umap_sub) <- c("#UMAP1", "#UMAP2")
peakVI_sub <- peakVI[common_cells, , drop = FALSE]
leiden_sub <- leiden[common_cells]
leiden_df <- data.frame(leiden_PeakVI = leiden_sub)
rownames(leiden_df) <- common_cells
proj@embeddings$UMAP <- SimpleList(df = umap_sub,  params = list())
rd <- SimpleList(matRD = peakVI_sub, params = list(corCutOff = NA), scaleDims = FALSE, corToDepth = NA,
  cellColData = proj@cellColData)
proj@reducedDims[["PeakVI"]] <- rd
proj <- addCellColData(ArchRProj = proj,data = leiden_df$leiden_PeakVI,
  name = "leiden_PeakVI",cells = rownames(leiden_df),force = TRUE)

# Transfer Multiome Vascular Cell Type & VSMC Cell Type Labels
rna_object <- "/COMMA/results/scrna/downstream_analysis/subclustering/vascular/vascular_annotated_small.h5ad"
seurat_object <- schard::h5ad2seurat(rna_object)
metadata <- seurat_object[[]] %>%
  filter(sequencing_method == "10X Multiome")
rna_cells <- data.frame(cell_name = rownames(metadata)) %>%
  mutate(rna_barcode = str_extract(cell_name, "[^_]+$"),fastq_name = str_remove(cell_name, "_[^_]+$")) 
archr_mapping <- proj@cellColData %>%
  as.data.frame() %>% rownames_to_column(var = "ATAC_cell") %>%
  mutate(ATAC_sample = sub("#.*", "", ATAC_cell)) %>% select(ATAC_sample, linked_RNA) %>%
  group_by(linked_RNA) %>% summarise(atac_sample_name = first(ATAC_sample), .groups = "drop")
rna_cells <- rna_cells %>% left_join(archr_mapping, by = c("fastq_name" = "linked_RNA"))
rna_cells <- rna_cells %>% filter(!is.na(atac_sample_name)) %>% mutate(new_cell_name = paste0(atac_sample_name, "#", rna_barcode, "-1"))
final_rna_metadata <- seurat_object[[]] %>% rownames_to_column("cell_name") %>% filter(sequencing_method == "10X Multiome")
new_cell_names <- merge(final_rna_metadata, rna_cells, by="cell_name")
archr_df <- proj@cellColData %>% as.data.frame() %>%
  rownames_to_column(var = "cell_name")
merged_metadata <- merge(archr_df, new_cell_names, by.x = "cell_name", by.y = "new_cell_name", all.x=T)
merged_metadata <- merged_metadata[match(rownames(proj@cellColData), merged_metadata$cell_name), ]
proj@cellColData$multiome_vascular_celltype <- merged_metadata$vascular_cell_type
proj@cellColData$multiome_vsmc_celltype <- merged_metadata$vsmc_cell_type

# Filter Multiome Only Cells
valid_cells <- proj$cellNames[!is.na(proj$multiome_vascular_celltype)]
proj_subset <- subsetArchRProject(ArchRProj = proj, cells = valid_cells,
  outputDirectory = "/COMMA/results/scatac/ArchR/vascular_paired", dropCells = TRUE, force = TRUE)

# Plot UMAPs of Clusters & Vascular Cell Type
p1 <- plotEmbedding(ArchRProj = proj_subset, colorBy = "cellColData", name = c("multiome_vascular_celltype"), embedding = "UMAP", size = 0.05)
ggsave(plot = p1, file="umap_multiome_vascular_celltype.png",  width =9, height =7)
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = c("leiden_PeakVI"), embedding = "UMAP", size = 0.05)
ggsave(plot = p2, file="umap_vascular_cluster.png")

# Define Differentially Expressed Genes
markersGS <- getMarkerFeatures(ArchRProj = proj, useMatrix = "GeneScoreMatrix", groupBy = "leiden_PeakVI",
    bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.5")
markerResult <- imap_dfr(markerList, ~ {as.data.frame(.x) %>%  
  tibble::rownames_to_column("Group") %>% mutate(group = .y)})
write.csv(markerResult, "leiden_cluster_markergenes.csv")

# Assign Clusters with 40% of Cells have a High Doublet Enrichment Score
doublet_score <- getCellColData(proj, "DoubletEnrichment")[, 1]
leiden <- getCellColData(proj, "leiden_PeakVI")[, 1]
top10_cutoff <- quantile(doublet_score, probs = 0.9, na.rm = TRUE)
cut_by_cluster <- split(doublet_score, leiden)
high_cut_frac <- sapply(cut_by_cluster,
  function(x) mean(x >= top10_cutoff, na.rm = TRUE))
doublet_cluster_ids <- names(high_cut_frac[high_cut_frac >= 0.4])
doublet_flag <- leiden %in% doublet_cluster_ids
names(doublet_flag) <- proj$cellNames
proj <- addCellColData(ArchRProj = proj, data = doublet_flag, name = "doublet_cluster", cells = names(doublet_flag), force = TRUE)

# Map Clusters to New Names
newClusterNames <- c("Pericyte",
  "Capillary Endothelial Cell", 
  "Capillary Endothelial Cell",
  "Arterial Endothelial Cell",
  "Pericyte",
  "Remove", # Remove Due to Absence of Vascular Marker Genes
  "Venous Endothelial Cell",
  "Vascular Smooth Muscle Cell",
  "Remove", # Removed Due to Absence of Strong Marker Genes
  "Pericyte",
  "Pericyte"
)
oldClusterNames <- as.character(sort(as.numeric(unique(proj$leiden_PeakVI))))
proj$vascular_cell_type <- mapLabels(proj$leiden_PeakVI,oldLabels = oldClusterNames,newLabels = newClusterNames)

# Remove Cells with High Doublet Scores & Assigned as Remove
keep_cells <- getCellNames(proj)[
  proj$vascular_cell_type != "Remove" & !as.logical(proj$doublet_cluster)]
proj_clustered <- subsetArchRProject(ArchRProj = proj, cells = keep_cells,
  outputDirectory = "/COMMA/results/scatac/ArchR/vascular_clustered", dropCells = TRUE, force = TRUE)

# Recompute UMAP after Cell Removal
proj_clustered <- addUMAP(ArchRProj = proj_clustered, reducedDims = "PeakVI", name = "UMAP", nNeighbors = 30, minDist = 0.2, force = TRUE)
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Study", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "vascular_cell_type", embedding = "UMAP")
p_combined <- p1 + p2
ggsave(filename = "UMAP_Study_and_vascularCellType.png", plot = p_combined,width = 12,height = 6,dpi = 300)

# Add Group Coverages
proj_clustered <- addGroupCoverages(ArchRProj = proj_clustered, groupBy = "vascular_cell_type", sampleLabels = "Sample", minReplicates = 3, maxReplicates = 3, minCells = 40, maxCells = 2000, force = TRUE)
proj_clustered <- addReproduciblePeakSet(
  ArchRProj = proj_clustered, 
  groupBy = "vascular_cell_type", 
  threads = 10,
  reproducibility = "2",
  pathToMacs2 = 'conda/skylake/envs/MACS2/bin/macs2')

# Save ArchR Project
saveArchRProject(proj_clustered, '/COMMA/results/scatac/ArchR/vascular_clustered', load = FALSE)

# Export to H5AD
obs <- getCellColData(proj_clustered)     # ArchRCellColData object
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

outfile = "/COMMA/results/scatac/vascular_clustered.h5ad"
adata$write_h5ad(outfile)
