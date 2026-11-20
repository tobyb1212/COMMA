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

proj <- loadArchRProject(path = '/COMMA/results/scatac/ArchR/initial', showLogo = FALSE)

# Plotting TSS Enrichment Scores and Number of Unique Fragments
p1 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Study", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE,
    baseSize = 10)

p2 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Study", 
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE,
    baseSize = 10)

p3 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Study", 
    name = "DoubletScore",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE,
    baseSize = 10)
p3 + ylim(0, 10)

# Saving Plots
ggsave(filename = "/COMMA/results/scatac/QC/QC_TSSE.pdf",plot = p1,width = 8,height = 5)
ggsave(filename = "/COMMA/results/scatac/QC/QC_TSSE.png",plot = p1,width = 8,height = 5)
ggsave(filename = "/COMMA/results/scatac/QC/QC_log10fragments.pdf",plot = p2,width = 8,height = 5)
ggsave(filename = "/COMMA/results/scatac/QC/QC_log10fragments.png",plot = p2,width = 8,height = 5)
ggsave(filename = "/COMMA/results/scatac/QC/QC_doublet_score.pdf",plot = p3+ylim(0,20),width = 8,height = 5)
ggsave(filename = "/COMMA/results/scatac/QC/QC_doublet_score.png",plot = p3+ylim(0,20),width = 8,height = 5)


df <- data.frame(getCellColData(proj)[, c("DoubletScore", "Study")])
colnames(df) <- c("DoubletScore", "Study")

# Plotting Doublet Scores for Each Study and Indicating Cut-Off Threshold
for (study_name in unique(df$Study)) {
  df_subset <- df %>% filter(Study == study_name)
  
  p <- ggplot(df_subset, aes(x = DoubletScore)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "black") + xlim(-1, 25) + 
    geom_vline(xintercept = 0.5, color = "red", linetype = "dotted", linewidth = 1) +
    labs(x = "Doublet score",
         y = "Cell count") +
    theme_bw()
  
  # Sanitize study_name for filenames (removing spaces) then saving
  clean_name <- gsub("[^A-Za-z0-9_]", "_", study_name)
  filename <- paste0(clean_name, "_doublet_score.png")
  outpath <- file.path("/COMMA/results/scatac/QC/", filename)
  ggsave(filename = outpath, plot = p, width = 6, height = 4, dpi = 300)
}

# Export to H5AD
obs <- getCellColData(proj)     # ArchRCellColData object
obs <- as.data.frame(obs)
var = as.data.frame("1")

adata <- AnnData(
  obs = obs,
  var = var)
outfile = "/COMMA/results/scatac/figures/atac_qc.h5ad"
adata$write_h5ad(outfile)

test <- readRDS("/COMMA/results/scatac/matrices/Peak_matrix_final_integrated.rds")
obs <- test@colData     # ArchRCellColData object
obs <- as.data.frame(obs)
umap <- reducedDim(test, "X_umap")
var = as.data.frame("1")
X <- matrix(0L, nrow = nrow(obs), ncol = length(var))

adata <- AnnData(
  X = X,
  obs = obs,
  var = var)
adata$obsm$X_umap <- umap

outfile = "/COMMA/results/scatac/figures/unfiltered_integrated_whole.h5ad"
adata$write_h5ad(outfile)

