library(SummarizedExperiment)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SingleCellExperiment)
library(anndata)
library(Seurat, lib.loc = "R/4.3/skylake/seurat_v5.3.0")
Sys.setenv(RETICULATE_PYTHON = "/conda/skylake/envs/r-reticulate/bin/python")
library(harmony)
library(reticulate)
library(optparse)
library(schard)

setwd("/COMMA/results/scrna/downstream_analysis/subclustering/vascular")

seurat_object <- schard::h5ad2seurat("vascular_hvg.h5ad")

seurat_object <- DietSeurat(seurat_object, dimreducs = NULL, layers = c("counts", "data"))
seurat_object <- ScaleData(seurat_object, vars.to.regress = c('pct_counts_mt','pct_counts_ribo'))
seurat_object <- RunPCA(seurat_object, npcs=30, reduction.name = "X_pca", features=rownames(seurat_object)) # Run PCA in Seurat
integrated_rds <- RunHarmony(seurat_object, c("Study", "Sample Name", "sample_type"), reduction="X_pca", reduction.save = "X_harmony")

# Addressing quirks of writing to anndata
if ("_index" %in% colnames(integrated_rds[["RNA"]]@meta.features)) {
  integrated_rds[["RNA"]]@meta.features <- integrated_rds[["RNA"]]@meta.features[, setdiff(colnames(integrated_rds[["RNA"]]@meta.features), "_index")]}
if ("_index" %in% colnames(integrated_rds@meta.data)) {
  integrated_rds@meta.data <- integrated_rds@meta.data[, setdiff(colnames(integrated_rds@meta.data), "_index")]}

sceasy::convertFormat(integrated_rds, from="seurat", to="anndata",
                       outFile='harmony_vascular_hvg.h5ad') # Convert Integrated object back to anndata format