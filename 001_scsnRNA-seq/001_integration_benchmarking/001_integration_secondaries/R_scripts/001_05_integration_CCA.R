library(SummarizedExperiment)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SingleCellExperiment)
library(anndata)
library(Seurat, lib.loc = "R/4.3/skylake/seurat_v5.3.0")
Sys.setenv(RETICULATE_PYTHON = "/conda/skylake/envs/r-reticulate/bin/python")
library(reticulate)
library(optparse)
library(schard)

option_list = list(
  make_option(c("-r", "--raw"), type="character", default=NULL, 
              help="Non-Integrated .h5ad object", metavar="character"),
    make_option(c("-v", "--hvg"), type="character", default=NULL, 
              help="HVG-selected .h5ad file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

dir.create("/COMMA/results/scrna/integrated/CCA", showWarnings = FALSE)
setwd("/COMMA/results/scrna/integrated/CCA")

set.seed(1234)

use_condaenv("r-reticulate")

print("Loaded Conda Environment")

ad_path <- opt$raw 
hvg_path <- opt$hvg

seurat_object <- schard::h5ad2seurat(hvg_path)
seurat_counts <- schard::h5ad2seurat("/COMMA/results/scrna/integrated/pre_integration/preprocessed_hvg_counts.h5ad")
seurat_object@assays$RNA@counts <- seurat_counts@assays$RNA@counts

seurat_object <- DietSeurat(seurat_object, dimreducs = NULL, layers = c("counts", "data")) # Slim down seurat to remove previous PCA but keep data
seurat_object[["RNA"]] <- split(seurat_object[["RNA"]], f = seurat_object$'Sample Name')
seurat_object <- ScaleData(seurat_object, vars.to.regress = c('pct_counts_mt','pct_counts_ribo')) # Scale data and regress out percentage mitochondrial and ribosomal genes
VariableFeatures(seurat_object) <- rownames(seurat_object)
seurat_object <- RunPCA(seurat_object, npcs=30, reduction.name = "X_pca", features=rownames(seurat_object)) # Run PCA in Seurat

integrated_embeddings <- IntegrateLayers(object = seurat_object, method = CCAIntegration, orig.reduction = "X_pca", normalization.method = "LogNormalize", new.reduction = "X_cca", verbose = FALSE) # Integrate layers with CCA

# Conversion back to anndata is difficult for v5 object, so transferring v5 embeddings to new v3 object with raw counts
integrated_rds <- schard::h5ad2seurat(hvg_path)
integrated_rds <- DietSeurat(integrated_rds, dimreducs = NULL, layers = c("counts"))
integrated_rds[["X_pca"]] <- integrated_embeddings[["X_pca"]]
integrated_rds[["X_cca"]] <- integrated_embeddings[["X_cca"]]

# Addressing quirks of writing to anndata
if ("_index" %in% colnames(integrated_rds[["RNA"]]@meta.features)) {
  integrated_rds[["RNA"]]@meta.features <- integrated_rds[["RNA"]]@meta.features[, setdiff(colnames(integrated_rds[["RNA"]]@meta.features), "_index")]}
if ("_index" %in% colnames(integrated_rds@meta.data)) {
  integrated_rds@meta.data <- integrated_rds@meta.data[, setdiff(colnames(integrated_rds@meta.data), "_index")]}

sceasy::convertFormat(integrated_rds, from="seurat", to="anndata",
                       outFile='cca_integrated_hvg.h5ad') # Convert Integrated object back to anndata format