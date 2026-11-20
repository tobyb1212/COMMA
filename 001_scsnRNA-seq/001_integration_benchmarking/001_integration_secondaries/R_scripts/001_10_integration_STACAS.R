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
library(STACAS)

option_list = list(
  make_option(c("-r", "--raw"), type="character", default=NULL, 
              help="Non-Integrated .h5ad object", metavar="character"),
    make_option(c("-v", "--hvg"), type="character", default=NULL, 
              help="HVG-selected .h5ad file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

dir.create("/COMMA/results/scrna/integrated/STACAS", showWarnings = FALSE)
setwd("/COMMA/results/scrna/integrated/STACAS")

set.seed(1234)

use_condaenv("r-reticulate")

print("Loaded Conda Environment")

ad_path <- opt$raw 
hvg_path <- opt$hvg

seurat_object <- schard::h5ad2seurat(hvg_path)
seurat_counts <- schard::h5ad2seurat("/COMMA/results/scrna/integrated/pre_integration/preprocessed_hvg_counts.h5ad")
seurat_object@assays$RNA@counts <- seurat_counts@assays$RNA@counts

# Add HCA Cell Types for Semi-Supervised Integration
df <- readr::read_csv("/COMMA/results/scrna/metadata/HCA_celltypes.csv")
metadata <- seurat_object@meta.data %>%
  rownames_to_column(var = "cell_id")
metadata <- metadata %>%
  left_join(df, by = c("cell_id" = "cell_name")) %>%
  mutate(HCA_cell_type = if_else(is.na(cell_type), "Unknown", cell_type)) %>%
  select(-cell_type)
rownames(metadata) <- metadata$cell_id
metadata <- metadata %>% select(-cell_id)
metadata <- metadata %>%
  mutate(HCA_cell_type = na_if(HCA_cell_type, "Unknown"))
seurat_object@meta.data <- metadata

# Formatting Seurat Object Correctly for ssSTACAS
seurat_object@assays$RNA@meta.features$highly_variable <- rownames(seurat_object[["RNA"]]) %in% rownames(seurat_object[["RNA"]]) 
seurat_object[["RNA"]]@var.features <- rownames(seurat_object[["RNA"]]) 
seurat_object[["pca"]] <- seurat_object[["Xpca_"]] 
seurat_object <- DietSeurat(seurat_object, dimreducs = c("pca"), layers = c("counts", "data"))

# Run semi-supervised STACAS for Sample Name Batches
obj.list <- SplitObject(seurat_object, split.by = "Sample Name")
integrated_rds <- obj.list %>%
  Run.STACAS(dims = 1:30, anchor.features = 2000, cell.labels = "HCA_cell_type", min.sample.size=1)
integrated_rds[["X_stacas"]] <- integrated_rds[["pca"]]
integrated_rds[["pca"]] <- NULL

# Addressing quirks of writing to anndata
if ("_index" %in% colnames(integrated_rds[["RNA"]]@meta.features)) {
  integrated_rds[["RNA"]]@meta.features <- integrated_rds[["RNA"]]@meta.features[, setdiff(colnames(integrated_rds[["RNA"]]@meta.features), "_index")]}
if ("_index" %in% colnames(integrated_rds@meta.data)) {
  integrated_rds@meta.data <- integrated_rds@meta.data[, setdiff(colnames(integrated_rds@meta.data), "_index")]}

sceasy::convertFormat(integrated_rds, from="seurat", to="anndata",
                       outFile='STACAS_integrated_hvg.h5ad') # Convert Integrated object back to anndata format