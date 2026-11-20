import pandas as pd
import anndata
import scanpy as sc
import os
import random

random.seed(1234)
os.chdir("/COMMA/results/scrna/integrated/pre_integration")

adata = sc.read_h5ad("/COMMA/results/scrna/count_matrix/filtered_count_matrix.h5ad")
adata.obs.rename(columns={'sequencinq_method': 'sequencing_method'}, inplace=True)

# Filter samples with less than 100 cells
sample_counts = adata.obs['Sample Name'].value_counts()
valid_samples = sample_counts[sample_counts >= 100].index
adata = adata[adata.obs['Sample Name'].isin(valid_samples)]

# Filter unneccessary adata.var
adata.var = adata.var[['mt', 'ribo', 'hb', 'n_cells']]

# Re-Filter Genes
sc.pp.filter_genes(adata, min_cells=50, inplace=True)

# Preserve Raw Counts (CellBender Counts) then normalise and logarithmise data
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Select Highly Variable Genes & Plot Graph of HVGs
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000, layer="counts", batch_key="Study")
sc.pl.highly_variable_genes(adata, save = "highlyvariablegenes.png")
adata_hvg = adata[:, adata.var.highly_variable]

# Perform PCA, Nearest Neighbour and UMAP
sc.pp.pca(adata_hvg)
sc.pl.pca_variance_ratio(adata_hvg, log=True, save = "elbowplot.png")
sc.pp.neighbors(adata_hvg, n_neighbors=30)
sc.tl.umap(adata_hvg)

# Transfer Embeddings & Connectivities from HVG Object to Full Adata Object
adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
adata.obsp = adata_hvg.obsp

# Plot Pre-Integration UMAPs
sc.pl.umap(adata, color=["Study", "Patient", "Sample Name", "sequencing_method", "sample_type", "Region", "total_counts", "pct_counts_mt", "pct_counts_ribo"], ncols=2, save="umaps.png")

# Save Both Full .h5ad Object and HVG .h5ad Object
adata.write_h5ad("preprocessed.h5ad")
adata_hvg.write_h5ad("preprocessed_hvg.h5ad")

adata
adata_hvg
