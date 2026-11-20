import pandas as pd
import anndata
import scanpy as sc
import os
import random
import scvi

random.seed(1234)

celltype = ["Endothelial Cells", "Mural Cells"]
celltype_dir = "vascular"
os.chdir("/COMMA/results/scrna/downstream_analysis/subclustering")
os.makedirs(celltype_dir, exist_ok=True)
os.chdir(celltype_dir)

# Subset anndata object for vascular cell types for Subclustering
adata = sc.read_h5ad("/COMMA/results/scrna/downstream_analysis/integrated/integrated_final.h5ad")
adata = adata[adata.obs["broad_cell_type"].isin(celltype)].copy()

# Filter samples with less than 100 cells
sample_counts = adata.obs['Sample Name'].value_counts()
valid_samples = sample_counts[sample_counts >= 100].index
adata = adata[adata.obs['Sample Name'].isin(valid_samples)]
adata = adata.copy()

# Re-Filter Genes
sc.pp.filter_genes(adata, min_cells=10, inplace=True)

# Normalise and logarithmise data
adata.X = adata.layers['counts'].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Select Highly Variable Genes & Plot Graph of HVGs
# Due to Small Size of Some Batches and Fragility of LOESS / Dispersion Model, using Cell_Ranger HVG Selection
sc.pp.highly_variable_genes(adata, flavor="cell_ranger", n_top_genes=4000, batch_key="Study")
sc.pl.highly_variable_genes(adata, save = "highlyvariablegenes.png")
adata_hvg = adata[:, adata.var.highly_variable]
adata_hvg = adata_hvg.copy()

# Perform PCA, Nearest Neighbour and UMAP
sc.pp.pca(adata_hvg)
sc.pl.pca_variance_ratio(adata_hvg, log=True, save = "elbowplot.png")
sc.pp.neighbors(adata_hvg, n_neighbors=30)
sc.tl.umap(adata_hvg)

# Transfer Embeddings & Connectivities from HVG Object to Full Adata Object
adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
adata.obsp = adata_hvg.obsp

adata.write_h5ad("vascular.h5ad")
adata_hvg.write_h5ad(f"vascular_hvg.h5ad")

