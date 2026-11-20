import pandas as pd
import anndata
import scanpy as sc
import os
import random
import scvi

random.seed(1234)
scvi.settings.seed = 1

vsmc_clusters=['Vascular Smooth Muscle Cell']
celltype_dir = "vsmc"
os.chdir("/COMMA/downstream_analysis/subclustering")
os.makedirs(celltype_dir, exist_ok=True)
os.chdir(celltype_dir)

# Subset anndata object for vsmc cell types for Subclustering
adata = sc.read_h5ad("//COMMA/results/scrna/downstream_analysis/subclustering/vascular/vascular_annotated_final.h5ad")
adata = adata[adata.obs["vascular_cell_type"].isin(vsmc_clusters)].copy()

# Filter samples with less than 20 cells
sample_counts = adata.obs['Sample Name'].value_counts()
valid_samples = sample_counts[sample_counts >= 20].index
adata = adata[adata.obs['Sample Name'].isin(valid_samples)]
adata = adata.copy()

# Re-Filter Genes
sc.pp.filter_genes(adata, min_cells=5, inplace=True)

# Normalise and logarithmise data
adata.X = adata.layers['counts'].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Select Highly Variable Genes & Plot Graph of HVGs
# Due to Small Size of Some Batches and Fragility of LOESS / Dispersion Model, using Cell_Ranger HVG Selection
sc.pp.highly_variable_genes(adata, flavor="cell_ranger", n_top_genes=2000, batch_key="Study")
sc.pl.highly_variable_genes(adata, save = "highlyvariablegenes.png")
adata_hvg = adata[:, adata.var.highly_variable]
adata_hvg = adata_hvg.copy()
adata_integrated = sc.read_h5ad("vsmc.h5ad")
adata_hvg = sc.read_h5ad("vsmc_hvg.h5ad")

# Construct scVI model for Subclustered Object on HVG Object & Transfer to Full Object
scvi.model.SCVI.setup_anndata(adata_hvg, layer="counts",batch_key="Study", categorical_covariate_keys=["Sample Name", "sample_type"],continuous_covariate_keys=['total_counts',"pct_counts_mt", "pct_counts_ribo"])
model = scvi.model.SCVI(adata_hvg, gene_likelihood="nb", dispersion ='gene-batch', n_latent=30, n_layers=2)
model.view_anndata_setup(adata_hvg)
model.train()
model.save(f'model.pt',overwrite=True)
adata_hvg.obsm["X_scVI"] = model.get_latent_representation()
adata_integrated.obsm["X_scVI"] = adata_hvg.obsm["X_scVI"]

adata_integrated.write_h5ad("scVI_vsmc_integrated.h5ad")


