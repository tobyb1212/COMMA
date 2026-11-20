import pandas as pd
import numpy as np
import scvi
import scanpy as sc
import os
import argparse
import gc
import random

os.chdir("/COMMA/results/scrna/downstream_analysis/epicardial_coronary/integrated")

adata_integrated = sc.read_h5ad("/COMMA/results/scrna/downstream_analysis/epicardial_coronary/integrated/preprocessed.h5ad")
adata_hvg = sc.read_h5ad("/COMMA/results/scrna/downstream_analysis/epicardial_coronary/integrated/preprocessed_hvg.h5ad")
adata_vascular = sc.read_h5ad("/COMMA/results/scrna/downstream_analysis/subclustering/vascular/vascular_annotated_final.h5ad")

# Set epicardial Study & Variables
adata_integrated.obs['Study'] = "Zhao 2025"
adata_hvg.obs['Study'] = "Zhao 2025"
adata_integrated.obs['sample_type'] = "Cells"
adata_hvg.obs['sample_type'] = "Cells"
adata_integrated.obs['broad_cell_type'] = "Unknown"
adata_hvg.obs['broad_cell_type'] = "Unknown"
adata_integrated.obs['vascular_cell_type'] = "Unknown"
adata_hvg.obs['vascular_cell_type'] = "Unknown"

# Intersections of Variables & Metadata between epicardial and Atlases
common_vascular_vars = adata_vascular.var_names.intersection(adata_integrated.var_names)
common_vascular_obs_cols = adata_vascular.obs.columns.intersection(adata_integrated.obs.columns)
common_vascular_var_cols = adata_vascular.var.columns.intersection(adata_integrated.var.columns)

# scVI & scANVI of Epicardial Coronary Cells with Microvascular Atlas
# Concatenate the Objects Together
a1 = adata_vascular[:, common_vascular_vars].copy()
a2 = adata_integrated[:, common_vascular_vars].copy()
a1.obs = a1.obs[common_vascular_obs_cols].copy(); a2.obs = a2.obs[common_vascular_obs_cols].copy()
a1.var = a1.var[common_vascular_var_cols].copy(); a2.var = a2.var[common_vascular_var_cols].copy()
adata_vascular_hvg = sc.concat([a1, a2], label="batch", keys=["atlas","epicardial"])
sc.pp.highly_variable_genes(adata_vascular_hvg, flavor="cell_ranger", n_top_genes=2000, layer="counts", batch_key="Study")
adata_vascular_hvg = adata_vascular_hvg[:, adata_vascular_hvg.var.highly_variable].copy()

#Construct scVI Model
scvi.model.SCVI.setup_anndata(adata_vascular_hvg, layer="counts",batch_key="Study", categorical_covariate_keys=["Sample Name", "sample_type"],continuous_covariate_keys=['total_counts',"pct_counts_mt", "pct_counts_ribo"])
model = scvi.model.SCVI(adata_vascular_hvg, gene_likelihood="nb", dispersion ='gene-batch', n_latent=30, n_layers=2)
model.view_anndata_setup(adata_vascular_hvg)
#model.train()
#model.save(f'vascular_model',overwrite=True)
model = scvi.model.SCVI.load("vascular_model", adata = adata_vascular_hvg)
adata_vascular_hvg.obsm["X_scVI"] = model.get_latent_representation()
sc.pp.neighbors(adata_vascular_hvg, use_rep="X_scVI")
sc.tl.umap(adata_vascular_hvg)
sc.pl.umap(adata_vascular_hvg, color=["Study", "Sample Name", "sample_type", "pct_counts_mt", "vascular_cell_type"], ncols=3,  save=f"scVI_vascular.png")
# Construct scANVI Model
lvae = scvi.model.SCANVI.from_scvi_model(model, adata=adata_vascular_hvg, labels_key="vascular_cell_type", unlabeled_category="Unknown")
#lvae.train(max_epochs=20, n_samples_per_label=100)
lvae = scvi.model.SCANVI.load("vascular_scanvi_model", adata = adata_vascular_hvg)
#lvae.save(f'vascular_scanvi_model',overwrite=True)
adata_vascular_hvg.obsm["X_scANVI"] = lvae.get_latent_representation()
epicardial_idx = adata_vascular_hvg.obs['batch'] == 'epicardial'
adata_vascular_cropped = adata_vascular_hvg[epicardial_idx].copy()
sc.pp.neighbors(adata_vascular_hvg, use_rep="X_scANVI")
sc.tl.umap(adata_vascular_hvg)
sc.pl.umap(adata_vascular_hvg, color=["Study", "Sample Name", "sample_type", "pct_counts_mt", "vascular_cell_type"], ncols=3,  save=f"scANVI_vascular.png")
# Transfer Embeddings
adata_integrated.obsm["X_scVI_vascular"] = adata_vascular_cropped.obsm["X_scVI"]
adata_integrated.obsm["X_scANVI_vascular"] = adata_vascular_cropped.obsm["X_scANVI"]

# Plot Embeddings of epicardial Data Alone
for emb in adata_integrated.obsm.keys():
    sc.pp.neighbors(adata_integrated, use_rep=emb)
    sc.tl.umap(adata_integrated)
    sc.pl.umap(adata_integrated, color=["Study", "Sample Name", "sample_type", "pct_counts_mt"], ncols=2,  save=f"{emb}_epicardial.png")

# Save Integrated Objects
adata_integrated.write_h5ad("epicardial_all_embeddings.h5ad")
adata_vascular_hvg.write_h5ad("vascular_epicardial.h5ad")
