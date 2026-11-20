import pandas as pd
import anndata
import scanpy as sc
import os
import random
import scvi

random.seed(1234)
scvi.settings.seed = 1

os.chdir("/COMMA/results/scrna/downstream_analysis/subclustering/vascular")

adata_integrated = sc.read_h5ad("vascular.h5ad")
adata_hvg = sc.read_h5ad("vascular_hvg.h5ad")

# Construct scVI model for Subclustered Object on HVG Object & Transfer to Full Object
scvi.model.SCVI.setup_anndata(adata_hvg, layer="counts",batch_key="Study", categorical_covariate_keys=["Sample Name", "sample_type"],continuous_covariate_keys=['total_counts',"pct_counts_mt", "pct_counts_ribo"])
model = scvi.model.SCVI(adata_hvg, gene_likelihood="nb", dispersion ='gene-batch', n_latent=30, n_layers=2)
model.view_anndata_setup(adata_hvg)
model.train()
model.save(f'model.pt',overwrite=True)
adata_hvg.obsm["X_scVI"] = model.get_latent_representation()
adata_integrated.obsm["X_scVI"] = adata_hvg.obsm["X_scVI"]

adata_integrated.write_h5ad("scVI_vascular_integrated.h5ad")

# Construct scANVI model from pre-generated scVI model
lvae = scvi.model.SCANVI.from_scvi_model(model, adata=adata_hvg, labels_key="broad_cell_type",unlabeled_category = "Unknown")
lvae.train(max_epochs=20, n_samples_per_label=100)
adata_hvg.obsm["X_scANVI"] = lvae.get_latent_representation()
adata_integrated.obsm["X_scANVI"] = adata_hvg.obsm["X_scANVI"]

adata_integrated.write_h5ad("scANVI_vascular_integrated.h5ad")


