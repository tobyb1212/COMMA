import numpy as np
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import classification_report
from scarches.models.scpoli import scPoli
import os
import argparse
import random

random.seed(1234)

os.chdir("/COMMA/results/scrna/downstream_analysis/subclustering/vascular")

# Subset anndata object for vascular cell types for Subclustering
adata_integrated = sc.read_h5ad("vascular.h5ad")
adata_hvg = sc.read_h5ad("vascular_hvg.h5ad")

# scPoli Integration
early_stopping_kwargs = { "early_stopping_metric": "val_prototype_loss","mode": "min","threshold": 0,"patience": 20, "reduce_lr": True,"lr_patience": 13,"lr_factor": 0.1}
condition_keys = ['Study', 'Sample Name', "sample_type"]
cell_type_key = 'broad_cell_type'
adata_ref = adata_hvg[adata_hvg.obs['broad_cell_type'] != 'Unknown'].copy()

# Setup scPoli Reference Model
scpoli_model = scPoli(adata=adata_ref,condition_keys=condition_keys,cell_type_keys=cell_type_key,embedding_dims=5,recon_loss='nb',)
scpoli_model.train(n_epochs=50,pretraining_epochs=40,early_stopping_kwargs=early_stopping_kwargs,eta=5)
scpoli_model.save("reference_model.pt")
labeled_indices = adata_hvg.obs.index[adata_hvg.obs['broad_cell_type'] != 'Unknown']
labeled_indices = [adata_hvg.obs.index.get_loc(idx) for idx in labeled_indices]

# Setup scPoli Query Model
scpoli_query = scPoli.load_query_data(adata=adata_hvg,reference_model=scpoli_model, labeled_indices=labeled_indices)
scpoli_query.train(n_epochs=50,pretraining_epochs=40,eta=10)
scpoli_query.save("query_model.pt")

# Return Latent Embeddings
for key in condition_keys:
    adata_hvg.obs[key] = np.array(adata_hvg.obs[key])
adata_hvg.obsm['X_scPoli'] = scpoli_query.get_latent(adata_hvg,mean=True)
adata_integrated.obsm["X_scPoli"] = adata_hvg.obsm["X_scPoli"]

# Save AnnData Object 
adata_integrated.write_h5ad("scPoli_vascular_integrated.h5ad")
