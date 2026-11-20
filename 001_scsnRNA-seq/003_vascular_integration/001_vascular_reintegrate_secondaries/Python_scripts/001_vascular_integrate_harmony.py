import pandas as pd
import anndata
import scanpy as sc
import os
import random
import scvi

os.chdir("/COMMA/results/scrna/downstream_analysis/subclustering/vascular")

# Subset anndata object for vascular cell types for Subclustering
adata_integrated = sc.read_h5ad("vascular.h5ad")
adata_hvg = sc.read_h5ad("harmony_vascular_hvg.h5ad")

# Transfer embeddings from HVG object to full object
adata_integrated.obsm["X_harmony"] = adata_hvg.obsm["X_x_harmony"]

adata_integrated.write_h5ad("harmony_vascular_integrated.h5ad")