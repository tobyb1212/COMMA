import pandas as pd
import anndata
import scanpy as sc
import os
import random

random.seed(1234)

adata = sc.read_h5ad("/COMMA/results/scrna/count_matrix/rna_count_matrix.h5ad")
# Supplemental Table 2
metadata = pd.read_excel("/COMMA/results/scrna/metadata/FINAL_METADATA.xlsx", sheet_name=0)

if 'fastq_name' not in adata.obs.columns:
    raise ValueError("'fastq_name' not found in adata.obs")
if 'fastq_name' not in metadata.columns:
    raise ValueError("'fastq_name' not found in metadata")

adata.obs['__orig_index'] = adata.obs.index
adata.obs = adata.obs.merge(metadata, on='fastq_name', how='left')
adata.obs.index = adata.obs['__orig_index']
adata.obs.drop(columns='__orig_index', inplace=True)

print(adata.obs.head(10))
print(adata.shape)
adata.write("/COMMA/results/scrna/count_matrix/rna_count_matrix_metadata.h5ad")