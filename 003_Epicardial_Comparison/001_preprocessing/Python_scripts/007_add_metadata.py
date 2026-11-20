import pandas as pd
import anndata as ad
import scanpy as sc
import os
import argparse
import random

random.seed(1234)

os.chdir("/COMMA/results/scrna/downstream_analysis/epicardial_coronary/count_matrix/")

adata = sc.read_h5ad("rna_count_matrix.h5ad")
# Supplemental Table 7
metadata = pd.read_excel("/COMMA/results/scrna/downstream_analysis/epicardial_coronary/count_matrix/Zhao_2025_metadata.xlsx", sheet_name=0)

adata.obs.index = adata.obs.index.str.replace(r'^processed_', '', regex=True)
adata.obs['fastq_name'] = adata.obs['fastq_name'].str.replace(r'^processed_', '', regex=True)

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
adata.write("rna_count_matrix_metadata.h5ad")

