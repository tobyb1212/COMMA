import pandas as pd
import anndata
import scanpy as sc
import os
import numpy as np
import random

random.seed(1234)

adata = sc.read_h5ad("/COMMA/results/scrna/count_matrix/pre_filter_combined_count_matrix.h5ad")

# Filter for Sequencing Saturation
saturation_threshold=0.2
filtered_cells = adata.obs[adata.obs['Sequencing Saturation'] <= saturation_threshold]
n_filtered = filtered_cells.shape[0]
samples_filtered = filtered_cells["fastq_name"].nunique()
filtered_cells[["fastq_name"]].drop_duplicates().to_csv("/COMMA/results/scrna/count_matrix/filtering/saturation_filter.csv", index=False) # Recording fastq_names that were Filtered
adata_filtered = adata[adata.obs['Sequencing Saturation']>saturation_threshold]
print(f"Cells filtered for sequencing saturation ≤ 0.2: {n_filtered}")
print(f"Samples filtered for sequencing saturation ≤ 0.2: {samples_filtered}")

# Filter for Q30 Bases in RNA Read
q30_threshold=0.7
filtered_cells = adata_filtered.obs[adata_filtered.obs['Q30 Bases in RNA read'] <= q30_threshold]
n_filtered = filtered_cells.shape[0]
samples_filtered = filtered_cells["fastq_name"].nunique()
filtered_cells[["fastq_name"]].drop_duplicates().to_csv("/COMMA/results/scrna/count_matrix/filtering/q30_filter.csv", index=False) # Recording fastq_names that were Filtered
adata_filtered = adata_filtered[adata_filtered.obs['Q30 Bases in RNA read']>q30_threshold]
print(f"Cells filtered for percentage Q30 bases in RNA read ≤ 0.7: {n_filtered}")
print(f"Samples filtered for percentage Q30 bases in RNA read ≤ 0.7: {samples_filtered}")

# Filter for Percentage of Reads Uniquely Mapped to Transcriptome
unique_threshold=0.3
filtered_cells = adata_filtered.obs[adata_filtered.obs['Reads Mapped to GeneFull_Ex50pAS: Unique GeneFull_Ex50pAS'] <= unique_threshold]
n_filtered = filtered_cells.shape[0]
samples_filtered = filtered_cells["fastq_name"].nunique()
filtered_cells[["fastq_name"]].drop_duplicates().to_csv("/COMMA/results/scrna/count_matrix/filtering/unique_filter.csv", index=False) # Recording fastq_names that were Filtered
adata_filtered = adata_filtered[adata_filtered.obs['Reads Mapped to GeneFull_Ex50pAS: Unique GeneFull_Ex50pAS']>unique_threshold]
print(f"Cells filtered for percentage reads uniquely mapped to the transcriptome ≤ 0.3: {n_filtered}")
print(f"Samples filtered for percentage reads uniquely mapped to the transcriptome ≤ 0.3: {samples_filtered}")

# Filter for Doublets
doublet_threshold=0.3
n_filtered = (adata_filtered.obs['scrublet_score'] >= doublet_threshold).sum()
proportion_doublets = np.sum(adata_filtered.obs['scrublet_score']>=doublet_threshold)/len(adata_filtered)
print(f'Proportion of doublet barcodes: {proportion_doublets}')
adata_filtered = adata_filtered[adata_filtered.obs['scrublet_score']<doublet_threshold]
print(f"Cells filtered for Doublet score >0.3: {n_filtered}")

# Filter for Mitochondrial Gene Percentage
cell_mask = adata_filtered.obs['sample_type'].isin(["CD45+ cells", "CD45- cells", "Cells"])
nuclei_mask = (adata_filtered.obs['sample_type'] == "Nuclei")
cell_filtered = ((adata_filtered.obs['pct_counts_mt'] > 20) & cell_mask).sum()
nuclei_filtered = ((adata_filtered.obs['pct_counts_mt'] > 5) & nuclei_mask).sum()
adata_filtered = adata_filtered[~((adata_filtered.obs['pct_counts_mt'] > 20) & cell_mask) & ~((adata_filtered.obs['pct_counts_mt'] > 5) & nuclei_mask)] 
print(f"Cells filtered for mitochondrial gene percentage >20%: {cell_filtered}")
print(f"Nuclei filtered for mitochondrial gene percentage >5%: {nuclei_filtered}")

# Filter for Ribosomal Gene Percentage
cell_mask = adata_filtered.obs['sample_type'].isin(["CD45+ cells", "CD45- cells", "Cells"])
nuclei_mask = (adata_filtered.obs['sample_type'] == "Nuclei")
cell_filtered = ((adata_filtered.obs['pct_counts_ribo'] > 20) & cell_mask).sum()
nuclei_filtered = ((adata_filtered.obs['pct_counts_ribo'] > 5) & nuclei_mask).sum()
adata_filtered = adata_filtered[~((adata_filtered.obs['pct_counts_ribo'] > 20) & cell_mask) & ~((adata_filtered.obs['pct_counts_ribo'] > 5) & nuclei_mask)] 
print(f"Cells filtered for ribosomal gene percentage >20%: {cell_filtered}")
print(f"Nuclei filtered for ribosomal gene percentage >5%: {nuclei_filtered}")

# Filter for Haemoglobin Gene Percentage
haem_threshold=3
n_filtered = (adata_filtered.obs['pct_counts_hb'] >= haem_threshold).sum()
adata_filtered = adata_filtered[adata_filtered.obs['pct_counts_hb']<haem_threshold]
print(f"Cells filtered for haemoglobin gene percentage >3%: {n_filtered}")

adata_filtered.write("/COMMA/results/scrna/count_matrix/filtered_count_matrix.h5ad")

print(f"Total number of cells retained after all filters: {adata_filtered.n_obs}")