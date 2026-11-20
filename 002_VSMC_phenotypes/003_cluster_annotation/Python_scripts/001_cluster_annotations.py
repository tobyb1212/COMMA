import pandas as pd
import numpy as np
import scanpy as sc
import os
import argparse
import random
import matplotlib.pyplot as plt
import matplotlib as mpl
import gc
import math
from collections import defaultdict
import seaborn as sns
import textwrap
from matplotlib.patches import Patch

random.seed(1234)

os.chdir("/COMMA/results/scrna/downstream_analysis/subclustering/vsmc/")

adata = sc.read_h5ad("/COMMA/results/scrna/downstream_analysis/subclustering/vsmc/integration/scVI_vsmc_integrated.h5ad")

sc.pp.neighbors(adata, n_neighbors=45, n_pcs=30, use_rep="X_scVI")
sc.tl.umap(adata, min_dist=0.2)
sc.tl.leiden(adata, resolution=0.3, key_added='leiden_0.3')

# Filtered & Unfiltered Gene Expression
sc.tl.rank_genes_groups(adata, groupby="leiden_0.3", method='wilcoxon')
sc.tl.filter_rank_genes_groups(adata, min_fold_change=0.5)
sc.pl.rank_genes_groups_dotplot(adata, key="rank_genes_groups_filtered", n_genes=5, save=f"_markers.png", dendrogram = False)
df_markers = sc.get.rank_genes_groups_df(adata, group=None, key="rank_genes_groups_filtered")
df_markers.to_csv(f"markers.csv", index=False)

# Perform Final Annotation of Cell States
vsmc_annotations = {
    '0': 'Arteriolar Pericyte-like VSMCs',
    '1': 'Microvascular Modulated VSMCs',
    '2': 'Contractile VSMCs',
    '3': 'Stressed VSMCs',
}
adata.obs['vsmc_cell_type'] = adata.obs['leiden_0.3'].map(vsmc_annotations)
adata.obs['vsmc_cell_type'] = adata.obs['vsmc_cell_type'].astype('category')
adata.obs['vsmc_cell_type'].value_counts()

# Figure 2A
adata.uns["vsmc_cell_type_colors"] = list(map(mpl.colors.to_hex, plt.get_cmap("Dark2").colors[:4]))
sc.pl.umap(adata, color="vsmc_cell_type", show=False, title=None)
plt.savefig("vsmc_umap.png", dpi=900, bbox_inches="tight")
plt.close()
sc.pl.umap(adata[adata.obs.groupby("vsmc_cell_type", group_keys=False).apply(lambda x: x.sample(1, random_state=0)).index], color="vsmc_cell_type", save="vsmc_labels.pdf")

# Final Plots of VSMC Subclusters
sc.tl.rank_genes_groups(adata, groupby="vsmc_cell_type", method='wilcoxon')
sc.tl.filter_rank_genes_groups(adata, min_fold_change=1)
# Figure 2B
sc.pl.rank_genes_groups_dotplot(adata, key="rank_genes_groups_filtered", n_genes=5, save=f"_final_markers.pdf", dendrogram = False)
sc.pl.umap(adata, color='vsmc_cell_type', ncols=3, save=f"_final_celltype.png")
# Supplemental Table 5
df_markers = sc.get.rank_genes_groups_df(adata, group=None, key="rank_genes_groups_filtered")
df_markers = df_markers.rename(columns={"pvals_adj": "FDR"})
df_markers["FDR"] = df_markers["FDR"].replace(0, "< 1x 10 ^ -308")
df_markers.to_csv(f"markers.csv", index=False)

# Supplemental Table 6
vsmc_cell = adata[adata.obs['vsmc_cell_type'].isin(['Microvascular Modulated VSMCs', 'Contractile VSMCs'])].copy()
sc.tl.rank_genes_groups(vsmc_cell, groupby="vsmc_cell_type", method='wilcoxon')
sc.tl.filter_rank_genes_groups(vsmc_cell, min_fold_change=1)
df_markers = sc.get.rank_genes_groups_df(vsmc_cell, group=None, key="rank_genes_groups_filtered")
df_markers.to_csv(f"markers_subcluster.csv", index=False)
df = sc.get.rank_genes_groups_df(vsmc_cell, group='Microvascular Modulated VSMCs', key='rank_genes_groups').dropna().copy()
x = df['logfoldchanges']
y = -np.log10(df['pvals_adj'].clip(lower=1e-308))
# Figure 2E
colors = np.where((df.pvals_adj < 0.05) & (x > 1), '#E66B2E',
         np.where((df.pvals_adj < 0.05) & (x < -1), '#807BB5', 'lightgrey'))
plt.figure(figsize=(5, 4))
plt.scatter(x, y, c=colors, s=12)
plt.axvline(0, ls='--', c='black')
plt.xlabel('log₂ Fold Change (Modulated vs Contractile)')
plt.ylabel('-log₁₀ Adjusted P-value')
for gene in ['KCNMA1', 'SLC38A11', 'IGFBP5']:
    d = df[df.names == gene]
    if len(d):
        plt.annotate(gene, (d.logfoldchanges.iloc[0],
                           -np.log10(max(d.pvals_adj.iloc[0], 1e-300))),
                     xytext=(5, 5), textcoords='offset points')
plt.tight_layout()
plt.savefig("modulated_contractile_volcano.pdf")

# Transfer Back to Initial AnnData
whole_adata = sc.read_h5ad("/COMMA/results/scrna/downstream_analysis/integrated/integrated_final.h5ad")
vascular_adata = sc.read_h5ad("/COMMA/results/scrna/downstream_analysis/subclustering/vascular/vascular_annotated_final.h5ad")
vascular_adata_small = sc.read_h5ad("/COMMA/results/scrna/downstream_analysis/subclustering/vascular/vascular_annotated_small.h5ad")
whole_adata.obs['vascular_cell_type'] = vascular_adata.obs['vascular_cell_type'].reindex(whole_adata.obs_names).astype(object).fillna(whole_adata.obs['broad_cell_type'].astype(object))
whole_adata.obs['vsmc_cell_type']      = adata.obs['vsmc_cell_type'].reindex(whole_adata.obs_names).astype(object).fillna(whole_adata.obs['vascular_cell_type']).fillna(whole_adata.obs['broad_cell_type'].astype(object))
vascular_adata.obs['vsmc_cell_type']  = adata.obs['vsmc_cell_type'].reindex(vascular_adata.obs_names).astype(object).fillna(vascular_adata.obs['vascular_cell_type'].astype(object))
vascular_adata_small.obs['vsmc_cell_type'] = adata.obs['vsmc_cell_type'].reindex(vascular_adata_small.obs_names).astype(object).fillna(vascular_adata_small.obs['vascular_cell_type'].astype(object))

# Write Updated Annotation Files to File & Filter Leblanc 2025 retrospectively from Whole Cardiac Atlas
whole_adata = whole_adata[whole_adata.obs["Study"] != "Leblanc 2025"].copy()
whole_adata.write_h5ad("/COMMA/results/scrna/downstream_analysis/integrated/integrated_final.h5ad")
adata.uns['rank_genes_groups_filtered'] = None
adata.uns['log1p'] = None
adata.write_h5ad("vsmc_annotated_final.h5ad")