import pandas as pd
import numpy as np
import scanpy as sc
import os
import argparse
import random
import celltypist
from celltypist import models
import re
import matplotlib.pyplot as plt
import matplotlib as mpl

os.chdir("/COMMA/results/scrna/downstream_analysis/integrated")

# Read in Original Integrated Adata Object
adata = sc.read_h5ad("/COMMA/results/scrna/integrated/scANVI/scANVI_integrated.h5ad")

# Filter Higher Resolution Clusters where 30% of Scrublet Scores >0.2
sc.tl.leiden(adata, resolution=0.5, key_added='leiden_scANVI_0.5')
adata.obs['doublet_cluster'] = adata.obs['leiden_scANVI_0.5'].map(adata.obs.groupby('leiden_scANVI_0.5')["scrublet_score"]
      .apply(lambda scores: (scores > .2).sum() / len(scores) > 0.3))

# Plot UMAP of Clusters, Doublet Score and Doublet Clusters
sc.pl.umap(adata, color=["doublet_cluster"], ncols=2, save="scANVI_doublet_umaps.png")

# Perform DGE of Leiden Resolution = 0.5
sc.tl.rank_genes_groups(adata, "leiden_scANVI_0.5", method='wilcoxon')
sc.tl.filter_rank_genes_groups(adata, min_fold_change=1)
sc.pl.rank_genes_groups_dotplot(adata, key="rank_genes_groups_filtered", n_genes=5, save="scANVI_predoublet_dotplot.png", dendrogram = False)

# Filter Selected Doublet Clusters 
adata_filtered = adata[adata.obs['doublet_cluster'] == False].copy()

# Recompute UMAP
sc.tl.umap(adata_filtered, min_dist=0.2)

# Read in pre-existing broad HCA annotations to train new CellTypist Model 
# New model is being trained as extra cell filtered as doublet clusters in the final integration
adata_hca = adata_filtered[adata_filtered.obs["Study"].isin(["Litvinukova 2020", "Kanemaru 2023"]), :]
df = pd.read_csv("/COMMA/results/scrna/metadata/HCA_celltypes.csv")
celltype_dict = dict(zip(df['cell_name'], df['cell_type']))
adata_hca.obs['HCA_cell_type'] = adata_hca.obs_names.map(celltype_dict)
adata_hca.obs['HCA_cell_type'] = adata_hca.obs['HCA_cell_type'].fillna('Unknown')
adata_hca = adata_hca[adata_hca.obs['HCA_cell_type'] != 'Unknown', :]

# Converting Separate Cardiomyocyte Labels into Single Label
adata_hca.obs['HCA_cell_type'] = adata_hca.obs['HCA_cell_type'].replace(
    {'Ventricular Cardiomyocyte': 'Cardiomyocyte',
    'Atrial Cardiomyocyte': 'Cardiomyocyte'})

# Train & save this new cell typist model
model = celltypist.train(adata_hca, labels = 'HCA_cell_type', n_jobs = 10, feature_selection = True)
model.write('HCA_celltypes_final.pkl')

# Predict cell types based on celltypist model from HCA data
predictions = celltypist.annotate(adata_filtered, model = 'HCA_celltypes_final.pkl', majority_voting = False)
adata_filtered = predictions.to_adata(prefix='HCA_celltypist_')

# Majority Cell Type Voting within Clusters
sc.tl.leiden(adata_filtered, resolution=0.3, key_added='leiden_scANVI_0.3')
cluster_labels = (adata_filtered.obs.groupby('leiden_scANVI_0.3')["HCA_celltypist_predicted_labels"].agg(lambda x: x.value_counts().idxmax()))
adata_filtered.obs["annotated_leiden_scANVI_0.3"] = adata_filtered.obs['leiden_scANVI_0.3'].map(cluster_labels)

# Identify clusters heavily biased by one study
study_props = (adata_filtered.obs.groupby(['leiden_scANVI_0.3', 'Study']).size().reset_index(name='count'))
study_props['prop'] = (study_props.groupby('leiden_scANVI_0.3')['count'].transform(lambda x: x / x.sum()))
biased_clusters = (study_props[study_props['prop'] >= 0.5].sort_values(['leiden_scANVI_0.3', 'prop'], ascending=[True, False]))
adata_filtered.obs['biased_cluster'] = adata_filtered.obs['leiden_scANVI_0.3'].isin(biased_clusters['leiden_scANVI_0.3'])
sc.pl.umap(adata_filtered, color=["biased_cluster"], save = "scANVI_biased_cluster.png")

# Perform Manual Cluster Annotation
manual_annotations = {
    '0': 'Cardiomyocytes',
    '1': 'Fibroblasts',
    '2': 'Cardiomyocytes',
    '3': 'Endothelial Cells',
    '4': 'Mural Cells',
    '5': 'Cardiomyocytes',
    '6': 'Myeloid Cells',
    '7': 'Endocardial Endothelial Cells',
    '8': 'Lymphoid Cells',
    '9': 'Mesothelial Cells',
    '10': 'Adipocytes',
    '11': 'Glial Cells',
    '12': 'Lymphatic Endothelial Cells',
    '13': 'Neural Cells'
}

adata_filtered.obs['broad_cell_type'] = adata_filtered.obs['leiden_scANVI_0.3'].map(manual_annotations)
adata_filtered.obs['broad_cell_type'] = adata_filtered.obs['broad_cell_type'].astype('category')

# Write Final UMAP & Integrated and Filtered .h5ad File - Small File for Multiomic Integration
sc.pl.umap(adata_filtered, color=["annotated_leiden_scANVI_0.3"], save = "scANVI_celltypist_umap.png")
sc.pl.umap(adata_filtered, color=["broad_cell_type"], save = "scANVI_celltype_umap.png")

# Print Summary Statistics for Integrated Object
adata_filtered
print(adata_filtered.obs["Study"].value_counts())
print(adata_filtered.obs["annotated_leiden_scANVI_0.3"].value_counts())
print(adata_filtered.obs["broad_cell_type"].value_counts())
print(adata_filtered.obs["sample_type"].value_counts())

# DGE Between Clusters - Supplemental Table 3
sc.tl.rank_genes_groups(adata_filtered, groupby='broad_cell_type', method='wilcoxon')
sc.tl.filter_rank_genes_groups(adata_filtered)
# Supplemental Figure 4
sc.pl.rank_genes_groups_dotplot(adata_filtered, groupby='broad_cell_type', dendrogram=False, key="rank_genes_groups_filtered", show=False, n_genes=3, save=f"whole_rna_dotplot.pdf")
df_markers = sc.get.rank_genes_groups_df(adata_filtered, group=None, key="rank_genes_groups_filtered")
df_markers = df_markers.rename(columns={"group":"Group", "names":"Gene", "logfoldchanges":"Log2(Fold Change in Expression)", "pvals_adj":"FDR"})[["Group","Gene","Log2(Fold Change in Expression)","FDR"]].dropna(subset=["Gene"])
df_markers = df_markers[df_markers["Gene"].astype(str).str.strip().ne("")]
df_markers["FDR"] = df_markers["FDR"].replace(0, "< 1x 10 ^ -308")
df_markers.to_csv(f"/COMMA/results/scrna/downstream_analysis/integrated/broad_cell_type_markers.csv", index=False)

# Figure 1B
col = "broad_cell_type"
adata_filtered.obs[col] = pd.Categorical(adata_filtered.obs[col])
adata_filtered_cats = list(adata_filtered.obs[col].cat.categories)
cmap_rna = mpl.cm.get_cmap("tab20", max(1, len(adata_filtered_cats)))
adata_filtered_colors = [mpl.colors.to_hex(cmap_rna(i)) for i in range(len(adata_filtered_cats))]
adata_filtered.uns[f"{col}_colors"] = adata_filtered_colors
plt.figure(figsize=(6, 4))
sc.pl.umap(adata_filtered, color=col, show=False, title=None)
plt.savefig("whole_rna_umap.png", dpi=900, bbox_inches="tight")
plt.close()
sc.pl.umap(adata_filtered[adata_filtered.obs.groupby(col, group_keys=False).apply(lambda x: x.sample(1, random_state=0)).index], color=col, save="whole_rna_labels.pdf")

adata_filtered.write_h5ad("integrated_final.h5ad")
# Leblanc et al 2025 subsequently removed from integrated_final.h5ad
