import pandas as pd
import numpy as np
import scanpy as sc
import os
import argparse
import random
import re
import scvi
import matplotlib.pyplot as plt
import gc
import matplotlib as mpl

random.seed(1234)

os.chdir("/COMMA/results/scrna/downstream_analysis/epicardial_coronary/analysis")

# Load Final Integration Model
adata_vascular = sc.read_h5ad("/COMMA/results/scrna/downstream_analysis/subclustering/vascular/vascular_annotated_final.h5ad")
adata = sc.read_h5ad("/COMMA/results/scrna/downstream_analysis/epicardial_coronary/integrated/preprocessed.h5ad")
adata.obs['Study'] = "Zhao 2025"
adata.obs['sample_type'] = "Cells"
adata.obs['broad_cell_type'] = "Unknown"
adata.obs['vascular_cell_type'] = "Unknown"
common_vascular_vars = adata_vascular.var_names.intersection(adata.var_names)
common_vascular_obs_cols = adata_vascular.obs.columns.intersection(adata.obs.columns)
common_vascular_var_cols = adata_vascular.var.columns.intersection(adata.var.columns)
a1 = adata_vascular[:, common_vascular_vars].copy()
a2 = adata[:, common_vascular_vars].copy()
a1.obs = a1.obs[common_vascular_obs_cols].copy(); a2.obs = a2.obs[common_vascular_obs_cols].copy()
a1.var = a1.var[common_vascular_var_cols].copy(); a2.var = a2.var[common_vascular_var_cols].copy()
adata_vascular = sc.concat([a1, a2], label="batch", keys=["atlas","epicardial"])
sc.pp.highly_variable_genes(adata_vascular, flavor="cell_ranger", n_top_genes=2000, layer="counts", batch_key="Study")
adata_vascular_hvg = adata_vascular[:, adata_vascular.var.highly_variable].copy()
lvae = scvi.model.SCANVI.load("/COMMA/results/scrna/downstream_analysis/epicardial_coronary/integrated/vascular_scanvi_model", adata = adata_vascular_hvg)
adata_vascular_hvg.obsm["X_scANVI"] = lvae.get_latent_representation()
adata_vascular.obsm['X_scANVI'] = adata_vascular_hvg.obsm["X_scANVI"]
probs = lvae.predict(soft=True)
adata_vascular.obs["scanvi_confidence"] = probs.max(axis=1)
adata_vascular.obs["scanvi_pred_label"] = probs.idxmax(axis=1)

# Filter Non-Vascular Clusters
sc.pp.neighbors(adata_vascular, use_rep="X_scANVI")
sc.tl.leiden(adata_vascular, resolution=0.3, key_added='leiden_scANVI_0.3')
sc.tl.umap(adata_vascular, min_dist = 0.2)
adata_vascular.write_h5ad("/COMMA/results/scrna/downstream_analysis/epicardial_all_cells.h5ad")
adata_vascular = sc.read_h5ad("/COMMA/results/scrna/downstream_analysis/epicardial_all_cells.h5ad")
sc.pl.umap(adata_vascular, color=["leiden_scANVI_0.3", "Study", "vascular_cell_type", "scanvi_pred_label", "scanvi_confidence", "scanvi_final_label"], ncols=3,save=f"_0.3_cluster.png")
mask = ~adata_vascular.obs['leiden_scANVI_0.3'].isin(['9', '10'])
adata_filtered = adata_vascular[mask].copy()
adata_filtered.write_h5ad("/COMMA/results/scrna/downstream_analysis/epicardial_coronary/analysis/vascular_only_cells.h5ad")
adata_filtered = sc.read_h5ad("/COMMA/results/scrna/downstream_analysis/epicardial_coronary/analysis/vascular_only_cells.h5ad")

# Preprocess Data
adata_filtered.X = adata_filtered.layers['counts'].copy()
sc.pp.normalize_total(adata_filtered, target_sum=1e4)
sc.pp.log1p(adata_filtered)
sc.pp.neighbors(adata_filtered, use_rep="X_scANVI")
sc.tl.umap(adata_filtered, min_dist = 0.2)
sc.pl.umap(adata_filtered, color = ["scanvi_pred_label", "scanvi_confidence"], save = "scanvi_predictions.png")

adata_filtered.obs["final_vascular_cell_type"] = adata_filtered.obs["vascular_cell_type"].astype(str).where(adata_filtered.obs["vascular_cell_type"] != "Unknown", adata_filtered.obs["scanvi_pred_label"].astype(str))
sc.pl.umap(adata_filtered, color = ["final_vascular_cell_type"], save = "final_vascular_cell_type.png")

# Perform DGE of Vascular Cell Types
sc.tl.rank_genes_groups(adata_filtered, "final_vascular_cell_type", method='wilcoxon')
sc.tl.filter_rank_genes_groups(adata_filtered, min_fold_change=1)

sc.pl.rank_genes_groups_dotplot(adata_filtered, key="rank_genes_groups_filtered", n_genes=5, save="scANVI_predoublet_dotplot.png", dendrogram = False)
df_markers = sc.get.rank_genes_groups_df(adata_filtered, group=None, key="rank_genes_groups_filtered")
df_markers.to_csv(f"markers.csv", index=False)

# Plot UMAP
order = ["Pericyte 1", "Capillary Endothelial Cell 1", "Capillary Endothelial Cell 2", "Arterial Endothelial Cell 1", "Vascular Smooth Muscle Cell", "Capillary Endothelial Cell 3", "Venous Endothelial Cell", "Pericyte 2", "Activated/Angiogenic Endothelial Cell", "Pericyte 3", "Atria-enriched Vascular Smooth Muscle Cell", "Arterial Endothelial Cell 2", "Pericyte 4"]
adata_filtered.obs["final_vascular_cell_type"] = pd.Categorical(adata_filtered.obs["final_vascular_cell_type"], categories=order, ordered=True)
rna_cats = list(adata_filtered.obs["final_vascular_cell_type"].cat.categories)
cmap_rna = mpl.cm.get_cmap("tab20", len(rna_cats))
rna_colors = [mpl.colors.to_hex(cmap_rna(i)) for i in range(len(rna_cats))]
adata_filtered.uns["final_vascular_cell_type_colors"] = rna_colors
# Figure 3B
plt.figure(figsize=(6, 4))
sc.pl.umap(adata_filtered, color="final_vascular_cell_type", show=False, title=None)
plt.savefig("vascular_umap.png", dpi=900, bbox_inches="tight")
plt.close()
sc.pl.umap(adata_filtered[adata_filtered.obs.groupby("final_vascular_cell_type", group_keys=False).apply(lambda x: x.sample(1, random_state=0)).index], color="final_vascular_cell_type", save="vascular_labels.pdf")

# Plot UMAPs by Origin
adata_filtered.obs["batch"] = pd.Categorical(adata_filtered.obs["batch"])
adata_filtered.uns["batch_colors"] = ["#2A9D8F", "#E9A23B"]
# Figure 3A
plt.figure(figsize=(6, 4))
sc.pl.umap(adata_filtered, color="batch", show=False, title=None)
plt.savefig("batch_umap.png", dpi=900, bbox_inches="tight")
plt.close()
sc.pl.umap(adata_filtered[adata_filtered.obs.groupby("batch", observed=True, group_keys=False).apply(lambda x: x.sample(1, random_state=0)).index],
    color="batch", save="batch_labels.pdf")

adata_filtered_small = adata_filtered[:, 0].copy()
adata_filtered_small.write_h5ad("epicardial_integrated_final_small.h5ad")
adata_filtered.write_h5ad("epicardial_integrated_final.h5ad")
