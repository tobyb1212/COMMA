import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np

# Figure 5C
rna = sc.read_h5ad("/COMMA/results/scrna/downstream_analysis/subclustering/vascular/vascular_annotated_final.h5ad")
atac = sc.read_h5ad("/COMMA/results/scatac/vascular_clustered.h5ad")
atac.obsm['X_umap'] = atac.obsm['umap']
atac.obsm["X_umap"] = atac.obsm["X_umap"].loc[atac.obs_names].to_numpy()
del atac.obsm["umap"]

col = "vascular_cell_type"

atac.obs["vascular_cell_type"] = atac.obs["vascular_cell_type"].replace({
    "Pericyte": "Pericyte 1", "Arterial Endothelial Cell": "Arterial Endothelial Cell 1", "Capillary Endothelial Cell": "Capillary Endothelial Cell 2", "Vascular Smooth Muscle Cell": "Vascular Smooth Muscle Cell"
})

rna.obs[col] = pd.Categorical(rna.obs[col])
atac.obs[col] = pd.Categorical(atac.obs[col])
rna_cats = list(rna.obs[col].cat.categories)
cmap_rna = mpl.cm.get_cmap("tab20", max(1, len(rna_cats)))
rna_colors = [mpl.colors.to_hex(cmap_rna(i)) for i in range(len(rna_cats))]
rna.uns[f"{col}_colors"] = rna_colors
plt.figure(figsize=(6, 4))
sc.pl.umap(rna, color=col, show=False, title=None)
plt.savefig("vascular_rna_umap.png", dpi=900, bbox_inches="tight")
plt.close()
sc.pl.umap(rna[rna.obs.groupby(col, group_keys=False).apply(lambda x: x.sample(1, random_state=0)).index], color=col, save="vascular_rna_labels.pdf")
rna_color_map = dict(zip(rna_cats, rna_colors))
atac_cats = list(atac.obs[col].cat.categories)
shared = [c for c in atac_cats if c in rna_color_map]
non_shared = [c for c in atac_cats if c not in rna_color_map]
cmap_other = mpl.cm.get_cmap("Set3", max(1, len(non_shared)))
other_colors = [mpl.colors.to_hex(cmap_other(i)) for i in range(len(non_shared))]
atac_colors = []
other_it = iter(other_colors)
for cat in atac_cats:
    if cat in rna_color_map:
        atac_colors.append(rna_color_map[cat])
    else:
        atac_colors.append(next(other_it))


atac.uns[f"{col}_colors"] = atac_colors
plt.figure(figsize=(6,4))
sc.pl.umap(atac, color=col, show=False, title=None)
plt.savefig("vascular_atac_umap.png", dpi=900, bbox_inches="tight")
plt.close()
sc.pl.umap(atac[atac.obs.groupby(col, group_keys=False).apply(lambda x: x.sample(1, random_state=0)).index], color=col, save="vascular_atac_labels.pdf")

# Figure 5D
ct = pd.crosstab(atac.obs['Region'], atac.obs["vascular_cell_type"], normalize="index")
ax = ct.plot(kind="bar", stacked=True, figsize=(10, 5))
ax.set_ylabel("Proportion within region")
ax.set_xlabel("Region")
ax.legend(title=f"Region", bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.)
plt.tight_layout()
plt.savefig("vascular_cell_type_cluster_proportions.png")
plt.close()