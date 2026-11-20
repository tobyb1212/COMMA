import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np

rna = sc.read_h5ad("/COMMA/results/scrna/downstream_analysis/integrated/integrated_final.h5ad")
atac = sc.read_h5ad("/COMMA/results/scatac/whole_clustered.h5ad")
atac.obsm['X_umap'] = atac.obsm['umap']
atac.obsm["X_umap"] = atac.obsm["X_umap"].loc[atac.obs_names].to_numpy()
del atac.obsm["umap"]

# Figure 5B
col = "broad_cell_type"
rna.obs[col] = pd.Categorical(rna.obs[col])
atac.obs[col] = pd.Categorical(atac.obs[col])
rna_cats = list(rna.obs[col].cat.categories)
cmap_rna = mpl.cm.get_cmap("tab20", max(1, len(rna_cats)))
rna_colors = [mpl.colors.to_hex(cmap_rna(i)) for i in range(len(rna_cats))]
rna.uns[f"{col}_colors"] = rna_colors
plt.figure(figsize=(6, 4))
sc.pl.umap(rna, color=col, show=False, title=None)
plt.savefig("whole_rna_umap.png", dpi=900, bbox_inches="tight")
plt.close()
sc.pl.umap(rna[rna.obs.groupby(col, group_keys=False).apply(lambda x: x.sample(1, random_state=0)).index], color=col, save="whole_rna_labels.pdf")
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
plt.figure(figsize=(6, 4))
sc.pl.umap(atac, color=col, show=False, title=None)
plt.savefig("whole_atac_umap.png", dpi=900, bbox_inches="tight")
plt.close()
sc.pl.umap(atac[atac.obs.groupby(col, group_keys=False).apply(lambda x: x.sample(1, random_state=0)).index], color=col, save="whole_atac_labels.pdf")
