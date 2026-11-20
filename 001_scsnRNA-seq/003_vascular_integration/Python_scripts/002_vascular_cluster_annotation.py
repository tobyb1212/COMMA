import pandas as pd
import numpy as np
import scanpy as sc
import os
import argparse
import random
import matplotlib.pyplot as plt
import gc

random.seed(1234)

os.chdir("/COMMA/results/scrna/downstream_analysis/subclustering/vascular/")
 
# Selected scANVI with Resolution =0.8
adata = sc.read_h5ad(f"integration/scANVI_vascular_integrated.h5ad")
sc.pp.neighbors(adata, use_rep=f"X_scANVI")
sc.tl.umap(adata, min_dist = 0.2)
sc.tl.leiden(adata, resolution=0.8, key_added=f'leiden_scANVI_0.8')

# Remove Clusters with Cardiomyocyte Markers 
cm_clusters = ['12', '10']
adata_filtered = adata[~adata.obs["leiden_scANVI_0.8"].isin(cm_clusters)].copy()

# Annotate Vascular Cells on Basis of Marker Genes
vascular_annotations = {
    '0': 'Pericyte 1',
    '1': 'Capillary Endothelial Cell 1',
    '2': 'Capillary Endothelial Cell 2',
    '3': 'Arterial Endothelial Cell 1',
    '4': 'Vascular Smooth Muscle Cell',
    '5': 'Capillary Endothelial Cell 3',
    '6': 'Venous Endothelial Cell',
    '7': 'Pericyte 2', 
    '8': 'Activated/Angiogenic Endothelial Cell',
    '9': 'Pericyte 3',
    '11': 'Atria-enriched Vascular Smooth Muscle Cell',
    '13': 'Arterial Endothelial Cell 2',
    '14': 'Pericyte 4'}

# Remove Leblanc 2025 Data as this was not included in scATAC-seq Meta-Analysis
adata_filtered = adata_filtered[adata_filtered.obs["Study"] != "Leblanc 2025"].copy()

adata_filtered.obs['vascular_cell_type'] = adata_filtered.obs['leiden_scANVI_0.8'].map(vascular_annotations)
adata_filtered.obs['vascular_cell_type'] = adata_filtered.obs['vascular_cell_type'].astype('category')
adata_filtered.obs['vascular_cell_type'].value_counts()

# Figure 1C - annotated RNA microvascular atlas UMAP
sc.pl.umap(adata_filtered, color="vascular_cell_type", show=False, title=None)
plt.savefig("vascular_cell_type_umap.pdf", bbox_inches="tight")
plt.close()

order = ["Pericyte 1", "Pericyte 2", "Pericyte 3", "Pericyte 4", "Atria-enriched Vascular Smooth Muscle Cell", "Vascular Smooth Muscle Cell", "Capillary Endothelial Cell 1", "Capillary Endothelial Cell 2", "Capillary Endothelial Cell 3", "Activated/Angiogenic Endothelial Cell", "Arterial Endothelial Cell 1", "Arterial Endothelial Cell 2", "Venous Endothelial Cell"]
adata_filtered.obs["vascular_cell_type"] = pd.Categorical(adata_filtered.obs["vascular_cell_type"], categories=order, ordered=True)

sc.tl.rank_genes_groups(adata_filtered, groupby="vascular_cell_type", method = "wilcoxon")
sc.tl.filter_rank_genes_groups(adata_filtered)
# Supplemental Table 4
df_markers = sc.get.rank_genes_groups_df(adata_filtered, group=None, key="rank_genes_groups_filtered")
df_markers = df_markers.rename(columns={"pvals_adj": "FDR"})
df_markers["FDR"] = df_markers["FDR"].replace(0, "< 1x 10 ^ -308")
df_markers.to_csv(f"vascular_markers.csv", index=False)

names = adata_filtered.uns["rank_genes_groups_filtered"]["names"]
gene_groups = {}
for g in order:
    for gene in names[g]:
        if pd.notna(gene):
            gene_groups.setdefault(gene, set()).add(g)

# Select top 2 + top 2 unique genes
markers = {}
for g in order:
    ranked = [gene for gene in names[g] if pd.notna(gene)]
    # Top 2 genes overall
    top2 = ranked[:2]
    # Highest-ranked genes found only in this cell type
    unique_genes = [
        gene for gene in ranked
        if len(gene_groups[gene]) == 1
        and gene not in top2
    ]
    # Top 2 unique genes
    top2_unique = unique_genes[:2]
    markers[g] = top2 + top2_unique
# Figure 1F
sc.pl.dotplot(adata_filtered, markers, groupby="vascular_cell_type", dendrogram=False, categories_order=order, show=False, save="vascular_dotplot.pdf")

# Vascular Proportions
# Figure 1D
ct = pd.crosstab(adata_filtered.obs['Region'], adata_filtered.obs["vascular_cell_type"], normalize="index")
ax = ct.plot(kind="bar", stacked=True, figsize=(10, 5))
ax.set_ylabel("Proportion within region")
ax.set_xlabel("Region")
ax.legend(title=f"Region", bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.)
plt.tight_layout()
plt.savefig("vascular_cell_type_cluster_proportions.png")
plt.close()

# Figure 1E - CHRM3 in mural cells
adata_filtered_mural = adata_filtered[adata_filtered.obs["vascular_cell_type"].isin(["Pericyte 1", "Pericyte 2", "Pericyte 3", "Pericyte 4", "Vascular Smooth Muscle Cell", "Atria-enriched Vascular Smooth Muscle Cell"])].copy()
adata_filtered_mural.obs["vascular_cell_type"] = adata_filtered_mural.obs["vascular_cell_type"].cat.remove_unused_categories()
old_order = adata_filtered_mural.obs["vascular_cell_type"].cat.categories
palette = adata_filtered.uns.get("vascular_cell_type_colors")
color_dict = dict(zip(adata_filtered.obs["vascular_cell_type"].cat.categories, palette)) if palette is not None else None
order = ["Pericyte 1", "Pericyte 2", "Pericyte 3", "Pericyte 4", "Vascular Smooth Muscle Cell", "Atria-enriched Vascular Smooth Muscle Cell"]
adata_filtered_mural.obs["vascular_cell_type"] = adata_filtered_mural.obs["vascular_cell_type"].cat.reorder_categories(order, ordered=True)
if color_dict is not None:
    adata_filtered_mural.uns["vascular_cell_type_colors"] = [color_dict[x] for x in order]
ax = sc.pl.violin(adata_filtered_mural, ["CHRM3"], groupby="vascular_cell_type", stripplot=False, show=False, rotation=90, bw=0.3, cut=0, width=0.95)
ax.figure.set_size_inches(14, 3.2)
plt.tight_layout()
plt.savefig("CHRM3_violin_mural_vascular_cell_type.pdf")
plt.close()

# Supplemental Figure 5A-D
for gene in ['SULF1', 'DKK2', 'NEBL', 'BTNL9']:
    plt.figure(figsize=(6, 4))
    sc.pl.umap(adata_filtered, color=gene, show=False, title=None)
    plt.savefig(f"vascular_cell_type_{gene}_umap.png", dpi=900, bbox_inches="tight")
    plt.close()
    sc.pl.umap(adata_filtered[adata_filtered.obs.groupby("vascular_cell_type", group_keys=False).apply(lambda x: x.sample(1, random_state=0)).index], color=gene, save=f"vascular_cell_type_{gene}_umap_legend.pdf")

adata_filtered.write_h5ad("/COMMA/results/scrna/downstream_analysis/subclustering/vascular/vascular_annotated_final.h5ad")
adata_filtered_small = adata_filtered[:, 0].copy()
adata_filtered_small.write_h5ad("/COMMA/results/scrna/downstream_analysis/subclustering/vascular/vascular_annotated_small.h5ad")