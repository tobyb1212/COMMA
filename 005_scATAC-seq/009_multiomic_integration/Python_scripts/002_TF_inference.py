import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
import networkx as nx
import os
from scipy import sparse
from ctxcore.genesig import GeneSignature
from pyscenic.aucell import aucell
import matplotlib as mpl

os.chdir("/COMMA/results/scatac/multiomic_integration/vsmc")
os.makedirs("aucell_analysis", exist_ok=True)
mpl.rcParams.update({
    "font.family": "Arial",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

# Read in RNA & Regulatory Network
adata_rna = sc.read_h5ad("rna_glue.h5ad")
grn = nx.read_graphml("vsmc_tf_network.graphml")

# Extract TF Regulons
regulons = []
for tf in grn.nodes:
    if grn.out_degree(tf) == 0:
        continue
    targets = list(set(grn.successors(tf)).intersection(adata_rna.var_names) - {tf})
    if len(targets) >= 3:
        regulons.append(GeneSignature(name=str(tf), gene2weight={gene: 1.0 for gene in targets}))

print("Number of TF regulons:", len(regulons))

# Calculate AUCell Activity
rna_expr = adata_rna.layers["norm"]
if sparse.issparse(rna_expr):
    rna_expr = rna_expr.toarray()
expr = pd.DataFrame(rna_expr, index=adata_rna.obs_names, columns=adata_rna.var_names)
auc = aucell(expr, regulons, num_workers=8)
auc.to_csv("aucell_analysis/regulon_activity.csv")

# Calculate Mean AUCell Activity Per VSMC Subtype
cluster_auc = auc.join(
    adata_rna.obs["predicted_vsmc_cell_type"]
).groupby("predicted_vsmc_cell_type", observed=True).mean()
cluster_auc.to_csv("aucell_analysis/mean_regulon_activity.csv")

# Z-score TF Activity Across VSMC Subtypes
auc_z = cluster_auc.T
auc_z = auc_z.sub(auc_z.mean(axis=1), axis=0).div(
    auc_z.std(axis=1).replace(0, np.nan), axis=0
).fillna(0)

# Select Top 5 TF Regulons Per VSMC Subtype
selected_tfs = []
for cluster in auc_z.columns:
    selected_tfs.extend(auc_z[cluster].nlargest(5).index.tolist())
selected_tfs = list(dict.fromkeys(selected_tfs))
heatmap = auc_z.loc[selected_tfs]

# Plot AUCell Heatmap
sns.set_style("white")
g = sns.clustermap(heatmap,cmap="RdBu_r",center=0,row_cluster=False,col_cluster=False,figsize=(12, max(8, len(selected_tfs) * 0.3)),xticklabels=True,yticklabels=True,cbar_kws={"label": "Mean AUCell activity (Z-score)"})

# Figure 5G
g.ax_heatmap.set_xlabel("VSMC subtype")
g.ax_heatmap.set_ylabel("TF regulon")
g.ax_heatmap.set_box_aspect(1)
plt.savefig("aucell_analysis/TF_regulon_heatmap.pdf", bbox_inches="tight")
plt.savefig("aucell_analysis/TF_regulon_heatmap.png", dpi=300, bbox_inches="tight")
plt.close()

print("AUCell analysis complete.")