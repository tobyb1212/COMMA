import numpy as np
import scanpy as sc
import pylemur
import os
import matplotlib.pyplot as plt
from scipy import sparse
import pandas as pd

os.chdir("/COMMA/results/scrna/downstream_analysis/DGE/lemur")

adata = sc.read_h5ad("/COMMA/results/scrna/downstream_analysis/subclustering/vsmc/vsmc_annotated_final.h5ad")

# Figure 4B
genes = pd.read_csv("/COMMA/results/scrna/downstream_analysis/DGE/consensus/outputs/vsmc/whole_vsmc_Cardiac_Disease_filtered.csv")
genes = genes["Gene"]
sc.pl.dotplot(adata, var_names=genes, groupby="Cardiac Disease", save="dge_dotplot.pdf")

# Supplemental Figure 9B
sc.pl.dotplot(adata, var_names=["TMEM51", "ADAMTS9"], groupby="Primary Cardiac Disease", save="tmem51_adamts9_dotplot.pdf")

# Figure 4C
for gene in ["TMEM51", "ADAMTS9"]:
    for status, label in [("No", "Control"), ("Yes", "CVD")]:
        selected = adata[adata.obs["Cardiac Disease"].astype(str) == status].copy()
        if selected.n_obs == 0:
            raise ValueError(f"No cells found for Figure 4C: {label}")
        sc.pl.umap(selected, color=gene, show=False, title=f"{gene} ({label})")
        plt.savefig(f"Figure_4C_{gene}_{label}.pdf", bbox_inches="tight")
        plt.close()
sc.pl.umap(adata, color="vsmc_cell_type", show=False, title=None)
plt.savefig("Figure_4C_VSMC_phenotypes.pdf", bbox_inches="tight")
plt.close()

# Preprocess Object
adata.X = adata.layers['counts'].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="cell_ranger", n_top_genes=2000, batch_key="Study")
sc.pl.highly_variable_genes(adata, save = "highlyvariablegenes.png")
adata_hvg = adata[:, adata.var.highly_variable]
adata_hvg = adata_hvg.copy()
sc.tl.pca(adata)
sc.pp.neighbors(adata)
adata.obs['Cardiac_Disease'] = adata.obs['Cardiac Disease']

# Train LEMUR Model
model = pylemur.tl.LEMUR(adata, design = "~ Study + `Sample Name` + sample_type + Cardiac_Disease", n_embedding=15)
model.fit()
print(model)

# Predict Control & Disease Differential Gene Expression
ctrl_pred = model.predict(new_condition=model.cond(Cardiac_Disease="No"))
disease_pred = model.predict(new_condition=model.cond(Cardiac_Disease="Yes"))
adata.layers["diff"] = disease_pred - ctrl_pred
adata.write_h5ad("lemur_vsmc_diff.h5ad")

# Investigate Specific Genes of Interest
# Supplemental Figure 9C
sel_genes = ["ADAMTS9", "TMEM51"]
for sel_gene in sel_genes:
    fsize = plt.rcParams['figure.figsize']
    fig = plt.figure(figsize=(fsize[0] * 4, fsize[1]))
    axs = [fig.add_subplot(1, 4, i+1) for i in range(4)]
    for ax in axs:
        ax.set_aspect('equal')
    sc.pl.umap(adata, layer="diff", color=[sel_gene], cmap = plt.get_cmap("seismic"), vcenter=0,
        vmin=-4, vmax=4, title="Predicted Differential Expression  (Disease / Control)", ax=axs[0], show=False)
    plt.savefig(f"{sel_gene}_dge.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{sel_gene}_dge.pdf", bbox_inches="tight")
    plt.close()
