import scanpy as sc
import pandas as pd
import os
import pandas as pd
import numpy as np
if not hasattr(pd.Series, "nonzero"):
    pd.Series.nonzero = lambda self: np.nonzero(self.to_numpy())
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from pathlib import Path

adata_atac = sc.read_h5ad("/COMMA/results/scatac/differential_analysis/genescores/vascular/vascular_genescore_matrix.h5ad")
adata_rna = sc.read_h5ad("/COMMA/results/scrna/downstream_analysis/subclustering/vascular/vascular_annotated_final.h5ad")
outdir = "/COMMA/results/scatac/differential_analysis/genescores/vascular/correlation"
os.chdir(outdir)

# Align Similar sc/snRNA-seq Cell Types
adata_rna.obs['vascular_cell_type_merged'] = adata_rna.obs['vascular_cell_type']
adata_rna.obs['vascular_cell_type_merged'] = adata_rna.obs['vascular_cell_type_merged'].replace({'Pericyte 1': 'Pericyte', 'Pericyte 2': 'Pericyte', 'Pericyte 3': 'Pericyte', 'Pericyte 4': 'Pericyte'})
adata_rna.obs['vascular_cell_type_merged'] = adata_rna.obs['vascular_cell_type_merged'].replace({'Arterial Endothelial Cell 1': 'Arterial Endothelial Cell', 'Arterial Endothelial Cell 2': 'Arterial Endothelial Cell'})
adata_rna.obs['vascular_cell_type_merged'] = adata_rna.obs['vascular_cell_type_merged'].replace({'Capillary Endothelial Cell 1': 'Capillary Endothelial Cell', 'Capillary Endothelial Cell 2': 'Capillary Endothelial Cell', 'Capillary Endothelial Cell 3': 'Capillary Endothelial Cell', 'Activated/Angiogenic Endothelial Cell': 'Capillary Endothelial Cell'})
adata_rna.obs['vascular_cell_type_merged'] = adata_rna.obs['vascular_cell_type_merged'].replace({'Vascular Smooth Muscle Cell': 'Vascular Smooth Muscle Cell', 'Atria-enriched Vascular Smooth Muscle Cell': 'Vascular Smooth Muscle Cell'})
adata_atac.obs['vascular_cell_type'] = adata_atac.obs['vascular_cell_type'].replace({'Vascular Smooth Muscle Cell': 'Vascular Smooth Muscle Cell'})

# Differential Gene Expression for Each Condition
sc.tl.rank_genes_groups(adata_rna, "vascular_cell_type_merged", method='wilcoxon')
sc.tl.rank_genes_groups(adata_atac, "vascular_cell_type", method='wilcoxon')
# Export to Dataframes
marker_genes_rna = sc.get.rank_genes_groups_df(adata_rna, group=None, key='rank_genes_groups')
marker_genes_atac = sc.get.rank_genes_groups_df(adata_atac, group=None, key='rank_genes_groups')

# Keep Common Genes from Dataframes
marker_genes_rna_filtered = marker_genes_rna[marker_genes_rna['names'].isin(marker_genes_atac['names'])].copy()
marker_genes_atac_filtered = marker_genes_atac[marker_genes_atac['names'].isin(marker_genes_rna['names'])].copy()

# Correlation of Marker Genes Per Cell-Type
results = []
for celltype in adata_rna.obs['vascular_cell_type_merged'].unique():
    celltype_safe = celltype.replace(" ", "_")
    marker_genes_rna_celltype = marker_genes_rna_filtered[marker_genes_rna_filtered['group'] == celltype].copy()
    marker_genes_atac_celltype = marker_genes_atac_filtered[marker_genes_atac_filtered['group'] == celltype].copy()
    rna_df = marker_genes_rna_celltype[["names", "scores"]].rename(columns={"scores": "scores_rna", "names": "gene"})
    atac_df = marker_genes_atac_celltype[["names", "scores"]].rename(columns={"scores": "scores_atac", "names": "gene"})
    merged = rna_df.merge(atac_df, on="gene", how="inner").dropna()
    if len(merged) >= 3:
        r, p = pearsonr(merged["scores_rna"], merged["scores_atac"])
        results.append({"celltype": celltype, "subset": "all_overlapping", "n_genes": len(merged), "r": r, "p": p})
        plt.figure(figsize=(6, 5))
        plt.scatter(merged["scores_rna"], merged["scores_atac"], alpha=0.7)
        plt.xlabel("RNA scores")
        plt.ylabel("ATAC scores")
        m, b = np.polyfit(merged["scores_rna"], merged["scores_atac"], 1)
        fit_x = np.linspace(merged["scores_rna"].min(), merged["scores_rna"].max(), 100)
        fit_y = m * fit_x + b
        plt.plot(fit_x, fit_y, linewidth=2, color='red')
        plt.title(f"{celltype}\nPearson r = {r:.3f}, p = {p:.2e}, n = {len(merged)}")
        plt.tight_layout()
        plt.savefig(f"{outdir}/{celltype_safe}_rna_atac_corr_all.pdf", dpi=300, bbox_inches="tight")
        plt.close()
    # Significant markers only, if p_adj exists in both tables
    sig_rna = marker_genes_rna_celltype[marker_genes_rna_celltype["pvals_adj"] < 0.05][["names", "scores"]]
    sig_atac = marker_genes_atac_celltype[marker_genes_atac_celltype["pvals_adj"] < 0.05][["names", "scores"]]
    sig_merged = sig_rna.rename(columns={"scores": "scores_rna", "names": "gene"}).merge(sig_atac.rename(columns={"scores": "scores_atac", "names": "gene"}), on="gene", how="inner").dropna()
    # Supplementary Figure 11
    if len(sig_merged) >= 3:
        r_sig, p_sig = pearsonr(sig_merged["scores_rna"], sig_merged["scores_atac"])
        results.append({"celltype": celltype, "subset": "significant_only", "n_genes": len(sig_merged), "r": r_sig, "p": p_sig})
        plt.figure(figsize=(6, 5))
        plt.scatter(sig_merged["scores_rna"], sig_merged["scores_atac"], alpha=0.7)
        plt.xlabel("RNA scores")
        plt.ylabel("ATAC scores")
        m, b = np.polyfit(sig_merged["scores_rna"], sig_merged["scores_atac"], 1)
        fit_x = np.linspace(sig_merged["scores_rna"].min(), sig_merged["scores_rna"].max(), 100)
        fit_y = m * fit_x + b
        plt.plot(fit_x, fit_y, linewidth=2, color='red')
        plt.title(f"{celltype} (pvals_adj < 0.05)\nPearson r = {r_sig:.3f}, p = {p_sig:.2e}, n = {len(sig_merged)}")
        plt.tight_layout()
        plt.savefig(f"{outdir}/{celltype_safe}_rna_atac_corr_sig.pdf", dpi=300, bbox_inches="tight")
        plt.close()

results_df = pd.DataFrame(results)
results_df.to_csv(f"{outdir}/correlation_results.csv", index=False)

# Figure 5E
sc.tl.filter_rank_genes_groups(adata_atac, min_fold_change=1)
sc.pl.rank_genes_groups_dotplot(adata_atac, key="rank_genes_groups_filtered", n_genes=3, save=f"atac_vascular_cell_type_dotplot.pdf", dendrogram = False)
