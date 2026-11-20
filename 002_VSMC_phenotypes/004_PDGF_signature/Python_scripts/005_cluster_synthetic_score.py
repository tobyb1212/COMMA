import scanpy as sc
import pandas as pd
import os

adata = sc.read_h5ad("/COMMA/results/scrna/downstream_analysis/subclustering/vsmc/vsmc_annotated_final.h5ad")

outdir = "/COMMA/results/scrna/downstream_analysis/subclustering/vsmc/synthetic_score"
os.chdir(outdir)

# Read in Synthetic Genes from Bickel et al 2024. PDGF-BB Treated microvascular SMCs
synthetic_genes = pd.read_csv("/COMMA/results/scrna/downstream_analysis/subclustering/vsmc/synthetic_score/synthetic_genes.txt", header=0).iloc[:, 0].dropna().tolist()
contractile_genes = pd.read_csv("/COMMA/results/scrna/downstream_analysis/subclustering/vsmc/synthetic_score/contractile_genes.txt", header=0).iloc[:, 0].dropna().tolist()

# Score Individual Cells Based on Gene Lists
sc.tl.score_genes(adata, gene_list=synthetic_genes, ctrl_size=len(synthetic_genes), score_name='synthetic_score', random_state = 1234)
sc.tl.score_genes(adata, gene_list=contractile_genes, ctrl_size=len(contractile_genes), score_name='contractile_score', random_state = 1234)

# Visualise Synthetic Score Predictions
sc.pl.umap(adata, color=["synthetic_score", "vsmc_cell_type"], save = "synthetic_cluster.png", cmap = "Reds")
avg_scores = (adata.obs.groupby("vsmc_cell_type")["consensus_score"].mean().sort_values(ascending=False))
print(avg_scores.head(10))

# Visualise contractile Score Predictions
sc.pl.umap(adata, color=["contractile_score", "vsmc_cell_type"], save = "contractile_cluster.png", cmap = "Blues")
avg_scores = (adata.obs.groupby("vsmc_cell_type")["contractile_score"].mean().sort_values(ascending=False))
print(avg_scores.head(10))

# Violin of Synthetic Score 
#Figure 2C
sc.pl.violin(adata, keys=["synthetic_score", "contractile_score"], stripplot = False, groupby="vsmc_cell_type", save="_synthetic_violin.pdf")
