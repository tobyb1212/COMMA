import matplotlib.pyplot as plt
import matplotlib as mpl
import anndata as ad
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
import torch
import glob
import os
import re
from scipy import sparse
import scglue
import networkx as nx
from itertools import chain

os.chdir("/COMMA/results/scatac/multiomic_integration/vsmc")

# Read in RNA & ATAC Files
adata_rna = sc.read_h5ad("/COMMA/results/scrna/downstream_analysis/subclustering/vsmc/vsmc_annotated_final.h5ad")
adata_rna.layers['norm'] = adata_rna.X
adata_rna.X = adata_rna.layers["counts"].copy()
adata_rna.var.modality = "RNA"
adata_rna.var['Modality'] = "RNA"
adatas_atac = {}
file_pattern = os.path.join("/COMMA/results/scatac/multiomic_integration/vsmc/peak_matrices", "**", "*_PeakMatrix.h5ad")
files = glob.glob(file_pattern, recursive=True)
for file_path in files:
    sample_adata = sc.read_h5ad(file_path)
    sample_id = os.path.basename(file_path).split("_PeakMatrix.h5ad")[0]
    adatas_atac[sample_id] = sample_adata
adata_atac = ad.concat(adatas_atac, label="fastq_name_new", merge="same")
adata_atac.obs_names_make_unique()
adata_atac.obs.rename(columns={'Sample.Name': 'Sample Name'}, inplace=True)
adata_atac.var.modality = "ATAC"
adata_atac.var['Modality'] = "ATAC"

# Rename ATAC Cells to Match RNA Barcode
archr_cells = adata_atac.obs_names.astype(str)
barcode_raw = [re.sub(r'.*#', '', s) for s in archr_cells]
barcode_clean = [re.sub(r'-1$', '', s) for s in barcode_raw]
linked = adata_atac.obs["linked_RNA"].astype(object)
cond = (linked.notna() &linked.astype(str).str.strip().ne("") &linked.astype(str).str.upper().ne("NA"))
rna_sample = linked.where(cond, adata_atac.obs["Sample Name"])
adata_atac.obs_names = (rna_sample.astype(str) + "_" + pd.Series(barcode_clean, index=adata_atac.obs.index)).tolist()
adata_atac.obs_names_make_unique()

# Trim Metadata to Clinical Data
rna_obs = adata_rna.obs.copy()
atac_obs = adata_atac.obs.copy()
common_cols = rna_obs.columns.intersection(atac_obs.columns)
rna_cols = common_cols.union(["sample_type"]).intersection(rna_obs.columns)
rna_obs = rna_obs[rna_cols].copy()
atac_obs = atac_obs[common_cols].copy()
adata_rna.obs = rna_obs
adata_atac.obs = atac_obs

# Filter to Nuclei Only
adata_rna = adata_rna[adata_rna.obs["sample_type"] == "Nuclei"].copy()
adata_atac = adata_atac[adata_atac.obs["sample_type"] == "Nuclei"].copy()
print(f"RNA nuclei: {adata_rna.n_obs}")
print(f"ATAC nuclei: {adata_atac.n_obs}")

# Select 2000 HVGs and 40000 Highly Variable Peaks
sc.pp.highly_variable_genes(adata_atac,n_top_genes=40000,flavor="seurat_v3", batch_key="Study")
sc.pp.normalize_total(adata_rna, target_sum=1e4)
sc.pp.highly_variable_genes(adata_rna, flavor="cell_ranger", n_top_genes=6000, batch_key="Study")

# GLUE Analysis
# Add Genomic Coordinates
scglue.data.get_gene_annotation(adata_rna, gtf="//COMMA/reference_genome/GENCODE/GRCh38.p14.GENCODE_release46/gencode.v46.basic.annotation_filtered.gtf", gtf_by="gene_name")
adata_atac.var["chrom"] = adata_atac.var['chr']
adata_atac.var["chromStart"] = adata_atac.var['start']
adata_atac.var["chromEnd"] = adata_atac.var['end']

# Preprocess RNA
adata_rna.X = adata_rna.layers['counts'].copy()
valid = adata_rna.var[["chrom", "chromStart", "chromEnd"]].notna().all(axis=1)
adata_rna = adata_rna[:, valid].copy()
adata_rna.var["chromStart"] = adata_rna.var["chromStart"].astype(int)
adata_rna.var["chromEnd"] = adata_rna.var["chromEnd"].astype(int)
sc.pp.normalize_total(adata_rna)
sc.pp.log1p(adata_rna)
sc.pp.scale(adata_rna, max_value=10)
sc.tl.pca(adata_rna, n_comps=100)

# Preprocess ATAC
scglue.data.lsi(adata_atac, n_components=100, n_iter=15)

# Construct Guidance Graph
guidance = scglue.genomics.rna_anchored_guidance_graph(adata_rna, adata_atac, promoter_len=2000, extend_range=150000, extend_fn=lambda x: np.exp(-x / 50000))
scglue.graph.check_graph(guidance, [adata_rna, adata_atac])
nx.write_graphml(guidance, "guidance.graphml.gz")

# Setup GLUE Model and Account for Batch Effects
scglue.models.configure_dataset(adata_rna, "NB", use_highly_variable=True, use_layer="counts", use_rep="X_pca", use_batch="Study", use_obs_names=True)
scglue.models.configure_dataset(adata_atac, "NB", use_highly_variable=True, use_rep="X_lsi", use_batch="Study", use_obs_names=True)

guidance_hvf = guidance.subgraph(chain(adata_rna.var.query("highly_variable").index,adata_atac.var.query("highly_variable").index)).copy()

model = scglue.models.fit_SCGLUE({"rna": adata_rna, "atac": adata_atac}, guidance_hvf, model=scglue.models.PairedSCGLUEModel)
model.save("glue_model.dill")

# Get Latent Representations
adata_rna.obsm["X_glue"] = model.encode_data("rna", adata_rna)
adata_atac.obsm["X_glue"] = model.encode_data("atac", adata_atac)
feature_embeddings = model.encode_graph(guidance_hvf)
feature_embeddings = pd.DataFrame(feature_embeddings, index=model.vertices)
adata_rna.varm["X_glue"] = feature_embeddings.reindex(adata_rna.var_names).to_numpy()
adata_atac.varm["X_glue"] = feature_embeddings.reindex(adata_atac.var_names).to_numpy()

# Combine RNA & ATAC Latent Spaces
scglue.data.transfer_labels(ref=adata_rna,query=adata_atac,field="vsmc_cell_type",use_rep="X_glue",n_neighbors=30,key_added="predicted_vsmc_cell_type")
adata_rna.obs['predicted_vsmc_cell_type'] = adata_rna.obs['vsmc_cell_type'].copy()
common = adata_rna.obs_names.intersection(adata_atac.obs_names)
rna_only, atac_only = adata_rna[~adata_rna.obs_names.isin(common)], adata_atac[~adata_atac.obs_names.isin(common)]
rna_only.obs["modality"], atac_only.obs["modality"], adata_rna.obs.loc[common, "modality"] = "expression", "accessibility", "paired"
obs = pd.concat([rna_only.obs, atac_only.obs, adata_rna[common].obs])
adata_glue = ad.AnnData(obs=obs)
adata_glue.obsm["X_glue"] = np.vstack([rna_only.obsm["X_glue"], atac_only.obsm["X_glue"], (adata_rna[common].obsm["X_glue"] + adata_atac[common].obsm["X_glue"]) / 2])

# Get UMAP
sc.pp.neighbors(adata_glue, use_rep="X_glue")
sc.tl.leiden(adata_glue, resolution=0.5, key_added="leiden_GLUE")
sc.tl.umap(adata_glue, min_dist = 1)
adata_glue.write_h5ad("combined_glue.h5ad")

# Plot UMAPs
sc.pl.umap(adata_glue, color=["leiden_GLUE", "Study"], show=False)
plt.savefig("umap_leiden_study.png")
plt.close()
sc.pl.umap(adata_glue, color=["vascular_cell_type", "modality"], show=False)
plt.savefig("umap_modality_celltype.png")
plt.close()

# Write RNA & ATAC Objects
adata_rna.write("rna_glue.h5ad")
adata_atac.write("atac_glue.h5ad")

# Plot Each Modality
study_colors = {
    "Amrute 2024": "#f57e20", "Hocker 2021": "#dbdb8c",
    "Kanemaru 2023": "#9372b2", "Leblanc 2025": "#9fdae5",
    "Read 2024": "#ae4a4b", "Selewa 2023": "#8c6d30",
    "Amrute 2023": "#3274A2", "Chaffin 2022": "#369067",
    "Hill 2024": "#C13D3E", "Koenig 2022": "#845C54",
    "Linna-Kuosmanen 2024": "#D484B8", "Litvinukova 2020": "#ABB26D",
    "Rao 2021": "#2CABB8", "Simonson 2023": "#B6C8E0",
    "Tucker 2020": "#EDBB88"}

mpl.rcParams.update({
    "font.family": "Arial",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

# Supplemental Figure 13
for mod in adata_glue.obs["modality"].unique():
    vsmc = adata_glue[adata_glue.obs["modality"] == mod].copy()
    name = f"glue_{mod}"
    for col in ["Study", "predicted_vsmc_cell_type"]:
        vsmc.obs[col] = vsmc.obs[col].astype("category")
        if col == "Study":
            vsmc.uns[f"{col}_colors"] = [study_colors[x] for x in vsmc.obs[col].cat.categories]
        elif col == "predicted_vsmc_cell_type":
            vsmc.uns[f"{col}_colors"] = list(map(mpl.colors.to_hex, plt.get_cmap("Dark2").colors[:len(vsmc.obs[col].cat.categories)]))
        sc.pl.umap(vsmc, color=col,show=False, title=None)
        plt.savefig(f"{name}_{col}_umap.png", dpi=900, bbox_inches="tight")
        plt.close()
        labels = vsmc.obs.groupby(col, observed=True).sample(1, random_state=0).index
        sc.pl.umap(vsmc[labels].copy(), color=col, size=5, show=False, title=None)
        plt.savefig(f"{name}_{col}_labels.pdf", bbox_inches="tight")
        plt.close()

# Predict Regulatory TFs
adata_rna.var["name"] = adata_rna.var_names
adata_atac.var["name"] = adata_atac.var_names
adata_rna.obs["cells"] = adata_rna.obs_names.astype(str)
adata_rna.var["genes"] = adata_rna.var_names.astype(str)
genes = adata_rna.var.query("highly_variable").index
peaks = adata_atac.var.query("highly_variable").index
features = pd.Index(np.concatenate([adata_rna.var_names, adata_atac.var_names]))
feature_embeddings = np.concatenate([adata_rna.varm["X_glue"], adata_atac.varm["X_glue"]])
skeleton = guidance_hvf.edge_subgraph(e for e, attr in dict(guidance_hvf.edges).items()if attr["type"] == "fwd").copy()
reginf = scglue.genomics.regulatory_inference(features, feature_embeddings,skeleton=skeleton, random_state=0)
gene2peak = reginf.edge_subgraph(e for e, attr in dict(reginf.edges).items()if attr["qval"] < 0.05).copy()
nx.write_graphml(gene2peak, "gene2peak.graphml.gz")
motif_bed = scglue.genomics.read_bed("/COMMA/reference_genome/GENCODE/unused/JASPAR2022-hg38.bed.gz")
tfs = pd.Index(motif_bed["name"].unique()).intersection(adata_rna.var_names)
grn = scglue.genomics.build_tf_target_network(adata_rna,adata_atac,gene2peak,genes,peaks,tfs,motif_bed)

print("GRN nodes:", grn.number_of_nodes(), "GRN edges:", grn.number_of_edges())
if grn.number_of_edges() == 0:
    raise ValueError("Empty GRN returned.")

nx.write_graphml(grn, "vsmc_tf_network.graphml")