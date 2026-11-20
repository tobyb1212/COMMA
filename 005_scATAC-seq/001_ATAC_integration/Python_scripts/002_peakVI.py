from scvi.model import PEAKVI as PeakVI
import anndata as ad
import scanpy as sc
import os
import glob
import random
import matplotlib.pyplot as plt

os.chdir("/COMMA/results/scatac/matrices")

# Read AnnData Object of Peak Matrices & Concatenate
adatas = {}
file_pattern = os.path.join("/COMMA/results/scatac/matrices/peaks/", "**", "*_PeakMatrix.h5ad")
files = glob.glob(file_pattern, recursive=True)
for file_path in files:
    sample_adata = sc.read_h5ad(file_path)
    sample_id = os.path.basename(file_path).split("_PeakMatrix.h5ad")[0]
    adatas[sample_id] = sample_adata

adata = ad.concat(adatas, label="fastq_name_new")
adata.obs_names_make_unique()
adata.obs.rename(columns={'Sample.Name': 'Sample_Name'}, inplace=True)

# Filter Peaks Present in Less than 1% Cells
min_cells = int(adata.shape[0] * 0.01)
sc.pp.filter_genes(adata, min_cells=min_cells)

# Select Highly Variable Peaks
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=80000, batch_key="Study")
adata_hvg = adata[:, adata.var.highly_variable].copy()

# Setup AnnData for PeakVI for Correcting for Sample Name and Sequencing Method
PeakVI.setup_anndata(adata_hvg, batch_key="Study", categorical_covariate_keys =["Sample_Name", "sequencing_method"], continuous_covariate_keys=["nFrags", "TSSEnrichment"])
model = PeakVI(adata_hvg)
model.train()

# Transfer Latent Space Embeddings and save PeakVI Model
model.save('peak_hvg_model_initial.pt', overwrite=True)
model = PeakVI.load("peak_hvg_model_initial.pt", adata=adata_hvg)
latent = model.get_latent_representation()
adata.obsm["X_PeakVI"] = latent

# Compute UMAP Coordinates for Cells
sc.pp.neighbors(adata, use_rep="X_PeakVI")
sc.tl.umap(adata, min_dist=0.2)
sc.tl.leiden(adata, resolution=0.8, key_added=f'leiden_0.8')
sc.pl.umap(adata, color=["Study", "sequencing_method", "leiden_0.8"], ncols=2, save="Peak_final.png")

# Supplemental Figure 10D - all-study PeakVI UMAP before Leblanc exclusion
study_colors = {"Amrute 2024": "#f57e20",
    "Hocker 2021": "#dbdb8c",
    "Kanemaru 2023": "#9372b2",
    "Leblanc 2025": "#9fdae5",
    "Read 2024": "#ae4a4b",
    "Selewa 2023": "#8c6d30",}
adata.uns["Study_colors"] = [study_colors[x] for x in adata.obs["Study"].cat.categories]
col = "Study"
plt.figure(figsize=(6, 4))
sc.pl.umap(adata, color=col, show=False, title=None)
plt.savefig(f"pre_umap.png", dpi=900, bbox_inches="tight")
plt.close()
sc.pl.umap(adata[adata.obs.groupby(col, group_keys=False).apply(lambda x: x.sample(1, random_state=0)).index], color=col, save=f"pre_labels.pdf")

# Remove Leblanc Data (Poorly Integrated)
adata_filtered = adata[~adata.obs['Study'].isin(['Leblanc 2025'])].copy()
sc.pp.highly_variable_genes(adata_filtered, flavor="seurat_v3", n_top_genes=80000, batch_key="Study")
adata_hvg = adata_filtered[:, adata_filtered.var.highly_variable].copy()
# Setup AnnData for PeakVI for Correcting for Sample Name and Sequencing Method
PeakVI.setup_anndata(adata_hvg, batch_key="Study", categorical_covariate_keys =["Sample_Name", "sequencing_method"], continuous_covariate_keys=["nFrags", "TSSEnrichment"])
model = PeakVI(adata_hvg)
model.train()
# Transfer Latent Space Embeddings and save PeakVI Model
model.save('peak_hvg_model_post_leblanc.pt', overwrite=True)
model = PeakVI.load("peak_hvg_model_post_leblanc.pt", adata=adata_hvg)
latent = model.get_latent_representation()
adata_filtered.obsm["X_PeakVI"] = latent
sc.pp.neighbors(adata_filtered, use_rep="X_PeakVI")
sc.tl.umap(adata_filtered, min_dist=0.2)
sc.tl.leiden(adata_filtered, resolution=0.8, key_added=f'leiden_0.8')
sc.pl.umap(adata_filtered, color=["Study", "sequencing_method", "leiden_0.8"], ncols=2, save="Peak_final_filtered.png")

# Save AnnData Objects
adata_filtered_save = adata_filtered[:, 0].copy() # Don't need .var, just .obsm and .obs
adata_filtered_save.write_h5ad("/COMMA/results/scatac/matrices/Peak_matrix_final_integrated_filtered.h5ad")


