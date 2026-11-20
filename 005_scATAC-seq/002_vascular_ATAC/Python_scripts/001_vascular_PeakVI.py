from scvi.model import PEAKVI as PeakVI
import anndata as ad
import scanpy as sc
import os
import glob
import matplotlib.pyplot as plt

os.chdir("/COMMA/results/scatac/matrices/vascular")

# Read AnnData Object of Peak Matrices & Concatenate
adatas = {}
file_pattern = os.path.join("/COMMA/results/scatac/matrices/vascular/peaks/", "**", "*_PeakMatrix.h5ad")
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
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=60000, batch_key="Study")
adata_hvg = adata[:, adata.var.highly_variable].copy()

# Setup AnnData for PeakVI for Correcting for Sample Name and Sequencing Method
PeakVI.setup_anndata(adata_hvg, batch_key="Study", categorical_covariate_keys =["Sample_Name", "sequencing_method"], continuous_covariate_keys=["nFrags", "TSSEnrichment"])
model = PeakVI(adata_hvg)
model.train()
model.save('peak_hvg_model.pt', overwrite=True)

# Transfer Latent Space Embeddings and save PeakVI Model
#model = PeakVI.load("peak_hvg_model.pt", adata=adata_hvg)
latent = model.get_latent_representation()
adata.obsm["X_PeakVI"] = latent

# Compute UMAP Coordinates for Cells
sc.pp.neighbors(adata, use_rep="X_PeakVI")

#Find Clusters & plot umaps
for dist in [0.05, 0.1, 0.2, 0.3, 0.4]:
    sc.tl.umap(adata, min_dist=dist)
    for res in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
            sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}')
            sc.pl.umap(adata, color=[f"leiden_{res}"], save=f"mindist_{dist}_leiden_{res}.png")

sc.tl.umap(adata, min_dist=0.2)
sc.tl.leiden(adata, resolution=0.8, key_added=f'leiden_0.8')
sc.pl.umap(adata, color=["Study", "sequencing_method", "leiden_0.8"], ncols=2, save="Peak_final.png")

# Supplemental Figure 10E
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
plt.savefig(f"peak_umap.png", dpi=900, bbox_inches="tight")
plt.close()
sc.pl.umap(adata[adata.obs.groupby(col, group_keys=False).apply(lambda x: x.sample(1, random_state=0)).index], color=col, save=f"peak_labels.pdf")

# Save AnnData Objects
adata_save = adata[:, 0].copy() # Don't need .var, just .obsm and .obs
adata_save.write_h5ad("/COMMA/results/scatac/matrices/vascular/Peak_matrix_integrated.h5ad")

# Filter Selewa et al Data
adata_filtered = adata[~adata.obs['Study'].isin(['Selewa 2023'])].copy()
sc.tl.umap(adata_filtered, min_dist=0.2)
sc.tl.leiden(adata_filtered, resolution=0.8, key_added=f'leiden_0.8')
sc.pl.umap(adata_filtered, color=["Study", "sequencing_method", "leiden_0.8"], ncols=2, save="filtered_Peak_final.png")

# Supplemental Figure 10F
study_colors = {"Amrute 2024": "#f57e20",
    "Hocker 2021": "#dbdb8c",
    "Kanemaru 2023": "#9372b2",
    "Leblanc 2025": "#9fdae5",
    "Read 2024": "#ae4a4b",
    "Selewa 2023": "#8c6d30",}
adata_filtered.uns["Study_colors"] = [study_colors[x] for x in adata_filtered.obs["Study"].cat.categories]
col = "Study"
plt.figure(figsize=(6, 4))
sc.pl.umap(adata_filtered, color=col, show=False, title=None)
plt.savefig(f"final_umap.png", dpi=900, bbox_inches="tight")
plt.close()
sc.pl.umap(adata_filtered[adata_filtered.obs.groupby(col, group_keys=False).apply(lambda x: x.sample(1, random_state=0)).index], color=col, save=f"final_labels.pdf")

# Save AnnData Objects
adata_save = adata_filtered[:, 0].copy() # Don't need .var, just .obsm and .obs
adata_save.write_h5ad("/COMMA/results/scatac/matrices/vascular/filtered_Peak_matrix_integrated.h5ad")