import scvi
import pandas as pd
import numpy as np
import scanpy as sc
import os
import argparse
import random

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process some directories.')
    parser.add_argument('-r', '--raw_file', required=True, help='Path to the raw file')
    parser.add_argument('-v', '--hvg_file', required=True, help='Path to the HVG file')
    return parser.parse_args()

def main():
    args = parse_arguments()
    adata_path = args.raw_file
    hvg_path = args.hvg_file

    directory="scANVI"
    if not os.path.isdir(directory):
        os.mkdir(directory)
    os.chdir(directory)

    random.seed(1234) 

    adata_integrated = sc.read_h5ad(adata_path)
    adata_hvg = sc.read_h5ad(hvg_path)
    # Load scVI Model
    vae = scvi.model.SCVI.load("/COMMA/results/scrna/integrated/scVI/model.pt/", adata = adata_hvg)
    # Load HCA Broad Cell Type Annotations
    df = pd.read_csv("/COMMA/results/scrna/metadata/HCA_celltypes.csv")
    celltype_dict = dict(zip(df['cell_name'], df['cell_type']))
    adata_hvg.obs['HCA_cell_type'] = adata_hvg.obs_names.map(celltype_dict)
    adata_hvg.obs['HCA_cell_type'] = adata_hvg.obs['HCA_cell_type'].fillna('Unknown')
    adata_hvg.obs['HCA_cell_type'] = adata_hvg.obs['HCA_cell_type'].replace({'Ventricular Cardiomyocyte': 'Cardiomyocyte','Atrial Cardiomyocyte': 'Cardiomyocyte'})

    # Establishing & Training Model with Partial Broad Cell Type Annotations from Heart Cell Atlas
    lvae = scvi.model.SCANVI.from_scvi_model(vae, adata=adata_hvg, labels_key="HCA_cell_type", unlabeled_category="Unknown")
    lvae.train(max_epochs=20, n_samples_per_label=100)
    adata_hvg.obsm["X_scANVI"] = lvae.get_latent_representation()

    # Transfer Annotations to Full Object
    adata_integrated.obsm["X_scANVI"] = adata_hvg.obsm["X_scANVI"]

    # Plotting Integrated UMAP to Assess Integration Success
    sc.pp.neighbors(adata_integrated, use_rep="X_scANVI")
    sc.tl.umap(adata_integrated)
    sc.pl.umap(adata_integrated, color=["Study", "Patient", "Sample Name", "sequencing_method", "sample_type", "Region", "total_counts", "pct_counts_mt", "pct_counts_ribo"], ncols=2,  save="newumaps.png")

    # Saving integrated adata object
    adata_integrated.write("scANVInew_integrated.h5ad")

if __name__ == "__main__":
    main()

