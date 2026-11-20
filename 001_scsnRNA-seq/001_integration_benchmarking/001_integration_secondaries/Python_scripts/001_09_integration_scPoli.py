import numpy as np
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import classification_report
from scarches.models.scpoli import scPoli
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

    directory="scPoli"
    if not os.path.isdir(directory):
        os.mkdir(directory)
    os.chdir(directory)

    random.seed(1234)

    adata_integrated = sc.read_h5ad(adata_path)
    adata_hvg = sc.read_h5ad(hvg_path)
    adata_hvg.X = adata_hvg.layers["counts"]

    # Set parameters for scPoli Model
    early_stopping_kwargs = { "early_stopping_metric": "val_prototype_loss","mode": "min","threshold": 0,"patience": 20,
    "reduce_lr": True,"lr_patience": 13,"lr_factor": 0.1}
    condition_keys = ['Study', 'Sample Name', "sample_type"]
    cell_type_key = 'HCA_cell_type'

    # Load HCA Annotations
    df = pd.read_csv("/COMMA/results/scrna/metadata/HCA_celltypes.csv")
    celltype_dict = dict(zip(df['cell_name'], df['cell_type']))
    adata_hvg.obs['HCA_cell_type'] = adata_hvg.obs_names.map(celltype_dict)
    adata_hvg.obs['HCA_cell_type'] = adata_hvg.obs['HCA_cell_type'].fillna('Unknown')
    adata_hvg.obs['HCA_cell_type'] = adata_hvg.obs['HCA_cell_type'].replace({'Ventricular Cardiomyocyte': 'Cardiomyocyte','Atrial Cardiomyocyte': 'Cardiomyocyte'})

    # Split adata into reference with known annotations
    adata_ref = adata_hvg[adata_hvg.obs['HCA_cell_type'] != 'Unknown'].copy()
    
    # Train model on Reference Anndata
    scpoli_model = scPoli(adata=adata_ref,condition_keys=condition_keys,cell_type_keys=cell_type_key,embedding_dims=5,recon_loss='nb',)
    scpoli_model.train(n_epochs=50,pretraining_epochs=40,early_stopping_kwargs=early_stopping_kwargs,eta=5)
    scpoli_model.save("reference_model.pt")

    # Transfer the learning function to all new data
    labeled_indices = adata_hvg.obs.index[adata_hvg.obs['HCA_cell_type'] != 'Unknown']
    labeled_indices = [adata_hvg.obs.index.get_loc(idx) for idx in labeled_indices]
    scpoli_query = scPoli.load_query_data(adata=adata_hvg,reference_model=scpoli_model, labeled_indices=labeled_indices)
    scpoli_query.train(n_epochs=50,pretraining_epochs=40,eta=10)
    scpoli_query.save("query_model.pt")
    
    # Get Latent Embeddings from scPoli Model
    for key in condition_keys:
        adata_hvg.obs[key] = np.array(adata_hvg.obs[key])
    adata_hvg.obsm['X_scPoli'] = scpoli_query.get_latent(adata_hvg,mean=True)

    # Transfer Annotations to Full Object
    adata_integrated.obsm["X_scPoli"] = adata_hvg.obsm["X_scPoli"]

    # Plotting Integrated UMAP to Assess Integration Success
    sc.pp.neighbors(adata_integrated, use_rep="X_scPoli")
    sc.tl.umap(adata_integrated)
    sc.pl.umap(adata_integrated, color=["Study", "Patient", "Sample Name", "sequencing_method", "sample_type", "Region", "total_counts", "pct_counts_mt", "pct_counts_ribo"], ncols=2,  save="newumaps.png")

    # Saving integrated adata object
    adata_integrated.write("scPolinew_integrated.h5ad")

if __name__ == "__main__":
    main()
