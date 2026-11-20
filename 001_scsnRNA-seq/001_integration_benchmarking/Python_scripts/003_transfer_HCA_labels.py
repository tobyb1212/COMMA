import celltypist
import scanpy as sc
from celltypist import models
import argparse
import pandas as pd
import os
import re

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process some directories.')
    parser.add_argument('-f', '--integrated_file', required=True, help='Path to the raw file')
    return parser.parse_args()

def main():
    args = parse_arguments()
    integrated_path = args.integrated_file
    adata = sc.read_h5ad(integrated_path)
    method = '_'.join(os.path.basename(integrated_path).split('_')[:-2]) # Extract integration method from filename
    
    # Read in pre-existing broad HCA annotations to train CellTypist model & Select HVGs
    adata_hca = adata[adata.obs["Study"].isin(["Litvinukova 2020", "Kanemaru 2023"]), :]
    df = pd.read_csv("/COMMA/results/scrna/metadata/HCA_celltypes.csv")
    celltype_dict = dict(zip(df['cell_name'], df['cell_type']))
    adata_hca.obs['HCA_cell_type'] = adata_hca.obs_names.map(celltype_dict)
    adata_hca.obs['HCA_cell_type'] = adata_hca.obs['HCA_cell_type'].fillna('Unknown')
    adata_hca.obs['HCA_cell_type'] = adata_hca.obs['HCA_cell_type'].replace({'Ventricular Cardiomyocyte': 'Cardiomyocyte','Atrial Cardiomyocyte': 'Cardiomyocyte'})
    adata_hca = adata_hca[adata_hca.obs['HCA_cell_type'] != 'Unknown', :]
    sc.pp.highly_variable_genes(adata_hca, flavor="seurat_v3", n_top_genes=2000, layer="counts", batch_key="Study")
    adata_hca_hvg = adata_hca[:, adata_hca.var.highly_variable].copy()
    adata_hca_hvg.X = adata_hca_hvg.layers['counts']
    sc.pp.normalize_total(adata_hca_hvg, target_sum=1e4)
    sc.pp.log1p(adata_hca_hvg)

    # Train & save cell typist model
    model = celltypist.train(adata_hca_hvg, labels = 'HCA_cell_type', n_jobs = 1, feature_selection = False)
    model_path = os.path.join(os.path.dirname(os.path.abspath(integrated_path)), f'HCA_celltypes_{method}.pkl')
    model.write(model_path)

    # Predict cell types based on celltypist model from HCA data
    adata_hvg = adata[:, adata_hca.var.highly_variable].copy()
    adata_hvg.X = adata_hvg.layers['counts']
    sc.pp.normalize_total(adata_hvg, target_sum=1e4)
    sc.pp.log1p(adata_hvg)
    predictions = celltypist.annotate(adata_hvg, model = f'/COMMA/003_rna_integration/models/HCA_celltypes_{method}.pkl', majority_voting = True)
    adata_hvg = predictions.to_adata(prefix='HCA_celltypist_')
    matching_cols = [col for col in adata_hvg.obs.columns if re.match(r"^HCA_celltypist_", col)]
    for col in matching_cols:
        adata.obs[col] = adata_hvg.obs[col].reindex(adata.obs.index)

    # Label Leiden clusters with majority predictions from the model
    cluster_labels = (adata.obs.groupby(f'leiden_{method}')["HCA_celltypist_predicted_labels"].agg(lambda x: x.value_counts().idxmax()))
    adata.obs[f"annotated_leiden_{method}"] = adata.obs[f'leiden_{method}'].map(cluster_labels)

    # Plot predicted cell types on UMAP
    out_dir = os.path.dirname(os.path.abspath(integrated_path))
    os.chdir(out_dir)
    sc.pl.umap(adata, color = ['HCA_celltypist_predicted_labels', 'HCA_celltypist_majority_voting', f"annotated_leiden_{method}"],  ncols=2, save = "celltypist_umap.png")

    # Write celltypist predicted .h5ad object
    outpath = os.path.join(out_dir, f"{method}_celltypist_predicted.h5ad")
    adata.write_h5ad(outpath)

if __name__ == "__main__":
    main()

