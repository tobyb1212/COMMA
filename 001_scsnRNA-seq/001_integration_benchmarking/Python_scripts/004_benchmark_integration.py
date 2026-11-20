import scib
import pandas as pd
import numpy as np
import scanpy as sc
import os
import argparse
import random

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process some directories.')
    parser.add_argument('-f', '--raw_path', required=True, help='Path to non-integrated .h5ad file')
    parser.add_argument('-i', '--integrated_file', required=True, help='Path to integrated file')
    return parser.parse_args()

def main():
    args = parse_arguments()
    random.seed(1234)
    adata_path = args.raw_path
    integrated_path = args.integrated_file

    integrated_adata=sc.read_h5ad(integrated_path)
    adata = sc.read_h5ad(adata_path)

    out_dir = os.path.dirname(os.path.abspath(integrated_path))
    os.chdir(out_dir)

    # Transfer HCA CellTypist cell type predictions to unintegrated object
    # Remove cells from unintegrated object that have been removed as doublet clusters
    adata = adata[integrated_adata.obs_names].copy()
    adata.obs['HCA_celltypist_predicted_labels'] = integrated_adata.obs['HCA_celltypist_predicted_labels'].values

    # Rename adata.obs to not break scib
    adata.obs.rename(columns={'Sample Name': 'Sample_Name'}, inplace=True)
    integrated_adata.obs.rename(columns={'Sample Name': 'Sample_Name'}, inplace=True)

    # Downsample to 500k cells
    common_cells = adata.obs_names.intersection(integrated_adata.obs_names)
    downsampled_cells = np.random.choice(common_cells, size=500000, replace=False)
    adata = adata[downsampled_cells].copy()
    integrated_adata = integrated_adata[downsampled_cells].copy()

    # Addressing some minor issues with code
    method = '_'.join(os.path.basename(integrated_path).split('_')[:-2])
    if "X_scanorama" in integrated_adata.obsm:
        integrated_adata.obsm[f"X_{method}"] = integrated_adata.obsm.pop("X_scanorama")
    elif "X_harmony" in integrated_adata.obsm:
        integrated_adata.obsm[f"X_{method}"] = integrated_adata.obsm.pop("X_harmony")
    elif "X_scANVI" in integrated_adata.obsm:
        integrated_adata.obsm[f"X_{method}"] = integrated_adata.obsm.pop("X_scANVI")
    elif "X_scPoli" in integrated_adata.obsm:
        integrated_adata.obsm[f"X_{method}"] = integrated_adata.obsm.pop("X_scPoli")

    # Change to embedding output to X_emb to match what scib requires and using knn for BBKNN
    if method == "BBKNN":
        type_ = "knn"
        embed = "X_pca"
    else:
        type_ = "embed"
        embed = f"X_{method}"
        integrated_adata.obsm["X_emb"] = integrated_adata.obsm[f"X_{method}"]

    metrics = scib.metrics.metrics(adata, integrated_adata, batch_key='Sample_Name',
    label_key='HCA_celltypist_predicted_labels', embed=embed, cluster_key = 'scib_cluster', cluster_nmi = 'cluster_nmi', 
    nmi_method= "arithmetic", si_metric ='euclidean', organism = 'human', subsample=30, type_=type_,
    trajectory_=False, ari_=True, nmi_=True, pcr_=True, cell_cycle_=True, hvg_score_=False, isolated_labels_=False, graph_conn_=True,
    kBET_=True, lisi_graph_=True, silhouette_=True)

    # Change to the correct directory and save metrics
    outpath = os.path.join(out_dir, f"{method}_metrics.csv")
    metrics.to_csv(outpath)

if __name__ == "__main__":
    main()