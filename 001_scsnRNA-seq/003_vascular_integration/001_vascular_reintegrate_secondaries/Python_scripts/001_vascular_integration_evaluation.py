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

out_dir = "/COMMA/results/scrna/downstream_analysis/subclustering/vascular"
os.chdir(out_dir)

def main():
    args = parse_arguments()
    random.seed(1234)
    adata_path = args.raw_path
    integrated_path = args.integrated_file

    adata_integrated=sc.read_h5ad(integrated_path)
    adata = sc.read_h5ad(adata_path)

    method = '_'.join(os.path.basename(integrated_path).split('_')[:-2])

    sc.pp.neighbors(adata_integrated, use_rep=f"X_{method}")
    sc.tl.umap(adata_integrated, min_dist = 0.2)
    sc.pl.umap(adata_integrated, color=["Study", "sequencing_method", "sample_type"], ncols=2, save=f"s_{method}.png") 
    for res in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
        sc.tl.leiden(adata_integrated, resolution=res, key_added=f'leiden_{method}_{res}')
        sc.pl.umap(adata_integrated, color=[f"leiden_{method}_{res}"], save=f"_{method}_leiden_{res}.png")
    
    adata.obs.rename(columns={'Sample Name': 'Sample_Name'}, inplace=True)
    adata_integrated.obs.rename(columns={'Sample Name': 'Sample_Name'}, inplace=True)

    embed = f"X_{method}"
    type_ = "embed"
    adata_integrated.obsm["X_emb"] = adata_integrated.obsm[f"X_{method}"]

    metrics = scib.metrics.metrics(adata, adata_integrated, batch_key='Sample_Name',
        label_key='broad_cell_type', embed=embed, cluster_key = 'scib_cluster', cluster_nmi = 'cluster_nmi', 
        nmi_method= "arithmetic", si_metric ='euclidean', organism = 'human', subsample=30, type_=type_,
        trajectory_=False, ari_=True, nmi_=True, pcr_=True, cell_cycle_=True, hvg_score_=False, isolated_labels_=False, graph_conn_=True,
        kBET_=True, lisi_graph_=True, silhouette_=True)

    outpath = os.path.join(out_dir, f"{method}_metrics.csv")
    metrics.to_csv(outpath)

if __name__ == "__main__":
    main()