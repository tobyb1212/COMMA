import pandas as pd
import numpy as np
import scanpy as sc
import os
import argparse
import random

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process some directories.')
    parser.add_argument('-f', '--integrated_file', required=True, help='Path to the raw file')
    return parser.parse_args()

def main():
    args = parse_arguments()
    integrated_path = args.integrated_file
    method = '_'.join(os.path.basename(integrated_path).split('_')[:-1]) # Extract integration method from filename
    out_dir = os.path.dirname(os.path.abspath(integrated_path))

    adata = sc.read_h5ad(integrated_path)
    sc.tl.leiden(adata, resolution=3, key_added=f'leiden_{method}') # Leiden clustering of cells

    # Finding clusters where 40% of the cells have a cluster score of >0.2 - marking as doublet cluster
    adata.obs['doublet_cluster'] = adata.obs[f'leiden_{method}'].map(adata.obs.groupby(f'leiden_{method}')["scrublet_score"]
        .apply(lambda scores: (scores > 0.2).sum() / len(scores) > 0.30))

    # Plot UMAPs with doublet score and list of doublets
    os.chdir(out_dir)
    sc.pl.umap(adata, color=["doublet_cluster", "scrublet_score"], ncols=2, save="doublet_umaps.png")

    # Filter doublet clusters from adata
    adata_new = adata[adata.obs['doublet_cluster'] == False].copy()

    # Plot UMAP for doublet clusters and doublet scores
    sc.pl.umap(adata_new, color=["doublet_cluster", "scrublet_score"], ncols=2, save="doublet_umaps_filtered.png")

    # Save filtered adata file
    outpath = os.path.join(out_dir, f"{method}_doublet_removed.h5ad")
    adata_new.write_h5ad(outpath)
    print(f"Saved file to {outpath}")

if __name__ == "__main__":
    main()
