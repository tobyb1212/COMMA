import scib
import pandas as pd
import numpy as np
import scanpy as sc
import scanpy.external as sce
import scanorama
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

    directory="scanorama"
    if not os.path.isdir(directory):
        os.mkdir(directory)
    os.chdir(directory)

    random.seed(1234)

    adata_integrated = sc.read_h5ad(adata_path)
    adata_hvg = sc.read_h5ad(hvg_path)

    # Scaling data and re-performing Dimensionality Reduction for adata_hvg
    # Scanorama seems to work best with scaled data - https://www.nature.com/articles/s41592-021-01336-8
    sc.pp.scale(adata_hvg)
    sc.tl.pca(adata_hvg)

    # Performing Integration on HVGs
    adata_hvg = adata_hvg[adata_hvg.obs['Sample Name'].argsort()]
    sce.pp.scanorama_integrate(adata_hvg, 'Sample Name', basis='X_pca', adjusted_basis='X_scanorama', batch_size=10000) 
    adata_hvg = adata_hvg[adata_integrated.obs_names]

    # Transferring Embeddings to Complete adata
    adata_integrated.obsm["X_scanorama"] = adata_hvg.obsm["X_scanorama"]

    # Plotting Integrated UMAP to Assess Integration Success
    sc.pp.neighbors(adata_integrated, use_rep="X_scanorama")
    sc.tl.umap(adata_integrated)
    sc.pl.umap(adata_integrated, color=["Study", "Patient", "Sample Name", "sequencing_method", "sample_type", "Region", "total_counts", "pct_counts_mt", "pct_counts_ribo"], ncols=2, save="umaps.png")

    # Saving integrated adata object
    adata_integrated.write("Scanorama_integrated.h5ad")

    print("Performed scanorama Integration")

if __name__ == "__main__":
    main()