import scib
import pandas as pd
import numpy as np
import scanpy as sc
import scanpy.external as sce
import random
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

    directory="BBKNN"
    if not os.path.isdir(directory):
        os.mkdir(directory)
    os.chdir(directory)

    random.seed(1234)

    adata_integrated = sc.read_h5ad(adata_path)
    adata_hvg = sc.read_h5ad(hvg_path)

    # Perform BBKNN on HVGs and transfer connectivities to full object
    sce.pp.bbknn(adata_hvg, batch_key='Sample Name', use_rep='X_pca', neighbors_within_batch=1)
    adata_integrated.obsp["connectivities"] = adata_hvg.obsp["connectivities"]
    adata_integrated.obsp["distances"] = adata_hvg.obsp["distances"]
    adata_integrated.uns['neighbors'] = adata_hvg.uns['neighbors'] 

    # Compute UMAP on full dataset
    sc.tl.umap(adata_integrated)
    sc.pl.umap(adata_integrated, color=["Study", "Patient", "Sample Name", "sequencing_method", "sample_type", "Region", "total_counts", "pct_counts_mt", "pct_counts_ribo"], ncols=2, save="umaps.png")

    # Saving integrated adata object
    adata_integrated.write("BBKNN_integrated.h5ad")

    print("Integrated successfully")

if __name__ == "__main__":
    main()