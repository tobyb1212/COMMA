import scib
import pandas as pd
import numpy as np
import scanpy as sc
import os
import random
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process some directories.')
    parser.add_argument('-f', '--raw_file', required=True, help='Path to the raw file')
    return parser.parse_args()

def main():
    args = parse_arguments()
    adata_path = args.raw_file

    directory="/COMMA/results/scrna/integrated/harmony"
    os.chdir(directory)

    random.seed(1234)

    adata_integrated = sc.read_h5ad(adata_path)
    adata_hvg = sc.read_h5ad("harmony_integrated_hvg.h5ad")

    # Transfer embeddings from HVG object to full object
    adata_integrated.obsm["X_harmony"] = adata_hvg.obsm["X_x_harmony"]

    # Plotting Integrated UMAP to Assess Integration Success
    sc.pp.neighbors(adata_integrated, use_rep="X_harmony")
    sc.tl.umap(adata_integrated)
    sc.pl.umap(adata_integrated, color=["Study", "Patient", "Sample Name", "sequencing_method", "sample_type", "Region", "total_counts", "pct_counts_mt", "pct_counts_ribo"], ncols=2, save = "umaps.png")

    # Save integrated adata object
    adata_integrated.write("harmony_integrated.h5ad")

    print("Evaluated Harmony Integration with Covariates Integration")

if __name__ == "__main__":
    main()

