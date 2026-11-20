import pandas as pd
import numpy as np
import scvi
import scanpy as sc
from scvi.external import SysVI
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

    directory="sysVI"
    if not os.path.isdir(directory):
        os.mkdir(directory)
    os.chdir(directory)

    random.seed(1234)

    adata_integrated = sc.read_h5ad(adata_path)
    adata_hvg = sc.read_h5ad(hvg_path)
    scvi.settings.seed = 0

    # Establish VAE Model with Primary Batch Effect of the Sample, but accounting for other categorical and clinical covariates
    SysVI.setup_anndata(adata=adata_hvg, batch_key="Study", categorical_covariate_keys=["Sample Name", "sample_type"])    # These options differ from default parameters to increase dimensionality of latent space and reflect differences in gene expression across batches
    model = SysVI(adata=adata_hvg, embed_categorical_covariates=True, n_latent=30, n_layers=2)
    model.train(plan_kwargs={"kl_weight": 0.2,"z_distance_cycle_weight": 10}, max_epochs=100)
    model.save(f'model.pt',overwrite=True)

    adata_hvg.obsm["X_sysVI"] = model.get_latent_representation()
    adata_integrated.obsm["X_sysVI"] = adata_hvg.obsm["X_sysVI"]

    # Plotting Integrated UMAP to Assess Integration Success
    sc.pp.neighbors(adata_integrated, use_rep="X_sysVI")
    sc.tl.umap(adata_integrated)
    sc.pl.umap(adata_integrated, color=["Study", "Patient", "Sample Name", "sequencing_method", "sample_type", "Region", "total_counts", "pct_counts_mt", "pct_counts_ribo"], ncols=2, save="umaps.png")

    # Saving integrated adata object
    adata_integrated.write("sysVI_integrated.h5ad")

    print("Performed sysVI Integration")

if __name__ == "__main__":
    main()
