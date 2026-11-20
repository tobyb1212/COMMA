import scib
import pandas as pd
import numpy as np
import scvi
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

    directory="scVI"
    if not os.path.isdir(directory):
        os.mkdir(directory)
    os.chdir(directory)

    random.seed(1234)

    adata_integrated = sc.read_h5ad(adata_path)
    adata_hvg = sc.read_h5ad(hvg_path)

    # Establish VAE Model with Primary Batch Effect of the Sample, but accounting for other categorical and clinical covariates
    scvi.model.SCVI.setup_anndata(adata_hvg, layer="counts",batch_key="Study", categorical_covariate_keys=["Sample Name", "sample_type"],continuous_covariate_keys=['total_counts',"pct_counts_mt", "pct_counts_ribo"])
    # These options differ from default parameters to increase dimensionality of latent space and reflect differences in gene expression across batches
    model = scvi.model.SCVI(adata_hvg, gene_likelihood="nb", dispersion ='gene-batch', n_latent=30, n_layers=2)
    model.view_anndata_setup(adata_hvg)
    model.train()
    model.save(f'model.pt',overwrite=True)

    adata_hvg.obsm["X_scVI"] = model.get_latent_representation()
    adata_integrated.obsm["X_scVI"] = adata_hvg.obsm["X_scVI"]

    # Plotting Integrated UMAP to Assess Integration Success
    sc.pp.neighbors(adata_integrated, use_rep="X_scVI")
    sc.tl.umap(adata_integrated)
    sc.pl.umap(adata_integrated, color=["Study", "Patient", "Sample Name", "sequencing_method", "sample_type", "Region", "total_counts", "pct_counts_mt", "pct_counts_ribo"], ncols=2, save="umaps.png")

    # Saving integrated adata object
    adata_integrated.write("scVI_integrated.h5ad")

    print("Performed SCVI Integration")

if __name__ == "__main__":
    main()
