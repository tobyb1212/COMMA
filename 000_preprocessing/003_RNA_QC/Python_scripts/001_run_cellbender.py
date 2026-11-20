import pandas as pd
import anndata as ad
import scanpy as sc
import os
import argparse
import random

random.seed(1234)

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process some directories.')
    parser.add_argument('-f', '--input_folder', type=str, required=True, help='Input directory containing matrix.mtz.gz, barcodes.tsv.gz, features.tsv.gz')
    return parser.parse_args()

def main():
    args = parse_arguments()
    input_folder = args.input_folder
    # Change Directory and Read gzipped output from STARsolo
    os.chdir(input_folder)
    adata = sc.read_10x_mtx(input_folder)
    # Make cellbender directory
    parent_dir = os.path.dirname(os.path.abspath(input_folder))
    cellbender_dir = os.path.join(parent_dir, 'cellbender')
    os.makedirs(cellbender_dir, exist_ok=True)
    # Write count matrix to .h5ad for cellbender
    outpath = os.path.join(cellbender_dir, 'raw_feature_barcode_matrix.h5ad')
    adata.write_h5ad(outpath)

if __name__ == "__main__":
    main()
