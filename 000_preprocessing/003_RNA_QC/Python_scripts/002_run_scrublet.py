import pandas as pd
import anndata as ad
import scanpy as sc
import os
import argparse
import random
import scrublet as scr

random.seed(1234)

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process some directories.')
    parser.add_argument('-f', '--input_file', type=str, required=True, help='Input .h5 output from cellbender')
    return parser.parse_args()

def main():
    args = parse_arguments()
    input_file = args.input_file

    # Read Cellbender Input
    adata = sc.read_10x_h5(input_file)

    # Make cellbender directory
    out_dir = os.path.dirname(os.path.abspath(input_file))

    # Run scrublet with standard settings - do not include automatically thresholded doublet prediction, manual threshold set later
    scrub = scr.Scrublet(adata.X)
    doublet_scores, automatic_predicted_doublets = scrub.scrub_doublets()
    adata.obs['scrublet_score'] = doublet_scores

    # Write count matrix to .h5ad for cellbender
    outpath = os.path.join(out_dir, 'scrublet.h5ad')
    adata.write_h5ad(outpath)

if __name__ == "__main__":
    main()


